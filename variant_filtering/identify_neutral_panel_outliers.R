#!/usr/bin/env Rscript
# Identify PCAdapt and OutFLANK outliers and merge their SNP positions.
# Run before build_neutral_panel.sh.
# Software: R 4.4.1, PLINK2 2.00a4.3, pcadapt 4.3.5, OutFLANK 0.2,
#           qvalue and vcfR.
# Usage: Rscript identify_neutral_panel_outliers.R [VCF] [K5_TABLE] [LOCATION_TABLE] [OUTPUT_DIR]
# Tables have two columns without a header: sample ID and group.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 4L) {
    stop("Usage: Rscript identify_neutral_panel_outliers.R [VCF] [K5_TABLE] [LOCATION_TABLE] [OUTPUT_DIR]")
}
defaults <- c("full.panel_5bpind.morePAN.HWE.DPqual.AC3.bi.SOR.vcf.gz",
              "K5_morePAN_n46_2.col.txt", "samp.locs_FLsep_morePAN_n46_2.col.txt",
              "outlier_lists")
defaults[seq_along(args)] <- args
vcf_path <- defaults[1]
k5_path <- defaults[2]
locations_path <- defaults[3]
out_dir <- defaults[4]

for (pkg in c("pcadapt", "qvalue", "vcfR", "OutFLANK")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
}
if (Sys.which("plink2") == "") stop("Load PLINK2 before running this script.")
if (!all(file.exists(defaults[1:3]))) stop("A required VCF or grouping table is missing.")
if (file.exists(out_dir)) stop("Use a new output directory: ", out_dir)
if (!dir.create(out_dir, recursive = TRUE)) stop("Cannot create output directory.")

read_groups <- function(path, expected_groups) {
    groups <- read.table(path, header = FALSE, sep = "\t", colClasses = "character",
                         comment.char = "", quote = "", stringsAsFactors = FALSE)
    if (ncol(groups) != 2L || anyNA(groups) || any(groups == "")) {
        stop("Grouping table must have two nonempty columns: ", path)
    }
    names(groups) <- c("sample", "group")
    if (anyDuplicated(groups$sample)) stop("Duplicate sample IDs in: ", path)
    if (nrow(groups) != 46L || length(unique(groups$group)) != expected_groups) {
        stop("Unexpected sample or group count in: ", path)
    }
    groups
}
align_groups <- function(groups, sample_ids) {
    if (anyDuplicated(sample_ids) || !setequal(groups$sample, sample_ids)) {
        stop("Grouping-table sample IDs do not match the genotype data.")
    }
    factor(groups$group[match(sample_ids, groups$sample)])
}
write_positions <- function(positions, filename) {
    write.table(positions, file.path(out_dir, filename), sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}
k5_groups <- read_groups(k5_path, 5L)
location_groups <- read_groups(locations_path, 10L)
if (!setequal(k5_groups$sample, location_groups$sample)) stop("Grouping tables differ in sample IDs.")

# 1. Convert the full SNP panel to PLINK format and run PCAdapt.
bed_prefix <- file.path(out_dir, "full_panel")
status <- system2("plink2", c("--vcf", shQuote(vcf_path), "--double-id",
                             "--make-bed", "--allow-extra-chr", "--out", shQuote(bed_prefix)))
if (status != 0L) stop("PLINK conversion failed; see the PLINK log.")
bim <- read.table(paste0(bed_prefix, ".bim"), header = FALSE,
                  colClasses = c("character", "character", "numeric", "numeric", "character", "character"))
fam <- read.table(paste0(bed_prefix, ".fam"), header = FALSE, colClasses = "character")
pca_groups <- align_groups(location_groups, fam[[2]])

pcadapt_input <- pcadapt::read.pcadapt(paste0(bed_prefix, ".bed"), type = "bed")
pcadapt_results <- pcadapt::pcadapt(pcadapt_input, K = 4, min.maf = 0.01)
save(pcadapt_results, file = file.path(out_dir, "PCAdapt.PCA_basic.results_full.panel.RData"))
if (length(pcadapt_results$pvalues) != nrow(bim)) stop("PCAdapt p-values and BIM loci differ in length.")
qval <- qvalue::qvalue(pcadapt_results$pvalues)$qvalues
pcadapt_positions <- data.frame(CHROM = bim[[1]], POS = bim[[4]], stringsAsFactors = FALSE)
pcadapt_positions <- pcadapt_positions[which(qval < 0.05), , drop = FALSE]
write_positions(pcadapt_positions, "outliers_p_0.05_from.PCAdapt_n46.tsv")

pdf(file.path(out_dir, "pcadapt_diagnostics.pdf"))
print(plot(pcadapt_results, option = "screeplot"))
print(plot(pcadapt_results, option = "scores", pop = pca_groups))
print(plot(pcadapt_results, option = "manhattan"))
print(plot(pcadapt_results, option = "qqplot"))
hist(pcadapt_results$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
print(plot(pcadapt_results, option = "stat.distribution"))
dev.off()

# 2. Prepare diploid allele counts for OutFLANK; exclude loci with any
# missing genotype. Keep the same locus mask for genotypes and coordinates.
vcf <- vcfR::read.vcfR(vcf_path)
gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)
vcf_samples <- colnames(gt)
k5_populations <- align_groups(k5_groups, vcf_samples)
location_populations <- align_groups(location_groups, vcf_samples)
coordinates <- data.frame(CHROM = vcf@fix[, "CHROM"],
                          POS = as.numeric(vcf@fix[, "POS"]), stringsAsFactors = FALSE)
if (anyNA(coordinates)) stop("Missing SNP coordinates.")
if (any(grepl(",", vcf@fix[, "ALT"], fixed = TRUE))) stop("Input must contain biallelic SNPs.")

dosage <- matrix(NA_integer_, nrow = nrow(gt), ncol = ncol(gt))
valid_gt <- c("0/0", "0/1", "1/0", "1/1", "0|0", "0|1", "1|0", "1|1")
for (j in seq_len(ncol(gt))) {
    g <- gt[, j]
    called <- !is.na(g) & !grepl(".", g, fixed = TRUE)
    if (any(called & !(g %in% valid_gt))) stop("Unexpected non-diploid/biallelic genotype in ", vcf_samples[j])
    dosage[g %in% c("0/0", "0|0"), j] <- 0L
    dosage[g %in% c("0/1", "1/0", "0|1", "1|0"), j] <- 1L
    dosage[g %in% c("1/1", "1|1"), j] <- 2L
}
valid_loci <- rowSums(is.na(dosage)) == 0L
if (!any(valid_loci)) stop("No complete loci remain for OutFLANK.")
SNPdata <- t(dosage[valid_loci, , drop = FALSE])
coordinates <- coordinates[valid_loci, , drop = FALSE]
rownames(SNPdata) <- vcf_samples
locus_names <- paste(coordinates$CHROM, coordinates$POS, sep = "_")
if (anyDuplicated(locus_names)) stop("Repeated CHROM/POS in the OutFLANK input.")
stopifnot(ncol(SNPdata) == nrow(coordinates),
          length(locus_names) == ncol(SNPdata), !anyNA(locus_names))
colnames(SNPdata) <- locus_names
rm(vcf, gt, dosage)
invisible(gc())

# 3. Run OutFLANK independently for K=5 and the ten sampling locations.
run_outflank <- function(populations, label, filename) {
    FstDataFrame <- OutFLANK::MakeDiploidFSTMat(SNPdata, locus_names, populations)
    save(FstDataFrame, file = file.path(out_dir, paste0("FstDataFrame_", label, ".RData")))
    fst_input <- FstDataFrame[!is.na(FstDataFrame$FST), , drop = FALSE]
    outliers <- OutFLANK::OutFLANK(fst_input,
                                  LeftTrimFraction = 0.05, RightTrimFraction = 0.05,
                                  Hmin = 0.10, NumberOfSamples = nlevels(populations),
                                  qthreshold = 0.05)
    save(outliers, file = file.path(out_dir, paste0("OutFLANK_", label, ".RData")))
    results <- outliers$results
    if (!any(!is.na(results$qvalues))) stop("No nonmissing OutFLANK q-values for ", label)
    cutoff <- quantile(results$qvalues, probs = 0.001, na.rm = TRUE)
    selected <- results[which(results$qvalues < cutoff), , drop = FALSE]
    indices <- match(as.character(selected$LocusName), locus_names)
    if (anyNA(indices)) stop("OutFLANK locus names cannot be matched to SNP coordinates.")
    positions <- coordinates[indices, , drop = FALSE]
    write_positions(positions, filename)
    write.csv(selected, file.path(out_dir, paste0("OutFLANK_", label, "_selected.csv")), row.names = FALSE)
    message(label, ": q-value cutoff = ", signif(cutoff, 8), "; selected loci = ", nrow(positions))
    positions
}
k5_positions <- run_outflank(k5_populations, "K5",
                             "outliers_new.panel_morePAN_btm.0.1.perc.qval_outlier.list_FOR.K5.tsv")
location_positions <- run_outflank(location_populations, "sampling_locations",
                                   "outliers_new.panel_morePAN_btm.0.1.perc.qval_outlier.list_samp.locs.tsv")

# 4. Merge all three lists, remove repeated positions, and sort.
merged <- unique(rbind(location_positions, k5_positions, pcadapt_positions))
merged <- merged[order(merged$CHROM, merged$POS), , drop = FALSE]
write_positions(merged, "merged.outlier.snps_sept.24.2025_n46.tsv")
message("Merged outlier positions: ", nrow(merged))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
