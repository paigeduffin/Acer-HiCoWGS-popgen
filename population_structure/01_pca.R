#!/usr/bin/env Rscript
# PCA of the neutral SNP panel using PCAdapt.
# Software: PLINK2 2.00a4.3, R 4.4.1, pcadapt 4.3.5, and ggplot2.
# Load the software before running; plink2 must be on PATH.
# Run from population_structure/:
# Rscript 01_pca.R [VCF] [LOCATION_TABLE] [NEW_OUTPUT_DIRECTORY]
# The location table has two columns without a header: sample ID and location.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L) {
    stop("Usage: Rscript 01_pca.R [VCF] [LOCATION_TABLE] [NEW_OUTPUT_DIRECTORY]")
}
defaults <- c("Acer_main.NEU.PANEL_n46.max10miss.MAF05.vcf.gz",
              "../variant_filtering/samp.locs_FLsep_morePAN_n46_2.col.txt",
              "pca_results")
defaults[seq_along(args)] <- args
vcf_path <- defaults[1]
location_path <- defaults[2]
out_dir <- defaults[3]

for (pkg in c("pcadapt", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
}
if (Sys.which("plink2") == "") stop("Load PLINK2 before running this script.")
if (!all(file.exists(defaults[1:2]))) stop("The VCF or location table is missing.")
if (file.exists(out_dir)) stop("Use a new output directory: ", out_dir)

locations <- read.table(location_path, header = FALSE, sep = "\t",
                        colClasses = "character", comment.char = "", quote = "")
if (ncol(locations) != 2L || anyNA(locations) || any(locations == "")) {
    stop("The location table must have two nonempty columns.")
}
if (anyDuplicated(locations[[1]])) stop("Duplicate sample IDs in the location table.")
if (!dir.create(out_dir, recursive = TRUE)) stop("Cannot create output directory.")

# 1. Convert the final neutral VCF to PLINK format.
bed_prefix <- file.path(out_dir, "Acer_main.NEU.PANEL_n46.max10miss.MAF05")
status <- system2("plink2", c("--vcf", shQuote(vcf_path), "--make-bed",
                             "--allow-extra-chr", "--double-id",
                             "--out", shQuote(bed_prefix)))
if (status != 0L) stop("PLINK conversion failed; see the PLINK log.")

# Match locations to the PLINK sample order by ID, not by table row order.
fam <- read.table(paste0(bed_prefix, ".fam"), header = FALSE, colClasses = "character")
sample_ids <- fam[[2]]
if (anyDuplicated(sample_ids) || !setequal(sample_ids, locations[[1]])) {
    stop("Sample IDs in the location table and PLINK data do not match.")
}
pop_table <- locations[[2]][match(sample_ids, locations[[1]])]

# 2. Main PCA and basic plots: K = 10, min.maf = 0.01.
pcadapt_file <- pcadapt::read.pcadapt(paste0(bed_prefix, ".bed"), type = "bed")
pcadapt_results <- pcadapt::pcadapt(input = pcadapt_file, K = 10, min.maf = 0.01)
stopifnot(nrow(pcadapt_results$scores) == length(sample_ids))

scores <- as.data.frame(pcadapt_results$scores)
names(scores) <- paste0("PC", seq_len(ncol(scores)))
write.table(data.frame(sample = sample_ids, location = pop_table, scores),
            file.path(out_dir, "pca_scores.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(data.frame(PC = seq_along(pcadapt_results$singular.values),
                       singular_value = pcadapt_results$singular.values),
            file.path(out_dir, "pca_singular_values.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

pdf(file.path(out_dir, "PCAdapt.PCA_plots.A.B_neu.pan_MAF05.pdf"))
print(plot(pcadapt_results, option = "screeplot"))
print(plot(pcadapt_results, option = "scores", pop = pop_table))
dev.off()

# 3. K = 2 grabbing coords for manual plotting
pcadapt_results_PC1.PC2 <- pcadapt::pcadapt(input = pcadapt_file, K = 2)
writeLines(capture.output(summary(pcadapt_results_PC1.PC2)),
           file.path(out_dir, "pcadapt_K2_summary.txt"))
save(pcadapt_results, pcadapt_results_PC1.PC2, sample_ids, pop_table,
     file = file.path(out_dir, "PCAdapt.PCA_neutral_panel.RData"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("PCA results saved in: ", out_dir)
