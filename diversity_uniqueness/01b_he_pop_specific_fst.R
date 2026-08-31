# Author: Paige Duffin
# Full-panel mean expected heterozygosity and population-specific FST.
# Sampling locations: 37 samples. K5: 46 samples.
# Usage: Rscript 01b_he_pop_specific_fst.R locations|K5 [INPUT.vcf.gz] [NEW_OUTPUT_DIR]
# Run from diversity_uniqueness/. Population tables are reused from
# ../differentiation_connectivity/; no duplicate helper files are needed.

library(vcfR)
library(adegenet)
library(hierfstat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 3L ||
    !args[1] %in% c("locations", "K5")) {
    stop("Usage: Rscript 01b_he_pop_specific_fst.R locations|K5 [INPUT.vcf.gz] [NEW_OUTPUT_DIR]")
}
grouping <- args[1]
if (grouping == "locations") {
    default_vcf <- "snps.morePAN.minDP10.maxDP50.qual.SOR.miss10.nosort_biallel_subsamp.4.FST.IBD_from.n46.vcf.gz"
    pop_file <- "../differentiation_connectivity/pop.IDs_samp.loc_n46_subsamp.4.evenness.txt"
    expected_n <- 37L
} else {
    default_vcf <- "snps.morePAN.minDP10.maxDP50.qual.SOR.miss10.nosort_biallel.vcf.gz"
    pop_file <- "../differentiation_connectivity/pop.IDs_n46_posthoc.K5.txt"
    expected_n <- 46L
}
vcf_file <- if (length(args) >= 2L) args[2] else default_vcf
out_dir <- if (length(args) >= 3L) args[3] else
    file.path("he_pop_specific_fst_results", grouping)
if (file.exists(out_dir)) stop("Output path already exists: ", out_dir)

# Shared preparation for the two calculations.
vcf <- read.vcfR(vcf_file)
genlte_obj <- vcfR2genlight(vcf)
rm(vcf)

pop_table <- read.table(pop_file, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE)
colnames(pop_table) <- c("sample", "pop")
stopifnot(nrow(pop_table) == expected_n, !anyDuplicated(pop_table$sample))
if (grouping == "locations") {
    genlte_obj <- genlte_obj[indNames(genlte_obj) %in% pop_table$sample, ]
}
stopifnot(setequal(pop_table$sample, indNames(genlte_obj)))
pop_table <- pop_table[match(indNames(genlte_obj), pop_table$sample), ]
stopifnot(all(pop_table$sample == indNames(genlte_obj)))
pop(genlte_obj) <- factor(pop_table$pop)
ploidy(genlte_obj) <- rep(2, nInd(genlte_obj))
pops <- levels(pop(genlte_obj))

if (!dir.create(out_dir, recursive = TRUE)) stop("Cannot create: ", out_dir)

# Mean He across loci, retaining the original glMean and missing-value handling.
He <- sapply(pops, function(p) {
    x <- genlte_obj[pop(genlte_obj) == p, ]
    allele_freq <- glMean(x)
    mean(2 * allele_freq * (1 - allele_freq), na.rm = TRUE)
})
print(He)
write.csv(data.frame(population = names(He), mean_He = as.numeric(He)),
          file.path(out_dir, "mean_he.csv"), row.names = FALSE)
saveRDS(He, file.path(out_dir, "mean_he.rds"))

# Population-specific FST from the original allele-dosage calculation.
dos <- as.matrix(genlte_obj)
psfst <- fs.dosage(dos = dos, pop = pop(genlte_obj))
print(psfst$Fs)

# Preserve the complete Fs table: its second row contains FST estimates;
# its first row contains FIS. Both include an overall summary.
write.csv(psfst$Fs, file.path(out_dir, "fs_statistics.csv"))
saveRDS(psfst, file.path(out_dir, "population_specific_fst.rds"))
