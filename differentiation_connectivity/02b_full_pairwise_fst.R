# Author: Paige Duffin
# Full-panel pairwise FST using StAMPP.
# Sampling locations: 37 samples, 3 bootstraps. K5: 46 samples, 5 bootstraps.
# The original 99% confidence-interval setting is retained. Only these small
# bootstrap counts completed; bootstrap CIs and p-values are poorly resolved.
# Usage: Rscript 02b_full_pairwise_fst.R locations|K5 [INPUT.vcf.gz] [NEW_OUTPUT_DIR] [NCORES]
# Run from this directory with the same two helper files used for neutral FST.

library(vcfR)
library(adegenet)
library(StAMPP)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 4L ||
    !args[1] %in% c("locations", "K5")) {
    stop("Usage: Rscript 02b_full_pairwise_fst.R locations|K5 [INPUT.vcf.gz] [NEW_OUTPUT_DIR] [NCORES]")
}
grouping <- args[1]
if (grouping == "locations") {
    default_vcf <- "snps.morePAN.minDP10.maxDP50.qual.SOR.miss10.nosort_biallel_subsamp.4.FST.IBD_from.n46.vcf.gz"
    pop_file <- "pop.IDs_samp.loc_n46_subsamp.4.evenness.txt"
    expected_n <- 37L
    nboots <- 3L
} else {
    default_vcf <- "snps.morePAN.minDP10.maxDP50.qual.SOR.miss10.nosort_biallel.vcf.gz"
    pop_file <- "pop.IDs_n46_posthoc.K5.txt"
    expected_n <- 46L
    nboots <- 5L
}
vcf_file <- if (length(args) >= 2L) args[2] else default_vcf
out_dir <- if (length(args) >= 3L) args[3] else
    file.path("full_pairwise_fst_results", grouping)
ncores <- as.integer(if (length(args) >= 4L) args[4] else
    Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (file.exists(out_dir)) stop("Output path already exists: ", out_dir)

vcf <- read.vcfR(vcf_file)
genlight_obj <- vcfR2genlight(vcf)
rm(vcf)

pop_table <- read.table(pop_file, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE)
colnames(pop_table) <- c("sample", "pop")
stopifnot(nrow(pop_table) == expected_n, !anyDuplicated(pop_table$sample))

# Retain the 37 listed individuals for locations, without filtering loci.
if (grouping == "locations") {
    genlight_obj <- genlight_obj[indNames(genlight_obj) %in% pop_table$sample, ]
}
stopifnot(setequal(pop_table$sample, indNames(genlight_obj)))
pop_table <- pop_table[match(indNames(genlight_obj), pop_table$sample), ]
stopifnot(all(pop_table$sample == indNames(genlight_obj)))
pop(genlight_obj) <- factor(pop_table$pop)
ploidy(genlight_obj) <- rep(2, nInd(genlight_obj))

if (!dir.create(out_dir, recursive = TRUE)) stop("Cannot create: ", out_dir)
res_stampp <- stamppFst(genlight_obj, nboots = nboots, percent = 99,
                        nclusters = ncores)

saveRDS(res_stampp, file.path(out_dir, "stampp_results.rds"))
write.csv(res_stampp$Fsts, file.path(out_dir, "pairwise_fst.csv"))
write.csv(res_stampp$Pvalues, file.path(out_dir, "pairwise_pvalues.csv"))
write.csv(res_stampp$Bootstraps, file.path(out_dir, "bootstraps.csv"))
