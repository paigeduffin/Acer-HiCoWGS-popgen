# Author: Paige Duffin
# Neutral-panel pairwise FST using StAMPP: 1,000 locus bootstraps, 95% CIs.
# Sampling locations use the 37 listed samples; K5 uses all 46 samples.
# The existing neutral-panel loci are retained without additional filtering.
# Usage: Rscript 01b_neutral_pairwise_fst.R locations|K5 [INPUT.vcf.gz] [NEW_OUTPUT_DIR] [NCORES]
# Run from this directory with both population-assignment helper files present.

library(vcfR)
library(adegenet)
library(StAMPP)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 4L ||
    !args[1] %in% c("locations", "K5")) {
    stop("Usage: Rscript 01b_neutral_pairwise_fst.R locations|K5 [INPUT.vcf.gz] [NEW_OUTPUT_DIR] [NCORES]")
}
grouping <- args[1]
vcf_file <- if (length(args) >= 2L) args[2] else
    "Acer_main.NEU.PANEL_n46.max10miss.MAF05.vcf.gz"
out_dir <- if (length(args) >= 3L) args[3] else
    file.path("neutral_pairwise_fst_results", grouping)
ncores <- as.integer(if (length(args) >= 4L) args[4] else
    Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (file.exists(out_dir)) stop("Output path already exists: ", out_dir)

pop_file <- if (grouping == "locations")
    "pop.IDs_samp.loc_n46_subsamp.4.evenness.txt" else
    "pop.IDs_n46_posthoc.K5.txt"
pop_table <- read.table(pop_file, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE)
colnames(pop_table) <- c("sample", "pop")
expected_n <- if (grouping == "locations") 37L else 46L
stopifnot(nrow(pop_table) == expected_n, !anyDuplicated(pop_table$sample))

vcf <- read.vcfR(vcf_file)
genlight_obj <- vcfR2genlight(vcf)
rm(vcf)
# Keep only the 37 listed individuals for the sampling-location analysis.
if (grouping == "locations") {
    genlight_obj <- genlight_obj[indNames(genlight_obj) %in% pop_table$sample, ]
}
stopifnot(setequal(pop_table$sample, indNames(genlight_obj)))
pop_table <- pop_table[match(indNames(genlight_obj), pop_table$sample), ]
stopifnot(all(pop_table$sample == indNames(genlight_obj)))
pop(genlight_obj) <- factor(pop_table$pop)
ploidy(genlight_obj) <- rep(2, nInd(genlight_obj))

if (!dir.create(out_dir, recursive = TRUE)) stop("Cannot create: ", out_dir)
res_stampp <- stamppFst(genlight_obj, nboots = 1000, percent = 95,
                        nclusters = ncores)

# Preserve the complete result and export the original three result components.
saveRDS(res_stampp, file.path(out_dir, "stampp_results.rds"))
write.csv(res_stampp$Fsts, file.path(out_dir, "pairwise_fst.csv"))
write.csv(res_stampp$Pvalues, file.path(out_dir, "pairwise_pvalues.csv"))
write.csv(res_stampp$Bootstraps, file.path(out_dir, "bootstraps.csv"))
