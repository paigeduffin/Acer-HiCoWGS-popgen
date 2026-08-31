#!/usr/bin/env Rscript
# Run sNMF on the neutral SNP panel and export every Q matrix for CLUMPAK.
# Software: R 4.4.1 and LEA.
# Prepare the PED file first with 02a_snmf_prepare.sh.
# Run from population_structure/:
# Rscript 02b_snmf.R [INPUT.ped] [NEW_OUTPUT_DIRECTORY]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
    stop("Usage: Rscript 02b_snmf.R [INPUT.ped] [NEW_OUTPUT_DIRECTORY]")
}
defaults <- c("snmf_input/Acer_main.NEU.PANEL_n46_plink_miss.MAF.filt.ped",
              "snmf_results")
defaults[seq_along(args)] <- args
ped_input_file <- defaults[1]
out_dir <- defaults[2]

if (!requireNamespace("LEA", quietly = TRUE)) stop("Install/load the LEA package.")
if (!file.exists(ped_input_file)) stop("Missing PED input: ", ped_input_file)
if (file.exists(out_dir)) stop("Use a new output directory: ", out_dir)
if (!dir.create(out_dir, recursive = TRUE)) stop("Cannot create output directory.")
ped_input_file <- normalizePath(ped_input_file, mustWork = TRUE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

# 1. Convert PED to GENO once, keeping the generated files in the results directory.
geno_file <- LEA::ped2geno(ped_input_file,
                          output.file = file.path(out_dir, "neutral_panel.geno"))

# 2. Estimate ancestry for K = 1-10, with 20 runs per K and cross-entropy.
obj.snmf_k.1.to.10_rep.20 <- LEA::snmf(geno_file, K = 1:10, entropy = TRUE,
                                     ploidy = 2, project = "new", repetitions = 20)
save(obj.snmf_k.1.to.10_rep.20,
     file = file.path(out_dir, "sNMF_K1_to_10_20_repetitions.RData"))

# 3. Export all 20 Q matrices for each K, grouped into CLUMPAK directories.
output_root <- file.path(out_dir, "clumpak")
for (K_val in 1:10) {
    dir_name <- file.path(output_root,
                          paste0("K", K_val, "_Acer.neu.pan_n46_snmf_admixture"))
    if (!dir.create(dir_name, recursive = TRUE)) stop("Cannot create: ", dir_name)

    for (run_num in 1:20) {
        Q_matrix <- LEA::Q(obj.snmf_k.1.to.10_rep.20, K = K_val, run = run_num)
        file_name <- paste0("Acer.neu.pan_n46_snmf_K", K_val,
                             "_run", run_num, "_admixture.txt")
        write.table(format(Q_matrix, scientific = FALSE, trim = TRUE),
                    file = file.path(dir_name, file_name), sep = "\t",
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("sNMF results and CLUMPAK files saved in: ", out_dir)

# After downloading the results directory, open a terminal inside that directory.
# Package the K-specific folders for CLUMPAK with:
# zip -r snmf_for_CLUMPAK.zip clumpak/K*_Acer.neu.pan_n46_snmf_admixture -x "*/.*"
