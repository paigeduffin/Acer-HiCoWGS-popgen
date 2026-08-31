# Variant filtering and neutral-panel construction

This directory contains scripts for filtering *Acropora cervicornis* SNPs
and constructing the neutral SNP panel for 46 samples. The workflow starts
with the merged SNP VCF produced during genotyping.

VCFs and intermediate analysis files are not stored in this repository.

## Workflow

| Order | File | Purpose |
|---|---|---|
| 1 | `01_ruth_hwe_filter.sh` | Estimate principal components, run RUTH, and filter SNPs for departures from Hardy–Weinberg equilibrium. |
| 2 | `02_quality_allele_filter.sh` | Apply depth, quality, and allele-count filters, then retain biallelic SNPs. |
| 3 | `03_sor_filter.sh` | Calculate strand odds ratio (SOR) and filter SNPs for strand bias. |
| 4 | `04_identify_neutral_panel_outliers.R` | Identify candidate outliers using PCAdapt and OutFLANK, then combine their SNP positions. |
| 5 | `05_build_neutral_panel.sh` | Remove candidate outliers, prune for linkage disequilibrium, and apply missingness and minor allele frequency filters. |

Run the scripts in numbered order. Steps 4 and 5 both use the full filtered
SNP panel from step 3; step 5 also uses the combined outlier list from step 4.
Input paths and command-line arguments are documented at the top of each script.

## Helper files

- `K5_morePAN_n46_2.col.txt`: sample assignments to the five inferred subpopulations.
- `samp.locs_FLsep_morePAN_n46_2.col.txt`: sample assignments to the ten sampling locations.

Both tables contain sample IDs and group assignments in two columns without a header.

## Software

The workflow uses bcftools 1.19, PLINK2 2.00a4.3, RUTH 3.29.4,
VCFtools 0.1.16, HTSlib 1.19.1, and R 4.4.1 with pcadapt 4.3.5,
OutFLANK 0.2, qvalue, and vcfR. Exact filtering options are recorded in the scripts.
