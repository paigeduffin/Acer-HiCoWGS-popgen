# Haplotype panel construction and imputation

This directory contains scripts used to construct the *Acropora cervicornis*
haplotype panel, impute the 330-sample low-coverage dataset, and assess accuracy
using five samples sequenced at both high and low coverage. The recorded
validation used the complete 44-sample reference panel.

Large input files, reference genomes, BAMs, VCFs, and software executables are
not stored in this repository.

## Workflow

| Order | File | Purpose |
|---|---|---|
| 1 | `01_prepare_whatshap_input.sh` | Prepare the 46-sample high-coverage VCF for read-backed phasing. |
| 2 | `02_whatshap_phase.sh` | Phase each sample with WhatsHap and merge the phased VCFs. |
| 3 | `03_between_phasing_pca_admixture.sh` | Examine population structure and remove two admixed samples from the reference panel. |
| 4 | `04_ruth_hwe_filter.sh` | Apply RUTH/HWE and repeat-mask filtering to the 44-sample panel. |
| 5 | `05_beagle_population_phasing.sh` | Perform population-based phasing with Beagle. |
| 6 | `06_prepare_imputation_reference.sh` | Prepare chromosome-specific Beagle reference files. |
| 7 | `07_two_round_loco_imputation.sh` | Recalibrate genotype probabilities and impute the 330-sample low-coverage VCF. |
| 8a | `08a_extract_validation_genotypes.sh` | Extract raw, recalibrated, and imputed genotypes for the five validation pairs. |
| 8b | `08b_calculate_validation_accuracy.R` | Calculate genotype concordance and Pearson correlation across validation stages. |

## Helper files

- `Acervicornis.repeatmask.bed` (not included): repeat-mask intervals used during
  RUTH/HWE filtering. Source: [Locatelli N. 2024. *Genomic insights into coral
  evolution and adaptation: A comparative study of Caribbean reef-builders.*
  Pennsylvania State University.](http://paperpile.com/b/BDpMYz/RKcV)
- `gen.map_sep.scaffs/`: chromosome-specific genetic maps used for Beagle phasing.
- `Acervicornis_SEXAVG.map`: genetic map used for low-coverage imputation.
- `../population_structure/main.15.chromo_rename.txt`: chromosome-name conversion file used in step 03.

## Software

WhatsHap, bcftools, HTSlib, VCFtools, PLINK2, ADMIXTURE, RUTH, Java, Beagle
4.1/5.5, and R with dplyr, tidyr, readr, and ggplot2. Exact options are
recorded in the scripts.
