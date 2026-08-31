# Population structure

This directory contains PCA, sNMF, DAPC, ADMIXTURE, and PopCluster analyses
of the neutral SNP panel for 46 *Acropora cervicornis* samples.

Genotype data and intermediate analysis files are not stored in this repository.

## Workflow

| Order | File | Purpose |
|---|---|---|
| 1a | `01a_pca_prepare.sh` | Convert the neutral VCF to PLINK files for PCA and PopCluster. |
| 1b | `01b_pca.R` | Run PCAdapt PCA and save scores, singular values, and plots. |
| 2a | `02a_snmf_prepare.sh` | Rename chromosomes and create PED files for sNMF. |
| 2b | `02b_snmf.R` | Run sNMF for K = 1–10 with 20 replicates and export Q matrices for CLUMPAK. |
| 3 | `03_dapc.R` | Summarize BIC across 50 runs and plot DAPC using the final run's automatically selected clusters. |
| 4 | `04_admixture.sh` | Retain the 15 main chromosomes, run ADMIXTURE, format Q matrices for CLUMPAK, and extract log summaries. |
| 5a | `05a_popcluster.sh` | Run PopCluster for K = 1–10 in five replicate jobs. |
| 5b | `05b_popcluster_collect.sh` | Format PopCluster Q matrices for CLUMPAK and collect K-file summaries. |

Run paired a/b scripts in that order. PopCluster reuses the PLINK files from
step 1a; run step 5b after all five jobs from step 5a finish. Default paths
assume execution from this directory. Script headers describe input arguments.

## Helper files

- `main.15.chromo_rename.txt`: chromosome-name mapping used by sNMF and ADMIXTURE.
- `Acer_popgen_n46_178k.PcPjt`: PopCluster project settings.
- `../variant_filtering/samp.locs_FLsep_morePAN_n46_2.col.txt`: shared sample-location table used by PCA and DAPC; retained in `variant_filtering/`.

## Software

Shell steps use bcftools, HTSlib, VCFtools, PLINK2, ADMIXTURE, and PopCluster.
R steps use pcadapt, ggplot2, LEA, vcfR, adegenet, and ade4.
Exact analysis settings are recorded in the scripts.
