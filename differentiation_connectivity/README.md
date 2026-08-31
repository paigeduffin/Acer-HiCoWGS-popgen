# Differentiation and connectivity

Pairwise FST analyses for the neutral and full SNP panels, plus isolation by
distance (IBD) using neutral-panel FST. Sampling-location analyses use 37 samples;
K5 analyses use 46 samples.

## Workflow

| Order | File | Purpose |
|---|---|---|
| 1a | `01a_neutral_pairwise_fst.sh` | Submit neutral-panel FST jobs for both groupings. |
| 1b | `01b_neutral_pairwise_fst.R` | Calculate neutral-panel pairwise FST and bootstrap results. |
| 2a | `02a_full_pairwise_fst.sh` | Submit full-panel FST jobs for both groupings. |
| 2b | `02b_full_pairwise_fst.R` | Calculate full-panel pairwise FST and bootstrap results. |
| 3 | `03_ibd.R` | Produce IBD plots and reproduce the original Mantel calculation. |

Submit shell scripts from this directory; each calls its paired R script.
IBD uses the separately formatted FST table listed below. Input paths and
analysis settings are documented in the scripts.

## Helper files

- `pop.IDs_samp.loc_n46_subsamp.4.evenness.txt`: sampling-location assignments.
- `pop.IDs_n46_posthoc.K5.txt`: K5 assignments.
- `Acer.haplo.locs_ABREV_w.panama.csv`: location coordinates for IBD.
- `STAMPP_fst.bootstrap_samploc.NEU.pan_178k_4.IBD.csv`: neutral-panel FST table for IBD, with `popA`, `popB`, and `fst` columns.

## Software

R with vcfR, adegenet, StAMPP, geosphere, readr, dplyr, ggplot2, reshape2,
and ade4; SLURM for job submission.
