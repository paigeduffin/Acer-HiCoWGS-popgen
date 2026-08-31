# Diversity and genetic uniqueness

Full-panel calculations of mean expected heterozygosity (He) and population-specific
FST for sampling locations (37 samples) and K5 groups (46 samples).

## Workflow

| Order | File | Purpose |
|---|---|---|
| 1a | `01a_he_pop_specific_fst.sh` | Submit the R calculations for both groupings. |
| 1b | `01b_he_pop_specific_fst.R` | Calculate and save mean He and population-specific FST. |

Submit the shell script from this directory; it calls the R script automatically.
Input paths and settings are documented in the scripts.

Plots and regressions comparing He with population-specific FST were produced
in Microsoft Excel.

## Helper files

Reused from `../differentiation_connectivity/`:

- `pop.IDs_samp.loc_n46_subsamp.4.evenness.txt`: sampling-location assignments.
- `pop.IDs_n46_posthoc.K5.txt`: K5 assignments.

## Software

R with vcfR, adegenet, and hierfstat; SLURM for job submission; Microsoft Excel
for plotting and regression.
