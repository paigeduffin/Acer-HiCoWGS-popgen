#!/usr/bin/env bash
#SBATCH --job-name=neutral_pairwise_fst
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --array=1-2
#SBATCH --output=neutral_pairwise_fst.%A_%a.out
#SBATCH --error=neutral_pairwise_fst.%A_%a.err

# Author: Paige Duffin
# Neutral-panel pairwise FST: task 1 = sampling locations (37 samples),
# task 2 = K5 groups (46 samples). R code is in 01b_neutral_pairwise_fst.R.
# Submit from this directory:
# sbatch 01a_neutral_pairwise_fst.sh [INPUT.vcf.gz] [NEW_OUTPUT_ROOT]
# Requires R packages vcfR, adegenet, and StAMPP.

set -euo pipefail
[[ $# -le 2 ]] || { echo "Usage: sbatch $0 [INPUT.vcf.gz] [NEW_OUTPUT_ROOT]" >&2; exit 1; }
cd -- "${SLURM_SUBMIT_DIR:-.}"

module load gcc/13.3.0 r/4.4.1 openblas/0.3.28 lapackpp/2023.11.05 openmpi/5.0.5 gdal/3.9.2 geos/3.12.2

case "${SLURM_ARRAY_TASK_ID:?Submit this script with sbatch}" in
    1) GROUP=locations ;;
    2) GROUP=K5 ;;
    *) echo "Expected array task 1 or 2." >&2; exit 1 ;;
esac

VCF=${1:-Acer_main.NEU.PANEL_n46.max10miss.MAF05.vcf.gz}
OUT_ROOT=${2:-neutral_pairwise_fst_results}
# Match StAMPP workers to the CPUs allocated to this R process.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript 01b_neutral_pairwise_fst.R "$GROUP" "$VCF" \
    "$OUT_ROOT/$GROUP" "${SLURM_CPUS_PER_TASK:-8}"
