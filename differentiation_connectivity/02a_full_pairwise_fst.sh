#!/usr/bin/env bash
#SBATCH --job-name=full_pairwise_fst
#SBATCH --partition=largemem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --array=1-2
#SBATCH --output=full_pairwise_fst.%A_%a.out
#SBATCH --error=full_pairwise_fst.%A_%a.err

# Author: Paige Duffin
# Full-panel pairwise FST: task 1 = sampling locations (37 samples, 3 bootstraps),
# task 2 = K5 groups (46 samples, 5 bootstraps).
# Submit from this directory:
# sbatch 02a_full_pairwise_fst.sh [INPUT.vcf.gz] [NEW_OUTPUT_ROOT]
# Without an input argument, each task uses its original VCF filename below.
# Requires R packages vcfR, adegenet, and StAMPP.

set -euo pipefail
cd -- "${SLURM_SUBMIT_DIR:-.}"

module load gcc/13.3.0 r/4.4.1 openblas/0.3.28 lapackpp/2023.11.05 openmpi/5.0.5 gdal/3.9.2 geos/3.12.2

case "${SLURM_ARRAY_TASK_ID:?Submit this script with sbatch}" in
    1)
        GROUP=locations
        DEFAULT_VCF=snps.morePAN.minDP10.maxDP50.qual.SOR.miss10.nosort_biallel_subsamp.4.FST.IBD_from.n46.vcf.gz
        ;;
    2)
        GROUP=K5
        DEFAULT_VCF=snps.morePAN.minDP10.maxDP50.qual.SOR.miss10.nosort_biallel.vcf.gz
        ;;
    *) echo "Expected array task 1 or 2." >&2; exit 1 ;;
esac

VCF=${1:-$DEFAULT_VCF}
OUT_ROOT=${2:-full_pairwise_fst_results}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript 02b_full_pairwise_fst.R "$GROUP" "$VCF" \
    "$OUT_ROOT/$GROUP" "${SLURM_CPUS_PER_TASK:-8}"
