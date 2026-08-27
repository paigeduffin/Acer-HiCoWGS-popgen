#!/bin/bash
#SBATCH --account=CEElab
#SBATCH --job-name=download_panama_sra
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=1-36
#SBATCH --output=download_panama_sra.%A_%a.out
#SBATCH --error=download_panama_sra.%A_%a.err

# Download the 36 Panama SRA accessions reprocessed for this study.
# The array index selects one accession marked "yes" in the
# reprocessed_2025 column of panama_sample_manifest.tsv.
#
# Usage:
#   sbatch 01_download_sra.sh [MANIFEST] [OUTPUT_DIRECTORY]

set -euo pipefail

MANIFEST=${1:-panama_sample_manifest.tsv}
OUTPUT_DIR=${2:-sra_downloads}

ml gcc/13.3.0 sratoolkit/3.2.0

RUN=$(awk -F '\t' 'NR > 1 && $3 == "yes" {print $1}' "$MANIFEST" \
    | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "$RUN" ]]; then
    echo "No reprocessed accession corresponds to array index ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "Downloading ${RUN}"
prefetch "$RUN"

