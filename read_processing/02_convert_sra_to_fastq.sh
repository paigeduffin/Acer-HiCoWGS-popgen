#!/bin/bash
#SBATCH --job-name=convert_panama_sra
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --array=1-36
#SBATCH --output=convert_panama_sra.%A_%a.out
#SBATCH --error=convert_panama_sra.%A_%a.err

# Convert the 36 downloaded SRA files to paired, gzipped FASTQ files.
#
# Usage:
#   sbatch 02_convert_sra_to_fastq.sh [MANIFEST] [SRA_DIRECTORY] [FASTQ_DIRECTORY]

set -euo pipefail

MANIFEST=${1:-panama_sample_manifest.tsv}
SRA_DIR=${2:-sra_downloads}
FASTQ_DIR=${3:-raw_reads}

ml gcc/11.3.0 fastq-dump/2.11.0

RUN=$(awk -F '\t' 'NR > 1 && $3 == "yes" {print $1}' "$MANIFEST" \
    | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "$RUN" ]]; then
    echo "No reprocessed accession corresponds to array index ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

SRA_FILE="${SRA_DIR}/${RUN}/${RUN}.sra"
if [[ ! -f "$SRA_FILE" ]]; then
    echo "SRA file not found: ${SRA_FILE}" >&2
    exit 1
fi

mkdir -p "$FASTQ_DIR"

echo "Converting ${RUN}"
fastq-dump --split-files --gzip --outdir "$FASTQ_DIR" "$SRA_FILE"
