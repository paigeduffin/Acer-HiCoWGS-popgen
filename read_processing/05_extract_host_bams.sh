#!/bin/bash
#SBATCH --job-name=extract_host_bams
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --array=1-36
#SBATCH --output=extract_host_bams.%A_%a.out
#SBATCH --error=extract_host_bams.%A_%a.err

# Remove alignments containing S. fitti reference names (prefix k127_)
# and create sorted host BAMs from the paired and unpaired SAM files.
# Symbiont BAM processing was conducted separately and is not included.
#
# Required program: samtools
# Usage:
#   sbatch 05_extract_host_bams.sh [WORK_DIRECTORY] [MANIFEST]

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORK_DIR=${1:-$PWD}
MANIFEST=${2:-${SCRIPT_DIR}/panama_sample_manifest.tsv}
THREADS=${SLURM_CPUS_PER_TASK:-32}

RUN=$(awk -F '\t' 'NR > 1 && $3 == "yes" {print $1}' "$MANIFEST" \
    | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "$RUN" ]]; then
    echo "No reprocessed accession corresponds to array index ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

PAIRED_SAM="${WORK_DIR}/mapped_reads/paired/${RUN}.paired.bwa.sam"
UNPAIRED_SAM="${WORK_DIR}/mapped_reads/unpaired/${RUN}.unpaired.bwa.sam"
PAIRED_HOST_DIR="${WORK_DIR}/host_bams/paired_host"
UNPAIRED_HOST_DIR="${WORK_DIR}/host_bams/unpaired_host"

mkdir -p "$PAIRED_HOST_DIR" "$UNPAIRED_HOST_DIR"

grep -v 'k127_' "$PAIRED_SAM" \
    | samtools view -b -S - \
    | samtools sort -@ "$THREADS" -O BAM \
        -o "${PAIRED_HOST_DIR}/${RUN}.paired_sorted_host.bam"

grep -v 'k127_' "$UNPAIRED_SAM" \
    | samtools view -b -S - \
    | samtools sort -@ "$THREADS" -O BAM \
        -o "${UNPAIRED_HOST_DIR}/${RUN}.unpaired_sorted_host.bam"
