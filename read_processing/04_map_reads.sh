#!/bin/bash
#SBATCH --job-name=map_panama_reads
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --array=1-36
#SBATCH --output=map_panama_reads.%A_%a.out
#SBATCH --error=map_panama_reads.%A_%a.err

# Map repaired paired and unpaired reads independently to the combined
# A. cervicornis-S. fitti reference.
#
# Required program: BWA 0.7.17
# Usage:
#   sbatch 04_map_reads.sh [WORK_DIRECTORY] [REFERENCE_FASTA] [MANIFEST]

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORK_DIR=${1:-$PWD}
REFERENCE=${2:-${WORK_DIR}/GCA_964034985.1_Acer_S.fitti.fa}
MANIFEST=${3:-${SCRIPT_DIR}/panama_sample_manifest.tsv}
THREADS=${SLURM_CPUS_PER_TASK:-32}

RUN=$(awk -F '\t' 'NR > 1 && $3 == "yes" {print $1}' "$MANIFEST" \
    | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "$RUN" ]]; then
    echo "No reprocessed accession corresponds to array index ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

PAIRED_INPUT="${WORK_DIR}/repaired_reads/paired"
UNPAIRED_INPUT="${WORK_DIR}/repaired_reads/unpaired"
PAIRED_OUTPUT="${WORK_DIR}/mapped_reads/paired"
UNPAIRED_OUTPUT="${WORK_DIR}/mapped_reads/unpaired"

mkdir -p "$PAIRED_OUTPUT" "$UNPAIRED_OUTPUT"

bwa mem -t "$THREADS" "$REFERENCE" \
    "${PAIRED_INPUT}/R1_${RUN}_R1.clean" \
    "${PAIRED_INPUT}/R2_${RUN}_R2.clean" \
    > "${PAIRED_OUTPUT}/${RUN}.paired.bwa.sam"

bwa mem -t "$THREADS" "$REFERENCE" \
    "${UNPAIRED_INPUT}/Unp_${RUN}_R1.clean_${RUN}_R2.clean" \
    > "${UNPAIRED_OUTPUT}/${RUN}.unpaired.bwa.sam"

