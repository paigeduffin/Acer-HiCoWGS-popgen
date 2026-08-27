#!/bin/bash
#SBATCH --job-name=finalize_host_bams
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --array=1-36
#SBATCH --output=finalize_host_bams.%A_%a.out
#SBATCH --error=finalize_host_bams.%A_%a.err

# Merge paired and unpaired host BAMs, add read groups, fix mate
# information, mark duplicates, and produce the sorted/indexed BAM used
# for genotyping. This combines the original array wrapper and per-sample
# Picard script.
#
# Required programs: samtools, Java, Picard 2.26.2
# Before submission, define the Picard JAR location, for example:
#   export PICARD_JAR=/path/to/picard.jar
#
# Usage:
#   sbatch 06_finalize_host_bams.sh [WORK_DIRECTORY] [MANIFEST]

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORK_DIR=${1:-$PWD}
MANIFEST=${2:-${SCRIPT_DIR}/panama_sample_manifest.tsv}
THREADS=${SLURM_CPUS_PER_TASK:-32}
: "${PICARD_JAR:?Set PICARD_JAR to the Picard 2.26.2 JAR location before submission}"

RUN=$(awk -F '\t' 'NR > 1 && $3 == "yes" {print $1}' "$MANIFEST" \
    | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "$RUN" ]]; then
    echo "No reprocessed accession corresponds to array index ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

HOST_DIR="${WORK_DIR}/host_bams"
PAIRED_BAM="${HOST_DIR}/paired_host/${RUN}.paired_sorted_host.bam"
UNPAIRED_BAM="${HOST_DIR}/unpaired_host/${RUN}.unpaired_sorted_host.bam"
OUTPUT_DIR="${HOST_DIR}/merged"
mkdir -p "$OUTPUT_DIR"

MERGED_BAM="${OUTPUT_DIR}/${RUN}_host_merged.bam"
MERGED_SORTED_BAM="${OUTPUT_DIR}/${RUN}_host_merged.sorted.bam"
RG_BAM="${OUTPUT_DIR}/${RUN}_host_RG.bam"
RG_SORTED_BAM="${OUTPUT_DIR}/${RUN}_host_RG.sort.bam"
FIXED_BAM="${OUTPUT_DIR}/${RUN}_host_RG.sort.fxm8s.bam"
DEDUP_BAM="${OUTPUT_DIR}/${RUN}_host_GATKd.bam"
FINAL_BAM="${OUTPUT_DIR}/${RUN}_GATKd_sort.bam"

samtools merge "$MERGED_BAM" "$UNPAIRED_BAM" "$PAIRED_BAM"
samtools sort -@ "$THREADS" -O BAM -o "$MERGED_SORTED_BAM" "$MERGED_BAM"

samtools depth "$MERGED_SORTED_BAM" \
    | awk '{sum += $3; sites++} END {if (sites > 0) print sum/sites; else print "NA"}' \
    > "${OUTPUT_DIR}/${RUN}_host_average_depth.txt"

java -Xmx4G -jar "$PICARD_JAR" AddOrReplaceReadGroups \
    I="$MERGED_SORTED_BAM" \
    O="$RG_BAM" \
    RGID="$RUN" \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM="$RUN"

samtools sort -@ "$THREADS" -O BAM -o "$RG_SORTED_BAM" "$RG_BAM"

java -Xmx4G -jar "$PICARD_JAR" FixMateInformation \
    I="$RG_SORTED_BAM" \
    O="$FIXED_BAM"

java -Xmx4G -jar "$PICARD_JAR" MarkDuplicates \
    I="$FIXED_BAM" \
    O="$DEDUP_BAM" \
    M="${OUTPUT_DIR}/${RUN}_host_metrics.txt"

# Maria's genotyping workflow expects the filename RUN_GATKd_sort.bam.
samtools sort -@ "$THREADS" -O BAM -o "$FINAL_BAM" "$DEDUP_BAM"
samtools index "$FINAL_BAM"

# Remove intermediate BAMs only after the final BAM and index succeed.
rm "$MERGED_BAM" "$RG_BAM" "$RG_SORTED_BAM" "$FIXED_BAM"
