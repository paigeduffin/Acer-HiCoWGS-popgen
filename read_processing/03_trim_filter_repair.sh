#!/bin/bash
#SBATCH --job-name=trim_filter_repair
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --array=1-36
#SBATCH --output=trim_filter_repair.%A_%a.out
#SBATCH --error=trim_filter_repair.%A_%a.err

# Trim, quality-filter, and restore pairing for one Panama sample per
# array task. This consolidates the original fastp, decompression,
# custom-filter, and rePair wrapper scripts.
#
# Required programs: fastp, Python 3, Perl
# Usage:
#   sbatch 03_trim_filter_repair.sh [WORK_DIRECTORY] [MANIFEST]

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORK_DIR=${1:-$PWD}
MANIFEST=${2:-${SCRIPT_DIR}/panama_sample_manifest.tsv}
THREADS=${SLURM_CPUS_PER_TASK:-4}

RUN=$(awk -F '\t' 'NR > 1 && $3 == "yes" {print $1}' "$MANIFEST" \
    | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "$RUN" ]]; then
    echo "No reprocessed accession corresponds to array index ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

RAW_DIR="${WORK_DIR}/raw_reads"
TRIMMED_DIR="${WORK_DIR}/trimmed_reads"
FILTERED_DIR="${WORK_DIR}/quality_filtered_reads"
PAIRED_DIR="${WORK_DIR}/repaired_reads/paired"
UNPAIRED_DIR="${WORK_DIR}/repaired_reads/unpaired"
REPORT_DIR="${WORK_DIR}/fastp_reports"

mkdir -p "$TRIMMED_DIR" "$FILTERED_DIR" "$PAIRED_DIR" \
    "$UNPAIRED_DIR" "$REPORT_DIR"

R1_RAW="${RAW_DIR}/${RUN}_1.fastq.gz"
R2_RAW="${RAW_DIR}/${RUN}_2.fastq.gz"
R1_TRIMMED_GZ="${TRIMMED_DIR}/${RUN}_1_fastp.trimmed.fastq.gz"
R2_TRIMMED_GZ="${TRIMMED_DIR}/${RUN}_2_fastp.trimmed.fastq.gz"
R1_TRIMMED="${TRIMMED_DIR}/${RUN}_1_fastp.trimmed.fastq"
R2_TRIMMED="${TRIMMED_DIR}/${RUN}_2_fastp.trimmed.fastq"

fastp \
    --in1 "$R1_RAW" \
    --in2 "$R2_RAW" \
    --out1 "$R1_TRIMMED_GZ" \
    --out2 "$R2_TRIMMED_GZ" \
    --thread "$THREADS" \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required 140 \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --detect_adapter_for_pe \
    --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20 \
    --cut_front --cut_front_window_size 1 --cut_front_mean_quality 20 \
    --cut_right --cut_right_window_size 10 --cut_right_mean_quality 20 \
    --trim_poly_g --poly_g_min_len 10 \
    --trim_poly_x \
    --correction \
    --html "${REPORT_DIR}/${RUN}.fastp.html" \
    --json "${REPORT_DIR}/${RUN}.fastp.json"

# The custom filter requires uncompressed FASTQ input. Keep the gzipped
# fastp outputs as well as the temporary uncompressed copies.
gzip -dc "$R1_TRIMMED_GZ" > "$R1_TRIMMED"
gzip -dc "$R2_TRIMMED_GZ" > "$R2_TRIMMED"

python3 "${SCRIPT_DIR}/custom_qual.filt_90perc.over.20.py" \
    "$R1_TRIMMED" "${FILTERED_DIR}/${RUN}_R1.clean"
python3 "${SCRIPT_DIR}/custom_qual.filt_90perc.over.20.py" \
    "$R2_TRIMMED" "${FILTERED_DIR}/${RUN}_R2.clean"

# rePair.pl writes its outputs into the current directory and constructs
# their names from the two input filenames, so run it from FILTERED_DIR.
cd "$FILTERED_DIR"
perl "${SCRIPT_DIR}/rePair.pl" "${RUN}_R1.clean" "${RUN}_R2.clean"

mv "R1_${RUN}_R1.clean" "$PAIRED_DIR/"
mv "R2_${RUN}_R2.clean" "$PAIRED_DIR/"
mv "Unp_${RUN}_R1.clean_${RUN}_R2.clean" "$UNPAIRED_DIR/"

