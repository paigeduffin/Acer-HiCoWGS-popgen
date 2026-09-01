#!/usr/bin/env bash
#SBATCH --job-name=prepare_loco_bams
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --array=1-330
#SBATCH --output=prepare_loco_bams.%A_%a.out
#SBATCH --error=prepare_loco_bams.%A_%a.err

# Author: Paige Duffin
# Process the 330 paired-end LoCo libraries from raw reads through the
# host/symbiont parsing and Picard steps used to create *_host_GATKd.bam.
# Only repaired paired reads were mapped in the original workflow; unpaired
# reads produced by rePair.pl are retained but are not mapped.
#
# Required programs: Trim Galore 0.6.6, Cutadapt 3.5, FastQC 0.11.9,
# Trimmomatic 0.39, seqtk 1.3, FASTX-Toolkit fastq_quality_filter,
# Perl, BWA 0.7.17, samtools, Java 17, and Picard 2.26.2.
# Set TRIMMOMATIC_JAR and PICARD_JAR before submission. REPAIR_PL and
# NEXTERA_ADAPTERS may also be set when their files are elsewhere.
#
# Usage:
#   sbatch 07a_prepare_loco_bams.sh [RAW_READ_DIRECTORY] [COMBINED_REFERENCE.fa] [NEW_OUTPUT_DIRECTORY] [MANIFEST]

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
RAW_DIR=${1:-raw_loco_reads}
REFERENCE=${2:-ref_genome_indexed/GCA_964034985.1_Acer_S.fitti.fa}
OUT_DIR=${3:-loco_read_processing}
MANIFEST=${4:-${SCRIPT_DIR}/loco_sample_manifest.tsv}
THREADS=${SLURM_CPUS_PER_TASK:-32}
TASK_ID=${SLURM_ARRAY_TASK_ID:?Submit as a 1-330 Slurm array}

: "${TRIMMOMATIC_JAR:?Set TRIMMOMATIC_JAR to the Trimmomatic 0.39 JAR}"
: "${PICARD_JAR:?Set PICARD_JAR to the Picard 2.26.2 JAR}"
REPAIR_PL=${REPAIR_PL:-${SCRIPT_DIR}/../read_processing/rePair.pl}
NEXTERA_ADAPTERS=${NEXTERA_ADAPTERS:-NexteraPE-PE.fa}

ROW=$(awk -F '\t' -v task="$TASK_ID" 'NR == task + 1 {print; exit}' "$MANIFEST")
if [[ -z "$ROW" ]]; then
    echo "No manifest entry corresponds to array task $TASK_ID." >&2
    exit 1
fi
IFS=$'\t' read -r SAMPLE RAW_R1 RAW_R2 <<< "$ROW"

SAMPLE_DIR="$OUT_DIR/samples/$SAMPLE"
TRIM_DIR="$SAMPLE_DIR/trimmed_reads"
QC_DIR="$SAMPLE_DIR/fastqc"
CLEAN_DIR="$SAMPLE_DIR/clean_reads"
MAP_DIR="$SAMPLE_DIR/mapping"
HOST_DIR="$OUT_DIR/bams/host_final"
SYM_DIR="$OUT_DIR/bams/symbiont"
METRICS_DIR="$OUT_DIR/metrics"
mkdir -p "$TRIM_DIR" "$QC_DIR" "$CLEAN_DIR" "$MAP_DIR" \
    "$HOST_DIR" "$SYM_DIR" "$METRICS_DIR"

R1_LINK="$SAMPLE_DIR/${SAMPLE}_R1.fastq.gz"
R2_LINK="$SAMPLE_DIR/${SAMPLE}_R2.fastq.gz"
ln -sf "$(cd "$RAW_DIR" && pwd)/$RAW_R1" "$R1_LINK"
ln -sf "$(cd "$RAW_DIR" && pwd)/$RAW_R2" "$R2_LINK"

# Adapter and quality trimming recorded in the original workflow.
trim_galore \
    --nextseq 20 \
    --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --clip_R1 2 \
    --clip_R2 2 \
    --fastqc \
    --paired \
    --output_dir "$TRIM_DIR" \
    "$R1_LINK" "$R2_LINK"

TRIMGALORE_R1="$TRIM_DIR/${SAMPLE}_R1_val_1.fq.gz"
TRIMGALORE_R2="$TRIM_DIR/${SAMPLE}_R2_val_2.fq.gz"
TRIMMED_R1="$TRIM_DIR/${SAMPLE}_1.trimmed.fastq"
TRIMMED_R2="$TRIM_DIR/${SAMPLE}_2.trimmed.fastq"

java -jar "$TRIMMOMATIC_JAR" PE -threads 4 \
    "$TRIMGALORE_R1" "$TRIMGALORE_R2" \
    "$TRIMMED_R1" "$TRIM_DIR/${SAMPLE}_1un.trimmed.fastq" \
    "$TRIMMED_R2" "$TRIM_DIR/${SAMPLE}_2un.trimmed.fastq" \
    "ILLUMINACLIP:${NEXTERA_ADAPTERS}:2:30:10" \
    SLIDINGWINDOW:4:20 \
    MINLEN:20

fastqc --outdir "$QC_DIR" "$TRIMMED_R1" "$TRIMMED_R2"

# Apply the additional trimming and Q20/90% read filter, then restore pairing.
TRIMMED2_R1="$TRIM_DIR/${SAMPLE}_1.trimmed2.fastq"
TRIMMED2_R2="$TRIM_DIR/${SAMPLE}_2.trimmed2.fastq"
seqtk trimfq "$TRIMMED_R1" > "$TRIMMED2_R1"
seqtk trimfq "$TRIMMED_R2" > "$TRIMMED2_R2"

CLEAN_R1="$CLEAN_DIR/${SAMPLE}_R1.clean"
CLEAN_R2="$CLEAN_DIR/${SAMPLE}_R2.clean"
fastq_quality_filter -q 20 -p 90 -i "$TRIMMED2_R1" -o "$CLEAN_R1"
fastq_quality_filter -q 20 -p 90 -i "$TRIMMED2_R2" -o "$CLEAN_R2"

(
    cd "$CLEAN_DIR"
    perl "$REPAIR_PL" "${SAMPLE}_R1.clean" "${SAMPLE}_R2.clean"
)

PAIRED_R1="$CLEAN_DIR/R1_${SAMPLE}_R1.clean"
PAIRED_R2="$CLEAN_DIR/R2_${SAMPLE}_R2.clean"
MAPPED_SAM="$MAP_DIR/${SAMPLE}.paired.bwa.sam"
bwa mem -t "$THREADS" "$REFERENCE" "$PAIRED_R1" "$PAIRED_R2" > "$MAPPED_SAM"

# Parse host and symbiont alignments from the combined-reference SAM.
HOST_SAM="$MAP_DIR/${SAMPLE}.paired_HOST.sam"
SYM_SAM="$MAP_DIR/${SAMPLE}.paired_SYM.sam"
grep -v 'k127_' "$MAPPED_SAM" > "$HOST_SAM"
grep -v 'OZ03\|CAXIUN' "$MAPPED_SAM" > "$SYM_SAM"

HOST_SORTED="$MAP_DIR/${SAMPLE}.paired_sorted_HOST.bam"
SYM_SORTED="$SYM_DIR/${SAMPLE}.paired_sorted_SYM.bam"
samtools view -b -S "$HOST_SAM" \
    | samtools sort -@ "$THREADS" -O BAM -o "$HOST_SORTED"
samtools view -b -S "$SYM_SAM" \
    | samtools sort -@ "$THREADS" -O BAM -o "$SYM_SORTED"
samtools index "$SYM_SORTED"

# Add read groups, sort, fix mate information, mark duplicates, and index.
RG_BAM="$MAP_DIR/${SAMPLE}_host_RG.bam"
RG_SORTED="$MAP_DIR/${SAMPLE}_host_RG.sort.bam"
FIXED_BAM="$MAP_DIR/${SAMPLE}_host_RG.sort.fxm8s.bam"
FINAL_BAM="$HOST_DIR/${SAMPLE}_host_GATKd.bam"

java -Xmx4G -jar "$PICARD_JAR" AddOrReplaceReadGroups \
    I="$HOST_SORTED" \
    O="$RG_BAM" \
    RGID="$SAMPLE" \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM="$SAMPLE"

samtools sort -@ "$THREADS" -O BAM -o "$RG_SORTED" "$RG_BAM"

java -Xmx4G -jar "$PICARD_JAR" FixMateInformation \
    I="$RG_SORTED" \
    O="$FIXED_BAM"

java -Xmx4G -jar "$PICARD_JAR" MarkDuplicates \
    I="$FIXED_BAM" \
    O="$FINAL_BAM" \
    M="$METRICS_DIR/${SAMPLE}_host_metrics.txt"

samtools index "$FINAL_BAM"

