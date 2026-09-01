#!/usr/bin/env bash
# Author: Paige Duffin
# Call one LoCo VCF per sample, merge the original 329-sample set, then recall
# ID_DRTO_102 at the retained LoCo sites and add it to produce the final
# 330-sample imputation input.
# Software: bcftools 1.19, HTSlib 1.19.1, and samtools.
#
# Call mode (submit as a 1-330 Slurm array):
#   sbatch --array=1-330 07b_call_and_merge_loco_vcfs.sh call [BAM_DIRECTORY] [COMBINED_REFERENCE.fa] [OUTPUT_DIRECTORY] [MANIFEST]
#
# Merge mode (run after all call tasks finish):
#   bash 07b_call_and_merge_loco_vcfs.sh merge [BAM_DIRECTORY] [COMBINED_REFERENCE.fa] [OUTPUT_DIRECTORY] [MANIFEST]

set -euo pipefail

if [[ $# -lt 1 || $# -gt 5 ]]; then
    echo "Usage: $0 {call|merge} [BAM_DIRECTORY] [COMBINED_REFERENCE.fa] [OUTPUT_DIRECTORY] [MANIFEST]" >&2
    exit 1
fi

MODE=$1
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BAM_DIR=${2:-loco_read_processing/bams/host_final}
REFERENCE=${3:-ref_genome_indexed/GCA_964034985.1_Acer_S.fitti.fa}
OUT_DIR=${4:-loco_variant_calling}
MANIFEST=${5:-${SCRIPT_DIR}/loco_sample_manifest.tsv}
VCF_DIR="$OUT_DIR/per_sample_vcfs"
mkdir -p "$VCF_DIR"

call_sample() {
    local TASK_ID=${SLURM_ARRAY_TASK_ID:?Submit call mode as a 1-330 Slurm array}
    local SAMPLE
    SAMPLE=$(awk -F '\t' -v task="$TASK_ID" 'NR == task + 1 {print $1; exit}' "$MANIFEST")
    if [[ -z "$SAMPLE" ]]; then
        echo "No manifest entry corresponds to array task $TASK_ID." >&2
        exit 1
    fi

    local BAM="$BAM_DIR/${SAMPLE}_host_GATKd.bam"
    local VCF="$VCF_DIR/${SAMPLE}.vcf.gz"
    if [[ ! -s "${BAM}.bai" ]]; then
        samtools index "$BAM"
    fi

    bcftools mpileup \
        -q 40 \
        -Q 20 \
        -f "$REFERENCE" \
        -Ou \
        -A \
        "$BAM" \
        | bcftools call -m -Oz -o "$VCF"
    tabix -p vcf "$VCF"
}

merge_samples() {
    local MERGE_DIR="$OUT_DIR/merge"
    local PART_DIR="$MERGE_DIR/parts"
    local N329_DIR="$MERGE_DIR/n329_chromosomes"
    local DRTO_DIR="$MERGE_DIR/DRTO_102"
    local FINAL_DIR="$MERGE_DIR/final_330_chromosomes"
    local N329_LIST="$MERGE_DIR/n329_vcfs.list"
    local PART_OUTPUT_LIST="$MERGE_DIR/part_outputs.list"
    local N329_CHROM_LIST="$MERGE_DIR/n329_chromosomes.list"
    local FINAL_CHROM_LIST="$MERGE_DIR/final_330_chromosomes.list"
    local MERGED_N329="$MERGE_DIR/Acer_LoCo_n329_merged.vcf.gz"
    local LOCO_SITES="$MERGE_DIR/LoCo_sites_all.vcf.gz"
    local DRTO_BAM="$BAM_DIR/ID_DRTO_102_host_GATKd.bam"
    local DRTO_ALL="$DRTO_DIR/DRTO_102_atLoCoSites.vcf.gz"

    mkdir -p "$PART_DIR" "$N329_DIR" "$DRTO_DIR" "$FINAL_DIR"

    # Reproduce the original 329-sample merge by excluding ID_DRTO_102.
    awk -F '\t' -v dir="$VCF_DIR" \
        'NR > 1 && $1 != "ID_DRTO_102" {print dir "/" $1 ".vcf.gz"}' \
        "$MANIFEST" \
        | LC_ALL=C sort \
        > "$N329_LIST"

    if [[ $(wc -l < "$N329_LIST") -ne 329 ]]; then
        echo "Expected 329 VCFs after excluding ID_DRTO_102." >&2
        exit 1
    fi
    while IFS= read -r VCF; do
        [[ -s "$VCF" ]] || { echo "Missing VCF: $VCF" >&2; exit 1; }
    done < "$N329_LIST"

    # The original workflow used four intermediate merges to reduce memory use.
    split -d -l 90 --additional-suffix=.list "$N329_LIST" "$PART_DIR/part_"
    : > "$PART_OUTPUT_LIST"
    for LIST in "$PART_DIR"/part_*.list; do
        PART=$(basename "$LIST" .list)
        PART_VCF="$PART_DIR/${PART}.vcf.gz"
        bcftools merge --file-list "$LIST" -Oz -o "$PART_VCF"
        tabix -p vcf "$PART_VCF"
        printf '%s\n' "$PART_VCF" >> "$PART_OUTPUT_LIST"
    done

    bcftools merge --force-samples --file-list "$PART_OUTPUT_LIST" \
        -Oz -o "$MERGED_N329"
    tabix -p vcf "$MERGED_N329"

    : > "$N329_CHROM_LIST"
    for CHROM_NUMBER in {966..979}; do
        CHROM="OZ035${CHROM_NUMBER}.1"
        READY="$N329_DIR/${CHROM}_noMulti.vcf.gz"

        bcftools view -r "$CHROM" "$MERGED_N329" -Ou \
            | bcftools norm -f "$REFERENCE" -c ws -Ou \
            | bcftools view -m2 -M2 -v snps -Ou \
            | bcftools annotate --set-id '%CHROM_%POS_%REF_%ALT' -Ou \
            | bcftools +fill-tags -Ou -- -t AC,AN,AF \
            | bcftools norm -d none -Ou \
            | bcftools view \
                -e 'INFO/AC<3 | (INFO/AN-INFO/AC)<3' \
                -Oz -o "$READY"
        tabix -p vcf "$READY"
        printf '%s\n' "$READY" >> "$N329_CHROM_LIST"
    done

    # Use the retained 329-sample sites to recall ID_DRTO_102, preserving 0/0 calls.
    bcftools concat --file-list "$N329_CHROM_LIST" -Oz -o "$LOCO_SITES"
    tabix -p vcf "$LOCO_SITES"
    if [[ ! -s "${DRTO_BAM}.bai" ]]; then
        samtools index "$DRTO_BAM"
    fi
    bcftools mpileup \
        -q 40 \
        -Q 20 \
        -f "$REFERENCE" \
        -R "$LOCO_SITES" \
        -Ou \
        "$DRTO_BAM" \
        | bcftools call -m -Oz -o "$DRTO_ALL"
    tabix -p vcf "$DRTO_ALL"

    : > "$FINAL_CHROM_LIST"
    for CHROM_NUMBER in {966..979}; do
        CHROM="OZ035${CHROM_NUMBER}.1"
        DRTO_CHROM="$DRTO_DIR/${CHROM}_DRTO_102_atSites.vcf.gz"
        MERGED_CHROM="$FINAL_DIR/${CHROM}_LoCo_plus_DRTO102.vcf.gz"
        FINAL_CHROM="$FINAL_DIR/${CHROM}_LoCo_plus_DRTO102_snps.biallelic.vcf.gz"

        bcftools view -r "$CHROM" -Oz -o "$DRTO_CHROM" "$DRTO_ALL"
        tabix -p vcf "$DRTO_CHROM"
        bcftools merge \
            "$DRTO_CHROM" \
            "$N329_DIR/${CHROM}_noMulti.vcf.gz" \
            -Oz -o "$MERGED_CHROM"
        tabix -p vcf "$MERGED_CHROM"

        bcftools view -m2 -M2 -v snps "$MERGED_CHROM" -Ou \
            | bcftools +fill-tags -Ou -- -t AC,AN,AF \
            | bcftools norm -d none -Oz -o "$FINAL_CHROM"
        tabix -p vcf "$FINAL_CHROM"
        printf '%s\n' "$FINAL_CHROM" >> "$FINAL_CHROM_LIST"
    done

    bcftools concat --file-list "$FINAL_CHROM_LIST" \
        -Oz -o "$OUT_DIR/LoCo_plus_DRTO102_merge.vcf.gz"
    tabix -p vcf "$OUT_DIR/LoCo_plus_DRTO102_merge.vcf.gz"
}

case "$MODE" in
    call)
        call_sample
        ;;
    merge)
        merge_samples
        ;;
    *)
        echo "MODE must be call or merge." >&2
        exit 1
        ;;
esac

