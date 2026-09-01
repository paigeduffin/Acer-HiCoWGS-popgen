#!/usr/bin/env bash
# Author: Paige Duffin
# Extract matching HiCo and LoCo genotype records for the five validation pairs.
# Raw, GP-recalibrated, and newly imputed LoCo datasets are handled separately.
# Newly imputed records are identified using INFO/IMP=1 before DR2 filtering.
# Software: bcftools 1.19 and HTSlib 1.19.1.
# Usage: bash 08a_extract_validation_genotypes.sh HICO.vcf.gz [RAW.vcf.gz] [RECAL.vcf.gz] [IMPUTED_PRE_DR2.vcf.gz] [NEW_OUTPUT_DIRECTORY]

set -euo pipefail

if [[ $# -lt 1 || $# -gt 5 ]]; then
    echo "Usage: $0 HICO.vcf.gz [RAW.vcf.gz] [RECAL.vcf.gz] [IMPUTED_PRE_DR2.vcf.gz] [NEW_OUTPUT_DIRECTORY]" >&2
    exit 1
fi

HICO_VCF=$1
RAW_VCF=${2:-loco_imputation/LoCo_RAW_with_GP_merged_chromosomes.vcf.gz}
RECAL_VCF=${3:-loco_imputation/LoCo_RECAL_merged_chromosomes.vcf.gz}
IMPUTED_VCF=${4:-loco_imputation/LoCo_imputed_n44_preDR2.vcf.gz}
OUT_DIR=${5:-validation_genotypes}

HICO_SAMPLES='CRF_Acer-059,Mote_AC75,Mote_AC76,Mote_AC80,DRTO_114'
LOCO_SAMPLES='ID_CRF_Acer59,ID_MML_ML75,ID_MML_ML76,ID_MML_ML80,ID_DRTO_114'

mkdir -- "$OUT_DIR"

printf 'HiCo\tLoCo\n' > "$OUT_DIR/validation_sample_pairs.tsv"
printf '%s\t%s\n' \
    'CRF_Acer-059' 'ID_CRF_Acer59' \
    'Mote_AC75' 'ID_MML_ML75' \
    'Mote_AC76' 'ID_MML_ML76' \
    'Mote_AC80' 'ID_MML_ML80' \
    'DRTO_114' 'ID_DRTO_114' \
    >> "$OUT_DIR/validation_sample_pairs.tsv"

HICO_SUBSET="$OUT_DIR/HiCo_validation_samples.vcf.gz"
bcftools view \
    --samples "$HICO_SAMPLES" \
    -Oz -o "$HICO_SUBSET" \
    "$HICO_VCF"
tabix -p vcf "$HICO_SUBSET"

bcftools query -f '%CHROM\t%POS\n' "$HICO_SUBSET" \
    | LC_ALL=C sort -u \
    > "$OUT_DIR/HiCo_sites.tsv"

extract_overlap_stage() {
    local STAGE=$1
    local LOCO_INPUT=$2
    local LOCO_SUBSET="$OUT_DIR/${STAGE}_LoCo_validation_samples.vcf.gz"
    local LOCO_SITES="$OUT_DIR/${STAGE}_LoCo_sites.tsv"
    local OVERLAP="$OUT_DIR/${STAGE}_overlap_sites.tsv"
    local HICO_OVERLAP="$OUT_DIR/${STAGE}_HiCo_overlap.vcf.gz"
    local LOCO_OVERLAP="$OUT_DIR/${STAGE}_LoCo_overlap.vcf.gz"

    bcftools view \
        --samples "$LOCO_SAMPLES" \
        -Oz -o "$LOCO_SUBSET" \
        "$LOCO_INPUT"

    bcftools query -f '%CHROM\t%POS\n' "$LOCO_SUBSET" \
        | LC_ALL=C sort -u \
        > "$LOCO_SITES"

    comm -12 "$OUT_DIR/HiCo_sites.tsv" "$LOCO_SITES" > "$OVERLAP"

    bcftools view -T "$OVERLAP" -Oz -o "$HICO_OVERLAP" "$HICO_SUBSET"
    bcftools view -T "$OVERLAP" -Oz -o "$LOCO_OVERLAP" "$LOCO_SUBSET"

    printf 'CHROM\tPOS\tID\tSAMPLE\tGT\tDP\n' \
        > "$OUT_DIR/${STAGE}_HiCo_genotypes.tsv"
    bcftools query \
        -f '[%CHROM\t%POS\t%ID\t%SAMPLE\t%GT\t%DP\n]' \
        "$HICO_OVERLAP" \
        >> "$OUT_DIR/${STAGE}_HiCo_genotypes.tsv"

    printf 'CHROM\tPOS\tID\tSAMPLE\tGT\n' \
        > "$OUT_DIR/${STAGE}_LoCo_genotypes.tsv"
    bcftools query \
        -f '[%CHROM\t%POS\t%ID\t%SAMPLE\t%GT\n]' \
        "$LOCO_OVERLAP" \
        >> "$OUT_DIR/${STAGE}_LoCo_genotypes.tsv"
}

extract_overlap_stage RAW "$RAW_VCF"
extract_overlap_stage RECAL "$RECAL_VCF"

# Restrict the imputed comparison to newly imputed sites (INFO/IMP=1).
IMPUTED_SUBSET="$OUT_DIR/IMP_LoCo_validation_samples.vcf.gz"
IMPUTED_SITES="$OUT_DIR/IMP_sites.tsv"
IMPUTED_HICO="$OUT_DIR/IMP_HiCo_overlap.vcf.gz"
IMPUTED_LOCO="$OUT_DIR/IMP_LoCo_overlap.vcf.gz"

bcftools view \
    --samples "$LOCO_SAMPLES" \
    -Oz -o "$IMPUTED_SUBSET" \
    "$IMPUTED_VCF"

bcftools view -i 'INFO/IMP=1' "$IMPUTED_SUBSET" -Ou \
    | bcftools query -f '%CHROM\t%POS\n' \
    | LC_ALL=C sort -u \
    > "$IMPUTED_SITES"

bcftools view -T "$IMPUTED_SITES" -Oz -o "$IMPUTED_HICO" "$HICO_SUBSET"
bcftools view -T "$IMPUTED_SITES" -i 'INFO/IMP=1' \
    -Oz -o "$IMPUTED_LOCO" "$IMPUTED_SUBSET"

printf 'CHROM\tPOS\tID\tSAMPLE\tGT\tDP\n' \
    > "$OUT_DIR/IMP_HiCo_genotypes.tsv"
bcftools query \
    -f '[%CHROM\t%POS\t%ID\t%SAMPLE\t%GT\t%DP\n]' \
    "$IMPUTED_HICO" \
    >> "$OUT_DIR/IMP_HiCo_genotypes.tsv"

printf 'CHROM\tPOS\tID\tDR2\tSAMPLE\tGT\n' \
    > "$OUT_DIR/IMP_LoCo_genotypes.tsv"
bcftools query \
    -f '[%CHROM\t%POS\t%ID\t%INFO/DR2\t%SAMPLE\t%GT\n]' \
    "$IMPUTED_LOCO" \
    >> "$OUT_DIR/IMP_LoCo_genotypes.tsv"
