#!/usr/bin/env bash
# Author: Paige Duffin
# Prepare the 330-sample LoCo VCF and perform the two-round imputation
# workflow used for haplotype-panel validation.
# Round 1: Beagle 4.1 genotype-probability estimation; set GP<0.99 calls
# to missing. Round 2: Beagle 5.5 imputation; retain sites with DR2>=0.99.
# This reproduces the analysis as performed with the complete n=44 panel.
# Software: bcftools 1.19, HTSlib 1.19.1, Java 17,
# Beagle 4.1 (beagle.27Jan18.7e1.jar), and Beagle 5.5
# (beagle.27Feb25.75f.jar).
# Set BEAGLE41_JAR and BEAGLE55_JAR if the jar paths differ.
# Usage: bash 07_two_round_loco_imputation.sh LOCO_330.vcf.gz COMBINED_REFERENCE.fa [REFERENCE_DIR] [GENETIC_MAP] [NEW_OUTPUT_DIRECTORY]

set -euo pipefail

if [[ $# -lt 2 || $# -gt 5 ]]; then
    echo "Usage: $0 LOCO_330.vcf.gz COMBINED_REFERENCE.fa [REFERENCE_DIR] [GENETIC_MAP] [NEW_OUTPUT_DIRECTORY]" >&2
    exit 1
fi

LOCO_VCF=$1
COMBINED_REFERENCE=$2
REFERENCE_DIR=${3:-imputation_reference}
GENETIC_MAP=${4:-Acervicornis_SEXAVG.map}
OUT_DIR=${5:-loco_imputation}
BEAGLE41_JAR=${BEAGLE41_JAR:-beagle.27Jan18.7e1.jar}
BEAGLE55_JAR=${BEAGLE55_JAR:-beagle.27Feb25.75f.jar}

INPUT_DIR="$OUT_DIR/input_chromosomes"
ROUND1_DIR="$OUT_DIR/round1_beagle41"
RECAL_DIR="$OUT_DIR/gp_recalibrated"
ROUND2_DIR="$OUT_DIR/round2_beagle55"
DR2_DIR="$OUT_DIR/dr2_filtered"

mkdir -- "$OUT_DIR"
mkdir -p -- "$INPUT_DIR" "$ROUND1_DIR" "$RECAL_DIR" "$ROUND2_DIR" "$DR2_DIR"

RAW_LIST="$OUT_DIR/round1_vcfs.list"
RECAL_LIST="$OUT_DIR/recalibrated_vcfs.list"
IMPUTED_LIST="$OUT_DIR/imputed_preDR2_vcfs.list"
DR2_LIST="$OUT_DIR/dr2_filtered_vcfs.list"
: > "$RAW_LIST"
: > "$RECAL_LIST"
: > "$IMPUTED_LIST"
: > "$DR2_LIST"

# Imputation validation used the 14 nuclear chromosomes below.
for CHROM_NUMBER in {966..979}; do
    CHROM="OZ035${CHROM_NUMBER}.1"
    INPUT_VCF="$INPUT_DIR/${CHROM}_LoCo_330_ready.vcf.gz"
    ROUND1_PREFIX="$ROUND1_DIR/${CHROM}_imputation1"
    ROUND1_VCF="${ROUND1_PREFIX}.vcf.gz"
    RECAL_VCF="$RECAL_DIR/${CHROM}_imputation1_GP0.99.vcf.gz"
    ROUND2_PREFIX="$ROUND2_DIR/${CHROM}_GP0.99_imputation2"
    ROUND2_VCF="${ROUND2_PREFIX}.vcf.gz"
    DR2_VCF="$DR2_DIR/${CHROM}_imputation2_DR2_0.99.vcf.gz"

    # Normalize to the mapping reference, retain biallelic SNPs, assign
    # standardized IDs, refresh allele tags, and exclude alleles with AC<3.
    bcftools view -r "$CHROM" "$LOCO_VCF" -Ou \
        | bcftools norm -f "$COMBINED_REFERENCE" -c ws -Ou \
        | bcftools view -m2 -M2 -v snps -Ou \
        | bcftools annotate --set-id '%CHROM_%POS_%REF_%ALT' -Ou \
        | bcftools +fill-tags -Ou -- -t AC,AN,AF \
        | bcftools view \
            -e 'INFO/AC<3 | (INFO/AN-INFO/AC)<3' \
            -Oz -o "$INPUT_VCF"
    tabix -p vcf "$INPUT_VCF"

    # Round 1: estimate genotype probabilities from LoCo likelihoods.
    java -Xss5m -Xmx16g \
        -jar "$BEAGLE41_JAR" \
        gl="$INPUT_VCF" \
        gprobs=true \
        ref="$REFERENCE_DIR/${CHROM}_reference_n44.bref" \
        map="$GENETIC_MAP" \
        window=2000 \
        overlap=200 \
        out="$ROUND1_PREFIX"
    tabix -p vcf "$ROUND1_VCF"
    printf '%s\n' "$ROUND1_VCF" >> "$RAW_LIST"

    # Set genotypes lacking a maximum GP of at least 0.99 to missing.
    bcftools +setGT "$ROUND1_VCF" \
        -Oz -o "$RECAL_VCF" -- \
        -t q -n . -e 'FORMAT/GP>=0.99'
    tabix -p vcf "$RECAL_VCF"
    printf '%s\n' "$RECAL_VCF" >> "$RECAL_LIST"

    # Round 2: impute missing genotypes from the GP-filtered calls.
    java -Xmx20g \
        -jar "$BEAGLE55_JAR" \
        gt="$RECAL_VCF" \
        ref="$REFERENCE_DIR/${CHROM}_reference_n44.bref3" \
        map="$GENETIC_MAP" \
        out="$ROUND2_PREFIX" \
        ne=10000 \
        impute=true \
        gp=true \
        seed=-99999
    tabix -p vcf "$ROUND2_VCF"
    printf '%s\n' "$ROUND2_VCF" >> "$IMPUTED_LIST"

    # Retain sites meeting the selected Beagle dosage-R-squared threshold.
    bcftools filter \
        -e 'INFO/DR2<0.99' \
        -Oz -o "$DR2_VCF" \
        "$ROUND2_VCF"
    tabix -p vcf "$DR2_VCF"
    printf '%s\n' "$DR2_VCF" >> "$DR2_LIST"
done

# Concatenate chromosome outputs for validation and downstream use.
bcftools concat --file-list "$RAW_LIST" \
    -Oz -o "$OUT_DIR/LoCo_RAW_with_GP_merged_chromosomes.vcf.gz"
tabix -p vcf "$OUT_DIR/LoCo_RAW_with_GP_merged_chromosomes.vcf.gz"

bcftools concat --file-list "$RECAL_LIST" \
    -Oz -o "$OUT_DIR/LoCo_RECAL_merged_chromosomes.vcf.gz"
tabix -p vcf "$OUT_DIR/LoCo_RECAL_merged_chromosomes.vcf.gz"

bcftools concat --file-list "$IMPUTED_LIST" \
    -Oz -o "$OUT_DIR/LoCo_imputed_n44_preDR2.vcf.gz"
tabix -p vcf "$OUT_DIR/LoCo_imputed_n44_preDR2.vcf.gz"

bcftools concat --file-list "$DR2_LIST" \
    -Oz -o "$OUT_DIR/LoCo_imputed_n44_DR2_0.99.vcf.gz"
tabix -p vcf "$OUT_DIR/LoCo_imputed_n44_DR2_0.99.vcf.gz"
