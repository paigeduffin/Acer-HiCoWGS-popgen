#!/usr/bin/env bash
# Author: Paige Duffin
# Perform population-based phasing of the 44-sample RUTH-filtered reference
# panel with Beagle, then concatenate the phased chromosome outputs.
# Inputs: the step-04 reference panel and chromosome-specific genetic maps.
# Software: bcftools 1.19, HTSlib 1.19.1, Java 17, and Beagle.
# Set BEAGLE_JAR to the Beagle .jar path if it is not named beagle.jar.
# Usage: bash 05_beagle_population_phasing.sh [REFPANEL.vcf.gz] [GENETIC_MAP_DIRECTORY] [NEW_OUTPUT_DIRECTORY]

set -euo pipefail

REFPANEL=${1:-ruth_hwe_filtering/refpanel_n44.vcf.gz}
MAP_DIR=${2:-gen.map_sep.scaffs}
OUT_DIR=${3:-beagle_phasing_output}
BEAGLE_JAR=${BEAGLE_JAR:-beagle.jar}

PARSED_DIR="$OUT_DIR/parsed_chromosomes"
PHASED_DIR="$OUT_DIR/phased_chromosomes"
PHASED_LIST="$OUT_DIR/phased_vcfs.list"
SKIPPED="$OUT_DIR/skipped_chromosomes.txt"

mkdir -- "$OUT_DIR"
mkdir -p -- "$PARSED_DIR" "$PHASED_DIR"
: > "$PHASED_LIST"
: > "$SKIPPED"

for CHROM_NUMBER in {966..980}; do
    CHROM="OZ035${CHROM_NUMBER}.1"
    INPUT_CHROM="$PARSED_DIR/${CHROM}_only_n44.vcf.gz"
    MAP_FILE="$MAP_DIR/gen.map_${CHROM}.map"
    OUTPUT_PREFIX="$PHASED_DIR/${CHROM}_n44_beagle_phased"

    bcftools view \
        -r "$CHROM" \
        -Oz \
        -o "$INPUT_CHROM" \
        "$REFPANEL"
    tabix -p vcf "$INPUT_CHROM"

    # Skip chromosomes without retained variants or a corresponding map.
    if [[ $(bcftools index -n "$INPUT_CHROM") -eq 0 || ! -s "$MAP_FILE" ]]; then
        printf '%s\n' "$CHROM" >> "$SKIPPED"
        continue
    fi

    java -Xmx32g -jar "$BEAGLE_JAR" \
        gt="$INPUT_CHROM" \
        map="$MAP_FILE" \
        out="$OUTPUT_PREFIX" \
        gp=true \
        ne=10000

    tabix -p vcf "${OUTPUT_PREFIX}.vcf.gz"
    printf '%s\n' "${OUTPUT_PREFIX}.vcf.gz" >> "$PHASED_LIST"
done

bcftools concat \
    --file-list "$PHASED_LIST" \
    -Oz \
    -o "$OUT_DIR/merged_beagle_phased_n44.vcf.gz"

tabix -p vcf "$OUT_DIR/merged_beagle_phased_n44.vcf.gz"
bcftools stats "$OUT_DIR/merged_beagle_phased_n44.vcf.gz" \
    > "$OUT_DIR/merged_beagle_phased_n44.stats.txt"
