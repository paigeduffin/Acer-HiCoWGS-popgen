#!/usr/bin/env bash
# Author: Paige Duffin
# Prepare the complete n=44 haplotype panel for the two-round LoCo
# imputation workflow used in the manuscript validation.
# This reproduces the analysis as performed; leave-one-out panels were not used.
# Inputs: the final phased panel and the host reference FASTA.
# Software: bcftools 1.19, HTSlib 1.19.1, Java 17,
# bref.27Jan18.7e1.jar, and bref3.27Feb25.75f.jar.
# Set BREF41_JAR and BREF3_JAR if the converter jars have different paths.
# Usage: bash 06_prepare_imputation_reference.sh [PANEL.vcf.gz] [HOST_REFERENCE.fa] [NEW_OUTPUT_DIRECTORY]

set -euo pipefail

PANEL=${1:-beagle_phasing_output/merged_beagle_phased_n44.vcf.gz}
HOST_REFERENCE=${2:-ref_genome_indexed/GCA_964034985.1_jaAcrCerv1.1_genomic.fna}
OUT_DIR=${3:-imputation_reference}
BREF41_JAR=${BREF41_JAR:-bref.27Jan18.7e1.jar}
BREF3_JAR=${BREF3_JAR:-bref3.27Feb25.75f.jar}

mkdir -- "$OUT_DIR"

# Validation and imputation used the 14 nuclear chromosomes below.
for CHROM_NUMBER in {966..979}; do
    CHROM="OZ035${CHROM_NUMBER}.1"
    NORMALIZED="$OUT_DIR/${CHROM}_reference_n44.normalized.vcf.gz"
    DEDUPLICATED="$OUT_DIR/${CHROM}_reference_n44.deduplicated.vcf.gz"
    FINAL_VCF="$OUT_DIR/${CHROM}_reference_n44.vcf.gz"
    DUPLICATE_IDS="$OUT_DIR/${CHROM}_duplicate_ids.txt"

    # Recalculate allele tags, retain alleles represented at least three
    # times, normalize, retain complete biallelic records, and standardize IDs.
    bcftools view -r "$CHROM" "$PANEL" -Ou \
        | bcftools +fill-tags -Ou -- -t AC,AN,AF \
        | bcftools view -e 'INFO/AC<3 | (INFO/AN-INFO/AC)<3' -Ou \
        | bcftools norm -m -any -Ou \
        | bcftools norm -f "$HOST_REFERENCE" -d none -Ou \
        | bcftools view -m2 -M2 -g ^miss -Ou \
        | bcftools annotate \
            --set-id '%CHROM_%POS_%REF_%ALT' \
            -Oz -o "$NORMALIZED"

    bcftools query -f '%ID\n' "$NORMALIZED" \
        | sort \
        | uniq -d \
        > "$DUPLICATE_IDS"

    if [[ -s "$DUPLICATE_IDS" ]]; then
        bcftools view \
            -e "ID=@${DUPLICATE_IDS}" \
            -Oz -o "$DEDUPLICATED" \
            "$NORMALIZED"
    else
        cp -- "$NORMALIZED" "$DEDUPLICATED"
    fi

    bcftools +fill-tags "$DEDUPLICATED" \
        -Oz -o "$FINAL_VCF" -- -t AC,AN,AF
    tabix -p vcf "$FINAL_VCF"

    # Beagle 4.1 converter writes the .bref file beside the input VCF.
    java -jar "$BREF41_JAR" "$FINAL_VCF"

    # Beagle 5.5 uses the bref3 reference format.
    java -jar "$BREF3_JAR" "$FINAL_VCF" \
        > "$OUT_DIR/${CHROM}_reference_n44.bref3"
done
