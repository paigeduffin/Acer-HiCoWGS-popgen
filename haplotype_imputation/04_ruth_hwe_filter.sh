#!/usr/bin/env bash
# Author: Paige Duffin
# Apply RUTH/HWE filtering to the 44-sample phased reference panel.
# Inputs: the merged WhatsHap VCF, the 44-sample keep list and PCA eigenvectors
# from step 03, and the A. cervicornis repeat-mask BED file.
# Software: bcftools 1.19, HTSlib 1.19.1, and RUTH 3.29.4.
# Set RUTH_BIN to the executable name or its path if needed.
# Usage: bash 04_ruth_hwe_filter.sh [PHASED_VCF.gz] [KEEP_FILE] [EIGENVEC] [REPEATMASK.bed] [NEW_OUTPUT_DIRECTORY]

set -euo pipefail

PHASED_VCF=${1:-whatshap_phasing_output/merged_step1phase_n46.vcf.gz}
KEEP_FILE=${2:-between_phasing_steps/after_remove_admixed/n44.keep.txt}
EIGENVEC=${3:-between_phasing_steps/after_remove_admixed/HWE_PD_n44_PCA_final.eigenvec}
REPEATMASK=${4:-Acervicornis.repeatmask.bed}
OUT_DIR=${5:-ruth_hwe_filtering}
RUTH_BIN=${RUTH_BIN:-ruth}

mkdir -- "$OUT_DIR"

# Prepare the 44-sample RUTH input: complete biallelic SNPs outside repeats,
# with minor allele count greater than three and refreshed INFO annotations.
bcftools view \
    -m2 -M2 -v snps \
    -i 'F_MISSING==0' \
    "$PHASED_VCF" \
    -Ou \
    | bcftools annotate --remove INFO -Ou \
    | bcftools view \
        --samples-file "$KEEP_FILE" \
        --targets-file "^${REPEATMASK}" \
        -Ou \
    | bcftools +fill-tags -Ou -- -t all \
    | bcftools view \
        --min-ac 3:minor \
        -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
        -Ov \
        -o "$OUT_DIR/ruth_input_n44.vcf"

# Reformat the PLINK eigenvectors for RUTH.
cut -f2- "$EIGENVEC" \
    | sed '1s/^IID/ID/' \
    > "$OUT_DIR/ruth_n44.eigenvec"

# Run RUTH using genotype likelihoods.
"$RUTH_BIN" \
    --evec "$OUT_DIR/ruth_n44.eigenvec" \
    --vcf "$OUT_DIR/ruth_input_n44.vcf" \
    --field PL \
    --site-only \
    --out "$OUT_DIR/ruth_output_n44.vcf"

# Retain positions with -5 < HWE_SLP_I < 5 and minor allele count > 3.
bcftools query \
    -f '%CHROM\t%POS\t%INFO/HWE_SLP_I\t%INFO/MAF\n' \
    "$OUT_DIR/ruth_output_n44.vcf" \
    > "$OUT_DIR/ruth_statistics_n44.tsv"

awk 'BEGIN {FS=OFS="\t"}
    $3 != "." && $4 != "." {
        logp=$3
        mac=$4*88
        if (logp > -5 && logp < 5 && mac > 3) print $1, $2
    }
' "$OUT_DIR/ruth_statistics_n44.tsv" \
    > "$OUT_DIR/ruth_pass_positions_n44.tsv"

# Apply the passing positions to the phased 44-sample panel.
bcftools view \
    -m2 -M2 -v snps \
    --samples-file "$KEEP_FILE" \
    -T "$OUT_DIR/ruth_pass_positions_n44.tsv" \
    "$PHASED_VCF" \
    -Ou \
    | bcftools annotate --remove INFO -Ou \
    | bcftools +fill-tags -Ou -- -t all \
    | bcftools view -Oz -o "$OUT_DIR/refpanel_n44.vcf.gz"

tabix -p vcf "$OUT_DIR/refpanel_n44.vcf.gz"
