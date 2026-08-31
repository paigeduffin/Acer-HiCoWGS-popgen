#!/usr/bin/env bash
# Author: Paige Duffin
# Prepare the expanded 46-sample VCF for read-backed phasing with WhatsHap.
# Input: the SNP VCF with SOR annotations, before the filters below.
# Sample names must match the corresponding BAM read-group sample names.
# Software: bcftools 1.19, including the tag2tag plugin.
# USC HPC modules: ml gcc/13.3.0 openblas/0.3.28 bcftools/1.19 htslib/1.19.1
# Usage: bash 01_prepare_whatshap_input.sh [INPUT.vcf.gz] [NEW_OUTPUT_DIRECTORY]

set -euo pipefail

INPUT=${1:-snps.5bpind.merged.morePAN.v2_SOR.add.vcf.gz}
OUT_DIR=${2:-whatshap_input}
mkdir -- "$OUT_DIR"
mkdir -- "$OUT_DIR/sep_by_samples"

PREP_VCF="$OUT_DIR/snps.5bpind.DPqual.AC3.nomiss.SOR.decomp2bi_n46.vcf.gz"
PREP_BCF="$OUT_DIR/snps.5bpind.DPqual.AC3.nomiss.SOR.decomp2bi.pl2gl_n46.bcf"

# Apply the recorded quality, allele-count, missingness, and SOR filters,
# then split multiallelic records into biallelic records.
bcftools view -i 'INFO/DP>10 && QUAL>20 && GQ>20 && MQ>40 && FS<10 && AC>3' \
    "$INPUT" -Ou \
    | bcftools view -i 'F_MISSING=0' -Ou \
    | bcftools view -i 'SOR<4' -Ou \
    | bcftools norm -m -any -Oz -o "$PREP_VCF"

# Convert PL to GL for the per-sample phasing inputs.
bcftools +tag2tag "$PREP_VCF" -- --pl-to-gl \
    | bcftools view -Ob -o "$PREP_BCF"

# Extract the sample list directly from the prepared VCF, then split and index.
bcftools query -l "$PREP_VCF" > "$OUT_DIR/samples.txt"
while IFS= read -r SAMPLE; do
    bcftools view -s "$SAMPLE" -Ob \
        -o "$OUT_DIR/sep_by_samples/${SAMPLE}_for.phasing.bcf" "$PREP_BCF"
    bcftools index "$OUT_DIR/sep_by_samples/${SAMPLE}_for.phasing.bcf"
done < "$OUT_DIR/samples.txt"
