#!/usr/bin/env bash
# Apply depth, quality, and allele-count filters, then retain biallelic SNPs.
# Run after ruth_hwe_filter.sh and before sor_filter.sh.
# Requirement: bcftools 1.19.
# Usage: bash quality_allele_filter.sh [INPUT.vcf.gz] [OUTPUT.vcf.gz]

set -euo pipefail

if [[ $# -gt 2 ]]; then
    echo "Usage: $0 [INPUT.vcf.gz] [OUTPUT.vcf.gz]" >&2
    exit 1
fi
INPUT=${1:-ruth_results/hwe_filtered.vcf.gz}
OUTPUT=${2:-full.panel_5bpind.morePAN.HWE.DPqual.AC3.bi.vcf.gz}

command -v bcftools >/dev/null || { echo "Missing program: bcftools" >&2; exit 1; }
[[ -r "$INPUT" ]] || { echo "Cannot read input: $INPUT" >&2; exit 1; }
for target in "$OUTPUT" "${OUTPUT}.tbi" "${OUTPUT}.csi"; do
    [[ ! -e "$target" ]] || { echo "Output already exists: $target" >&2; exit 1; }
done

# Apply the site-selection expression without changing individual genotypes.
# DP is the INFO field; AC > 3 is a strict alternate-allele count threshold.
bcftools view -i 'INFO/DP>10 && QUAL>20 && GQ>20 && MQ>40 && AC>3' \
    -Ou "$INPUT" \
    | bcftools view -m2 -M2 -v snps -Oz -o "$OUTPUT"
bcftools index -t "$OUTPUT"

echo "Input for sor_filter.sh: $OUTPUT"
