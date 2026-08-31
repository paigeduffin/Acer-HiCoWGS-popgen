#!/usr/bin/env bash
# Build the neutral SNP panel by removing outliers, pruning LD, and
# retaining sites with <=10% missing genotypes and MAF >=0.05.
#
# Input: the full filtered SNP panel and merged PCAdapt/OutFLANK positions.
# Outlier positions: two tab-separated columns, CHROM and POS (1-based).
# Software: bcftools 1.19, VCFtools 0.1.16, HTSlib 1.19.1,
#           PLINK2 2.00a4.3. Load the software before running this script.
#
# Usage: bash build_neutral_panel.sh [INPUT.vcf.gz] [OUTLIERS.tsv] [OUTPUT_DIR]
# With no arguments, the filenames below are used in the current directory.

set -euo pipefail

if [[ $# -gt 3 ]]; then
    echo "Usage: $0 [INPUT.vcf.gz] [OUTLIERS.tsv] [OUTPUT_DIR]" >&2
    exit 1
fi

INPUT=${1:-full.panel_5bpind.morePAN.HWE.DPqual.AC3.bi.SOR.vcf.gz}
OUTLIERS=${2:-merged.outlier.snps_sept.24.2025_n46.tsv}
OUT_DIR=${3:-neutral_panel}

for program in bcftools vcftools bgzip tabix plink2; do
    command -v "$program" >/dev/null || {
        echo "Missing program: $program" >&2
        exit 1
    }
done
for input_file in "$INPUT" "$OUTLIERS"; do
    [[ -s "$input_file" ]] || { echo "Missing or empty input: $input_file" >&2; exit 1; }
done
# Use a new directory so previous results are not overwritten.
mkdir -- "$OUT_DIR"

ONLY_OUTLIERS="$OUT_DIR/Acer_popgen_ONLY.outliers_n46"
NO_OUTLIERS="$OUT_DIR/Acer_popgen_rm.outliers_n46"
LD_PREFIX="$OUT_DIR/Acer_popgen_rm.outliers_n46_LD.prune_50_5_0.2"
NEUTRAL="$OUT_DIR/Acer_main.NEU.PANEL_n46"

# 1. Export outlier-only SNPs and exclude those positions from the full panel.
vcftools --gzvcf "$INPUT" --positions "$OUTLIERS" \
    --recode --recode-INFO-all --out "$ONLY_OUTLIERS"
vcftools --gzvcf "$INPUT" --exclude-positions "$OUTLIERS" \
    --recode --recode-INFO-all --out "$NO_OUTLIERS"

for prefix in "$ONLY_OUTLIERS" "$NO_OUTLIERS"; do
    bgzip -c "${prefix}.recode.vcf" > "${prefix}.vcf.gz"
    tabix -p vcf "${prefix}.vcf.gz"
done

# 2. Assign chromosome_position IDs and identify SNPs to retain after LD
# pruning. Window: 50 SNPs; step: 5 SNPs; r-squared threshold: 0.2.
# No additional physical-distance thinning is applied.
bcftools annotate --set-id '%CHROM\_%POS' \
    -Ov -o "${NO_OUTLIERS}_unique.ids.vcf" "${NO_OUTLIERS}.vcf.gz"
plink2 --vcf "${NO_OUTLIERS}_unique.ids.vcf" \
    --indep-pairwise 50 5 0.2 --allow-extra-chr --bad-ld \
    --out "$LD_PREFIX"

bcftools view --include "ID=@${LD_PREFIX}.prune.in" \
    -Oz -o "${NEUTRAL}.vcf.gz" "${NO_OUTLIERS}_unique.ids.vcf"
tabix -p vcf "${NEUTRAL}.vcf.gz"

# 3. Retain sites with at least 90% called genotypes.
vcftools --gzvcf "${NEUTRAL}.vcf.gz" --max-missing 0.9 \
    --recode --recode-INFO-all --out "${NEUTRAL}.max10miss"
bgzip -c "${NEUTRAL}.max10miss.recode.vcf" > "${NEUTRAL}.max10miss.vcf.gz"
tabix -p vcf "${NEUTRAL}.max10miss.vcf.gz"

# 4. Recalculate MAF and retain SNPs with MAF >=0.05.
bcftools +fill-tags "${NEUTRAL}.max10miss.vcf.gz" -Ou -- -t MAF \
    | bcftools view -i 'MAF>=0.05' -Oz \
        -o "${NEUTRAL}.max10miss.MAF05.vcf.gz"
tabix -p vcf "${NEUTRAL}.max10miss.MAF05.vcf.gz"

echo "Final neutral panel: ${NEUTRAL}.max10miss.MAF05.vcf.gz"
