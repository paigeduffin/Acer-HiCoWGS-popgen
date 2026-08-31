#!/usr/bin/env bash
# Rename chromosomes and create PED/MAP files for 02b_snmf.R.
# Software: bcftools 1.19, HTSlib 1.19.1, VCFtools 0.1.16.
# Load the software before running.
# Usage: bash 02a_snmf_prepare.sh [INPUT.vcf.gz] [CHROMOSOME_MAP] [NEW_OUTPUT_DIRECTORY]
# The chromosome map has two columns without a header: original name and number.

set -euo pipefail

if [[ $# -gt 3 ]]; then
    echo "Usage: $0 [INPUT.vcf.gz] [CHROMOSOME_MAP] [NEW_OUTPUT_DIRECTORY]" >&2
    exit 1
fi
INPUT=${1:-Acer_main.NEU.PANEL_n46.max10miss.MAF05.vcf.gz}
CHROMOSOME_MAP=${2:-main.15.chromo_rename.txt}
OUT_DIR=${3:-snmf_input}

for program in bcftools tabix vcftools; do
    command -v "$program" >/dev/null || { echo "Missing program: $program" >&2; exit 1; }
done
for input_file in "$INPUT" "$CHROMOSOME_MAP"; do
    [[ -s "$input_file" ]] || { echo "Missing or empty input: $input_file" >&2; exit 1; }
done
mkdir -- "$OUT_DIR"

# 1. Rename chromosomes using the supplied mapping table.
RENAMED="$OUT_DIR/Acer_main.NEU.PANEL_n46_renm.chr.for.snmf_miss.MAF.filt.vcf.gz"
bcftools annotate --rename-chrs "$CHROMOSOME_MAP" -Oz -o "$RENAMED" "$INPUT"
tabix -p vcf "$RENAMED"

# 2. Convert the renamed VCF to PED/MAP format.
PREFIX="$OUT_DIR/Acer_main.NEU.PANEL_n46_plink_miss.MAF.filt"
vcftools --gzvcf "$RENAMED" --plink --out "$PREFIX"
echo "PED input for 02b_snmf.R: ${PREFIX}.ped"
