#!/usr/bin/env bash
# Convert the neutral SNP panel to PLINK format for 01b_pca.R.
# Software: PLINK2 2.00a4.3.
# USC HPC modules: ml gcc/13.3.0 plink2/2.00a4.3
# Usage: bash 01a_pca_prepare.sh [INPUT.vcf.gz] [NEW_OUTPUT_DIRECTORY]

set -euo pipefail

if [[ $# -gt 2 ]]; then
    echo "Usage: $0 [INPUT.vcf.gz] [NEW_OUTPUT_DIRECTORY]" >&2
    exit 1
fi
INPUT=${1:-Acer_main.NEU.PANEL_n46.max10miss.MAF05.vcf.gz}
OUT_DIR=${2:-pca_input}
command -v plink2 >/dev/null || { echo "Load PLINK2 before running this script." >&2; exit 1; }
[[ -r "$INPUT" ]] || { echo "Cannot read input: $INPUT" >&2; exit 1; }
mkdir -- "$OUT_DIR"

PREFIX="$OUT_DIR/Acer_main.NEU.PANEL_n46.max10miss.MAF05"
plink2 --vcf "$INPUT" --make-bed --allow-extra-chr --double-id --out "$PREFIX"

echo "PLINK files prepared. BED_PREFIX for 01b_pca.R: $PREFIX"
