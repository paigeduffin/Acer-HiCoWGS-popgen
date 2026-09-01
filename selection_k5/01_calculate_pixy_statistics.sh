#!/usr/bin/env bash
#SBATCH --job-name=Acer_K5_pixy
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=Acer_K5_pixy.%j.out
#SBATCH --error=Acer_K5_pixy.%j.err

# Calculate nucleotide diversity (pi) and pairwise FST for the five post hoc
# population groupings in non-overlapping 10 kb windows. Only windows with at
# least 5 kb of callable sites are retained.
#
# Required inputs:
# - Combined variant/invariant VCF used for the final Pixy analysis.
# - Two-column callable-site file used with Pixy's --sites_file option.
# - Two-column sample/population file. The repository's
#   ../variant_filtering/K5_morePAN_n46_2.col.txt contains the assignments used.
#
# Software: Pixy 2.0.0 and HTSlib. Activate the Pixy environment before
# submitting this script.
#
# Usage:
#   sbatch 01_calculate_pixy_statistics.sh [VCF] [SITES] [POPULATIONS] [OUT_DIR]

set -euo pipefail

if [[ $# -gt 4 ]]; then
    echo "Usage: $0 [VCF] [SITES] [POPULATIONS] [OUT_DIR]" >&2
    exit 1
fi

VCF=${1:-snps.invariant.morePAN.minDP10.maxDP50.qual.SOR.miss10.nosort.vcf.gz}
SITES=${2:-sites_rmHWE.txt}
POPULATIONS=${3:-../variant_filtering/K5_morePAN_n46_2.col.txt}
OUT_DIR=${4:-pixy_results}

command -v pixy >/dev/null || {
    echo "Pixy is not available. Activate the Pixy 2.0.0 environment first." >&2
    exit 1
}

for input_file in "$VCF" "$SITES" "$POPULATIONS"; do
    [[ -s "$input_file" ]] || {
        echo "Missing or empty input: $input_file" >&2
        exit 1
    }
done

mkdir -- "$OUT_DIR"
WINDOWS="$OUT_DIR/10kb_windows_min5kb_present.bed"

# Build a BED file containing 10 kb windows represented by at least 5,000
# callable positions in the sites file.
awk 'BEGIN {OFS="\t"}
     !/^#/ && NF >= 2 && $2 ~ /^[0-9]+$/ {
         chromosome=$1
         start=int(($2-1)/10000)*10000
         end=start+10000
         count[chromosome SUBSEP start SUBSEP end]++
     }
     END {
         for (window in count) {
             if (count[window] >= 5000) {
                 split(window, fields, SUBSEP)
                 print fields[1], fields[2], fields[3]
             }
         }
     }' "$SITES" | LC_ALL=C sort -k1,1 -k2,2n > "$WINDOWS"

[[ -s "$WINDOWS" ]] || {
    echo "No windows passed the 5 kb callable-site requirement." >&2
    exit 1
}

pixy --stats pi fst \
    --vcf "$VCF" \
    --populations "$POPULATIONS" \
    --include_multiallelic_snps \
    --sites_file "$SITES" \
    --bed_file "$WINDOWS" \
    --n_cores 4 \
    --output_prefix 10kb_min5kb \
    --output_folder "$OUT_DIR"

echo "Pixy outputs saved in: $OUT_DIR"
