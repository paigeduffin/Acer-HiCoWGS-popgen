#!/usr/bin/env bash
# RUTH/HWE filtering for the expanded manuscript dataset.
# Reconstructed from Paige Duffin's August 2025 workflow, subsequently
# rerun on the expanded dataset. This is not the haplotype-panel filter.
#
# Usage: bash ruth_hwe_filter.sh INPUT.vcf[.gz] NEW_OUTPUT_DIRECTORY
# Requirements: bcftools 1.19 (+fill-tags), PLINK2, RUTH (study: v3.29.4).
# Activate the appropriate software environment before running/submitting.
# If needed: export RUTH_BIN=/absolute/path/to/ruth

set -euo pipefail
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 INPUT.vcf[.gz] NEW_OUTPUT_DIRECTORY" >&2
    exit 1
fi
INPUT=$1
OUT_DIR=$2
RUTH_BIN=${RUTH_BIN:-ruth}
for program in bcftools plink2 "$RUTH_BIN"; do
    command -v "$program" >/dev/null || { echo "Missing program: $program" >&2; exit 1; }
done
[[ -r "$INPUT" ]] || { echo "Cannot read input: $INPUT" >&2; exit 1; }
bcftools view -h "$INPUT" | awk '
    /^##FORMAT=<ID=PL,/ {found=1}
    END {exit !found}
' || { echo "RUTH input must contain FORMAT/PL genotype likelihoods." >&2; exit 1; }
# Require a new directory to avoid overwriting a previous run.
mkdir -- "$OUT_DIR"

# 1. Use biallelic SNPs to estimate PCs. The full SNP input is retained
# separately for RUTH; multiallelic sites are not removed from that input.
# No additional MAF, MAC, thinning, or LD pruning is applied here.
bcftools view -m2 -M2 -v snps -Ou "$INPUT" \
    | bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
        -Ov -o "$OUT_DIR/pca_input.vcf"

# --double-id explicitly preserves full VCF sample names in the IID field.
plink2 --vcf "$OUT_DIR/pca_input.vcf" --double-id \
    --allow-extra-chr --make-bed --out "$OUT_DIR/pca_input"
plink2 --bfile "$OUT_DIR/pca_input" --allow-extra-chr \
    --freq --out "$OUT_DIR/pca_frequencies"
plink2 --bfile "$OUT_DIR/pca_input" --allow-extra-chr \
    --read-freq "$OUT_DIR/pca_frequencies.afreq" \
    --pca --out "$OUT_DIR/pca"

# Select IID and PC columns by header name, whether or not FID is present.
awk 'BEGIN {OFS="\t"}
    NR==1 {
        for (i=1; i<=NF; i++) {
            if ($i=="IID" || $i=="#IID") iid=i
            if ($i ~ /^PC[0-9]+$/) pc[++n]=i
        }
        if (!iid || !n) {print "Unrecognized PLINK eigenvec header" > "/dev/stderr"; exit 1}
        printf "ID"; for (j=1; j<=n; j++) printf "\t%s", $pc[j]; print ""
        next
    }
    {printf "%s", $iid; for (j=1; j<=n; j++) printf "\t%s", $pc[j]; print ""}
' "$OUT_DIR/pca.eigenvec" > "$OUT_DIR/ruth.eigenvec"

# 2. Recalculate INFO tags in a temporary RUTH input; FORMAT/PL remains.
bcftools annotate --remove INFO -Ou "$INPUT" \
    | bcftools +fill-tags -Ov -o "$OUT_DIR/ruth_input.vcf" -- -t all
"$RUTH_BIN" --evec "$OUT_DIR/ruth.eigenvec" \
    --vcf "$OUT_DIR/ruth_input.vcf" --field PL --site-only \
    --out "$OUT_DIR/ruth_results.vcf"

# 3. Retain positions strictly between -5 and 5, excluding missing scores.
bcftools query -f '%CHROM\t%POS\t%INFO/HWE_SLP_I\t%INFO/MAF\n' \
    "$OUT_DIR/ruth_results.vcf" > "$OUT_DIR/ruth_statistics.tsv"
awk 'BEGIN {FS=OFS="\t"}
    $3 != "." && $3 > -5 && $3 < 5 {print $1, $2}
' "$OUT_DIR/ruth_statistics.tsv" > "$OUT_DIR/ruth_pass_positions.tsv"

if [[ ! -s "$OUT_DIR/ruth_pass_positions.tsv" ]]; then
    echo "No positions passed; inspect ruth_results.vcf before continuing." >&2
    exit 1
fi

# Apply the passing positions to the original, genotyped VCF. In contrast
# to the notebook's final INFO-removal step, retain the original INFO
# annotations (including DP4/SOR, if present) in this deliverable.
bcftools view -T "$OUT_DIR/ruth_pass_positions.tsv" \
    -Oz -o "$OUT_DIR/hwe_filtered.vcf.gz" "$INPUT"
bcftools index -t "$OUT_DIR/hwe_filtered.vcf.gz"

printf 'dataset\tvariant_records\n' > "$OUT_DIR/record_counts.tsv"
for item in input retained; do
    VCF=$INPUT
    [[ "$item" != retained ]] || VCF="$OUT_DIR/hwe_filtered.vcf.gz"
    count=$(bcftools view -H "$VCF" | awk 'END {print NR+0}')
    printf '%s\t%s\n' "$item" "$count" >> "$OUT_DIR/record_counts.tsv"
done
