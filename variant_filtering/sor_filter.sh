#!/usr/bin/env bash
# Calculate site-level SOR from INFO/DP4, annotate the VCF, and retain SOR<4.
# Based on Paige Duffin's SOR_calc_10.16.25.sh and supplied filter command.
#
# Usage: bash sor_filter.sh INPUT.vcf[.gz] NEW_OUTPUT_DIRECTORY
# Requirement: bcftools 1.19; activate the environment before running.
# INPUT must retain INFO/DP4 (ref-forward, ref-reverse, alt-forward,
# alt-reverse counts). Existing SOR values are recalculated from DP4.

set -euo pipefail
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 INPUT.vcf[.gz] NEW_OUTPUT_DIRECTORY" >&2
    exit 1
fi
INPUT=$1
OUT_DIR=$2
command -v bcftools >/dev/null || { echo "Missing program: bcftools" >&2; exit 1; }
[[ -r "$INPUT" ]] || { echo "Cannot read input: $INPUT" >&2; exit 1; }
bcftools view -h "$INPUT" | awk '
    /^##INFO=<ID=DP4,/ {found=1}
    END {exit !found}
' || { echo "Input requires an INFO/DP4 header and strand-count data." >&2; exit 1; }
mkdir -- "$OUT_DIR"

# Process each VCF record directly, avoiding coordinate-only joins that
# can misassign annotations when multiple records share CHROM and POS.
# Install the header before filtering. Missing DP4 produces SOR=.; malformed
# DP4 stops the workflow instead of silently interpreting invalid counts.
bcftools view -Ov "$INPUT" | awk -v table="$OUT_DIR/sor_values.tsv" '
BEGIN {
    FS=OFS="\t"
    print "CHROM", "POS", "REF", "ALT", "SOR" > table
}
/^##INFO=<ID=SOR,/ {next}
/^#CHROM/ {
    print "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Strand Odds Ratio calculated from INFO/DP4 with a pseudocount of one; smaller is better\">"
    print; next
}
/^#/ {print; next}
{
    n=split($8, info, ";"); rebuilt=""; raw="."
    for (i=1; i<=n; i++) {
        if (info[i] ~ /^DP4=/) raw=substr(info[i], 5)
        if (info[i] !~ /^SOR=/ && info[i] != "." && info[i] != "")
            rebuilt=rebuilt (rebuilt=="" ? "" : ";") info[i]
    }
    value="."
    if (raw != ".") {
        k=split(raw, dp4, ",")
        valid=(k==4)
        for (i=1; i<=k; i++) if (dp4[i] !~ /^[0-9]+$/) valid=0
        if (!valid) {
            print "Invalid DP4 at " $1 ":" $2 > "/dev/stderr"
            exit 2
        }
        ref_fw=dp4[1]+1; ref_rv=dp4[2]+1
        alt_fw=dp4[3]+1; alt_rv=dp4[4]+1
        term1=(ref_fw*alt_rv)/(alt_fw*ref_rv)
        ln_ref=log((ref_fw<ref_rv ? ref_fw : ref_rv)/(ref_fw>ref_rv ? ref_fw : ref_rv))
        ln_alt=log((alt_fw<alt_rv ? alt_fw : alt_rv)/(alt_fw>alt_rv ? alt_fw : alt_rv))
        sor=log(term1+1/term1)+ln_ref-ln_alt
        # Match the original awk table precision before applying SOR<4.
        value=sprintf("%.6g", sor)
    }
    $8=rebuilt (rebuilt=="" ? "" : ";") "SOR=" value
    print
    print $1, $2, $4, $5, value > table
}' | bcftools view -Oz -o "$OUT_DIR/sor_annotated.vcf.gz"
bcftools index -t "$OUT_DIR/sor_annotated.vcf.gz"

# Strict inequality, exactly as in the supplied filtering command.
bcftools view -i 'INFO/SOR<4' -Oz \
    -o "$OUT_DIR/sor_filtered.vcf.gz" "$OUT_DIR/sor_annotated.vcf.gz"
bcftools index -t "$OUT_DIR/sor_filtered.vcf.gz"

printf 'dataset\tvariant_records\n' > "$OUT_DIR/record_counts.tsv"
for item in annotated filtered; do
    count=$(bcftools view -H "$OUT_DIR/sor_${item}.vcf.gz" | awk 'END {print NR+0}')
    printf '%s\t%s\n' "$item" "$count" >> "$OUT_DIR/record_counts.tsv"
done
