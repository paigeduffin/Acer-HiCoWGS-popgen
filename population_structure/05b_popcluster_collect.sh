#!/usr/bin/env bash
# Collect PopCluster Q matrices and K-file summaries after all five jobs finish.
# Usage: bash 05b_popcluster_collect.sh [INPUT_ROOT] [NEW_OUTPUT_DIRECTORY]
# DLK2 and FST/FIS are extracted when reported; this script does not calculate them.

set -euo pipefail

if [[ $# -gt 2 ]]; then
    echo "Usage: $0 [INPUT_ROOT] [NEW_OUTPUT_DIRECTORY]" >&2
    exit 1
fi
INPUT_ROOT=${1:-popcluster_results}
OUT_DIR=${2:-popcluster_collected}

# Check for every expected K/replicate output before collecting.
for K in {1..10}; do
    for replicate in {1..5}; do
        prefix="$INPUT_ROOT/PopClust_rand.seed_run${replicate}/long_K${K}_r${replicate}"
        for input_file in "${prefix}.K" "${prefix}.K.${K}.R.1_Q"; do
            [[ -s "$input_file" ]] || { echo "Missing or empty output: $input_file" >&2; exit 1; }
        done
    done
done
mkdir -- "$OUT_DIR"
mkdir "$OUT_DIR/PopCluster_CLUMPAK_ready" "$OUT_DIR/K_files_full"
summary="$OUT_DIR/summary_output.txt"
dlk2="$OUT_DIR/bestK_DLK2_output.txt"
fstfis="$OUT_DIR/bestK_FSTFIS_output.txt"
: > "$summary"
: > "$dlk2"
: > "$fstfis"

for K in {1..10}; do
    destination="$OUT_DIR/PopCluster_CLUMPAK_ready/K${K}_Acer.neu.pan_PopCluster_admixture"
    mkdir "$destination"
    for replicate in {1..5}; do
        prefix="$INPUT_ROOT/PopClust_rand.seed_run${replicate}/long_K${K}_r${replicate}"

        # Remove the individual-ID column and write the K ancestry columns as tabs.
        awk -F',' -v k="$K" 'BEGIN {OFS="\t"}
            {
                sub(/\r$/, "")
                if (NF != k+1) {
                    print "Unexpected column count in " FILENAME > "/dev/stderr"
                    exit 1
                }
                for (i=2; i<=NF; i++) printf "%s%s", $i, (i==NF ? ORS : OFS)
            }
            END {
                if (NR != 46) {
                    print "Expected 46 individuals in " FILENAME > "/dev/stderr"
                    exit 1
                }
            }
        ' "${prefix}.K.${K}.R.1_Q" > "$destination/K${K}_r${replicate}.Q"

        cp "${prefix}.K" "$OUT_DIR/K_files_full/long_K${K}_r${replicate}.K"
        awk '$1 ~ /^[0-9]+$/ && $2 ~ /^"long_K/ {print; exit}' \
            "${prefix}.K" >> "$summary"
        awk '$1 == "DLK2" {print; exit}' "${prefix}.K" >> "$dlk2"
        awk '$1 == "FST/FIS" {print; exit}' "${prefix}.K" >> "$fstfis"
    done
done

[[ -s "$dlk2" ]] || echo "No DLK2 lines were present in the K files."
[[ -s "$fstfis" ]] || echo "No FST/FIS lines were present in the K files."
echo "Collected PopCluster outputs: $OUT_DIR"

# After downloading the collected output directory, open a terminal inside it:
# zip -r PopCluster_CLUMPAK_ready.zip PopCluster_CLUMPAK_ready/K*_Acer.neu.pan_PopCluster_admixture -x "*/.*"
