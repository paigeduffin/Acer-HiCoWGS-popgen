#!/usr/bin/env bash
#SBATCH --job-name=Acer_PopCluster
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=24:00:00
#SBATCH --array=1-5
#SBATCH --output=Acer_PopCluster.%A_%a.out
#SBATCH --error=Acer_PopCluster.%A_%a.err

# Run five independent replicates, each evaluating K = 1-10.
# Software: PopCluster (PopClusterLnx); set POPCLUSTER_BIN to its full path
# if it is not on PATH. Adjust Slurm resources for your cluster as needed.
# Input: full neutral-panel PLINK files, prepared by 01a_pca_prepare.sh.
# Usage: sbatch 05a_popcluster.sh [BED_PREFIX] [PROJECT_FILE] [OUTPUT_ROOT]
# Without Slurm, run each replicate using REPLICATE=1 bash 05a_popcluster.sh
# (repeat for REPLICATE=2 through 5). BED_PREFIX excludes .bed/.bim/.fam.

set -euo pipefail

if [[ $# -gt 3 ]]; then
    echo "Usage: $0 [BED_PREFIX] [PROJECT_FILE] [OUTPUT_ROOT]" >&2
    exit 1
fi
BED_PREFIX=${1:-pca_input/Acer_main.NEU.PANEL_n46.max10miss.MAF05}
PROJECT_FILE=${2:-Acer_popgen_n46_178k.PcPjt}
OUTPUT_ROOT=${3:-popcluster_results}
POPCLUSTER_BIN=${POPCLUSTER_BIN:-PopClusterLnx}
replicate=${SLURM_ARRAY_TASK_ID:-${REPLICATE:-1}}
[[ "$replicate" =~ ^[1-5]$ ]] || { echo "Replicate must be 1-5." >&2; exit 1; }

command -v "$POPCLUSTER_BIN" >/dev/null || { echo "Cannot find PopClusterLnx." >&2; exit 1; }
POPCLUSTER_BIN=$(command -v "$POPCLUSTER_BIN")
[[ "$POPCLUSTER_BIN" = /* ]] || POPCLUSTER_BIN="$PWD/$POPCLUSTER_BIN"
for input_file in "$PROJECT_FILE" "${BED_PREFIX}.bed" "${BED_PREFIX}.bim" "${BED_PREFIX}.fam"; do
    [[ -s "$input_file" ]] || { echo "Missing or empty input: $input_file" >&2; exit 1; }
done

# The supplied project specifies 46 individuals and 177,516 loci.
[[ $(awk 'END {print NR}' "${BED_PREFIX}.fam") -eq 46 ]] || { echo "Expected 46 individuals." >&2; exit 1; }
[[ $(awk 'END {print NR}' "${BED_PREFIX}.bim") -eq 177516 ]] || { echo "Expected 177516 loci." >&2; exit 1; }
BED_PREFIX="$(cd "$(dirname "$BED_PREFIX")" && pwd)/$(basename "$BED_PREFIX")"
project_bed=$(awk 'NR==6 {sub(/\r$/, ""); print; exit}' "$PROJECT_FILE")
[[ "$project_bed" == "Acer_main.NEU.PANEL_n46.max10miss.MAF05.bed" ]] || {
    echo "Unexpected BED filename in the project file." >&2; exit 1;
}

# Each array task has its own directory; existing replicate results are protected.
mkdir -p "$OUTPUT_ROOT"
OUTPUT_ROOT=$(cd "$OUTPUT_ROOT" && pwd)
run_dir="$OUTPUT_ROOT/PopClust_rand.seed_run${replicate}"
mkdir -- "$run_dir"
cp "$PROJECT_FILE" "$run_dir/Acer_popgen_n46_178k.PcPjt"
for extension in bed bim fam; do
    ln -s "${BED_PREFIX}.${extension}" "$run_dir/${project_bed%.bed}.${extension}"
done
cd "$run_dir"

printf 'K\treplicate\tseed\n' > seeds.tsv
for K in {1..10}; do
    seed=$(( (RANDOM << 16) ^ (RANDOM << 1) ^ replicate ))
    printf '%s\t%s\t%s\n' "$K" "$replicate" "$seed" >> seeds.tsv
    "$POPCLUSTER_BIN" INP:Acer_popgen_n46_178k.PcPjt \
        "K_1:$K" "K_2:$K" RUN:1 SKI:0 "RAN:$seed" \
        "OUT:long_K${K}_r${replicate}" >> popcluster.log 2>&1
done
echo "PopCluster replicate $replicate finished: $run_dir"
