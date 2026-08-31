#!/usr/bin/env bash
#SBATCH --job-name=Acer_ADMIXTURE
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --output=Acer_ADMIXTURE.%j.out
#SBATCH --error=Acer_ADMIXTURE.%j.err

# Prepare the 15-chromosome neutral panel, run ADMIXTURE, and collect outputs.
# Software: ADMIXTURE 1.3.0, bcftools 1.19, HTSlib 1.19.1,
#           VCFtools 0.1.16, PLINK2 2.00a4.3.
# Load the software before running/submitting; adjust Slurm resources as needed.
# Usage: bash 04_admixture.sh [INPUT.vcf.gz] [CHROMOSOME_MAP] [NEW_OUTPUT_DIRECTORY]
# Or: sbatch 04_admixture.sh [INPUT.vcf.gz] [CHROMOSOME_MAP] [NEW_OUTPUT_DIRECTORY]
# If ADMIXTURE is not on PATH, set ADMIXTURE_BIN to the executable's full path.

set -euo pipefail

if [[ $# -gt 3 ]]; then
    echo "Usage: $0 [INPUT.vcf.gz] [CHROMOSOME_MAP] [NEW_OUTPUT_DIRECTORY]" >&2
    exit 1
fi
INPUT=${1:-Acer_main.NEU.PANEL_n46.max10miss.MAF05.vcf.gz}
CHROMOSOME_MAP=${2:-main.15.chromo_rename.txt}
OUT_DIR=${3:-admixture_results}
ADMIXTURE_BIN=${ADMIXTURE_BIN:-admixture}

for program in vcftools bcftools tabix plink2 "$ADMIXTURE_BIN"; do
    command -v "$program" >/dev/null || { echo "Missing program: $program" >&2; exit 1; }
done
ADMIXTURE_BIN=$(command -v "$ADMIXTURE_BIN")
[[ "$ADMIXTURE_BIN" = /* ]] || ADMIXTURE_BIN="$PWD/$ADMIXTURE_BIN"
for input_file in "$INPUT" "$CHROMOSOME_MAP"; do
    [[ -s "$input_file" ]] || { echo "Missing or empty input: $input_file" >&2; exit 1; }
done
mkdir -- "$OUT_DIR"
OUT_DIR=$(cd "$OUT_DIR" && pwd)
BED_PREFIX="Acer_main.NEU.PANEL_n46_miss.MAF.filt_15.chromo"
PREFIX="$OUT_DIR/$BED_PREFIX"

# 1. Retain only the 15 main chromosomes.
vcftools --gzvcf "$INPUT" \
    --chr OZ035966.1 --chr OZ035967.1 --chr OZ035968.1 \
    --chr OZ035969.1 --chr OZ035970.1 --chr OZ035971.1 \
    --chr OZ035972.1 --chr OZ035973.1 --chr OZ035974.1 \
    --chr OZ035975.1 --chr OZ035976.1 --chr OZ035977.1 \
    --chr OZ035978.1 --chr OZ035979.1 --chr OZ035980.1 \
    --recode --recode-INFO-all --out "$PREFIX"

# 2. Rename chromosomes and convert to binary PLINK format.
bcftools annotate --rename-chrs "$CHROMOSOME_MAP" \
    -Oz -o "${PREFIX}.vcf.gz" "${PREFIX}.recode.vcf"
tabix -p vcf "${PREFIX}.vcf.gz"
plink2 --vcf "${PREFIX}.vcf.gz" --make-bed --allow-extra-chr --out "$PREFIX"

# 3. Run K = 1-10 with 20 replicates each, preserving the per-run seed formula.
# Runs are sequential; move each Q/P pair before starting the next replicate.
cd "$OUT_DIR"
for K in {1..10}; do
    for run in {1..20}; do
        run_dir="output_files/K${K}_run${run}"
        mkdir -p "$run_dir"
        seed=$((K * 100000 + run))
        log_file="$run_dir/run.log"
        "$ADMIXTURE_BIN" --cv=10 --seed="$seed" "${BED_PREFIX}.bed" "$K" \
            > "$log_file" 2>&1
        for extension in Q P; do
            result="${BED_PREFIX}.${K}.${extension}"
            [[ -s "$result" ]] || { echo "Missing ADMIXTURE output: $result" >&2; exit 1; }
            mv "$result" "$run_dir/Acer_neu.pan_n46_K${K}_run${run}.${extension}"
        done
        printf 'K=%s run=%s seed=%s\n' "$K" "$run" "$seed" >> "$log_file"
    done
done

# 4. Group all Q matrices by K for CLUMPAK.
for K in {1..10}; do
    k_dir="clumpak/K${K}_Acer_neu.pan_n46_admixture"
    mkdir -p "$k_dir"
    for run in {1..20}; do
        cp "output_files/K${K}_run${run}/Acer_neu.pan_n46_K${K}_run${run}.Q" \
            "$k_dir/Acer_neu.pan_n46_K${K}_run${run}_admixture.txt"
    done
done

# 5. Extract CV errors, log-likelihoods, and convergence lines with their sources.
awk '/CV error/ {print FILENAME ":" FNR ":" $0}' \
    output_files/K*_run*/run.log > CV_values.txt
awk '/^Loglikelihood:/ {print FILENAME ":" FNR ":" $0}' \
    output_files/K*_run*/run.log > Loglikelihood.txt
awk '/^Converged in/ {print FILENAME ":" FNR ":" $0}' \
    output_files/K*_run*/run.log > Converged_in.txt

echo "ADMIXTURE results, CLUMPAK files, and log summaries saved in: $OUT_DIR"

# After downloading the results directory, open a terminal inside its clumpak folder:
# zip -r ../admixture_for_CLUMPAK.zip K*_Acer_neu.pan_n46_admixture -x "*/.*"
