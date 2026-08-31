#!/usr/bin/env bash
# Author: Paige Duffin
# Thin and LD-prune the phased 46-sample panel, inspect structure with PCA and
# ADMIXTURE, remove two admixed samples, and generate the 44-sample PCA input
# used for RUTH filtering.
# Software: VCFtools 0.1.16, bcftools 1.19, HTSlib 1.19.1,
# PLINK2 2.00a4.3, and ADMIXTURE 1.3.0.
# Set ADMIXTURE_BIN to the executable name or its absolute path if needed.
# Usage: bash 03_between_phasing_pca_admixture.sh [PHASED_VCF.gz] [CHROMOSOME_RENAME_FILE] [NEW_OUTPUT_DIRECTORY]

set -euo pipefail

PHASED_VCF=${1:-whatshap_phasing_output/merged_step1phase_n46.vcf.gz}
CHROM_RENAME=${2:-../population_structure/main.15.chromo_rename.txt}
OUT_DIR=${3:-between_phasing_steps}
ADMIXTURE_BIN=${ADMIXTURE_BIN:-admixture}

if [[ "$OUT_DIR" != /* ]]; then
    OUT_DIR="$(pwd)/$OUT_DIR"
fi
if [[ "$ADMIXTURE_BIN" == */* && "$ADMIXTURE_BIN" != /* ]]; then
    ADMIXTURE_BIN="$(pwd)/$ADMIXTURE_BIN"
fi

INTERMEDIATE="$OUT_DIR/intermediate_vcfs"
ADMIX_DIR="$OUT_DIR/admixture"
ADMIX_RUNS="$ADMIX_DIR/output_files"
CLUMPAK_DIR="$ADMIX_DIR/admix4phasing_clumpak_ready"
N44_DIR="$OUT_DIR/after_remove_admixed"

mkdir -- "$OUT_DIR"
mkdir -p -- "$INTERMEDIATE" "$ADMIX_DIR" "$ADMIX_RUNS" "$CLUMPAK_DIR" "$N44_DIR"

# Retain biallelic SNPs with MAC >= 3 and thin to one site per 20 kb.
vcftools \
    --gzvcf "$PHASED_VCF" \
    --mac 3 \
    --max-alleles 2 \
    --min-alleles 2 \
    --remove-indels \
    --thin 20000 \
    --recode \
    --out "$INTERMEDIATE/HWE_thinned_PD_n46"

# Assign chromosome-position-allele IDs to the retained variants.
bcftools annotate \
    --set-id +'%CHROM_%POS_%REF_%FIRST_ALT' \
    "$INTERMEDIATE/HWE_thinned_PD_n46.recode.vcf" \
    -Ov \
    -o "$INTERMEDIATE/HWE_PD_n46.vcf"

# Build the PLINK files, identify LD-pruned sites, and run the 46-sample PCA.
plink2 \
    --vcf "$INTERMEDIATE/HWE_PD_n46.vcf" \
    --allow-extra-chr \
    --make-bed \
    --out "$INTERMEDIATE/HWE_PD_n46_clean"

plink2 \
    --bfile "$INTERMEDIATE/HWE_PD_n46_clean" \
    --bad-ld \
    --allow-extra-chr \
    --indep-pairwise 500kb 0.2 \
    --out "$INTERMEDIATE/HWE_pca_PD_n46_clean"

plink2 \
    --bfile "$INTERMEDIATE/HWE_PD_n46_clean" \
    --allow-extra-chr \
    --freq \
    --out "$INTERMEDIATE/HWE_PD_freq_n46_clean"

plink2 \
    --bfile "$INTERMEDIATE/HWE_PD_n46_clean" \
    --allow-extra-chr \
    --extract "$INTERMEDIATE/HWE_pca_PD_n46_clean.prune.in" \
    --read-freq "$INTERMEDIATE/HWE_PD_freq_n46_clean.afreq" \
    --pca \
    --out "$INTERMEDIATE/HWE_pca_PD_n46_final"

vcftools \
    --vcf "$INTERMEDIATE/HWE_PD_n46.vcf" \
    --snps "$INTERMEDIATE/HWE_pca_PD_n46_clean.prune.in" \
    --recode \
    --out "$INTERMEDIATE/HWE_pca_PD_n46_vcftools"

# Restrict the LD-pruned panel to the 15 main chromosomes for ADMIXTURE.
CHR_ARGS=()
for CHROM_NUMBER in {966..980}; do
    CHR_ARGS+=(--chr "OZ035${CHROM_NUMBER}.1")
done

vcftools \
    --vcf "$INTERMEDIATE/HWE_pca_PD_n46_vcftools.recode.vcf" \
    "${CHR_ARGS[@]}" \
    --recode \
    --recode-INFO-all \
    --out "$ADMIX_DIR/N.Loca_filtering_n46_vcf1"

bcftools annotate \
    --rename-chrs "$CHROM_RENAME" \
    "$ADMIX_DIR/N.Loca_filtering_n46_vcf1.recode.vcf" \
    -Oz \
    -o "$ADMIX_DIR/N.Loca_filtering_n46_vcf2.vcf.gz"

tabix -p vcf "$ADMIX_DIR/N.Loca_filtering_n46_vcf2.vcf.gz"

plink2 \
    --vcf "$ADMIX_DIR/N.Loca_filtering_n46_vcf2.vcf.gz" \
    --make-bed \
    --allow-extra-chr \
    --out "$ADMIX_DIR/N.Loca_filtering_n46_vcf2"

# Run 20 reproducibly seeded ADMIXTURE replicates for each K from 1 through 10.
ADMIX_BASENAME="N.Loca_filtering_n46_vcf2"
for K in {1..10}; do
    K_CLUMPAK_DIR="$CLUMPAK_DIR/K${K}_Acer_admix4phasing"
    mkdir -p -- "$K_CLUMPAK_DIR"

    for RUN in {1..20}; do
        RUN_DIR="$ADMIX_RUNS/K${K}_run${RUN}_Acer_admix4phasing"
        LOG_FILE="$RUN_DIR/CVflag_admix4phasing_K${K}_run${RUN}.log"
        SEED=$((K * 100000 + RUN))
        mkdir -- "$RUN_DIR"

        (
            cd "$ADMIX_DIR"
            "$ADMIXTURE_BIN" \
                --cv=10 \
                --seed="$SEED" \
                "${ADMIX_BASENAME}.bed" \
                "$K" \
                > "$LOG_FILE" 2>&1
        )

        mv -- \
            "$ADMIX_DIR/${ADMIX_BASENAME}.${K}.Q" \
            "$RUN_DIR/${ADMIX_BASENAME}.${K}.Q"
        mv -- \
            "$ADMIX_DIR/${ADMIX_BASENAME}.${K}.P" \
            "$RUN_DIR/${ADMIX_BASENAME}.${K}.P"

        cp -- \
            "$RUN_DIR/${ADMIX_BASENAME}.${K}.Q" \
            "$K_CLUMPAK_DIR/Acer_admix4phasing_K${K}_run${RUN}_admixture.txt"
    done
done

# Collect ADMIXTURE diagnostics across all runs.
grep -R --include='*.log' -nH 'CV error' "$ADMIX_RUNS" \
    > "$ADMIX_DIR/CV_values_Acer_admix4phasing.txt" || true
grep -R --include='*.log' -nH '^Loglikelihood:' "$ADMIX_RUNS" \
    > "$ADMIX_DIR/Loglikelihood_Acer_admix4phasing.txt" || true
grep -R --include='*.log' -nH '^Converged in' "$ADMIX_RUNS" \
    > "$ADMIX_DIR/Converged_in_Acer_admix4phasing.txt" || true

# Remove the two admixed samples identified above and generate the 44-sample PCA.
bcftools query -l "$INTERMEDIATE/HWE_PD_n46.vcf" \
    | awk '$0 != "CRF_Acer-099" && $0 != "SRR24007608"' \
    > "$N44_DIR/n44.keep.txt"

vcftools \
    --vcf "$INTERMEDIATE/HWE_PD_n46.vcf" \
    --keep "$N44_DIR/n44.keep.txt" \
    --recode \
    --out "$N44_DIR/rm.admixed_n44"

plink2 \
    --vcf "$N44_DIR/rm.admixed_n44.recode.vcf" \
    --allow-extra-chr \
    --make-bed \
    --out "$N44_DIR/HWE_PD_n44_PCA"

plink2 \
    --bfile "$N44_DIR/HWE_PD_n44_PCA" \
    --allow-extra-chr \
    --freq \
    --out "$N44_DIR/HWE_PD_n44_freq_PCA"

plink2 \
    --bfile "$N44_DIR/HWE_PD_n44_PCA" \
    --allow-extra-chr \
    --extract "$INTERMEDIATE/HWE_pca_PD_n46_clean.prune.in" \
    --read-freq "$N44_DIR/HWE_PD_n44_freq_PCA.afreq" \
    --pca \
    --out "$N44_DIR/HWE_PD_n44_PCA_final"

# Optional: package the CLUMPAK input directories after downloading them.
# zip -r admix4phasing_clumpak_ready.zip admix4phasing_clumpak_ready/K*_Acer_admix4phasing -x '*/.*'
