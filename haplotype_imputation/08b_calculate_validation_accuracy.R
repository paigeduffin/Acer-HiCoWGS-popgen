#!/usr/bin/env Rscript
# Author: Paige Duffin
# Calculate genotype concordance and mean sample-specific Pearson correlation
# for raw, GP-recalibrated, and newly imputed LoCo genotypes.
# HiCo truth calls are required to have depth >=10.
# Usage: Rscript 08b_calculate_validation_accuracy.R [EXTRACTION_DIRECTORY] [NEW_OUTPUT_DIRECTORY]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "validation_genotypes"
output_dir <- if (length(args) >= 2) args[[2]] else "validation_accuracy"
dir.create(output_dir, showWarnings = FALSE)

parse_gt <- function(x) {
  case_when(
    x %in% c("0/0", "0|0") ~ 0,
    x %in% c("0/1", "1/0", "0|1", "1|0") ~ 1,
    x %in% c("1/1", "1|1") ~ 2,
    TRUE ~ NA_real_
  )
}

standardize_hico_names <- function(x) {
  recode(
    x,
    "CRF_Acer-059" = "CRF_Acer59",
    "Mote_AC75" = "MML_ML75",
    "Mote_AC76" = "MML_ML76",
    "Mote_AC80" = "MML_ML80",
    "DRTO_114" = "DRTO_114"
  )
}

standardize_loco_names <- function(x) {
  sub("^ID_", "", x)
}

read_stage <- function(stage) {
  hico <- read_tsv(
    file.path(input_dir, paste0(stage, "_HiCo_genotypes.tsv")),
    show_col_types = FALSE
  ) %>%
    transmute(
      CHROM, POS,
      Genet = standardize_hico_names(SAMPLE),
      GT.Hi = parse_gt(GT),
      DP = as.numeric(DP)
    )

  loco <- read_tsv(
    file.path(input_dir, paste0(stage, "_LoCo_genotypes.tsv")),
    show_col_types = FALSE
  )

  if (!"DR2" %in% names(loco)) {
    loco$DR2 <- NA_real_
  }

  loco <- loco %>%
    transmute(
      CHROM, POS,
      Genet = standardize_loco_names(SAMPLE),
      GT.Lo = parse_gt(GT),
      DR2 = as.numeric(DR2)
    )

  inner_join(loco, hico, by = c("CHROM", "POS", "Genet"))
}

summarize_samples <- function(data, stage, threshold = NA_real_) {
  data %>%
    filter(DP >= 10, !is.na(GT.Lo), !is.na(GT.Hi)) %>%
    group_by(Genet) %>%
    summarize(
      Stage = stage,
      DR2_threshold = threshold,
      Nsites = n(),
      HomRefSites = sum(GT.Hi == 0),
      HomRefConcordant = sum(GT.Hi == GT.Lo & GT.Hi == 0),
      HetSites = sum(GT.Hi == 1),
      HetConcordant = sum(GT.Hi == GT.Lo & GT.Hi == 1),
      HomAltSites = sum(GT.Hi == 2),
      HomAltConcordant = sum(GT.Hi == GT.Lo & GT.Hi == 2),
      ConcordanceHomRef = HomRefConcordant / HomRefSites,
      ConcordanceHet = HetConcordant / HetSites,
      ConcordanceHomAlt = HomAltConcordant / HomAltSites,
      ConcordanceOverall = mean(GT.Hi == GT.Lo),
      Pearson_r = cor(GT.Lo, GT.Hi, method = "pearson", use = "complete.obs"),
      .groups = "drop"
    )
}

summarize_means <- function(per_sample) {
  per_sample %>%
    group_by(Stage, DR2_threshold) %>%
    summarize(
      MeanNsites = mean(Nsites),
      MeanConcordanceHomRef = mean(ConcordanceHomRef, na.rm = TRUE),
      MeanConcordanceHet = mean(ConcordanceHet, na.rm = TRUE),
      MeanConcordanceHomAlt = mean(ConcordanceHomAlt, na.rm = TRUE),
      MeanConcordanceOverall = mean(ConcordanceOverall, na.rm = TRUE),
      MeanPearson_r = mean(Pearson_r, na.rm = TRUE),
      .groups = "drop"
    )
}

raw_data <- read_stage("RAW")
recal_data <- read_stage("RECAL")
imputed_data <- read_stage("IMP")

raw_per_sample <- summarize_samples(raw_data, "RAW")
recal_per_sample <- summarize_samples(recal_data, "RECAL")

dr2_thresholds <- c(0, seq(0.10, 0.90, by = 0.10), 0.95, 0.99)
imputed_per_sample <- bind_rows(lapply(dr2_thresholds, function(threshold) {
  summarize_samples(
    filter(imputed_data, DR2 >= threshold),
    "IMPUTED",
    threshold
  )
}))

all_per_sample <- bind_rows(raw_per_sample, recal_per_sample, imputed_per_sample)
all_means <- summarize_means(all_per_sample)

write_tsv(
  all_per_sample,
  file.path(output_dir, "validation_accuracy_per_sample.tsv")
)
write_tsv(
  all_means,
  file.path(output_dir, "validation_accuracy_means.tsv")
)
write_tsv(
  filter(all_means, Stage == "IMPUTED"),
  file.path(output_dir, "imputed_accuracy_by_DR2_threshold.tsv")
)

plot_data <- filter(all_means, Stage == "IMPUTED") %>%
  select(
    DR2_threshold,
    MeanConcordanceHomRef,
    MeanConcordanceHet,
    MeanConcordanceHomAlt,
    MeanPearson_r
  ) %>%
  pivot_longer(
    -DR2_threshold,
    names_to = "Metric",
    values_to = "Value"
  )

accuracy_plot <- ggplot(
  plot_data,
  aes(x = DR2_threshold, y = Value, color = Metric)
) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0.99, color = "red") +
  scale_y_continuous(limits = c(0.5, 1.0)) +
  labs(
    x = "Beagle DR2 threshold",
    y = "Accuracy",
    color = NULL
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "imputed_accuracy_by_DR2_threshold.pdf"),
  accuracy_plot,
  width = 8,
  height = 5
)
