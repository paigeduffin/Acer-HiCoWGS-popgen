# Author: Paige Duffin
# Merge K = 5 meadow-plot outlier windows into candidate peaks and summarize
# the number of outlier windows per peak.
#
# Peaks may span up to three consecutive retained non-outlier windows. Windows
# removed upstream for insufficient callable data are ignored, and intervening
# non-outlier windows are excluded from the adjusted peak-size score. This
# reproduces the rule used in the manuscript without manual counting in Excel.
#
# Usage:
#   Rscript 03_count_candidate_peaks.R [MEADOW_DATA.csv] [OUT_DIR]

library(dplyr)
library(readr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2) {
  stop("Usage: Rscript 03_count_candidate_peaks.R [MEADOW_DATA.csv] [OUT_DIR]")
}

input_file <- if (length(args) >= 1) {
  args[1]
} else {
  "meadow_plot_results/meadow_window_statistics.csv"
}
out_dir <- if (length(args) >= 2) args[2] else "candidate_peak_results"

if (!file.exists(input_file)) stop("Missing input: ", input_file)
if (dir.exists(out_dir)) stop("Output directory already exists: ", out_dir)
dir.create(out_dir, recursive = TRUE)

meadow_data <- read_csv(input_file, show_col_types = FALSE)
required_columns <- c(
  "comparison", "chromosome", "window_pos_1", "window_pos_2", "outlier"
)
if (!all(required_columns %in% names(meadow_data))) {
  stop(
    "Input is missing: ",
    paste(setdiff(required_columns, names(meadow_data)), collapse = ", ")
  )
}

# Give every retained window a sequential index within its chromosome and
# comparison. Outlier indices differing by no more than four belong to the same
# peak because at most three retained non-outlier windows separate them.
indexed_windows <- meadow_data %>%
  arrange(comparison, chromosome, window_pos_1) %>%
  group_by(comparison, chromosome) %>%
  mutate(retained_window_index = row_number()) %>%
  ungroup()

candidate_peaks <- indexed_windows %>%
  filter(outlier) %>%
  group_by(comparison, chromosome) %>%
  arrange(retained_window_index, .by_group = TRUE) %>%
  mutate(
    new_peak = is.na(lag(retained_window_index)) |
      retained_window_index - lag(retained_window_index) > 4,
    peak_number = cumsum(new_peak)
  ) %>%
  group_by(comparison, chromosome, peak_number) %>%
  summarise(
    first_window_start = min(window_pos_1),
    last_window_end = max(window_pos_2),
    n_outlier_windows = n(),
    n_intervening_retained_nonoutliers =
      max(retained_window_index) - min(retained_window_index) + 1 - n(),
    total_genomic_span_kb =
      (max(window_pos_2) - min(window_pos_1) + 1) / 1000,
    adjusted_outlier_span_kb = n_outlier_windows * 10,
    .groups = "drop"
  ) %>%
  arrange(desc(n_outlier_windows), comparison, chromosome, first_window_start)

write_csv(candidate_peaks, file.path(out_dir, "candidate_peaks.csv"))

outlier_window_counts <- indexed_windows %>%
  filter(outlier) %>%
  count(comparison, name = "n_outlier_windows")
write_csv(
  outlier_window_counts,
  file.path(out_dir, "outlier_window_counts_by_comparison.csv")
)

comparison_order <- c(
  "ARCU.v.BEME", "ARCU.v.PANM", "DRJA.v.ARCU", "DRJA.v.BEME",
  "DRJA.v.PANM", "FLDT.v.ARCU", "FLDT.v.BEME", "FLDT.v.DRJA",
  "FLDT.v.PANM", "PANM.v.BEME"
)

plot_data <- candidate_peaks %>%
  mutate(comparison = factor(comparison, levels = comparison_order))

peak_length_plot <- ggplot(
  plot_data,
  aes(comparison, n_outlier_windows, group = comparison)
) +
  geom_violin(width = 1, color = "#2D5564", fill = "#2D5564", alpha = 0.25) +
  geom_jitter(
    aes(size = n_outlier_windows),
    width = 0.25, shape = 21, color = "#2D5564", fill = "#2D5564",
    alpha = 0.85
  ) +
  scale_size(range = c(1.5, 2.5), guide = "none") +
  labs(
    x = "K = 5 population comparison",
    y = "Outlier 10 kb windows per candidate peak"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out_dir, "candidate_peak_length_distributions.png"),
  peak_length_plot,
  width = 16, height = 8, units = "cm", dpi = 600, bg = "transparent"
)
ggsave(
  file.path(out_dir, "candidate_peak_length_distributions.pdf"),
  peak_length_plot,
  width = 16, height = 8, units = "cm"
)

message("Candidate-peak results saved in: ", out_dir)
