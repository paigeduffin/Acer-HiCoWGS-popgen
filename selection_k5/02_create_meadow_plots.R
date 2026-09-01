# Author: Paige Duffin
# Identify comparison-specific outlier windows and create K = 5 meadow plots.
#
# Outliers are defined independently for each of the 10 pairwise comparisons
# as windows in the top 1% of FST and either tail (lower or upper 1%) of the
# log2 pi-ratio distribution. Negative FST estimates are excluded first.
#
# Usage:
#   Rscript 02_create_meadow_plots.R [FST_FILE] [PI_FILE] [OUT_DIR]

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(scales)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3) {
  stop("Usage: Rscript 02_create_meadow_plots.R [FST_FILE] [PI_FILE] [OUT_DIR]")
}

fst_file <- if (length(args) >= 1) args[1] else "pixy_results/10kb_min5kb_fst.txt"
pi_file <- if (length(args) >= 2) args[2] else "pixy_results/10kb_min5kb_pi.txt"
out_dir <- if (length(args) >= 3) args[3] else "meadow_plot_results"

for (input_file in c(fst_file, pi_file)) {
  if (!file.exists(input_file)) stop("Missing input: ", input_file)
}
if (dir.exists(out_dir)) stop("Output directory already exists: ", out_dir)
dir.create(out_dir, recursive = TRUE)
dir.create(file.path(out_dir, "plots"))

chromosome_order <- paste0("OZ035", 966:979, ".1")

# The first population is the numerator of each pi ratio; the second is the
# denominator. These orientations match the original analysis.
comparison_map <- tibble::tribble(
  ~comparison,    ~population_1,          ~population_2,
  "FLDT.v.ARCU", "Florida_DryTort",      "Aruba_Curacao",
  "FLDT.v.BEME", "Florida_DryTort",      "Belize_Mexico",
  "FLDT.v.DRJA", "Florida_DryTort",      "DominRepub_Jamaica",
  "FLDT.v.PANM", "Florida_DryTort",      "Panama",
  "DRJA.v.ARCU", "DominRepub_Jamaica",  "Aruba_Curacao",
  "DRJA.v.BEME", "DominRepub_Jamaica",  "Belize_Mexico",
  "DRJA.v.PANM", "DominRepub_Jamaica",  "Panama",
  "ARCU.v.BEME", "Aruba_Curacao",       "Belize_Mexico",
  "ARCU.v.PANM", "Aruba_Curacao",       "Panama",
  "PANM.v.BEME", "Panama",              "Belize_Mexico"
)

pair_key <- function(pop_a, pop_b) {
  ifelse(pop_a <= pop_b,
         paste(pop_a, pop_b, sep = "__"),
         paste(pop_b, pop_a, sep = "__"))
}

comparison_map <- comparison_map %>%
  mutate(pair_key = pair_key(population_1, population_2))

fst <- read_tsv(fst_file, show_col_types = FALSE)
pi <- read_tsv(pi_file, show_col_types = FALSE)

required_fst <- c("pop1", "pop2", "chromosome", "window_pos_1",
                  "window_pos_2", "avg_wc_fst")
required_pi <- c("pop", "chromosome", "window_pos_1", "window_pos_2", "avg_pi")
if (!all(required_fst %in% names(fst))) {
  stop("FST file is missing: ", paste(setdiff(required_fst, names(fst)), collapse = ", "))
}
if (!all(required_pi %in% names(pi))) {
  stop("Pi file is missing: ", paste(setdiff(required_pi, names(pi)), collapse = ", "))
}

fst_prepared <- fst %>%
  filter(chromosome %in% chromosome_order) %>%
  mutate(pair_key = pair_key(pop1, pop2)) %>%
  inner_join(comparison_map %>% select(comparison, pair_key), by = "pair_key") %>%
  select(chromosome, window_pos_1, window_pos_2, comparison, avg_wc_fst)

pi_wide <- pi %>%
  filter(chromosome %in% chromosome_order) %>%
  select(pop, chromosome, window_pos_1, window_pos_2, avg_pi) %>%
  pivot_wider(names_from = pop, values_from = avg_pi)

make_pi_comparison <- function(comparison, population_1, population_2) {
  pi_wide %>%
    transmute(
      chromosome,
      window_pos_1,
      window_pos_2,
      comparison = comparison,
      population_1 = population_1,
      population_2 = population_2,
      pi_population_1 = .data[[population_1]],
      pi_population_2 = .data[[population_2]],
      pi_ratio = pi_population_1 / pi_population_2,
      log2_pi_ratio = log2(pi_ratio)
    )
}

pi_comparisons <- Map(
  make_pi_comparison,
  comparison_map$comparison,
  comparison_map$population_1,
  comparison_map$population_2
) %>% bind_rows()

meadow_data <- pi_comparisons %>%
  inner_join(
    fst_prepared,
    by = c("chromosome", "window_pos_1", "window_pos_2", "comparison")
  ) %>%
  filter(!is.na(avg_wc_fst), !is.na(log2_pi_ratio), avg_wc_fst >= 0) %>%
  group_by(comparison) %>%
  mutate(
    fst_99 = quantile(avg_wc_fst, 0.99, na.rm = TRUE),
    pi_ratio_01 = quantile(log2_pi_ratio, 0.01, na.rm = TRUE),
    pi_ratio_99 = quantile(log2_pi_ratio, 0.99, na.rm = TRUE),
    outlier = avg_wc_fst >= fst_99 &
      (log2_pi_ratio <= pi_ratio_01 | log2_pi_ratio >= pi_ratio_99),
    outlier_driver = case_when(
      outlier & log2_pi_ratio < 0 ~ population_1,
      outlier & log2_pi_ratio > 0 ~ population_2,
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(
    chromosome_number = match(chromosome, chromosome_order),
    chromosome_number = factor(chromosome_number, levels = 1:14)
  ) %>%
  arrange(comparison, chromosome_number, window_pos_1)

write_csv(meadow_data, file.path(out_dir, "meadow_window_statistics.csv"))
write_csv(filter(meadow_data, outlier),
          file.path(out_dir, "candidate_outlier_windows.csv"))

outlier_counts <- meadow_data %>%
  filter(outlier) %>%
  count(comparison, outlier_driver, name = "n_outlier_windows")
write_csv(outlier_counts, file.path(out_dir, "outlier_window_counts_by_driver.csv"))

auto_pi_scale <- function(values, legend_title = expression(log[2](pi[1] / pi[2]))) {
  finite_values <- values[is.finite(values)]
  limits <- if (length(finite_values) == 0) c(-1, 1) else range(c(finite_values, 0))
  if (diff(limits) == 0) limits <- limits + c(-1e-9, 1e-9)

  negative_colors <- c(
    "#3E207B", "#4C3194", "#5A42AD", "#5E63CC", "#5970D2",
    "#5C89DC", "#94BBEC", "#A3CCD0", "#9CCD9E", "#92D050"
  )
  positive_colors <- c(
    "#92D050", "#CAD961", "#FFE171", "#FACA42", "#F3A532",
    "#EB8E2D", "#E2662A", "#BE3D28", "#B12F27", "#752A24"
  )
  zero_position <- scales::rescale(0, from = limits)
  n_colors <- 50
  colors <- c(
    colorRampPalette(negative_colors)(ceiling(n_colors / 2)),
    colorRampPalette(positive_colors)(ceiling(n_colors / 2))[-1]
  )
  positions <- c(
    seq(0, zero_position, length.out = ceiling(n_colors / 2)),
    seq(zero_position, 1, length.out = ceiling(n_colors / 2))[-1]
  )

  scale_color_gradientn(
    name = legend_title,
    colors = colors,
    values = positions,
    limits = limits,
    oob = scales::squish
  )
}

make_meadow_plot <- function(comparison_name) {
  plot_data <- meadow_data %>%
    filter(comparison == comparison_name) %>%
    mutate(
      point_size = if_else(outlier, 1.6, 0.45),
      point_alpha = if_else(outlier, 1, 0.6)
    )
  finite_data <- filter(plot_data, is.finite(log2_pi_ratio))
  negative_infinite <- filter(plot_data, log2_pi_ratio == -Inf)
  positive_infinite <- filter(plot_data, log2_pi_ratio == Inf)

  ggplot() +
    geom_point(
      data = finite_data,
      aes(window_pos_1, avg_wc_fst, color = log2_pi_ratio,
          size = point_size, alpha = point_alpha)
    ) +
    geom_point(
      data = negative_infinite,
      aes(window_pos_1, avg_wc_fst, size = point_size, alpha = point_alpha),
      shape = 8, color = "#6579D7"
    ) +
    geom_point(
      data = positive_infinite,
      aes(window_pos_1, avg_wc_fst, size = point_size, alpha = point_alpha),
      shape = 8, color = "#D55C2C"
    ) +
    auto_pi_scale(finite_data$log2_pi_ratio) +
    scale_size_identity() +
    scale_alpha_identity() +
    coord_cartesian(ylim = c(0, 1)) +
    facet_grid(
      ~ chromosome_number,
      scales = "free_x", space = "free_x", switch = "x"
    ) +
    labs(title = comparison_name, x = "Chromosome", y = expression(F[ST])) +
    theme_classic() +
    theme(
      panel.spacing.x = grid::unit(0.03, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      legend.position = "right"
    )
}

comparison_order <- c(
  "DRJA.v.ARCU", "FLDT.v.BEME", "FLDT.v.DRJA", "PANM.v.BEME",
  "ARCU.v.PANM", "DRJA.v.PANM", "DRJA.v.BEME", "FLDT.v.PANM",
  "FLDT.v.ARCU", "ARCU.v.BEME"
)

meadow_plots <- lapply(comparison_order, make_meadow_plot)
names(meadow_plots) <- comparison_order

for (comparison_name in names(meadow_plots)) {
  ggsave(
    file.path(out_dir, "plots", paste0("meadow_", comparison_name, ".png")),
    meadow_plots[[comparison_name]],
    width = 21, height = 8, units = "cm", dpi = 600, bg = "transparent"
  )
}

combined_plot <- arrangeGrob(grobs = meadow_plots, ncol = 1)
ggsave(
  file.path(out_dir, "plots", "meadow_plots_all_comparisons.png"),
  combined_plot,
  width = 21, height = 70, units = "cm", dpi = 600, bg = "transparent"
)
ggsave(
  file.path(out_dir, "plots", "meadow_plots_all_comparisons.pdf"),
  combined_plot,
  width = 21, height = 70, units = "cm"
)

message("Meadow-plot results saved in: ", out_dir)
