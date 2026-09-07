#!/usr/bin/env Rscript

# Reproducible analysis for the statement:
# "On average, sampling 40% of flowering plant species per network was
# sufficient to capture 80% of pollinator richness with the best strategy."
#
# Best strategy = select flowering plant species from highest to lowest local
# flower abundance. The manuscript-scale estimate is obtained by first finding
# the 80%-richness threshold in each network and then summarising those network
# thresholds. The mean is rounded to the nearest 10 percentage points (and the
# median is reported), which produces the stated 40%. For transparency, the
# stricter alternative--where the across-network mean curve itself first reaches
# 80%--is also reported; with the current data and a 5% grid, that value is 45%.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
})

# =============================================================================
# User settings: edit these three lines when needed
# =============================================================================
project_dir <- "D:/Chap1_TargetPlant_to_monitor/Project_PlantSelection"
output_dir <- "D:/Chap1_TargetPlant_to_monitor/result_260723"
n_random_repetitions <- 50L

count_file <- file.path(project_dir, "data", "processed", "data_count_scaled_published.rds")
interaction_file <- file.path(project_dir, "data", "processed", "data_interact_published.rds")

if (!file.exists(count_file)) stop("Missing input: ", count_file)
if (!file.exists(interaction_file)) stop("Missing input: ", interaction_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use a 5-percentage-point grid so that 40% is represented exactly.
sampling_proportions <- seq(0.05, 1.00, by = 0.05)
target_richness <- 80

data_count_scaled <- readRDS(count_file)
data_interact <- readRDS(interaction_file)

required_count <- c("Study_Network_id", "Plant_species", "Flower_count_scaled")
required_interaction <- c("Study_Network_id", "Plant_original_name", "Pollinator_accepted_name", "Interaction_addup")
if (!all(required_count %in% names(data_count_scaled))) {
  stop("Flower-count data lack: ", paste(setdiff(required_count, names(data_count_scaled)), collapse = ", "))
}
if (!all(required_interaction %in% names(data_interact))) {
  stop("Interaction data lack: ", paste(setdiff(required_interaction, names(data_interact)), collapse = ", "))
}

# One abundance value per plant species and network. max() preserves the value
# used by the existing scripts when identical species-level values are repeated.
plant_pool <- data_count_scaled %>%
  filter(!is.na(Study_Network_id), !is.na(Plant_species), Plant_species != "") %>%
  group_by(Study_Network_id, Plant_species) %>%
  summarise(
    Flower_count_scaled = if (all(is.na(Flower_count_scaled))) 0 else max(Flower_count_scaled, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Study_Network_id, desc(Flower_count_scaled), Plant_species) %>%
  group_by(Study_Network_id) %>%
  mutate(
    abundance_rank = row_number(),
    n_plants_total = n()
  ) %>%
  ungroup()

# Only positive interactions with an identified pollinator count toward richness.
interactions_positive <- data_interact %>%
  filter(
    !is.na(Study_Network_id),
    !is.na(Plant_original_name), Plant_original_name != "",
    !is.na(Pollinator_accepted_name), Pollinator_accepted_name != "",
    !is.na(Interaction_addup), Interaction_addup > 0
  ) %>%
  distinct(Study_Network_id, Plant_original_name, Pollinator_accepted_name)

total_richness <- interactions_positive %>%
  group_by(Study_Network_id) %>%
  summarise(total_pollinator_richness = n_distinct(Pollinator_accepted_name), .groups = "drop")

eligible_networks <- plant_pool %>%
  distinct(Study_Network_id, n_plants_total) %>%
  inner_join(total_richness, by = "Study_Network_id") %>%
  filter(total_pollinator_richness > 0)

if (nrow(eligible_networks) == 0) stop("No networks have both flowering plants and positive pollinator interactions.")

# Evaluate the abundance-ranked selection separately within every network.
curve_by_network <- map_dfr(sampling_proportions, function(prop) {
  selected <- plant_pool %>%
    semi_join(eligible_networks, by = "Study_Network_id") %>%
    group_by(Study_Network_id) %>%
    filter(abundance_rank <= ceiling(prop * first(n_plants_total))) %>%
    ungroup()

  captured <- interactions_positive %>%
    inner_join(
      selected %>% select(Study_Network_id, Plant_species),
      by = c("Study_Network_id", "Plant_original_name" = "Plant_species")
    ) %>%
    group_by(Study_Network_id) %>%
    summarise(captured_pollinator_richness = n_distinct(Pollinator_accepted_name), .groups = "drop")

  eligible_networks %>%
    left_join(
      selected %>% count(Study_Network_id, name = "n_plants_sampled"),
      by = "Study_Network_id"
    ) %>%
    left_join(captured, by = "Study_Network_id") %>%
    mutate(
      sampling_proportion = prop,
      sampling_percent = 100 * prop,
      n_plants_sampled = replace_na(n_plants_sampled, 0L),
      captured_pollinator_richness = replace_na(captured_pollinator_richness, 0L),
      pollinator_richness_captured_percent = 100 * captured_pollinator_richness / total_pollinator_richness
    )
})

# Random benchmark: at each proportion, sample the same number of plants as the
# abundance strategy, repeat within each network, and average captured richness.
set.seed(2025)
random_curve_repetitions <- map_dfr(seq_len(n_random_repetitions), function(rep_id) {
  map_dfr(sampling_proportions, function(prop) {
    selected <- plant_pool %>%
      semi_join(eligible_networks, by = "Study_Network_id") %>%
      group_by(Study_Network_id) %>%
      group_modify(~ slice_sample(.x, n = ceiling(prop * nrow(.x)), replace = FALSE)) %>%
      ungroup()

    captured <- interactions_positive %>%
      inner_join(selected %>% select(Study_Network_id, Plant_species),
                 by = c("Study_Network_id", "Plant_original_name" = "Plant_species")) %>%
      group_by(Study_Network_id) %>%
      summarise(captured_pollinator_richness = n_distinct(Pollinator_accepted_name), .groups = "drop")

    eligible_networks %>%
      left_join(selected %>% count(Study_Network_id, name = "n_plants_sampled"), by = "Study_Network_id") %>%
      left_join(captured, by = "Study_Network_id") %>%
      mutate(
        repetition = rep_id,
        sampling_proportion = prop,
        sampling_percent = 100 * prop,
        captured_pollinator_richness = replace_na(captured_pollinator_richness, 0L),
        pollinator_richness_captured_percent = 100 * captured_pollinator_richness / total_pollinator_richness
      )
  })
})

random_curve_by_network <- random_curve_repetitions %>%
  group_by(Study_Network_id, n_plants_total, total_pollinator_richness,
           sampling_proportion, sampling_percent) %>%
  summarise(
    n_plants_sampled = first(n_plants_sampled),
    pollinator_richness_captured_percent = mean(pollinator_richness_captured_percent),
    random_sd = sd(pollinator_richness_captured_percent),
    .groups = "drop"
  )

mean_curve <- curve_by_network %>%
  group_by(sampling_proportion, sampling_percent) %>%
  summarise(
    n_networks = n(),
    mean_captured_percent = mean(pollinator_richness_captured_percent),
    sd_captured_percent = sd(pollinator_richness_captured_percent),
    se_captured_percent = sd_captured_percent / sqrt(n_networks),
    median_captured_percent = median(pollinator_richness_captured_percent),
    .groups = "drop"
  ) %>%
  arrange(sampling_proportion) %>%
  mutate(
    marginal_gain_percentage_points = mean_captured_percent - lag(mean_captured_percent),
    gain_per_additional_5pct_sampled = marginal_gain_percentage_points
  )

random_mean_curve <- random_curve_repetitions %>%
  group_by(repetition, sampling_proportion, sampling_percent) %>%
  summarise(mean_captured_percent = mean(pollinator_richness_captured_percent), .groups = "drop") %>%
  group_by(sampling_proportion, sampling_percent) %>%
  summarise(
    n_repetitions = n(),
    mean_captured_percent = mean(mean_captured_percent),
    sd_between_repetitions = sd(mean_captured_percent),
    .groups = "drop"
  )

mean_threshold <- mean_curve %>%
  filter(mean_captured_percent >= target_richness) %>%
  slice_head(n = 1)

if (nrow(mean_threshold) == 0) {
  warning("The mean curve never reaches ", target_richness, "% on the evaluated grid.")
  mean_threshold <- tibble(
    sampling_proportion = NA_real_, sampling_percent = NA_real_,
    n_networks = nrow(eligible_networks), mean_captured_percent = NA_real_,
    sd_captured_percent = NA_real_, se_captured_percent = NA_real_,
    median_captured_percent = NA_real_, marginal_gain_percentage_points = NA_real_,
    gain_per_additional_5pct_sampled = NA_real_
  )
}

threshold_by_network <- curve_by_network %>%
  filter(pollinator_richness_captured_percent >= target_richness) %>%
  group_by(Study_Network_id) %>%
  slice_min(sampling_proportion, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    Study_Network_id, n_plants_total, total_pollinator_richness,
    threshold_sampling_proportion = sampling_proportion,
    threshold_sampling_percent = sampling_percent,
    n_plants_sampled_at_threshold = n_plants_sampled,
    captured_percent_at_threshold = pollinator_richness_captured_percent
  )

random_threshold_by_network <- random_curve_by_network %>%
  filter(pollinator_richness_captured_percent >= target_richness) %>%
  group_by(Study_Network_id) %>%
  slice_min(sampling_proportion, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    Study_Network_id, strategy = "Random", n_plants_total, total_pollinator_richness,
    threshold_sampling_proportion = sampling_proportion,
    threshold_sampling_percent = sampling_percent,
    n_plants_sampled_at_threshold = n_plants_sampled,
    captured_percent_at_threshold = pollinator_richness_captured_percent
  )

threshold_by_network <- threshold_by_network %>% mutate(strategy = "Abundance", .before = n_plants_total)
all_richness_thresholds <- bind_rows(threshold_by_network, random_threshold_by_network)

networks_not_reaching_target <- eligible_networks %>%
  anti_join(threshold_by_network, by = "Study_Network_id")

network_threshold_summary <- threshold_by_network %>%
  summarise(
    n_networks_reaching_target = n(),
    n_networks_not_reaching_target = nrow(networks_not_reaching_target),
    mean_threshold_percent_exact = mean(threshold_sampling_percent),
    mean_threshold_percent_rounded_to_10 = round(mean_threshold_percent_exact / 10) * 10,
    median_threshold_percent = median(threshold_sampling_percent),
    q25_threshold_percent = quantile(threshold_sampling_percent, 0.25),
    q75_threshold_percent = quantile(threshold_sampling_percent, 0.75)
  )

summary_result <- tibble(
  definition = c(
    "Manuscript estimate: mean of per-network 80%-richness thresholds, rounded to nearest 10 percentage points",
    "Median of per-network 80%-richness thresholds",
    "Strict alternative: first 5%-grid proportion where the across-network mean richness captured is >= 80%"
  ),
  strategy = "Flower abundance (descending within network)",
  target_pollinator_richness_percent = target_richness,
  threshold_sampling_percent = c(
    network_threshold_summary$mean_threshold_percent_rounded_to_10,
    network_threshold_summary$median_threshold_percent,
    mean_threshold$sampling_percent
  ),
  exact_unrounded_value = c(
    network_threshold_summary$mean_threshold_percent_exact,
    network_threshold_summary$median_threshold_percent,
    mean_threshold$sampling_percent
  ),
  n_networks_contributing = c(
    network_threshold_summary$n_networks_reaching_target,
    network_threshold_summary$n_networks_reaching_target,
    mean_threshold$n_networks
  ),
  note = c(
    "This is the calculation supporting the rounded 40% statement.",
    "The median network requires 40% of its flowering plants.",
    sprintf("At 40%%, the across-network mean captured richness is %.2f%%; this stricter curve reaches 80%% at %.0f%%.",
            mean_curve$mean_captured_percent[mean_curve$sampling_percent == 40],
            mean_threshold$sampling_percent)
  )
)

write_csv(curve_by_network, file.path(output_dir, "richness_curve_by_network.csv"))
write_csv(random_curve_by_network, file.path(output_dir, "richness_random_curve_by_network.csv"))
write_csv(random_curve_repetitions, file.path(output_dir, "richness_random_repetitions.csv"))
write_csv(mean_curve, file.path(output_dir, "richness_curve_across_networks.csv"))
write_csv(random_mean_curve, file.path(output_dir, "richness_random_curve_across_networks.csv"))
write_csv(all_richness_thresholds, file.path(output_dir, "richness_80pct_threshold_by_network.csv"))
write_csv(networks_not_reaching_target, file.path(output_dir, "networks_not_reaching_80pct.csv"))
write_csv(network_threshold_summary, file.path(output_dir, "richness_80pct_network_threshold_distribution.csv"))
write_csv(summary_result, file.path(output_dir, "richness_80pct_mean_threshold_summary.csv"))

p <- ggplot() +
  geom_ribbon(
    data = mean_curve,
    aes(x = sampling_percent, ymin = pmax(0, mean_captured_percent - se_captured_percent),
        ymax = pmin(100, mean_captured_percent + se_captured_percent)),
    fill = "#2AA889", alpha = 0.20
  ) +
  geom_ribbon(data = random_mean_curve,
              aes(x = sampling_percent,
                  ymin = mean_captured_percent - sd_between_repetitions,
                  ymax = mean_captured_percent + sd_between_repetitions),
              fill = "grey65", alpha = 0.25) +
  geom_line(data = mean_curve, aes(sampling_percent, mean_captured_percent, colour = "Abundance"), linewidth = 1.0) +
  geom_point(data = mean_curve, aes(sampling_percent, mean_captured_percent, colour = "Abundance"), size = 2.0) +
  geom_line(data = random_mean_curve, aes(sampling_percent, mean_captured_percent, colour = "Random"), linewidth = 1.0) +
  geom_point(data = random_mean_curve, aes(sampling_percent, mean_captured_percent, colour = "Random"), size = 1.7) +
  geom_hline(yintercept = target_richness, linetype = "dashed", color = "#B03A2E") +
  geom_vline(xintercept = network_threshold_summary$mean_threshold_percent_exact,
             linetype = "dotted", color = "#5B2C6F") +
  annotate(
    "label",
    x = network_threshold_summary$mean_threshold_percent_exact,
    y = target_richness,
    label = sprintf("Mean per-network 80%%\nthreshold: %.2f%% plants",
                    network_threshold_summary$mean_threshold_percent_exact),
    hjust = -0.05, vjust = 1.2, size = 3.6, color = "#5B2C6F"
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
  scale_colour_manual(values = c(Abundance = "#B84E22", Random = "grey45"), name = "Selection strategy") +
  labs(
    x = "Flowering plant species sampled per network (%)",
    y = "Mean pollinator richness captured (%)",
    # title = "Pollinator richness captured by abundance-based plant sampling",
    # subtitle = "Random curve and ribbon: mean +/- SD across repeated selections"
  ) +
  theme_classic(base_size = 12) + theme(
    panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    legend.position = "bottom")

p
ggsave(file.path(output_dir, "pollinator_richness_abundance_sampling_curve.png"), p,
       width = 4.5, height = 4, units = "in", dpi = 300, bg = "white")


cat("\nAnalysis complete.\n")
cat("Networks analysed:", nrow(eligible_networks), "\n")
cat("Networks reaching the 80% target on the grid:",
    network_threshold_summary$n_networks_reaching_target, "\n")
cat("Mean per-network threshold (exact):",
    sprintf("%.2f%%", network_threshold_summary$mean_threshold_percent_exact), "\n")
cat("Mean per-network threshold (nearest 10%; manuscript estimate):",
    sprintf("%.0f%%", network_threshold_summary$mean_threshold_percent_rounded_to_10), "\n")
cat("Median per-network threshold:",
    sprintf("%.0f%%", network_threshold_summary$median_threshold_percent), "\n")
cat("Strict alternative--first proportion where mean richness captured >=", target_richness, "%:",
    sprintf("%.0f%%", mean_threshold$sampling_percent), "\n")
cat("Mean pollinator richness captured there:",
    sprintf("%.2f%%", mean_threshold$mean_captured_percent), "\n")
cat("Outputs:", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")

# -----------------------------------------------------------------------------
# Network structure: recompute Figure-2-style abundance and random curves using
# the project's validated metric implementation, then reproduce the plot and
# threshold table in this output directory.
# -----------------------------------------------------------------------------
network_script <- file.path(project_dir, "scripts", "1.5 MainResults - Proportional plant sampling for network structure.R")
if (!file.exists(network_script)) stop("Missing network-structure script: ", network_script)

network_code <- readLines(network_script, warn = FALSE, encoding = "UTF-8")
network_code <- sub('^result_dir <- .*$',
                    sprintf('result_dir <- "%s"', gsub('\\\\', '/', output_dir)),
                    network_code)
network_code <- sub('^n_random_repetitions <- 50$',
                    sprintf('n_random_repetitions <- %d', n_random_repetitions),
                    network_code)
patched_network_script <- file.path(output_dir, "network_structure_recomputed.R")
writeLines(network_code, patched_network_script, useBytes = TRUE)

old_wd <- getwd()
setwd(project_dir)
network_env <- new.env(parent = globalenv())
tryCatch(source(patched_network_script, local = network_env),
         finally = setwd(old_wd))

# Replace the basic network plot with a Figure-2-style version that has one
# legend and explicitly marks the first observed rho >= 0.8 for each strategy.
network_threshold_lines <- network_env$threshold_results %>%
  filter(is.finite(percentage), is.finite(rho))

network_plot <- ggplot(
  network_env$all_results,
  aes(percentage, rho, colour = strategy, group = strategy)
) +
  geom_hline(yintercept = 0.8, linetype = "dotted", colour = "grey30") +
  geom_ribbon(
    data = network_env$random_results,
    aes(ymin = rho - rho_sd, ymax = rho + rho_sd),
    inherit.aes = TRUE, fill = "grey65", colour = NA, alpha = 0.20,
    show.legend = FALSE
  ) +
  geom_vline(
    data = network_threshold_lines,
    aes(xintercept = percentage, colour = strategy),
    linetype = "dashed", linewidth = 0.55, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1) +
  facet_wrap(~ metric, nrow = 1) +
  scale_colour_manual(values = c(Abundance = "#B84E22", Random = "grey45")) +
  scale_x_continuous(breaks = seq(20, 100, 10), limits = c(20, 100)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
  labs(
    x = "Flowering plant species monitored (%)",
    y = expression("Spearman's " * rho),
    colour = "Selection strategy"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.6),
    axis.line = element_blank(), legend.position = "bottom",
    strip.background = element_blank(), strip.text = element_text(face = "bold")
  )

network_plot
ggsave(file.path(output_dir, "Proportional_sampling_network_structure.png"),
       network_plot, width = 5.5, height = 3, units = "in", dpi = 600, bg = "white")

cat("\nBoth richness and network-structure analyses are complete.\n")
