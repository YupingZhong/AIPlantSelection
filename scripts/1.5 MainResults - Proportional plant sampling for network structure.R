################################################################################
# Proportional plant sampling and representation of network structure
#
# Purpose
#   Test what proportion of flowering plant species must be monitored for
#   subnetwork metrics to show a strong association with full-network metrics.
#
# Strategies
#   1. Abundance: select the locally most abundant proportion of plant species.
#   2. Random: select the same proportion of plant species at random.
#
# Metrics
#   Connectance, NODF and network-level specialisation (H2').
#
# Threshold
#   The smallest observed proportion at which Spearman's rho reaches or exceeds
#   0.8. The threshold is calculated separately for every strategy and metric.
#
# Important
#   This is an independent proportional-sampling analysis. It does not replace
#   or overwrite the fixed-number analysis used for Figure 2.
################################################################################


# ==============================================================================
# 1. Packages and user settings
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bipartite)

set.seed(2025)

data_dir <- "data/processed"
result_dir <- "D:/Chap1_TargetPlant_to_monitor/result_260723"

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Starting at 20% ensures that networks with the minimum of 10 eligible plants
# contribute at least two selected plant species, which is required for NODF.
sampling_proportions <- seq(0.20, 1.00, by = 0.05)

# Number of random selections performed at every sampling proportion.
n_random_repetitions <- 50

rho_threshold <- 0.80


# ==============================================================================
# 2. Load processed data and full-network metrics
# ==============================================================================

data_count_scaled <- readRDS(
  file.path(data_dir, "data_count_scaled_published.rds")
)

data_interact <- readRDS(
  file.path(data_dir, "data_interact_published.rds")
)

all_nestedness_df <- readRDS(
  file.path(data_dir, "all_nestedness_df.rds")
)

all_connectance_df <- readRDS(
  file.path(data_dir, "all_connectance_df.rds")
)

all_H2_df <- readRDS(
  file.path(data_dir, "all_H2_df.rds")
)

full_metrics <- all_nestedness_df %>%
  select(
    Study_Network_id,
    NODF_full = Nestedness_orig
  ) %>%
  left_join(
    all_connectance_df %>%
      select(
        Study_Network_id,
        Connectance_full = Connectance_orig
      ),
    by = "Study_Network_id"
  ) %>%
  left_join(
    all_H2_df %>%
      select(
        Study_Network_id,
        H2_full = H2_orig
      ),
    by = "Study_Network_id"
  )


# ==============================================================================
# 3. Prepare the eligible plant pool
# ==============================================================================

plant_pool <- data_count_scaled %>%
  filter(
    !is.na(Study_Network_id),
    !is.na(Plant_species),
    Plant_species != "",
    !is.na(Flower_count_scaled),
    Flower_count_scaled > 0
  ) %>%
  distinct(
    Study_Network_id,
    Plant_species,
    .keep_all = TRUE
  ) %>%
  group_by(Study_Network_id) %>%
  mutate(
    total_plant_richness = n()
  ) %>%
  ungroup()

plant_pool_check <- plant_pool %>%
  distinct(
    Study_Network_id,
    total_plant_richness
  )

if (any(plant_pool_check$total_plant_richness < 10)) {
  stop(
    paste0(
      "At least one network contains fewer than 10 eligible plants. ",
      "Rerun the data-preparation filters before this analysis."
    )
  )
}

# Restrict every calculation to networks present in the plant survey,
# interaction data and full-network metric table.
analysis_networks <- Reduce(
  intersect,
  list(
    unique(plant_pool$Study_Network_id),
    unique(data_interact$Study_Network_id),
    unique(full_metrics$Study_Network_id)
  )
)

plant_pool <- plant_pool %>%
  filter(Study_Network_id %in% analysis_networks)

data_interact <- data_interact %>%
  filter(Study_Network_id %in% analysis_networks)

full_metrics <- full_metrics %>%
  filter(Study_Network_id %in% analysis_networks)

cat("Networks in proportional-sampling analysis:",
    length(analysis_networks), "\n")


# ==============================================================================
# 4. Metric functions
# ==============================================================================

make_interaction_matrix <- function(df) {

  df_matrix <- df %>%
    filter(
      !is.na(Plant_original_name),
      Plant_original_name != "",
      !is.na(Pollinator_accepted_name),
      Pollinator_accepted_name != "",
      !is.na(Interaction_addup),
      Interaction_addup > 0
    ) %>%
    group_by(
      Plant_original_name,
      Pollinator_accepted_name
    ) %>%
    summarise(
      Interaction_addup = sum(Interaction_addup),
      .groups = "drop"
    )

  if (nrow(df_matrix) == 0) {
    return(NULL)
  }

  mat <- xtabs(
    Interaction_addup ~
      Plant_original_name + Pollinator_accepted_name,
    data = df_matrix
  )

  mat <- as.matrix(mat)
  mat <- mat[rowSums(mat) > 0, colSums(mat) > 0, drop = FALSE]

  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(NULL)
  }

  mat
}


calculate_network_metrics <- function(df) {

  mat <- make_interaction_matrix(df)

  if (is.null(mat)) {
    return(
      tibble(
        Connectance_sub = NA_real_,
        NODF_sub = NA_real_,
        H2_sub = NA_real_,
        n_interacting_plants = 0L
      )
    )
  }

  mat_binary <- (mat > 0) * 1

  connectance <- tryCatch(
    as.numeric(
      bipartite::networklevel(
        mat_binary,
        index = "connectance"
      )
    ),
    error = function(e) NA_real_
  )

  nodf <- if (nrow(mat_binary) > 1 && ncol(mat_binary) > 1) {
    tryCatch(
      as.numeric(
        bipartite::nested(
          mat_binary,
          method = "NODF"
        )
      ),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }

  h2 <- if (nrow(mat) > 1 && ncol(mat) > 1) {
    tryCatch(
      as.numeric(
        bipartite::networklevel(
          mat,
          index = "H2"
        )
      ),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }

  tibble(
    Connectance_sub = connectance,
    NODF_sub = nodf,
    H2_sub = h2,
    n_interacting_plants = nrow(mat)
  )
}


calculate_spearman <- function(data, sub_metric, full_metric) {

  x <- data[[sub_metric]]
  y <- data[[full_metric]]

  valid <- is.finite(x) & is.finite(y)

  if (sum(valid) < 3 ||
      length(unique(x[valid])) < 2 ||
      length(unique(y[valid])) < 2) {
    return(
      tibble(
        rho = NA_real_,
        p_value = NA_real_,
        n_networks = sum(valid)
      )
    )
  }

  test <- suppressWarnings(
    cor.test(
      x[valid],
      y[valid],
      method = "spearman",
      alternative = "two.sided",
      exact = FALSE
    )
  )

  tibble(
    rho = unname(test$estimate),
    p_value = test$p.value,
    n_networks = sum(valid)
  )
}


# ==============================================================================
# 5. Plant selection functions
# ==============================================================================

select_abundant_plants <- function(proportion) {

  plant_pool %>%
    group_by(Study_Network_id) %>%
    group_modify(function(.x, .y) {

      total_richness <- nrow(.x)
      n_to_select <- max(
        2L,
        ceiling(proportion * total_richness)
      )

      .x %>%
        arrange(
          desc(Flower_count_scaled),
          Plant_species
        ) %>%
        slice_head(n = n_to_select) %>%
        mutate(
          total_plant_richness = total_richness,
          n_selected = n_to_select
        )
    }) %>%
    ungroup() %>%
    select(
      Study_Network_id,
      Plant_species,
      total_plant_richness,
      n_selected
    )
}


select_random_plants <- function(proportion) {

  plant_pool %>%
    group_by(Study_Network_id) %>%
    group_modify(function(.x, .y) {

      total_richness <- nrow(.x)
      n_to_select <- max(
        2L,
        ceiling(proportion * total_richness)
      )

      .x %>%
        slice_sample(
          n = n_to_select,
          replace = FALSE
        ) %>%
        mutate(
          total_plant_richness = total_richness,
          n_selected = n_to_select
        )
    }) %>%
    ungroup() %>%
    select(
      Study_Network_id,
      Plant_species,
      total_plant_richness,
      n_selected
    )
}


build_subnetwork_metrics <- function(selected_plants) {

  selected_interactions <- data_interact %>%
    inner_join(
      selected_plants,
      by = c(
        "Study_Network_id",
        "Plant_original_name" = "Plant_species"
      )
    )

  selected_interactions %>%
    group_by(Study_Network_id) %>%
    group_modify(
      ~ calculate_network_metrics(.x)
    ) %>%
    ungroup() %>%
    right_join(
      selected_plants %>%
        distinct(
          Study_Network_id,
          total_plant_richness,
          n_selected
        ),
      by = "Study_Network_id"
    ) %>%
    left_join(
      full_metrics,
      by = "Study_Network_id"
    )
}


summarise_correlations <- function(metric_data) {

  bind_rows(
    calculate_spearman(
      metric_data,
      "Connectance_sub",
      "Connectance_full"
    ) %>%
      mutate(metric = "Connectance"),

    calculate_spearman(
      metric_data,
      "NODF_sub",
      "NODF_full"
    ) %>%
      mutate(metric = "NODF"),

    calculate_spearman(
      metric_data,
      "H2_sub",
      "H2_full"
    ) %>%
      mutate(metric = "H2'" )
  )
}


# ==============================================================================
# 6. Abundance-based proportional sampling
# ==============================================================================

abundance_results <- map_dfr(
  sampling_proportions,
  function(proportion) {

    selected <- select_abundant_plants(proportion)
    metrics <- build_subnetwork_metrics(selected)

    summarise_correlations(metrics) %>%
      mutate(
        strategy = "Abundance",
        proportion = proportion,
        percentage = proportion * 100,
        median_n_selected = median(
          metrics$n_selected,
          na.rm = TRUE
        ),
        min_n_selected = min(
          metrics$n_selected,
          na.rm = TRUE
        ),
        max_n_selected = max(
          metrics$n_selected,
          na.rm = TRUE
        )
      )
  }
)


# ==============================================================================
# 7. Random proportional sampling
# ==============================================================================

random_by_repetition <- map_dfr(
  seq_len(n_random_repetitions),
  function(repetition) {

    map_dfr(
      sampling_proportions,
      function(proportion) {

        selected <- select_random_plants(proportion)
        metrics <- build_subnetwork_metrics(selected)

        summarise_correlations(metrics) %>%
          mutate(
            repetition = repetition,
            strategy = "Random",
            proportion = proportion,
            percentage = proportion * 100,
            median_n_selected = median(
              metrics$n_selected,
              na.rm = TRUE
            ),
            min_n_selected = min(
              metrics$n_selected,
              na.rm = TRUE
            ),
            max_n_selected = max(
              metrics$n_selected,
              na.rm = TRUE
            )
          )
      }
    )
  }
)

random_results <- random_by_repetition %>%
  group_by(
    strategy,
    proportion,
    percentage,
    metric
  ) %>%
  summarise(
    rho = mean(rho, na.rm = TRUE),
    rho_sd = sd(rho, na.rm = TRUE),
    p_value = NA_real_,
    n_networks = round(mean(n_networks, na.rm = TRUE)),
    median_n_selected = median(median_n_selected, na.rm = TRUE),
    min_n_selected = min(min_n_selected, na.rm = TRUE),
    max_n_selected = max(max_n_selected, na.rm = TRUE),
    .groups = "drop"
  )

abundance_results <- abundance_results %>%
  mutate(rho_sd = NA_real_)

all_results <- bind_rows(
  abundance_results,
  random_results
) %>%
  arrange(metric, strategy, proportion)


# ==============================================================================
# 8. First observed proportion reaching rho >= 0.8
# ==============================================================================

threshold_results <- all_results %>%
  group_by(strategy, metric) %>%
  arrange(proportion, .by_group = TRUE) %>%
  filter(
    is.finite(rho),
    rho >= rho_threshold
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(
    strategy,
    metric,
    proportion,
    percentage,
    rho,
    rho_sd,
    n_networks,
    median_n_selected,
    min_n_selected,
    max_n_selected
  )

expected_threshold_rows <- tidyr::expand_grid(
  strategy = c("Abundance", "Random"),
  metric = c("Connectance", "NODF", "H2'")
)

threshold_results <- expected_threshold_rows %>%
  left_join(
    threshold_results,
    by = c("strategy", "metric")
  ) %>%
  arrange(metric, strategy)

print(threshold_results)


# ==============================================================================
# 9. Plot
# ==============================================================================

proportional_plot <- ggplot(
  all_results,
  aes(
    x = percentage,
    y = rho,
    colour = strategy,
    group = strategy
  )
) +
  geom_hline(
    yintercept = rho_threshold,
    linetype = "dotted",
    colour = "grey35"
  ) +
  geom_ribbon(
    data = random_results,
    aes(
      ymin = rho - rho_sd,
      ymax = rho + rho_sd,
      fill = strategy
    ),
    colour = NA,
    alpha = 0.20,
    inherit.aes = TRUE
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(
    ~ metric,
    nrow = 1
  ) +
  scale_colour_manual(
    values = c(
      "Abundance" = "#B84E22",
      "Random" = "grey45"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Random" = "grey65"
    )
  ) +
  scale_x_continuous(
    breaks = seq(20, 100, by = 10),
    limits = c(20, 100)
  ) +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.2)
  ) +
  labs(
    x = "Flowering plant species monitored (%)",
    y = "Spearman's rho",
    colour = "Selection strategy",
    fill = "Selection strategy"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(
      colour = "grey50",
      fill = NA,
      linewidth = 0.6
    ),
    axis.line = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(proportional_plot)


# ==============================================================================
# 10. Save results
# ==============================================================================

saveRDS(
  all_results,
  file.path(
    data_dir,
    "proportional_sampling_correlations.rds"
  )
)

saveRDS(
  random_by_repetition,
  file.path(
    data_dir,
    "proportional_sampling_random_repetitions.rds"
  )
)

write.csv(
  all_results,
  file.path(
    result_dir,
    "Proportional_sampling_correlations.csv"
  ),
  row.names = FALSE
)

write.csv(
  threshold_results,
  file.path(
    result_dir,
    "Proportional_sampling_thresholds.csv"
  ),
  row.names = FALSE
)

ggsave(
  filename = file.path(
    result_dir,
    "Proportional_sampling_network_structure.png"
  ),
  plot = proportional_plot,
  width = 10,
  height = 4.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

cat("\nAnalysis complete.\n")
cat("Threshold table:\n")
print(threshold_results)
