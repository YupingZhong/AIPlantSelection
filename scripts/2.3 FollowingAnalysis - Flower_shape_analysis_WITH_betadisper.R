################################################################################
# Following analysis: Do different flower shapes attract different pollinators?
#
# Main questions
#   1. FLOWER-SHAPE EFFECT
#      After accounting for study identity and restricting permutations within
#      networks, is flower shape associated with pollinator-genus composition?
#
#   2. WITHIN-SHAPE HETEROGENEITY
#      Within the same network, do plant species assigned to the same flower-
#      shape category nevertheless support different pollinator assemblages?
#
# Interpretation
#   - Question 1 is tested with a network-adjusted, blocked PERMANOVA.
#   - Question 2 is tested using pairwise Sorensen dissimilarity calculated
#     separately within every network x flower-shape combination.
#   - The former global betadisper analysis is not used to claim within-network
#     variation because its distances are measured from global shape centroids.
#
# Cache workflow
#   FIRST COMPLETE RUN:
#       RUN_MULTIVARIATE      <- TRUE
#       RUN_BETADISPER         <- TRUE
#       RUN_WITHIN_SHAPE_BETA  <- TRUE
#
#   LATER PLOTTING RUNS:
#       Change all three settings to FALSE. Saved RDS/CSV files will be loaded,
#       so plots can be edited without rerunning PERMANOVA, betadisper, or beta
#       diversity.
################################################################################


# ==============================================================================
# 1. Packages and user settings
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(vegan)
  library(permute)
  library(betapart)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(multcomp)
  library(multcompView)
  library(ggplot2)
  library(viridis)
  library(cowplot)
  library(readr)
  # Reload dplyr last because MASS/multcomp dependencies export functions with
  # overlapping names such as select().
  library(dplyr)
})

project_dir <- "D:/Chap1_TargetPlant_to_monitor/Project_PlantSelection"
result_dir  <- "D:/Chap1_TargetPlant_to_monitor/result_260723"
cache_dir   <- file.path(project_dir, "data", "processed", "flower_shape_cache")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Set TRUE only when the analysis must be recalculated.
RUN_MULTIVARIATE <- FALSE
RUN_BETADISPER <- FALSE
RUN_WITHIN_SHAPE_BETA <- FALSE

# Default = 999 for the manuscript. FLOWER_SHAPE_NPERM is only an optional
# command-line/environment override for quick diagnostic runs.
N_PERMUTATIONS <- as.integer(
  Sys.getenv("FLOWER_SHAPE_NPERM", unset = "999")
)
RANDOM_SEED <- 2025
set.seed(RANDOM_SEED)

# Excluded because these categories contain too few usable observations in the
# current analysis. Keep this definition identical across all sections.
excluded_shapes <- c("trap flowers", "brush flowers")


# ==============================================================================
# 2. Load and prepare data
# ==============================================================================

data_merge <- readRDS(
  file.path(project_dir, "data", "processed", "data_merge.rds")
)

all_data <- data_merge %>%
  filter(
    !is.na(Study_id),
    !is.na(Study_Network_id),
    !is.na(Plant_accepted_name),
    !is.na(Pollinator_genus),
    !is.na(flw_shape_revised),
    !flw_shape_revised %in% excluded_shapes
  )

# Analysis unit: one plant species x network combination.
# Pollinator composition is recorded at genus level using presence/absence.
all_data_unique <- all_data %>%
  distinct(
    Plant_accepted_name,
    Study_Network_id,
    Pollinator_genus,
    .keep_all = TRUE
  )

community_long <- all_data_unique %>%
  transmute(
    Study_id,
    Study_Network_id,
    Plant_accepted_name,
    FlowerShape = flw_shape_revised,
    Pollinator_genus,
    presence = 1L,
    SampleID = paste(
      Plant_accepted_name,
      Study_Network_id,
      sep = "__NETWORK__"
    )
  ) %>%
  group_by(
    Study_id,
    Study_Network_id,
    Plant_accepted_name,
    FlowerShape,
    Pollinator_genus,
    SampleID
  ) %>%
  summarise(presence = max(presence), .groups = "drop")

meta <- community_long %>%
  distinct(
    SampleID,
    Study_id,
    Study_Network_id,
    Plant_accepted_name,
    FlowerShape
  ) %>%
  arrange(SampleID) %>%
  mutate(
    Study_id = factor(Study_id),
    Study_Network_id = factor(Study_Network_id),
    FlowerShape = factor(FlowerShape)
  )

# A SampleID must correspond to exactly one plant species in one network.
if (anyDuplicated(meta$SampleID)) {
  stop("SampleID is not unique. Check plant names and network identifiers.")
}

community_wide <- community_long %>%
  dplyr::select(SampleID, Pollinator_genus, presence) %>%
  pivot_wider(
    names_from = Pollinator_genus,
    values_from = presence,
    values_fill = 0
  )

comm_pa <- community_wide %>%
  arrange(match(SampleID, meta$SampleID)) %>%
  tibble::column_to_rownames("SampleID") %>%
  as.matrix()

storage.mode(comm_pa) <- "numeric"
comm_pa <- (comm_pa > 0) * 1

stopifnot(identical(rownames(comm_pa), meta$SampleID))

cat("Plant species x network observations:", nrow(meta), "\n")
cat("Networks:", n_distinct(meta$Study_Network_id), "\n")
cat("Flower-shape categories:", n_distinct(meta$FlowerShape), "\n")
cat("Pollinator genera:", ncol(comm_pa), "\n")


# ==============================================================================
# 3. Question 1: flower-shape effect on pollinator composition
# ==============================================================================
#
# Study_id is included as a covariate to account for differences among studies.
# Network identity is handled through restricted permutations: observations may
# only be permuted within their own network. This is a blocked PERMANOVA, not a
# separate PERMANOVA for every network. FlowerShape is the predictor of interest.

permanova_cache <- file.path(cache_dir, "permanova_network_blocked.rds")
pairwise_cache  <- file.path(cache_dir, "pairwise_permanova_network_blocked.rds")

if (RUN_MULTIVARIATE ||
    !file.exists(permanova_cache) ||
    !file.exists(pairwise_cache)) {

  permutation_control <- how(nperm = N_PERMUTATIONS)
  setBlocks(permutation_control) <- meta$Study_Network_id

  adonis_res <- adonis2(
    comm_pa ~ Study_id + FlowerShape,
    data = meta,
    method = "jaccard",
    binary = TRUE,
    permutations = permutation_control,
    by = "margin"
  )

  saveRDS(adonis_res, permanova_cache)

  # ---------------------------------------------------------------------------
  # Pairwise comparisons between flower-shape categories
  # ---------------------------------------------------------------------------
  # Each comparison uses the same study adjustment and within-network
  # permutation restriction as the overall PERMANOVA.

  shape_pairs <- combn(
    levels(droplevels(meta$FlowerShape)),
    2,
    simplify = FALSE
  )

  pairwise_permanova <- map_dfr(shape_pairs, function(pair) {

    keep <- meta$FlowerShape %in% pair
    meta_pair <- droplevels(meta[keep, , drop = FALSE])
    mat_pair <- comm_pa[keep, , drop = FALSE]

    # A valid blocked comparison requires at least one network containing both
    # flower-shape categories.
    paired_networks <- meta_pair %>%
      distinct(Study_Network_id, FlowerShape) %>%
      count(Study_Network_id, name = "n_shapes") %>%
      filter(n_shapes == 2) %>%
      pull(Study_Network_id)

    keep_network <- meta_pair$Study_Network_id %in% paired_networks
    meta_pair <- droplevels(meta_pair[keep_network, , drop = FALSE])
    mat_pair <- mat_pair[keep_network, , drop = FALSE]

    if (nrow(meta_pair) < 4 || length(paired_networks) < 2) {
      return(tibble(
        Shape_1 = pair[1], Shape_2 = pair[2],
        N_networks = length(paired_networks),
        F_value = NA_real_, R2 = NA_real_, P_value = NA_real_
      ))
    }

    pair_control <- how(nperm = N_PERMUTATIONS)
    setBlocks(pair_control) <- meta_pair$Study_Network_id

    fit <- tryCatch(
      adonis2(
        mat_pair ~ Study_id + FlowerShape,
        data = meta_pair,
        method = "jaccard",
        binary = TRUE,
        permutations = pair_control,
        by = "margin"
      ),
      error = function(e) NULL
    )

    if (is.null(fit) || !"FlowerShape" %in% rownames(fit)) {
      return(tibble(
        Shape_1 = pair[1], Shape_2 = pair[2],
        N_networks = length(paired_networks),
        F_value = NA_real_, R2 = NA_real_, P_value = NA_real_
      ))
    }

    shape_row <- which(rownames(fit) == "FlowerShape")

    tibble(
      Shape_1 = pair[1],
      Shape_2 = pair[2],
      N_networks = length(paired_networks),
      F_value = fit$F[shape_row],
      R2 = fit$R2[shape_row],
      P_value = fit$`Pr(>F)`[shape_row]
    )
  }) %>%
    mutate(
      P_adjusted = p.adjust(P_value, method = "BH"),
      Significance = case_when(
        is.na(P_adjusted) ~ NA_character_,
        P_adjusted < 0.001 ~ "***",
        P_adjusted < 0.01  ~ "**",
        P_adjusted < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    arrange(P_adjusted)

  saveRDS(pairwise_permanova, pairwise_cache)

} else {

  adonis_res <- readRDS(permanova_cache)
  pairwise_permanova <- readRDS(pairwise_cache)
}

print(adonis_res)
print(pairwise_permanova)

write_csv(
  as.data.frame(adonis_res) %>% tibble::rownames_to_column("Term"),
  file.path(result_dir, "Flower_shape_PERMANOVA_network_blocked.csv")
)

write_csv(
  pairwise_permanova,
  file.path(result_dir, "Flower_shape_pairwise_PERMANOVA_network_blocked.csv")
)

# Extract the flower-shape row for manuscript reporting. The model accounts
# for study-level variation, while inference for FlowerShape is based on
# permutations restricted within network.
shape_test <- as.data.frame(adonis_res) %>%
  tibble::rownames_to_column("Term") %>%
  filter(Term == "FlowerShape")

if (nrow(shape_test) == 1) {
  cat("\nFlower-shape PERMANOVA result:\n")
  cat("R2 =", round(shape_test$R2, 4), "\n")
  cat("F  =", round(shape_test$F, 3), "\n")
  cat("P  =", shape_test$`Pr(>F)`, "\n")
}


# ==============================================================================
# 3b. PERMANOVA diagnostic: multivariate homogeneity of dispersion
# ==============================================================================
#
# Purpose:
#   PERMANOVA can be significant because group centroids differ, because group
#   dispersions differ, or both. betadisper is therefore retained as a formal
#   diagnostic for interpreting the PERMANOVA.
#
# Important interpretation:
#   Each distance is measured from a plant species x network observation to the
#   GLOBAL centroid of its flower-shape category. These distances must not be
#   described as direct within-network dissimilarities. The within-network
#   ecological question is addressed separately in Sections 4-5.
#
# Design:
#   - The same Jaccard distance matrix used conceptually by PERMANOVA is used.
#   - bias.adjust = TRUE reduces small-sample bias caused by unequal group size.
#   - Permutations are restricted within Study_Network_id, matching the blocked
#     inference used for the FlowerShape term in the PERMANOVA.

betadisper_cache <- file.path(cache_dir, "betadisper_network_blocked.rds")

if (RUN_BETADISPER || !file.exists(betadisper_cache)) {

  jaccard_dist <- vegan::vegdist(
    comm_pa,
    method = "jaccard",
    binary = TRUE
  )

  dispersion_fit <- vegan::betadisper(
    jaccard_dist,
    group = meta$FlowerShape,
    type = "centroid",
    bias.adjust = TRUE
  )

  dispersion_control <- permute::how(nperm = N_PERMUTATIONS)
  setBlocks(dispersion_control) <- meta$Study_Network_id

  dispersion_test <- vegan::permutest(
    dispersion_fit,
    permutations = dispersion_control,
    pairwise = TRUE
  )

  dispersion_data <- meta %>%
    mutate(
      Distance_to_global_centroid = as.numeric(dispersion_fit$distances)
    )

  saveRDS(
    list(
      fit = dispersion_fit,
      test = dispersion_test,
      data = dispersion_data
    ),
    betadisper_cache
  )

} else {

  dispersion_objects <- readRDS(betadisper_cache)
  dispersion_fit <- dispersion_objects$fit
  dispersion_test <- dispersion_objects$test
  dispersion_data <- dispersion_objects$data
}

cat("\nBlocked betadisper diagnostic:\n")
print(dispersion_test)

# Overall permutation test. The Groups row is the flower-shape test to report.
dispersion_omnibus <- as.data.frame(dispersion_test$tab) %>%
  tibble::rownames_to_column("Term")

write_csv(
  dispersion_omnibus,
  file.path(result_dir, "Supplement_betadisper_omnibus_network_blocked.csv")
)

# Save every plant species x network distance to its GLOBAL flower-shape
# centroid. This file supports the supplementary diagnostic figure.
write_csv(
  dispersion_data,
  file.path(result_dir, "Supplement_betadisper_distances_to_global_centroid.csv")
)

# Descriptive summary for each flower-shape category.
dispersion_summary <- dispersion_data %>%
  group_by(FlowerShape) %>%
  summarise(
    N_observations = n(),
    Mean_distance = mean(Distance_to_global_centroid, na.rm = TRUE),
    SD_distance = sd(Distance_to_global_centroid, na.rm = TRUE),
    Median_distance = median(Distance_to_global_centroid, na.rm = TRUE),
    Q1_distance = quantile(
      Distance_to_global_centroid, 0.25, na.rm = TRUE
    ),
    Q3_distance = quantile(
      Distance_to_global_centroid, 0.75, na.rm = TRUE
    ),
    .groups = "drop"
  )

write_csv(
  dispersion_summary,
  file.path(result_dir, "Supplement_betadisper_flower_shape_summary.csv")
)

# Pairwise results returned by permutest.betadisper are named vectors. Build the
# two flower-shape columns explicitly; coercing these vectors to tables would
# incorrectly place the complete comparison name in one column.
shape_comparisons <- t(combn(levels(meta$FlowerShape), 2))

dispersion_pairwise <- tibble(
  Shape_1 = shape_comparisons[, 1],
  Shape_2 = shape_comparisons[, 2],
  Observed_P_value = unname(dispersion_test$pairwise$observed),
  Permutation_P_value = unname(dispersion_test$pairwise$permuted)
) %>%
  mutate(
    P_adjusted_BH = p.adjust(Permutation_P_value, method = "BH")
  ) %>%
  arrange(P_adjusted_BH)

write_csv(
  dispersion_pairwise,
  file.path(result_dir, "Supplement_betadisper_pairwise_network_blocked.csv")
)

# Compact-letter display for the supplementary figure. Flower shapes sharing a
# letter do not differ at alpha = 0.05 according to BH-adjusted, blocked
# pairwise PERMDISP comparisons.
dispersion_letter_p <- dispersion_pairwise$P_adjusted_BH
names(dispersion_letter_p) <- paste(
  dispersion_pairwise$Shape_1,
  dispersion_pairwise$Shape_2,
  sep = "-"
)

dispersion_letters <- multcompView::multcompLetters(
  dispersion_letter_p,
  threshold = 0.05
)$Letters %>%
  tibble::enframe(
    name = "FlowerShape",
    value = "Letters"
  )

# ==============================================================================
# 4. Question 2: variation within network x flower-shape combinations
# ==============================================================================
#
# This section directly answers:
#   Within one network, how dissimilar are the pollinator assemblages of plant
#   species assigned to the same flower-shape category?
#
# Each output row represents one network x flower-shape combination.
# Combinations containing fewer than two plant species cannot yield a pairwise
# distance and are excluded.
#
# beta.sor = total Sorensen dissimilarity
# beta.sim = turnover component
# beta.sne = nestedness-resultant component
# beta.sor = beta.sim + beta.sne

within_beta_cache <- file.path(cache_dir, "within_network_shape_beta.rds")

if (RUN_WITHIN_SHAPE_BETA || !file.exists(within_beta_cache)) {

  network_shape_groups <- meta %>%
    count(
      Study_id,
      Study_Network_id,
      FlowerShape,
      name = "N_plants"
    ) %>%
    filter(N_plants >= 2)

  within_shape_beta <- pmap_dfr(
    network_shape_groups,
    function(Study_id, Study_Network_id, FlowerShape, N_plants) {

      ids <- meta %>%
        filter(
          .data$Study_Network_id == .env$Study_Network_id,
          .data$FlowerShape == .env$FlowerShape
        ) %>%
        pull(SampleID)

      sub_comm <- comm_pa[ids, , drop = FALSE]

      # Remove pollinator columns absent from this network x shape subset.
      sub_comm <- sub_comm[
        ,
        colSums(sub_comm) > 0,
        drop = FALSE
      ]

      if (nrow(sub_comm) < 2 || ncol(sub_comm) == 0) {
        return(tibble())
      }

      beta_parts <- tryCatch(
        betapart::beta.pair(
          sub_comm,
          index.family = "sorensen"
        ),
        error = function(e) NULL
      )

      if (is.null(beta_parts)) return(tibble())

      total_values <- as.matrix(beta_parts$beta.sor)[upper.tri(as.matrix(beta_parts$beta.sor))]
      turnover_values <- as.matrix(beta_parts$beta.sim)[upper.tri(as.matrix(beta_parts$beta.sim))]
      nestedness_values <- as.matrix(beta_parts$beta.sne)[upper.tri(as.matrix(beta_parts$beta.sne))]

      tibble(
        Study_id = as.character(Study_id),
        Study_Network_id = as.character(Study_Network_id),
        FlowerShape = as.character(FlowerShape),
        N_plants = N_plants,
        N_pairs = length(total_values),
        Mean_Sorensen = mean(total_values, na.rm = TRUE),
        Mean_turnover = mean(turnover_values, na.rm = TRUE),
        Mean_nestedness = mean(nestedness_values, na.rm = TRUE)
      )
    }
  ) %>%
    filter(
      is.finite(Mean_Sorensen),
      is.finite(Mean_turnover),
      is.finite(Mean_nestedness)
    )

  saveRDS(within_shape_beta, within_beta_cache)

} else {

  within_shape_beta <- readRDS(within_beta_cache)
}

write_csv(
  within_shape_beta,
  file.path(result_dir, "Within_network_flower_shape_Sorensen.csv")
)

cat("\nNetwork x flower-shape combinations:", nrow(within_shape_beta), "\n")
cat("Networks represented:", n_distinct(within_shape_beta$Study_Network_id), "\n")


# ==============================================================================
# 5. Mixed model for within-shape Sorensen dissimilarity
# ==============================================================================
#
# Response: mean pairwise dissimilarity among plant species sharing the same
# flower shape within the same network.
#
# Random effects account for the hierarchical grouping of observations by study
# and network. No weights are used, so each network x shape combination has equal
# influence. N_pairs is retained in the output for sensitivity checks.

within_shape_beta <- within_shape_beta %>%
  mutate(
    Study_id = factor(Study_id),
    Study_Network_id = factor(Study_Network_id),
    FlowerShape = factor(FlowerShape)
  )

within_shape_model <- lmer(
  Mean_Sorensen ~ FlowerShape +
    (1 | Study_id) +
    (1 | Study_Network_id),
  data = within_shape_beta,
  REML = TRUE
)

print(anova(within_shape_model))

within_shape_emm <- emmeans(
  within_shape_model,
  ~ FlowerShape
)

within_shape_pairs <- pairs(
  within_shape_emm,
  adjust = "sidak"
)

within_shape_letters <- multcomp::cld(
  within_shape_emm,
  adjust = "sidak",
  Letters = letters,
  sort = FALSE
) %>%
  as.data.frame() %>%
  transmute(
    FlowerShape,
    Letters = trimws(.group)
  )

within_shape_emm_df <- as.data.frame(within_shape_emm) %>%
  left_join(within_shape_letters, by = "FlowerShape")

write_csv(
  as.data.frame(within_shape_pairs),
  file.path(result_dir, "Within_network_shape_Sorensen_pairwise.csv")
)

write_csv(
  within_shape_emm_df,
  file.path(result_dir, "Within_network_shape_Sorensen_emmeans.csv")
)


# ==============================================================================
# 6. Pollinator functional-group composition by flower shape
# ==============================================================================
#
# This panel provides an intuitive guild-level description. Proportions are
# calculated within each network x flower-shape combination and then averaged
# across networks, preventing interaction-rich networks from dominating.

data_functional <- data_merge %>%
  filter(
    !is.na(flw_shape_revised),
    !is.na(Pollinator_accepted_name),
    !is.na(Plant_accepted_name),
    !is.na(Study_Network_id),
    !flw_shape_revised %in% excluded_shapes,
    Pollinator_genus != "Apis"
  ) %>%
  mutate(
    functional_group = case_when(
      Pollinator_family == "Syrphidae" ~ "Syrphidae",
      Pollinator_family %in% c(
        "Apidae", "Halictidae", "Andrenidae", "Megachilidae",
        "Colletidae", "Melittidae", "Stenotritidae"
      ) ~ "Bees",
      Pollinator_order == "Hymenoptera" ~ "Non-bee Hymenoptera",
      Pollinator_order == "Diptera" ~ "Non-syrphid Diptera",
      Pollinator_order == "Lepidoptera" ~ "Lepidoptera",
      Pollinator_order == "Coleoptera" ~ "Coleoptera",
      TRUE ~ "Other"
    )
  )

group_order <- c(
  "Bees", "Non-bee Hymenoptera", "Syrphidae",
  "Non-syrphid Diptera", "Lepidoptera", "Coleoptera", "Other"
)

group_cols <- c(
  "Bees" = "#F5E066",
  "Syrphidae" = "#85CCAE",
  "Non-bee Hymenoptera" = "#F5A88A",
  "Coleoptera" = "#E8A3D1",
  "Lepidoptera" = "#A5B5D9",
  "Non-syrphid Diptera" = "#B8DD7F",
  "Other" = "#AAAAAA"
)

network_shape <- data_functional %>%
  distinct(Study_Network_id, flw_shape_revised)

network_fun <- data_functional %>%
  group_by(
    Study_Network_id,
    flw_shape_revised,
    functional_group
  ) %>%
  summarise(
    n_interactions = sum(Interaction_addup, na.rm = TRUE),
    .groups = "drop"
  )

network_props <- network_shape %>%
  crossing(functional_group = group_order) %>%
  left_join(
    network_fun,
    by = c(
      "Study_Network_id",
      "flw_shape_revised",
      "functional_group"
    )
  ) %>%
  mutate(n_interactions = replace_na(n_interactions, 0)) %>%
  group_by(Study_Network_id, flw_shape_revised) %>%
  mutate(
    total_interactions = sum(n_interactions),
    proportion = if_else(
      total_interactions > 0,
      n_interactions / total_interactions,
      NA_real_
    )
  ) %>%
  ungroup() %>%
  filter(is.finite(proportion))

plot_df <- network_props %>%
  group_by(flw_shape_revised, functional_group) %>%
  summarise(
    mean_prop = mean(proportion),
    n_networks = n_distinct(Study_Network_id),
    .groups = "drop"
  )

shape_order <- plot_df %>%
  filter(functional_group == "Bees") %>%
  arrange(mean_prop) %>%
  pull(flw_shape_revised)

plot_df <- plot_df %>%
  mutate(
    flw_shape_revised = factor(
      flw_shape_revised,
      levels = shape_order
    ),
    functional_group = factor(
      functional_group,
      levels = group_order
    )
  )


# ==============================================================================
# 6b. Supplementary Figure: global multivariate dispersion diagnostic
# ==============================================================================
#
# This is retained only as a PERMANOVA diagnostic. Unlike Figure 4b, it does not
# measure dissimilarity among plant species within the same network.

dispersion_plot_data <- dispersion_data %>%
  mutate(
    FlowerShape = factor(
      as.character(FlowerShape),
      levels = shape_order
    )
  ) %>%
  filter(!is.na(FlowerShape))

dispersion_letter_data <- dispersion_plot_data %>%
  group_by(FlowerShape) %>%
  summarise(
    label_position = max(Distance_to_global_centroid, na.rm = TRUE) + 0.025,
    .groups = "drop"
  ) %>%
  mutate(FlowerShape = as.character(FlowerShape)) %>%
  left_join(dispersion_letters, by = "FlowerShape") %>%
  mutate(
    FlowerShape = factor(FlowerShape, levels = shape_order)
  )

supp_betadisper_plot <- ggplot(
  dispersion_plot_data,
  aes(
    x = FlowerShape,
    y = Distance_to_global_centroid
  )
) +
  geom_boxplot(
    fill = "grey75",
    width = 0.65,
    outlier.shape = NA,
    linewidth = 0.35
  ) +
  geom_jitter(
    width = 0.12,
    size = 0.55,
    alpha = 0.22,
    colour = "#4D4D4D"
  ) +
  geom_text(
    data = dispersion_letter_data,
    aes(
      x = FlowerShape,
      y = label_position,
      label = Letters
    ),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Distance to global flower-shape centroid"
  ) +
  theme_classic(base_size = 10) +
  theme(
    panel.border = element_rect(
      colour = "grey55",
      fill = NA,
      linewidth = 0.6
    ),
    axis.line = element_blank(),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.text = element_text(colour = "black"),
    legend.position = "none"
  )

print(supp_betadisper_plot)

ggsave(
  file.path(result_dir, "Supplementary_betadisper_global_centroids.png"),
  supp_betadisper_plot,
  width = 5,
  height = 3,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(result_dir, "Supplementary_betadisper_global_centroids.pdf"),
  supp_betadisper_plot,
  width = 7.5,
  height = 4.5,
  units = "in",
  device = cairo_pdf
)


# ==============================================================================
# 7. Figure 4a: functional-group composition
# ==============================================================================

# Omit Other from display only; retain the original proportion denominator.
plot_group_order <- setdiff(group_order, "Other")
plot_df_display <- plot_df %>%
  filter(functional_group != "Other") %>%
  mutate(functional_group = factor(functional_group, levels = plot_group_order))

visit_group <- ggplot(
  plot_df_display,
  aes(
    x = flw_shape_revised,
    y = mean_prop,
    fill = functional_group
  )
) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_fill_manual(
    values = group_cols[plot_group_order],
    breaks = plot_group_order,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Mean proportion of pollinator visits",
    fill = "Pollinator group"
  ) +
  theme_classic(base_size = 10) +
  theme(
    panel.border = element_rect(
      colour = "grey55",
      fill = NA,
      linewidth = 0.6
    ),
    axis.line = element_blank(),
    axis.title = element_text(face = "bold", size = 11),
    axis.text.y = element_text(colour = "black", size = 9),
    axis.text.x = element_text(colour = "black", size = 8),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8)
  ) +
  guides(
    fill = guide_legend(
      nrow = 2,
      ncol = 3,
      byrow = TRUE,
      title.position = "top"
    )
  )


# ==============================================================================
# 8. Figure 4b: within-network, within-shape dissimilarity
# ==============================================================================
#
# Each point is one network x flower-shape combination. The y-axis measures the
# mean pairwise Sorensen dissimilarity among plant species sharing that flower
# shape within that network. Higher values mean greater within-shape variation.

within_shape_beta <- within_shape_beta %>%
  mutate(
    FlowerShape = factor(
      as.character(FlowerShape),
      levels = shape_order
    )
  ) %>%
  filter(!is.na(FlowerShape))

within_shape_emm_df <- within_shape_emm_df %>%
  mutate(
    FlowerShape = factor(
      as.character(FlowerShape),
      levels = shape_order
    )
  ) %>%
  filter(!is.na(FlowerShape))

p_within_shape <- ggplot(
  within_shape_beta,
  aes(
    x = FlowerShape,
    y = Mean_Sorensen
  )
) +
  geom_boxplot(
    fill = "#7B95C6",
    alpha = 0.45,
    width = 0.65,
    outlier.shape = NA,
    linewidth = 0.35
  ) +
  geom_jitter(
    aes(colour = N_plants),
    height = 0, # Jitter categories only; keep dissimilarity values within their bounds.
    width = 0.12,
    size = 1.2,
    alpha = 0.40
  ) +
  scale_colour_viridis_c(
    option = "viridis",
    direction = -1,
    name = "Within-network richness"
  ) +
  geom_errorbar(
    data = within_shape_emm_df,
    aes(
      x = FlowerShape,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.45
  ) +
  geom_point(
    data = within_shape_emm_df,
    aes(
      x = FlowerShape,
      y = emmean
    ),
    inherit.aes = FALSE,
    size = 2.3,
    shape = 21,
    fill = "white",
    stroke = 0.8
  ) +
  geom_text(
    data = within_shape_emm_df,
    aes(
      x = FlowerShape,
      y = 1.07,
      label = Letters
    ),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold",
    hjust = 0.5,
    colour = "black"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(0, 1.12),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = NULL,
    y = "Mean within-network dissimilarity"
  ) +
  theme_classic(base_size = 10) +
  theme(
    panel.border = element_rect(
      colour = "grey55",
      fill = NA,
      linewidth = 0.6
    ),
    axis.line = element_blank(),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.text.x = element_text(colour = "black", size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )


# ==============================================================================
# 9. Combine and save Figure 4
# ==============================================================================

get_bottom_legend <- function(p) {
  g <- ggplot2::ggplotGrob(p + theme(legend.position = "bottom"))
  ids <- grep("^guide-box", g$layout$name)
  ids <- ids[!vapply(g$grobs[ids], inherits, logical(1), "zeroGrob")]
  if (!length(ids)) stop("No visible legend was generated.")
  g$grobs[[ids[1]]]
}

legend_functional <- get_bottom_legend(visit_group)
legend_richness <- get_bottom_legend(p_within_shape)

panel_a <- visit_group + theme(legend.position = "none")
panel_b <- p_within_shape + theme(legend.position = "none")

# Add panel labels in a separate space above each plotting panel.
# draw_plot(height = 0.96) reserves the upper 4% for the label.
panel_a_labeled <- cowplot::ggdraw() +
  cowplot::draw_plot(
    panel_a,
    x = 0,
    y = 0,
    width = 1,
    height = 0.96
  ) +
  cowplot::draw_label(
    "(a)",
    x = 0.005,
    y = 1,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 13
  )

panel_b_labeled <- cowplot::ggdraw() +
  cowplot::draw_plot(
    panel_b,
    x = 0,
    y = 0,
    width = 1,
    height = 0.96
  ) +
  cowplot::draw_label(
    "(b)",
    x = 0.005,
    y = 1,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 13
  )

main_panels <- cowplot::plot_grid(
  panel_a_labeled,
  panel_b_labeled,
  nrow = 1,
  rel_widths = c(1.55, 1),
  align = "h",
  axis = "tb"
)

combined_legends <- cowplot::plot_grid(
  legend_functional,
  legend_richness,
  nrow = 1,
  rel_widths = c(1.7, 1)
)

figure_4 <- cowplot::plot_grid(
  main_panels,
  combined_legends,
  ncol = 1,
  rel_heights = c(1, 0.22)
)

print(figure_4)

ggsave(
  file.path(result_dir, "Figure_4_flower_shape_REVISED.png"),
  figure_4,
  width = 7.5,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(result_dir, "Figure_4_flower_shape_REVISED.pdf"),
  figure_4,
  width = 9,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)


