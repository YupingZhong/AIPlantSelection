################################################################################
# ==============================================================================
# HYPOTHESIS 2
# Favored plant species by pollinators are rare
#
# ANALYSIS SCRIPT ONLY
#
# This script:
#   1. Prepares plant-level abundance / interaction / richness data
#   2. Fits the mixed model:
#        Proportion_of_richness ~ Abundance_scaled +
#        (1 | Study_id / Study_Network_id)
#   3. Calculates marginal and conditional R2
#   4. Calculates predicted richness capture and residual capture
#   5. Identifies rare attractive plants
#   6. Determines network-level Present / Absent groups
#   7. Performs Wilcoxon test
#   8. Performs sensitivity analysis
#   9. Saves all analysis results for later plotting
#
# IMPORTANT:
#   This script does NOT create the final three-panel figure.
# ==============================================================================

################################################################################
##################################################

library(dplyr)
library(lme4)
library(lmerTest)
library(purrr)

# ==============================================================================
# 1. Paths
# ==============================================================================

data_dir <- "data/processed"
result_dir <- "/Chap1_TargetPlant_to_monitor/result_260723"

if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)


# ==============================================================================
# 2. Load data
# ==============================================================================

data_count_scaled <- readRDS(
  file.path(data_dir, "data_count_scaled_published.rds")
)

data_interact <- readRDS(
  file.path(data_dir, "data_interact_published.rds")
)

data_merge <- readRDS(
  file.path(data_dir, "data_merge.rds")
)

percent_10 <- readRDS(
  file.path(data_dir, "percent_10.rds")
)

traits <- read.csv(
  file.path(data_dir, "merge.trait.csv"),
  header = TRUE,
  fileEncoding = "UTF-8"
)


# ==============================================================================
# 3. Top 10 abundant plants in each network
# ==============================================================================

top_10_species <- data_count_scaled %>%
  distinct(Study_Network_id, Plant_species, Flower_count_scaled) %>%
  group_by(Study_Network_id) %>%
  slice_max(Flower_count_scaled, n = 10, with_ties = FALSE) %>%
  filter(n_distinct(Plant_species) > 9) %>%
  ungroup() %>%
  select(Study_Network_id, Plant_species)


# ==============================================================================
# 4. Network-level total interaction and pollinator richness
# ==============================================================================

interaction_all <- data_interact %>%
  filter(Interaction_addup != 0) %>%
  group_by(Study_Network_id) %>%
  summarise(
    Interaction_sum = sum(Interaction_addup, na.rm = TRUE),
    richness_sum = n_distinct(Pollinator_accepted_name),
    .groups = "drop"
  )


# ==============================================================================
# 5. Plant name mapping
# ==============================================================================

plant_name_map <- data_interact %>%
  distinct(Plant_accepted_name, Plant_original_name) %>%
  filter(!is.na(Plant_accepted_name), !is.na(Plant_original_name)) %>%
  rename(Plant_species = Plant_original_name)


# ==============================================================================
# 6. Plant-level analysis data
# ==============================================================================

data_combined_all <- data_merge %>%
  left_join(interaction_all, by = "Study_Network_id") %>%
  group_by(Study_id, Study_Network_id, Plant_species) %>%
  summarise(
    Plant_accepted_name = first(Plant_accepted_name),
    Abundance = first(Flower_count_scaled),
    Pollinator_richness = n_distinct(Pollinator_accepted_name, na.rm = TRUE),
    Interaction_times = sum(Interaction_addup, na.rm = TRUE),
    Proportion_of_interactions =
      Interaction_times / first(Interaction_sum) * 100,
    Proportion_of_richness =
      Pollinator_richness / first(richness_sum) * 100,
    .groups = "drop"
  ) %>%
  filter(
    is.finite(Abundance),
    is.finite(Interaction_times),
    is.finite(Proportion_of_richness)
  ) %>%
  group_by(Study_Network_id) %>%
  mutate(
    Abundance_scaled = scale(log(Abundance + 1))[, 1]
  ) %>%
  ungroup() %>%
  left_join(plant_name_map, by = "Plant_species",
            suffix = c("", "_map")) %>%
  mutate(
    Plant_accepted_name = coalesce(
      Plant_accepted_name,
      Plant_accepted_name_map
    )
  ) %>%
  select(-Plant_accepted_name_map)


# ==============================================================================
# 7. Mixed model
# ==============================================================================

model_h2 <- lmer(
  Proportion_of_richness ~ Abundance_scaled +
    (1 | Study_id / Study_Network_id),
  data = data_combined_all
)


# ==============================================================================
# 8. Model statistics
# ==============================================================================

coef_tab <- as.data.frame(summary(model_h2)$coefficients)
coef_tab$term <- rownames(coef_tab)

beta_abundance <- coef_tab$Estimate[
  coef_tab$term == "Abundance_scaled"
]

se_abundance <- coef_tab$`Std. Error`[
  coef_tab$term == "Abundance_scaled"
]

t_abundance <- coef_tab$`t value`[
  coef_tab$term == "Abundance_scaled"
]

p_abundance <- coef_tab$`Pr(>|t|)`[
  coef_tab$term == "Abundance_scaled"
]


# ==============================================================================
# 9. R2
# ==============================================================================

fixed_pred <- predict(model_h2, re.form = NA)

var_fixed <- var(fixed_pred, na.rm = TRUE)
var_rand <- sum(as.data.frame(VarCorr(model_h2))$vcov)
var_resid <- attr(VarCorr(model_h2), "sc")^2

r2_marginal <- var_fixed / (var_fixed + var_rand + var_resid)

r2_conditional <- (var_fixed + var_rand) /
  (var_fixed + var_rand + var_resid)


model_results <- tibble(
  beta = beta_abundance,
  SE = se_abundance,
  t = t_abundance,
  p = p_abundance,
  Marginal_R2 = r2_marginal,
  Conditional_R2 = r2_conditional
)

print(model_results)


# ==============================================================================
# 10. Predicted richness capture and residual
# ==============================================================================

data_combined_all <- data_combined_all %>%
  mutate(
    Predicted_capture = predict(
      model_h2,
      newdata = .,
      re.form = NA,
      allow.new.levels = TRUE
    ),
    Residual_capture =
      Proportion_of_richness - Predicted_capture
  )


# ==============================================================================
# 11. Rare attractive plants
# ==============================================================================

get_rare_attr <- function(data, rare_q = 0.05, attr_q = 0.90) {
  
  rare_thr <- quantile(
    data$Abundance_scaled,
    rare_q,
    na.rm = TRUE
  )
  
  attr_thr <- quantile(
    data$Residual_capture,
    attr_q,
    na.rm = TRUE
  )
  
  data %>%
    filter(
      Abundance_scaled < rare_thr,
      Residual_capture > attr_thr
    )
}

Rare_favour_plants <- get_rare_attr(
  data_combined_all,
  0.05,
  0.90
)

rare_plant_list <- Rare_favour_plants %>%
  distinct(Plant_species, Plant_accepted_name)


rare_species_accepted <- rare_plant_list %>%
  filter(!is.na(Plant_accepted_name)) %>%
  distinct(Plant_accepted_name) %>%
  pull(Plant_accepted_name)


# ==============================================================================
# 12. Network-level Present / Absent classification
# ==============================================================================

network_highattr <- data_interact %>%
  distinct(Study_Network_id, Plant_accepted_name) %>%
  mutate(
    is_highattr =
      Plant_accepted_name %in% rare_species_accepted
  )

top10_tbl <- top_10_species %>%
  left_join(plant_name_map, by = "Plant_species") %>%
  filter(!is.na(Plant_accepted_name)) %>%
  group_by(Study_Network_id) %>%
  summarise(
    top10_list = list(unique(Plant_accepted_name)),
    .groups = "drop"
  )

network_group <- network_highattr %>%
  left_join(top10_tbl, by = "Study_Network_id") %>%
  group_split(Study_Network_id) %>%
  map_dfr(function(df) {
    
    highattr <- df$Plant_accepted_name[df$is_highattr]
    top10 <- unique(unlist(df$top10_list[[1]]))
    
    tibble(
      Study_Network_id = unique(df$Study_Network_id),
      present = any(df$is_highattr),
      fully_covered =
        length(highattr) > 0 &&
        all(highattr %in% top10)
    )
    
  }) %>%
  mutate(
    group = ifelse(
      present & !fully_covered,
      "Present",
      "Absent"
    )
  )


# ==============================================================================
# 13. Wilcoxon test
# ==============================================================================

attract_percent_10 <- percent_10 %>%
  rename(percentage_Abun10 = percentage) %>%
  left_join(
    network_group %>% select(Study_Network_id, group),
    by = "Study_Network_id"
  ) %>%
  rename(present_group = group) %>%
  filter(
    !is.na(percentage_Abun10),
    !is.na(present_group)
  )

attract_percent_10$present_group <- factor(
  attract_percent_10$present_group,
  levels = c("Present", "Absent")
)

wilcox_res <- wilcox.test(
  percentage_Abun10 ~ present_group,
  data = attract_percent_10
)

print(wilcox_res)


# ==============================================================================
# 14. Sensitivity analysis
# ==============================================================================

run_present_absent_test <- function(data, rare_q, attr_q) {
  
  rare_thr <- quantile(
    data$Abundance_scaled,
    rare_q,
    na.rm = TRUE
  )
  
  attr_thr <- quantile(
    data$Residual_capture,
    attr_q,
    na.rm = TRUE
  )
  
  rare_species <- data %>%
    filter(
      Abundance_scaled < rare_thr,
      Residual_capture > attr_thr
    ) %>%
    distinct(Plant_species) %>%
    left_join(plant_name_map, by = "Plant_species") %>%
    filter(!is.na(Plant_accepted_name)) %>%
    distinct(Plant_accepted_name) %>%
    pull(Plant_accepted_name)
  
  network_highattr <- data_interact %>%
    distinct(Study_Network_id, Plant_accepted_name) %>%
    mutate(
      is_highattr =
        Plant_accepted_name %in% rare_species
    )
  
  network_group <- network_highattr %>%
    left_join(top10_tbl, by = "Study_Network_id") %>%
    group_split(Study_Network_id) %>%
    map_dfr(function(df) {
      
      highattr <- df$Plant_accepted_name[df$is_highattr]
      top10 <- unique(unlist(df$top10_list[[1]]))
      
      tibble(
        Study_Network_id = unique(df$Study_Network_id),
        present = any(df$is_highattr),
        fully_covered =
          length(highattr) > 0 &&
          all(highattr %in% top10)
      )
      
    }) %>%
    mutate(
      group = ifelse(
        present & !fully_covered,
        "Present",
        "Absent"
      )
    )
  
  df_test <- percent_10 %>%
    rename(percentage_Abun10 = percentage) %>%
    left_join(
      network_group %>% select(Study_Network_id, group),
      by = "Study_Network_id"
    ) %>%
    filter(
      !is.na(percentage_Abun10),
      !is.na(group)
    )
  
  wt <- wilcox.test(
    percentage_Abun10 ~ group,
    data = df_test
  )
  
  median_present <- median(
    df_test$percentage_Abun10[df_test$group == "Present"],
    na.rm = TRUE
  )
  
  median_absent <- median(
    df_test$percentage_Abun10[df_test$group == "Absent"],
    na.rm = TRUE
  )
  
  tibble(
    rare_q = rare_q,
    attr_q = attr_q,
    p_value = wt$p.value,
    effect_size = median_present - median_absent,
    n_species = length(rare_species)
  )
}


rare_cut <- c(0.05, 0.10, 0.15)
attr_cut <- c(0.90, 0.85, 0.80)

sensitivity_results <- expand.grid(
  rare_q = rare_cut,
  attr_q = attr_cut
) %>%
  split(seq_len(nrow(.))) %>%
  map_dfr(function(x) {
    run_present_absent_test(
      data_combined_all,
      x$rare_q,
      x$attr_q
    )
  })


print(sensitivity_results)


# ==============================================================================
# 15. Save results
# ==============================================================================

saveRDS(
  data_combined_all,
  file.path(data_dir, "H2_data_combined_all.rds")
)

saveRDS(
  model_h2,
  file.path(data_dir, "H2_model.rds")
)

saveRDS(
  model_results,
  file.path(data_dir, "H2_model_summary.rds")
)

saveRDS(
  Rare_favour_plants,
  file.path(data_dir, "H2_Rare_favour_plants.rds")
)

saveRDS(
  rare_species_accepted,
  file.path(data_dir, "H2_rare_species_accepted.rds")
)

saveRDS(
  top_10_species,
  file.path(data_dir, "H2_top_10_species.rds")
)

saveRDS(
  network_group,
  file.path(data_dir, "H2_network_group.rds")
)

saveRDS(
  attract_percent_10,
  file.path(data_dir, "H2_attract_percent_10.rds")
)

saveRDS(
  wilcox_res,
  file.path(data_dir, "H2_wilcox_res.rds")
)

saveRDS(
  sensitivity_results,
  file.path(data_dir, "H2_sensitivity_results.rds")
)

write.csv(
  model_results,
  file.path(result_dir, "H2_model_summary.csv"),
  row.names = FALSE
)

write.csv(
  sensitivity_results,
  file.path(result_dir, "H2_sensitivity_results.csv"),
  row.names = FALSE
)

write.csv(
  rare_plant_list,
  file.path(result_dir, "H2_rare_attractive_plants.csv"),
  row.names = FALSE
)

