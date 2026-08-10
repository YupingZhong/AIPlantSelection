
############################################################################
# 3. Network Metrics 
# calculation of connectance and nestedness for complete networks and sub-networks
###################################################################################

# ==========================================================
# Required input files
# ==========================================================

# - data/processed/data_count_scaled_published.rds
# - data/processed/data_interact_published.rds
# - data/processed/Subsampled_networks.rds

# ==========================================================
# Load libraries
# ==========================================================
library(bipartite)
library(vegan)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)

# ==========================================================
# Load data
# ==========================================================
data_count_scaled <- readRDS("data/processed/data_count_scaled_published.rds")
data_interact <- readRDS("data/processed/data_interact_published.rds")
datasets <- readRDS("data/processed/Subsampled_networks.rds" )

# ==========================================================
# Network IDs
# ==========================================================

result_orig <- data_interact%>%
  filter(!Interaction_addup == 0)%>%
  group_by(Study_Network_id) %>%
  ungroup()

# ==========================================================
# Function 1: calculate nestedness for one network
# ==========================================================
calculate_nestedness <- function(df) {
  
  dt <- as.data.table(df)
  
  interaction_matrix <- dcast(
    dt,
    Plant_accepted_name ~ Pollinator_accepted_name,
    value.var = "Interaction_addup",
    fill = 0
  )
  
  mat <- as.matrix(interaction_matrix[, -1, with = FALSE])
  mat[is.na(mat)] <- 0
  
  # 去掉空行列
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NA_real_)
  
  tryCatch(
    bipartite::nested(mat, method = "NODF"),
    error = function(e) NA_real_
  )
}

################
# ID skeleton

all_ids <- data.frame(
  Study_Network_id = sort(
    unique(
      data_interact %>%
        filter(Interaction_addup != 0) %>%
        pull(Study_Network_id))))

# ==========================================================
# Function 2: calculate nestedness for all networks
# ==========================================================

nested_df_from_data <- function(data, colname, id_skeleton) {
  
  df <- data %>%
    filter(
      !is.na(Pollinator_accepted_name),
      !is.na(Plant_accepted_name),
      !is.na(Interaction_addup)
    ) %>%
    dplyr::select(
      Plant_accepted_name,
      Pollinator_accepted_name,
      Interaction_addup,
      Study_Network_id
    ) %>%
    group_by(Study_Network_id) %>%
    summarise(
      value = calculate_nestedness(cur_data()),
      .groups = "drop"
    )
  
  # align skeleton（关键步骤）
  df_aligned <- id_skeleton %>%
    left_join(df, by = "Study_Network_id") %>%
    rename(!!colname := value)
  
  return(df_aligned)
}

# ==========================================================
# Calculate nestedness for all datasets
# ==========================================================

nested_list <- imap(
  datasets,
  ~ nested_df_from_data(.x, .y, all_ids))

all_nestedness_df <- reduce(nested_list, left_join, by = "Study_Network_id") %>%
  rename(
    Nestedness_orig = orig,
    Nestedness_Abun10 = Abun10,
    Nestedness_Abun5 = Abun5,
    Nestedness_Abun3 = Abun3,
    Nestedness_FlwShape5 = FlwShape5,
    Nestedness_FlwShape3 = FlwShape3,
    Nestedness_Pylo10 = Pylo10,
    Nestedness_Pylo5 = Pylo5,
    Nestedness_Pylo3 = Pylo3
  )%>%
  mutate(across(
    starts_with("Nestedness"),
    ~ as.numeric(unname(.))))

# ==========================================================
# Quality checks
# ==========================================================

# Check that network IDs are aligned across all datasets
id_check <- sapply(
  nested_list,
  function(df) identical(df$Study_Network_id, all_ids$Study_Network_id))
stopifnot(all(id_check))

# Check for duplicated network IDs
stopifnot(!anyDuplicated(all_nestedness_df$Study_Network_id))

# Summarise missing values
colSums(is.na(all_nestedness_df))

# =========================
# Correlation between subsampled and complete-network NODF
# =========================

nodf_cols <- c(
  "Nestedness_orig",
  "Nestedness_Abun10",
  "Nestedness_Abun5",
  "Nestedness_Abun3",
  "Nestedness_FlwShape5",
  "Nestedness_FlwShape3",
  "Nestedness_Pylo10",
  "Nestedness_Pylo5",
  "Nestedness_Pylo3"
)

# ==========================================================
# Calculate Spearman correlations
# ==========================================================

cor_results <- lapply(nodf_cols[-1], function(x) {
  
  x_vec <- as.numeric(all_nestedness_df[[x]])
  y_vec <- as.numeric(all_nestedness_df$Nestedness_orig)
  
  # exclude NA
  valid_idx <- !is.na(x_vec) & !is.na(y_vec)
  
  x_clean <- x_vec[valid_idx]
  y_clean <- y_vec[valid_idx]
  
  if(length(x_clean) < 3 || length(unique(x_clean)) < 2){
    return(c(r = NA, p_value = NA, n = length(x_clean)))
  }
  
  test <- cor.test(x_clean, y_clean, method = "pearson")
  
  c(r = unname(test$estimate), p_value = test$p.value, n = length(x_clean))
})


# ==========================================================
# Combine results
# ==========================================================

cor_results_df <- as.data.frame(do.call(rbind, cor_results)) %>%
  mutate(
    Comparison = nodf_cols[-1],
    Metric = "Nestedness",
    r = as.numeric(r),
    p_value = as.numeric(p_value)
  )

rownames(cor_results_df) <- nodf_cols[-1]
# 
# print("Pearson correlation with full network:")
print(cor_results_df)

cor(
  all_nestedness_df$Nestedness_FlwShape5,
  all_nestedness_df$Nestedness_orig,
  use = "complete.obs")

# =========================
# Save results
# =========================

saveRDS(all_nestedness_df, "data/processed/all_nestedness_df.rds")
saveRDS(cor_results_df, "data/processed/cor_results_df.rds")

