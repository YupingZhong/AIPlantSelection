############################################################################
# 3. Network Metrics 
# calculation of H2 for complete networks and sub-networks
###################################################################################
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

################
# ID skeleton

all_ids <- data.frame(
  Study_Network_id = sort(
    unique(
      data_interact %>%
        filter(Interaction_addup != 0) %>%
        pull(Study_Network_id))))
#### function
calculate_H2 <- function(df) {
  
  dt <- as.data.table(df)
  
  interaction_matrix <- dcast(
    dt,
    Plant_accepted_name ~ Pollinator_accepted_name,
    value.var = "Interaction_addup",
    fill = 0
  )
  
  mat <- as.matrix(interaction_matrix[, -1, with = FALSE])
  mat[is.na(mat)] <- 0
  
  # remove empty rows/columns
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  
  if (nrow(mat) < 2 || ncol(mat) < 2)
    return(NA_real_)
  
  tryCatch(
    bipartite::networklevel(
      mat,
      index = "H2"
    ),
    error = function(e) NA_real_
  )
}


H2_df_from_data <- function(data, colname, id_skeleton) {
  
  df <- data %>%
    filter(
      !is.na(Pollinator_accepted_name),
      !is.na(Plant_accepted_name),
      !is.na(Interaction_addup)
    ) %>%
    select(
      Plant_accepted_name,
      Pollinator_accepted_name,
      Interaction_addup,
      Study_Network_id
    ) %>%
    group_by(Study_Network_id) %>%
    summarise(
      value = calculate_H2(cur_data()),
      .groups = "drop"
    )
  
  df_aligned <- id_skeleton %>%
    left_join(df, by = "Study_Network_id") %>%
    rename(!!colname := value)
  
  return(df_aligned)
}

H2_list <- imap(
  datasets,
  ~ H2_df_from_data(.x, .y, all_ids)
)

all_H2_df <- reduce(
  H2_list,
  left_join,
  by = "Study_Network_id"
) %>%
  rename(
    H2_orig = orig,
    H2_Abun10 = Abun10,
    H2_Abun5 = Abun5,
    H2_Abun3 = Abun3,
    H2_FlwShape5 = FlwShape5,
    H2_FlwShape3 = FlwShape3,
    H2_Pylo10 = Pylo10,
    H2_Pylo5 = Pylo5,
    H2_Pylo3 = Pylo3
  ) %>%
  mutate(
    across(
      starts_with("H2"),
      ~ as.numeric(unname(.))
    )
  )


summary(all_H2_df)
colSums(is.na(all_H2_df))

H2_cols <- c(
  "H2_orig",
  "H2_Abun10",
  "H2_Abun5",
  "H2_Abun3",
  "H2_FlwShape5",
  "H2_FlwShape3",
  "H2_Pylo10",
  "H2_Pylo5",
  "H2_Pylo3"
)

cor_results_H2 <- lapply(H2_cols[-1], function(x){
  
  x_vec <- as.numeric(all_H2_df[[x]])
  y_vec <- as.numeric(all_H2_df$H2_orig)
  
  valid_idx <- !is.na(x_vec) & !is.na(y_vec)
  
  x_clean <- x_vec[valid_idx]
  y_clean <- y_vec[valid_idx]
  
  if(length(x_clean) < 3 ||
     length(unique(x_clean)) < 2){
    
    return(
      c(
        rho = NA,
        p_value = NA,
        n = length(x_clean)
      )
    )
  }
  
  test <- cor.test(
    x_clean,
    y_clean,
    method = "spearman",
    exact = FALSE
  )
  
  c(
    rho = unname(test$estimate),
    p_value = test$p.value,
    n = length(x_clean)
  )
})

cor_results_H2_df <- as.data.frame(
  do.call(rbind, cor_results_H2)
) %>%
  mutate(
    Comparison = H2_cols[-1],
    Metric = "H2",
    rho = as.numeric(rho),
    p_value = as.numeric(p_value)
  )

rownames(cor_results_H2_df) <- H2_cols[-1]

saveRDS(
  all_H2_df,
  "data/processed/all_H2_df.rds"
)

saveRDS(
  cor_results_H2_df,
  "data/processed/cor_results_H2_df.rds"
)
