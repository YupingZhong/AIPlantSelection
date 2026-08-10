###############################################################################
# 4. Combine network metric correlation results
###############################################################################

# ==========================================================
# Load libraries
# ==========================================================

library(dplyr)
library(tibble)

# ==========================================================
# Load data
# ==========================================================

nestedness_cor <- readRDS(
  "data/processed/cor_results_df.rds")

connectance_cor <- readRDS(
  "data/processed/cor_results_df2.rds")

# ==========================================================
# Combine correlation results
# ==========================================================

cor_results_all <- bind_rows(
  nestedness_cor %>%
    mutate(Metric = "Nestedness"),
  
  connectance_cor %>%
    mutate(Metric = "Connectance")
)

cor_results_all


###############################################################################
# 3.3 Network metrics for random plant selection
# Note: randomisation may take several minutes
###############################################################################

# ==========================================================
# Load libraries
# ==========================================================

library(dplyr)
library(tidyr)
library(pbapply)
library(ggplot2)
library(data.table)
library(bipartite)

# ==========================================================
# Reproducibility
# ==========================================================

set.seed(2025)

# ==========================================================
# Prepare plant pool
# ==========================================================

plant_pool <- data_count_scaled %>%
  distinct(
    Study_Network_id,
    Plant_species,
    Flower_count_scaled,
    .keep_all = TRUE
  )

plant_pool_split <- split(
  plant_pool,
  plant_pool$Study_Network_id
)

network_list <- unique(
  plant_pool$Study_Network_id
)

# ==========================================================
# Function 1: randomly select plant species
# ==========================================================

get_random_interactions <- function(n_sp) {
  
  sampled_plants <- do.call(
    rbind,
    lapply(
      plant_pool_split,
      function(df) {
        
        n_select <- min(
          n_sp,
          nrow(df)
        )
        
        df[
          sample(
            nrow(df),
            n_select
          ),
        ]
      }
    )
  )
  
  data_interact %>%
    inner_join(
      sampled_plants,
      by = c(
        "Study_Network_id",
        "Plant_original_name" = "Plant_species"
      )
    )
}

# ==========================================================
# Function 2: calculate NODF
# ==========================================================

calc_nodf <- function(df) {
  
  dt <- as.data.table(df)
  
  mat_df <- dcast(
    dt,
    Plant_accepted_name ~ Pollinator_accepted_name,
    value.var = "Interaction_addup",
    fill = 0
  )
  
  mat <- as.matrix(
    mat_df[, -1, with = FALSE]
  )
  
  mat[is.na(mat)] <- 0
  
  # Remove empty rows and columns
  mat <- mat[
    rowSums(mat) > 0,
    ,
    drop = FALSE
  ]
  
  mat <- mat[
    ,
    colSums(mat) > 0,
    drop = FALSE
  ]
  
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    return(NA_real_)
  }
  
  tryCatch(
    bipartite::nested(
      mat,
      method = "NODF"
    ),
    error = function(e) NA_real_
  )
}

# ==========================================================
# Function 3: calculate connectance
# ==========================================================

calc_connectance <- function(df) {
  
  dt <- as.data.table(df)
  
  mat_df <- dcast(
    dt,
    Plant_accepted_name ~ Pollinator_accepted_name,
    value.var = "Interaction_addup",
    fill = 0
  )
  
  mat <- as.matrix(
    mat_df[, -1, with = FALSE]
  )
  
  mat[is.na(mat)] <- 0
  
  # Remove empty rows and columns
  mat <- mat[
    rowSums(mat) > 0,
    ,
    drop = FALSE
  ]
  
  mat <- mat[
    ,
    colSums(mat) > 0,
    drop = FALSE
  ]
  
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    return(NA_real_)
  }
  
  sum(mat > 0) /
    (nrow(mat) * ncol(mat))
}

# ==========================================================
# Function 4: run random network metrics
# ==========================================================

run_random_metrics <- function(
    n_sp,
    n_iter = 1000
) {
  
  pboptions(
    type = "timer"
  )
  
  res <- pblapply(
    seq_len(n_iter),
    function(i) {
      
      random_df <- get_random_interactions(
        n_sp
      )
      
      random_df %>%
        group_by(
          Study_Network_id
        ) %>%
        group_modify(
          ~ tibble(
            NODF = calc_nodf(.x),
            Connectance = calc_connectance(.x)
          )
        ) %>%
        ungroup()
    }
  )
  
  bind_rows(res)
}

# ==========================================================
# Run random subsampling
# ==========================================================

random_10_all <- run_random_metrics(
  10,
  1000
)

random_5_all <- run_random_metrics(
  5,
  1000
)

random_3_all <- run_random_metrics(
  3,
  1000
)

# ==========================================================
# Summarise random metrics
# ==========================================================

random_10 <- random_10_all %>%
  group_by(
    Study_Network_id
  ) %>%
  summarise(
    Random_NODF = mean(
      NODF,
      na.rm = TRUE
    ),
    Random_Connectance = mean(
      Connectance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Method = "random 10"
  )

random_5 <- random_5_all %>%
  group_by(
    Study_Network_id
  ) %>%
  summarise(
    Random_NODF = mean(
      NODF,
      na.rm = TRUE
    ),
    Random_Connectance = mean(
      Connectance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Method = "random 5"
  )

random_3 <- random_3_all %>%
  group_by(
    Study_Network_id
  ) %>%
  summarise(
    Random_NODF = mean(
      NODF,
      na.rm = TRUE
    ),
    Random_Connectance = mean(
      Connectance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Method = "random 3"
  )

# ==========================================================
# Combine random results
# ==========================================================

random_all <- bind_rows(
  random_10,
  random_5,
  random_3
)

# ==========================================================
# Save random network metrics
# ==========================================================

saveRDS(
  random_all,
  "data/processed/random_all_metrics.rds"
)

