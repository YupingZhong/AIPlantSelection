###############################################################################
# Sampling curves for NODF, Connectance and H2
#
# Purpose:
#   Evaluate how many plant species are required for the sub-network
#   to reproduce the structure of the complete network.
#
# IMPORTANT:
#   Full-network metrics are NOT recalculated here.
#   They are read from previously generated RDS files.
#
#   Sub-network metrics are calculated for each sampling effort n.
#   Spearman's rank correlation (rho) is then calculated between
#   sub-network and complete-network metrics.
###############################################################################


# =============================================================================
# 0. Load libraries
# =============================================================================

library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(bipartite)
library(minpack.lm)
library(segmented)
library(cowplot)


# =============================================================================
# 1. Load original data
# =============================================================================

data_count_scaled <- readRDS(
  "data/processed/data_count_scaled_published.rds"
)

data_interact <- readRDS(
  "data/processed/data_interact_published.rds"
)

traits <- read.csv(
  "data/processed/merge.trait.csv",
  header = TRUE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 2. Load previously calculated full-network metrics
#
# IMPORTANT:
#   These metrics were already calculated in independent scripts.
#   DO NOT recalculate them here.
# =============================================================================

all_nestedness_df <- readRDS(
  "data/processed/all_nestedness_df.rds"
)

all_connectance_df <- readRDS(
  "data/processed/all_connectance_df.rds"
)

all_H2_df <- readRDS(
  "data/processed/all_H2_df.rds"
)


# =============================================================================
# 3. Prepare full-network metric table
# =============================================================================

full_metrics <- all_nestedness_df %>%
  dplyr::select(
    Study_Network_id,
    NODF_full = Nestedness_orig
  ) %>%
  left_join(
    all_connectance_df %>%
      dplyr::select(
        Study_Network_id,
        Connectance_full = Connectance_orig
      ),
    by = "Study_Network_id"
  ) %>%
  left_join(
    all_H2_df %>%
      dplyr::select(
        Study_Network_id,
        H2_full = H2_orig
      ),
    by = "Study_Network_id"
  )


# =============================================================================
# 4. Quality check
# =============================================================================

cat("\n============================================\n")
cat("Full-network metric check\n")
cat("============================================\n")

print(summary(full_metrics))

cat("\nMissing values:\n")
print(colSums(is.na(full_metrics)))

cat("\nNumber of networks:\n")
print(nrow(full_metrics))


# =============================================================================
# 5. Prepare data
# =============================================================================

data_count_scaled <- data_count_scaled %>%
  mutate(
    Plant_species = trimws(Plant_species)
  )

data_interact <- data_interact %>%
  mutate(
    Plant_original_name = trimws(Plant_original_name),
    Interaction_addup = replace_na(
      Interaction_addup,
      0
    )
  )

networks <- unique(
  data_interact$Study_Network_id
)


# =============================================================================
# 6. Function: calculate NODF for a sub-network
# =============================================================================

calc_nodf <- function(df) {
  
  if (nrow(df) == 0)
    return(NA_real_)
  
  dt <- as.data.table(df)
  
  mat_df <- tryCatch(
    
    dcast(
      dt,
      Plant_accepted_name ~ Pollinator_accepted_name,
      value.var = "Interaction_addup",
      fun.aggregate = sum,
      fill = 0
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(mat_df))
    return(NA_real_)
  
  if (ncol(mat_df) < 2)
    return(NA_real_)
  
  mat <- as.matrix(
    mat_df[, -1, with = FALSE]
  )
  
  mat[is.na(mat)] <- 0
  
  # Remove empty rows
  mat <- mat[
    rowSums(mat) > 0,
    ,
    drop = FALSE
  ]
  
  # Remove empty columns
  mat <- mat[
    ,
    colSums(mat) > 0,
    drop = FALSE
  ]
  
  if (
    nrow(mat) < 2 ||
    ncol(mat) < 2 ||
    sum(mat) == 0
  ) {
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


# =============================================================================
# 7. Function: calculate Connectance for a sub-network
# =============================================================================

calc_connectance <- function(df) {
  
  if (nrow(df) == 0)
    return(NA_real_)
  
  dt <- as.data.table(df)
  
  mat_df <- tryCatch(
    
    dcast(
      dt,
      Plant_accepted_name ~ Pollinator_accepted_name,
      value.var = "Interaction_addup",
      fun.aggregate = sum,
      fill = 0
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(mat_df))
    return(NA_real_)
  
  if (ncol(mat_df) < 2)
    return(NA_real_)
  
  mat <- as.matrix(
    mat_df[, -1, with = FALSE]
  )
  
  mat[is.na(mat)] <- 0
  
  # Remove empty rows
  mat <- mat[
    rowSums(mat) > 0,
    ,
    drop = FALSE
  ]
  
  # Remove empty columns
  mat <- mat[
    ,
    colSums(mat) > 0,
    drop = FALSE
  ]
  
  if (
    nrow(mat) < 2 ||
    ncol(mat) < 2
  ) {
    return(NA_real_)
  }
  
  num_links <- sum(mat > 0)
  
  num_possible <- (
    nrow(mat) *
      ncol(mat)
  )
  
  if (num_possible == 0)
    return(NA_real_)
  
  num_links / num_possible
}


# =============================================================================
# 8. Function: calculate H2 for a sub-network
# =============================================================================

calc_H2 <- function(df) {
  
  if (nrow(df) == 0)
    return(NA_real_)
  
  dt <- as.data.table(df)
  
  mat_df <- tryCatch(
    
    dcast(
      dt,
      Plant_accepted_name ~ Pollinator_accepted_name,
      value.var = "Interaction_addup",
      fun.aggregate = sum,
      fill = 0
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(mat_df))
    return(NA_real_)
  
  if (ncol(mat_df) < 2)
    return(NA_real_)
  
  mat <- as.matrix(
    mat_df[, -1, with = FALSE]
  )
  
  mat[is.na(mat)] <- 0
  
  # Remove empty rows
  mat <- mat[
    rowSums(mat) > 0,
    ,
    drop = FALSE
  ]
  
  # Remove empty columns
  mat <- mat[
    ,
    colSums(mat) > 0,
    drop = FALSE
  ]
  
  if (
    nrow(mat) < 2 ||
    ncol(mat) < 2
  ) {
    return(NA_real_)
  }
  
  tryCatch(
    
    bipartite::networklevel(
      mat,
      index = "H2"
    ),
    
    error = function(e) NA_real_
  )
}


# =============================================================================
# 9. Helper function: Spearman's rho
# =============================================================================

calculate_spearman <- function(
    x,
    y,
    min_n = 3
) {
  
  valid <- (
    !is.na(x) &
      !is.na(y)
  )
  
  x <- x[valid]
  y <- y[valid]
  
  if (length(x) < min_n)
    return(NA_real_)
  
  if (length(unique(x)) < 2)
    return(NA_real_)
  
  if (length(unique(y)) < 2)
    return(NA_real_)
  
  suppressWarnings(
    cor(
      x,
      y,
      method = "spearman"
    )
  )
}


# =============================================================================
# 10. Abundance-based sampling curve
#
# Select the top n plant species according to Flower_count_scaled.
#
# For each n:
#   1. calculate sub-network NODF
#   2. calculate sub-network Connectance
#   3. calculate sub-network H2
#   4. compare with PRE-CALCULATED full-network metrics
#   5. calculate Spearman's rho
# =============================================================================

get_structure_curve <- function(
    max_n,
    data_count_scaled,
    data_interact,
    full_metrics,
    verbose = TRUE
) {
  
  results <- data.frame(
    n = integer(),
    rho_nodf = numeric(),
    rho_conn = numeric(),
    rho_H2 = numeric(),
    n_nodf = integer(),
    n_conn = integer(),
    n_H2 = integer()
  )
  
  networks <- unique(
    data_interact$Study_Network_id
  )
  
  for (i in 1:max_n) {
    
    nodf_sub <- c()
    nodf_full <- c()
    
    conn_sub <- c()
    conn_full <- c()
    
    H2_sub <- c()
    H2_full <- c()
    
    
    # -------------------------------------------------------------------------
    # Loop through networks
    # -------------------------------------------------------------------------
    
    for (net in networks) {
      
      # Get already calculated full-network metrics
      
      full_metric <- full_metrics %>%
        filter(
          Study_Network_id == net
        )
      
      if (nrow(full_metric) == 0)
        next
      
      
      # -----------------------------------------------------------------------
      # Select top-abundant plant species
      # -----------------------------------------------------------------------
      
      top_species <- data_count_scaled %>%
        filter(
          Study_Network_id == net
        ) %>%
        arrange(
          desc(Flower_count_scaled)
        ) %>%
        slice_head(
          n = i
        ) %>%
        pull(
          Plant_species
        )
      
      if (length(top_species) < i)
        next
      
      
      # -----------------------------------------------------------------------
      # Construct sub-network
      # -----------------------------------------------------------------------
      
      sub_df <- data_interact %>%
        filter(
          Study_Network_id == net,
          Plant_original_name %in% top_species
        )
      
      if (nrow(sub_df) == 0)
        next
      
      
      # -----------------------------------------------------------------------
      # Calculate SUBNETWORK metrics
      #
      # These are the metrics that need to be recalculated for each n.
      # -----------------------------------------------------------------------
      
      sub_nodf <- calc_nodf(
        sub_df
      )
      
      sub_conn <- calc_connectance(
        sub_df
      )
      
      sub_H2 <- calc_H2(
        sub_df
      )
      
      
      # -----------------------------------------------------------------------
      # Pair sub-network and full-network NODF
      # -----------------------------------------------------------------------
      
      if (
        !is.na(sub_nodf) &&
        !is.na(full_metric$NODF_full)
      ) {
        
        nodf_sub <- c(
          nodf_sub,
          sub_nodf
        )
        
        nodf_full <- c(
          nodf_full,
          full_metric$NODF_full
        )
      }
      
      
      # -----------------------------------------------------------------------
      # Pair sub-network and full-network Connectance
      # -----------------------------------------------------------------------
      
      if (
        !is.na(sub_conn) &&
        !is.na(full_metric$Connectance_full)
      ) {
        
        conn_sub <- c(
          conn_sub,
          sub_conn
        )
        
        conn_full <- c(
          conn_full,
          full_metric$Connectance_full
        )
      }
      
      
      # -----------------------------------------------------------------------
      # Pair sub-network and full-network H2
      # -----------------------------------------------------------------------
      
      if (
        !is.na(sub_H2) &&
        !is.na(full_metric$H2_full)
      ) {
        
        H2_sub <- c(
          H2_sub,
          sub_H2
        )
        
        H2_full <- c(
          H2_full,
          full_metric$H2_full
        )
      }
    }
    
    
    # -------------------------------------------------------------------------
    # Spearman's rank correlations
    # -------------------------------------------------------------------------
    
    rho_nodf <- calculate_spearman(
      nodf_sub,
      nodf_full
    )
    
    rho_conn <- calculate_spearman(
      conn_sub,
      conn_full
    )
    
    rho_H2 <- calculate_spearman(
      H2_sub,
      H2_full
    )
    
    
    # -------------------------------------------------------------------------
    # Save results
    # -------------------------------------------------------------------------
    
    results <- rbind(
      
      results,
      
      data.frame(
        n = i,
        
        rho_nodf = rho_nodf,
        
        rho_conn = rho_conn,
        
        rho_H2 = rho_H2,
        
        n_nodf = length(nodf_sub),
        
        n_conn = length(conn_sub),
        
        n_H2 = length(H2_sub)
      )
    )
    
    
    # -------------------------------------------------------------------------
    # Progress
    # -------------------------------------------------------------------------
    
    if (verbose) {
      
      cat(
        "n =", i,
        "| NODF n =", length(nodf_sub),
        "| Connectance n =", length(conn_sub),
        "| H2 n =", length(H2_sub),
        "\n"
      )
    }
  }
  
  return(results)
}


# =============================================================================
# 11. Random sampling curve
#
# For every n:
#   randomly select n plant species within each network
#
# Repeat n_rep times.
#
# For each repetition:
#   calculate Spearman's rho between sub-network and full-network metrics.
#
# Final output:
#   mean rho
#   SD of rho
# =============================================================================

get_structure_curve_random <- function(
    max_n,
    data_count_scaled,
    data_interact,
    full_metrics,
    n_rep = 50,
    verbose = TRUE
) {
  
  results <- data.frame(
    n = integer(),
    
    rho_nodf_mean = numeric(),
    rho_nodf_sd = numeric(),
    
    rho_conn_mean = numeric(),
    rho_conn_sd = numeric(),
    
    rho_H2_mean = numeric(),
    rho_H2_sd = numeric()
  )
  
  
  networks <- unique(
    data_interact$Study_Network_id
  )
  
  
  for (i in 1:max_n) {
    
    rho_nodf_rep <- c()
    rho_conn_rep <- c()
    rho_H2_rep <- c()
    
    
    # ========================================================================
    # Repeat random sampling
    # ========================================================================
    
    for (rep in 1:n_rep) {
      
      nodf_sub <- c()
      nodf_full <- c()
      
      conn_sub <- c()
      conn_full <- c()
      
      H2_sub <- c()
      H2_full <- c()
      
      
      # ======================================================================
      # Loop through networks
      # ======================================================================
      
      for (net in networks) {
        
        # Get full-network metrics already calculated
        
        full_metric <- full_metrics %>%
          filter(
            Study_Network_id == net
          )
        
        if (nrow(full_metric) == 0)
          next
        
        
        # --------------------------------------------------------------------
        # Get plant species pool
        # --------------------------------------------------------------------
        
        all_species <- data_count_scaled %>%
          filter(
            Study_Network_id == net
          ) %>%
          distinct(
            Plant_species
          ) %>%
          pull(
            Plant_species
          )
        
        
        if (length(all_species) < i)
          next
        
        
        # --------------------------------------------------------------------
        # Randomly select i plant species
        # --------------------------------------------------------------------
        
        random_species <- sample(
          all_species,
          size = i,
          replace = FALSE
        )
        
        
        # --------------------------------------------------------------------
        # Construct sub-network
        # --------------------------------------------------------------------
        
        sub_df <- data_interact %>%
          filter(
            Study_Network_id == net,
            Plant_original_name %in% random_species
          )
        
        
        if (nrow(sub_df) == 0)
          next
        
        
        # --------------------------------------------------------------------
        # Calculate SUBNETWORK metrics
        # --------------------------------------------------------------------
        
        sub_nodf <- calc_nodf(
          sub_df
        )
        
        sub_conn <- calc_connectance(
          sub_df
        )
        
        sub_H2 <- calc_H2(
          sub_df
        )
        
        
        # --------------------------------------------------------------------
        # Pair NODF
        # --------------------------------------------------------------------
        
        if (
          !is.na(sub_nodf) &&
          !is.na(full_metric$NODF_full)
        ) {
          
          nodf_sub <- c(
            nodf_sub,
            sub_nodf
          )
          
          nodf_full <- c(
            nodf_full,
            full_metric$NODF_full
          )
        }
        
        
        # --------------------------------------------------------------------
        # Pair Connectance
        # --------------------------------------------------------------------
        
        if (
          !is.na(sub_conn) &&
          !is.na(full_metric$Connectance_full)
        ) {
          
          conn_sub <- c(
            conn_sub,
            sub_conn
          )
          
          conn_full <- c(
            conn_full,
            full_metric$Connectance_full
          )
        }
        
        
        # --------------------------------------------------------------------
        # Pair H2
        # --------------------------------------------------------------------
        
        if (
          !is.na(sub_H2) &&
          !is.na(full_metric$H2_full)
        ) {
          
          H2_sub <- c(
            H2_sub,
            sub_H2
          )
          
          H2_full <- c(
            H2_full,
            full_metric$H2_full
          )
        }
      }
      
      
      # ======================================================================
      # Calculate Spearman's rho for this repetition
      # ======================================================================
      
      rho_nodf <- calculate_spearman(
        nodf_sub,
        nodf_full
      )
      
      rho_conn <- calculate_spearman(
        conn_sub,
        conn_full
      )
      
      rho_H2 <- calculate_spearman(
        H2_sub,
        H2_full
      )
      
      
      # ----------------------------------------------------------------------
      # Save repetition results
      # ----------------------------------------------------------------------
      
      if (!is.na(rho_nodf))
        rho_nodf_rep <- c(
          rho_nodf_rep,
          rho_nodf
        )
      
      if (!is.na(rho_conn))
        rho_conn_rep <- c(
          rho_conn_rep,
          rho_conn
        )
      
      if (!is.na(rho_H2))
        rho_H2_rep <- c(
          rho_H2_rep,
          rho_H2
        )
    }
    
    
    # ========================================================================
    # Summarise random repetitions
    # ========================================================================
    
    results <- rbind(
      
      results,
      
      data.frame(
        
        n = i,
        
        rho_nodf_mean =
          ifelse(
            length(rho_nodf_rep) > 0,
            mean(
              rho_nodf_rep,
              na.rm = TRUE
            ),
            NA_real_
          ),
        
        rho_nodf_sd =
          ifelse(
            length(rho_nodf_rep) > 1,
            sd(
              rho_nodf_rep,
              na.rm = TRUE
            ),
            NA_real_
          ),
        
        rho_conn_mean =
          ifelse(
            length(rho_conn_rep) > 0,
            mean(
              rho_conn_rep,
              na.rm = TRUE
            ),
            NA_real_
          ),
        
        rho_conn_sd =
          ifelse(
            length(rho_conn_rep) > 1,
            sd(
              rho_conn_rep,
              na.rm = TRUE
            ),
            NA_real_
          ),
        
        rho_H2_mean =
          ifelse(
            length(rho_H2_rep) > 0,
            mean(
              rho_H2_rep,
              na.rm = TRUE
            ),
            NA_real_
          ),
        
        rho_H2_sd =
          ifelse(
            length(rho_H2_rep) > 1,
            sd(
              rho_H2_rep,
              na.rm = TRUE
            ),
            NA_real_
          )
      )
    )
    
    
    if (verbose) {
      
      cat(
        "Random n =", i,
        "| repetitions:",
        n_rep,
        "\n"
      )
    }
  }
  return(results)
}


# =============================================================================
# 12. Run abundance-based sampling
# =============================================================================

set.seed(2025)

r_curve <- get_structure_curve(
  max_n = 50,
  
  data_count_scaled =
    data_count_scaled,
  
  data_interact =
    data_interact,
  
  full_metrics =
    full_metrics,
  
  verbose = TRUE
)


# =============================================================================
# 13. Run random sampling
# =============================================================================

set.seed(2025)

r_curve_random <- get_structure_curve_random(
  max_n = 50,
  
  data_count_scaled =
    data_count_scaled,
  
  data_interact =
    data_interact,
  
  full_metrics =
    full_metrics,
  
  n_rep = 50,
  
  verbose = TRUE
)


# =============================================================================
# 14. Save sampling curve results
# =============================================================================

saveRDS(
  r_curve,
  "data/processed/spearman_r_curve.rds"
)

saveRDS(
  r_curve_random,
  "data/processed/spearman_r_curve_random.rds"
)


# =============================================================================
# 15. Inspect results
# =============================================================================

print(r_curve)

print(r_curve_random)



######################### read and plot. dont have to run radom process every time

r_curve <- readRDS("data/processed/spearman_r_curve.rds")
r_curve_random <- readRDS("data/processed/spearman_r_curve_random.rds")
# =============================================================================
# 16. Saturation model
# =============================================================================

fit_saturation_model <- function(
    df,
    response
) {
  
  df <- df %>%
    filter(
      !is.na(.data[[response]])
    )
  
  if (nrow(df) < 4)
    return(NULL)
  
  nlsLM(
    as.formula(
      paste0(
        response,
        " ~ a * n / (b + n)"
      )
    ),
    
    data = df,
    
    start = list(
      a = max(
        df[[response]],
        na.rm = TRUE
      ),
      
      b = median(
        df$n,
        na.rm = TRUE
      )
    )
  )
}


# =============================================================================
# 17. Breakpoint model
# =============================================================================

fit_breakpoint <- function(
    df,
    response
) {
  
  df <- df %>%
    filter(
      !is.na(.data[[response]])
    )
  
  if (nrow(df) < 5)
    return(NULL)
  
  lm_base <- lm(
    as.formula(
      paste0(
        response,
        " ~ n"
      )
    ),
    
    data = df
  )
  
  tryCatch(
    
    segmented(
      lm_base,
      seg.Z = ~ n
    ),
    
    error = function(e) NULL
  )
}


# =============================================================================
# 18. Plateau detection
#
# Minimum sampling effort is defined as the smallest n at which
# Spearman's rho reaches rho >= 0.8.
# =============================================================================

detect_plateau <- function(
    df,
    response,
    threshold = 0.8
) {
  
  df <- df %>%
    filter(
      !is.na(.data[[response]])
    )
  
  idx <- which(
    df[[response]] >= threshold
  )
  
  if (length(idx) == 0)
    return(NA_real_)
  
  df$n[
    min(idx)
  ]
}


# =============================================================================
# 19. Fit models
# =============================================================================

threshold_val <- 0.8


# -----------------------------------------------------------------------------
# NODF
# -----------------------------------------------------------------------------

sat_model_nodf <- fit_saturation_model(
  r_curve,
  "rho_nodf"
)

sat_model_nodf_rand <- fit_saturation_model(
  r_curve_random,
  "rho_nodf_mean"
)

bp_model_nodf <- fit_breakpoint(
  r_curve,
  "rho_nodf"
)

plateau_nodf <- detect_plateau(
  r_curve,
  "rho_nodf",
  threshold_val
)

plateau_nodf_rand <- detect_plateau(
  r_curve_random,
  "rho_nodf_mean",
  threshold_val
)


# -----------------------------------------------------------------------------
# Connectance
# -----------------------------------------------------------------------------

sat_model_conn <- fit_saturation_model(
  r_curve,
  "rho_conn"
)

sat_model_conn_rand <- fit_saturation_model(
  r_curve_random,
  "rho_conn_mean"
)

bp_model_conn <- fit_breakpoint(
  r_curve,
  "rho_conn"
)

plateau_conn <- detect_plateau(
  r_curve,
  "rho_conn",
  threshold_val
)

plateau_conn_rand <- detect_plateau(
  r_curve_random,
  "rho_conn_mean",
  threshold_val
)


# -----------------------------------------------------------------------------
# H2
# -----------------------------------------------------------------------------

sat_model_H2 <- fit_saturation_model(
  r_curve,
  "rho_H2"
)

sat_model_H2_rand <- fit_saturation_model(
  r_curve_random,
  "rho_H2_mean"
)

bp_model_H2 <- fit_breakpoint(
  r_curve,
  "rho_H2"
)

plateau_H2 <- detect_plateau(
  r_curve,
  "rho_H2",
  threshold_val
)

plateau_H2_rand <- detect_plateau(
  r_curve_random,
  "rho_H2_mean",
  threshold_val
)


# =============================================================================
# 21. Prediction data
# =============================================================================

pred_df <- data.frame(
  n = r_curve$n
)


# -----------------------------------------------------------------------------
# NODF predictions
# -----------------------------------------------------------------------------

if (!is.null(sat_model_nodf)) {
  
  coef_nodf <- coef(
    sat_model_nodf
  )
  
  pred_df$pred_nodf_top <-
    (
      coef_nodf["a"] *
        pred_df$n
    ) /
    (
      coef_nodf["b"] +
        pred_df$n
    )
}

if (!is.null(sat_model_nodf_rand)) {
  
  coef_nodf_rand <- coef(
    sat_model_nodf_rand
  )
  
  pred_df$pred_nodf_rand <-
    (
      coef_nodf_rand["a"] *
        pred_df$n
    ) /
    (
      coef_nodf_rand["b"] +
        pred_df$n
    )
}


# -----------------------------------------------------------------------------
# Connectance predictions
# -----------------------------------------------------------------------------

if (!is.null(sat_model_conn)) {
  
  coef_conn <- coef(
    sat_model_conn
  )
  
  pred_df$pred_conn_top <-
    (
      coef_conn["a"] *
        pred_df$n
    ) /
    (
      coef_conn["b"] +
        pred_df$n
    )
}

if (!is.null(sat_model_conn_rand)) {
  
  coef_conn_rand <- coef(
    sat_model_conn_rand
  )
  
  pred_df$pred_conn_rand <-
    (
      coef_conn_rand["a"] *
        pred_df$n
    ) /
    (
      coef_conn_rand["b"] +
        pred_df$n
    )
}


# -----------------------------------------------------------------------------
# H2 predictions
# -----------------------------------------------------------------------------

if (!is.null(sat_model_H2)) {
  
  coef_H2 <- coef(
    sat_model_H2
  )
  
  pred_df$pred_H2_top <-
    (
      coef_H2["a"] *
        pred_df$n
    ) /
    (
      coef_H2["b"] +
        pred_df$n
    )
}

if (!is.null(sat_model_H2_rand)) {
  
  coef_H2_rand <- coef(
    sat_model_H2_rand
  )
  
  pred_df$pred_H2_rand <-
    (
      coef_H2_rand["a"] *
        pred_df$n
    ) /
    (
      coef_H2_rand["b"] +
        pred_df$n
    )
}




# =============================================================================
# 22. Connectance plot
# =============================================================================

p_conn <- ggplot() +
  
  geom_ribbon(
    data = r_curve_random,
    aes(
      x = n,
      ymin = rho_conn_mean - rho_conn_sd,
      ymax = rho_conn_mean + rho_conn_sd
    ),
    fill = "grey80",
    alpha = 0.5
  ) +
  
  # geom_line(
  #   data = r_curve_random,
  #   aes(
  #     x = n,
  #     y = rho_conn_mean
  #   ),
  #   color = "grey50",
  #   linetype = "dashed",
  #   linewidth = 1,
  #   alpha = 0.7
  # ) +
  
  geom_point(
    data = r_curve,
    aes(
      x = n,
      y = rho_conn
    ),
    color = "#B84E22",
    size = 2
  ) +
  
  geom_vline(
    xintercept = plateau_conn,
    color = "#F5A88A",
    linewidth = 1.2,
    alpha = 0.8
  ) +
  
  geom_vline(
    xintercept = plateau_conn_rand,
    color = "grey50",
    linetype = "dashed",
    linewidth = 1.2,
    alpha = 0.8
  ) +
  
  annotate(
    "text",
    x = plateau_conn - 12,
    y = max(r_curve$rho_conn, na.rm = TRUE) - 0.15,
    label = paste0(
      "Abundant:\nn = ",
      plateau_conn
    ),
    color = "#B84E22",
    vjust = 1,
    hjust = -0.1,
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = plateau_conn_rand + 5,
    y = max(
      r_curve_random$rho_conn_mean,
      na.rm = TRUE
    ) - 0.15,
    label = paste0(
      "Random:\nn = ",
      plateau_conn_rand
    ),
    color = "grey50",
    vjust = 1,
    hjust = 0,
    size = 3.5
  ) +
  
  labs(
    title = "Connectance",
    x = "Sampled plant species",
    y = NULL
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 12)
  )

# =============================================================================
# 23. NODF plot
# =============================================================================

p_nodf <- ggplot() +
  
  geom_ribbon(
    data = r_curve_random,
    aes(
      x = n,
      ymin = rho_nodf_mean - rho_nodf_sd,
      ymax = rho_nodf_mean + rho_nodf_sd
    ),
    fill = "grey80",
    alpha = 0.5
  ) +
  
  # geom_line(
  #   data = r_curve_random,
  #   aes(
  #     x = n,
  #     y = rho_nodf_mean
  #   ),
  #   color = "grey50",
  #   linetype = "dashed",
  #   linewidth = 1,
  #   alpha = 0.7
# ) +

geom_point(
  data = r_curve,
  aes(
    x = n,
    y = rho_nodf
  ),
  color = "#355A9A",
  size = 2
) +
  
  geom_vline(
    xintercept = plateau_nodf,
    color = "#355A9A",
    linewidth = 1.2,
    alpha = 0.8
  ) +
  
  geom_vline(
    xintercept = plateau_nodf_rand,
    color = "grey50",
    linetype = "dashed",
    linewidth = 1.2,
    alpha = 0.8
  ) +
  
  annotate(
    "text",
    x = plateau_nodf - 12,
    y = max(r_curve$rho_nodf, na.rm = TRUE) - 0.15,
    label = paste0(
      "Abundant:\nn = ",
      plateau_nodf
    ),
    color = "#355A9A",
    vjust = 1,
    hjust = -0.1,
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = plateau_nodf_rand - 12,
    y = max(
      r_curve_random$rho_nodf_mean,
      na.rm = TRUE
    ) - 0.15,
    label = paste0(
      "Random:\nn = ",
      plateau_nodf_rand
    ),
    color = "grey50",
    vjust = 1,
    hjust = -0.1,
    size = 3.5
  ) +
  
  labs(
    title = "Nestedness",
    x = "Sampled plant species",
    y = NULL
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 12)
  )

# =============================================================================
# 24. H2 plot
# =============================================================================

p_H2 <- ggplot() +
  
  geom_ribbon(
    data = r_curve_random,
    aes(
      x = n,
      ymin = rho_H2_mean - rho_H2_sd,
      ymax = rho_H2_mean + rho_H2_sd
    ),
    fill = "grey80",
    alpha = 0.5
  ) +
  
  # geom_line(
  #   data = r_curve_random,
  #   aes(
  #     x = n,
  #     y = rho_H2_mean
  #   ),
  #   color = "grey50",
  #   linetype = "dashed",
  #   linewidth = 1,
  #   alpha = 0.7
  # ) +
  
  geom_point(
    data = r_curve,
    aes(
      x = n,
      y = rho_H2
    ),
    color = "#6A4C93",
    size = 2
  ) +
  
  geom_vline(
    xintercept = plateau_H2,
    color = "#6A4C93",
    linewidth = 1.2,
    alpha = 0.8
  ) +
  
  geom_vline(
    xintercept = plateau_H2_rand,
    color = "grey50",
    linetype = "dashed",
    linewidth = 1.2,
    alpha = 0.8
  ) +
  
  annotate(
    "text",
    x = plateau_H2 - 12,
    y = max(r_curve$rho_H2, na.rm = TRUE) - 0.15,
    label = paste0(
      "Abundant:\nn = ",
      plateau_H2
    ),
    color = "#6A4C93",
    vjust = 1,
    hjust = -0.1,
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = plateau_H2_rand - 12,
    y = max(
      r_curve_random$rho_H2_mean,
      na.rm = TRUE
    ) - 0.15,
    label = paste0(
      "Random:\nn = ",
      plateau_H2_rand
    ),
    color = "grey50",
    vjust = 1,
    hjust = -0.1,
    size = 3.5
  ) +
  
  labs(
    title = expression(bold("Selectivity (" * H[2] * "\u2032)")),
    x = "Sampled plant species",
    y = NULL
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 12)
  )


# =============================================================================
# 25. Display plots
# =============================================================================

p_nodf

p_conn

p_H2


# =============================================================================
# 26. Combine plots
# =============================================================================

combined_plot_suf <- plot_grid(
  p_conn,
  p_nodf,
  p_H2,
  ncol = 3,
  align = "hv",
  
  labels = c("(a)", "(b)", "(c)"),
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  label_size = 18,
  label_fontface = "bold"
)

y_title <- ggdraw() +
  draw_label(
    "Across-network Spearman's ρ",
    angle = 90,
    fontface = "bold",
    size = 14
  )

combined_plot_suf <- plot_grid(
  y_title,
  combined_plot_suf,
  ncol = 2,
  rel_widths = c(0.05, 1)
)

combined_plot_suf



# =============================================================================
# 27. Save final figure
# =============================================================================

ggsave(
  "D:/Chap1_TargetPlant_to_monitor/result_260723/spearman_sufficient_n_plant_network_metrics.png",
  combined_plot_suf,
  width = 11.5,
  height = 3.8,
  units = "in",
  dpi = 600
)

# =============================================================================
# 28. Save model objects
# =============================================================================

saveRDS(
  list(
    NODF = sat_model_nodf,
    NODF_random = sat_model_nodf_rand,
    
    Connectance = sat_model_conn,
    Connectance_random = sat_model_conn_rand,
    
    H2 = sat_model_H2,
    H2_random = sat_model_H2_rand
  ),
  "data/processed/sampling_curve_saturation_models.rds"
)


saveRDS(
  list(
    NODF = bp_model_nodf,
    Connectance = bp_model_conn,
    H2 = bp_model_H2
  ),
  "data/processed/sampling_curve_breakpoint_models.rds"
)


# =============================================================================
# END
# =============================================================================