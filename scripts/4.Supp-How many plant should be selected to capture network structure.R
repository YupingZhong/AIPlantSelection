
##1.Data overview
library(dplyr)
library(ggplot2) #For plotting
library(giscoR) #For plotting
library(patchwork)#For plotting (binding plots)
library(ggstar) #For plotting (cool shapes)
library(scales) #For plotting (decimals on axes)
library(tidyr)
library(viridis)
library(stringr)
library(vegan)
setwd("E:/Chap1_TargetPlant_to_monitor")
##################plant and pollinator species distribution
data_count_scaled<-readRDS("data_count_scaled_published_0526.rds")
data_interact<-readRDS("data_interact_published_0526.rds")
traits<-read.csv("test.merge.trait.csv", header = TRUE, fileEncoding = "UTF-8")

data_interact<-data_interact%>%
  mutate(
    Plant_accepted_name = str_replace_all(Plant_accepted_name, "×", "") %>%  
      str_squish()  
  )

data_merge<- merge(data_interact, data_count_scaled[,c("Flower_data_merger","Flower_count_scaled","Plant_species","Study_Network_id")], 
                   by = "Flower_data_merger",all = TRUE)%>%
  filter(!is.na(Flower_data_merger))%>%
  mutate(
    Study_Network_id = coalesce(Study_Network_id.x, Study_Network_id.y)
  ) %>%
  dplyr::select(-Study_Network_id.x, -Study_Network_id.y)%>%
  mutate(Interaction_addup = ifelse(is.na(Interaction_addup), 0, Interaction_addup))%>% # 替换 `Interaction_addup` 为 NA 的值为 0
  mutate(
    Plant_accepted_name = str_squish(str_replace_all(replace_na(Plant_accepted_name, ""), "×", ""))
  ) 


data_merge_trait <- left_join(data_merge, traits, by = "Plant_accepted_name")


########################################################
#########################
library(dplyr)
library(data.table)
library(minpack.lm)
library(segmented)
library(tidyr)
library(ggplot2)
library(bipartite)

# =========================
# 1. NODF calculator
# =========================
calc_nodf <- function(df) {
  
  if (nrow(df) == 0) return(NA_real_)
  
  dt <- as.data.table(df)
  
  mat_df <- tryCatch({
    dcast(
      dt,
      Plant_accepted_name ~ Pollinator_accepted_name,
      value.var = "Interaction_addup",
      #fun.aggregate = sum,
      fill = 0
    )
  }, error = function(e) return(NULL))
  
  if (is.null(mat_df) || ncol(mat_df) < 2) return(NA_real_)
  
  mat <- as.matrix(mat_df[, -1, with = FALSE])
  mat[is.na(mat)] <- 0
  
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  
  if (nrow(mat) < 2 || ncol(mat) < 2 || sum(mat) == 0) {
    return(NA_real_)
  }
  
  tryCatch(
    bipartite::nested(mat, method = "NODF"),
    error = function(e) NA_real_
  )
}

#====================
# Connectance
#=======================

calculate_connectance <- function(df) {
  
  if (nrow(df) == 0) return(NA_real_)
  
  dt <- as.data.table(df)
  
  mat_df <- tryCatch({
    dcast(dt,
          Plant_accepted_name ~ Pollinator_accepted_name,
          value.var = "Interaction_addup",
          fill = 0)
  }, error = function(e) return(NULL))
  
  if (is.null(mat_df)) return(NA_real_)
  
  mat <- as.matrix(mat_df[, -1, with = FALSE])
  mat[is.na(mat)] <- 0
  
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NA_real_)
  
  num_links <- sum(mat > 0)
  num_possible <- nrow(mat) * ncol(mat)
  
  if (num_possible == 0) return(NA_real_)
  
  num_links / num_possible
}


# =========================
# 2. Sampling curve 
# =========================
get_structure_curve <- function(max_n,
                                data_count_scaled,
                                data_interact,
                                min_networks = 5,
                                verbose = TRUE) {
  
  results <- data.frame(
    n = integer(),
    r_nodf = numeric(),
    r_conn = numeric()
  )
  
  data_count_scaled <- data_count_scaled %>%
    mutate(Plant_species = trimws(Plant_species))
  
  data_interact <- data_interact %>%
    mutate(
      Plant_original_name = trimws(Plant_original_name),
      Interaction_addup = replace_na(Interaction_addup, 0)
    )
  
  networks <- unique(data_interact$Study_Network_id)
  
  full_nodf <- list()
  full_conn <- list()
  full_mat_list <- list()
  
  # =========================
  # FULL MATRICES (once only)
  # =========================
  for (net in networks) {
    
    full_df_net <- data_interact %>%
      filter(Study_Network_id == net)
    
    dt <- as.data.table(full_df_net)
    
    mat_df <- dcast(
      dt,
      Plant_accepted_name ~ Pollinator_accepted_name,
      value.var = "Interaction_addup",
      fun.aggregate = sum,
      fill = 0
    )
    
    mat <- as.matrix(mat_df[, -1, with = FALSE])
    rownames(mat) <- mat_df[[1]]
    
    mat[is.na(mat)] <- 0
    
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
    
    full_mat_list[[net]] <- mat
    
    full_nodf[[net]] <- tryCatch(
      bipartite::nested(mat, method = "NODF"),
      error = function(e) NA_real_
    )
    
    full_conn[[net]] <- sum(mat > 0) / (nrow(mat) * ncol(mat))
  }
  
  # =========================
  # LOOP n
  # =========================
  for (i in 1:max_n) {
    
    nodf_sub <- c()
    nodf_full <- c()
    
    conn_sub <- c()
    conn_full <- c()
    
    for (net in networks) {
      
      full_mat <- full_mat_list[[net]]
      
      # etwork-specific full data
      full_df_net <- data_interact %>%
        filter(Study_Network_id == net)
      
      top_species <- data_count_scaled %>%
        filter(Study_Network_id == net) %>%
        arrange(desc(Flower_count_scaled)) %>%
        slice_head(n = i) %>%
        pull(Plant_species)
      

      sub_df <- full_df_net %>%
        filter(Plant_original_name %in% top_species)
      
      dt_sub <- as.data.table(sub_df)
      
      sub_mat_df <- dcast(
        dt_sub,
        Plant_accepted_name ~ Pollinator_accepted_name,
        value.var = "Interaction_addup",
        fun.aggregate = sum,
        fill = 0
      )
      
      sub_mat <- as.matrix(sub_mat_df[, -1, with = FALSE])
      rownames(sub_mat) <- sub_mat_df[[1]]
      
      sub_mat[is.na(sub_mat)] <- 0
      
      sub_mat <- sub_mat[rowSums(sub_mat) > 0, , drop = FALSE]
      sub_mat <- sub_mat[, colSums(sub_mat) > 0, drop = FALSE]
      
      if (nrow(sub_mat) < 2 || ncol(sub_mat) < 2) next
      
      sub_nodf <- tryCatch(
        bipartite::nested(sub_mat, method = "NODF"),
        error = function(e) NA_real_
      )
      
      sub_conn <- sum(sub_mat > 0) / (nrow(sub_mat) * ncol(sub_mat))
      
      if (!is.na(sub_nodf) && !is.na(full_nodf[[net]])) {
        nodf_sub <- c(nodf_sub, sub_nodf)
        nodf_full <- c(nodf_full, full_nodf[[net]])
      }
      
      if (!is.na(sub_conn) && !is.na(full_conn[[net]])) {
        conn_sub <- c(conn_sub, sub_conn)
        conn_full <- c(conn_full, full_conn[[net]])
      }
    }
    
    r_nodf <- if (length(nodf_sub) > 2)
      cor(nodf_sub, nodf_full) else NA_real_
    
    r_conn <- if (length(conn_sub) > 2)
      cor(conn_sub, conn_full) else NA_real_
    
    results <- rbind(
      results,
      data.frame(
        n = i,
        r_nodf = r_nodf,
        r_conn = r_conn
      )
    )
    
    if (verbose) {
      cat("n =", i,
          "| NODF n =", length(nodf_sub),
          "| Conn n =", length(conn_sub), "\n")
    }
  }
  
  return(results)
}
#===================
#-----random sampling-----------
#======================
get_structure_curve_random <- function(max_n,
                                       data_count_scaled,
                                       data_interact,
                                       n_rep = 100,
                                       verbose = TRUE) {
  
  library(data.table)
  library(dplyr)
  library(bipartite)
  
  results <- data.frame(
    n = integer(),
    r_nodf = numeric(),
    r_conn = numeric()
  )
  
  networks <- unique(data_interact$Study_Network_id)
  
  # ===== FULL NETWORK =====
  full_nodf <- list()
  full_conn <- list()
  
  for (net in networks) {
    
    full_df_net <- data_interact %>%
      filter(Study_Network_id == net)
    
    dt <- as.data.table(full_df_net)
    
    mat_df <- dcast(
      dt,
      Plant_accepted_name ~ Pollinator_accepted_name,
      value.var = "Interaction_addup",
      fun.aggregate = sum,
      fill = 0
    )
    
    mat <- as.matrix(mat_df[, -1, with = FALSE])
    mat[is.na(mat)] <- 0
    
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
    
    if (nrow(mat) < 2 || ncol(mat) < 2) next
    
    full_nodf[[net]] <- bipartite::nested(mat, method = "NODF")
    full_conn[[net]] <- sum(mat > 0) / (nrow(mat) * ncol(mat))
  }
  
  # ===== LOOP n =====
  for (i in 1:max_n) {
    
    r_nodf_rep <- c()
    r_conn_rep <- c()
    
    for (rep in 1:n_rep) {
      
      nodf_sub <- c()
      nodf_full <- c()
      
      conn_sub <- c()
      conn_full <- c()
      
      for (net in networks) {
        
        full_df_net <- data_interact %>%
          filter(Study_Network_id == net)
        
        # RANDOM sampling
        all_species <- unique(
          data_count_scaled$Plant_species[
            data_count_scaled$Study_Network_id == net
          ]
        )
        
        if (length(all_species) < i) next
        
        random_species <- sample(all_species, i)
        
        sub_df <- full_df_net %>%
          filter(Plant_original_name %in% random_species)
        
        dt_sub <- as.data.table(sub_df)
        
        sub_mat_df <- dcast(
          dt_sub,
          Plant_accepted_name ~ Pollinator_accepted_name,
          value.var = "Interaction_addup",
          fun.aggregate = sum,
          fill = 0
        )
        
        sub_mat <- as.matrix(sub_mat_df[, -1, with = FALSE])
        sub_mat[is.na(sub_mat)] <- 0
        
        sub_mat <- sub_mat[rowSums(sub_mat) > 0, , drop = FALSE]
        sub_mat <- sub_mat[, colSums(sub_mat) > 0, drop = FALSE]
        
        if (nrow(sub_mat) < 2 || ncol(sub_mat) < 2) next
        
        sub_nodf <- bipartite::nested(sub_mat, method = "NODF")
        sub_conn <- sum(sub_mat > 0) / (nrow(sub_mat) * ncol(sub_mat))
        
        if (!is.na(sub_nodf) && !is.na(full_nodf[[net]])) {
          nodf_sub <- c(nodf_sub, sub_nodf)
          nodf_full <- c(nodf_full, full_nodf[[net]])
        }
        
        if (!is.na(sub_conn) && !is.na(full_conn[[net]])) {
          conn_sub <- c(conn_sub, sub_conn)
          conn_full <- c(conn_full, full_conn[[net]])
        }
      }
      
      if (length(nodf_sub) > 2) {
        r_nodf_rep <- c(r_nodf_rep, cor(nodf_sub, nodf_full))
      }
      
      if (length(conn_sub) > 2) {
        r_conn_rep <- c(r_conn_rep, cor(conn_sub, conn_full))
      }
    }
    
    results <- rbind(
      results,
      data.frame(
        n = i,
        
        r_nodf_mean = mean(r_nodf_rep, na.rm = TRUE),
        r_nodf_sd   = sd(r_nodf_rep, na.rm = TRUE),
        
        r_conn_mean = mean(r_conn_rep, na.rm = TRUE),
        r_conn_sd   = sd(r_conn_rep, na.rm = TRUE)
      )
    )
    
    
    if (verbose) {
      cat("Random n =", i, "\n")
    }
  }
  
  return(results)
}

# =========================
# 3. Saturation model
# =========================
fit_saturation_model <- function(df, response) {
  
  df <- df %>% filter(!is.na(.data[[response]]))
  
  nlsLM(
    as.formula(paste0(response, " ~ a * n / (b + n)")),
    data = df,
    start = list(
      a = max(df[[response]], na.rm = TRUE),
      b = median(df$n)
    )
  )
}



# =========================
# 4. Breakpoint model
fit_breakpoint <- function(df, response) {
  
  df <- df %>% filter(!is.na(.data[[response]]))
  
  lm_base <- lm(as.formula(paste0(response, " ~ n")), data = df)
  
  segmented(lm_base, seg.Z = ~ n)
}


# =========================
# 5. Plateau detection 
#We defined the minimum sampling effort as the smallest number of top-abundant plant species at 
#which the correlation between subsampled and full-network nestedness exceeded 0.6.
# =========================
detect_plateau <- function(df, response, threshold = 0.6) {
  
  df <- df %>% filter(!is.na(.data[[response]]))
  
  idx <- which(df[[response]] >= threshold)
  
  if (length(idx) == 0) return(NA)
  
  return(df$n[min(idx)])
}

# =========================
# 6. RUN
# =========================
r_curve <- get_structure_curve(
  max_n = 50,
  data_count_scaled = data_count_scaled,
  data_interact = data_interact
)

r_curve_random <- get_structure_curve_random(
  max_n = 50,
  data_count_scaled = data_count_scaled,
  data_interact = data_interact,
  n_rep = 50
)

# save
# saveRDS(r_curve_random,"r_curve_random.rds")
# saveRDS(r_curve,"r_curve.rds")
r_curve_random<-readRDS("r_curve_random.rds")
r_curve<-readRDS("r_curve.rds")

`# =========================
# 7. MODELS (NODF + CONNECTANCE)
# =========================
threshold_val <- 0.8
# ---- NODF ----
sat_model_nodf <- fit_saturation_model(r_curve, "r_nodf")
sat_model_nodf_rand <- fit_saturation_model(r_curve_random, "r_nodf_mean")

bp_model_nodf  <- fit_breakpoint(r_curve, "r_nodf")

plateau_nodf <- detect_plateau(r_curve, "r_nodf", threshold_val)
plateau_nodf_rand <- detect_plateau(r_curve_random, "r_nodf_mean", threshold_val)


cat("NODF Plateau (r ≥ 0.8) at n =", plateau_nodf, "\n")


# ---- CONNECTANCE ----
sat_model_conn <- fit_saturation_model(r_curve, "r_conn")
sat_model_conn_rand <- fit_saturation_model(r_curve_random, "r_conn_mean")
bp_model_conn  <- fit_breakpoint(r_curve, "r_conn")

plateau_conn <- detect_plateau(r_curve, "r_conn", threshold_val)
plateau_conn_rand <- detect_plateau(r_curve_random, "r_conn_mean", threshold_val)

cat("Connectance Plateau (r ≥ 0.8) at n =", plateau_conn, "\n")


# =========================
# 8. PREDICTION
# =========================

# ---------- NODF ----------
pred_df_nodf <- data.frame(n = r_curve$n)

# top
coef_nodf <- coef(sat_model_nodf)
pred_df_nodf$pred_top <- (coef_nodf["a"] * pred_df_nodf$n) /
  (coef_nodf["b"] + pred_df_nodf$n)

# random
coef_nodf_rand <- coef(sat_model_nodf_rand)
pred_df_nodf$pred_rand <- (coef_nodf_rand["a"] * pred_df_nodf$n) /
  (coef_nodf_rand["b"] + pred_df_nodf$n)


# ---------- CONNECTANCE ----------
pred_df_conn <- data.frame(n = r_curve$n)

# top
coef_conn <- coef(sat_model_conn)
pred_df_conn$pred_top <- (coef_conn["a"] * pred_df_conn$n) /
  (coef_conn["b"] + pred_df_conn$n)

# random
coef_conn_rand <- coef(sat_model_conn_rand)
pred_df_conn$pred_rand <- (coef_conn_rand["a"] * pred_df_conn$n) /
  (coef_conn_rand["b"] + pred_df_conn$n)
# =========================
# 9. FINAL FIGURE 
# =========================

library(ggplot2)
#======NODF
# ===== NODF PLOT =====
p_nodf <- ggplot() +
  
  # ===== RANDOM CI =====
geom_ribbon(
  data = r_curve_random,
  aes(x = n, ymin = r_nodf_mean - r_nodf_sd, ymax = r_nodf_mean + r_nodf_sd),
  fill = "grey80",
  alpha = 0.5
) +
  
  # ===== RANDOM OBSERVED LINE =====
geom_line(
  data = r_curve_random,
  aes(x = n, y = r_nodf_mean),
  color = "grey50",           # 改为灰色
  linetype = "dashed",
  linewidth = 1,
  alpha = 0.7
) +
  
  # ===== TOP POINTS =====
geom_point(
  data = r_curve,
  aes(x = n, y = r_nodf),
  color = "#D55E00",
  size = 2
) +
  
  # ===== TOP FIT =====
geom_line(
  data = pred_df_nodf,
  aes(x = n, y = pred_top),
  color = "#D55E00",
  linewidth = 1.3
) +
  
  # ===== RANDOM PREDICTION LINE =====
geom_line(
  data = pred_df_nodf,
  aes(x = n, y = pred_rand),
  color = "grey50",           # 改为灰色
  linetype = "dotted",        # 改为点线区分
  linewidth = 0.9
) +
  
  # ===== PLATEAU LINES =====
geom_vline(
  xintercept = plateau_nodf,
  color = "#D55E00",
  linewidth = 1.2,
  alpha = 0.8
) +
  
  geom_vline(
    xintercept = plateau_nodf_rand,
    linetype = "dashed",
    color = "grey50",           # 改为灰色
    linewidth = 1.2,
    alpha = 0.8
  ) +
  
  # ===== LABELS（统一高度设置）=====
annotate(
  "text",
  x = plateau_nodf-12,
  y = max(r_curve$r_nodf, na.rm = TRUE) - 0.15,  # 统一偏移
  label = paste0("Abundant:\nn = ", plateau_nodf),
  color = "#D55E00",
  vjust = 1,
  hjust = -0.1,
  size = 3.5
) +
  
  annotate(
    "text",
    x = plateau_nodf_rand,
    y = max(r_curve$r_nodf, na.rm = TRUE) - 0.15,  # 统一偏移
    label = paste0("Random:\nn = ", plateau_nodf_rand),
    color = "grey50",           # 改为灰色
    vjust = 1,
    hjust = 1.1,
    size = 3.5
  ) +
  
  labs(
    x = "Number of plant species (n)",
    y = "Correlation (R) — NODF"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )


# ===== CONNECTANCE PLOT =====
p_conn <- ggplot() +
  
  geom_ribbon(
    data = r_curve_random,
    aes(x = n, ymin = r_conn_mean - r_conn_sd, ymax = r_conn_mean + r_conn_sd),
    fill = "grey80",
    alpha = 0.5
  ) +
  
  # ===== RANDOM OBSERVED LINE =====
geom_line(
  data = r_curve_random,
  aes(x = n, y = r_conn_mean),
  color = "grey50",           # 改为灰色
  linetype = "dashed",
  linewidth = 1,              # 加上linewidth
  alpha = 0.7
) +
  
  geom_point(
    data = r_curve,
    aes(x = n, y = r_conn),
    color = "#0072B2",
    size = 2
  ) +
  
  geom_line(
    data = pred_df_conn,
    aes(x = n, y = pred_top),
    color = "#0072B2",
    linewidth = 1.3
  ) +
  
  # ===== RANDOM PREDICTION LINE =====
geom_line(
  data = pred_df_conn,
  aes(x = n, y = pred_rand),
  color = "grey50",           # 改为灰色
  linetype = "dotted",        # 改为点线区分
  linewidth = 0.9
) +
  
  geom_vline(
    xintercept = plateau_conn,
    color = "#0072B2",
    linewidth = 1.2,            # 加上linewidth
    alpha = 0.8
  ) +
  
  geom_vline(
    xintercept = plateau_conn_rand,
    linetype = "dashed",
    color = "grey50",           # 改为灰色
    linewidth = 1.2,            # 加上linewidth
    alpha = 0.8
  ) +
  
  # ===== LABELS（统一高度设置）=====
annotate(
  "text",
  x = plateau_conn-12,
  y = max(r_curve$r_conn, na.rm = TRUE) - 0.15,  # 统一偏移
  label = paste0("Abundant:\nn = ", plateau_conn),
  color = "#0072B2",
  vjust = 1,
  hjust = -0.1,
  size = 3.5
) +
  
  annotate(
    "text",
    x = plateau_conn_rand+5,
    y = max(r_curve$r_conn, na.rm = TRUE) - 0.15,  # 统一偏移
    label = paste0("Random:\nn = ", plateau_conn_rand),
    color = "grey50",           # 改为灰色
    vjust = 1,
    hjust = 1.1,
    size = 3.5
  ) +
  
  labs(
    x = "Number of plant species (n)",
    y = "Correlation (R) — Connectance"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )
# =========================
# SHOW BOTH PLOTS
# =========================
library(cowplot)
p_nodf
p_conn

combined_plot_suf <- plot_grid(
  p_nodf, p_conn, # 右边
  ncol = 2                       
)

combined_plot_suf
ggsave("/Chap1_TargetPlant_to_monitor/result_260723/sufficient_n_plant.png", combined_plot_suf, width = 8, 
       height = 3.5, units = "in", dpi = 300)
