
############################################################
## Visualization of network metrics
##
## Metrics:
## 1. Connectance
## 2. Nestedness (NODF)
## 3. H2 specialization
##
## Four subsampling strategies:
## a. Abundance-based
## b. Abundance + flower shape
## c. Phylogenetic
## d. Random
##
## Lines represent linear regression fits for visualization.
##
## rho (ρ) and p indicate Spearman rank correlations.
##
## Solid lines indicate significant Spearman correlations
## (p < 0.05).
##
## Dashed lines indicate non-significant correlations.
############################################################


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)


############################################################
# 1. Read data
############################################################

all_nestedness_df <- readRDS(
  "data/processed/all_nestedness_df.rds"
)

all_connectance_df <- readRDS(
  "data/processed/all_connectance_df.rds"
)

all_H2_df <- readRDS(
  "data/processed/all_H2_df.rds"
)

random_all <- readRDS(
  "data/processed/random_all_metrics.rds"
)


############################################################
# 2. Colors
############################################################

my_colors <- c(
  "10" = "#66C2A5",
  "5"  = "#B07AA1",
  "3"  = "#E9C46A"
)


############################################################
# 3. General theme
############################################################

base_theme <- theme_classic(
  base_size = 12
) +
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 11
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text = element_text(
      color = "black",
      size = 11
    ),
    
    legend.position = "none",
    
    plot.margin = margin(
      t = 2,
      r = 8,
      b = 2,
      l = 2
    )
  )


############################################################
# 4. Function:
#    calculate Spearman rho and p
############################################################

get_spearman_stats <- function(
    df,
    x_col,
    y_col,
    group_col
) {
  
  df %>%
    group_by(
      across(all_of(group_col))
    ) %>%
    summarise(
      
      rho = {
        
        tmp <- pick(
          all_of(c(x_col, y_col))
        )
        
        tmp <- tmp[
          complete.cases(tmp),
          ,
          drop = FALSE
        ]
        
        if (nrow(tmp) < 3) {
          
          NA_real_
          
        } else {
          
          cor.test(
            tmp[[x_col]],
            tmp[[y_col]],
            method = "spearman",
            exact = FALSE
          )$estimate
        }
      },
      
      p = {
        
        tmp <- pick(
          all_of(c(x_col, y_col))
        )
        
        tmp <- tmp[
          complete.cases(tmp),
          ,
          drop = FALSE
        ]
        
        if (nrow(tmp) < 3) {
          
          NA_real_
          
        } else {
          
          cor.test(
            tmp[[x_col]],
            tmp[[y_col]],
            method = "spearman",
            exact = FALSE
          )$p.value
        }
      },
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      linetype = case_when(
        
        is.na(p) ~ "dashed",
        
        p < 0.05 ~ "solid",
        
        TRUE ~ "dashed"
      )
    )
}


############################################################
# 5. General plotting function
############################################################

plot_metric <- function(
    df,
    x_col,
    y_col,
    group_col,
    linetype_df,
    title,
    x_label,
    y_label,
    y_limits,
    jitter_width = 0.02,
    jitter_height = 0.2,
    cor_label_x = "left",
    cor_label_y = NULL
) {
  
  # --------------------------------------------------------
  # Number of groups
  # --------------------------------------------------------
  
  n_group <- length(
    unique(
      df[[group_col]]
    )
  )
  
  
  # --------------------------------------------------------
  # Correlation label positions
  # --------------------------------------------------------
  
  if (is.null(cor_label_y)) {
    
    cor_label_y <- seq(
      0.98,
      0.98 - 0.04 * (n_group - 1),
      length.out = n_group
    )
    
  }
  
  
  # --------------------------------------------------------
  # Plot
  # --------------------------------------------------------
  
  ggplot(
    
    df,
    
    aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = .data[[group_col]]
    )
    
  ) +
    
    
    ########################################################
  # Points
  ########################################################
  
  geom_point(
    
    shape = 1,
    
    size = 2,
    
    position = position_jitter(
      width = jitter_width,
      height = jitter_height
    )
  ) +
    
    
    ########################################################
  # Linear regression line
  #
  # IMPORTANT:
  # This line is only for visualization.
  #
  # Statistical test = Spearman correlation.
  ########################################################
  
  geom_smooth(
    
    method = "lm",
    
    se = FALSE,
    
    linewidth = 0.8,
    
    aes(
      linetype = .data[[group_col]]
    ),
    
    show.legend = FALSE
  ) +
    
    
    ########################################################
  # Significant = solid
  # Non-significant = dashed
  #
  # Based on Spearman p-value
  ########################################################
  
  scale_linetype_manual(
    
    values = setNames(
      linetype_df$linetype,
      linetype_df[[group_col]]
    )
  ) +
    
    
    ########################################################
  # Spearman rho and p
  #
  # cor.coef.name = "rho"
  #
  # This makes ggpubr display:
  #
  # ρ = 0.XX, p = 0.XXX
  ########################################################
  
  stat_cor(
    
    aes(
      group = .data[[group_col]],
      color = .data[[group_col]]
    ),
    
    method = "spearman",
    
    cor.coef.name = "rho",
    
    size = 3.5,
    
    label.x.npc = cor_label_x,
    
    label.y.npc = cor_label_y,
    
    output.type = "expression"
  ) +
    
    
    ########################################################
  # Colors
  ########################################################
  
  scale_color_manual(
    
    values = my_colors
  ) +
    
    
    ########################################################
  # Y axis
  ########################################################
  
  scale_y_continuous(
    
    limits = y_limits
  ) +
    
    
    ########################################################
  # Labels
  ########################################################
  
  labs(
    
    title = title,
    
    x = x_label,
    
    y = y_label
  ) +
    
    
    ########################################################
  # Theme
  ########################################################
  
  base_theme
}


############################################################
# 6. Prepare RANDOM data
############################################################

random_all <- random_all %>%
  
  mutate(
    
    Method = trimws(Method),
    
    Method = case_when(
      
      Method %in% c(
        "random10",
        "Random10",
        "random 10"
      ) ~ "10",
      
      Method %in% c(
        "random5",
        "Random5",
        "random 5"
      ) ~ "5",
      
      Method %in% c(
        "random3",
        "Random3",
        "random 3"
      ) ~ "3",
      
      TRUE ~ Method
    )
  )


random_all$Method <- factor(
  
  random_all$Method,
  
  levels = c(
    "10",
    "5",
    "3"
  )
)


############################################################
# ==========================================================
# NESTEDNESS
# ==========================================================
############################################################


############################################################
# 7. Nestedness - Abundance
############################################################

nested_abun <- all_nestedness_df %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    Nestedness_orig,
    
    Nestedness_Abun10,
    Nestedness_Abun5,
    Nestedness_Abun3
  ) %>%
  
  pivot_longer(
    
    cols = starts_with(
      "Nestedness_Abun"
    ),
    
    names_to = "Level",
    
    values_to = "Value"
  ) %>%
  
  mutate(
    
    Level = factor(
      
      sub(
        "Nestedness_Abun",
        "",
        Level
      ),
      
      levels = c(
        "10",
        "5",
        "3"
      )
    )
  ) %>%
  
  filter(
    !is.na(Value)
  )


nested_abun_sig <- get_spearman_stats(
  
  nested_abun,
  
  x_col = "Nestedness_orig",
  
  y_col = "Value",
  
  group_col = "Level"
)


p_nested_abun <- plot_metric(
  
  nested_abun,
  
  x_col = "Nestedness_orig",
  
  y_col = "Value",
  
  group_col = "Level",
  
  linetype_df = nested_abun_sig,
  
  title = "a. Abundance-based",
  
  x_label = "NODF of Full Network",
  
  y_label = "NODF in Subnetwork",
  
  y_limits = c(
    0,
    60
  ),
  
  jitter_width = 0.03,
  
  jitter_height = 0.5
)


############################################################
# 8. Nestedness - Flower shape
############################################################

nested_flw <- all_nestedness_df %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    Nestedness_orig,
    
    Nestedness_FlwShape5,
    Nestedness_FlwShape3
  ) %>%
  
  pivot_longer(
    
    cols = starts_with(
      "Nestedness_FlwShape"
    ),
    
    names_to = "Level",
    
    values_to = "Value"
  ) %>%
  
  mutate(
    
    Level = case_when(
      
      grepl(
        "5$",
        Level
      ) ~ "5",
      
      grepl(
        "3$",
        Level
      ) ~ "3"
    ),
    
    Level = factor(
      
      Level,
      
      levels = c(
        "5",
        "3"
      )
    )
  ) %>%
  
  filter(
    !is.na(Value)
  )


nested_flw_sig <- get_spearman_stats(
  
  nested_flw,
  
  x_col = "Nestedness_orig",
  
  y_col = "Value",
  
  group_col = "Level"
)


p_nested_flw <- plot_metric(
  
  nested_flw,
  
  x_col = "Nestedness_orig",
  
  y_col = "Value",
  
  group_col = "Level",
  
  linetype_df = nested_flw_sig,
  
  title = "b. Abundance + flower shape",
  
  x_label = "NODF of Full Network",
  
  y_label = "NODF in Subnetwork",
  
  y_limits = c(
    0,
    60
  ),
  
  jitter_width = 0.03,
  
  jitter_height = 0.5
)


############################################################
# 9. Nestedness - Phylogenetic
############################################################

nested_phylo <- all_nestedness_df %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    Nestedness_orig,
    
    Nestedness_Pylo10,
    Nestedness_Pylo5,
    Nestedness_Pylo3
  ) %>%
  
  pivot_longer(
    
    cols = starts_with(
      "Nestedness_Pylo"
    ),
    
    names_to = "Level",
    
    values_to = "Value"
  ) %>%
  
  mutate(
    
    Level = factor(
      
      sub(
        "Nestedness_Pylo",
        "",
        Level
      ),
      
      levels = c(
        "10",
        "5",
        "3"
      )
    )
  ) %>%
  
  filter(
    !is.na(Value)
  )


nested_phylo_sig <- get_spearman_stats(
  
  nested_phylo,
  
  x_col = "Nestedness_orig",
  
  y_col = "Value",
  
  group_col = "Level"
)


p_nested_phylo <- plot_metric(
  
  nested_phylo,
  
  x_col = "Nestedness_orig",
  
  y_col = "Value",
  
  group_col = "Level",
  
  linetype_df = nested_phylo_sig,
  
  title = "c. Phylogenetic",
  
  x_label = "NODF of Full Network",
  
  y_label = "NODF in Subnetwork",
  
  y_limits = c(
    0,
    60
  ),
  
  jitter_width = 0.03,
  
  jitter_height = 0.5
)


############################################################
# 10. Nestedness - Random
############################################################

random_nodf <- random_all %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    Method,
    
    Random_NODF
  ) %>%
  
  group_by(
    
    Study_Network_id,
    
    Method
  ) %>%
  
  summarise(
    
    Value = mean(
      Random_NODF,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  left_join(
    
    all_nestedness_df %>%
      
      dplyr::select(
        
        Study_Network_id,
        
        Nestedness_orig
      ),
    
    by = "Study_Network_id"
  ) %>%
  
  filter(
    !is.na(Value)
  )


random_nodf_sig <- get_spearman_stats(
  
  random_nodf,
  
  x_col = "Nestedness_orig",
  
  y_col = "Value",
  
  group_col = "Method"
)


p_nested_random <- plot_metric(
  
  random_nodf,
  
  x_col = "Nestedness_orig",
  
  y_col = "Value",
  
  group_col = "Method",
  
  linetype_df = random_nodf_sig,
  
  title = "d. Random",
  
  x_label = "NODF of Full Network",
  
  y_label = "NODF in Subnetwork",
  
  y_limits = c(
    0,
    60
  ),
  
  jitter_width = 0.02,
  
  jitter_height = 0.2,
  
  cor_label_x = 0.33
)


############################################################
# ==========================================================
# CONNECTANCE
# ==========================================================
############################################################


############################################################
# 11. Connectance - Abundance
############################################################

conn_abun <- all_connectance_df %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    Connectance_orig,
    
    Connectance_Abun10,
    Connectance_Abun5,
    Connectance_Abun3
  ) %>%
  
  pivot_longer(
    
    cols = starts_with(
      "Connectance_Abun"
    ),
    
    names_to = "Level",
    
    values_to = "Value"
  ) %>%
  
  mutate(
    
    Level = factor(
      
      sub(
        "Connectance_Abun",
        "",
        Level
      ),
      
      levels = c(
        "10",
        "5",
        "3"
      )
    )
  ) %>%
  
  filter(
    !is.na(Value)
  )


conn_abun_sig <- get_spearman_stats(
  
  conn_abun,
  
  x_col = "Connectance_orig",
  
  y_col = "Value",
  
  group_col = "Level"
)


p_conn_abun <- plot_metric(
  
  conn_abun,
  
  x_col = "Connectance_orig",
  
  y_col = "Value",
  
  group_col = "Level",
  
  linetype_df = conn_abun_sig,
  
  title = "a. Abundance-based",
  
  x_label = "Connectance of Full Network",
  
  y_label = "Connectance in Subnetwork",
  
  y_limits = c(
    0,
    1
  ),
  
  jitter_width = 0.005,
  
  jitter_height = 0.005
)


############################################################
# 12. Connectance - Flower shape
############################################################

conn_flw <- all_connectance_df %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    Connectance_orig,
    
    Connectance_FlwShape5,
    Connectance_FlwShape3
  ) %>%
  
  pivot_longer(
    
    cols = starts_with(
      "Connectance_FlwShape"
    ),
    
    names_to = "Level",
    
    values_to = "Value"
  ) %>%
  
  mutate(
    
    Level = case_when(
      
      grepl(
        "5$",
        Level
      ) ~ "5",
      
      grepl(
        "3$",
        Level
      ) ~ "3"
    ),
    
    Level = factor(
      
      Level,
      
      levels = c(
        "5",
        "3"
      )
    )
  ) %>%
  
  filter(
    !is.na(Value)
  )


conn_flw_sig <- get_spearman_stats(
  
  conn_flw,
  
  x_col = "Connectance_orig",
  
  y_col = "Value",
  
  group_col = "Level"
)


p_conn_flw <- plot_metric(
  
  conn_flw,
  
  x_col = "Connectance_orig",
  
  y_col = "Value",
  
  group_col = "Level",
  
  linetype_df = conn_flw_sig,
  
  title = "b. Abundance + flower shape",
  
  x_label = "Connectance of Full Network",
  
  y_label = "Connectance in Subnetwork",
  
  y_limits = c(
    0,
    1
  ),
  
  jitter_width = 0.005,
  
  jitter_height = 0.005
)


############################################################
# 13. Connectance - Phylogenetic
############################################################

conn_phylo <- all_connectance_df %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    Connectance_orig,
    
    Connectance_Pylo10,
    Connectance_Pylo5,
    Connectance_Pylo3
  ) %>%
  
  pivot_longer(
    
    cols = starts_with(
      "Connectance_Pylo"
    ),
    
    names_to = "Level",
    
    values_to = "Value"
  ) %>%
  
  mutate(
    
    Level = factor(
      
      sub(
        "Connectance_Pylo",
        "",
        Level
      ),
      
      levels = c(
        "10",
        "5",
        "3"
      )
    )
  ) %>%
  
  filter(
    !is.na(Value)
  )


conn_phylo_sig <- get_spearman_stats(
  
  conn_phylo,
  
  x_col = "Connectance_orig",
  
  y_col = "Value",
  
  group_col = "Level"
)


p_conn_phylo <- plot_metric(
  
  conn_phylo,
  
  x_col = "Connectance_orig",
  
  y_col = "Value",
  
  group_col = "Level",
  
  linetype_df = conn_phylo_sig,
  
  title = "c. Phylogenetic",
  
  x_label = "Connectance of Full Network",
  
  y_label = "Connectance in Subnetwork",
  
  y_limits = c(
    0,
    1
  ),
  
  jitter_width = 0.005,
  
  jitter_height = 0.005
)


############################################################
# 14. Connectance - Random
############################################################

random_conn <- random_all %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    Method,
    
    Random_Connectance
  ) %>%
  
  group_by(
    
    Study_Network_id,
    
    Method
  ) %>%
  
  summarise(
    
    Value = mean(
      Random_Connectance,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  left_join(
    
    all_connectance_df %>%
      
      dplyr::select(
        
        Study_Network_id,
        
        Connectance_orig
      ),
    
    by = "Study_Network_id"
  ) %>%
  
  filter(
    !is.na(Value)
  )


random_conn_sig <- get_spearman_stats(
  
  random_conn,
  
  x_col = "Connectance_orig",
  
  y_col = "Value",
  
  group_col = "Method"
)


p_conn_random <- plot_metric(
  
  random_conn,
  
  x_col = "Connectance_orig",
  
  y_col = "Value",
  
  group_col = "Method",
  
  linetype_df = random_conn_sig,
  
  title = "d. Random",
  
  x_label = "Connectance of Full Network",
  
  y_label = "Connectance in Subnetwork",
  
  y_limits = c(
    0,
    1
  ),
  
  jitter_width = 0.005,
  
  jitter_height = 0.005,
  
  cor_label_x = 0.33
)


############################################################
# ==========================================================
# H2
# ==========================================================
############################################################


############################################################
# 15. H2 - Abundance
############################################################

H2_abun <- all_H2_df %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    H2_orig,
    
    H2_Abun10,
    H2_Abun5,
    H2_Abun3
  ) %>%
  
  pivot_longer(
    
    cols = starts_with(
      "H2_Abun"
    ),
    
    names_to = "Level",
    
    values_to = "Value"
  ) %>%
  
  mutate(
    
    Level = factor(
      
      sub(
        "H2_Abun",
        "",
        Level
      ),
      
      levels = c(
        "10",
        "5",
        "3"
      )
    )
  ) %>%
  
  filter(
    !is.na(Value)
  )


H2_abun_sig <- get_spearman_stats(
  
  H2_abun,
  
  x_col = "H2_orig",
  
  y_col = "Value",
  
  group_col = "Level"
)


p_H2_abun <- plot_metric(
  
  H2_abun,
  
  x_col = "H2_orig",
  
  y_col = "Value",
  
  group_col = "Level",
  
  linetype_df = H2_abun_sig,
  
  title = "a. Abundance-based",
  
  x_label = "H2 of Full Network",
  
  y_label = "H2 in Subnetwork",
  
  y_limits = c(
    0,
    1
  ),
  
  jitter_width = 0.005,
  
  jitter_height = 0.005
)


############################################################
# 16. H2 - Flower shape
############################################################

H2_flw <- all_H2_df %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    H2_orig,
    
    H2_FlwShape5,
    H2_FlwShape3
  ) %>%
  
  pivot_longer(
    
    cols = starts_with(
      "H2_FlwShape"
    ),
    
    names_to = "Level",
    
    values_to = "Value"
  ) %>%
  
  mutate(
    
    Level = case_when(
      
      grepl(
        "5$",
        Level
      ) ~ "5",
      
      grepl(
        "3$",
        Level
      ) ~ "3"
    ),
    
    Level = factor(
      
      Level,
      
      levels = c(
        "5",
        "3"
      )
    )
  ) %>%
  
  filter(
    !is.na(Value)
  )


H2_flw_sig <- get_spearman_stats(
  
  H2_flw,
  
  x_col = "H2_orig",
  
  y_col = "Value",
  
  group_col = "Level"
)


p_H2_flw <- plot_metric(
  
  H2_flw,
  
  x_col = "H2_orig",
  
  y_col = "Value",
  
  group_col = "Level",
  
  linetype_df = H2_flw_sig,
  
  title = "b. Abundance + flower shape",
  
  x_label = "H2 of Full Network",
  
  y_label = "H2 in Subnetwork",
  
  y_limits = c(
    0,
    1
  ),
  
  jitter_width = 0.005,
  
  jitter_height = 0.005
)


############################################################
# 17. H2 - Phylogenetic
############################################################

H2_phylo <- all_H2_df %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    H2_orig,
    
    H2_Pylo10,
    H2_Pylo5,
    H2_Pylo3
  ) %>%
  
  pivot_longer(
    
    cols = starts_with(
      "H2_Pylo"
    ),
    
    names_to = "Level",
    
    values_to = "Value"
  ) %>%
  
  mutate(
    
    Level = factor(
      
      sub(
        "H2_Pylo",
        "",
        Level
      ),
      
      levels = c(
        "10",
        "5",
        "3"
      )
    )
  ) %>%
  
  filter(
    !is.na(Value)
  )


H2_phylo_sig <- get_spearman_stats(
  
  H2_phylo,
  
  x_col = "H2_orig",
  
  y_col = "Value",
  
  group_col = "Level"
)


p_H2_phylo <- plot_metric(
  
  H2_phylo,
  
  x_col = "H2_orig",
  
  y_col = "Value",
  
  group_col = "Level",
  
  linetype_df = H2_phylo_sig,
  
  title = "c. Phylogenetic",
  
  x_label = "H2 of Full Network",
  
  y_label = "H2 in Subnetwork",
  
  y_limits = c(
    0,
    1
  ),
  
  jitter_width = 0.005,
  
  jitter_height = 0.005
)


############################################################
# 18. H2 - Random
############################################################

random_H2 <- random_all %>%
  
  dplyr::select(
    
    Study_Network_id,
    
    Method,
    
    Random_H2
  ) %>%
  
  group_by(
    
    Study_Network_id,
    
    Method
  ) %>%
  
  summarise(
    
    Value = mean(
      Random_H2,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  left_join(
    
    all_H2_df %>%
      
      dplyr::select(
        
        Study_Network_id,
        
        H2_orig
      ),
    
    by = "Study_Network_id"
  ) %>%
  
  filter(
    !is.na(Value)
  )


random_H2_sig <- get_spearman_stats(
  
  random_H2,
  
  x_col = "H2_orig",
  
  y_col = "Value",
  
  group_col = "Method"
)


p_H2_random <- plot_metric(
  
  random_H2,
  
  x_col = "H2_orig",
  
  y_col = "Value",
  
  group_col = "Method",
  
  linetype_df = random_H2_sig,
  
  title = "d. Random",
  
  x_label = "H2 of Full Network",
  
  y_label = "H2 in Subnetwork",
  
  y_limits = c(
    0,
    1
  ),
  
  jitter_width = 0.005,
  
  jitter_height = 0.005,
  
  cor_label_x = 0.33
)


############################################################
# 19. Arrange NODF row
############################################################

nestedness_row <- plot_grid(
  
  p_nested_abun,
  
  p_nested_flw,
  
  p_nested_phylo,
  
  p_nested_random,
  
  ncol = 4,
  
  align = "hv",
  
  axis = "tblr",
  
  rel_widths = c(
    1,
    1,
    1,
    1
  )
)


############################################################
# 20. Arrange Connectance row
############################################################

connectance_row <- plot_grid(
  
  p_conn_abun,
  
  p_conn_flw,
  
  p_conn_phylo,
  
  p_conn_random,
  
  ncol = 4,
  
  align = "hv",
  
  axis = "tblr",
  
  rel_widths = c(
    1,
    1,
    1,
    1
  )
)


############################################################
# 21. Arrange H2 row
############################################################

H2_row <- plot_grid(
  
  p_H2_abun,
  
  p_H2_flw,
  
  p_H2_phylo,
  
  p_H2_random,
  
  ncol = 4,
  
  align = "hv",
  
  axis = "tblr",
  
  rel_widths = c(
    1,
    1,
    1,
    1
  )
)


############################################################
# 22. Create unified legend
############################################################

legend_df <- data.frame(
  
  Level = factor(
    
    c(
      "10",
      "5",
      "3"
    ),
    
    levels = c(
      "10",
      "5",
      "3"
    )
  ),
  
  x = 1:3,
  
  y = 1:3
)


legend_plot <- ggplot(
  
  legend_df,
  
  aes(
    x = x,
    y = y,
    color = Level,
    linetype = Level
  )
  
) +
  
  geom_point(
    size = 3
  ) +
  
  geom_line(
    linewidth = 0.8
  ) +
  
  scale_color_manual(
    
    values = my_colors,
    
    labels = c(
      "10" = "10",
      "5" = "5",
      "3" = "3"
    )
  ) +
  
  scale_linetype_manual(
    
    values = c(
      "10" = "solid",
      "5" = "solid",
      "3" = "solid"
    )
  ) +
  
  labs(
    
    color = "Number of plant species subsampled"
  ) +
  
  theme_void() +
  
  theme(
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    )
  )


legend <- cowplot::get_legend(
  legend_plot
)


############################################################
# 23. Combine NODF + Connectance + H2
############################################################

combined_plot <- plot_grid(
  
  connectance_row,
  
  nestedness_row,
  
  H2_row,
  
  legend,
  
  ncol = 1,
  
  rel_heights = c(
    1,
    1,
    1,
    0.12
  )
)


############################################################
# 24. Display
############################################################

combined_plot


############################################################
# 25. Save
############################################################

ggsave(
  
  "D:/Chap1_TargetPlant_to_monitor/result_260723/network.matrix_H2.png",
  
  combined_plot,
  
  width = 18,
  
  height = 13,
  
  units = "in",
  
  dpi = 300
)

