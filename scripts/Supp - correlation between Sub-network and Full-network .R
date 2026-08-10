
############################################################

##----3.5 plot Relationship of the richness of pollinators of subsampled networks ---------------
#         with the richness of the full network for each subsampling strategy


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)

# ============================================================
# 1. 读取数据
# ============================================================

result_all <- read.csv(
  "data/processed/result_all_published_PD.csv",
  header = TRUE
)


# ============================================================
# 2. 基础颜色和主题
# ============================================================

my_colors <- c(
  "Top 10" = "#66C2A5",
  "Top 5"  = "#B07AA1",
  "Top 3"  = "#E9C46A",
  "flw 5"  = "#B07AA1",
  "flw 3"  = "#E9C46A",
  "Phylo 10" = "#66C2A5",
  "Phylo 5"  = "#B07AA1",
  "Phylo 3"  = "#E9C46A",
  "Random 10" = "#66C2A5",
  "Random 5"  = "#B07AA1",
  "Random 3" = "#E9C46A"
)

base_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15
    ),
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    legend.position = "none"
  )


# ============================================================
# 3. 计算每个 subsampling strategy 的 Spearman p-value
# ============================================================

get_linetype <- function(df) {
  
  df %>%
    group_by(Subsampling) %>%
    summarise(
      p = {
        tmp <- cur_data()
        
        tmp <- tmp[
          complete.cases(
            tmp$Full_Richness,
            tmp$Subsampled_Richness
          ),
        ]
        
        if (nrow(tmp) < 3) {
          NA_real_
        } else {
          cor.test(
            tmp$Full_Richness,
            tmp$Subsampled_Richness,
            method = "spearman"
          )$p.value
        }
      },
      .groups = "drop"
    ) %>%
    mutate(
      linetype = ifelse(
        p < 0.05,
        "solid",
        "dashed"
      )
    )
}


# ============================================================
# 4. 把四种方法全部整理成统一的 long format
# ============================================================

rich_abun_long <- result_all %>%
  select(
    Study_Network_id,
    Full_Richness = total_pollinator_count_Abun10,
    percentage_Abun10,
    percentage_Abun5,
    percentage_Abun3
  ) %>%
  pivot_longer(
    cols = starts_with("percentage_Abun"),
    names_to = "Subsampling",
    values_to = "Subsampled_Richness"
  ) %>%
  mutate(
    Subsampling = factor(
      Subsampling,
      levels = c(
        "percentage_Abun10",
        "percentage_Abun5",
        "percentage_Abun3"
      ),
      labels = c(
        "Top 10",
        "Top 5",
        "Top 3"
      )
    )
  )


rich_shape_long <- result_all %>%
  select(
    Study_Network_id,
    Full_Richness = total_pollinator_count_Abun10,
    percentage_FlwShape5,
    percentage_FlwShape3
  ) %>%
  pivot_longer(
    cols = starts_with("percentage_FlwShape"),
    names_to = "Subsampling",
    values_to = "Subsampled_Richness"
  ) %>%
  mutate(
    Subsampling = factor(
      Subsampling,
      levels = c(
        "percentage_FlwShape5",
        "percentage_FlwShape3"
      ),
      labels = c(
        "flw 5",
        "flw 3"
      )
    )
  )


rich_phylo_long <- result_all %>%
  select(
    Study_Network_id,
    Full_Richness = total_pollinator_count_Abun10,
    percentage_Pylo10,
    percentage_Pylo5,
    percentage_Pylo3
  ) %>%
  pivot_longer(
    cols = starts_with("percentage_Pylo"),
    names_to = "Subsampling",
    values_to = "Subsampled_Richness"
  ) %>%
  mutate(
    Subsampling = factor(
      Subsampling,
      levels = c(
        "percentage_Pylo10",
        "percentage_Pylo5",
        "percentage_Pylo3"
      ),
      labels = c(
        "Phylo 10",
        "Phylo 5",
        "Phylo 3"
      )
    )
  )


rich_random_long <- result_all %>%
  select(
    Study_Network_id,
    Full_Richness = total_pollinator_count_Abun10,
    random_mean_percentage_Random10,
    random_mean_percentage_Random5,
    random_mean_percentage_Random3
  ) %>%
  pivot_longer(
    cols = starts_with("random_mean_percentage"),
    names_to = "Subsampling",
    values_to = "Subsampled_Richness"
  ) %>%
  mutate(
    Subsampling = factor(
      Subsampling,
      levels = c(
        "random_mean_percentage_Random10",
        "random_mean_percentage_Random5",
        "random_mean_percentage_Random3"
      ),
      labels = c(
        "Random 10",
        "Random 5",
        "Random 3"
      )
    )
  )


# ============================================================
# 5. 分别计算四组 regression line 的 linetype
# ============================================================

linetype_abun <- get_linetype(rich_abun_long)
linetype_shape <- get_linetype(rich_shape_long)
linetype_phylo <- get_linetype(rich_phylo_long)
linetype_random <- get_linetype(rich_random_long)


# ============================================================
# 6. 统一绘图函数
# ============================================================

plot_rich <- function(
    df,
    linetype_df,
    title,
    cor_label_x = "left",
    cor_label_y = NULL
) {
  
  n_group <- length(unique(df$Subsampling))
  
  if (is.null(cor_label_y)) {
    cor_label_y <- seq(
      0.98,
      0.98 - 0.02 * (n_group - 1),
      length.out = n_group
    )
  }
  
  ggplot(
    df,
    aes(
      x = Full_Richness,
      y = Subsampled_Richness,
      color = Subsampling
    )
  ) +
    
    # points
    geom_point(
      shape = 1,
      size = 2,
      position = position_jitter(
        width = 0.03,
        height = 0.5
      )
    ) +
    
    # regression
    geom_smooth(
      method = "lm",
      se = FALSE,
      aes(
        linetype = Subsampling
      ),
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    
    # solid / dashed according to Spearman p
    scale_linetype_manual(
      values = setNames(
        linetype_df$linetype,
        linetype_df$Subsampling
      )
    ) +
    
    # Spearman correlation
    stat_cor(
      aes(
        group = Subsampling,
        color = Subsampling
      ),
      method = "spearman",
      size = 3.5,
      label.x.npc = cor_label_x,
      label.y.npc = cor_label_y
    ) +
    
    scale_color_manual(
      values = my_colors
    ) +
    
    labs(
      x = "Pollinator richness in full network",
      y = NULL,
      title = title
    ) +
    
    base_theme +
    
    theme(
      plot.title = element_text(
        size = 11,
        hjust = 0.5
      ),
      axis.title.x = element_text(
        size = 12
      ),
      axis.text = element_text(
        size = 11
      ),
      plot.margin = margin(
        t = 0,
        r = 10,
        b = 0,
        l = 0
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 150)
    )
}


# ============================================================
# 7. 四张图
# ============================================================

p_abun <- plot_rich(
  rich_abun_long,
  linetype_abun,
  title = "a. Abundance-based"
)

p_shape <- plot_rich(
  rich_shape_long,
  linetype_shape,
  title = "b. Abundance + flower shape"
)

p_phylo <- plot_rich(
  rich_phylo_long,
  linetype_phylo,
  title = "c. Phylogenetic"
)

p_random <- plot_rich(
  rich_random_long,
  linetype_random,
  title = "d. Random",
  cor_label_x = 0.33
)


# ============================================================
# 8. 单独提取 legend
# ============================================================

legend_plot <- p_abun +
  aes(color = Subsampling) +
  labs(
    colour = "Number of plant species subsampled"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

legend <- cowplot::get_legend(
  legend_plot
)


# ============================================================
# 9. 四张图横向排列
# ============================================================

plots_row <- plot_grid(
  p_abun,
  p_shape,
  p_phylo,
  p_random,
  ncol = 4,
  align = "h",
  axis = "tb",
  rel_widths = c(1, 1, 1, 1)
)


# ============================================================
# 10. 加 legend
# ============================================================

final_plot <- plot_grid(
  plots_row,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)


# ============================================================
# 11. 添加共同 Y 轴标题
# ============================================================

final_plot <- ggdraw(final_plot) +
  draw_label(
    "Pollinator richness in subsampled network",
    x = 0,
    y = 0.5,
    angle = 90,
    fontface = "bold",
    size = 12,
    hjust = 0.41,
    vjust = 1.5
  )


# ============================================================
# 12. 显示
# ============================================================

final_plot


# ============================================================
# 13. 保存
# ============================================================

ggsave(
  "/Chap1_TargetPlant_to_monitor/result_260526/rich_4.png",
  final_plot,
  width = 12,
  height = 5,
  units = "in",
  dpi = 300
)







############################################################
## Visualization of nestedness results
############################################################

#Lines represent linear regression fits for visualization.
#ρ and p indicate Spearman rank correlations and their significance.
#Solid lines indicate significant Spearman correlations (p < 0.05), whereas dashed lines indicate non-significant correlations. 


############################################################
## Network metrics: Nestedness + Connectance
## Four subsampling strategies:
## a. Abundance-based
## b. Abundance + flower shape
## c. Phylogenetic
## d. Random
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

base_theme <- theme_classic(base_size = 12) +
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
# 4. Function: calculate significance
############################################################

get_linetype <- function(
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
      
      p = {
        
        tmp <- cur_data()
        
        tmp <- tmp[
          complete.cases(
            tmp[, c(x_col, y_col)]
          ),
        ]
        
        if (nrow(tmp) < 3) {
          
          NA_real_
          
        } else {
          
          cor.test(
            tmp[[x_col]],
            tmp[[y_col]],
            method = "spearman"
          )$p.value
          
        }
      },
      
      .groups = "drop"
      
    ) %>%
    mutate(
      
      linetype = ifelse(
        p < 0.05,
        "solid",
        "dashed"
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
  
  # number of groups
  n_group <- length(
    unique(df[[group_col]])
  )
  
  # correlation label positions
  if (is.null(cor_label_y)) {
    
    cor_label_y <- seq(
      0.98,
      0.98 - 0.04 * (n_group - 1),
      length.out = n_group
    )
    
  }
  
  
  ggplot(
    df,
    aes(
      x = .data[[x_col]],
      y = .data[[y_col]],
      color = .data[[group_col]]
    )
  ) +
    
    # points
    geom_point(
      shape = 1,
      size = 2,
      position = position_jitter(
        width = jitter_width,
        height = jitter_height
      )
    ) +
    
    # regression
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.8,
      aes(
        linetype = .data[[group_col]]
      ),
      show.legend = FALSE
    ) +
    
    # significant = solid
    # non-significant = dashed
    scale_linetype_manual(
      values = setNames(
        linetype_df$linetype,
        linetype_df[[group_col]]
      )
    ) +
    
    # Spearman correlation
    stat_cor(
      aes(
        group = .data[[group_col]],
        color = .data[[group_col]]
      ),
      method = "spearman",
      size = 3.5,
      label.x.npc = cor_label_x,
      label.y.npc = cor_label_y
    ) +
    
    scale_color_manual(
      values = my_colors
    ) +
    
    scale_y_continuous(
      limits = y_limits
    ) +
    
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    
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
  
  select(
    Study_Network_id,
    Nestedness_orig,
    Nestedness_Abun10,
    Nestedness_Abun5,
    Nestedness_Abun3
  ) %>%
  
  pivot_longer(
    cols = starts_with("Nestedness_Abun"),
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


nested_abun_sig <- get_linetype(
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
  y_limits = c(0, 60),
  jitter_width = 0.03,
  jitter_height = 0.5
)


############################################################
# 8. Nestedness - Flower shape
############################################################

nested_flw <- all_nestedness_df %>%
  
  select(
    Study_Network_id,
    Nestedness_orig,
    Nestedness_FlwShape5,
    Nestedness_FlwShape3
  ) %>%
  
  pivot_longer(
    cols = starts_with("Nestedness_FlwShape"),
    names_to = "Level",
    values_to = "Value"
  ) %>%
  
  mutate(
    Level = case_when(
      grepl("5$", Level) ~ "5",
      grepl("3$", Level) ~ "3"
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


nested_flw_sig <- get_linetype(
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
  y_limits = c(0, 60),
  jitter_width = 0.03,
  jitter_height = 0.5
)


############################################################
# 9. Nestedness - Phylogenetic
############################################################

nested_phylo <- all_nestedness_df %>%
  
  select(
    Study_Network_id,
    Nestedness_orig,
    Nestedness_Pylo10,
    Nestedness_Pylo5,
    Nestedness_Pylo3
  ) %>%
  
  pivot_longer(
    cols = starts_with("Nestedness_Pylo"),
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


nested_phylo_sig <- get_linetype(
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
  y_limits = c(0, 60),
  jitter_width = 0.03,
  jitter_height = 0.5
)


############################################################
# 10. Nestedness - Random
############################################################

random_nodf <- random_all %>%
  
  select(
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
      select(
        Study_Network_id,
        Nestedness_orig
      ),
    by = "Study_Network_id"
  ) %>%
  
  filter(
    !is.na(Value)
  )


random_nodf_sig <- get_linetype(
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
  y_limits = c(0, 60),
  jitter_width = 0.02,
  jitter_height = 0.2,
  cor_label_x = 0.33
)


############################################################

# CONNECTANCE

############################################################


############################################################
# 11. Connectance - Abundance
############################################################

conn_abun <- all_connectance_df %>%
  
  select(
    Study_Network_id,
    Connectance_orig,
    Connectance_Abun10,
    Connectance_Abun5,
    Connectance_Abun3
  ) %>%
  
  pivot_longer(
    cols = starts_with("Connectance_Abun"),
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


conn_abun_sig <- get_linetype(
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
  y_limits = c(0, 1),
  jitter_width = 0.005,
  jitter_height = 0.005
)


############################################################
# 12. Connectance - Flower shape
############################################################

conn_flw <- all_connectance_df %>%
  
  select(
    Study_Network_id,
    Connectance_orig,
    Connectance_FlwShape5,
    Connectance_FlwShape3
  ) %>%
  
  pivot_longer(
    cols = starts_with("Connectance_FlwShape"),
    names_to = "Level",
    values_to = "Value"
  ) %>%
  
  mutate(
    Level = case_when(
      grepl("5$", Level) ~ "5",
      grepl("3$", Level) ~ "3"
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


conn_flw_sig <- get_linetype(
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
  y_limits = c(0, 1),
  jitter_width = 0.005,
  jitter_height = 0.005
)


############################################################
# 13. Connectance - Phylogenetic
############################################################

conn_phylo <- all_connectance_df %>%
  
  select(
    Study_Network_id,
    Connectance_orig,
    Connectance_Pylo10,
    Connectance_Pylo5,
    Connectance_Pylo3
  ) %>%
  
  pivot_longer(
    cols = starts_with("Connectance_Pylo"),
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


conn_phylo_sig <- get_linetype(
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
  y_limits = c(0, 1),
  jitter_width = 0.005,
  jitter_height = 0.005
)


############################################################
# 14. Connectance - Random
############################################################

random_conn <- random_all %>%
  
  select(
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
      select(
        Study_Network_id,
        Connectance_orig
      ),
    by = "Study_Network_id"
  ) %>%
  
  filter(
    !is.na(Value)
  )


random_conn_sig <- get_linetype(
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
  y_limits = c(0, 1),
  jitter_width = 0.005,
  jitter_height = 0.005,
  cor_label_x = 0.33
)


############################################################
# ==========================================================
# 15. Remove legends
# ==========================================================
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
# 16. Create a clean unified legend
############################################################

legend_df <- data.frame(
  Level = factor(
    c("10", "5", "3"),
    levels = c("10", "5", "3")
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
# 17. Combine Nestedness + Connectance
############################################################

combined_plot <- plot_grid(
  
  connectance_row,
  nestedness_row,
  legend,
  
  ncol = 1,
  
  rel_heights = c(
    1,
    1,
    0.12
  )
)


############################################################
# 18. Display
############################################################

combined_plot


############################################################
# 19. Save
############################################################

ggsave(
  "D:/Chap1_TargetPlant_to_monitor/result_260723/network.matrix.png",
  combined_plot,
  width = 18,
  height = 9.5,
  units = "in",
  dpi = 300
)

