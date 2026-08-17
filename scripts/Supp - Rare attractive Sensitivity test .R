################################################################################
# H2 SENSITIVITY ANALYSIS
#
# This script ONLY:
#   1. Loads H2 analysis results
#   2. Creates sensitivity-analysis figures
#   3. Exports sensitivity-analysis table
#
# All statistical analyses are performed in Script 1.
################################################################################


# ==============================================================================
# 1. Libraries
# ==============================================================================

library(dplyr)
library(ggplot2)
library(cowplot)


# ==============================================================================
# 2. Paths
# ==============================================================================

data_dir <- "data/processed"

result_dir <- "/Chap1_TargetPlant_to_monitor/result_260723"

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. Load H2 results
# ==============================================================================

data_combined_all <- readRDS(
  file.path(
    data_dir,
    "H2_data_combined_all.rds"
  )
)

model_summary <- readRDS(
  file.path(
    data_dir,
    "H2_model_summary.rds"
  )
)

sensitivity_results <- readRDS(
  file.path(
    data_dir,
    "H2_sensitivity_results.rds"
  )
)


# ==============================================================================
# 4. Common theme
# ==============================================================================

theme_sensitivity <- theme_classic(
  base_size = 12
) +
  theme(
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    ),
    legend.position = "bottom",
    plot.margin = margin(
      5,
      5,
      5,
      5
    )
  )


# ==============================================================================
# 5. FIGURE 1
#    Abundance vs observed richness capture
# ==============================================================================

r2_label <- paste0(
  "Marginal R² = ",
  round(
    model_summary$Marginal_R2,
    2
  ),
  "\nConditional R² = ",
  round(
    model_summary$Conditional_R2,
    2
  )
)


fig1_relation <- ggplot(
  data_combined_all,
  aes(
    x = Abundance_scaled,
    y = Proportion_of_richness
  )
) +
  
  geom_point(
    alpha = 0.50,
    size = 1.8,
    color = "#8EC9D8"
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1,
    color = "#5C9FB2"
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = r2_label,
    hjust = 1.1,
    vjust = 1.5,
    size = 3.8
  ) +
  
  labs(
    x = "Standardized log flower abundance",
    y = "Percentage of pollinator richness captured (%)"
  ) +
  
  theme_sensitivity


# ==============================================================================
# 6. FIGURE 2
#    Abundance vs deviation from expected richness capture
# ==============================================================================

data_combined_all <- data_combined_all %>%
  mutate(
    
    performance_class = case_when(
      
      Abundance_scaled <
        quantile(
          Abundance_scaled,
          0.05,
          na.rm = TRUE
        ) &
        
        Residual_capture >
        quantile(
          Residual_capture,
          0.90,
          na.rm = TRUE
        )
      
      ~ "Rare and attractive",
      
      Residual_capture > 0
      ~ "Above expected",
      
      TRUE
      ~ "Below expected"
    )
  )


fig2_relation <- ggplot(
  data_combined_all,
  aes(
    x = Abundance_scaled,
    y = Residual_capture
  )
) +
  
  geom_point(
    aes(
      color = performance_class
    ),
    alpha = 0.70,
    size = 1.8
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  scale_color_manual(
    values = c(
      "Rare and attractive" = "#D97A6B",
      "Above expected" = "#C9A66B",
      "Below expected" = "#B8C4B2"
    )
  ) +
  
  labs(
    x = "Standardized log flower abundance",
    y = "Deviation from expected richness capture (%)",
    color = NULL
  ) +
  
  theme_sensitivity


# ==============================================================================
# 7. FIGURE 3
#    Sensitivity across rarity × attractiveness thresholds
#
#    Tile colour = effect size
#    Text         = P value
#    Black border = P < 0.05
# ==============================================================================

sensitivity_plot_data <- sensitivity_results %>%
  mutate(
    
    rare_q_label = factor(
      paste0(
        "Rare < ",
        rare_q * 100,
        "%"
      ),
      levels = c(
        "Rare < 5%",
        "Rare < 10%",
        "Rare < 15%"
      )
    ),
    
    attr_q_label = factor(
      paste0(
        "Attractive > ",
        attr_q * 100,
        "%"
      ),
      levels = c(
        "Attractive > 90%",
        "Attractive > 85%",
        "Attractive > 80%"
      )
    ),
    
    significant = p_value < 0.05,
    
    p_label = format.pval(
      p_value,
      digits = 2,
      eps = 0.001
    )
  )


sens_test <- ggplot(
  sensitivity_plot_data,
  aes(
    x = attr_q_label,
    y = rare_q_label
  )
) +
  
  # Effect size
  geom_tile(
    aes(
      fill = effect_size
    ),
    color = "white",
    linewidth = 0.5
  ) +
  
  # P value
  geom_text(
    aes(
      label = p_label
    ),
    size = 3.5,
    color = "black"
  ) +
  
  # Significant results
  geom_tile(
    data = sensitivity_plot_data %>%
      filter(significant),
    aes(
      x = attr_q_label,
      y = rare_q_label
    ),
    fill = NA,
    color = "black",
    linewidth = 1.0,
    inherit.aes = FALSE
  ) +
  
  scale_fill_gradient2(
    low = "#90D4A4",
    mid = "#FFFFFF",
    high = "#D97A6B",
    midpoint = 0,
    name = "Effect size\n(median difference)"
  ) +
  
  labs(
    x = "Attractiveness threshold",
    y = "Rarity threshold"
  ) +
  
  theme_sensitivity


# ==============================================================================
# 8. Combined sensitivity figure
# ==============================================================================

final_sensitivity_fig <- plot_grid(
  
  fig1_relation,
  fig2_relation,
  sens_test,
  
  ncol = 3,
  
  labels = c(
    "a",
    "b",
    "c"
  ),
  
  label_size = 14,
  label_fontface = "bold",
  
  align = "h"
)


print(final_sensitivity_fig)


# ==============================================================================
# 9. Export figures
# ==============================================================================

ggsave(
  file.path(
    result_dir,
    "Fig_H2_sensitivity_three_panel.png"
  ),
  final_sensitivity_fig,
  width = 13,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ==============================================================================
# 10. Export sensitivity table
# ==============================================================================

write.csv(
  sensitivity_results,
  file.path(
    result_dir,
    "H2_sensitivity_results.csv"
  ),
  row.names = FALSE
)


