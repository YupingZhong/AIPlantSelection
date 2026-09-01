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


library(dplyr)
library(ggplot2)
library(cowplot)

# Paths and data
data_dir <- "data/processed"
result_dir <- "/Chap1_TargetPlant_to_monitor/result_260723"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

data_combined_all <- readRDS(file.path(data_dir, "H2_data_combined_all.rds"))
model_summary <- readRDS(file.path(data_dir, "H2_model_summary.rds"))
sensitivity_results <- readRDS(file.path(data_dir, "H2_sensitivity_results.rds"))

#Holm矫正
sensitivity_results <- sensitivity_results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "holm")
  )

# Shared panel style
theme_sensitivity <- theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "grey55", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(colour = "black", size = 10),
    legend.position = "bottom",
    plot.tag = element_text(
      face = "bold",
      size = 14,
      hjust = 0,
      vjust = 1,
      margin = margin(b = 4)
    ),
    plot.tag.position = "topleft",
    plot.margin = margin(8, 8, 5, 5)
  )

# (a) Abundance versus observed richness capture
r2_label <- paste0(
  "Marginal R² = ", round(model_summary$Marginal_R2, 2),
  "\nConditional R² = ", round(model_summary$Conditional_R2, 2)
)

fig1_relation <- ggplot(
  data_combined_all,
  aes(Abundance_scaled, Proportion_of_richness)
) +
  geom_point(alpha = 0.50, size = 1.8, colour = "#B07AA1") +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    level = 0.95,
    linewidth = 1,
    colour = "#76506F"
  ) +
  
  labs(
    tag = "(a)",
    x = "Standardized log flower abundance",
    y = "Pollinator richness captured (%)"
  ) +
  theme_sensitivity

# Classify performance for panel (b)
data_combined_all <- data_combined_all %>%
  mutate(
    performance_class = case_when(
      Abundance_scaled < quantile(Abundance_scaled, 0.05, na.rm = TRUE) &
        Residual_capture > quantile(Residual_capture, 0.90, na.rm = TRUE) ~
        "Rare and attractive",
      Residual_capture > 0 ~ "Above expected",
      TRUE ~ "Below expected"
    )
  )

# (b) Abundance versus deviation from expected richness capture
fig2_relation <- ggplot(
  data_combined_all,
  aes(Abundance_scaled, Residual_capture)
) +
  geom_point(aes(colour = performance_class), alpha = 0.70, size = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c(
      "Rare and attractive" = "#D97A6B",
      "Above expected" = "#C9A66B",
      "Below expected" = "#B8C4B2"
    )
  ) +
  labs(
    tag = "(b)",
    x = "Standardized log flower abundance",
    y = "Deviation from expected richness (%)",
    colour = NULL
  ) +
  theme_sensitivity


fig2_relation <- fig2_relation +
  guides(
    colour = guide_legend(ncol = 1)
  ) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.height = grid::unit(0.35, "cm")
  )

# Prepare threshold-sensitivity results for panel (c)
sensitivity_plot_data <- sensitivity_results %>%
  mutate(
    rare_q_label = factor(
      paste0("Rare < ", rare_q * 100, "%"),
      levels = c(
        "Rare < 5%",
        "Rare < 10%",
        "Rare < 15%"
      )
    ),
    
    attr_q_label = factor(
      paste0("Attractive > ", attr_q * 100, "%"),
      levels = c(
        "Attractive > 90%",
        "Attractive > 85%",
        "Attractive > 80%"
      )
    ),
    
    significant = p_adjusted < 0.05,
    
    p_label = case_when(
      is.na(p_adjusted) ~ "NA",
      p_adjusted < 0.001 ~ "<0.001",
      TRUE ~ paste0(
        "P = ",
        formatC(
          p_adjusted,
          format = "f",
          digits = 3
        )
      )
    )
  )

# (c) Sensitivity across rarity and attractiveness thresholds
sens_test <- ggplot(sensitivity_plot_data, aes(attr_q_label, rare_q_label)) +
  geom_tile(aes(fill = effect_size), colour = "white", linewidth = 0.5) +
  geom_text(aes(label = p_label), size = 3.5, colour = "black") +
  geom_tile(
    data = filter(sensitivity_plot_data, significant),
    aes(attr_q_label, rare_q_label),
    fill = NA,
    colour = "black",
    linewidth = 1,
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
    tag = "(c)",
    x = "Attractiveness threshold",
    y = "Rarity threshold"
  ) +
  theme_sensitivity

sens_test <- sens_test +
  scale_x_discrete(
    labels = c(
      "Attractive > 90%" = ">90th",
      "Attractive > 85%" = ">85th",
      "Attractive > 80%" = ">80th"
    )
  ) +
  scale_y_discrete(
    labels = c(
      "Rare < 5%" = "<5th",
      "Rare < 10%" = "<10th",
      "Rare < 15%" = "<15th"
    )
  ) +
  labs(
    x = "Attractiveness threshold\n(residual percentile)",
    y = "Rarity threshold\n(abundance percentile)"
  ) +
  theme(
    axis.text.x = element_text(size = 10, angle = 0),
    axis.title = element_text(size = 10, face = "bold")
  )

# Combine, display, and export
final_sensitivity_fig <- plot_grid(
  fig1_relation,
  fig2_relation,
  sens_test,
  ncol = 3,
  align = "h",
  axis = "tblr"
)

final_sensitivity_fig

ggsave(
  file.path(result_dir, "Fig_H2_sensitivity_three_panel.png"),
  final_sensitivity_fig,
  width = 13,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

write.csv(
  sensitivity_results,
  file.path(result_dir, "H2_sensitivity_results.csv"),
  row.names = FALSE
)

