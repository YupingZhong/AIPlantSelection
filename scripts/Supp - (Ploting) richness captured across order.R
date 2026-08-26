#  Supplement richness captured across order

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)
library(ggpubr)

# ==========================================================
# Load corrected results
# ==========================================================

result_all_order <- readRDS("data/processed/result_all_by_pollinator_order.rds")


# ==========================================================
# Network size
# ==========================================================

network_size <- data_count_scaled %>%
  group_by(Study_Network_id) %>%
  summarise(
    n_plants_total = n_distinct(Plant_species),
    .groups = "drop"
  )

# ==========================================================
# Prepare plotting data - CRITICAL: Set factor levels early
# ==========================================================

plot_data <- result_all_order %>%
  mutate(
    
    Strategy = case_when(
      strategy == "Abundance" ~
        "Flower abundance",
      
      strategy == "Flower shape" ~
        "Flower abundance + shapes",
      
      strategy == "Phylogenetic" ~
        "Phylogenetic distance",
      
      strategy == "Random" ~
        "Random"
    ),
    
    # 立即转为因子并设置顺序
    Strategy = factor(
      Strategy,
      levels = c(
        "Flower abundance",
        "Flower abundance + shapes",
        "Phylogenetic distance",
        "Random"
      )
    ),
    
    Plant_Number = factor(
      paste0("Top", n_plants),
      levels = c("Top10", "Top5", "Top3")
    ),
    
    pollinator_order = factor(
      pollinator_order,
      levels = c(
        "Diptera",
        "Hymenoptera",
        "Lepidoptera",
        "Coleoptera"
      )
    )
  ) %>%
  
  filter(!is.na(percentage)) %>%
  
  # Flower shape has no Top10
  filter(
    !(
      Plant_Number == "Top10" &
        Strategy == "Flower abundance + shapes"
    )
  ) %>%
  
  left_join(network_size, by = "Study_Network_id")


# 验证因子顺序
print("Plant_Number levels:")
print(levels(plot_data$Plant_Number))
print("pollinator_order levels:")
print(levels(plot_data$pollinator_order))


# ==========================================================
# Calculate summary statistics
# ==========================================================

plot_summary <- plot_data %>%
  group_by(
    pollinator_order,
    Plant_Number,
    Strategy
  ) %>%
  summarise(
    mean = mean(percentage, na.rm = TRUE),
    median = median(percentage, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  # 重新确保因子顺序
  mutate(
    Strategy = factor(Strategy, 
                      levels = c("Flower abundance", 
                                 "Flower abundance + shapes",
                                 "Phylogenetic distance", 
                                 "Random")),
    Plant_Number = factor(Plant_Number, 
                          levels = c("Top10", "Top5", "Top3")),
    pollinator_order = factor(pollinator_order,
                              levels = c("Diptera", "Hymenoptera", 
                                         "Lepidoptera", "Coleoptera"))
  )


# ==========================================================
# Prepare paired data for statistical tests
# ==========================================================

paired_data <- plot_data %>%
  dplyr::select(
    Study_Network_id,
    pollinator_order,
    Plant_Number,
    Strategy,
    percentage
  ) %>%
  pivot_wider(
    names_from = Strategy,
    values_from = percentage
  ) %>%
  # 重新确保因子顺序
  mutate(
    Plant_Number = factor(Plant_Number, 
                          levels = c("Top10", "Top5", "Top3")),
    pollinator_order = factor(pollinator_order,
                              levels = c("Diptera", "Hymenoptera", 
                                         "Lepidoptera", "Coleoptera"))
  )


# ==========================================================
# Function: paired Wilcoxon test
# ==========================================================

run_paired_test <- function(df, strategy1, strategy2) {
  
  if (!strategy1 %in% names(df) | !strategy2 %in% names(df)) {
    return(NA_real_)
  }
  
  x <- df[[strategy1]]
  y <- df[[strategy2]]
  
  valid <- !is.na(x) & !is.na(y)
  
  if (sum(valid) < 3) {
    return(NA_real_)
  }
  
  wilcox.test(
    x[valid],
    y[valid],
    paired = TRUE,
    exact = FALSE
  )$p.value
}


# ==========================================================
# Run statistical tests
# ==========================================================

sig_results <- list()
counter <- 1

for (order_level in levels(plot_data$pollinator_order)) {
  
  for (n_level in levels(plot_data$Plant_Number)) {
    
    df <- paired_data %>%
      filter(
        pollinator_order == order_level,
        Plant_Number == n_level
      )
    
    p1 <- run_paired_test(
      df,
      "Flower abundance",
      "Random"
    )
    
    p2 <- run_paired_test(
      df,
      "Flower abundance",
      "Phylogenetic distance"
    )
    
    if (n_level == "Top10") {
      
      sig_results[[counter]] <- tibble(
        pollinator_order = order_level,
        Plant_Number = n_level,
        comparison = c(
          "Abundance vs Random",
          "Abundance vs Phylogenetic"
        ),
        p = c(p1, p2)
      )
      
    } else {
      
      p3 <- run_paired_test(
        df,
        "Flower abundance + shapes",
        "Flower abundance"
      )
      
      sig_results[[counter]] <- tibble(
        pollinator_order = order_level,
        Plant_Number = n_level,
        comparison = c(
          "Abundance vs Random",
          "Abundance vs Phylogenetic",
          "Abundance + shape vs Abundance"
        ),
        p = c(p1, p2, p3)
      )
    }
    
    counter <- counter + 1
  }
}

sig_results <- bind_rows(sig_results) %>%
  # 重新确保因子顺序
  mutate(
    Plant_Number = factor(Plant_Number, 
                          levels = c("Top10", "Top5", "Top3")),
    pollinator_order = factor(pollinator_order,
                              levels = c("Diptera", "Hymenoptera", 
                                         "Lepidoptera", "Coleoptera"))
  )


# ==========================================================
# Format significance results for plotting
# ==========================================================

sig_results <- sig_results %>%
  mutate(
    
    group1 = case_when(
      comparison == "Abundance vs Random" ~
        "Flower abundance",
      
      comparison == "Abundance vs Phylogenetic" ~
        "Flower abundance",
      
      comparison == "Abundance + shape vs Abundance" ~
        "Flower abundance"
    ),
    
    group2 = case_when(
      comparison == "Abundance vs Random" ~
        "Random",
      
      comparison == "Abundance vs Phylogenetic" ~
        "Phylogenetic distance",
      
      comparison == "Abundance + shape vs Abundance" ~
        "Flower abundance + shapes"
    ),
    
    label = case_when(
      is.na(p) ~ "NA",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )


# ==========================================================
# Y positions for significance brackets
# ==========================================================

panel_max <- plot_data %>%
  group_by(
    pollinator_order,
    Plant_Number
  ) %>%
  summarise(
    ymax = max(percentage, na.rm = TRUE),
    .groups = "drop"
  )

sig_results <- sig_results %>%
  left_join(
    panel_max,
    by = c("pollinator_order", "Plant_Number")
  ) %>%
  mutate(
    y.position = case_when(
      comparison == "Abundance vs Random" ~
        pmin(ymax + 8, 130),
      
      comparison == "Abundance vs Phylogenetic" ~
        pmin(ymax + 15, 130),
      
      comparison == "Abundance + shape vs Abundance" ~
        pmin(ymax + 22, 130)
    )
  )


# ==========================================================
# Strategy colours
# ==========================================================

strategy_colors <- c(
  "Flower abundance" = "#2AA889",
  "Flower abundance + shapes" = "#8E6BBE",
  "Phylogenetic distance" = "#E69F00",
  "Random" = "#B8B8B8"
)


# ==========================================================
# Supplementary Figure
# Flower-selection strategies by pollinator order
# ==========================================================

plot_supp_order <- ggplot(
  plot_data,
  aes(x = Strategy, y = percentage)
) +
  
  # ========================================================
# Violin
# ========================================================
geom_violin(
  aes(fill = Strategy),
  width = 0.80,
  alpha = 0.55,
  color = "grey35",
  linewidth = 0.35,
  trim = TRUE
) +
  
  # ========================================================
# Connect networks
# ========================================================
geom_line(
  aes(group = Study_Network_id),
  color = "grey60",
  alpha = 0.08,
  linewidth = 0.20
) +
  
  # ========================================================
# Individual networks
# ========================================================
geom_point(
  aes(color = n_plants_total),
  alpha = 0.45,
  size = 0.9
) +
  
  # ========================================================
# Median
# ========================================================
geom_crossbar(
  data = plot_summary,
  aes(
    x = Strategy,
    y = median,
    ymin = median,
    ymax = median
  ),
  inherit.aes = FALSE,
  width = 0.28,
  color = "black",
  linewidth = 0.40
) +
  
  # ========================================================
# Mean
# ========================================================
geom_point(
  data = plot_summary,
  aes(
    x = Strategy,
    y = mean
  ),
  inherit.aes = FALSE,
  color = "black",
  shape = 18,
  size = 3.0
) +
  
  # ========================================================
# Facets - 注意：行是pollinator_order，列是Plant_Number
# ========================================================
facet_grid(
  pollinator_order ~ Plant_Number,
  scales = "free_x",
  space = "free_x"
) +
  
  # ========================================================
# Strategy colours
# ========================================================
scale_fill_manual(
  values = strategy_colors,
  name = "Subsampling strategy"
) +
  
  # ========================================================
# Network size
# ========================================================
scale_color_viridis_c(
  option = "turbo",
  begin = 0.18,
  end = 0.90,
  name = "Flowering species\nper network",
  breaks = breaks_pretty(n = 5)
) +
  
  # ========================================================
# X labels
# ========================================================
scale_x_discrete(
  labels = c(
    "Flower abundance" = "Abundance",
    "Flower abundance + shapes" = "Abundance + shape",
    "Phylogenetic distance" = "Phylogenetic",
    "Random" = "Random"
  ),
  drop = TRUE
) +
  
  # ========================================================
# Y axis
# ========================================================
scale_y_continuous(
  limits = c(0, 135),
  breaks = seq(0, 100, by = 20),
  expand = expansion(mult = c(0.01, 0.01))
) +
  
  # ========================================================
# Labels
# ========================================================
labs(
  x = "Subsampling strategy",
  y = "Percent of pollinator richness captured"
) +
  
  # ========================================================
# Theme
# ========================================================
theme_classic(base_size = 11) +
  
  theme(
    panel.border = element_rect(
      colour = "grey55",
      fill = NA,
      linewidth = 0.6
    ),
    
    # Avoid overlapping theme_classic axis lines
    axis.line = element_blank(),
    
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text.y = element_text(
      color = "black",
      size = 9
    ),
    
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      color = "black",
      size = 8
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.position = "bottom",
    
    legend.box = "vertical",
    
    legend.title = element_text(
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  # ========================================================
# Legends
# ========================================================
guides(
  fill = guide_legend(
    order = 1,
    override.aes = list(alpha = 0.7)
  ),
  
  color = guide_colorbar(
    order = 2,
    barwidth = 8,
    barheight = 0.7
  )
) +
  
  # ========================================================
# Significance annotations
# ========================================================
stat_pvalue_manual(
  sig_results,
  label = "label",
  xmin = "group1",
  xmax = "group2",
  y.position = "y.position",
  tip.length = 0.01,
  bracket.size = 0.3,
  size = 3,
  hide.ns = FALSE
)

plot_supp_order

ggsave(
  filename = "/Chap1_TargetPlant_to_monitor/result_260723/Supplementary_Figure_pollinator_order.png",
  plot = plot_supp_order,
  width = 10,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "/Chap1_TargetPlant_to_monitor/result_260723/Supplementary_Figure_pollinator_order.pdf",
  plot = plot_supp_order,
  width = 10,
  height = 9,
  units = "in",
  device = cairo_pdf
)
