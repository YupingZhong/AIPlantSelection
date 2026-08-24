#######################################################################
# 2.4 Plot main results: violin + paired networks + network size
#
# violin fill = subsampling strategy
# point colour = network size
# black diamond = mean
# black bar = median
#
# Statistical comparisons:
# paired Wilcoxon tests
#
# Top10:
#   Abundance vs Random
#   Abundance vs Phylogenetic
#
# Top5 / Top3:
#   Abundance vs Random
#   Abundance vs Phylogenetic
#   Abundance + shape vs Abundance
#######################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)
library(ggpubr)
library(flextable)
library(officer)


# ==========================================================
# Load data
# ==========================================================

data_count_scaled <- readRDS(
  "data/processed/data_count_scaled_published.rds"
)


result_all <- read.csv(
  "data/processed/result_all_published_PD.csv",
  header = TRUE
) 
# ==========================================================
# Calculate network size
# ==========================================================

network_size <- data_count_scaled %>%
  group_by(Study_Network_id) %>%
  summarise(
    n_plants_total = n_distinct(Plant_species),
    .groups = "drop"
  )

summary(network_size$n_plants_total)

range(
  network_size$n_plants_total,
  na.rm = TRUE
)


# ==========================================================
# Prepare plotting data
# ==========================================================

plot_data <- result_all %>%
  dplyr::select(
    Study_Network_id,
    
    # Abundance
    percentage_Abun10,
    percentage_Abun5,
    percentage_Abun3,
    
    # Abundance + shape
    percentage_FlwShape5,
    percentage_FlwShape3,
    
    # Phylogenetic
    percentage_Pylo10,
    percentage_Pylo5,
    percentage_Pylo3,
    
    # Random
    random_mean_percentage_Random10,
    random_mean_percentage_Random5,
    random_mean_percentage_Random3
  ) %>%
  
  pivot_longer(
    cols = -Study_Network_id,
    names_to = "Option",
    values_to = "Percentage"
  ) %>%
  
  filter(
    !is.na(Percentage)
  ) %>%
  
  mutate(
    
    # ------------------------------------------------------
    # Strategy
    # ------------------------------------------------------
    
    Strategy = case_when(
      
      grepl("Abun", Option) ~
        "Flower abundance",
      
      grepl("FlwShape", Option) ~
        "Flower abundance + shapes",
      
      grepl("Pylo", Option) ~
        "Phylogenetic distance",
      
      grepl("Random", Option) ~
        "Random"
    ),
    
    
    # ------------------------------------------------------
    # Sampling effort
    # ------------------------------------------------------
    
    Plant_Number = factor(
      case_when(
      grepl("10$", Option) ~ "Top10",
      grepl("5$", Option) ~ "Top5",
      grepl("3$", Option) ~ "Top3"
    ),
    levels = c("Top10", "Top5", "Top3")
  )
)%>%
  
  # --------------------------------------------------------
# Remove abundance + shape from Top10
# --------------------------------------------------------

filter(
  !(
    Plant_Number == "Top10" &
      Strategy == "Flower abundance + shapes"
  )
) %>%
  
  # --------------------------------------------------------
# Add network size
# --------------------------------------------------------

left_join(
  network_size,
  by = "Study_Network_id"
)


# ==========================================================
# Factor order
# ==========================================================

plot_data <- plot_data %>%
  mutate(
    
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
      case_when(
        grepl("10$", Option) ~ "Top10",
        grepl("5$", Option) ~ "Top5",
        grepl("3$", Option) ~ "Top3"
      ),
      levels = c("Top10", "Top5", "Top3")
    )
  )



# Check
print(levels(plot_data$Plant_Number))


# ==========================================================
# Calculate mean and median
# ==========================================================

plot_summary <- plot_data %>%
  group_by(
    Plant_Number,
    Strategy
  ) %>%
  summarise(
    
    mean = mean(
      Percentage,
      na.rm = TRUE
    ),
    
    median = median(
      Percentage,
      na.rm = TRUE
    ),
    
    n = n(),
    
    .groups = "drop"
  )


plot_summary$Strategy <- factor(
  plot_summary$Strategy,
  levels = levels(plot_data$Strategy)
)

plot_summary$Plant_Number <- factor(
  plot_summary$Plant_Number,
  levels = levels(plot_data$Plant_Number)
)


# ==========================================================
# Paired data
# ==========================================================

paired_data <- plot_data %>%
  dplyr::select(
    Study_Network_id,
    Plant_Number,
    Strategy,
    Percentage
  ) %>%
  
  pivot_wider(
    names_from = Strategy,
    values_from = Percentage
  )


# ==========================================================
# Function: paired Wilcoxon test
# ==========================================================

run_paired_test <- function(
    df,
    strategy1,
    strategy2
) {
  
  # If one strategy does not exist
  if (
    !strategy1 %in% names(df) |
    !strategy2 %in% names(df)
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  x <- df[[strategy1]]
  y <- df[[strategy2]]
  
  
  # Only complete pairs
  valid <- !is.na(x) & !is.na(y)
  
  
  if (sum(valid) < 3) {
    
    return(
      NA_real_
    )
  }
  
  
  test <- wilcox.test(
    x[valid],
    y[valid],
    paired = TRUE,
    exact = FALSE
  )
  
  
  return(
    test$p.value
  )
}


# ==========================================================
# Run paired tests
# ==========================================================

sig_results <- list()


for (
  n_level in levels(plot_data$Plant_Number)
) {
  
  df <- paired_data %>%
    filter(
      Plant_Number == n_level
    )
  
  
  # --------------------------------------------------------
  # Abundance vs Random
  # --------------------------------------------------------
  
  p1 <- run_paired_test(
    df,
    "Flower abundance",
    "Random"
  )
  
  
  # --------------------------------------------------------
  # Abundance vs Phylogenetic
  # --------------------------------------------------------
  
  p2 <- run_paired_test(
    df,
    "Flower abundance",
    "Phylogenetic distance"
  )
  
  
  # --------------------------------------------------------
  # Store results
  #
  # Top10 does NOT have abundance + shape
  # --------------------------------------------------------
  
  if (n_level == "Top10") {
    
    sig_results[[n_level]] <- data.frame(
      
      Plant_Number = n_level,
      
      comparison = c(
        "Abundance vs Random",
        "Abundance vs Phylogenetic"
      ),
      
      p = c(
        p1,
        p2
      )
    )
    
  } else {
    
    # ------------------------------------------------------
    # Abundance + shape vs Abundance
    # ------------------------------------------------------
    
    p3 <- run_paired_test(
      df,
      "Flower abundance + shapes",
      "Flower abundance"
    )
    
    
    sig_results[[n_level]] <- data.frame(
      
      Plant_Number = n_level,
      
      comparison = c(
        "Abundance vs Random",
        "Abundance vs Phylogenetic",
        "Abundance + shape vs Abundance"
      ),
      
      p = c(
        p1,
        p2,
        p3
      )
    )
  }
}


# ==========================================================
# Combine significance results
# ==========================================================

sig_results <- bind_rows(
  sig_results
)

sig_results$Plant_Number <- factor(
  sig_results$Plant_Number,
  levels = c("Top10", "Top5", "Top3")
)

# ==========================================================
# Define comparison groups
# ==========================================================

sig_results <- sig_results %>%
  mutate(
    
    # ------------------------------------------------------
    # Group 1
    # ------------------------------------------------------
    
    group1 = case_when(
      
      comparison ==
        "Abundance vs Random" ~
        "Random",
      
      comparison ==
        "Abundance vs Phylogenetic" ~
        "Phylogenetic distance",
      
      comparison ==
        "Abundance + shape vs Abundance" ~
        "Flower abundance"
    ),
    
    
    # ------------------------------------------------------
    # Group 2
    # ------------------------------------------------------
    
    group2 = case_when(
      
      comparison ==
        "Abundance vs Random" ~
        "Flower abundance",
      
      comparison ==
        "Abundance vs Phylogenetic" ~
        "Flower abundance",
      
      comparison ==
        "Abundance + shape vs Abundance" ~
        "Flower abundance + shapes"
    ),
    
    
    # ------------------------------------------------------
    # Significance label
    # ------------------------------------------------------
    
    label = case_when(
      
      is.na(p) ~
        "NA",
      
      p < 0.001 ~
        "***",
      
      p < 0.01 ~
        "**",
      
      p < 0.05 ~
        "*",
      
      TRUE ~
        "ns"
    )
  )


# ==========================================================
# Y positions for significance brackets
# ==========================================================
#
# Keep the brackets inside the 0-120 plotting range.
#
# Top10:
#   105 = Abundance vs Random
#   112 = Abundance vs Phylogenetic
#
# Top5 / Top3:
#   105 = Abundance vs Random
#   112 = Abundance vs Phylogenetic
#   119 = Abundance + shape vs Abundance
#
# ==========================================================

sig_results <- sig_results %>%
  mutate(
    
    y.position = case_when(
      
      comparison ==
        "Abundance vs Random" ~
        105,
      
      comparison ==
        "Abundance vs Phylogenetic" ~
        112,
      
      comparison ==
        "Abundance + shape vs Abundance" ~
        119
    )
  )


# ==========================================================
# Check significance results
# ==========================================================

print(sig_results)


# ==========================================================
# Strategy colours
# ==========================================================

strategy_colors <- c(
  
  "Flower abundance" =
    "#2AA889",
  
  "Flower abundance + shapes" =
    "#8E6BBE",
  
  "Phylogenetic distance" =
    "#E69F00",
  
  "Random" =
    "#B8B8B8"
)



# ==========================================================
# Plot
# ==========================================================

# Panel labels

panel_labels <- data.frame(
  Plant_Number = factor(
    c("Top10", "Top5", "Top3"),
    levels = c("Top10", "Top5", "Top3")
  ),
  label = c("(a)", "(b)", "(c)")
)


plot_violin_size <- ggplot(
  
  plot_data,
  
  aes(
    x = Strategy,
    y = Percentage
  )
) +
  
  
  # ========================================================
# Violin distributions
# ========================================================

geom_violin(
  
  aes(
    fill = Strategy
  ),
  
  width = 0.85,
  alpha = 0.55,
  color = "grey35",
  linewidth = 0.4,
  trim = TRUE
) +
  
  
  # ========================================================
# Connect the same network across strategies
# ========================================================

geom_line(
  
  aes(
    group = Study_Network_id
  ),
  
  color = "grey60",
  alpha = 0.08,
  linewidth = 0.2
) +
  
  
  # ========================================================
# Raw observations
# ========================================================

geom_point(
  
  aes(
    color = n_plants_total
  ),
  
  alpha = 0.50,
  size = 1.1
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
  
  width = 0.32,
  color = "black",
  linewidth = 0.45
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
  size = 4
) +
  

  # ========================================================
# Sampling-effort panels
# ========================================================
facet_grid(
  . ~ Plant_Number,
  scales = "free_x",
  space = "free_x"
) +
  
  # ========================================================
# Panel labels
# ========================================================

geom_text(
  data = panel_labels,
  aes(
    x = -Inf,
    y = Inf,
    label = label
  ),
  inherit.aes = FALSE,
  hjust = -0.3,
  vjust = -0.8,       
  fontface = "bold",
  size = 5
)  +
  

  
  # ========================================================
# Fill scale
# ========================================================

scale_fill_manual(values = strategy_colors, name = "Subsampling strategy") +
  
  # ========================================================
# Point colour = network size
# ========================================================

scale_color_viridis_c(
  option = "turbo",
  begin = 0.18,
  end = 0.90,
  name = "Flowering species per network",
  breaks = breaks_pretty(
    n = 5
  )
) +
  
  
  # ========================================================
# X-axis labels
# ========================================================

scale_x_discrete(
  labels = c(
    
    "Flower abundance" =
      "Abundance",
    
    "Flower abundance + shapes" =
      "Abundance + shape",
    
    "Phylogenetic distance" =
      "Phylogenetic",
    
    "Random" =
      "Random"
  ),
  drop = TRUE
) +
  
  # ========================================================
# Y-axis
# ========================================================

scale_y_continuous(
  
  limits = c(
    0,
    125
  ),
  
  breaks = seq(
    0,
    100,
    by = 20
  ),
  
  expand = expansion(
    mult = c(
      0.02,
      0.02
    )
  )
) +

  coord_cartesian(
    clip = "off"
  ) +
  
  # ========================================================
# Labels
# ========================================================

labs(
  
  x = "Subsampling strategy",
  
  y = "pollinator richness captured (%)"
) +
  
  
  # ========================================================
# Theme
# ========================================================

theme_classic(
  base_size = 12
) +
  
  theme(
    # Border around each facet panel
    panel.border = element_rect(
      colour = "grey75",
      fill = NA,
      linewidth = 0.7
    ),
    
    # Avoid overlapping theme_classic axis lines
    axis.line = element_blank(),
    
    # Space between facet panels
    panel.spacing = unit(
      0.8,
      "lines"
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.text.y = element_text(
      color = "black",
      size = 11
    ),
    
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      color = "black",
      size = 10
    ),
    
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "center",
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    # Extra top space for a, b and c
    plot.margin = margin(
      t = 22,
      r = 10,
      b = 10,
      l = 10
    )
  ) +
  
  
  # ========================================================
# Legends
# ========================================================

guides(
  
  fill = guide_legend(
    
    order = 1,
    
    override.aes = list(
      alpha = 0.7
    )
  ),
  
  color = guide_colorbar(
    
    order = 2,
    
    barwidth = 8,
    
    barheight = 0.8
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
  
  bracket.size = 0.4,
  
  size = 4,
  
  hide.ns = FALSE
)


# ==========================================================
# Display
# ==========================================================

plot_violin_size


ggsave(
  "/Chap1_TargetPlant_to_monitor/result_260723/main_violin_network_size.png",
  plot = plot_violin_size,
  width = 9,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ==========================================================
# Supplementary Table S1
# Pairwise Wilcoxon comparisons among subsampling strategies
# ==========================================================


# ==========================================================
# 1. Prepare paired comparison data
# ==========================================================

table_data <- result_all %>%
  dplyr::select(
    Study_Network_id,
    
    # Abundance
    Top10_Abundance = percentage_Abun10,
    Top5_Abundance  = percentage_Abun5,
    Top3_Abundance  = percentage_Abun3,
    
    # Abundance + flower shape
    Top5_Shape = percentage_FlwShape5,
    Top3_Shape = percentage_FlwShape3,
    
    # Phylogenetic
    Top10_Phylo = percentage_Pylo10,
    Top5_Phylo  = percentage_Pylo5,
    Top3_Phylo  = percentage_Pylo3,
    
    # Random
    Random10 = random_mean_percentage_Random10,
    Random5  = random_mean_percentage_Random5,
    Random3  = random_mean_percentage_Random3
  )


# ==========================================================
# 2. Function for paired Wilcoxon comparison
# ==========================================================

run_table_test <- function(
    data,
    var1,
    var2,
    comparison
) {
  
  x <- data[[var1]]
  y <- data[[var2]]
  
  valid <- !is.na(x) & !is.na(y)
  
  x <- x[valid]
  y <- y[valid]
  
  n <- length(x)
  
  if (n < 3) {
    
    return(
      tibble(
        Comparison = comparison,
        n = n,
        Median1 = NA_real_,
        Median2 = NA_real_,
        P_value = NA_real_
      )
    )
    
  }
  
  test <- wilcox.test(
    x,
    y,
    paired = TRUE,
    exact = FALSE
  )
  
  tibble(
    Comparison = comparison,
    n = n,
    Median1 = median(x, na.rm = TRUE),
    Median2 = median(y, na.rm = TRUE),
    P_value = test$p.value
  )
}


# ==========================================================
# 3. Run all Table S1 comparisons
# ==========================================================

table_s1 <- bind_rows(
  
  # --------------------------------------------------------
  # Within-strategy comparisons
  # --------------------------------------------------------
  
  run_table_test(
    table_data,
    "Top10_Abundance",
    "Top5_Abundance",
    "Top 10 abundant vs Top 5 abundant"
  ),
  
  run_table_test(
    table_data,
    "Top5_Abundance",
    "Top3_Abundance",
    "Top 5 abundant vs Top 3 abundant"
  ),
  
  run_table_test(
    table_data,
    "Top5_Shape",
    "Top3_Shape",
    "Top 5 flower-shape vs Top 3 flower-shape"
  ),
  
  run_table_test(
    table_data,
    "Random10",
    "Random5",
    "Random 10 vs Random 5"
  ),
  
  run_table_test(
    table_data,
    "Random5",
    "Random3",
    "Random 5 vs Random 3"
  ),
  
  # --------------------------------------------------------
  # Abundance vs Random
  # --------------------------------------------------------
  
  run_table_test(
    table_data,
    "Top10_Abundance",
    "Random10",
    "Top 10 abundant vs Random 10"
  ),
  
  run_table_test(
    table_data,
    "Top5_Abundance",
    "Random5",
    "Top 5 abundant vs Random 5"
  ),
  
  run_table_test(
    table_data,
    "Top3_Abundance",
    "Random3",
    "Top 3 abundant vs Random 3"
  ),
  
  # --------------------------------------------------------
  # Flower shape vs Random
  # --------------------------------------------------------
  
  run_table_test(
    table_data,
    "Top5_Shape",
    "Random5",
    "Top 5 flower-shape vs Random 5"
  ),
  
  run_table_test(
    table_data,
    "Top3_Shape",
    "Random3",
    "Top 3 flower-shape vs Random 3"
  ),
  
  # --------------------------------------------------------
  # Abundance vs Abundance + flower shape
  # --------------------------------------------------------
  
  run_table_test(
    table_data,
    "Top5_Abundance",
    "Top5_Shape",
    "Top 5 abundant vs Top 5 flower-shape"
  ),
  
  run_table_test(
    table_data,
    "Top3_Abundance",
    "Top3_Shape",
    "Top 3 abundant vs Top 3 flower-shape"
  )
)


# ==========================================================
# 4. Format P-values
# ==========================================================

table_s1 <- table_s1 %>%
  mutate(
    
    `Median1` = round(
      Median1,
      2
    ),
    
    `Median2` = round(
      Median2,
      2
    ),
    
    `P-value` = case_when(
      
      is.na(P_value) ~
        "NA",
      
      P_value < 0.001 ~
        "<0.001",
      
      TRUE ~
        format.pval(
          P_value,
          digits = 3,
          eps = 0.001
        )
    )
  ) %>%
  select(
    Comparison,
    n,
    Median1,
    Median2,
    `P-value`
  )


# ==========================================================
# 5. Print table in R
# ==========================================================

print(table_s1)


# ==========================================================
# 6. Export CSV
# ==========================================================

write.csv(
  table_s1,
  "/Chap1_TargetPlant_to_monitor/result_260723/Table_S1_Pairwise_Wilcoxon.csv",
  row.names = FALSE
)


# ==========================================================
# 7. Create Word table
# ==========================================================

ft_table_s1 <- flextable(
  table_s1
) %>%
  
  set_header_labels(
    Comparison = "Comparison",
    n = "n",
    Median1 = "Median1",
    Median2 = "Median2",
    `P-value` = "P-value"
  ) %>%
  
  theme_booktabs() %>%
  
  fontsize(
    size = 10,
    part = "all"
  ) %>%
  
  bold(
    part = "header"
  ) %>%
  
  align(
    j = "Comparison",
    align = "left",
    part = "all"
  ) %>%
  
  align(
    j = c(
      "n",
      "Median1",
      "Median2",
      "P-value"
    ),
    align = "center",
    part = "all"
  ) %>%
  
  width(
    j = "Comparison",
    width = 3.8
  ) %>%
  
  width(
    j = "n",
    width = 0.8
  ) %>%
  
  width(
    j = "Median1",
    width = 1.2
  ) %>%
  
  width(
    j = "Median2",
    width = 1.2
  ) %>%
  
  width(
    j = "P-value",
    width = 1.2
  )


# ==========================================================
# 8. Create Word document
# ==========================================================

doc_table_s1 <- read_docx() %>%
  
  body_add_par(
    "Table S1",
    style = "heading 2"
  ) %>%
  
  body_add_par(
    paste0(
      "Pairwise statistical comparisons among subsampling ",
      "strategies using paired Wilcoxon signed-rank tests. ",
      "Sample sizes and medians of each subsampling strategy ",
      "being compared are given, as well as p-values."
    )
  ) %>%
  
  body_add_flextable(
    ft_table_s1
  )


# ==========================================================
# 9. Save Word document
# ==========================================================

print(
  doc_table_s1,
  target =
    "/Chap1_TargetPlant_to_monitor/result_260723/Table_S1_Pairwise_Wilcoxon.docx"
)

