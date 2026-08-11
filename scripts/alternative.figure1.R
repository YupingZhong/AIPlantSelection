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


# ==========================================================
# Load data
# ==========================================================

data_count_scaled <- readRDS(
  "data/processed/data_count_scaled_published.rds"
)

Kal_id_list <- data_count_scaled %>%
  select(
    Study_id,
    Study_Network_id
  ) %>%
  filter(
    Study_id == "22_Kallnik"
  ) %>%
  distinct()

result_all <- read.csv(
  "data/processed/result_all_published_PD.csv",
  header = TRUE
) %>%
  filter(
    !Study_Network_id %in% Kal_id_list$Study_Network_id
  )

# 检查是否还有 Kallnik
result_all %>%
  filter(
    Study_Network_id %in% Kal_id_list$Study_Network_id
  )

result_all %>%
  summarise(
    n_sample = n_distinct(Study_Network_id),
    .groups = "drop"
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
  select(
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
    
    Plant_Number = case_when(
      grepl("10$", Option) ~ "Top10",
      grepl("5$", Option) ~ "Top5",
      grepl("3$", Option) ~ "Top3"
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
      Plant_Number,
      levels = c(
        "Top10",
        "Top5",
        "Top3"
      )
    )
  )



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
  select(
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

facet_wrap(
  ~ Plant_Number,
  nrow = 1,
  ncol = 3,
  scales = "free_x",
  drop = FALSE
) +
  
  
  # ========================================================
# Fill scale
# ========================================================

scale_fill_manual(
  
  values = strategy_colors,
  
  name = "Subsampling strategy"
) +
  
  
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

theme_classic(
  base_size = 12
) +
  
  theme(
    
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
    
    plot.margin = margin(
      10,
      10,
      10,
      10
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

plot_data %>%
  group_by(Plant_Number) %>%
  summarise(
    n_sample = n_distinct(Study_Network_id),
    .groups = "drop"
  )

