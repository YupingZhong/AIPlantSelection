#######################################################################
# 2.4 Plot main results: violin + paired networks + network size
#     violin fill = subsampling strategy
#     point colour = network size
#     black diamond = mean
#     black bar = median
#######################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(scales)

# ==========================================================
# Load data
# ==========================================================

result_all <- read.csv(
  "data/processed/result_all_published_PD.csv",
  header = TRUE
)

data_count_scaled <- readRDS(
  "data/processed/data_count_scaled_published.rds"
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
range(network_size$n_plants_total, na.rm = TRUE)

# ==========================================================
# Prepare plotting data
# ==========================================================

plot_data <- result_all %>%
  select(
    Study_Network_id,
    percentage_Abun10,
    percentage_Abun5,
    percentage_Abun3,
    percentage_FlwShape5,
    percentage_FlwShape3,
    percentage_Pylo10,
    percentage_Pylo5,
    percentage_Pylo3,
    random_mean_percentage_Random10,
    random_mean_percentage_Random5,
    random_mean_percentage_Random3
  ) %>%
  pivot_longer(
    cols = -Study_Network_id,
    names_to = "Option",
    values_to = "Percentage"
  ) %>%
  filter(!is.na(Percentage)) %>%
  mutate(
    Strategy = case_when(
      grepl("Abun", Option) ~ "Flower abundance",
      grepl("FlwShape", Option) ~ "Flower abundance + shapes",
      grepl("Pylo", Option) ~ "Phylogenetic distance",
      grepl("Random", Option) ~ "Random"
    ),
    
    Plant_Number = case_when(
      grepl("10$", Option) ~ "Top10",
      grepl("5$", Option) ~ "Top5",
      grepl("3$", Option) ~ "Top3"
    )
  ) %>%
  left_join(
    network_size,
    by = "Study_Network_id"
  )

# ==========================================================
# Check join
# ==========================================================

plot_data %>%
  filter(is.na(n_plants_total)) %>%
  distinct(Study_Network_id)

# ==========================================================
# Set factor order
# ==========================================================

plot_data$Strategy <- factor(
  plot_data$Strategy,
  levels = c(
    "Random",
    "Phylogenetic distance",
    "Flower abundance + shapes",
    "Flower abundance"
  )
)

plot_data$Plant_Number <- factor(
  plot_data$Plant_Number,
  levels = c(
    "Top10",
    "Top5",
    "Top3"
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
    mean = mean(Percentage, na.rm = TRUE),
    median = median(Percentage, na.rm = TRUE),
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
# Strategy colours for violin fills
# ==========================================================

strategy_colors <- c(
  "Random" = "#B8B8B8",
  "Phylogenetic distance" = "#E69F00",
  "Flower abundance + shapes" = "#8E6BBE",
  "Flower abundance" = "#2AA889"
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
  
  # --------------------------------------------------------
# Violin distributions
# --------------------------------------------------------

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
  
  # --------------------------------------------------------
# Connect the same network across strategies
# --------------------------------------------------------

geom_line(
  aes(
    group = Study_Network_id
  ),
  color = "grey60",
  alpha = 0.08,
  linewidth = 0.2
) +
  
  # --------------------------------------------------------
# Raw observations
# Colour = total flowering-plant richness
# --------------------------------------------------------

geom_point(
  aes(
    color = n_plants_total
  ),
  alpha = 0.50,
  size = 1.1
) +
  
  # --------------------------------------------------------
# Median
# --------------------------------------------------------

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
  
  # --------------------------------------------------------
# Mean
# --------------------------------------------------------

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
  
  # --------------------------------------------------------
# Sampling-effort panels
# --------------------------------------------------------

facet_wrap(
  ~ Plant_Number,
  nrow = 1,
  scales = "free_x"
) +
  
  # ========================================================
# Scales
# ========================================================

scale_fill_manual(
  values = strategy_colors,
  name = "Subsampling strategy"
) +
  
  scale_color_viridis_c(
    option = "turbo",
    begin = 0.18,
    end = 0.90,
    name = "Flowering species per network",
    breaks = breaks_pretty(n = 5)
  ) +
  
  scale_x_discrete(
    labels = c(
      "Random" = "Random",
      "Phylogenetic distance" = "Phylogenetic",
      "Flower abundance + shapes" = "Abundance + shape",
      "Flower abundance" = "Abundance"
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
      size = 12
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
      10, 10, 10, 10
    )
  ) +
  
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
  )

plot_violin_size