###############################################################################
# 3.4 Shannon diversity
###############################################################################

# ==========================================================
# Load data for plotting
# ==========================================================

result_all <- read.csv(
  "result_all_published.csv",
  header = TRUE
)

plotshannon <- result_all %>%
  select(
    Study_Network_id,
    shannon_index_abun5,
    shannon_index_abun3,
    shannon_index_flw5,
    shannon_index_flw3
  )

# If only networks in Op1better should be retained, use:
#
# plotshannon <- result_all %>%
#   select(
#     Study_Network_id,
#     shannon_index_abun5,
#     shannon_index_abun3,
#     shannon_index_flw5,
#     shannon_index_flw3
#   ) %>%
#   filter(
#     Study_Network_id %in% Op1better$Study_Network_id
#   )

# ==========================================================
# Reshape data
# ==========================================================

shannon_long <- plotshannon %>%
  pivot_longer(
    cols = -Study_Network_id,
    names_to = "Option",
    values_to = "Shannon_index"
  )

shannon_long_clean <- shannon_long %>%
  drop_na()

# ==========================================================
# Summarise Shannon diversity
# ==========================================================

data_summary <- shannon_long_clean %>%
  group_by(
    Option
  ) %>%
  summarise(
    Mean_Shannon_index = mean(
      Shannon_index,
      na.rm = TRUE
    ),
    SD_Shannon_index = sd(
      Shannon_index,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

custom_order <- c(
  "shannon_index_abun5",
  "shannon_index_flw5",
  "shannon_index_abun3",
  "shannon_index_flw3"
)

shannon_long_clean$Option <- factor(
  shannon_long_clean$Option,
  levels = custom_order
)

non_missing_count <- shannon_long_clean %>%
  group_by(
    Option
  ) %>%
  summarise(
    n_non_missing = sum(
      !is.na(Shannon_index)
    ),
    .groups = "drop"
  )

data_summary <- data_summary %>%
  left_join(
    non_missing_count,
    by = "Option"
  ) %>%
  mutate(
    `Focus Species` = case_when(
      Option %in% c(
        "shannon_index_abun5",
        "shannon_index_abun3"
      ) ~ "Flower abundance",
      
      Option %in% c(
        "shannon_index_flw5",
        "shannon_index_flw3"
      ) ~ "Flower shapes"
    ),
    
    `Plant Number` = case_when(
      Option %in% c(
        "shannon_index_abun5",
        "shannon_index_flw5"
      ) ~ "Top5",
      
      TRUE ~ "Top3"
    )
  )

data_summary$`Plant Number` <- factor(
  data_summary$`Plant Number`,
  levels = c(
    "Top5",
    "Top3"
  )
)

option_info <- data_summary %>%
  select(
    Option,
    `Focus Species`,
    `Plant Number`
  )

shannon_long_clean <- shannon_long_clean %>%
  left_join(
    option_info,
    by = "Option"
  )

shannon_long_clean$`Plant Number` <- factor(
  shannon_long_clean$`Plant Number`,
  levels = c(
    "Top5",
    "Top3"
  )
)

# ==========================================================
# Colours
# ==========================================================

focus_species_colors <- c(
  "Flower abundance" = "#1B9E77",
  "Flower shapes" = "#9E9AC8"
)

data_summary$`Focus Species` <- factor(
  data_summary$`Focus Species`,
  levels = c(
    "Flower abundance",
    "Flower shapes"
  )
)

shannon_long_clean$`Focus Species` <- factor(
  shannon_long_clean$`Focus Species`,
  levels = c(
    "Flower abundance",
    "Flower shapes"
  )
)

# ==========================================================
# Plot Shannon diversity
# ==========================================================

set.seed(2025)

plot_shannon <- ggplot(
  data_summary,
  aes(
    x = Option,
    y = Mean_Shannon_index,
    color = `Focus Species`
  )
) +
  
  geom_jitter(
    data = shannon_long_clean,
    aes(
      x = Option,
      y = Shannon_index,
      color = `Focus Species`
    ),
    width = 0.2,
    size = 1,
    alpha = 0.5
  ) +
  
  geom_crossbar(
    aes(
      ymin = Mean_Shannon_index - SD_Shannon_index,
      ymax = Mean_Shannon_index + SD_Shannon_index
    ),
    position = position_dodge(
      width = 0.3
    ),
    width = 0.5,
    color = "grey40"
  ) +
  
  geom_text(
    aes(
      label = paste(
        "n =",
        n_non_missing
      )
    ),
    position = position_dodge(
      width = 0.5
    ),
    vjust = -1,
    size = 4,
    color = "black"
  ) +
  
  geom_text(
    aes(
      label = round(
        Mean_Shannon_index,
        2
      )
    ),
    position = position_dodge(
      width = 0.5
    ),
    vjust = 1.5,
    size = 4,
    color = "black"
  ) +
  
  facet_wrap(
    ~ `Plant Number`,
    nrow = 1,
    scales = "free_x"
  ) +
  
  scale_color_manual(
    values = focus_species_colors
  ) +
  
  labs(
    x = "Subsampling option",
    y = "Shannon index",
    color = "Plants selected based on"
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "grey80",
      fill = NA,
      linewidth = 0.5
    ),
    plot.margin = margin(
      10,
      10,
      10,
      10
    ),
    strip.background = element_rect(
      fill = "white",
      color = NA
    ),
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 10
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 5
      )
    )
  )

plot_shannon

# ==========================================================
# Save Shannon plot
# ==========================================================

ggsave(
  "result251105_published/shannon.png",
  plot_shannon,
  width = 7.5,
  height = 4.5,
  units = "in",
  dpi = 300
)

