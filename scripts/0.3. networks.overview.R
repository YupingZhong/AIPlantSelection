############################################################
# Optimising plant species selection for automated monitoring
# Networks overview
############################################################

library(dplyr)
library(ggplot2)
library(cowplot)

# Load data
data_count_scaled <- readRDS("data/processed/data_count_scaled_published.rds")
data_interact <- readRDS("data/processed/data_interact_published.rds")

# Prepare summary data
plant_diversity <- data_count_scaled %>%
  filter(Flower_count != 0) %>%
  group_by(Study_Network_id) %>%
  summarise(plant_sp_number = n_distinct(Plant_species, na.rm = TRUE))

pollinator_diversity <- data_interact %>%
  group_by(Study_Network_id) %>%
  summarise(pollinator_sp_number = n_distinct(Pollinator_accepted_name, na.rm = TRUE))

shape_freq <- data_count_scaled %>%
  filter(Flower_count_scaled != 0, !is.na(Plant_species)) %>%
  group_by(Study_Network_id) %>%
  summarise(shape_number = n_distinct(flw_shape_revised, na.rm = TRUE)) %>%
  count(shape_number, name = "freq") %>%
  mutate(shape_number = factor(shape_number))

shape_class_freq <- data_count_scaled %>%
  filter(
    Flower_count_scaled != 0,
    !is.na(Plant_species),
    !is.na(flw_shape_revised)
  ) %>%
  distinct(Study_Network_id, flw_shape_revised) %>%
  count(flw_shape_revised, name = "network_number") %>%
  arrange(desc(network_number))

shape_order <- shape_class_freq$flw_shape_revised
x_max <- max(
  plant_diversity$plant_sp_number,
  pollinator_diversity$pollinator_sp_number,
  na.rm = TRUE
)
breaks_all <- seq(0, x_max, by = 10)

# Shared theme for the first three plots
distribution_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(size = 12, face = "bold")
  )

# Individual plots
histogram_1 <- ggplot(plant_diversity, aes(plant_sp_number, fill = "#66C2A5")) +
  geom_histogram(binwidth = 2, colour = "black", boundary = 0) +
  labs(x = "Flower species", y = "Number of Networks") +
  scale_fill_identity() +
  scale_x_continuous(
    breaks = breaks_all,
    limits = c(0, x_max),
    expand = c(0, 0)
  ) +
  distribution_theme

histogram_2 <- ggplot(pollinator_diversity, aes(pollinator_sp_number, fill = "#E9C46A")) +
  geom_histogram(binwidth = 2, colour = "black", boundary = 0) +
  labs(x = "Pollinator species in interaction data", y = "Number of Networks") +
  scale_fill_identity() +
  scale_x_continuous(
    breaks = breaks_all,
    limits = c(0, x_max),
    expand = c(0, 0)
  ) +
  distribution_theme

shape_barplot <- ggplot(shape_freq, aes(shape_number, freq, fill = "#A5B5D9")) +
  geom_col(colour = "black") +
  labs(x = "Number of flower shapes", y = "Number of Networks") +
  scale_fill_identity() +
  distribution_theme

shape_class <- ggplot(
  shape_class_freq,
  aes(factor(flw_shape_revised, levels = shape_order), network_number)
) +
  geom_col(
    width = 0.8,
    fill = "#B07AA1",
    colour = "black",
    linewidth = 0.35
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Flower shape categories", y = "Number of networks") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "grey55", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.text.x = element_text(colour = "black", size = 9, angle = 45, hjust = 1),
    plot.margin = margin(8, 8, 5, 5)
  )

# Shared supplementary-figure styling
supp_theme <- theme(
  panel.border = element_rect(colour = "grey55", fill = NA, linewidth = 0.6),
  axis.line = element_blank(),
  axis.title = element_text(face = "bold", size = 12),
  axis.text = element_text(colour = "black", size = 10),
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

histogram_1_supp <- histogram_1 +
  labs(tag = "(a)", y = "Number of networks") +
  supp_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

histogram_2_supp <- histogram_2 +
  labs(tag = "(b)", y = "Number of networks") +
  supp_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

shape_barplot_supp <- shape_barplot +
  labs(tag = "(c)", y = "Number of networks") +
  supp_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

shape_class_supp <- shape_class +
  labs(tag = "(d)") +
  supp_theme +
  theme(axis.text.y = element_text(colour = "black", size = 9))

# Arrange and export the figure
top_row <- plot_grid(
  histogram_1_supp,
  histogram_2_supp,
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1)
)

bottom_row <- plot_grid(
  shape_barplot_supp,
  shape_class_supp,
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1)
)

supp_fig <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1.1)
)

supp_fig

ggsave(
  "/Chap1_TargetPlant_to_monitor/result_260723/Supplementary_Figure.png",
  supp_fig,
  width = 10,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "/Chap1_TargetPlant_to_monitor/result_260723/Supplementary_Figure.pdf",
  supp_fig,
  width = 12,
  height = 9,
  units = "in",
  bg = "white"
)

