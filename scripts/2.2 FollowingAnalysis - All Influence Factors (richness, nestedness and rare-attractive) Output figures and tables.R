################################################################################
# Factors associated with abundance-based subsampling effectiveness
#
# This script:
#   1. Loads the processed analysis results
#   2. Calculates richness and NODF associations
#   3. Creates the main three-panel figure
#   4. Creates the associated result tables
################################################################################


# ==============================================================================
# 1. Libraries
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(vegan)
library(bipartite)
library(stringr)
library(lme4)
library(flextable)
library(officer)


# ==============================================================================
# 2. Paths
# ==============================================================================

data_dir <- file.path("data", "processed")
result_dir <- "/Chap1_TargetPlant_to_monitor/result_submission"

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. Load original data required for plotting
# ==============================================================================

data_count_scaled <- readRDS(
  file.path(data_dir, "data_count_scaled_published.rds")
)

data_interact <- readRDS(
  file.path(data_dir, "data_interact_published.rds")
)

result_all <- read.csv(
  file.path(data_dir, "result_all_published_PD.csv")
)


# ==============================================================================
# 4. Load H2 analysis results
# ==============================================================================

Rare_favour_plants <- readRDS(
  file.path(data_dir, "H2_Rare_favour_plants.rds")
)

network_group <- readRDS(
  file.path(data_dir, "H2_network_group.rds")
)

attract_percent_10 <- readRDS(
  file.path(data_dir, "H2_attract_percent_10.rds")
)

wilcox_res <- readRDS(
  file.path(data_dir, "H2_wilcox_res.rds")
)


# ==============================================================================
# 5. Load trait data for Rare–Attractive species table
# ==============================================================================

traits <- read.csv("data/processed/merge.trait.csv")

# ==============================================================================
# 6. Summary statistics for Panel C
# ==============================================================================

data_summary <- attract_percent_10 %>%
  group_by(present_group) %>%
  summarise(
    n = n(),
    Mean_Percentage = mean(
      percentage_Abun10,
      na.rm = TRUE
    ),
    SD_Percentage = sd(
      percentage_Abun10,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ==============================================================================
# 7. Common settings
# ==============================================================================

common_ylim <- c(0, 100)
common_breaks <- seq(0, 100, 20)

richness_colors <- c(
  "Plant richness" = "#69B89C",
  "Pollinator richness" = "#E6B84A"
)

presence_colors <- c(
  "Present" = "#8064A2",
  "Absent" = "#AFA9B8"
)

presence_point_colors <- c(
  "Present" = "#604477",
  "Absent" = "#77717D"
)

nestedness_color <- "#5B8FD9"


# ==============================================================================
# 8. Plant richness
# ==============================================================================

plant_rich <- data_count_scaled %>%
  filter(Flower_count_scaled != 0) %>%
  group_by(Study_Network_id) %>%
  summarise(
    Plant_Richness = n_distinct(
      Plant_species,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ==============================================================================
# 9. Pollinator richness
# ==============================================================================

poll_rich <- result_all %>%
  dplyr::select(
    Study_Network_id,
    total_pollinator_count_Abun10
  )


# ==============================================================================
# 10. Richness comparison data
# ==============================================================================

richness_compare <- bind_rows(
  
  result_all %>%
    dplyr::select(
      Study_Network_id,
      percentage_Abun10
    ) %>%
    left_join(
      plant_rich,
      by = "Study_Network_id"
    ) %>%
    filter(
      !is.na(percentage_Abun10),
      !is.na(Plant_Richness)
    ) %>%
    transmute(
      Study_Network_id,
      percentage_Abun10,
      Predictor = "Plant richness",
      Richness = Plant_Richness
    ),
  
  result_all %>%
    dplyr::select(
      Study_Network_id,
      percentage_Abun10
    ) %>%
    left_join(
      poll_rich,
      by = "Study_Network_id"
    ) %>%
    filter(
      !is.na(percentage_Abun10),
      !is.na(total_pollinator_count_Abun10)
    ) %>%
    transmute(
      Study_Network_id,
      percentage_Abun10,
      Predictor = "Pollinator richness",
      Richness = total_pollinator_count_Abun10
    )
)


# ==============================================================================
# 11. Richness statistics
# ==============================================================================

# Correlations and linear models are fitted using richness on its original
# scale. Panel (a) subsequently applies a log10 coordinate transformation for
# display. Consequently, the fitted lines may appear curved on the plotted
# log10 x-axis even though the fitted models are linear on the original scale.

stats_richness <- richness_compare %>%
  group_by(Predictor) %>%
  summarise(
    r = cor(
      Richness,
      percentage_Abun10,
      method = "pearson"
    ),
    p = summary(
      lm(
        percentage_Abun10 ~ Richness
      )
    )$coefficients[2, 4],
    .groups = "drop"
  )


# ==============================================================================
# 12. Network nestedness
# ==============================================================================

nestedness_df <- data_interact %>%
  filter(
    !is.na(Interaction_addup),
    Interaction_addup > 0
  ) %>%
  group_by(Study_Network_id) %>%
  summarise(
    
    Nestedness = {
      
      mat <- table(
        Pollinator_accepted_name,
        Plant_accepted_name
      )
      
      mat_bin <- ifelse(mat > 0, 1, 0)
      
      if (
        nrow(mat_bin) > 1 &&
        ncol(mat_bin) > 1
      ) {
        bipartite::networklevel(
          mat_bin,
          index = "NODF"
        )
      } else {
        NA_real_
      }
      
    },
    
    .groups = "drop"
  )


df_nested <- result_all %>%
  dplyr::select(
    Study_Network_id,
    percentage_Abun10
  ) %>%
  left_join(
    nestedness_df,
    by = "Study_Network_id"
  ) %>%
  filter(
    !is.na(percentage_Abun10),
    !is.na(Nestedness)
  )


cor_nested <- cor.test(
  df_nested$Nestedness,
  df_nested$percentage_Abun10,
  method = "pearson"
)

lm_nested <- lm(
  percentage_Abun10 ~ Nestedness,
  data = df_nested
)

p_nested <- summary(
  lm_nested
)$coefficients[2, 4]

format_p <- function(p) {
  ifelse(
    p < 0.001,
    "p < 0.001",
    paste0("p = ", format.pval(p, digits = 2))
  )
}
# ==============================================================================
# 13. PANEL A
# ==============================================================================
plant_stat_label <- paste0(
  "Plant: R = ",
  round(
    stats_richness$r[
      stats_richness$Predictor == "Plant richness"
    ],
    2
  ),
  ", ",
  format_p(
    stats_richness$p[
      stats_richness$Predictor == "Plant richness"
    ]
  )
)

pollinator_stat_label <- paste0(
  "Pollinator: R = ",
  round(
    stats_richness$r[
      stats_richness$Predictor == "Pollinator richness"
    ],
    2
  ),
  ", ",
  format_p(
    stats_richness$p[
      stats_richness$Predictor == "Pollinator richness"
    ]
  )
)

p_richness <- ggplot(
  richness_compare,
  aes(
    x = Richness,
    y = percentage_Abun10,
    color = Predictor,
    shape = Predictor
  )
) +
  geom_point(
    aes(size = Predictor),
    alpha = 0.75,
    position = position_jitter(
      width = 0.03,
      height = 0.4
    )
  ) +
  
  scale_size_manual(
    values = c(
      "Plant richness" = 1.6,
      "Pollinator richness" = 1.2
    ),
    guide = "none"
  )  +
  
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.9
  )  +
  
  scale_color_manual(
    values = richness_colors,
    breaks = c(
      "Plant richness",
      "Pollinator richness"
    ),
    labels = c(
      "Plant",
      "Pollinator"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Plant richness" = 18,
      "Pollinator richness" = 15
    ),
    breaks = c(
      "Plant richness",
      "Pollinator richness"
    ),
    labels = c(
      "Plant",
      "Pollinator"
    )
  ) +
  scale_x_continuous(
    breaks = scales::breaks_log(n = 5),
    labels = scales::label_number()
  )  +
  
  scale_y_continuous(
    limits = common_ylim,
    breaks = common_breaks
  ) +
  
  labs(
    x = expression("Species richness (" * log[10] * " scale)"),
    y = "Pollinator richness captured (%)",
    color = NULL,
    shape = NULL
  ) +
  annotate(
    "text",
    x = 10.5,
    y = 13,
    label = plant_stat_label,
    color = richness_colors["Plant richness"],
    hjust = 0,
    vjust = 0,
    size = 5.5
  ) +
  
  annotate(
    "text",
    x = 10.5,
    y = 6,
    label = pollinator_stat_label,
    color = richness_colors["Pollinator richness"],
    hjust = 0,
    vjust = 0,
    size = 5.5
  )  +
  
  # Transform the displayed x-axis after fitting the linear models above.
  # This retains Pearson correlations and lm fits on raw richness while
  # presenting plant and pollinator richness on a log10 scale.
  coord_trans(
    x = "log10"
  ) +
  
  theme_classic(base_size = 17) +
  theme(
    axis.title = element_text(
      face = "bold",
      size = 17
    ),
    axis.text.y = element_text(
      color = "black",
      size = 13
    ),
    axis.ticks.y = element_line(
      color = "black"
    ),
    legend.position = "none"
  )


# ==============================================================================
# 14. PANEL B
# ==============================================================================
nested_annotation_x <- min(
  df_nested$Nestedness,
  na.rm = TRUE
) + 0.05 * diff(
  range(df_nested$Nestedness, na.rm = TRUE)
)

p_nestedness <- ggplot(
  df_nested,
  aes(
    x = Nestedness,
    y = percentage_Abun10
  )
) +
  
  geom_point(
    size = 1.4,
    alpha = 0.70,
    color = nestedness_color,
    position = position_jitter(
      width = 0.01,
      height = 0.4
    )
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.9,
    color = nestedness_color,
    linetype = ifelse(
      p_nested < 0.05,
      "solid",
      "dashed"
    )
  ) +
  
  scale_y_continuous(
    limits = common_ylim,
    breaks = common_breaks
  ) +
  annotate(
    "text",
    x = nested_annotation_x,
    y = 5,
    hjust = 0,
    vjust = 0,
    label = paste0(
      "R = ",
      round(cor_nested$estimate, 2),
      "\n",
      format_p(cor_nested$p.value)
    ),
    size = 5.5
  )  +
  
  labs(
    x = "Network nestedness",
    y = NULL
  ) +
  
  theme_classic(base_size = 17) +
  theme(
    axis.title.x = element_text(
      face = "bold",
      size = 17
    ),
    axis.text.x = element_text(
      color = "black",
      size = 13
    ),
    axis.text.y = element_text(
      color = "black",
      size = 13
    ),
    axis.ticks.y = element_line(
      color = "black"
    )
  )


# ==============================================================================
# 15. PANEL C
# ==============================================================================

attract_percent_10 <- attract_percent_10 %>%
  mutate(
    present_group = factor(
      str_to_title(
        str_squish(present_group)
      ),
      levels = c(
        "Present",
        "Absent"
      )
    )
  )


sig_label <- case_when(
  wilcox_res$p.value < 0.001 ~ "***",
  wilcox_res$p.value < 0.01 ~ "**",
  wilcox_res$p.value < 0.05 ~ "*",
  TRUE ~ "ns"
)


p_presence <- ggplot(
  attract_percent_10,
  aes(
    x = present_group,
    y = percentage_Abun10,
    fill = present_group
  )
) +
  
  geom_violin(
    trim = FALSE,
    alpha = 0.25,
    color = NA
  ) +
  
  geom_boxplot(
    width = 0.16,
    alpha = 0.65,
    outlier.shape = NA,
    color = "grey75"
  ) +
  
  geom_jitter(
    aes(color = present_group),
    width = 0.10,
    size = 1.5,
    alpha = 0.55
  ) +
  
  geom_text(
    data = data_summary,
    aes(
      x = present_group,
      y = 8,
      label = paste0("n = ", n)
    ),
    inherit.aes = FALSE,
    size = 5.5,
    vjust = 0
  ) +
  
  geom_text(
    data = data_summary,
    aes(
      x = present_group,
      y = 2,
      label = sprintf(
        "Mean = %.1f%%",
        Mean_Percentage
      )
    ),
    inherit.aes = FALSE,
    size = 5.5,
    vjust = 0
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = 88,
    yend = 88
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = 88,
    yend = 86
  ) +
  
  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = 88,
    yend = 86
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = 90,
    label = sig_label,
    size = 7,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = presence_colors
  ) +
  
  scale_color_manual(
    values = presence_point_colors
  ) +
  
  scale_y_continuous(
    limits = common_ylim,
    breaks = common_breaks
  ) +
  
  labs(
    x = "Rare attractive plants",
    y = NULL
  ) +
  
  theme_classic(base_size = 17) +
  theme(
    axis.title.x = element_text(
      face = "bold",
      size = 17
    ),
    axis.text.x = element_text(
      color = "black",
      size = 13
    ),
    axis.text.y = element_text(
      color = "black",
      size = 13
    ),
    axis.ticks.y = element_line(
      color = "black"
    ),
    legend.position = "none"
  )

common_panel_theme <- theme(
  # Make the plotting regions approximately square
  
  
  # Border around all four sides
  panel.border = element_rect(
    colour = "grey55",
    fill = NA,
    linewidth = 2.2
  ),
  
  # Prevent theme_classic axis lines from overlapping the border
  axis.line = element_blank(),
  
  # Retain ticks
  axis.ticks.x = element_line(
    colour = "black",
    linewidth = 0.6
  ),
  axis.ticks.y = element_line(
    colour = "black",
    linewidth = 0.6
  ),
  
  axis.text.x = element_text(
    colour = "black",
    size = 15
  ),
  axis.text.y = element_text(
    colour = "black",
    size = 15
  ),
  
  panel.grid = element_blank(),
  
  plot.margin = margin(
    t = 24,
    r = 8,
    b = 8,
    l = 8
  )
)

p_richness <- p_richness + common_panel_theme
p_nestedness <- p_nestedness + common_panel_theme
p_presence <- p_presence + common_panel_theme
# ==============================================================================
# 16. Main three-panel figure
# ==============================================================================

# Extract a horizontal legend
shared_legend <- cowplot::get_legend(
  p_richness +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "left",
      legend.box.just = "left",
      legend.margin = margin(2, 0, 2, 0),
      legend.text = element_text(size = 15),
      legend.key.size = unit(1.1, "lines")
    ) +
    guides(
      color = guide_legend(nrow = 1),
      shape = guide_legend(nrow = 1)
    )
)

# Three square panels without individual legends

panels <- plot_grid(
  p_richness + theme(legend.position = "none"),
  p_nestedness,
  p_presence,
  ncol = 3,
  rel_widths = c(1, 1, 1),
  labels = c("(a)", "(b)", "(c)"),
  label_size = 25,
  label_fontface = "bold",
  label_x = 0.01,
  label_y = 0.997,
  hjust = 0,
  vjust = 0.95,
  align = "hv",
  axis = "tblr"
)

# Left-aligned legend
legend_left <- plot_grid(
  NULL,
  shared_legend,
  NULL,
  ncol = 3,
  rel_widths = c(0.07, 0.48, 0.45)
)

# Final figure
final_factor_fig <- plot_grid(
  panels,
  legend_left,
  ncol = 1,
  rel_heights = c(8.9, 0.43)
)

print(final_factor_fig)

ggsave(
  filename = file.path(
    result_dir,
    "Figure_3.tiff"
  ),
  plot = final_factor_fig,
  device = "tiff",
  width = 14,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ==============================================================================
# 17. Rare–Attractive plant species table
# ==============================================================================

rarelist <- Rare_favour_plants %>%
  distinct(Plant_accepted_name)


sp_list_rare <- rarelist %>%
  left_join(
    traits,
    by = "Plant_accepted_name"
  ) %>%
  distinct(
    Plant_accepted_name,
    flw_shape_revised
  )


sp_list_grouped <- sp_list_rare %>%
  filter(
    !is.na(flw_shape_revised)
  ) %>%
  group_by(
    flw_shape_revised
  ) %>%
  summarise(
    Plant_species = paste(
      sort(unique(Plant_accepted_name)),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  arrange(
    flw_shape_revised
  ) %>%
  mutate(
    Illustration = ""
  ) %>%
  dplyr::select(
    Illustration,
    flw_shape_revised,
    Plant_species
  )


sp_list_grouped$Plant_species <- stringr::str_wrap(
  sp_list_grouped$Plant_species,
  width = 60
)


ft_rare <- flextable(
  sp_list_grouped
) %>%
  set_header_labels(
    Illustration = "",
    flw_shape_revised = "Flower shape",
    Plant_species = "Plant species"
  ) %>%
  theme_booktabs() %>%
  fontsize(
    size = 11,
    part = "all"
  ) %>%
  italic(
    j = "Plant_species",
    part = "body"
  ) %>%
  width(
    j = "Illustration",
    width = 1
  ) %>%
  width(
    j = "flw_shape_revised",
    width = 2.5
  ) %>%
  width(
    j = "Plant_species",
    width = 4.5
  ) %>%
  align(
    j = c(
      "Illustration",
      "flw_shape_revised",
      "Plant_species"
    ),
    align = "left"
  )


doc_rare <- read_docx() %>%
  body_add_par(
    "Rare but highly attractive plant species requiring additional sampling. ",
    style = "heading 2"
  ) %>%
  body_add_flextable(
    ft_rare
  )


print(
  doc_rare,
  target = file.path(
    result_dir,
    "sp_list_rare.docx"
  )
)


# ==============================================================================
# 18. Wilcoxon comparison table
# ==============================================================================

table2 <- attract_percent_10 %>%
  group_by(
    present_group
  ) %>%
  summarise(
    Mean = mean(
      percentage_Abun10,
      na.rm = TRUE
    ),
    SD = sd(
      percentage_Abun10,
      na.rm = TRUE
    ),
    n = sum(
      !is.na(percentage_Abun10)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    `Richness captured (%)` = sprintf(
      "%.2f ± %.2f",
      Mean,
      SD
    ),
    `P-value` = signif(
      wilcox_res$p.value,
      3
    )
  ) %>%
  dplyr::select(
    Group = present_group,
    `Richness captured (%)`,
    n,
    `P-value`
  )


ft_table2 <- flextable(
  table2
) %>%
  theme_booktabs() %>%
  autofit() %>%
  fontsize(
    size = 10,
    part = "all"
  ) %>%
  align(
    align = "center",
    part = "all"
  ) %>%
  bold(
    part = "header"
  ) %>%
  align(
    j = "Group",
    align = "left",
    part = "all"
  )


doc_table2 <- read_docx() %>%
  body_add_par(
    paste0(
      "Table 2. Comparison of pollinator richness captured by the ",
      "top 10 abundant plants between networks with and without ",
      "highly attractive plant species (Wilcoxon rank-sum test)."
    ),
    style = "heading 2"
  ) %>%
  body_add_flextable(
    ft_table2
  )


print(
  doc_table2,
  target = file.path(
    result_dir,
    "Table2_Wilcoxon.docx"
  )
)


# ==============================================================================
# 19. Export CSV summaries
# ==============================================================================

write.csv(
  table2,
  file.path(
    result_dir,
    "H2_Wilcoxon_group_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  rarelist,
  file.path(
    result_dir,
    "H2_rare_attractive_plants.csv"
  ),
  row.names = FALSE
)

