################################################################################
# RESULT 2
# Figures and tables for H2
#
# This script ONLY:
#   1. Loads H2 analysis results
#   2. Creates the main three-panel figure
#   3. Creates the main H2 result tables
#
# All statistical analyses are performed in Script 1.
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

data_dir <- "data/processed"
result_dir <- "/Chap1_TargetPlant_to_monitor/result_260723"

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
  filter(!is.na(Interaction_addup)) %>%
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
          index = "nestedness"
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


# ==============================================================================
# 13. PANEL A
# ==============================================================================

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
    size = 2.4,
    alpha = 0.75,
    position = position_jitter(
      width = 0.03,
      height = 0.4
    )
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.9
  ) +
  
  scale_color_manual(
    values = richness_colors
  ) +
  
  scale_shape_manual(
    values = c(
      "Plant richness" = 18,
      "Pollinator richness" = 15
    )
  ) +
  
  scale_y_continuous(
    limits = common_ylim,
    breaks = common_breaks
  ) +
  
  labs(
    x = "Richness",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = 97,
    hjust = 1.05,
    vjust = 1,
    label = paste0(
      "Plant: R = ",
      round(
        stats_richness$r[
          stats_richness$Predictor == "Plant richness"
        ],
        2
      ),
      ", p = ",
      format.pval(
        stats_richness$p[
          stats_richness$Predictor == "Plant richness"
        ],
        digits = 2
      ),
      "\nPollinator: R = ",
      round(
        stats_richness$r[
          stats_richness$Predictor == "Pollinator richness"
        ],
        2
      ),
      ", p = ",
      format.pval(
        stats_richness$p[
          stats_richness$Predictor == "Pollinator richness"
        ],
        digits = 2
      )
    ),
    size = 3.1
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.y = element_text(
      color = "black",
      size = 11
    ),
    axis.ticks.y = element_line(
      color = "black"
    ),
    legend.position = "top"
  )


# ==============================================================================
# 14. PANEL B
# ==============================================================================

p_nestedness <- ggplot(
  df_nested,
  aes(
    x = Nestedness,
    y = percentage_Abun10
  )
) +
  
  geom_point(
    size = 2.4,
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
    x = Inf,
    y = 97,
    hjust = 1.05,
    vjust = 1,
    label = paste0(
      "R = ",
      round(cor_nested$estimate, 2),
      "\np = ",
      format.pval(
        cor_nested$p.value,
        digits = 2
      )
    ),
    size = 3.1
  ) +
  
  labs(
    x = "Network nestedness",
    y = NULL
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.title.x = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.x = element_text(
      color = "black",
      size = 11
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
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
    color = "grey25"
  ) +
  
  geom_jitter(
    aes(color = present_group),
    width = 0.10,
    size = 1.6,
    alpha = 0.55
  ) +
  
  geom_text(
    data = data_summary,
    aes(
      x = present_group,
      y = 72,
      label = sprintf(
        "Mean = %.1f%%",
        Mean_Percentage
      )
    ),
    inherit.aes = FALSE,
    size = 3.1,
    fontface = "bold"
  ) +
  
  geom_text(
    data = data_summary,
    aes(
      x = present_group,
      y = 78,
      label = paste0("n = ", n)
    ),
    inherit.aes = FALSE,
    size = 3
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
    size = 5,
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
    x = "Rare attractive plants in full network",
    y = NULL
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.title.x = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.x = element_text(
      color = "black",
      size = 11
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none"
  )


# ==============================================================================
# 16. Main three-panel figure
# ==============================================================================

final_factor_fig <- plot_grid(
  p_richness,
  p_nestedness,
  p_presence,
  ncol = 3,
  labels = c("a", "b", "c"),
  label_size = 16,
  label_fontface = "bold",
  align = "h"
)

y_title <- ggdraw() +
  draw_label(
    "Pollinator richness captured (%)",
    angle = 90,
    fontface = "bold",
    size = 13
  )

final_factor_fig <- plot_grid(
  y_title,
  final_factor_fig,
  ncol = 2,
  rel_widths = c(0.06, 1)
)

print(final_factor_fig)

ggsave(
  file.path(
    result_dir,
    "Fig_Factors_subsampling_effectiveness.png"
  ),
  final_factor_fig,
  width = 10.5,
  height = 4,
  units = "in",
  dpi = 600,
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

