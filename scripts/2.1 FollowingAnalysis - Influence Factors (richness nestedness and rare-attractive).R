################################################################################
# ==============================================================================
# RESULT 2 Exploring Op1
# Factors influencing the effectiveness of plant subsampling
#
# Final Figure:
# A. Plant & pollinator richness
# B. Network nestedness
# C. Presence/absence of rare attractive plants
#
# All panels share one common Y-axis:
# Percent of pollinator richness captured (%)
#
# Y-axis is fixed from 0 to 100.
# Only Panel A displays the Y-axis.
# ==============================================================================

################################################################################


# ==============================================================================
# 0. Libraries
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(vegan)
library(bipartite)
library(flextable)
library(officer)
library(stargazer)
library(sandwich)
library(lmtest)
library(MASS)
library(broom)
library(stringr)
library(lme4)
library(lmerTest)
library(purrr)
library(ggsignif)


# ==============================================================================
# 1. Load data
# ==============================================================================

data_count_scaled <- readRDS(
  "data/processed/data_count_scaled_published.rds"
)

data_interact <- readRDS(
  "data/processed/data_interact_published.rds"
)

result_all <- read.csv(
  "data/processed/result_all_published_PD.csv"
)

data_merge <- readRDS(
  "data/processed/data_merge.rds"
)

percent_10 <- readRDS(
  "data/processed/percent_10.rds"
)


# ==============================================================================
# 2. Plant richness
# ==============================================================================

plant_rich <- data_count_scaled %>%
  filter(
    Flower_count_scaled != 0
  ) %>%
  group_by(
    Study_Network_id
  ) %>%
  summarise(
    Plant_Richness = n_distinct(
      Plant_species,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ==============================================================================
# 3. Pollinator richness
# ==============================================================================

poll_rich <- result_all[
  ,
  c(
    "Study_Network_id",
    "total_pollinator_count_Abun10"
  )
]


# ==============================================================================
# 4. Network nestedness
# ==============================================================================

nestedness_df <- data_interact %>%
  filter(
    !is.na(Interaction_addup)
  ) %>%
  group_by(
    Study_Network_id
  ) %>%
  summarise(
    
    Nestedness = {
      
      mat <- table(
        Pollinator_accepted_name,
        Plant_accepted_name
      )
      
      mat_bin <- ifelse(
        mat > 0,
        1,
        0
      )
      
      if (
        nrow(mat_bin) > 1 &
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


# ==============================================================================
# 5. Prepare richness comparison data
# ==============================================================================

plant_df <- result_all[
  ,
  c(
    "Study_Network_id",
    "percentage_Abun10"
  )
] %>%
  
  left_join(
    plant_rich,
    by = "Study_Network_id"
  ) %>%
  
  filter(
    !is.na(percentage_Abun10),
    !is.na(Plant_Richness)
  ) %>%
  
  mutate(
    Predictor = "Plant richness",
    Richness = Plant_Richness
  )


pollinator_df <- result_all[
  ,
  c(
    "Study_Network_id",
    "percentage_Abun10"
  )
] %>%
  
  left_join(
    poll_rich,
    by = "Study_Network_id"
  ) %>%
  
  filter(
    !is.na(percentage_Abun10),
    !is.na(total_pollinator_count_Abun10)
  ) %>%
  
  mutate(
    Predictor = "Pollinator richness",
    Richness = total_pollinator_count_Abun10
  )


richness_compare <- bind_rows(
  
  plant_df[
    ,
    c(
      "Study_Network_id",
      "percentage_Abun10",
      "Predictor",
      "Richness"
    )
  ],
  
  pollinator_df[
    ,
    c(
      "Study_Network_id",
      "percentage_Abun10",
      "Predictor",
      "Richness"
    )
  ]
  
)


# ==============================================================================
# 6. Richness statistics
# ==============================================================================

stats_richness <- richness_compare %>%
  
  group_by(
    Predictor
  ) %>%
  
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

print(stats_richness)


# ==============================================================================
# 7. Nestedness data + statistics
# ==============================================================================

df_nested <- result_all[
  ,
  c(
    "Study_Network_id",
    "percentage_Abun10"
  )
] %>%
  
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


nested_line <- ifelse(
  p_nested < 0.05,
  "solid",
  "dashed"
)


nested_label <- paste0(
  "R = ",
  round(
    cor_nested$estimate,
    2
  ),
  "\n",
  "p = ",
  format.pval(
    cor_nested$p.value,
    digits = 2
  )
)


# ==============================================================================
# 8. Prepare rare attractive plant analysis
# ==============================================================================

traits <- read.csv(
  "data/processed/merge.trait.csv",
  header = TRUE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------------------------
# Top 10 abundant plant species in each network
# ------------------------------------------------------------------------------

top_10_species <- data_count_scaled %>%
  
  distinct(
    Study_Network_id,
    Plant_species,
    Flower_count_scaled,
    .keep_all = FALSE
  ) %>%
  
  group_by(
    Study_Network_id
  ) %>%
  
  slice_max(
    order_by = Flower_count_scaled,
    n = 10,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  group_by(
    Study_Network_id
  ) %>%
  
  filter(
    n_distinct(Plant_species) > 9
  ) %>%
  
  ungroup() %>%
  
  dplyr::select(
    Study_Network_id,
    Plant_species
  )


# ------------------------------------------------------------------------------
# Total interactions in each network
# ------------------------------------------------------------------------------

interaction_all <- data_interact %>%
  
  group_by(
    Study_Network_id
  ) %>%
  
  filter(
    Interaction_addup != 0
  ) %>%
  
  summarise(
    
    Interaction_sum =
      sum(
        Interaction_addup
      ),
    
    richness_sum =
      n_distinct(
        Pollinator_accepted_name
      ),
    
    .groups = "drop"
    
  )


# ------------------------------------------------------------------------------
# Plant name matching
# ------------------------------------------------------------------------------

match <- data_interact %>%
  
  ungroup() %>%
  
  dplyr::select(
    Plant_accepted_name,
    Plant_original_name
  ) %>%
  
  distinct() %>%
  
  rename(
    Plant_species = Plant_original_name
  )


# ------------------------------------------------------------------------------
# Combined plant-level data
# ------------------------------------------------------------------------------

data_combined_all <- data_merge %>%
  
  ungroup() %>%
  
  left_join(
    interaction_all,
    by = "Study_Network_id"
  ) %>%
  
  left_join(
    traits,
    by = "Plant_accepted_name"
  ) %>%
  
  group_by(
    Study_id,
    Study_Network_id,
    Plant_species
  ) %>%
  
  summarise(
    
    Abundance =
      first(
        Flower_count_scaled
      ),
    
    Pollinator_richness =
      n_distinct(
        Pollinator_accepted_name,
        na.rm = TRUE
      ),
    
    Interaction_times =
      sum(
        Interaction_addup,
        na.rm = TRUE
      ),
    
    Proportion_of_interactions =
      (
        Interaction_times /
          first(Interaction_sum)
      ) * 100,
    
    Proportion_of_richness =
      (
        Pollinator_richness /
          first(richness_sum)
      ) * 100,
    
    .groups = "drop"
    
  ) %>%
  
  filter(
    is.finite(Abundance),
    is.finite(Interaction_times)
  ) %>%
  
  group_by(
    Study_Network_id
  ) %>%
  
  mutate(
    
    Abundance_scaled =
      scale(
        log(
          Abundance + 1
        )
      )[, 1]
    
  ) %>%
  
  ungroup() %>%
  
  left_join(
    match,
    by = "Plant_species"
  )


# ==============================================================================
# 9. Mixed model for plant abundance
# ==============================================================================

model <- lmer(
  
  Proportion_of_richness ~
    Abundance_scaled +
    (1 | Study_id / Study_Network_id),
  
  data = data_combined_all
  
)

summary(model)


# ------------------------------------------------------------------------------
# Marginal and conditional R2
# ------------------------------------------------------------------------------

fixed_pred <- predict(
  model,
  re.form = NA
)

var_fixed <- var(
  fixed_pred,
  na.rm = TRUE
)

var_rand <- sum(
  as.data.frame(
    VarCorr(model)
  )$vcov
)

var_resid <- attr(
  VarCorr(model),
  "sc"
)^2


r2_marginal <-
  var_fixed /
  (
    var_fixed +
      var_rand +
      var_resid
  )


r2_conditional <-
  (
    var_fixed +
      var_rand
  ) /
  (
    var_fixed +
      var_rand +
      var_resid
  )


# ==============================================================================
# 10. Predicted capture + residual
# ==============================================================================

data_combined_all$Predicted_capture <-
  predict(
    model,
    newdata = data_combined_all,
    re.form = NA,
    allow.new.levels = TRUE
  )


data_combined_all$Residual_capture <-
  data_combined_all$Proportion_of_richness -
  data_combined_all$Predicted_capture


# ==============================================================================
# 11. Identify rare attractive plants
# ==============================================================================

get_rare_attr <- function(
    data,
    rare_q = 0.05,
    attr_q = 0.90
) {
  
  rare_thr <- quantile(
    data$Abundance_scaled,
    rare_q,
    na.rm = TRUE
  )
  
  attr_thr <- quantile(
    data$Residual_capture,
    attr_q,
    na.rm = TRUE
  )
  
  data %>%
    
    filter(
      
      Abundance_scaled <
        rare_thr,
      
      Residual_capture >
        attr_thr
      
    )
  
}


Rare_favour_plants <-
  get_rare_attr(
    data_combined_all,
    0.05,
    0.90
  )


# ==============================================================================
# 12. Network-level presence / absence
# ==============================================================================

rare_species <- unique(
  Rare_favour_plants$Plant_species
)


network_highattr <- data_interact %>%
  
  distinct(
    Study_Network_id,
    Plant_accepted_name
  ) %>%
  
  mutate(
    
    is_highattr =
      Plant_accepted_name %in%
      rare_species
    
  )


top10_tbl <- top_10_species %>%
  
  rename(
    Plant_accepted_name = Plant_species
  ) %>%
  
  group_by(
    Study_Network_id
  ) %>%
  
  summarise(
    
    top10_list =
      list(
        unique(
          Plant_accepted_name
        )
      ),
    
    .groups = "drop"
    
  )


network_group <- network_highattr %>%
  
  left_join(
    top10_tbl,
    by = "Study_Network_id"
  ) %>%
  
  group_split(
    Study_Network_id
  ) %>%
  
  purrr::map_dfr(
    
    function(df) {
      
      id <- unique(
        df$Study_Network_id
      )
      
      highattr_species <-
        df$Plant_accepted_name[
          df$is_highattr
        ]
      
      top10 <-
        unique(
          unlist(
            df$top10_list[1]
          )
        )
      
      tibble(
        
        Study_Network_id = id,
        
        present =
          any(
            df$is_highattr,
            na.rm = TRUE
          ),
        
        fully_covered =
          all(
            highattr_species %in%
              top10
          )
        
      )
      
    }
    
  ) %>%
  
  mutate(
    
    group = ifelse(
      present &
        !fully_covered,
      "Present",
      "Absent"
    )
    
  )


# ==============================================================================
# 13. Merge with Top10 capture
# ==============================================================================

attract_percent_10 <-
  percent_10 %>%
  
  rename(
    percentage_Abun10 = percentage
  ) %>%
  
  left_join(
    
    network_group %>%
      dplyr::select(
        Study_Network_id,
        group
      ),
    
    by = "Study_Network_id"
    
  ) %>%
  
  rename(
    present_group = group
  ) %>%
  
  mutate(
    
    present_group =
      str_squish(
        present_group
      ),
    
    present_group =
      str_to_title(
        present_group
      ),
    
    present_group =
      factor(
        present_group,
        levels = c(
          "Present",
          "Absent"
        )
      )
    
  )


# ==============================================================================
# 14. Presence / absence statistics
# ==============================================================================

wilcox_res <- wilcox.test(
  
  percentage_Abun10 ~ present_group,
  
  data = attract_percent_10
  
)


sig_label <- if (
  wilcox_res$p.value < 0.001
) {
  
  "***"
  
} else if (
  wilcox_res$p.value < 0.01
) {
  
  "**"
  
} else if (
  wilcox_res$p.value < 0.05
) {
  
  "*"
  
} else {
  
  "ns"
  
}


data_summary <- attract_percent_10 %>%
  
  group_by(
    present_group
  ) %>%
  
  summarise(
    
    Mean_Percentage =
      mean(
        percentage_Abun10,
        na.rm = TRUE
      ),
    
    SD_Percentage =
      sd(
        percentage_Abun10,
        na.rm = TRUE
      ),
    
    n =
      sum(
        !is.na(
          percentage_Abun10
        )
      ),
    
    .groups = "drop"
    
  )


# ==============================================================================
# 15. Common Y-axis
# ==============================================================================

# All three panels use the same Y-axis:
# Percent of pollinator richness captured (%)
#
# Fixed range: 0–100
# Only Panel A displays the Y-axis.

common_ylim <- c(0, 100)

common_breaks <- seq(
  0,
  100,
  by = 20
)


# ==============================================================================
# 16. Harmonized colors
# ==============================================================================

# ------------------------------------------------------------------------------

# Plant richness: harmonious green

# Pollinator richness: harmonious yellow

# Nestedness: warm orange

# Presence / absence: purple–grey

# ------------------------------------------------------------------------------

richness_colors <- c(
  
  "Plant richness" =
    "#69B89C",
  
  "Pollinator richness" =
    "#E6B84A"
  
)

nestedness_color <-
  "#5B8FD9"


# ------------------------------------------------------------------------------

# Presence / absence colors

# ------------------------------------------------------------------------------

presence_colors <- c(
  
  "Present" =
    "#8064A2",
  
  "Absent" =
    "#AFA9B8"
  
)

# Slightly darker colors for individual points
# This keeps points visible but still within the same color family.

presence_point_colors <- c(
  
  "Present" =
    "#604477",
  
  "Absent" =
    "#77717D"
  
)


# ==============================================================================
# 17. PANEL A
# Plant + Pollinator richness
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
  
  # --------------------------------------------------------------------------
# Points
# --------------------------------------------------------------------------

geom_point(
  
  size = 2.4,
  
  alpha = 0.75,
  
  position =
    position_jitter(
      width = 0.03,
      height = 0.4
    )
  
) +
  
  # --------------------------------------------------------------------------
# Regression lines
# --------------------------------------------------------------------------

geom_smooth(
  
  aes(
    color = Predictor
  ),
  
  method = "lm",
  
  se = FALSE,
  
  linewidth = 0.9
  
) +
  
  # --------------------------------------------------------------------------
# Colors
# --------------------------------------------------------------------------

scale_color_manual(
  
  values =
    richness_colors
  
) +
  
  # --------------------------------------------------------------------------
# New symbols
#
# Avoid the previously used:
# 1 = open circle
# 2 = open triangle
#
# Here:
# 18 = filled diamond
# 15 = filled square
# --------------------------------------------------------------------------

scale_shape_manual(
  
  values = c(
    
    "Plant richness" =
      18,
    
    "Pollinator richness" =
      15
    
  )
  
) +
  
  # --------------------------------------------------------------------------
# Y-axis
# --------------------------------------------------------------------------

scale_y_continuous(
  
  limits =
    common_ylim,
  
  breaks =
    common_breaks,
  
  expand =
    expansion(
      mult = c(
        0,
        0.02
      )
    )
  
) +
  
  # --------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------

labs(
  
  x = "Richness",
  
  y =
    "Percent of pollinator richness captured (%)",
  
  color = NULL,
  
  shape = NULL
  
) +
  
  # --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------

annotate(
  
  "text",
  
  x = Inf,
  
  y = 97,
  
  hjust = 1.05,
  
  vjust = 1,
  
  label =
    paste0(
      
      "Plant: R = ",
      
      round(
        stats_richness$r[
          stats_richness$Predictor ==
            "Plant richness"
        ],
        2
      ),
      
      ", p = ",
      
      format.pval(
        stats_richness$p[
          stats_richness$Predictor ==
            "Plant richness"
        ],
        digits = 2
      ),
      
      "\n",
      
      "Pollinator: R = ",
      
      round(
        stats_richness$r[
          stats_richness$Predictor ==
            "Pollinator richness"
        ],
        2
      ),
      
      ", p = ",
      
      format.pval(
        stats_richness$p[
          stats_richness$Predictor ==
            "Pollinator richness"
        ],
        digits = 2
      )
      
    ),
  
  size = 3.1
  
) +
  
  # --------------------------------------------------------------------------
# Theme
# --------------------------------------------------------------------------

theme_classic(
  
  base_size = 12
  
) +
  
  theme(
    
    axis.title.x =
      element_text(
        face = "bold",
        size = 13
      ),
    
    axis.title.y =
      element_text(
        face = "bold",
        size = 13
      ),
    
    axis.text =
      element_text(
        color = "black",
        size = 11
      ),
    
    legend.position =
      "top",
    
    legend.text =
      element_text(
        size = 10
      ),
    
    legend.key.width =
      unit(
        1.1,
        "cm"
      ),
    
    plot.margin =
      margin(
        5,
        5,
        5,
        5
      )
    
  )


# ==============================================================================
# 18. PANEL B
# Network nestedness
# ==============================================================================

p_nestedness <- ggplot(
  
  df_nested,
  
  aes(
    
    x = Nestedness,
    
    y = percentage_Abun10
    
  )
  
) +
  
  # --------------------------------------------------------------------------
# Points
# --------------------------------------------------------------------------

geom_point(
  
  shape = 16,
  
  size = 2.4,
  
  alpha = 0.70,
  
  color =
    nestedness_color,
  
  position =
    position_jitter(
      width = 0.01,
      height = 0.4
    )
  
) +
  
  # --------------------------------------------------------------------------
# Regression
# --------------------------------------------------------------------------

geom_smooth(
  
  method = "lm",
  
  se = FALSE,
  
  linewidth = 0.9,
  
  color =
    nestedness_color,
  
  linetype =
    nested_line
  
) +
  
  # --------------------------------------------------------------------------
# Y-axis
#
# Keep the scale identical to Panel A,
# but hide the Y-axis itself.
# --------------------------------------------------------------------------

scale_y_continuous(
  
  limits =
    common_ylim,
  
  breaks =
    common_breaks,
  
  expand =
    expansion(
      mult = c(
        0,
        0.02
      )
    )
  
) +
  
  # --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------

annotate(
  
  "text",
  
  x = Inf,
  
  y = 97,
  
  hjust = 1.05,
  
  vjust = 1,
  
  label =
    nested_label,
  
  size = 3.1
  
) +
  
  # --------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------

labs(
  
  x =
    "Network nestedness",
  
  y = NULL
  
) +
  
  # --------------------------------------------------------------------------
# Theme
# --------------------------------------------------------------------------

theme_classic(
  
  base_size = 12
  
) +
  
  theme(
    
    axis.title.x =
      element_text(
        face = "bold",
        size = 13
      ),
    
    # Hide Y-axis completely
    
    axis.title.y =
      element_blank(),
    
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    axis.text.x =
      element_text(
        color = "black",
        size = 11
      ),
    
    legend.position =
      "none",
    
    plot.margin =
      margin(
        5,
        5,
        5,
        5
      )
    
  )


# ==============================================================================
# 19. PANEL C
# Presence / absence of rare attractive plants
# ==============================================================================

# ------------------------------------------------------------------------------

# Significance bracket

# ------------------------------------------------------------------------------

y_pos_c <- 88


p_presence <- ggplot(
  
  attract_percent_10,
  
  aes(
    
    x = present_group,
    
    y = percentage_Abun10,
    
    fill = present_group
    
  )
  
) +
  
  # --------------------------------------------------------------------------
# Violin
# --------------------------------------------------------------------------

geom_violin(
  
  trim = FALSE,
  
  alpha = 0.25,
  
  color = NA,
  
  width = 0.8
  
) +
  
  # --------------------------------------------------------------------------
# Boxplot
# --------------------------------------------------------------------------

geom_boxplot(
  
  width = 0.16,
  
  alpha = 0.65,
  
  outlier.shape = NA,
  
  color = "grey25"
  
) +
  
  # --------------------------------------------------------------------------
# Points
#
# IMPORTANT:
# The points now use the same purple-grey family
# as the violin.
#
# Present = dark purple
# Absent  = dark grey-purple
# --------------------------------------------------------------------------

geom_jitter(
  
  aes(
    color = present_group
  ),
  
  width = 0.10,
  
  size = 1.6,
  
  alpha = 0.55
  
) +
  
  # --------------------------------------------------------------------------
# Mean
#
# FIXED POSITION:
# No longer attached to the actual mean.
# Therefore it will not be covered by points.
# --------------------------------------------------------------------------

geom_text(
  
  data =
    data_summary,
  
  aes(
    
    x = present_group,
    
    y = 72,
    
    label =
      sprintf(
        "Mean = %.1f%%",
        Mean_Percentage
      )
    
  ),
  
  inherit.aes = FALSE,
  
  size = 3.1,
  
  fontface = "bold"
  
) +
  
  # --------------------------------------------------------------------------
# n
#
# Also fixed at a clean position.
# --------------------------------------------------------------------------

geom_text(
  
  data =
    data_summary,
  
  aes(
    
    x = present_group,
    
    y = 78,
    
    label =
      paste0(
        "n = ",
        n
      )
    
  ),
  
  inherit.aes = FALSE,
  
  size = 3.0
  
) +
  
  # --------------------------------------------------------------------------
# Significance bracket
# --------------------------------------------------------------------------

geom_segment(
  
  aes(
    
    x = 1,
    
    xend = 2,
    
    y = y_pos_c,
    
    yend = y_pos_c
    
  ),
  
  inherit.aes = FALSE
  
) +
  
  geom_segment(
    
    aes(
      
      x = 1,
      
      xend = 1,
      
      y = y_pos_c,
      
      yend =
        y_pos_c - 2
      
    ),
    
    inherit.aes = FALSE
    
  ) +
  
  geom_segment(
    
    aes(
      
      x = 2,
      
      xend = 2,
      
      y = y_pos_c,
      
      yend =
        y_pos_c - 2
      
    ),
    
    inherit.aes = FALSE
    
  ) +
  
  annotate(
    
    "text",
    
    x = 1.5,
    
    y =
      y_pos_c + 2,
    
    label =
      sig_label,
    
    size = 5,
    
    fontface = "bold"
    
  ) +
  
  # --------------------------------------------------------------------------
# Fill colors
# --------------------------------------------------------------------------

scale_fill_manual(
  
  values =
    presence_colors
  
) +
  
  # --------------------------------------------------------------------------
# Point colors
# --------------------------------------------------------------------------

scale_color_manual(
  
  values =
    presence_point_colors
  
) +
  
  # --------------------------------------------------------------------------
# Y-axis
#
# Same 0–100 range.
# Completely hidden in Panel C.
# --------------------------------------------------------------------------

scale_y_continuous(
  
  limits =
    common_ylim,
  
  breaks =
    common_breaks,
  
  expand =
    expansion(
      mult = c(
        0,
        0.02
      )
    )
  
) +
  
  # --------------------------------------------------------------------------
# Labels
# --------------------------------------------------------------------------

labs(
  
  x =
    "Rare attractive plants in full network",
  
  y = NULL,
  
  fill = NULL,
  
  color = NULL
  
) +
  
  # --------------------------------------------------------------------------
# Theme
# --------------------------------------------------------------------------

theme_classic(
  
  base_size = 12
  
) +
  
  theme(
    
    axis.title.x =
      element_text(
        face = "bold",
        size = 13
      ),
    
    # Hide Y-axis completely
    
    axis.title.y =
      element_blank(),
    
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    axis.text.x =
      element_text(
        color = "black",
        size = 11
      ),
    
    legend.position =
      "none",
    
    plot.margin =
      margin(
        5,
        5,
        5,
        5
      )
    
  )


# ==============================================================================
# 20. Inspect individual panels
# ==============================================================================

p_richness

p_nestedness

p_presence


# ==============================================================================
# 21. Combine into three-panel figure
# ==============================================================================

final_factor_fig <- plot_grid(
  
  p_richness,
  
  p_nestedness,
  
  p_presence,
  
  ncol = 3,
  
  labels =
    c(
      "a",
      "b",
      "c"
    ),
  
  label_size = 14,
  
  label_fontface =
    "bold",
  
  align = "h",
  
  axis = "tb",
  
  rel_widths =
    c(
      1.05,
      1,
      0.95
    )
  
)


# ==============================================================================
# 22. Display final figure
# ==============================================================================

final_factor_fig


# ==============================================================================
# 23. Save final figure
# ==============================================================================

ggsave(
  
  "/Chap1_TargetPlant_to_monitor/result_260723/Fig_Factors_subsampling_effectiveness.png",
  
  final_factor_fig,
  
  width = 13,
  
  height = 4.6,
  
  units = "in",
  
  dpi = 600,
  
  bg = "white"
  
)


# ==============================================================================
#  Display final figure
# ==============================================================================

final_factor_fig

