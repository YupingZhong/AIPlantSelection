
##-----------(3) hypothesis 3 : Different flower shapes attract different pollinators groups 
# 花型对传粉者组成的分析（268 sites not pooled）genus级别

##########################################
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(betapart)
library(multcompView)
library(flextable)
library(officer)
library(cowplot)

data_merge<-readRDS("data/processed/data_merge.rds")

all_data <- data_merge %>%
  filter(
    !is.na(flw_shape_revised),
    !is.na(Pollinator_genus),
    !flw_shape_revised %in% c("trap flowers", "brush flowers")
  )


# ---------------------------
# 4.1 Matrix construction
## 4.1.1 Plant × Site × Pollinator genus Matrix（presence/absence）
# ---------------------------（选项一）
# 每个 Plant × Site × Pollinator genus 只算一次
all_data_unique <- all_data %>%
  distinct(Plant_accepted_name, Study_Network_id, Pollinator_genus, .keep_all = TRUE)


mat <- all_data_unique %>%
  mutate(presence = 1) %>%   # 转为 presence/absence
  group_by(Plant_accepted_name, Study_Network_id, Pollinator_genus) %>%
  summarise(presence = max(presence), .groups = "drop") %>%
  pivot_wider(names_from = Pollinator_genus,
              values_from = presence,
              values_fill = 0)



plant_flw <- all_data_unique %>%
  dplyr::select(Plant_accepted_name, flw_shape_revised, Study_Network_id,Study_id) %>%
  distinct()

mat_numeric <- mat %>%
  mutate(SampleID = paste(Plant_accepted_name,
                          Study_Network_id, sep = "_")) %>%
  tibble::column_to_rownames("SampleID") %>%
  dplyr::select(-Plant_accepted_name, -Study_Network_id)

meta <- plant_flw %>%
  mutate(
    SampleID = paste(Plant_accepted_name, Study_Network_id, sep = "_")
  )

mat_numeric2 <- mat_numeric[meta$SampleID, ]


# ---------------------------
# 4.2.1 PERMANOVA
# ---------------------------
d <- vegdist(mat_numeric2, method = "jaccard")  # frequency用method = "bray"

adonis_res <- adonis2(d ~ Study_id +flw_shape_revised,
                                 data = meta,
                                 permutations = 999,
                                 strata = meta$Study_Network_id,
                                 by = "margin")

head(meta)
print(adonis_res)
saveRDS(adonis_res,"data/processed/adonis_res.rds")

#PERMANOVA结果显著，但R方极小
#这说明globally, 花型对传粉者群落组成的影响显著，但不是主要影响因素
#但如果单独对每一个网络进行分析，结果可能会不同

# ==========================================================
# Pairwise PERMANOVA between flower-shape categories
# ==========================================================

library(permute)

shape_pairs <- combn(
  unique(meta$flw_shape_revised),
  2,
  simplify = FALSE
)

pairwise_permanova <- lapply(shape_pairs, function(pair) {
  
  # Select the two flower shapes
  keep <- meta$flw_shape_revised %in% pair
  
  meta_pair <- droplevels(meta[keep, ])
  mat_pair  <- mat_numeric2[keep, , drop = FALSE]
  
  # Restricted permutations within networks
  permutation_control <- how(nperm = 999)
  setBlocks(permutation_control) <- meta_pair$Study_Network_id
  
  fit <- tryCatch(
    adonis2(
      mat_pair ~ Study_id + flw_shape_revised,
      data = meta_pair,
      method = "jaccard",
      binary = TRUE,
      permutations = permutation_control,
      by = "margin"
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(
      data.frame(
        Shape_1 = pair[1],
        Shape_2 = pair[2],
        F_value = NA_real_,
        R2 = NA_real_,
        P_value = NA_real_
      )
    )
  }
  
  shape_row <- which(
    rownames(fit) == "flw_shape_revised"
  )
  
  data.frame(
    Shape_1 = pair[1],
    Shape_2 = pair[2],
    F_value = fit$F[shape_row],
    R2 = fit$R2[shape_row],
    P_value = fit$`Pr(>F)`[shape_row]
  )
})

pairwise_permanova <- bind_rows(pairwise_permanova) %>%
  mutate(
    P_adjusted = p.adjust(
      P_value,
      method = "BH"
    ),
    Significance = case_when(
      is.na(P_adjusted) ~ NA_character_,
      P_adjusted < 0.001 ~ "***",
      P_adjusted < 0.01  ~ "**",
      P_adjusted < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(P_adjusted)

print(pairwise_permanova)

#############

# (2)组间差异主要由大部分还是少部分网络驱动

###########
adonis_res<-readRDS("data/processed/adonis_res.rds")

networks <- unique(meta$Study_Network_id)

res_list <- list()

skip_list <- data.frame(
  Network = character(),
  Reason = character()
)


for(net in networks){
  
  sub_meta <- meta %>%
    filter(Study_Network_id == net)
  
  keep_ids <- sub_meta$SampleID
  
  sub_mat <- mat_numeric2[keep_ids, , drop = FALSE]
  
  
  # 少于两个花型
  if(length(unique(sub_meta$flw_shape_revised)) < 2){
    skip_list <- rbind(
      skip_list,
      data.frame(
        Network = net,
        Reason = "Less than two flower shapes"
      )
    )
    next
  }
  
  
  # 样本太少
  if(nrow(sub_mat) < 5){
    skip_list <- rbind(
      skip_list,
      data.frame(
        Network = net,
        Reason = "Too few samples"
      )
    )
    next
  }
  
  
  # 没有interaction
  if(sum(sub_mat) == 0){
    skip_list <- rbind(
      skip_list,
      data.frame(
        Network = net,
        Reason = "No pollinator interactions"
      )
    )
    next
  }
  
  
  sub_d <- vegdist(sub_mat,
                   method = "jaccard")
  
  
  # 所有样本组成完全一样
  if(sum(sub_d) == 0){
    skip_list <- rbind(
      skip_list,
      data.frame(
        Network = net,
        Reason = "No variation in pollinator composition"
      )
    )
    next
  }
  
  
  ad <- tryCatch(
    adonis2(
      sub_d ~ flw_shape_revised,
      data = sub_meta,
      permutations = 999
    ),
    error = function(e) NULL
  )
  
  
  if(is.null(ad)){
    skip_list <- rbind(
      skip_list,
      data.frame(
        Network = net,
        Reason = "PERMANOVA failed"
      )
    )
    next
  }
  
  
  res_list[[as.character(net)]] <- data.frame(
    
    Network = net,
    
    # 样本数量
    N_samples = nrow(sub_meta),
    
    # 植物物种数量
    N_plants = n_distinct(sub_meta$Plant_accepted_name),
    
    # 花型数量
    N_shapes = n_distinct(sub_meta$flw_shape_revised),
    
    # interaction richness
    N_interactions = sum(sub_mat),
    
    # pollinator genus richness
    N_pollinators = sum(colSums(sub_mat) > 0),
    
    # PERMANOVA结果
    R2 = ad$R2[1],
    F = ad$F[1],
    P = ad$`Pr(>F)`[1]
  )
}


network_perm <- bind_rows(res_list)


network_perm %>%
  summarise(
    total_network = n(),
    significant = sum(P < 0.05),
    proportion = mean(P < 0.05)
  )



#--------
# 4.2.2 dispersion difference test
# If：P > 0.05，PERMANOVA result reliable
#-------
disp <- betadisper(d, meta$flw_shape_revised)

saveRDS(disp, "data/processed/disp.rds")

R2_strata <- round(adonis_res$R2[1], 3)
P_strata <- adonis_res$`Pr(>F)`[1]
P_strata_text <- ifelse(P_strata < 0.001, "P < 0.001", paste0("P = ", round(P_strata, 3)))

anova(disp)

# ############################################
# Figure A
# Variation in pollinator community composition among flower shapes
# betadisper + sidak
# ############################################

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)

disp <- readRDS("data/processed/disp.rds")

disp_df <- data.frame(
  Distance = disp$distances,
  FlowerShape = meta$flw_shape_revised,
  Study_id = meta$Study_id,
  Study_Network_id = meta$Study_Network_id
) %>%
  filter(!is.na(FlowerShape)) %>%
  mutate(
    Study_id = factor(Study_id),
    Study_Network_id = factor(Study_Network_id),
    FlowerShape = factor(FlowerShape)
  )

disp_model <- lmer(
  Distance ~ Study_id + FlowerShape + (1 | Study_Network_id),
  data = disp_df
)

anova(disp_model)

emm <- emmeans(disp_model, ~ FlowerShape)

# 先定义绘图排序：由小到大
shape_order <- as.data.frame(emm) %>%
  arrange(emmean) %>%
  pull(FlowerShape)

# sidak 字母分组
letters_df <- multcomp::cld(
  emm,
  adjust = "sidak",
  Letters = letters,
  sort = FALSE
) %>%
  as.data.frame() %>%
  dplyr::select(FlowerShape, .group) %>%
  rename(Letters = .group) %>%
  mutate(
    Letters = trimws(Letters),
    FlowerShape = factor(FlowerShape, levels = shape_order)
  )

# 用于作图的 estimated marginal means
emm_df <- as.data.frame(emm) %>%
  mutate(
    FlowerShape = factor(FlowerShape, levels = shape_order)
  ) %>%
  left_join(letters_df, by = "FlowerShape")

# 两两比较结果
pairs_res <- pairs(emm, adjust = "sidak")

# 若后续箱线图/散点图使用这些对象，再统一排序
disp_df <- disp_df %>%
  mutate(FlowerShape = factor(FlowerShape, levels = shape_order))


############################################
# B. Beta diversity partition
# Study → Network → Plant species
############################################

# -----------------------------------------
# 1. Presence / absence matrix
# -----------------------------------------

comm_pa <- as.data.frame(mat_numeric2)

comm_pa[] <- lapply(
  comm_pa,
  function(x) as.integer(x > 0)
)

comm_pa <- as.matrix(comm_pa)


# -----------------------------------------
# 2. Metadata + unique SampleID
# -----------------------------------------

meta_beta <- meta

# -----------------------------------------
# 3. Match community matrix and metadata
# -----------------------------------------

common_ids <- intersect(
  rownames(comm_pa),
  meta_beta$SampleID
)

comm_pa <- comm_pa[
  common_ids,
  ,
  drop = FALSE
]

meta_beta <- meta_beta %>%
  filter(SampleID %in% common_ids) %>%
  slice(match(common_ids, SampleID))


# 检查是否完全对应
stopifnot(
  identical(
    rownames(comm_pa),
    meta_beta$SampleID
  )
)


# -----------------------------------------
# 4. Beta diversity within each Network
# -----------------------------------------

network_beta_list <- list()

for(net in unique(meta_beta$Study_Network_id)){
  
  # 当前 Network 的 plant species
  ids <- meta_beta$SampleID[
    meta_beta$Study_Network_id == net
  ]
  
  if(length(ids) < 2) next
  
  sub_comm <- comm_pa[
    ids,
    ,
    drop = FALSE
  ]
  
  sub_meta <- meta_beta[
    match(ids, meta_beta$SampleID),
    ,
    drop = FALSE
  ]
  
  
  # ---------------------------------------
  # Sorensen beta diversity
  # ---------------------------------------
  
  bc <- beta.pair(
    sub_comm,
    index.family = "sorensen"
  )
  
  turnover_mat <- as.matrix(bc$beta.sim)
  nestedness_mat <- as.matrix(bc$beta.sne)
  
  
  # ---------------------------------------
  # Within each flower shape
  # ---------------------------------------
  
  for(shp in unique(sub_meta$flw_shape_revised)){
    
    shp_ids <- sub_meta$SampleID[
      sub_meta$flw_shape_revised == shp
    ]
    
    # 至少两个 plant species
    if(length(shp_ids) < 2) next
    
    idx <- match(
      shp_ids,
      rownames(turnover_mat)
    )
    
    if(any(is.na(idx))) next
    
    
    # 同一 flower shape 内的 pairwise beta
    turnover_sub <- turnover_mat[
      idx,
      idx,
      drop = FALSE
    ]
    
    nestedness_sub <- nestedness_mat[
      idx,
      idx,
      drop = FALSE
    ]
    
    turnover_values <- turnover_sub[
      upper.tri(turnover_sub)
    ]
    
    nestedness_values <- nestedness_sub[
      upper.tri(nestedness_sub)
    ]
    
    
    # -------------------------------------
    # Network × Flower shape
    # -------------------------------------
    
    network_beta_list[[length(network_beta_list) + 1]] <-
      data.frame(
        
        Study_id = unique(sub_meta$Study_id),
        
        Study_Network_id = net,
        
        FlowerShape = shp,
        
        turnover = mean(
          turnover_values,
          na.rm = TRUE
        ),
        
        nestedness = mean(
          nestedness_values,
          na.rm = TRUE
        ),
        
        N_plants = length(shp_ids),
        
        N_pairs = length(turnover_values)
      )
  }
}


# -----------------------------------------
# 5. Combine Network × FlowerShape results
# -----------------------------------------

network_beta <- bind_rows(
  network_beta_list
)


# -----------------------------------------
# 6. Average across Networks
# -----------------------------------------

df_beta <- network_beta %>%
  group_by(FlowerShape) %>%
  summarise(
    
    turnover = mean(
      turnover,
      na.rm = TRUE
    ),
    
    nestedness = mean(
      nestedness,
      na.rm = TRUE
    ),
    
    N_networks = n_distinct(
      Study_Network_id
    ),
    
    N_study = n_distinct(
      Study_id
    ),
    
    .groups = "drop"
  )


# -----------------------------------------
# 7. Relative contribution
# -----------------------------------------

df_beta <- df_beta %>%
  mutate(
    
    total = turnover + nestedness,
    
    Turnover = turnover / total,
    
    Nestedness = nestedness / total
    
  ) %>%
  dplyr::select(
    FlowerShape,
    Turnover,
    Nestedness,
    N_networks,
    N_study
  )


# -----------------------------------------
# 8. Long format for plotting
# -----------------------------------------

df_beta_long <- df_beta %>%
  dplyr::select(
    FlowerShape,
    Turnover,
    Nestedness
  ) %>%
  pivot_longer(
    cols = c(
      Turnover,
      Nestedness
    ),
    names_to = "Component",
    values_to = "Proportion"
  )


df_beta$FlowerShape <- factor(df_beta$FlowerShape, levels = shape_order)

df_beta_long$FlowerShape <- factor(
  df_beta_long$FlowerShape,
  levels = shape_order
)

#=========================================

# 4.3 Visialize Pollinator composition across different flower shapes 

#========================================
data_functional <- data_merge %>%
  filter(!is.na(flw_shape_revised),
         !is.na(Pollinator_accepted_name),
         !is.na(Plant_accepted_name),
         !is.na(Study_Network_id),
         Pollinator_genus != "Apis") %>%
  filter(!flw_shape_revised %in% c("brush flowers", "trap flowers")) %>%
  mutate(functional_group = case_when(
    Pollinator_family == "Syrphidae"                          ~ "Syrphidae",
    Pollinator_family %in% c("Apidae", "Halictidae",
                             "Andrenidae", "Megachilidae",
                             "Colletidae", "Melittidae",
                             "Stenotritidae")                ~ "Bees",
    Pollinator_order == "Hymenoptera"                         ~ "Non-bee Hymenoptera",
    Pollinator_order == "Diptera"                             ~ "Non-syrphid Diptera",
    Pollinator_order == "Lepidoptera"                         ~ "Lepidoptera",
    Pollinator_order == "Coleoptera"                          ~ "Coleoptera",
    TRUE                                                      ~ "Other"
  ))



# 检查 Other 里还有什么
data_functional %>%
  filter(functional_group == "Other") %>%
  distinct(Pollinator_order, Pollinator_genus) %>%
  arrange(Pollinator_order)


# ⭐ Step 3: 每个网络内 unique interaction，统计各功能群数量
# 每个真实存在的 network × flower shape 组合
network_shape <- data_functional %>%
  distinct(
    Study_Network_id,
    flw_shape_revised
  )


# 每个 network × flower shape × functional group 的interaction数量
network_fun <- data_functional %>%
  group_by(
    Study_Network_id,
    flw_shape_revised,
    functional_group
  ) %>%
  summarise(
    n_interactions = sum(Interaction_addup,
                         na.rm = TRUE),
    .groups = "drop"
  )


# 只补充真实存在的 flower shape 中缺失的 functional group
network_props <- network_shape %>%
  tidyr::crossing(
    functional_group = unique(data_functional$functional_group)
  ) %>%
  left_join(
    network_fun,
    by = c(
      "Study_Network_id",
      "flw_shape_revised",
      "functional_group"
    )
  ) %>%
  mutate(
    n_interactions = replace_na(n_interactions, 0)
  )

# (option2 )Step 3: 每个网络内 interaction sum，统计各功能群数量
network_props <- network_props %>%
  group_by(
    Study_Network_id,
    flw_shape_revised
  ) %>%
  mutate(
    prop=n_interactions/sum(n_interactions)
  ) %>%
  ungroup()

# step 4: average across network
plot_df <- network_props %>%
  group_by(
    flw_shape_revised,
    functional_group
  ) %>%
  summarise(
    mean_prop=mean(prop),
    .groups="drop"
  )


shape_order <- plot_df %>%
  filter(functional_group=="Bees") %>%
  arrange(mean_prop) %>%
  pull(flw_shape_revised)

###or : visit sum


plot_df$flw_shape_revised <- factor(plot_df$flw_shape_revised, levels = shape_order)



group_cols <- c(
  "Bees"                 = "#F5E066",   
  "Syrphidae"            = "#85CCAE",   
  "Non-bee Hymenoptera"  = "#F5A88A",   
  "Coleoptera"           = "#E8A3D1",   
  "Lepidoptera"          = "#A5B5D9",   
  "Non-syrphid Diptera"  = "#B8DD7F",   
  "Other"                = "#AAAAAA"    
)

# 定义堆叠顺序（从下到上）
group_order <- c(
  "Bees", 
  "Non-bee Hymenoptera", 
  "Syrphidae", 
  "Non-syrphid Diptera", 
  "Lepidoptera", 
  "Coleoptera", 
  "Other"
)

# 设置因子水平
plot_df$functional_group <- factor(plot_df$functional_group, levels = group_order)

# 设置因子顺序
#plot_df$flw_shape_revised <- factor(plot_df$flw_shape_revised,
#                                          levels = shape_order_final)

# plot
visit_group <- ggplot(plot_df,
                      aes(x = flw_shape_revised, y = mean_prop, fill = functional_group)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0,1),
    breaks = c(0,0.25,0.5,0.75,1),
    labels = c("0%","25%","50%","75%","100%"),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  labs(x = "Flower shape categories", y = "Mean proportion of pollinator visits", fill = "Pollinator group") +
  theme_classic(base_size = 10) +
  theme(
    panel.border = element_rect(
      colour = "grey55",
      fill = NA,
      linewidth = 0.9
    ),
    axis.line = element_blank(),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.text.y = element_text(color = "black", size = 9),
    axis.text.x = element_text(color = "black", size = 8),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 9),
    legend.position = "right",
    legend.direction = "vertical",
    plot.margin = margin(3, 5, 3, 3)
  )

visit_group

#ggsave("/Chap1_TargetPlant_to_monitor/result_260723/visit_larger_group_sum.png", visit_group, width = 5, 
#       height = 3.0, units = "in", dpi = 600, bg = "white")  

# 不带图例的主体图
visit_group_no_legend <- visit_group +
  theme(legend.position = "none")

# 单独提取图例：底部、两行
visit_group_legend <- cowplot::get_legend(
  visit_group +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical"
    ) +
    guides(
      fill = guide_legend(
        nrow = 2,
        byrow = TRUE,
        title.position = "top"
      )
    )
)

# 分别保存
ggsave(
  "/Chap1_TargetPlant_to_monitor/result_260723/visit_group_no_legend.png",
  visit_group_no_legend,
  width = 5, height = 3,
  units = "in", dpi = 600, bg = "white"
)

ggsave(
  "/Chap1_TargetPlant_to_monitor/result_260723/visit_group_legend.png",
  visit_group_legend,
  width = 8, height = 1.1,
  units = "in", dpi = 600, bg = "white"
)


###################################################################
# Panel A
############################################
# ==========================================================
# 让 Panel A 的 flower-shape 顺序跟随 visit_group
# ==========================================================

shared_shape_order <- levels(plot_df$flw_shape_revised)

disp_df <- disp_df %>%
  filter(FlowerShape %in% shared_shape_order) %>%
  mutate(
    FlowerShape = factor(
      as.character(FlowerShape),
      levels = shared_shape_order
    )
  )

emm_df <- emm_df %>%
  filter(FlowerShape %in% shared_shape_order) %>%
  mutate(
    FlowerShape = factor(
      as.character(FlowerShape),
      levels = shared_shape_order
    )
  )

letters_df <- letters_df %>%
  filter(FlowerShape %in% shared_shape_order) %>%
  mutate(
    FlowerShape = factor(
      as.character(FlowerShape),
      levels = shared_shape_order
    )
  )

label_pos <- disp_df %>%
  group_by(FlowerShape) %>%
  summarise(
    xpos = max(Distance),
    .groups = "drop"
  ) %>%
  left_join(
    letters_df,
    by = "FlowerShape"
  ) %>%
  mutate(
    xpos = xpos * 1.08
  )


pA <- ggplot(
  disp_df,
  aes(x = FlowerShape, y = Distance)
) +
  
  # 原始数据
  geom_boxplot(
    fill = "#7b95c6",
    alpha = 0.45,
    width = 0.65,
    outlier.shape = NA,
    linewidth = 0.35
  ) +
  
  geom_jitter(
    width = 0.12,
    size = 0.55,
    color = "grey60",
    alpha = 0.18
  ) +
  
  # adjusted mean 的 95% CI
  geom_errorbar(
    data = emm_df,
    aes(
      x = FlowerShape,
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.45
  ) +
  
  # adjusted mean
  geom_point(
    data = emm_df,
    aes(
      x = FlowerShape,
      y = emmean
    ),
    inherit.aes = FALSE,
    size = 2.3,
    shape = 21,
    fill = "white",
    stroke = 0.8
  ) +
  
  # Tukey letters
  geom_text(
    data = emm_df,
    aes(
      x = FlowerShape,
      y = asymp.UCL + 0.025,
      label = Letters
    ),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  
  coord_flip(clip = "off") +
  
  labs(
    x = NULL,
    y = "Pollinator community dispersion"
  ) +
  
  theme_classic(base_size = 10) +
  theme(
    panel.border = element_rect(
      colour = "grey55",
      fill = NA,
      linewidth = 0.5
    ),
    axis.line = element_blank(),
    axis.title.x = element_text(
      face = "bold",
      size = 11
    ),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 9
    ),
    plot.margin = margin(3, 5, 3, 3)
  )

pA


ggsave(
  "/Chap1_TargetPlant_to_monitor/result_260723/flower_shape_beta_partition_A.png",
  pA,
  width = 2.8,
  height = 3.0,
  dpi = 600,
  bg = "white"
)



# ############################################
# # Panel B
# ############################################
# 
# pB <- ggplot(
#   df_beta_long,
#   aes(
#     x = FlowerShape,
#     y = Proportion,
#     fill = Component
#   )
# ) +
#   geom_col(width = 0.65) +
#   coord_flip() +
#   scale_y_continuous(
#     limits = c(0, 1),
#     breaks = c(0, 0.25, 0.5, 0.75, 1),
#     labels = c("0%", "25%", "50%", "75%", "100%"),
#     expand = c(0, 0)
#   ) +
#   scale_fill_manual(
#     values = c(
#       "Turnover" = "#a1d8e8",
#       "Nestedness" = "#a2c986"
#     )
#   ) +
#   labs(
#     x = NULL,
#     y = "Relative contribution to β-diversity",
#     fill = NULL
#   ) +
#   theme_classic(base_size = 10) +
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.line.y = element_blank(),
#     axis.title.x = element_text(face = "bold", size = 11),
#     axis.text.x = element_text(color = "black", size = 8),
#     axis.ticks.x = element_line(linewidth = 0.3),
#     axis.line.x = element_line(linewidth = 0.4),
#     legend.position = "right",
#     legend.direction = "vertical",
#     legend.text = element_text(size = 9),
#     plot.margin = margin(3, 5, 3, 3)
#   )
# ############################################
# # Combine
# ############################################
# figure_all <- cowplot::plot_grid(
#   pA,
#   pB,
#   nrow = 1,
#   rel_widths = c(1.8, 1),
#   align = "h"
# )
# 
# figure_all
# 
# ggsave(
#   "/Chap1_TargetPlant_to_monitor/result_260723/flower_shape_beta_partition_AB.png",
#   figure_all,
#   width = 8.5,
#   height = 3.0,
#   dpi = 600,
#   bg = "white"
# )
# 


###############

#Table output

###############

sp_list_rare_info<-sp_list_rare%>%
  left_join(data_interact%>%dplyr::select(Plant_accepted_name,Plant_order,Plant_family,Plant_genus)%>%
              distinct(), by = "Plant_accepted_name")

sp_count_shape <- sp_list_rare_info %>%
  group_by(flw_shape_revised) %>%
  summarise(
    n_plants = n_distinct(Plant_accepted_name)
  )


sp_count_family <- sp_list_rare_info %>%
  group_by(Plant_family) %>%
  summarise(
    n_plants = n_distinct(Plant_accepted_name)
  ) %>%
  arrange(desc(n_plants))


# PERMANOVA
perm_table <- data.frame(
  Test = "PERMANOVA",
  Factor = c("Flower shape", "Residual"),
  Df = adonis_res$Df[1:2],
  F = c(round(adonis_res$F[1], 2), NA),
  R2 = round(adonis_res$R2[1:2], 3),
  P_value = c(
    ifelse(adonis_res$`Pr(>F)`[1] < .001, "<0.001",
           sprintf("%.3f", adonis_res$`Pr(>F)`[1])),
    NA
  )
)

# Betadisper
disp_aov <- anova(disp)

disp_table <- data.frame(
  Test = "Betadisper",
  Factor = rownames(disp_aov),
  Df = disp_aov$Df,
  F = round(disp_aov$`F value`, 2),
  R2 = NA,
  P_value = ifelse(
    disp_aov$`Pr(>F)` < .001, "<0.001",
    sprintf("%.3f", disp_aov$`Pr(>F)`)
  )
)

# 合并
final_table <- bind_rows(perm_table, disp_table)

# Flextable
ft <- flextable(final_table) |>
  theme_booktabs() |>
  set_header_labels(
    Test = "Test", Factor = "Factor", Df = "Df",
    F = "F", R2 = "R²", P_value = "P-value"
  ) |>
  font(fontname = "Arial", part = "all") |>
  fontsize(size = 10, part = "all") |>
  bold(part = "header") |>
  align(j = c("Test", "Factor"), align = "left", part = "all") |>
  align(j = c("Df", "F", "R2", "P_value"), align = "center", part = "all") |>
  autofit()

# Word
read_docx() |>
  body_add_par(
    "Table X. PERMANOVA and multivariate dispersion (betadisper) results for pollinator community composition across flower shapes.",
    style = "heading 2"
  ) |>
  body_add_flextable(ft) |>
  print(target = "./result_260526/community_shape_multivariate_test.docx")