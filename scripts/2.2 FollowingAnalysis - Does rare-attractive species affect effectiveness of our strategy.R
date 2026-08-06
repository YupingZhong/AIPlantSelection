###############################################################################

##-----------(3) hypothesis 3 : Favored plant sp by pollinators is rare (plant abundance-interaction times) 

################################
library(ggplot2)
library(dplyr)
library(broom)
library(purrr)

traits<-read.csv("data/processed/merge.trait.csv", header = TRUE, fileEncoding = "UTF-8")
colnames(data_count_scaled)

top_10_species <- data_count_scaled %>%
  distinct(Study_Network_id, Plant_species, Flower_count_scaled, .keep_all = FALSE) %>%
  group_by(Study_Network_id) %>%
  slice_max(order_by = Flower_count_scaled, n = 10, with_ties = FALSE) %>%  # Keep only the top 10 species per network (no ties allowed)
  ungroup()%>%
  group_by(Study_Network_id) %>%
  filter(n_distinct(Plant_species) > 9) %>%
  ungroup()%>%
  dplyr::select(Study_Network_id, Plant_species)


######################全部网络
# 准备数据
# 先算出每个网络总互动次数
interaction_all <- data_interact %>%
  group_by(Study_Network_id) %>%
  filter(Interaction_addup!= 0)%>%
  summarize(Interaction_sum = sum(Interaction_addup),
            richness_sum = n_distinct(Pollinator_accepted_name), .groups = "drop")

# 合并后继续处理
match<-data_interact%>%
  ungroup%>%
  dplyr::select(Plant_accepted_name,Plant_original_name)%>%
  distinct()%>%
  rename(Plant_species = Plant_original_name)

unique(data_combined_all$Study_Network_id[!is.na(data_combined_all$Study_Network_id)])
unique(interaction_all$Study_Network_id)
unique(data_merge$Study_Network_id)

data_combined_all <- data_merge %>% 
  ungroup%>%#
  left_join(interaction_all, by = "Study_Network_id") %>%
  left_join(traits,by = "Plant_accepted_name")%>%
  group_by(Study_id,Study_Network_id, Plant_species) %>%
  summarize(
    Abundance = first(Flower_count_scaled),
    Pollinator_richness = n_distinct(Pollinator_accepted_name, na.rm = TRUE), #对每一种植物，计算累积的pollinator物种数
    Interaction_times = sum(Interaction_addup, na.rm = TRUE), #对每一种植物，计算累积的pollinator交互次数
    Proportion_of_interactions = (Interaction_times / first(Interaction_sum))*100,
    Proportion_of_richness = (Pollinator_richness / first(richness_sum))*100,
    .groups = "drop"
  ) %>%
  filter(is.finite(Abundance) & is.finite(Interaction_times)) %>%
  group_by(Study_Network_id) %>%
  mutate(
    #Abundance_normalized = (Abundance - min(Abundance)) / (max(Abundance) - min(Abundance)) * 100, ##每个植物在本网络中的相对多度
    Abundance_scaled = scale(log(Abundance + 1))[,1]
  ) %>%
  ungroup()%>%
  left_join(match,by="Plant_species")

data_combined_all%>%filter(is.na(Abundance))


# fit mix model
library(lme4)

model <- lmer(
  Proportion_of_richness ~ Abundance_scaled +
    (1 | Study_id/Study_Network_id),
  data = data_combined_all
)

summary(model)


data_combined_all$Predicted_capture <- predict(
  model,
  newdata = data_combined_all,
  re.form = NA,
  allow.new.levels = TRUE
)


############################################
###  3.2   Plotting
##3.2.1 总图
#########（1）richness
colnames(data_combined_all)

get_rare_attr <- function(data, rare_q = 0.05, attr_q = 0.9) {
  
  rare_thr <- quantile(data$Abundance_scaled, rare_q, na.rm = TRUE)
  attr_thr <- quantile(data$Residual_capture, attr_q, na.rm = TRUE)
  
  data %>%
    filter(
      Abundance_scaled < rare_thr,
      Residual_capture > attr_thr
    )
}



data_combined_all <- data_combined_all %>%
  mutate(
    Residual_capture =
      Proportion_of_richness - Predicted_capture
  )%>%
  mutate(
    performance_class = case_when(
      Abundance_scaled < quantile(Abundance_scaled, 0.05, na.rm = TRUE) &
        Residual_capture > quantile(Residual_capture, 0.9, na.rm = TRUE)
      ~ "Rare and attractive",
      
      Residual_capture > 0 ~ "Above expected",
      TRUE ~ "Below expected"
    )
  )



Rare_favour_plants <- get_rare_attr(data_combined_all, 0.05, 0.9)

unique(Rare_favour_plants$Plant_species)


## calculate R2 and p
# 固定效应预测
library(lmerTest)

model <- lmer(
  Proportion_of_richness ~ Abundance_scaled +
    (1 | Study_id/Study_Network_id),
  data = data_combined_all
)

summary(model)
#Floral abundance was positively associated with pollinator richness capture (β = 3.03, SE = 0.12, t = 25.53, P < 0.001).

fixed_pred <- predict(model, re.form = NA)

# 总方差
var_fixed <- var(fixed_pred, na.rm = TRUE)

# 随机效应方差
var_rand <- sum(as.data.frame(VarCorr(model))$vcov)

# 残差方差
var_resid <- attr(VarCorr(model), "sc")^2

# Marginal R2（固定效应）
r2_marginal <- var_fixed / (var_fixed + var_rand + var_resid)

# Conditional R2（固定 + 随机）
r2_conditional <- (var_fixed + var_rand) / (var_fixed + var_rand + var_resid)

label <- paste0(
  "Marginal R² = ", round(r2_marginal, 2),
  "\nConditional R² = ", round(r2_conditional, 2)
)


fig1_relation <- ggplot(data_combined_all,
                        aes(x = Abundance_scaled,
                            y = Proportion_of_richness)) +
  
  geom_point(alpha = 0.5, size = 1.8, color = "#9E9AC8") +
  
  geom_smooth(method = "lm",
              color = "#6A51A3",
              se = TRUE,
              linewidth = 1) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom"
  )+
  
  
  labs(
    x = "Standardized log Flower Abundance",
    y = "Percentage of Pollinator richness captured"
  ) +
  
  annotate("text",
           x = Inf, y = Inf,
           label = label,
           hjust = 1.1, vjust = 1.5,
           size = 4)

data_combined_all$Predicted_capture <- predict(
  model,
  newdata = data_combined_all,
  re.form = NA,
  allow.new.levels = TRUE
)

data_combined_all$Residual_capture <-
  data_combined_all$Proportion_of_richness -
  data_combined_all$Predicted_capture

fig2_relation <- ggplot(data_combined_all,
                        aes(x = Abundance_scaled,
                            y = Residual_capture)) +
  
  geom_point(aes(color = performance_class),
             alpha = 0.7, size = 1.8) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  scale_color_manual(values = c(
    "Rare and attractive" = "#D7301F",
    "Above expected" = "#FC8D59",
    "Below expected" = "grey"
  )) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom"
  )+
  
  
  labs(
    x = "Standardized log Flower Abundance",
    y = "Deviation from expected richness capture (%)"
  )




##----------------------
##-- Sensitive test 
##----------------------
run_present_absent_test <- function(data, rare_q, attr_q) {
  
  # cutoff
  rare_thr <- quantile(data$Abundance_scaled, rare_q, na.rm = TRUE)
  attr_thr <- quantile(data$Residual_capture, attr_q, na.rm = TRUE)
  
  # rare-attractive species
  rare_attr <- data %>%
    filter(
      Abundance_scaled < rare_thr,
      Residual_capture > attr_thr
    )
  
  rare_species <- unique(rare_attr$Plant_species)
  
  # 每个网络是否包含这些 species
  network_highattr <- data_interact %>%
    distinct(Study_Network_id, Plant_accepted_name) %>%
    mutate(is_highattr = Plant_accepted_name %in% rare_species)
  
  top10_tbl <- top_10_species %>%
    rename(Plant_accepted_name = Plant_species) %>%
    group_by(Study_Network_id) %>%
    summarise(top10_list = list(unique(Plant_accepted_name)), .groups = "drop")
  
  network_group <- network_highattr %>%
    left_join(top10_tbl, by = "Study_Network_id") %>%
    group_split(Study_Network_id) %>%
    purrr::map_dfr(function(df){
      
      id <- unique(df$Study_Network_id)
      
      highattr_species <- df$Plant_accepted_name[df$is_highattr]
      top10 <- unique(unlist(df$top10_list[1]))
      
      tibble(
        Study_Network_id = id,
        present = any(df$is_highattr, na.rm = TRUE),
        fully_covered = all(highattr_species %in% top10)
      )
    }) %>%
    mutate(
      group = ifelse(present & !fully_covered, "Present", "Absent")
    )
  
  # merge response
  df_test <- percent_10 %>%
    left_join(network_group %>%
                select(Study_Network_id, group),
              by = "Study_Network_id")
  
  # Wilcoxon
  wt <- wilcox.test(percentage_Abun10 ~ group, data = df_test)
  
  # effect size (Cliff’s delta approx via median diff)
  eff <- with(df_test,
              median(percentage_Abun10[group=="Present"], na.rm=TRUE) -
                median(percentage_Abun10[group=="Absent"], na.rm=TRUE))
  
  tibble(
    rare_q = rare_q,
    attr_q = attr_q,
    p_value = wt$p.value,
    effect_size = eff,
    n_species = length(rare_species)
  )
}

rare_cut <- c(0.05, 0.10, 0.15)
attr_cut <- c(0.90, 0.85, 0.80)

results_list <- list()

for (r in rare_cut) {
  for (a in attr_cut) {
    
    results_list[[paste(r, a, sep="_")]] <- 
      run_present_absent_test(data_combined_all, r, a)
  }
}

results <- bind_rows(results_list)

sens_test <- ggplot(results, aes(x = attr_q,
                                 y = rare_q)) +
  
  # 1️⃣ effect size surface
  geom_tile(aes(fill = effect_size)) +
  
  scale_fill_gradient2(
    low = "#D7301F",
    mid = "white",
    high = "#1A9850",
    name = "Effect size (slope)"
  ) +
  
  # 2️⃣ significance border
  geom_tile(aes(color = p_value < 0.05),
            fill = NA,
            linewidth = 0.8) +
  
  scale_color_manual(values = c(
    "TRUE" = "black",
    "FALSE" = NA
  ), guide = "none") +
  
  # 3️⃣ original grid points
  geom_point(color = "white", size = 2) +
  
  labs(
    x = "Attractiveness threshold (quantile)",
    y = "Rarity threshold (quantile)"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "bottom"
  )


library(cowplot)

blank <- ggplot() + theme_void()

# 统一风格（很重要，避免大小/边距不一致）
theme_all <- theme(
  plot.margin = margin(5, 5, 5, 5),
  axis.title = element_text(size = 11),
  axis.text = element_text(size = 10)
)

fig1_relation <- fig1_relation + theme_all
fig2_relation <- fig2_relation + theme_all
sens_test <- sens_test + theme_all

# 1️⃣ 第一行
col1 <- plot_grid(
  fig1_relation,
  fig2_relation,
  #labels = c("A", "B"),
  label_size = 14,
  ncol = 1,
  rel_heights = c(1, 1),
  align = "hv"
)

# 2️⃣ 第二行
col2 <- plot_grid(
  sens_test,
  blank,
  #labels = c("C", ""),
  label_size = 14,
  ncol = 1,
  rel_heights = c(1, 1),
  align = "hv"
)

# 3️⃣ 总图
final_fig <- plot_grid(
  col1,
  col2,
  ncol = 2,
  rel_widths = c(1.2, 1)
)

final_fig

ggsave(
  "./result_260526/rare_attrac_sens_test.png",
  final_fig,
  width = 10,
  height = 9,
  dpi = 300
)
##########################################

###  Rare favored plant list


###########################################################################

#------Rare attractive plant list Word output

library(dplyr)
# rarelist: 高吸引力植物
rarelist <- Rare_favour_plants%>%distinct(Plant_accepted_name)
#%>% filter(p.value > 0.05)

library(dplyr)
library(flextable)
library(officer)
library(stringr)

sp_list_rare <- rarelist %>%
  left_join(traits, by = "Plant_accepted_name") %>%
  distinct(Plant_accepted_name, flw_shape_revised)

# sp_list_rare$flw_shape_revised <- stringr::str_wrap(
#   sp_list_rare$flw_shape_revised,
#   width = 40
# )

# # 创建 flextable
# ft <- flextable(sp_list_rare) %>%
#   set_header_labels(
#     Plant_accepted_name = "Plant species",
#     flw_shape_revised = "Flower shape"
#   ) %>%
#   theme_booktabs() %>%
#   fontsize(size = 11, part = "all") %>%
#   italic(j = "Plant_accepted_name", part = "body") %>%
#   width(j = "Plant_accepted_name", width = 3) %>%
#   width(j = "flw_shape_revised", width = 5) %>%
#   align(j = c("Plant_accepted_name", "flw_shape_revised"), align = "left")

# ⭐按花型分组，把物种名拼接成一行（斜体在表格里再处理）
sp_list_grouped <- sp_list_rare %>%
  filter(!is.na(flw_shape_revised)) %>%
  group_by(flw_shape_revised) %>%
  summarise(
    Plant_species = paste(sort(unique(Plant_accepted_name)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(flw_shape_revised) %>%
  mutate(Illustration = "") %>%               # ⭐留空列，方便后续手动插图
  select(Illustration, flw_shape_revised, Plant_species)

# 换行处理（如果花型名称或物种列表过长）
sp_list_grouped$Plant_species <- stringr::str_wrap(sp_list_grouped$Plant_species, width = 60)

# 创建 flextable
ft <- flextable(sp_list_grouped) %>%
  set_header_labels(
    Illustration = "",
    flw_shape_revised = "Flower shape",
    Plant_species = "Plant species"
  ) %>%
  theme_booktabs() %>%
  fontsize(size = 11, part = "all") %>%
  italic(j = "Plant_species", part = "body") %>%
  width(j = "Illustration", width = 1) %>%
  width(j = "flw_shape_revised", width = 2.5) %>%
  width(j = "Plant_species", width = 4.5) %>%
  align(j = c("Illustration", "flw_shape_revised", "Plant_species"), align = "left")

# 导出 Word
# 导出 Word
doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = "./result_260526/sp_list_rare_new.docx")



#  每个网络哪些植物是高吸引力
network_highattr <- data_interact %>%
  distinct(Study_Network_id, Plant_accepted_name) %>%
  mutate(is_highattr = Plant_accepted_name %in% rarelist$Plant_accepted_name)

# 每个网络 top 10 植物
top10_tbl <- top_10_species %>%
  rename(Plant_accepted_name = Plant_species) %>%
  group_by(Study_Network_id) %>%
  summarise(top10_list = list(unique(Plant_accepted_name)), .groups = "drop")


# 每个网络判断 Present / Absent
network_group <- network_highattr %>%
  left_join(top10_tbl, by = "Study_Network_id") %>%
  group_split(Study_Network_id) %>%
  map_dfr(function(df){
    
    id <- unique(df$Study_Network_id)
    
    highattr_species <- df$Plant_accepted_name[df$is_highattr]
    top10 <- unique(unlist(df$top10_list[1]))
    
    tibble(
      Study_Network_id = id,
      highattr_in_network = any(df$is_highattr, na.rm = TRUE),
      highattr_all_in_top10 = all(highattr_species %in% top10)
    )
  }) %>%
  mutate(
    present_group = case_when(
      highattr_in_network & !highattr_all_in_top10 ~ "Present",
      TRUE ~ "Absent"
    )
  )

network_group

# network_group 有列: Study_Network_id, group (Present/Absent)
attract_percent_10 <- percent_10 %>%
  left_join(network_group %>% dplyr::select(Study_Network_id, present_group), by = "Study_Network_id")

shapiro.test(attract_percent_10$percentage_Abun10[attract_percent_10$present_group=="Present"])
shapiro.test(attract_percent_10$percentage_Abun10[attract_percent_10$present_group=="Absent"])

wilcox.test(percentage_Abun10 ~ present_group, data = attract_percent_10)


library(ggplot2)
library(dplyr)
library(stringr)
library(ggsignif)

# 1️⃣ 分组
attract_percent_10 <- attract_percent_10 %>%
  mutate(
    present_group = str_squish(present_group) %>% str_to_title(),
    present_group = factor(present_group, levels = c("Present", "Absent"))
  )

# 2️⃣ Wilcoxon test
wilcox_res <- wilcox.test(percentage_Abun10 ~ present_group, data = attract_percent_10)

sig_label <- if(wilcox_res$p.value < 0.001) {
  "***"
} else if(wilcox_res$p.value < 0.01) {
  "**"
} else if(wilcox_res$p.value < 0.05) {
  "*"
} else {
  "ns"
}

# 3️⃣ summary（n）
data_summary <- attract_percent_10 %>%
  group_by(present_group) %>%
  summarise(
    Mean_Percentage = mean(percentage_Abun10, na.rm = TRUE),
    SD_Percentage = sd(percentage_Abun10, na.rm = TRUE),
    n_non_missing = sum(!is.na(percentage_Abun10)),
    .groups = "drop"
  )

group_colors <- c("Present" = "#E41A1C", "Absent" = "#377EB8")

# 4️⃣ y 位置
y_max <- max(attract_percent_10$percentage_Abun10, na.rm = TRUE)
y_pos <- y_max * 1.1

# =========================
# 5️⃣ FINAL PLOT
# =========================
plot <- ggplot(attract_percent_10,
               aes(x = present_group,
                   y = percentage_Abun10,
                   fill = present_group)) +
  
  geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.15, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.3, alpha = 0.4) +
  
  # mean
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = sprintf("%.1f%%", stat(y))),
    vjust = -1,
    position = position_nudge(x = -0.4)
  ) +
  
  # n
  geom_text(
    data = data_summary,
    aes(x = present_group,
        y = Mean_Percentage + SD_Percentage,
        label = paste0("n = ", n_non_missing)),
    vjust = -0.5,
    size = 4,
    position = position_nudge(x = -0.4)
  ) +
  
  
  #号显著性（重点！）
  geom_signif(
    comparisons = list(c("Present", "Absent")),
    annotations = sig_label,
    y_position = y_pos,
    tip_length = 0.01,
    textsize = 6
  ) +
  
  scale_fill_manual(values = group_colors) +
  
  labs(
    x = "Presence of highly attractive plants in full network",
    y = "Percent of pollinator richness captured (%)"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.margin = margin(5, 15, 5, 5)
  )

plot

ggsave("./result_260526/rare_atract.png", plot, width = 6.5, height = 4, units = "in", dpi = 300)

#########################################################################
# 
library(dplyr)
library(flextable)
library(officer)

# 1️⃣ 汇总 + 直接生成论文格式列
table2 <- attract_percent_10 %>%
  group_by(present_group) %>%
  summarise(
    Mean = mean(percentage_Abun10, na.rm = TRUE),
    SD = sd(percentage_Abun10, na.rm = TRUE),
    n = sum(!is.na(percentage_Abun10)),
    .groups = "drop"
  ) %>%
  
  # ⭐ 合并 Mean ± SD（论文标准写法）
  mutate(
    `Richness captured (%)` = sprintf("%.2f ± %.2f", Mean, SD),
    `P-value` = signif(wilcox_res$p.value, 3)
  ) %>%
  
  # ⭐ 改列名（关键）
  dplyr::select(
    Group = present_group,
    `Richness captured (%)`,
    n,
    `P-value`
  )

# 2️⃣ flextable 美化
ft_table2 <- flextable(table2) %>%
  theme_booktabs() %>%
  autofit() %>%
  fontsize(size = 10, part = "all") %>%
  align(align = "center", part = "all") %>%
  
  # ⭐ 表头加粗
  bold(part = "header") %>%
  
  # ⭐ Group 左对齐（更像论文）
  align(j = "Group", align = "left", part = "all")

# 3️⃣ Word 输出（优化标题）
doc <- read_docx() %>%
  body_add_par(
    "Table 2. Comparison of pollinator richness captured by the top 10 abundant plants between networks with and without highly attractive plant species (Wilcoxon rank-sum test).",
    style = "heading 2"
  ) %>%
  body_add_flextable(ft_table2)

# 4️⃣ 保存
print(doc, target = "./result_260526/Table2_Wilcoxon.docx")

##########################################

# Supp abundance-richness in each network
# effect size forest plot
library(dplyr)
library(purrr)
library(broom)
library(ggplot2)

effect_table <- data_combined_all %>%
  group_by(Study_Network_id) %>%
  group_split() %>%
  map_dfr(function(df){
    
    # 每个network一个独立模型 ✔
    model <- lm(Proportion_of_richness ~ Abundance_scaled, data = df)
    
    tt <- broom::tidy(model, conf.int = TRUE)
    
    tibble(
      Study_Network_id = unique(df$Study_Network_id),
      slope = tt$estimate[2],
      CI_low = tt$conf.low[2],
      CI_high = tt$conf.high[2],
      p = tt$p.value[2],
      n = nrow(df)
    )
  })


library(ggplot2)

forest<-ggplot(effect_table,
               aes(x = slope,
                   y = reorder(Study_Network_id, slope))) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
                 height = 0.2,
                 color = "grey40") +
  
  geom_point(aes(size = n),
             color = "#2C7FB8",
             alpha = 0.8) +
  
  theme_classic(base_size = 3) +
  
  labs(
    x = "Effect size: Abundance → richness capture",
    y = "Network"
  ) +
  
  theme(legend.position = "none")

ggsave("./result_260526/forest.png", forest, width = 9, height = 15, units = "in", dpi = 300)


##-----------(4) hypothesis 4 : Different flower shapes attract different pollinators groups 
# 花型对传粉者组成的分析（268 sites not pooled）genus级别

##########################################
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(ggpattern)


all_data <- data_merge %>%
  filter(!is.na(flw_shape_revised),
         !is.na(Pollinator_genus))%>%
  filter(! flw_shape_revised %in% c("trap flowers","brush flowers"))


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

#####--------------------（选项2开始）
# 4.1.2 visit frequency
# #  Plant × Site × Pollinator Matrix (Interaction Frequency considered)
# all_data_unique <- all_data %>%
#   group_by(Plant_accepted_name, Study_Network_id, Pollinator_genus) %>%
#   summarise(Visit_frequency = sum(Interaction_addup, na.rm = TRUE), .groups = "drop")
# 
# # 构建矩阵
# mat <- all_data_unique %>%
#   pivot_wider(names_from = Pollinator_genus,
#               values_from = Visit_frequency,
#               values_fill = 0)


##########------------------（2结束）
# 接2

plant_flw <- all_data_unique %>%
  dplyr::select(Plant_accepted_name, flw_shape_revised, Study_Network_id,Study_id) %>%
  distinct()

mat_df <- as.data.frame(mat)

mat_df$SampleID <- paste(mat_df$Plant_accepted_name,
                         mat_df$Study_Network_id,
                         sep = "_")

rownames(mat_df) <- mat_df$SampleID

mat_numeric <- mat_df[, !(names(mat_df) %in%
                            c("Plant_accepted_name",
                              "Study_Network_id", 
                              "SampleID"))]

meta <- plant_flw
table(meta$flw_shape_revised)

meta$SampleID <- paste(meta$Plant_accepted_name,
                       meta$Study_Network_id,
                       sep = "_")

meta2 <- meta %>%
  filter(!flw_shape_revised %in% c("trap flowers","brush flowers"))

mat_numeric2 <- mat_numeric[meta2$SampleID, ]


# ---------------------------
# 4.2   NMDS analysis
# ---------------------------
# set.seed(123)
# 
# # presence/absence 数据可用 Jaccard 或 Bray-Curtis
d <- vegdist(mat_numeric2, method = "jaccard")  # frequency用method = "bray"
# 
# ord_nmds <- metaMDS(mat_numeric2,
#                     distance = "jaccard", #0/1数据用Jaccard,直接使用 interaction frequency则用 Bray-Curtis 会更合适
#                     k = 2,
#                     trymax = 50,
#                     autotransform = FALSE)
# 
# saveRDS(ord_nmds,"ord_nmds260527")
# ord_nmds<-readRDS("ord_nmds260527")
# 
# scores_df <- as.data.frame(scores(ord_nmds)$sites)
# scores_df$SampleID <- rownames(scores_df)
# 
# scores_df <- left_join(scores_df, meta2, by = "SampleID")
# colnames(scores_df)[1:2] <- c("NMDS1", "NMDS2")


# ---------------------------
# 4.2.1 PERMANOVA
# ---------------------------
# adonis_res <- adonis2(d ~ flw_shape_revised,
#                       data = meta2,
#                       permutations = 999,
#                       strata = meta2$Study_Network_id)

adonis_res_study_ctrl <- adonis2(d ~ Study_id +flw_shape_revised,
                                 data = meta2,
                                 permutations = 999,
                                 strata = meta2$Study_Network_id,
                                 by = "margin")
# 是否
head(meta2)
print(adonis_res)
print(adonis_res_study_ctrl)
saveRDS(adonis_res,"data/processed/adonis_res.rds")

adonis_res<-readRDS("data/processed/adonis_res.rds")
#PERMANOVA结果显著，但R方极小
#这说明globally, 花型对传粉者群落组成的影响显著，但不是主要影响因素
#但如果单独对每一个网络进行分析，结果可能会不同

#############

# 组间差异主要由大部分还是少部分网络驱动

###########
library(dplyr)
library(vegan)

networks <- unique(meta2$Study_Network_id)

res_list <- list()

skip_list <- data.frame(
  Network = character(),
  Reason = character()
)


for(net in networks){
  
  sub_meta <- meta2 %>%
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
nrow(network_perm)
head(network_perm)

network_perm %>%
  summarise(
    total_network=n(),
    significant=sum(P < 0.05),
    proportion=mean(P < 0.05)
  )

table(skip_list$Reason)

network_perm %>%
  summarise(
    total_network = n(),
    significant = sum(P < 0.05),
    proportion = mean(P < 0.05)
  )

## 检查R2是否受到sampling effeort 的影响
library(ggplot2)


ggplot(network_perm,
       aes(x=N_samples,
           y=R2))+
  geom_point()+
  geom_smooth(method="lm",
              se=TRUE)+
  theme_bw()+
  labs(
    x="Number of plant samples",
    y=expression(PERMANOVA~R^2)
  )

summary(
  lm(R2 ~ N_samples,
     data=network_perm)
)

summary(
  lm(R2 ~ N_samples + N_shapes,
     data=network_perm)
)

library(car)

vif(
  lm(R2 ~ N_samples + N_shapes,
     data=network_perm)
)

############33
model <- lm(
  R2 ~ N_samples + N_shapes,
  data=network_perm
)

network_perm$R2_adjusted <- residuals(model)
summary(network_perm$R2_adjusted)

summary(network_perm$R2)

ggplot(network_perm,
       aes(R2_adjusted))+
  geom_histogram()

#########哪些网络的花型影响最大
strong_network <- network_perm %>%
  filter(
    P < 0.05,
    R2_adjusted > 0.1
  ) %>%
  arrange(desc(R2_adjusted))


strong_network

strongdata<-plottest%>%filter(Study_Network_id %in% strong_network$Network)%>%select(Study_Network_id,percentage_Abun5,FlwShape_Top5)

t.test(strongdata$percentage_Abun5,
       strongdata$FlwShape_Top5,
       paired = TRUE)
wilcox.test(strongdata$percentage_Abun5,
            strongdata$FlwShape_Top5,
            paired = TRUE)


#
#--------
# 4.2.2 dispersion difference test
# If：P > 0.05，PERMANOVA result reliable
#-------
disp <- betadisper(d, meta2$flw_shape_revised)
anova(disp) #P<0.05

#PERMANOVA 显著 + betadisper 显著 ⇒ 差异可能由 dispersion 驱动，而不是 centroid difference
#因此，无法证明花型显著性影响传粉者群落组成

R2_strata <- round(adonis_res_study_ctrl$R2[1], 3)
P_strata <- adonis_res_study_ctrl$`Pr(>F)`[1]
P_strata_text <- ifelse(P_strata < 0.001, "P < 0.001", paste0("P = ", round(P_strata, 3)))

############################################
# Figure A
# Variation in pollinator community composition among flower shapes
# betadisper + Tukey HSD
############################################


library(dplyr)
library(ggplot2)
library(betapart)
library(vegan)
library(multcompView)
library(patchwork)

############################################
# A. Betadisper + Tukey
############################################

# dataframe

disp_df <- data.frame(
  FlowerShape = meta2$flw_shape_revised,
  Distance = disp$distances
)



# Tukey test

disp_aov <- aov(
  Distance ~ FlowerShape,
  data = disp_df
)


tukey_res <- TukeyHSD(disp_aov)


tukey_letters <- multcompLetters4(
  disp_aov,
  tukey_res
)


letters_df <- data.frame(
  FlowerShape = names(tukey_letters[[1]]$Letters),
  Letters = tukey_letters[[1]]$Letters
)



############################################
# 统一flower shape排序
# 按 median dispersion
############################################


shape_order <- disp_df %>%
  group_by(FlowerShape) %>%
  summarise(
    median_disp = median(Distance)
  ) %>%
  arrange(median_disp) %>%
  pull(FlowerShape)



disp_df$FlowerShape <- factor(
  disp_df$FlowerShape,
  levels = shape_order
)



letters_df$FlowerShape <- factor(
  letters_df$FlowerShape,
  levels = shape_order
)



############################################
# Panel A
############################################


label_pos <- disp_df %>%
  group_by(FlowerShape) %>%
  summarise(
    xpos=max(Distance)
  ) %>%
  left_join(
    letters_df,
    by="FlowerShape"
  )


label_pos$xpos <- label_pos$xpos*1.08



pA <- ggplot(
  disp_df,
  aes(
    x=FlowerShape,
    y=Distance
  )
)+
  
  geom_boxplot(
    fill="#A23B72",
    alpha=0.65,
    width=0.65,
    outlier.shape=NA
  )+
  
  geom_jitter(
    width=0.15,
    size=1,
    alpha=0.25
  )+
  
  geom_text(
    data=label_pos,
    aes(
      x=FlowerShape,
      y=xpos,
      label=Letters
    ),
    size=5,
    fontface="bold"
  )+
  
  coord_flip()+
  
  labs(
    x=NULL,
    y="Distance to centroid"
  )+
  
  theme_classic(base_size=12) +
  theme(
    axis.title    = element_text(face = "bold", size = 15),
    axis.text     = element_text(color = "black", size = 12),
  )



############################################
# B. beta diversity partition
############################################


comm_pa <- decostand(
  mat_numeric2,
  "pa"
)


bc <- beta.pair(
  comm_pa,
  index.family="sorensen"
)


turnover_mat <- as.matrix(bc$beta.sim)

nestedness_mat <- as.matrix(bc$beta.sne)



flower_shape_vector <- meta2$flw_shape_revised



get_shape_mean <- function(mat, group_vec){
  
  sapply(shape_order,function(shp){
    
    sites <- which(group_vec==shp)
    
    mean(
      mat[sites,sites][upper.tri(mat[sites,sites])],
      na.rm=TRUE
    )
    
  })
}



df_beta <- data.frame(
  
  FlowerShape = shape_order,
  
  turnover =
    get_shape_mean(
      turnover_mat,
      flower_shape_vector
    ),
  
  nestedness =
    get_shape_mean(
      nestedness_mat,
      flower_shape_vector
    )
  
)

df_beta$FlowerShape <- factor(
  df_beta$FlowerShape,
  levels = shape_order
)



############################################
# Panel B1 turnover
############################################
pB <- ggplot(
  df_beta,
  aes(
    x=FlowerShape,
    y=turnover
  )
)+
  
  geom_point(
    size=3,
    color="tomato"
  )+
  
  coord_flip()+
  
  labs(
    x=NULL,
    y="Turnover (βsim)"
  )+
  scale_y_continuous(
    expand = expansion(mult=c(0.1,0.1))
  )+
  
  theme_classic(base_size=12)+
  
  theme(
    axis.title    = element_text(face = "bold", size = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank()
  )  

############################################
# Panel B2 nestedness
############################################
pC <- ggplot(
  df_beta,
  aes(
    x=FlowerShape,
    y=nestedness
  )
)+
  
  geom_point(
    size=3,
    color="steelblue"
  )+
  
  coord_flip()+
  
  labs(
    x=NULL,
    y="Nestedness (βsne)"
  )+
  scale_y_continuous(
    expand = expansion(mult=c(0.1,0.1))
  )+
  
  theme_classic(base_size=12)+
  
  theme(
    axis.title    = element_text(face = "bold", size = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank()
  )  
############################################
# Combine
############################################


figure_beta1 <- cowplot::plot_grid(
  pB,
  pC,
  ncol = 2,
  rel_widths = c(0.8,0.8),
  align = "h"
)

figure_beta1

pA

ggsave(
  "result_260526/flower_shape_beta_partitionA.png",
  pA,
  width=6.5,
  height=5,
  dpi=300
)

ggsave(
  "result_260526/flower_shape_beta_partitionBC.png",
  figure_beta1,
  width=3.7,
  height=5,
  dpi=300
)


# disp_aov <- aov(
#   disp$distances ~ meta2$flw_shape_revised
# )
# 
# tukey_res <- TukeyHSD(disp_aov)
# tukey_res
# library(multcompView)
# tukey_letters <- multcompLetters4(
#   disp_aov,
#   tukey_res
# )
# 
# tukey_letters
# 
# letters_df <- data.frame(
#   FlowerShape = names(tukey_letters[[1]]$Letters),
#   Letters = tukey_letters[[1]]$Letters
# )
# 
# letters_df
# 
# disp_df <- data.frame(
#   FlowerShape = meta2$flw_shape_revised,
#   Distance = disp$distances
# )
# 
# label_pos <- disp_df %>%
#   group_by(FlowerShape) %>%
#   summarise(
#     y = max(Distance, na.rm = TRUE)
#   ) %>%
#   left_join(
#     letters_df,
#     by="FlowerShape"
#   )
# 
# label_pos$y <- label_pos$y * 1.1
# 
# 
# 
# p_disp <- ggplot(
#   disp_df,
#   aes(
#     x = FlowerShape,
#     y = Distance
#   )
# )+
#   
#   geom_boxplot(
#     width=0.6,
#     fill="#A23B72",
#     alpha=0.6,
#     outlier.shape=NA
#   )+
#   
#   geom_jitter(
#     width=0.15,
#     size=1.5,
#     alpha=0.35
#   )+
#   
#   geom_text(
#     data=label_pos,
#     aes(
#       x=FlowerShape,
#       y=y,
#       label=Letters
#     ),
#     size=5,
#     fontface="bold"
#   )+
#   
#   labs(
#     x="Flower shape",
#     y="Distance to centroid"
#   )+
#   
#   theme_classic(
#     base_size=14
#   )+
#   
#   theme(
#     axis.text.x = element_text(
#       angle=45,
#       hjust=1
#     )
#   )
# 
# 
# p_disp
#结合以上结果，无法证明花型对传粉者群落产生显著性影响。

##########################

# Which flower shapes drive dispersion? (组内差异由哪个花型带来)

#################################

library(dplyr)

disp_df <- data.frame(
  shape = meta2$flw_shape_revised,
  dist = disp$distances
)

shape_var <- disp_df %>%
  group_by(shape) %>%
  summarise(
    mean_dispersion = mean(dist),
    sd_dispersion = sd(dist),
    n = n()
  )

shape_var

#mean_dispersion = 这个花型的“传粉者组成不稳定程度”

#越大 = 越“乱 / heterogenous”
#越小 = 越“稳定 / consistent”

head(meta2)

shape_richness <- meta2 %>%
  group_by(flw_shape_revised) %>%
  summarise(
    plant_species = n_distinct(Plant_accepted_name)
  )

df <- left_join(shape_var, shape_richness,
                by=c("shape"="flw_shape_revised"))

# saveRDS(df, "data/processed/mean_dis_plot.rds")
df<-readRDS("data/processed/mean_dis_plot.rds")

library(ggrepel)  # 更好地处理标签位置

# 先计算相关系数
corr_result <- cor.test(df$mean_dispersion,
                        df$plant_species,
                        method = "spearman")
line_type <- ifelse(corr_result$p.value < 0.05,
                    "solid",
                    "dashed")

disper <- ggplot(df, aes(x = plant_species,
                         y = mean_dispersion)) +
  
  geom_point(size = 4.5, 
             color = "#2E86AB",
             alpha = 0.75) +
  
  geom_smooth(method = "lm", 
              se = TRUE,
              color = "#A23B72",
              fill = "#A23B72",
              alpha = 0.1,
              linewidth = 0.9,
              linetype = line_type) +
  
  ggrepel::geom_text_repel(
    aes(label = shape),
    size = 3.5,
    fontface = "italic",
    segment.color = "gray50",
    segment.size = 0.3,
    segment.alpha = 0.5,
    max.overlaps = 15
  ) +
  
  # 添加相关系数信息
  annotate("text", 
           x = Inf, y = -Inf,
           label = paste0("ρ = ", round(corr_result$estimate, 3), 
                          "\np = ", round(corr_result$p.value, 4)),
           hjust = 1.1, vjust = -0.5,
           size = 3.8,
           fontface = "italic",
           color = "#333333") +
  
  labs(
    x = "Plant species richness within flower shape category",
    y = "Mean dispersion of pollinator community composition"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black", size = 0.6),
    plot.margin = margin(10, 10, 10, 10)
  )

disper

ggsave("/Chap1_TargetPlant_to_monitor/result_260526/shape disper.png", disper, width = 5, 
       height = 5, units = "in", dpi = 300)  

######计算turnover和嵌套
# 
# library(betapart)
# 
# comm_pa <- decostand(mat_numeric2, "pa")
# 
# bc <- beta.pair(comm_pa, index.family = "sorensen")
# str(bc)
# turnover_mat <- as.matrix(bc$beta.sim)
# nestedness_mat <- as.matrix(bc$beta.sne)
# total_mat <- as.matrix(bc$beta.sor)
# 
# flower_shape_vector <- meta2$flw_shape_revised
# 
# library(dplyr)
# 
# turnover_df <- data.frame(
#   shape = flower_shape_vector
# )
# 
# 
# 
# get_shape_mean <- function(mat, group_vec) {
#   sapply(unique(group_vec), function(shp) {
#     sites <- which(group_vec == shp)
#     mean(mat[sites, sites][upper.tri(mat[sites, sites])], na.rm = TRUE)
#   })
# }
# 
# shape_turnover <- get_shape_mean(turnover_mat, flower_shape_vector)
# shape_nestedness <- get_shape_mean(nestedness_mat, flower_shape_vector)
# 
# df_beta <- data.frame(
#   shape = unique(flower_shape_vector),
#   turnover = shape_turnover,
#   nestedness = shape_nestedness
# )
# df_beta[order(-df_beta$turnover), ]
# df_all <- left_join(shape_var, df_beta,
#                     by = "shape")
# df_all$shape <- factor(df_all$shape,
#                        levels = df_all$shape[order(df_all$turnover)])
# 
# library(ggplot2)
# library(patchwork)
# 
# p1 <- ggplot(df_all, aes(x = shape, y = mean_dispersion)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(
#     ymin = mean_dispersion - sd_dispersion,
#     ymax = mean_dispersion + sd_dispersion
#   ), width = 0.2) +
#   coord_flip() +
#   theme_classic(base_size = 13)
# 
# p2 <- ggplot(df_all, aes(x = shape, y = turnover)) +
#   geom_point(size = 3, color = "tomato") +
#   coord_flip() +
#   labs(y = "Turnover (βsim)", x = NULL) +
#   theme_classic(base_size = 13)
# 
# p3 <- ggplot(df_all, aes(x = shape, y = nestedness)) +
#   geom_point(size = 3, color = "steelblue") +
#   coord_flip() +
#   labs(y = "Nestedness (βsne)", x = NULL) +
#   theme_classic(base_size = 13)
# 
# library(gridExtra)
# 
# plot_turn_over<-grid.arrange(p1, p2, p3, ncol = 3)
# 
# ggsave("/Chap1_TargetPlant_to_monitor/result_260526/plot_turn_over.png", plot_turn_over, width = 12, 
#        height = 5, units = "in", dpi = 300)  

###############
library(dplyr)
sp_list_rare_info<-sp_list_rare%>%
  left_join(data_interact%>%select(Plant_accepted_name,Plant_order,Plant_family,Plant_genus)%>%
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

##----------- PERMANOVA result output--------------

library(dplyr)
library(flextable)
library(officer)

perm_table <- data.frame(
  Test = c("PERMANOVA", "PERMANOVA"),
  Factor = c("Flower shape", "Residual"),
  Df = c(adonis_res$Df[1], adonis_res$Df[2]),
  F = c(round(adonis_res$F[1], 2), NA),
  R2 = c(round(adonis_res$R2[1], 3), round(adonis_res$R2[2], 3)),
  P_value = c(as.character(
    ifelse(adonis_res$`Pr(>F)`[1] < 0.001, "<0.001",
           round(adonis_res$`Pr(>F)`[1], 3))
  ),
  NA)
  
)

disp_aov <- anova(disp)

disp_table <- data.frame(
  Test = "Betadisper",
  Factor = rownames(disp_aov),
  Df = disp_aov$Df,
  F = round(disp_aov$`F value`, 2),
  R2 = NA,
  P_value = as.character(
    ifelse(disp_aov$`Pr(>F)` < 0.001, "<0.001",
           round(disp_aov$`Pr(>F)`, 3))
  )
  
)

final_table <- bind_rows(perm_table, disp_table)


ft <- flextable(final_table) %>%
  theme_booktabs() %>%
  autofit() %>%
  fontsize(size = 10, part = "all") %>%
  
  flextable::font(fontname = "Arial", part = "all") %>%   # ⭐关键修复
  
  bold(part = "header") %>%
  
  align(j = c("Test", "Factor"), align = "left", part = "all") %>%
  align(j = c("Df", "F", "R2", "P_value"), align = "center", part = "all") %>%
  
  set_header_labels(
    Test = "Test",
    Factor = "Factor",
    Df = "Df",
    F = "F",
    R2 = "R²",
    P_value = "P-value"
  )


# =========================
#  Word
# =========================
doc <- read_docx()

doc <- body_add_par(
  doc,
  "Table X. PERMANOVA and multivariate dispersion (betadisper) results for pollinator community composition across flower shapes.",
  style = "heading 2"
)

doc <- body_add_flextable(doc, value = ft)

print(doc, target = "./result_260526/community_shape_multivariate_test.docx")

#-------NMDS plot
ggplot(scores_df, aes(x = NMDS1, y = NMDS2, color = flw_shape_revised)) +
  geom_point(alpha = 0.6) +
  stat_ellipse(aes(group = flw_shape_revised), linetype = 2) +
  theme_classic()




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

saveRDS(data_merge,"data_merge0713.rds")
a<-data_functional%>%filter(Pollinator_order == "Coleoptera")

# 检查 Other 里还有什么
data_functional %>%
  filter(functional_group == "Other") %>%
  distinct(Pollinator_order, Pollinator_genus) %>%
  arrange(Pollinator_order)

all_groups <- data_functional %>%
  distinct(functional_group)


network_fun <- data_functional %>%
  group_by(
    Study_Network_id,
    flw_shape_revised,
    functional_group
  ) %>%
  summarise(
    n_interactions = sum(Interaction_addup,
                         na.rm = TRUE),
    .groups="drop"
  )

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

# ⭐ Step 4: 跨网络取均值（合作者建议）
plot_df <- network_props %>%
  group_by(
    flw_shape_revised,
    functional_group
  ) %>%
  summarise(
    mean_prop=mean(prop),
    .groups="drop"
  )

# ⭐ Step 5: 花按 Bees 比例排序：
shape_order <- plot_df %>%
  filter(functional_group=="Bees") %>%
  arrange(mean_prop) %>%
  pull(flw_shape_revised)

###or : visit sum


plot_df$flw_shape_revised <- factor(plot_df$flw_shape_revised, levels = shape_order)

# ⭐ Step 6: 功能群颜色

group_cols <- c(
  "Bees"                 = "#F5E066",   # 更柔和的黄
  "Syrphidae"            = "#85CCAE",   # 更柔和的青绿
  "Non-bee Hymenoptera"  = "#F5A88A",   # 更柔和的橙
  "Coleoptera"           = "#E8A3D1",   # 更柔和的粉
  "Lepidoptera"          = "#A5B5D9",   # 更柔和的蓝
  "Non-syrphid Diptera"  = "#B8DD7F",   # 更柔和的黄绿
  "Other"                = "#AAAAAA"    # 更柔和的灰
)
# ⭐ 关键：定义堆叠顺序（从下到上）
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
# ⭐ Step 7: 作图
visit_group <- ggplot(plot_df,
                      aes(x = flw_shape_revised,
                          y = mean_prop,
                          fill = functional_group)) +
  geom_bar(stat="identity",
           position="fill") +
  scale_fill_manual(values=group_cols) +
  labs(
    x="Flower shape",
    y="Mean proportion of pollinator visits",
    fill="Pollinator group"
  ) +
  coord_flip() +
  theme_classic(base_size=12) +
  theme(
    axis.title    = element_text(face = "bold", size = 15),
    axis.text     = element_text(color = "black", size = 12),
    legend.title  = element_text(face = "bold", size = 11),
    legend.text   = element_text(size = 10),
    legend.position = "right"
  )

visit_group

ggsave("/Chap1_TargetPlant_to_monitor/result_260526/visit_larger_group_sum.png", visit_group, width = 8.5, 
       height = 4, units = "in", dpi = 300)  

######################################################################3

all_data <- data_merge %>%
  filter(!is.na(flw_shape_revised),
         !is.na(Pollinator_genus))%>%
  filter(! flw_shape_revised %in% c("brush flowers","trap flowers"))#保留"trap flowers"吗

unique(all_data$Pollinator_family)
unique(shape_genus$Pollinator_genus)


# shape_genus <- all_data %>%
#   group_by(flw_shape_revised, Pollinator_genus) %>%
#   summarise(visits = sum(Interaction_addup, na.rm = TRUE),
#             .groups = "drop") #interaction sum

# shape_genus <- data_merge %>%
#   filter(!is.na(flw_shape_revised),
#          !is.na(Pollinator_genus),
#          !is.na(Pollinator_accepted_name),
#          !is.na(Plant_accepted_name),
#          !is.na(Study_Network_id)) %>%
#   filter(!flw_shape_revised %in% c("brush flowers", "trap flowers")) %>%
#   # ⭐ 每个网络内每对 plant × pollinator 只计一次
#   distinct(Study_Network_id, flw_shape_revised,
#            Plant_accepted_name, Pollinator_accepted_name, Pollinator_genus) %>%
#   group_by(flw_shape_revised, Pollinator_genus) %>%
#   summarise(visits = n(), .groups = "drop")# unique interactions

all_data%>%distinct(flw_shape_revised)

top10_genus <- shape_genus %>%
  group_by(Pollinator_genus) %>%
  summarise(total = sum(visits)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 10) %>%
  pull(Pollinator_genus)


shape_genus2 <- shape_genus %>%
  mutate(genus_group = ifelse(Pollinator_genus %in% top10_genus,
                              Pollinator_genus,
                              "Other"))

plot_df <- shape_genus2 %>%
  group_by(flw_shape_revised, genus_group) %>%
  summarise(visits = sum(visits), .groups = "drop") %>%
  group_by(flw_shape_revised) %>%
  mutate(prop = visits / sum(visits))

shape_order <- shape_genus %>%
  group_by(flw_shape_revised) %>%
  summarise(total = sum(visits)) %>%
  arrange(desc(total)) %>%
  pull(flw_shape_revised)

plot_df$flw_shape_revised <- factor(plot_df$flw_shape_revised,
                                    levels = shape_order)

library(RColorBrewer)

cols <- c(
  brewer.pal(8, "Set2"),
  brewer.pal(8, "Set3"),
  "#999999"
)

names(cols) <- c(top10_genus, "Other")
unique(plot_df$genus_group)

genus_order <- plot_df %>%
  group_by(genus_group) %>%
  summarise(total = sum(prop), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(genus_group)

genus_order <- c("Other", setdiff(genus_order, "Other"))

plot_df$genus_group <- factor(plot_df$genus_group, levels = genus_order)

#bruch:3+1. trap: 1 exclude them
plot_df_clean<-plot_df%>%filter(!flw_shape_revised %in% c("brush flowers"))# 这个"trap flowers"保留是因为绘图统一性


shape_order_visits <- plot_df_clean %>%
  group_by(flw_shape_revised) %>%
  summarise(total_visits = sum(visits)) %>%
  arrange(total_visits) %>%  # 升序
  pull(flw_shape_revised) #参考Overvieiw中的排列顺序


# 设置因子顺序
plot_df_clean$flw_shape_revised <- factor(plot_df_clean$flw_shape_revised,
                                          levels = shape_order_final)
# 绘图
visit_group <- ggplot(plot_df_clean,
                      aes(x = flw_shape_revised,
                          y = prop,
                          fill = genus_group)) +
  geom_bar(stat = "identity", position = "fill", alpha = 0.95) +
  scale_fill_manual(values = cols) +
  labs(
    x = "Flower shape",
    y = "Relative proportion of pollinator visits",
    fill = "Pollinator genus"
  ) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )


visit_group

library(gridExtra)

combined <- grid.arrange( visit_group, shape_class,ncol = 2)
combined

ggsave("/Chap1_TargetPlant_to_monitor/result_260526/visit_group.png", visit_group, width = 8, 
       height = 4, units = "in", dpi = 300)  






###########

#Factors influencing the effectiveness of subsampling

#################################################################################

##-----------(1) hypothesis 1 : 10 plant not enough for large networks------------

# 2.绘制对传粉者交互多度的影响 （possible supplement）

# Interaction_abun <- result10 %>%
#   filter(!is.na(Interaction_addup)) %>%  # 删除 Interaction_addup 是 NA 的行
#   group_by(Study_Network_id) %>%
#   summarise(Interaction_addup = sum(Interaction_addup), .groups = "drop")
# 
# total_Interaction<-data_interact %>%
#   group_by(Study_Network_id) %>%
#   summarise(Interaction_total = sum(Interaction_addup), .groups = "drop")

# Interaction_captured<-total_Interaction%>%
#   left_join(Interaction_abun, by = "Study_Network_id") %>%
#   mutate(percent = (Interaction_addup/Interaction_total) * 100 )
# 
# Interaction_captured<-merge(Interaction_captured, plant_percent,by.x = "Study_Network_id", by.y = "Study_Network_id", all = TRUE)%>% drop_na()

#model4 <- lm(percent ~ n_all_plant, data = Interaction_captured)
#summary(model4)

#############################
# 添加分组标签，并确保列名一致
# df1 <- plotplotplot %>%
#   mutate(
#     percent = percentage_Abun10,
#     Group = "Pollinator Richness"
#   ) %>%
#   select(plant_percent, percent, Group)

# df2 <- Interaction_captured %>%
#   mutate(Group = "Pollinator Interactions") %>%
#   select(plant_percent, percent, Group)
# 
# # 合并两个数据框
# combined_df <- bind_rows(df1, df2)
# 
# # 自定义颜色
# my_colors <- c(
#   "Pollinator Richness" = "darkred",     # 红色
#   "Pollinator Interactions" = "#9E9AC8"     # 蓝色
# )

# 绘图
# combined_plot <- ggplot(combined_df, aes(x = plant_percent, y = percent, color = Group)) +
#   geom_point(shape = 1, size = 2.5, position = position_jitter(width = 0.03, height = 0.5)) +
#   geom_smooth(method = "lm", se = FALSE, size = 0.8) +
#   scale_color_manual(values = my_colors) +
#   labs(
#     x = "Plant Species Coverage by Top 10 Abundant Plants",
#     y = "% Captured",
#     color = NULL
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
#     plot.margin = margin(10, 10, 10, 10),
#     strip.background = element_rect(fill = "white", color = NA),
#     strip.text = element_text(size = 12, face = "bold"),
#     legend.position = "top",
#     legend.text = element_text(size = 11),
#     legend.title = element_blank()
#   )
# 
# print(combined_plot)#500*500
#####################################
# 数据准备
# df1 <- plotplotplot %>%
#   mutate(
#     percent = percentage_Abun10,
#     Group = "Pollinator Richness"
#   ) %>%
#   dplyr::select(n_all_plant, percent, Group)

# df2 <- Interaction_captured %>%
#   mutate(Group = "Pollinator Interactions") %>%
#   select(n_all_plant, percent, Group)
# 
# combined_df <- bind_rows(df1, df2)
# 
# my_colors <- c(
#   "Pollinator Richness" = "darkred",
#   "Pollinator Interactions" = "#9E9AC8"
# )


# combined_plot <- ggplot(combined_df, aes(x = n_all_plant, y = percent, color = Group)) +
#   geom_point(shape = 1, size = 2.5, position = position_jitter(width = 0.03, height = 0.5)) +
#   geom_smooth(method = "lm", se = FALSE, size = 0.8, fullrange = FALSE) + # 拟合线只在数据范围内
#   scale_color_manual(values = my_colors) +
#   labs(
#     x = "Number of Plant Species in the Network",
#     y = "Percent of Pollinator Richness Captured by Subsampling",
#     color = NULL
#   ) +
#   coord_cartesian(ylim = c(0, 100)) +   # 固定 y 轴范围在 0-100
#   theme_bw(base_size = 12) +
#   theme(
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
#     plot.margin = margin(10, 10, 10, 10),
#     strip.background = element_rect(fill = "white", color = NA),
#     strip.text = element_text(size = 12, face = "bold"),
#     legend.position = "top",
#     legend.text = element_text(size = 11),
#     legend.title = element_blank()
#   )

# library(ggpubr)
# set.seed(2025)
# df1_plot <- ggplot(df1, aes(x = n_all_plant, y = percent)) +
#   geom_point(shape = 1, size = 2.5, color = "#9E9AC8", 
#              position = position_jitter(width = 0.03, height = 0.5)) +
#   geom_smooth(method = "lm", se = FALSE, size = 0.8, color = "#9E9AC8") +
#   stat_cor(
#     method = "pearson",
#     label.x.npc = 0.05,
#     label.y.npc = 0.05,
#     aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
#     r.accuracy = 0.01,
#     p.accuracy = 0.001,
#     size = 5
#   ) +
#   labs(
#     x = "Number of Plant Species in the Network",
#     y = "Percent of Pollinator Richness\nCaptured by Subsampling"
#   ) +
#   scale_x_continuous(breaks = seq(floor(min(df1$n_all_plant)),
#                                   ceiling(max(df1$n_all_plant)),
#                                   by = 10))+
#   #coord_cartesian(ylim = c(0, 100)) +
#   theme_classic(base_size = 12) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
#     axis.title = element_text(face = "bold", size = 15),
#     axis.text = element_text(color = "black", size = 12),
#     legend.position = "none",
#     plot.margin = margin(5, 15, 5, 5)
#   )
# 
# print(df1_plot)#500*500
# ggsave("./result_260526/n_plant_inter_rich_remove_grass.png", 
#        df1_plot, width = 5, height = 4.5, units = "in", dpi = 300)
# 
# summary(model2)
# 




########### Test the covirence betweenpollinator richness, plant richness and network size, and the effect of them on susampling effectiveness ######
# ########## marginal effect
# library(ggeffects)
# library(ggplot2)
# 
# plant_eff <- ggpredict(
#   m2,
#   terms = "Plant_Richness"
# )
# 
# p1 <- ggplot(plant_eff,
#              aes(x=x,
#                  y=predicted)) +
#   geom_ribbon(
#     aes(ymin=conf.low,
#         ymax=conf.high),
#     alpha=0.2
#   ) +
#   geom_line(
#     color = "#66C2A5",
#     linewidth=1
#   ) +
#   labs(
#     x="Plant richness",
#     y="Pollinator richness captured (%)"
#   ) +
#   theme_classic(base_size=13)
# 
# p1
# 
# nested_eff <- ggpredict(
#   m2,
#   terms="Nestedness"
# )
# 
# 
# p2 <- ggplot(nested_eff,
#              aes(x=x,
#                  y=predicted)) +
#   geom_ribbon(
#     aes(ymin=conf.low,
#         ymax=conf.high),
#     alpha=0.2
#   ) +
#   geom_line(
#     color = "#66C2A5",
#     linewidth=1
#   ) +
#   labs(
#     x="Nestedness",
#     y="Pollinator richness captured (%)"
#   ) +
#   theme_classic(base_size=13)
# 
# p2