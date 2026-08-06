###############################################################################

##-----------(2) hypothesis 2 : Favored plant sp by pollinators is rare (plant abundance-interaction times) 

################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(broom)
library(purrr)
library(cowplot)
library(flextable)
library(officer)
library(stringr)
library(ggsignif)

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


# rarelist: 高吸引力植物
rarelist <- Rare_favour_plants%>%distinct(Plant_accepted_name)
#%>% filter(p.value > 0.05)


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



# network_group 有列: Study_Network_id, group (Present/Absent)
attract_percent_10 <- percent_10 %>%
  left_join(network_group %>% dplyr::select(Study_Network_id, present_group), by = "Study_Network_id")

shapiro.test(attract_percent_10$percentage_Abun10[attract_percent_10$present_group=="Present"])
shapiro.test(attract_percent_10$percentage_Abun10[attract_percent_10$present_group=="Absent"])

wilcox.test(percentage_Abun10 ~ present_group, data = attract_percent_10)



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
# 
# # Supp abundance-richness in each network
# # effect size forest plot
# library(dplyr)
# library(purrr)
# library(broom)
# library(ggplot2)
# 
# effect_table <- data_combined_all %>%
#   group_by(Study_Network_id) %>%
#   group_split() %>%
#   map_dfr(function(df){
#     
#     # 每个network一个独立模型 ✔
#     model <- lm(Proportion_of_richness ~ Abundance_scaled, data = df)
#     
#     tt <- broom::tidy(model, conf.int = TRUE)
#     
#     tibble(
#       Study_Network_id = unique(df$Study_Network_id),
#       slope = tt$estimate[2],
#       CI_low = tt$conf.low[2],
#       CI_high = tt$conf.high[2],
#       p = tt$p.value[2],
#       n = nrow(df)
#     )
#   })
# 
# 
# library(ggplot2)
# 
# forest<-ggplot(effect_table,
#                aes(x = slope,
#                    y = reorder(Study_Network_id, slope))) +
#   
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#   
#   geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
#                  height = 0.2,
#                  color = "grey40") +
#   
#   geom_point(aes(size = n),
#              color = "#2C7FB8",
#              alpha = 0.8) +
#   
#   theme_classic(base_size = 3) +
#   
#   labs(
#     x = "Effect size: Abundance → richness capture",
#     y = "Network"
#   ) +
#   
#   theme(legend.position = "none")
# 
# ggsave("./result_260526/forest.png", forest, width = 9, height = 15, units = "in", dpi = 300)

