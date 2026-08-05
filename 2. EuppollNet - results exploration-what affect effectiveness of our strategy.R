
setwd("D:/Chap1_TargetPlant_to_mo-nitor")
library(dplyr)
library(ggplot2) #For plotting
library(giscoR) #For plotting
library(sf) #For handling coordinates and plotting
library(patchwork)#For plotting (binding plots)
library(ggstar) #For plotting (cool shapes)
library(scales) #For plotting (decimals on axes)
library(tidyr)
library(viridis)
library(stringr)
library(vegan)

####################################################################################

#=============================RESULT 2 Exploring Op1==================================

#Factors influencing the effectiveness of subsampling

#################################################################################

##-----------(1) hypothesis 1 : 10 plant not enough for large networks------------


data_count_scaled<-readRDS("data_count_scaled_published_0526.rds")
data_interact<-readRDS("data_interact_published_0526.rds")

data_interact<-data_interact%>%
  mutate(
    Plant_accepted_name = str_replace_all(Plant_accepted_name, "×", "") %>%  # 去掉 ×
      str_squish()  # 去掉多余空格
  )

data_count_scaled_species <- data_count_scaled %>%
  filter(!is.na(Plant_species)) %>%       # 排除 NA
  filter(str_detect(Plant_species, " "))  # 保留含空格的名字（双名）

data_merge<- merge(data_interact, data_count_scaled[,c("Flower_data_merger","Flower_count_scaled","Plant_species","Study_Network_id")], 
                   by = "Flower_data_merger",all = TRUE)%>%
  filter(!is.na(Flower_data_merger))%>%
  mutate(
    Study_Network_id = coalesce(Study_Network_id.x, Study_Network_id.y)
  ) %>%
  dplyr::select(-Study_Network_id.x, -Study_Network_id.y)%>%
  mutate(Interaction_addup = ifelse(is.na(Interaction_addup), 0, Interaction_addup))%>% # 替换 `Interaction_addup` 为 NA 的值为 0
  mutate(
    Plant_accepted_name = str_squish(str_replace_all(replace_na(Plant_accepted_name, ""), "×", ""))
  ) 


unique(data_interact$Study_Network_id)
unique(data_count_scaled_species$Study_Network_id)
unique(data_merge$Study_Network_id)

percent_10<-readRDS("percent_10.rds")

#Top 10 common
data_count_scaled_species %>%
  filter(is.na(Flower_count_scaled!=0))

top_10_species <- data_count_scaled_species %>%
  filter(Flower_count_scaled!=0)%>%
  distinct(Study_Network_id, Plant_species,Flower_count_scaled, .keep_all = TRUE)%>%
  group_by(Study_Network_id) %>%
  arrange(Study_Network_id, desc(Flower_count_scaled)) %>%
  top_n(10, wt = Flower_count_scaled)

n_top10<-top_10_species%>%
  group_by(Study_Network_id)%>%
  summarise(n_top10=n_distinct(Plant_species,na.rm = TRUE))
print(n_top10)

n_all_plant<-data_count_scaled %>%
  filter(Flower_count_scaled!=0)%>%
  group_by(Study_Network_id)%>%
  summarise(n_all_plant=n_distinct(Plant_species,na.rm = TRUE))


plant_percent<-merge(n_top10, n_all_plant,by.x = "Study_Network_id", by.y = "Study_Network_id", all = TRUE)%>%
  mutate(plant_percent=n_top10/n_all_plant)
#merge
plotplotplot <- merge(percent_10,plant_percent,by.x = "Study_Network_id", by.y = "Study_Network_id", 
                      all = TRUE)
plotplotplot_clean <- plotplotplot %>% drop_na()


# check number of sites
num_points_before <- nrow(plotplotplot_clean)
print(paste("Number of points before jitter:", num_points_before))#389 networks

# regression
model1<- lm(Abundant_Top10 ~ plant_percent, data = plotplotplot_clean)
model2 <- lm(Abundant_Top10 ~ n_all_plant, data = plotplotplot_clean)
summary(model1)
summary(model2)

###########
#传粉者交互多度
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
# 
# # regression
# model3 <- lm(percent ~ plant_percent, data = Interaction_captured)
#model4 <- lm(percent ~ n_all_plant, data = Interaction_captured)
# 
# summary(model3)
#summary(model4)
#############################3
# 添加分组标签，并确保列名一致
# df1 <- plotplotplot_clean %>%
#   mutate(
#     percent = Abundant_Top10,
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
df1 <- plotplotplot_clean %>%
  mutate(
    percent = Abundant_Top10,
    Group = "Pollinator Richness"
  ) %>%
  dplyr::select(n_all_plant, percent, Group)

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


#####################################

##-----------(2) hypothesis 2 : total plant richness total pollinator richness , pollinator evenness, NODF affect the effectiveness of our strategies------------
# 
# ################################

# # for all strategies
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(vegan)
library(bipartite)

# 
# # -------------------- function  --------------------
# make_long <- function(df, cols, labels, xvar) {
#   df %>%
#     pivot_longer(
#       cols = all_of(cols),
#       names_to = "Top_Plant_Level",
#       values_to = "Richness_Value"
#     ) %>%
#     filter(!is.na(Richness_Value), !is.na(.data[[xvar]])) %>%
#     mutate(
#       Top_Plant_Level = factor(Top_Plant_Level,
#                                levels = cols,
#                                labels = labels)
#     )
# }
# 
# get_reg_lines <- function(df, formula) {
#   if(nrow(df) == 0) {
#     warning("No data available for regression!")
#     return(tibble(Top_Plant_Level = character(), signif_line = character()))
#   }
#   df %>%
#     group_by(Top_Plant_Level) %>%
#     summarise(lm_res = list(lm(formula, data = cur_data_all())), .groups = "drop") %>%
#     rowwise() %>%
#     mutate(
#       p_value = summary(lm_res)$coefficients[2,4],
#       signif_line = ifelse(p_value < 0.05, "solid", "dashed")
#     ) %>%
#     dplyr::select(Top_Plant_Level, signif_line)
# }
# 
# make_plot <- function(df, xvar, ylab, xlab, reg_lines, colors, base_theme, label_y_offsets) {
#   if(nrow(df) == 0) return(NULL)
#   ggplot(df, aes_string(x = xvar, y = "Richness_Value", color = "Top_Plant_Level")) +
#     geom_point(shape = 1, size = 2, position = position_jitter(width = 0.03, height = 0.5)) +
#     geom_smooth(method = "lm", se = FALSE, aes(linetype = Top_Plant_Level), linewidth = 0.8) +
#     scale_linetype_manual(values = setNames(reg_lines$signif_line, reg_lines$Top_Plant_Level)) +
#     stat_cor(aes(color = Top_Plant_Level), method = "pearson",
#              size = 3, label.x.npc = 0.6, label.y.npc = 0.05) +
#     scale_color_manual(values = colors, breaks = names(colors)) +
#     labs(y = ylab, x = xlab) +
#     base_theme +
#     theme(plot.margin = margin(5, 15, 5, 5)) 
# }
# 
# # -------------------- color --------------------
# my_colors <- c("Top 10"="#66C2A5", "Top 5"="#4E79A7", "Top 3"="#F28E2B",
#                "flw 5"="#4E79A7", "flw 3"="#F28E2B",
#                "Random 10"="#66C2A5", "Random 5"="#4E79A7", "Random 3"="#F28E2B")
# 
# base_theme <- theme_classic(base_size = 12) +
#   theme(plot.title = element_text(hjust=0.5, face="bold", size=15),
#         axis.title = element_text(face="bold", size=15),
#         axis.text = element_text(color="black", size=12),
#         legend.position = "none",
#         panel.background = element_blank(),
#         panel.grid = element_blank())
# 
# -------------------- Evenness calculation--------------------
evenness_df <- data_interact %>%
  filter(!is.na(Interaction_addup)) %>%
  group_by(Study_Network_id, Pollinator_accepted_name) %>%
  summarise(Total_interaction = sum(as.numeric(Interaction_addup)), .groups = "drop") %>%
  pivot_wider(names_from = Pollinator_accepted_name,
              values_from = Total_interaction,
              values_fill = 0) %>%
  rowwise() %>%
  mutate(
    H = diversity(c_across(-Study_Network_id), index = "shannon"),
    S = specnumber(c_across(-Study_Network_id)),
    Evenness = H / log(S)
  ) %>%
  dplyr::select(Study_Network_id, Evenness)

# --------------------  Nestedness calculation--------------------
library(bipartite)
nestedness_df <- data_interact %>%
  group_by(Study_Network_id) %>%
  summarise(
    Nestedness = {
      mat <- table(Pollinator_accepted_name, Plant_accepted_name)
      mat_bin <- ifelse(mat > 0, 1, 0)
      if(nrow(mat_bin) > 1 & ncol(mat_bin) > 1) {
        bipartite::networklevel(mat_bin, index = "nestedness")
      } else NA_real_
    },
    .groups = "drop"
  )


result_all<-read.csv("result_all_published.csv",header=TRUE)%>%
  rename(Abun_Top5_Rich = pollinator_count.x,
         Abun_Top3_Rich = pollinator_count.y,
         Abun_Top10_Rich = pollinator_count.x.x,
         FlwShape_Top5_Rich = pollinator_count.y.y,
         FlwShape_Top3_Rich = pollinator_count.x.x.x,
         Random_10_Rich = pollinator_count.y.y.y,
         Random_5_Rich = pollinator_count.x.x.x.x,
         Random_3_Rich = pollinator_count.y.y.y.y,
         Total_Rich = total_pollinator_count.x)

colnames(result_all)
# # -------------------- plotting --------------------
# plot_richness_relation <- function(metric_df, metric_name, metric_label, save_path_prefix) {
#   
#   # abundance-based
#   plot_abun <- result_all[, c("Study_Network_id","Abundant_Top10","Abundant_Top5","Abundant_Top3")] %>%
#     left_join(metric_df, by = "Study_Network_id")
#   
#   abun_long <- make_long(plot_abun,
#                          cols = c("Abundant_Top10","Abundant_Top5","Abundant_Top3"),
#                          labels = c("Top 10","Top 5","Top 3"),
#                          xvar = metric_name)
#   
#   reg_abun <- get_reg_lines(abun_long, formula = as.formula(paste0("Richness_Value ~ ", metric_name)))
#   
#   p_abun <- make_plot(abun_long, metric_name,
#                       "",#Percent of Pollinator Richness Captured (Top-N)
#                       metric_label, reg_abun, my_colors, base_theme, c(0.98,0.96,0.94))
#   
#   # flower-shape-based
#   plot_shape <- result_all[, c("Study_Network_id","FlwShape_Top5","FlwShape_Top3")] %>%
#     left_join(metric_df, by = "Study_Network_id")
#   
#   shape_long <- make_long(plot_shape,
#                           cols = c("FlwShape_Top5","FlwShape_Top3"),
#                           labels = c("flw 5","flw 3"),
#                           xvar = metric_name)
#   
#   reg_shape <- get_reg_lines(shape_long, formula = as.formula(paste0("Richness_Value ~ ", metric_name)))
#   
#   p_shape <- make_plot(shape_long, metric_name,
#                        "",#Percent of Pollinator Richness Captured (Flw-N)
#                        metric_label, reg_shape, my_colors, base_theme, c(0.98,0.96))
#   
#   # ---------- random-based ----------
#   plot_random <- result_all[, c("Study_Network_id",
#                                 "Random_10",
#                                 "Random_5",
#                                 "Random_3")] %>%
#     left_join(metric_df, by = "Study_Network_id")
#   
#   random_long <- make_long(
#     plot_random,
#     cols   = c("Random_10","Random_5","Random_3"),
#     labels = c("Random 10","Random 5","Random 3"),
#     xvar   = metric_name
#   )
#   
#   reg_random <- get_reg_lines(
#     random_long,
#     formula = as.formula(paste0("Richness_Value ~ ", metric_name))
#   )
#   
#   p_random <- make_plot(
#     random_long,
#     metric_name,
#     "",
#     metric_label,
#     reg_random, my_colors, base_theme, c(0.98,0.96,0.94)
#   )
#   
#   
#   # 拼图并保存
#   final_plot <- plot_grid(p_abun, p_shape, p_random, ncol = 3)
#   ggsave(paste0(save_path_prefix, "_", tolower(metric_name), ".png"),
#          final_plot, width = 9.8, height = 4.5, units = "in", dpi = 300)
#   
#   return(final_plot)
# }
# 
# # -------------------- 生成图像 --------------------
# getwd()
# plot_evenness <- plot_richness_relation(
#   metric_df = evenness_df,
#   metric_name = "Evenness",
#   metric_label = "Pollinator Evenness (Within-Network)",
#   save_path_prefix = "result251105_published/evenness_topN"
# )
# 
# plot_nestedness <- plot_richness_relation(
#   metric_df = nestedness_df,
#   metric_name = "Nestedness",
#   metric_label = "Network Nestedness",
#   save_path_prefix = "result251105_published/nestedness_topN"
# )
# 

plot_rich_site<-result_all[,c("Study_Network_id","Total_Rich")]
head(plot_rich_site)

# plot_totalrich <- plot_richness_relation(
#   metric_df = plot_rich_site,   # 这里是你包含 Total_Rich 的表
#   metric_name = "Total_Rich",       # 横轴变量
#   metric_label = "Total Pollinator Richness",
#   save_path_prefix = "result251105_published/totalrich_topN"
# )
# 
# 
# print(plot_nestedness)
# plot_evenness
# 
# library(cowplot)
# 
# final_three_plot <- plot_grid(
#   plot_totalrich,
#   plot_evenness,
#   plot_nestedness,
#   ncol = 1,             # 竖排
#   align = "v",
#   rel_heights = c(1, 1, 1)
# )
# 
# 
# final_three_plot_with_ylabel <- plot_grid(
#   ggdraw() + draw_label(
#     "Percent of Pollinator Richness Captured in Subnetwork",
#     angle = 90,
#     fontface = "bold",
#     size = 14
#   ),
#   final_three_plot,
#   ncol = 2,
#   rel_widths = c(0.12, 1)  # increase left margin for long label
# )
# 
# 
# print(final_three_plot_with_ylabel)
# 
# 
# # 导出高分辨率 PNG 文件
# ggsave("/Chap1_TargetPlant_to_monitor/result_260526/three_panel_totalrich_evenness_nestedness.png",
#        final_three_plot_with_ylabel,
#        width = 13.5, height = 11, dpi = 300)

###########################

##-------------Show solely Abundant-top10-----------

##########################################

  
  # ---------- Pearson ----------
  get_cor_info <- function(df, x, y) {
    test <- cor.test(df[[x]], df[[y]], method = "pearson")
    
    list(
      r = unname(test$estimate),
      p = test$p.value,
      line_type = ifelse(test$p.value < 0.05, "solid", "dashed")
    )
  }
#-----------linner regression result--------
get_lm_info <- function(df, x, y) {
  
  model <- lm(df[[y]] ~ df[[x]])
  
  p <- summary(model)$coefficients[2, 4]
  
  list(
    p = p,
    line_type = ifelse(p < 0.05, "solid", "dashed")
  )
}

#-------------------------------
plot_top10_relation <- function(metric_df, metric_name, metric_label) {
  
  df <- result_all[, c("Study_Network_id", "Abundant_Top10")] %>%
    left_join(metric_df, by = "Study_Network_id") %>%
    filter(!is.na(Abundant_Top10), !is.na(.data[[metric_name]])) %>%
    mutate(
      Predictor = metric_label
    )
  
  # ---------- Pearson  ----------
  cor_info <- get_cor_info(df, metric_name, "Abundant_Top10")
  lm_info <- get_lm_info(df, metric_name, "Abundant_Top10")
  
  
  label = paste0(
    "R = ", round(cor_info$r, 2),
    "\np = ", format.pval(cor_info$p, digits = 2)
  )
  

  ggplot(
    df,
    aes_string(
      x = metric_name,
      y = "Abundant_Top10",
      color = "Predictor",
      shape = "Predictor"
    )
  ) +
    
    geom_point(
      size = 2,
      alpha = 0.7,
      position = position_jitter(
        width = 0.03,
        height = 0.5
      )
    ) +
    
    geom_smooth(
      aes(color = Predictor),
      method = "lm",
      se = FALSE,
      linewidth = 1,
      linetype = lm_info$line_type
    ) +
    
    scale_color_manual(
      values=setNames("#9E9AC8", metric_label)
    )+
    
    scale_shape_manual(
      values=setNames(1, metric_label)
    )+
    
    annotate(
      "text",
      x=Inf,
      y=-Inf,
      label=label,
      hjust=1.1,
      vjust=-0.5,
      size=3
    )+
    
    labs(
      x=stringr::str_wrap(metric_label,width=25),
      y="Percent of Pollinator Richness Captured",# or set as NULL
      color=NULL,
      shape=NULL
    )+
    
    theme_classic(base_size = 12) +
  
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.width = unit(1.5,"cm")
  )
}

# -------------plant data
n_all_plant <- data_count_scaled %>%
  filter(Flower_count_scaled != 0) %>%
  group_by(Study_Network_id) %>%
  summarise(
    Plant_Richness = n_distinct(Plant_species, na.rm = TRUE)
  )


## covariance test

network_size <- data_interact %>%
  group_by(Study_Network_id) %>%
  summarise(
    Plant_interaction_Richness = n_distinct(Plant_accepted_name),
    Pollinator_Richness = n_distinct(Pollinator_accepted_name),
    Network_size = Plant_interaction_Richness + Pollinator_Richness
  )

# ==============================================================================
# 一键生成论文级报告：Abundant_Top10 回归分析
# 作者：你的名字
# 日期：2025年4月5日
# 功能：自动完成数据合并、建模、诊断、绘图、表格输出
# ==============================================================================

# 1. 加载必要包
library(dplyr)
library(flextable)
library(officer)
library(stargazer)
library(sandwich)
library(lmtest)
library(MASS)
library(broom)

# 2. 读取数据并合并
test <- n_all_plant %>%
  left_join(plot_rich_site, by = "Study_Network_id") %>%
  rename(Pollinator_Richness = Total_Rich) %>%
  left_join(evenness_df, by = "Study_Network_id") %>%
  left_join(nestedness_df, by = "Study_Network_id") %>%
  left_join(result_all[, c("Study_Network_id", "Abundant_Top10")], by = "Study_Network_id")%>%
  left_join(network_size[, c("Study_Network_id", "Network_size")],by = "Study_Network_id" )

# 3. 检查缺失值
sum(is.na(test$Abundant_Top10))  # 确保无缺失

m1 <- lm(
  Abundant_Top10 ~ Plant_Richness + Pollinator_Richness,
  data=test
)

# 4. 普通线性回归
m2 <- lm(Abundant_Top10 ~ Plant_Richness + Pollinator_Richness + Nestedness, data = test)


# 共线性检验
library(car)
anova(m1, m2)
summary(m2)
R2_m1 <- summary(m1)$r.squared
R2_m2 <- summary(m2)$r.squared

R2_m2 - R2_m1


vif_values <- vif(m2)
vif_values
vif_table <- data.frame(
  Variable=names(vif_values),
  VIF=round(vif_values,2)
)

# 5. 稳健回归
m2_robust <- rlm(Abundant_Top10 ~ Plant_Richness + Pollinator_Richness + Nestedness, data = test)
vif(m2_robust)
summary(m2_robust)

# 6. 模型诊断
# 残差图
par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))

# 残差正态性
shapiro_test <- shapiro.test(residuals(m2))
bptest_test <- bptest(m2)

# Cook's distance
cooks_d <- cooks.distance(m2)
n <- length(cooks_d)
threshold <- 4 / n
influential <- sum(cooks_d > threshold)

#The effectiveness of abundance-based plant selection was primarily constrained 
# by plant community size, with larger plant communities reducing the proportion of 
# pollinator richness captured by a fixed number of monitored species. However, 
# network structure also played an important role: more nested networks showed higher
# capture efficiency, suggesting that abundant plants represent pollinator communities 
# more effectively when interactions are concentrated around core plant species.

# 7. 提取稳健回归结果（用于表格）
# 使用 tidy() + 手动计算 p 值
coef_tab <- as.data.frame(
  summary(m2_robust)$coefficients
)

coef_tab <- coef_tab[-1,]


tidy_robust <- data.frame(
  Variable = rownames(coef_tab),
  Estimate_SE = paste0(
    round(coef_tab$Value,3),
    " ± ",
    round(coef_tab$`Std. Error`,3)
  ),
  t = round(coef_tab$`t value`,2),
  P = 2*(1-pt(abs(coef_tab$`t value`),df=264))
)


tidy_robust$P <- ifelse(
  tidy_robust$P <0.001,
  "<0.001",
  sprintf("%.3f",tidy_robust$P)
)

rownames(tidy_robust)<-NULL

tidy_robust <- left_join(
  tidy_robust,
  vif_table,
  by=c("Variable"="Variable")
)

# =========================
# 8. Word报告
# =========================

doc <- read_docx()

doc <- doc %>%
  body_add_par(
    "Abundant Top 10 Species Analysis",
    style="heading 1"
  ) %>%
  body_add_par(
    "This report presents the results of linear and robust regression models predicting abundant species capture.",
    style="Normal"
  )


doc <- doc %>%
  body_add_par(
    paste0(
      "Model diagnostics showed heteroscedasticity ",
      "(Breusch-Pagan test: p < 0.001) and deviations ",
      "from residual normality (Shapiro-Wilk test: p = 0.001). ",
      "All predictors showed low multicollinearity ",
      "(VIF < 2). Robust regression was therefore used ",
      "for final inference."
    ),
    style="Normal"
  )


# flextable

ft <- flextable(tidy_robust)

ft <- set_header_labels(
  ft,
  Variable="Variable",
  Estimate_SE="Estimate ± SE",
  t="t",
  P="p-value",
  VIF="VIF"
)


ft <- ft %>%
  align(
    align="center",
    part="all"
  ) %>%
  bold(
    part="header"
  ) %>%
  font(
    fontname="Times New Roman",
    part="all"
  ) %>%
  fontsize(
    size=10,
    part="all"
  ) %>%
  theme_vanilla() %>%
  autofit()



doc <- doc %>%
  body_add_flextable(ft) %>%
  body_add_par(
    "Table 1. Robust regression results: Predictors of Abundant Top 10 Species (%)",
    style="Normal"
  ) %>%
  body_add_par(
    "Note: SE = standard error; p-values are two-tailed; VIF values indicate multicollinearity diagnostics.",
    style="Normal"
  )


print(
  doc,
  target="result_260526/Influence_factor_Top10_Report.docx"
)


########## marginal effect
library(ggeffects)
library(ggplot2)

plant_eff <- ggpredict(
  m2,
  terms = "Plant_Richness"
)

p1 <- ggplot(plant_eff,
             aes(x=x,
                 y=predicted)) +
  geom_ribbon(
    aes(ymin=conf.low,
        ymax=conf.high),
    alpha=0.2
  ) +
  geom_line(
    color = "#66C2A5",
    linewidth=1
  ) +
  labs(
    x="Plant richness",
    y="Pollinator richness captured (%)"
  ) +
  theme_classic(base_size=13)

p1

nested_eff <- ggpredict(
  m2,
  terms="Nestedness"
)


p2 <- ggplot(nested_eff,
             aes(x=x,
                 y=predicted)) +
  geom_ribbon(
    aes(ymin=conf.low,
        ymax=conf.high),
    alpha=0.2
  ) +
  geom_line(
    color = "#66C2A5",
    linewidth=1
  ) +
  labs(
    x="Nestedness",
    y="Pollinator richness captured (%)"
  ) +
  theme_classic(base_size=13)

p2
# -------------------- plot--------------------
p_top10_plant <- plot_top10_relation(
  n_all_plant,
  "Plant_Richness",
  "Total Plant Richness"
)


p_top10_rich <- plot_top10_relation(
  plot_rich_site,
  "Total_Rich",
  "Total Pollinator Richness"
)

p_top10_even <- plot_top10_relation(
  evenness_df,
  "Evenness",
  "Pollinator Evenness (Within-Network)"
)

p_top10_nested <- plot_top10_relation(
  nestedness_df,
  "Nestedness",
  "Network Nestedness"
)


top10_four_plot <- plot_grid(
  p_top10_plant,
  p_top10_rich,
  p_top10_even,
  p_top10_nested,
  ncol = 2,
  nrow = 2,
  align = "hv"
)

# -------------------- 添加共享Y轴标题 --------------------
top10_with_ylabel <- plot_grid(
  ggdraw() +
    draw_label(
      "Percent of Pollinator Richness Captured",
      angle = 90,
      fontface = "bold",
      size = 14
    ),
  
  top10_four_plot,
  
  ncol = 2,
  rel_widths = c(0.12, 1)
)

# 显示
print(top10_with_ylabel)


ggsave(
  "result_260526/top10_influence_factor_4penal.png",
  top10_with_ylabel,
  width = 6.5, height = 6,
  dpi = 300
)

############新图
plant_df <- result_all[, c("Study_Network_id", "Abundant_Top10")] %>%
  left_join(
    n_all_plant,
    by = "Study_Network_id"
  ) %>%
  filter(
    !is.na(Abundant_Top10),
    !is.na(Plant_Richness)
  ) %>%
  mutate(
    Predictor = "Plant richness",
    Richness = Plant_Richness
  )

pollinator_df <- result_all[, c("Study_Network_id", "Abundant_Top10")] %>%
  left_join(
    plot_rich_site,
    by="Study_Network_id"
  ) %>%
  filter(
    !is.na(Abundant_Top10),
    !is.na(Total_Rich)
  ) %>%
  mutate(
    Predictor="Pollinator richness",
    Richness=Total_Rich
  )

richness_compare <- bind_rows(
  plant_df[,c("Study_Network_id",
              "Abundant_Top10",
              "Predictor",
              "Richness")],
  
  pollinator_df[,c("Study_Network_id",
                   "Abundant_Top10",
                   "Predictor",
                   "Richness")]
)

stats_richness <- richness_compare %>%
  group_by(Predictor) %>%
  summarise(
    r = cor(Richness,
            Abundant_Top10,
            method="pearson"),
    p = summary(
      lm(
        Abundant_Top10 ~ Richness
      )
    )$coefficients[2,4]
  )

stats_richness

p_top10_richness_compare <- ggplot(
  richness_compare,
  aes(
    x = Richness,
    y = Abundant_Top10,
    group = Predictor,
    color = Predictor,
    shape = Predictor
  )
) +
  
  geom_point(
    size = 2,
    alpha = 0.7,
    position = position_jitter(
      width = 0.03,
      height = 0.5
    )
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  
  scale_color_manual(
    values = c(
      "Plant richness" = "#66C2A5",
      "Pollinator richness" = "#E69F00"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Plant richness" = 1,        # 空心圆
      "Pollinator richness" = 2    # 空心三角
    )
  ) +
  
  labs(
    x = "Richness",
    y = "Percent of Pollinator Richness Captured",
    color = NULL,
    shape = NULL
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.1,
    vjust = -0.5,
    label =
      paste0(
        "Plant:\nR=",
        round(stats_richness$r[1],2),
        ", p=",
        format.pval(stats_richness$p[1],digits=2),
        "\n",
        "Pollinator:\nR=",
        round(stats_richness$r[2],2),
        ", p=",
        format.pval(stats_richness$p[2],digits=2)
      ),
    size = 3
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.width = unit(1.5,"cm")
  )


p_top10_richness_compare

ggsave(
  "result_260526/top10_influence_factor_plantpoll.png",
  p_top10_richness_compare,
  width = 6, height = 5.5,
  dpi = 300
)

p_top10_nested <- plot_top10_relation(
  nestedness_df,
  "Nestedness",
  "Network Nestedness"
)
p_top10_nested
ggsave(
  "result_260526/top10_influence_factor_nest.png",
  p_top10_nested,
  width = 6, height = 5.5,
  dpi = 300
)
###############################################################################

##-----------(3) hypothesis 3 : Favored plant sp by pollinators is rare (plant abundance-interaction times) 

################################
library(ggplot2)
library(dplyr)
library(broom)
library(purrr)

traits<-read.csv("test.merge.trait.csv", header = TRUE, fileEncoding = "UTF-8")
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
  wt <- wilcox.test(Abundant_Top10 ~ group, data = df_test)
  
  # effect size (Cliff’s delta approx via median diff)
  eff <- with(df_test,
              median(Abundant_Top10[group=="Present"], na.rm=TRUE) -
                median(Abundant_Top10[group=="Absent"], na.rm=TRUE))
  
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

shapiro.test(attract_percent_10$Abundant_Top10[attract_percent_10$present_group=="Present"])
shapiro.test(attract_percent_10$Abundant_Top10[attract_percent_10$present_group=="Absent"])

wilcox.test(Abundant_Top10 ~ present_group, data = attract_percent_10)
unique(attract_percent_10$present_group)

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
wilcox_res <- wilcox.test(Abundant_Top10 ~ present_group, data = attract_percent_10)

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
    Mean_Percentage = mean(Abundant_Top10, na.rm = TRUE),
    SD_Percentage = sd(Abundant_Top10, na.rm = TRUE),
    n_non_missing = sum(!is.na(Abundant_Top10)),
    .groups = "drop"
  )

group_colors <- c("Present" = "#E41A1C", "Absent" = "#377EB8")

# 4️⃣ y 位置
y_max <- max(attract_percent_10$Abundant_Top10, na.rm = TRUE)
y_pos <- y_max * 1.1

# =========================
# 5️⃣ FINAL PLOT
# =========================
plot <- ggplot(attract_percent_10,
               aes(x = present_group,
                   y = Abundant_Top10,
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
    Mean = mean(Abundant_Top10, na.rm = TRUE),
    SD = sd(Abundant_Top10, na.rm = TRUE),
    n = sum(!is.na(Abundant_Top10)),
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
adonis_res <- adonis2(d ~ flw_shape_revised,
                      data = meta2,
                      permutations = 999,
                      strata = meta2$Study_Network_id)

adonis_res_study_ctrl <- adonis2(d ~ Study_id +flw_shape_revised,
                      data = meta2,
                      permutations = 999,
                      strata = meta2$Study_Network_id,
                      by = "margin")
# 是否
head(meta2)
print(adonis_res)
print(adonis_res_study_ctrl)
saveRDS(adonis_res,"adonis_res.rds")

adonis_res<-readRDS("adonis_res.rds")
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

strongdata<-plottest%>%filter(Study_Network_id %in% strong_network$Network)%>%select(Study_Network_id,Abundant_Top5,FlwShape_Top5)

t.test(strongdata$Abundant_Top5,
       strongdata$FlwShape_Top5,
       paired = TRUE)
wilcox.test(strongdata$Abundant_Top5,
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

# saveRDS(df, "mean_dis_plot.rds")
df<-readRDS("mean_dis_plot.rds")

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

##----------- NMDS result output--------------

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





# ####################################################################################
# 
# #####===================== exploring OP2 ====================================
# 
# ############################################################################
# 
# 
# ### (1) For OP1 and OP2, which performs better and why?
# library(ggplot2)
# library(dplyr)
# #提取结果
# result_all<-read.csv("result_all_published.csv",header=TRUE)
# plottest<-result_all[,c("Study_Network_id","Abundant_Top10","Abundant_Top5","Abundant_Top3","FlwShape_Top5","FlwShape_Top3")]
# plot_OP1_OP2<-plottest[,]
# # 进行一元线性回归
# model <- lm(Abundant_Top5 ~ FlwShape_Top5, data = plot_OP1_OP2)
# summary(model)
# 
# ####新的筛选标准：R=|a-b|/(a+b) 用来比较二者差异的大小，筛选标准R>0.3(文献？)
# # Outliner overview with threshold R = 0.3
# 
# # Calculate R for all entries
# plot_OP1_OP2 <- plot_OP1_OP2 %>%
#   mutate(R = abs(Abundant_Top5 - FlwShape_Top5) / (Abundant_Top5 + FlwShape_Top5),
#          op2minusop1 = FlwShape_Top5 - Abundant_Top5)
# 
# # OP1 better (Abundant_Top5 significantly higher than FlwShape_Top5)
# Op1better <- plot_OP1_OP2 %>%
#   filter(Abundant_Top5 > FlwShape_Top5 & R > 0.3 & !is.na(R))
# 
# # OP2 better (FlwShape significantly higher than Abundant_Top10)
# Op2better <- plot_OP1_OP2 %>%
#   filter(FlwShape_Top5 > Abundant_Top5 & R > 0.3 & !is.na(R))
# 
# # 对 Op1 / Op2 表现更好的网络筛选的 top 5 物种，补充 Plant_accepted_name 和花型
# 
# op2_sp_choosen_op1 <- top_5_species %>%
#   filter(Study_Network_id %in% Op1better$Study_Network_id |
#            Study_Network_id %in% Op2better$Study_Network_id) %>%
#   distinct(Study_Network_id, Plant_species, .keep_all = TRUE)%>%
#   left_join(
#     traits,
#     by = c("Plant_species" = "Plant_accepted_name")
#   ) %>%
#   group_by(Study_Network_id) %>%
#   summarise(
#     op1_sp_choosen = n_distinct(Plant_species),
#     Plant_species = paste(unique(Plant_species), collapse = "; "),
#     flw_shapes = paste(unique(flw_shape_revised), collapse = "; "),
#     .groups = "drop"
#   )
# 
# 
# 
# op2_sp_choosen_op2 <- abun_species_top5 %>%  # *此数据见 rscript.op1-3
#   filter(Study_Network_id %in% Op1better$Study_Network_id | 
#            Study_Network_id %in% Op2better$Study_Network_id) %>%
#   distinct(Study_Network_id, Plant_species, .keep_all = TRUE)%>%
#   left_join(
#     traits,
#     by = c("Plant_species" = "Plant_accepted_name")
#   ) %>%
#   group_by(Study_Network_id) %>%
#   summarise(
#     op2_sp_choosen = n_distinct(Plant_species),
#     flw_shape_count = n_distinct(flw_shape_revised)   # 新增列
#   )
# 
# 
# # Optional: Combine results for easier interpretation
# # 为每种策略选择的物种数添加分类标签
# op2_df <- op2_sp_choosen_op2 %>%
#   mutate(op2_sp_choosen = op2_sp_choosen) %>%
#   dplyr::select(Study_Network_id, op2_sp_choosen)
# 
# op1_df <- op2_sp_choosen_op1 %>%
#   mutate(op1_sp_choosen = op1_sp_choosen) %>%
#   dplyr::select(Study_Network_id, op1_sp_choosen)
# 
# # 合并两个策略物种数据
# species_combined <- full_join(op2_df, op1_df, by = "Study_Network_id")
# 
# # 加入性能指标（如 R、Top5 等）
# op2_sp_choosen_combined <- species_combined %>%
#   left_join(plot_OP1_OP2, by = "Study_Network_id")%>%
#   mutate(
#     Category = case_when(
#       Abundant_Top5 > FlwShape_Top5 & R > 0.3 ~ "OP1better",
#       FlwShape_Top5 > Abundant_Top5 & R > 0.3 ~ "OP2better",
#       TRUE ~ "Similar"
#     )
#   )
# 
# op2_sp_choosen_combined
# write.csv(op2_sp_choosen_combined,"op1_vs_op2.250807.csv")
# 
# ##optional:看看这些op2表现更好的network里选了哪些种?
# 
# match<-na.omit(unique(data_merge[,c("Plant_accepted_name", "Plant_species")]))
# 
# op1_species <- top_5_species %>%
#   left_join(match,by = "Plant_species")%>%
#   dplyr::select(Study_Network_id, Plant_accepted_name,Plant_species ) %>%
#   mutate(method = "Op1")
# 
# op2_species <- abun_species_top5 %>%
#   left_join(match,by = "Plant_species")%>%
#   dplyr::select(Study_Network_id, Plant_accepted_name,Plant_species) %>%
#   mutate(method = "Op2")
# 
# op1_species%>%filter(is.na(Plant_species ))
# op2_species%>%filter(is.na(Plant_species ))
# 
# # 合并
# combined_species <- bind_rows(op1_species, op2_species)
# target_networks <- op2_sp_choosen_combined$Study_Network_id
# 
# interaction_sum <- data_merge %>%
#   ungroup() %>%
#   filter(Study_Network_id %in% target_networks) %>%
#   dplyr::select(Study_Network_id, Plant_accepted_name, Plant_species, Pollinator_accepted_name, Interaction_addup) %>%
#   group_by(Study_Network_id, Plant_accepted_name, Plant_species) %>%
#   summarise(Interaction_addup_sum = sum(Interaction_addup, na.rm = TRUE), .groups = "drop")
# 
# interaction_sum%>%filter(is.na(Interaction_addup_sum))
# unique(interaction_sum$Study_Network_id)
# 
# # 统计每个物种在不同方法中的出现情况
# species_comparison <- combined_species %>%
#   filter(Study_Network_id %in% target_networks) %>%
#   group_by(Study_Network_id, Plant_accepted_name, Plant_species) %>%
#   summarise(methods = paste(unique(method), collapse = ", "), .groups = "drop") %>%
#   mutate(category = case_when(
#     methods == "Op1, Op2" ~ "Both",
#     methods == "Op1" ~ "Op1 only",
#     methods == "Op2" ~ "Op2 only"
#   )) %>%
#   left_join(op2_sp_choosen_combined %>% dplyr::select(Study_Network_id, Category), by = "Study_Network_id")%>%
#   left_join(interaction_sum, by = c("Study_Network_id", "Plant_species")) %>%
#   mutate(Interaction_addup_sum = ifelse(is.na(Interaction_addup_sum), 0, Interaction_addup_sum))%>% # 替换 `Interaction_addup` 为 NA 的值为 0
#   left_join(traits %>% rename(Plant_species = "Plant_accepted_name"), by = "Plant_species") %>%
#   left_join(flower_count_match, by = c("Study_Network_id", "Plant_species")) %>%
#   distinct()
# 
# 
# unique(species_comparison$Study_Network_id)
# species_comparison
# 
# 
# a<-species_comparison 
# unique(a$Study_Network_id)
# 
# # 1️⃣ 统计每个 network 中不同方法的物种数
# network_counts_in_sc <- species_comparison %>%
#   filter(!is.na(Interaction_addup_sum)) %>%
#   separate_rows(methods, sep = ", ") %>%
#   mutate(methods = trimws(methods)) %>%  # 🧹 去掉多余空格
#   group_by(Study_Network_id, methods) %>%
#   summarise(num_unique_plants = n_distinct(Plant_accepted_name), .groups = "drop") %>%
#   pivot_wider(
#     names_from = methods,
#     values_from = num_unique_plants,
#     values_fill = 0,
#     names_prefix = "Num_"
#   )
# 
# # 2️⃣ 把统计结果合并回 species_comparison
# species_comparison2 <- species_comparison %>%
#   left_join(network_counts_in_sc, by = "Study_Network_id") %>%
#   mutate(Interaction_addup_sum = replace_na(Interaction_addup_sum, 0))
# 
# 
#  write.csv(species_comparison2,"E:/Chap1_TargetPlant_to_monitor/result250917_published/species_comparison0911.csv")
# 
#  
#  data_merge%>%filter(Study_Network_id=="30_Smith_lsk12")
#  
#  
#  
#  
# # 绘制散点图和回归线
# scatterplot <- ggplot(plot_OP1_OP2, aes(x = Abundant_Top5, y = FlwShape_Top5)) +
#   geom_point(color = "gray", size = 2, alpha = 0.7, position = position_jitter(width = 0.5, height = 0.5)) +  # 修改点的颜色和大小
#   geom_point(data = Op1better, aes(x = Abundant_Top5, y = FlwShape_Top5), color = "#1B9E77", size = 3, shape = 24, fill = "#1B9E77", alpha = 1) +  # 标记 OP1better 的点
#   geom_point(data = Op2better, aes(x = Abundant_Top5, y = FlwShape_Top5), color = "#9E9AC8", size = 3, shape = 21, fill = "#9E9AC8", alpha = 1) +  # 标记 OP2better 的点
#   geom_smooth(method = "lm", se = FALSE, color = "gray") +  # 回归线，不显示标准误差带
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray", size = 1) +  # 添加 y=x 线
#   labs(x = "Percent of Richness Captured by Option 1", y = "Percent of Richness Captured by Option 2") +
#   theme_classic(base_size=12) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
#     axis.title = element_text(face = "bold", size = 15),
#     axis.text = element_text(color = "black", size = 12),
#     legend.title = element_text(face = "bold", size = 11),
#     legend.text = element_text(size = 10),
#     legend.position = "right",
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     plot.margin = margin(5, 15, 5, 5)
#   )
# # 打印散点图和回归线
# print(scatterplot)
# ggsave("result251105_published/op2-op1.png", scatterplot, width = 5, height = 5, units = "in", dpi = 300)
# 
# 
# # 打印回归模型的摘要
# model <- lm(FlwShape_Top5 ~ Abundant_Top5, data = plot_OP1_OP2)
# summary(model)
# 
# ######################################################################
# 
# #分析op1更好的这8个网络，或op2更好的这七个网络，是因为什么? 怎么分析呢？
# 
# #1. 哪些植物？
# #hypothesis: op1 better because flower shape 
# data_merge%>%filter(Study_Network_id == "20_Hoiss_HO1")
# 
# library(dplyr)
# library(tidyr)
# 
# top_5_species
# abun_species_top5
# traits
# 
# op1better_plant<-top_5_species%>%
#   filter(Study_Network_id%in%Op1better$Study_Network_id)%>%left_join(traits %>% dplyr::select(Plant_accepted_name, flw_shape_revised),
#                                                                      by = c("Plant_species" = "Plant_accepted_name")) 
# 
# 
# op2better_plant<-abun_species_top5%>%
#   filter(Study_Network_id%in%Op2better$Study_Network_id)%>%left_join(traits %>% dplyr::select(Plant_accepted_name, flw_shape_revised),
#                                                                      by = c("Plant_species" = "Plant_accepted_name")) 
# #能否将两组网络都列出植物数量，花型数量，交互数量，传粉者多样性，交互次数。用data_merge计算。
# 
# library(dplyr)
# 
# # 确认 OP1better 和 OP2better 的网络
# op1_ids <- unique(Op1better$Study_Network_id)
# op2_ids <- unique(Op2better$Study_Network_id)
# 
# # 基于 data_merge 计算指标
# network_summary <- data_merge %>%
#   # 添加花型信息（若 traits 不在 data_merge 中）
#   left_join(traits %>% dplyr::select(Plant_accepted_name, flw_shape_revised),
#             by = c("Plant_accepted_name" = "Plant_accepted_name")) %>%
#   group_by(Study_Network_id) %>%
#   summarise(
#     Plant_richness = n_distinct(Plant_accepted_name),
#     FlwShape_richness = n_distinct(flw_shape_revised),
#     Pollinator_richness = n_distinct(Pollinator_accepted_name),
#     Interaction_richness = n_distinct(paste(Plant_accepted_name, Pollinator_accepted_name, sep = "_")),
#     Interaction_sum = sum(Interaction_addup, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     Category = case_when(
#       Study_Network_id %in% op1_ids ~ "OP1better",
#       Study_Network_id %in% op2_ids ~ "OP2better",
#       TRUE ~ "Other"
#     )
#   )
# 
# # 查看结果
# network_summary
# 
# # 可选：分别列出 OP1better 和 OP2better 的子表
# op1_summary <- network_summary %>% filter(Category == "OP1better")
# op2_summary <- network_summary %>% filter(Category == "OP2better")
# 
# ########################
# #预测：在op1表现更好的网络中植物多度与传粉者多度成正比，而在2中由于某类花型更attractive，考虑花型更能捕获
# 
# library(dplyr)
# library(ggplot2)
# 
# # 计算每个网络的回归显著性
# reg_lines <- abun_richness %>%
#   group_by(Study_Network_id) %>%
#   summarise(
#     lm_res = list(lm(Proportion_of_richness ~ Abundance_normalized, data = cur_data_all())),
#     Group = first(Group),  # 保留 Group 信息用于上色
#     .groups = "drop"
#   ) %>%
#   rowwise() %>%
#   mutate(
#     p_value = summary(lm_res)$coefficients[2,4],
#     signif_line = ifelse(p_value < 0.05, "solid", "dashed")
#   ) %>%
#   dplyr::select(Study_Network_id, Group, signif_line)
# 
# # 合并回原始数据
# abun_richness_plot <- abun_richness %>%
#   left_join(reg_lines, by = c("Study_Network_id", "Group"))
# 
# # 绘图
# set.seed(2025)
# plot1.2<-ggplot(abun_richness_plot, 
#        aes(x = Abundance_normalized, y = Proportion_of_richness, 
#            shape = Group, color = Group, group = Study_Network_id)) +
#   geom_point(size = 2, position = position_jitter(width = 0.03, height = 0.03)) +
#   geom_smooth(method = "lm", se = FALSE, aes(linetype = signif_line, color = Group), size = 0.8, show.legend = FALSE) +
#   scale_color_manual(values = c("op1better" = "#1B9E77",  
#                                 "op2better" = "#9E9AC8")) +  
#   scale_shape_manual(values = c("op1better" = 16,   # 实心圆
#                                 "op2better" = 17)) + # 实心三角形
#   scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) +
#   labs(x = "Flower Abundance (Normalized)",
#        y = "Percent of Richness Captured in Subnetwork") +
#   theme_classic(base_size=12) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
#     axis.title = element_text(face = "bold", size = 12),
#     axis.text = element_text(color = "black", size = 12),
#     legend.title = element_text(face = "bold", size = 11),
#     legend.text = element_text(size = 10),
#     legend.position = "right",
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     plot.margin = margin(5, 15, 5, 5)
#   )
# 
# ggsave("result251105_published/op1-2 explore2.png", plot1.2, width = 5.5, height = 4, units = "in", dpi = 300)
# 
# 
# ####################
# #假设1：网络中某类花型更attractive因此导致在这些网络中考虑花型能更好地捕捉传粉者
# #结论：无显著差异
# #假设2：网络中某类花型包含更多植物种因此导致在这些网络中考虑花型能更好地捕捉传粉者
# #结论：无显著差异
# 
# library(dplyr)
# 
# # 确保花型信息存在
# data_flw_poll <- data_merge %>%
#   left_join(traits %>% dplyr::select(Plant_accepted_name, flw_shape_revised),
#             by = "Plant_accepted_name") %>%
#   filter(!is.na(flw_shape_revised), !is.na(Pollinator_accepted_name))
# 
# # 计算：每个网络 × 每种花型 的传粉者多样性
# poll_div_by_flw <- data_flw_poll %>%
#   group_by(Study_Network_id, flw_shape_revised) %>%
#   summarise(
#     Pollinator_richness = n_distinct(Pollinator_accepted_name),
#     Interaction_sum = sum(Interaction_addup, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # 可选：计算花型在该网络内的植物数量（用于补充理解）
# flw_plant_count <- data_flw_poll %>%
#   group_by(Study_Network_id, flw_shape_revised) %>%
#   summarise(
#     Plant_richness = n_distinct(Plant_accepted_name),
#     .groups = "drop"
#   )
# 
# # 合并结果
# poll_div_by_flw_full <- left_join(poll_div_by_flw, flw_plant_count,
#                                   by = c("Study_Network_id", "flw_shape_revised"))%>%
#   group_by(Study_Network_id) %>%
#   mutate(
#     Total_poll_richness = sum(Pollinator_richness, na.rm = TRUE),
#     Pollinator_richness_percent = Pollinator_richness / Total_poll_richness * 100
#   ) %>%
#   ungroup()%>%
#   group_by(Study_Network_id) %>%
#   mutate(
#     Total_plant_richness = sum(Plant_richness, na.rm = TRUE),
#     Plant_richness_percent = Plant_richness / Total_plant_richness * 100
#   ) %>%
#   ungroup()
# 
# # 查看结果
# head(poll_div_by_flw_full)
# 
# # 保存结果
# write.csv(poll_div_by_flw_full, "pollinator_diversity_by_flower_shape.csv", row.names = FALSE)
# 
# 
# poll_div_by_flw_full <- poll_div_by_flw_full %>%
#   mutate(Group = case_when(
#     Study_Network_id %in% Op1better$Study_Network_id ~ "op1better",
#     Study_Network_id %in% Op2better$Study_Network_id ~ "op2better",
#     TRUE ~ "others"
#   )) %>%
#   filter(Group %in% c("op1better", "op2better"))
# 
# 
# library(dplyr)
# 
# poll_div_summary <- poll_div_by_flw_full %>%
#   group_by(Group, flw_shape_revised) %>%
#   summarise(
#     mean_richness = mean(Pollinator_richness, na.rm = TRUE),
#     sd_richness = sd(Pollinator_richness, na.rm = TRUE),
#     n_networks = n(),
#     .groups = "drop"
#   )
# 
# poll_div_summary
# 
# library(ggplot2)
# 
# ppoll<-ggplot(poll_div_by_flw_full, 
#        aes(x = flw_shape_revised, y = Pollinator_richness_percent, fill = Group, shape = Group, color = Group)) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
#   geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
#   scale_fill_manual(values = c("op1better" = "#1B9E77",  
#                                 "op2better" = "#9E9AC8"))+
#   scale_color_manual(values = c("op1better" = "#1B9E77",  
#                                 "op2better" = "#9E9AC8")) +  
#   scale_shape_manual(values = c("op1better" = 16,   # 实心圆
#                                 "op2better" = 17)) + # 实心三角形
#   stat_compare_means(aes(group = Group), 
#                      method = "wilcox.test", 
#                      label = "p.signif", 
#                      hide.ns = TRUE) +
#   labs(
#     x = "Flower Shape",
#     y = "Percent of Pollinator Richness Captured (%)"
#   ) +
#   theme_classic(base_size=12) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "top",
#      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
#     axis.title = element_text(face = "bold", size = 15),
#     axis.text = element_text(color = "black", size = 12),
#     legend.title = element_text(face = "bold", size = 11),
#     legend.text = element_text(size = 10),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     plot.margin = margin(5, 15, 5, 5)
#   )
# 
# pplant<-ggplot(poll_div_by_flw_full, 
#        aes(x = flw_shape_revised, y = Plant_richness_percent, fill = Group, shape = Group, color = Group)) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
#   geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
#   scale_fill_manual(values = c("op1better" = "#1B9E77",  
#                                "op2better" = "#9E9AC8"))+
#   scale_color_manual(values = c("op1better" = "#1B9E77",  
#                                 "op2better" = "#9E9AC8")) +  
#   scale_shape_manual(values = c("op1better" = 16,   # 实心圆
#                                 "op2better" = 17)) + # 实心三角形
#   stat_compare_means(aes(group = Group), 
#                      method = "wilcox.test", 
#                      label = "p.signif", 
#                      hide.ns = TRUE) +
#   labs(
#     #title = "Relative Plant Richness by Flower Shape",
#     x = "Flower Shape",
#     y = "Percent of Plant Richness (%)"
#   ) +
#   theme_classic(base_size=12) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "top",
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
#     axis.title = element_text(face = "bold", size = 15),
#     axis.text = element_text(color = "black", size = 12),
#     legend.title = element_text(face = "bold", size = 11),
#     legend.text = element_text(size = 10),
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     plot.margin = margin(5, 15, 5, 5)
#   )
# final_plot <- plot_grid(ppoll, pplant, labels = c("A", "B"), ncol = 2, align = "h", label_size = 14)
# 
# ggsave("result251105_published/op1-2 explore.png", final_plot, width = 9.8, height = 6, units = "in", dpi = 300)
# 
# 
# ###################
# 
# 
# #############################
# # library(dplyr)
# # library(tidyr)
# # library(vegan)
# # library(ggplot2)
# # library(RColorBrewer)
# # library(cowplot)
# # 
# # # 1. 只保留 Op2better 的数据，并补充花型
# # op2_data <- data_merge %>% 
# #   filter(Study_Network_id %in% Op2better$Study_Network_id) %>%
# #   left_join(traits %>% 
# #               select(Plant_accepted_name, flw_shape_revised),
# #             by = "Plant_accepted_name") %>%
# #   filter(!is.na(flw_shape_revised),
# #          !is.na(Pollinator_accepted_name))
# # 
# # plot_list <- list()
# # adonis_list <- list()
# 
# 
# ###################
# 
# ######################################################################
# 
# ###(2) Is the number of flower shape classes in each site related to the performance of OP2?
# 
# #################################
# 
# ##计算flora shape
# 
# #0visit的,只出现在植物调查之中的植物也应该参与此分析
# data_merge<- merge(data_interact, data_count_scaled_species[,c("Flower_data_merger","Flower_count_scaled","Plant_species")], 
#                    by = "Flower_data_merger",all = TRUE)%>%
#   filter(!is.na(Study_Network_id))%>%
#   mutate(Interaction_addup = ifelse(is.na(Interaction_addup), 0, Interaction_addup))%>% # 替换 `Interaction_addup` 为 NA 的值为 0
#   mutate(
#     Plant_accepted_name = str_replace_all(Plant_accepted_name, "×", "") %>%  # 去掉 ×
#       str_squish()  # 去掉多余空格
#   ) 
# data_merge_trait_with_0 <- data_merge%>%
#   left_join( traits%>%rename(Plant_species = Plant_accepted_name), by = "Plant_species") %>%
#   filter(str_count(Plant_species, " ") >= 1)
# 
# shape_per_site<-data_merge_trait_with_0%>%
#   group_by(Study_Network_id)%>%
#   summarise(flw_shape=n_distinct(flw_shape_revised))
# 
# 
# plot_flw_OP2<-merge(plottest,shape_per_site,by="Study_Network_id")
# # 进行一元线性回归
# model <- lm(flw_shape ~ FlwShape_Top5, data = plot_flw_OP2)
# model_summary <- summary(model)
# p_value <- model_summary$coefficients[2, 4]
# line_type <- ifelse(p_value < 0.05, "solid", "dashed")
# # 绘制散点图和回归线
# scatterplot <- ggplot(plot_flw_OP2, aes(x = flw_shape, y = FlwShape_Top5)) +
#   geom_point(position = position_jitter(width = 0, height = 0.5), color = "#9E9AC8") +
#   geom_smooth(method = "lm", se = FALSE, linetype = line_type, color = "#7B61FF") +
#   scale_x_continuous(breaks = seq(min(plot_flw_OP2$flw_shape),
#                                   max(plot_flw_OP2$flw_shape), by = 1)) +  # 横轴整数
#   labs(x = "Number of Flower Shapes within Networks", y = "Percent of Richness Captured by Subsampling op2") +
#   theme_minimal(base_size = 12) +
#   theme(
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
#     plot.margin = margin(10, 10, 10, 10),
#     strip.background = element_rect(fill = "white", color = NA),
#     strip.text = element_text(size = 12, face = "bold"),
#     axis.ticks.x = element_blank(),
#     legend.position = "right",
#     legend.title = element_text(size = 11),
#     legend.text = element_text(size = 10)
#   )
# 
# print(scatterplot)
# #######
# ##计算flora shape
# plot_flw_OP2<-merge(plot_OP1_OP2,shape_per_site,by="Study_Network_id")
# # 进行一元线性回归
# model2 <- lm(flw_shape ~ op2minusop1, data = plot_flw_OP2)
# model_summary2 <- summary(model2)
# p_value2 <- model_summary2$coefficients[2, 4]
# line_type2 <- ifelse(p_value2 < 0.05, "solid", "dashed")
# # 绘制散点图和回归线
# scatterplot2 <- ggplot(plot_flw_OP2, aes(x = flw_shape, y = op2minusop1)) +
#   geom_point(position = position_jitter(width = 0, height = 0.5), color = "#9E9AC8") +
#   geom_smooth(method = "lm", se = FALSE, linetype = line_type2, color = "#7B61FF") +
#   scale_x_continuous(breaks = seq(min(plot_flw_OP2$flw_shape),
#                                   max(plot_flw_OP2$flw_shape), by = 1)) +  # 横轴整数
#   labs(x = "Number of Flower Shapes within Networks", y = "Δ Percentage of Captured Richness") +
#   theme_minimal(base_size = 12) +
#   theme(
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
#     plot.margin = margin(10, 10, 10, 10),
#     strip.background = element_rect(fill = "white", color = NA),
#     strip.text = element_text(size = 12, face = "bold"),
#     axis.ticks.x = element_blank(),
#     legend.position = "right",
#     legend.title = element_text(size = 11),
#     legend.text = element_text(size = 10)
#   )
# 
# print(scatterplot2)
# ########
# # 进行一元线性回归
# model3 <- lm(flw_shape ~ Abundant_Top5, data = plot_flw_OP2)
# model_summary3 <- summary(model3)
# p_value3 <- model_summary3$coefficients[2, 4]
# line_type3 <- ifelse(p_value < 0.05, "solid", "dashed")
# # 绘制散点图和回归线
# scatterplot3 <- ggplot(plot_flw_OP2, aes(x = flw_shape, y = Abundant_Top5)) +
#   geom_point(position = position_jitter(width = 0, height = 0.5), color = "#9E9AC8") +
#   geom_smooth(method = "lm", se = FALSE, linetype = line_type3, color = "#7B61FF") +
#   scale_x_continuous(breaks = seq(min(plot_flw_OP2$flw_shape),
#                                   max(plot_flw_OP2$flw_shape), by = 1)) +  # 横轴整数
#   labs(x = "Number of Flower Shapes within Networks", y = "Percent of Richness Captured by Subsampling op1") +
#   theme_minimal(base_size = 12) +
#   theme(
#     panel.background = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
#     plot.margin = margin(10, 10, 10, 10),
#     strip.background = element_rect(fill = "white", color = NA),
#     strip.text = element_text(size = 12, face = "bold"),
#     axis.ticks.x = element_blank(),
#     legend.position = "right",
#     legend.title = element_text(size = 11),
#     legend.text = element_text(size = 10)
#   )
# 
# print(scatterplot3)
# 
# ###(3) Overview for flowershape
# #wind sp gets visit
# wind_sp<-data_merge_trait_with_0 %>%
#   filter(flw_shape_revised == "wind flowers") %>%
#   distinct(flw_shape_revised, Plant_accepted_name,Flower_count_scaled)
# 
# unique(wind_sp$Plant_accepted_name)
# unique(wind_sp$Study_Network_id)
# 
# # 计算每种花形状的Study_Network_id比例
# flower_shape_proportion <- data_merge_trait_with_0 %>%
#   filter(!is.na(flw_shape_revised) & flw_shape_revised != "" & flw_shape_revised != "?")%>%
#   group_by(flw_shape_revised) %>%
#   summarise(count = n_distinct(Study_Network_id)) %>%
#   mutate(proportion = count / sum(count)) %>%
#   arrange(desc(proportion))
# 
# # 创建条形图
# barplot<-ggplot(flower_shape_proportion, aes(x = reorder(flw_shape_revised, -proportion), y = proportion)) +
#   geom_bar(stat = "identity", fill = "#9E9AC8", color = "black") +
#   labs(x = "Flower Shapes", y = "Proportion of Networks") +
#   theme_minimal(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转x轴标签，防止重叠+
#         panel.background = element_blank(),     # 去掉灰色背景
#         panel.grid = element_blank(),           # 去掉网格线
#         panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),  # 给每个facet加浅灰边框
#         plot.margin = margin(10, 10, 10, 10),   # 增加整体留白
#         strip.background = element_rect(fill = "white", color = NA), # 去掉标题底色
#         strip.text = element_text(size = 12, face = "bold"),
#         axis.ticks.x = element_blank(),
#         legend.position = "right",
#         legend.title = element_text(size = 11),
#         legend.text = element_text(size = 10)
#   )
# 
# library(cowplot)
# # combine
# final_plot <- plot_grid(
#   barplot, scatterplot2,
#   # labels = c("A", "B"),   # 可选：给子图加标签
#   ncol = 2                # 横向拼接；改成 nrow=2 可以纵向拼接
# )
# final_plot #650*400
# 
# 
# 
# # 保存图像，大小 850*400 像素
# ggsave("E:/Chap1_TargetPlant_to_monitor/result250917_published/flw shape class_capture difference.png", final_plot, width = 8.5, height = 4, units = "in", dpi = 300)
# 
