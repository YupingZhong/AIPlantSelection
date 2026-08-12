
####################################################################################

#=============================RESULT 2 Exploring Op1==================================

# Factors influencing the effectiveness of subsampling

#################################################################################

##-----------(1) hypothesis 1 : plant/pollinator richness, nestedness affect the effectiveness of our strategies (10 plant species not enough)

###
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

# load data
data_count_scaled<-readRDS("data/processed/data_count_scaled_published.rds")
data_interact<-readRDS("data/processed/data_interact_published.rds")
result_all<-read.csv("data/processed/result_all_published_PD.csv")


##------------- Shows solely Abundant-top10-----------

poll_rich<-result_all[,c("Study_Network_id","total_pollinator_count_Abun10")]
# plant data
plant_rich <- data_count_scaled %>%
  filter(Flower_count_scaled != 0) %>%
  group_by(Study_Network_id) %>%
  summarise(
    Plant_Richness = n_distinct(Plant_species, na.rm = TRUE)
  )


# --------------------  Nestedness calculation--------------------

nestedness_df <- data_interact %>%
  filter(!is.na(Interaction_addup)) %>%
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
  
  df <- result_all[, c("Study_Network_id", "percentage_Abun10")] %>%
    left_join(metric_df, by = "Study_Network_id") %>%
    filter(!is.na(percentage_Abun10), !is.na(.data[[metric_name]])) %>%
    mutate(
      Predictor = metric_label
    )
  
  # ---------- Pearson  ----------
  cor_info <- get_cor_info(df, metric_name, "percentage_Abun10")
  lm_info <- get_lm_info(df, metric_name, "percentage_Abun10")
  
  
  label = paste0(
    "R = ", round(cor_info$r, 2),
    "\np = ", format.pval(cor_info$p, digits = 2)
  )
  

  ggplot(
    df,
    aes_string(
      x = metric_name,
      y = "percentage_Abun10",
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
      values=setNames("#E9A77B", metric_label)
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

###############################################
# -------------------- plot--------------------

############

plant_df <- result_all[, c("Study_Network_id", "percentage_Abun10")] %>%
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

pollinator_df <- result_all[, c("Study_Network_id", "percentage_Abun10")] %>%
  left_join(
    poll_rich,
    by="Study_Network_id"
  ) %>%
  filter(
    !is.na(percentage_Abun10),
    !is.na(total_pollinator_count_Abun10)
  ) %>%
  mutate(
    Predictor="Pollinator richness",
    Richness=total_pollinator_count_Abun10
  )

richness_compare <- bind_rows(
  plant_df[,c("Study_Network_id",
              "percentage_Abun10",
              "Predictor",
              "Richness")],
  
  pollinator_df[,c("Study_Network_id",
                   "percentage_Abun10",
                   "Predictor",
                   "Richness")]
)

stats_richness <- richness_compare %>%
  group_by(Predictor) %>%
  summarise(
    r = cor(Richness,
            percentage_Abun10,
            method="pearson"),
    p = summary(
      lm(
        percentage_Abun10 ~ Richness
      )
    )$coefficients[2,4]
  )



p_top10_richness_compare <- ggplot(
  richness_compare,
  aes(
    x = Richness,
    y = percentage_Abun10,
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





# =========================

# Table output

# Abundant_Top10 回归分析

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




#####################################
# # for all strategies



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
# # -------------------- Evenness calculation--------------------
# evenness_df <- data_interact %>%
#   filter(!is.na(Interaction_addup)) %>%
#   group_by(Study_Network_id, Pollinator_accepted_name) %>%
#   summarise(Total_interaction = sum(as.numeric(Interaction_addup)), .groups = "drop") %>%
#   pivot_wider(names_from = Pollinator_accepted_name,
#               values_from = Total_interaction,
#               values_fill = 0) %>%
#   rowwise() %>%
#   mutate(
#     H = diversity(c_across(-Study_Network_id), index = "shannon"),
#     S = specnumber(c_across(-Study_Network_id)),
#     Evenness = H / log(S)
#   ) %>%
#   dplyr::select(Study_Network_id, Evenness)




# # -------------------- plotting (figures for all strategies, all sampling effort) --------------------
# plot_richness_relation <- function(metric_df, metric_name, metric_label, save_path_prefix) {
#   
#   # abundance-based
#   plot_abun <- result_all[, c("Study_Network_id","percentage_Abun10","percentage_Abun5","percentage_Abun3")] %>%
#     left_join(metric_df, by = "Study_Network_id")
#   
#   abun_long <- make_long(plot_abun,
#                          cols = c("percentage_Abun10","percentage_Abun5","percentage_Abun3"),
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

#poll_rich<-result_all[,c("Study_Network_id","total_pollinator_count_Abun10")]

# plot_totalrich <- plot_richness_relation(
#   metric_df = poll_rich,   # 这里是你包含 Total_Rich 的表
#   metric_name = "total_pollinator_count_Abun10",       # 横轴变量
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