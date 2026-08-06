
#######################################################################

#------2.4 Plot main results---------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
result_all<-read.csv("result_all_published_PD.csv",header=TRUE)

#result_all<-read.csv("unic_inter_result_all_published.csv",header=TRUE)#when calculate unique interaction coverage
head(result_all)
plottest<-result_all[,c("Study_Network_id","percentage_Abun10","percentage_Abun5","percentage_Abun3","percentage_FlwShape5","percentage_FlwShape3",
                        "percentage_Pylo10","percentage_Pylo5","percentage_Pylo3","random_mean_percentage_Random10","random_mean_percentage_Random5","random_mean_percentage_Random3")]
n_distinct(result_all$Study_Network_id)
data_long <- pivot_longer(plottest, -Study_Network_id, names_to = "Option", values_to = "Percentage")
data_long_clean <- data_long[complete.cases(data_long), ]
identical(data_long, data_long_clean)

print(data_long_clean)
length(unique(data_long$Study_Network_id))  

# sample size
data_long %>%
  group_by(Option) %>%
  summarise(n_non_missing = sum(!is.na(Percentage)))

#data_long_clean$Percentage <- data_long_clean$Percentage * 100
# summary
data_summary <- data_long_clean %>%
  group_by(Option) %>%
  summarise(
    Mean_Percentage = mean(Percentage, na.rm = TRUE),
    SD_Percentage = sd(Percentage, na.rm = TRUE)
  )
print(data_summary)


custom_order <- c("percentage_Abun10","percentage_Pylo10","random_mean_percentage_Random10",
                  "percentage_Abun5","percentage_FlwShape5","percentage_Pylo5","random_mean_percentage_Random5",
                  "percentage_Abun3","percentage_FlwShape3","percentage_Pylo3","random_mean_percentage_Random3")

data_long_clean$Option <- factor(data_long_clean$Option, levels = custom_order)

non_missing_count <- data_long_clean %>%
  group_by(Option) %>%
  summarise(
    n_non_missing = sum(!is.na(Percentage))
  )

data_summary <- left_join(data_summary, non_missing_count, by = "Option") %>%
  mutate("Focus Species" = case_when(
    Option %in% c("percentage_Abun10",
                  "percentage_Abun5",
                  "percentage_Abun3") ~ "Flower abundance",
    
    Option %in% c("percentage_FlwShape5",
                  "percentage_FlwShape3") ~ "Flower abundance + shapes",
    
    Option %in% c("percentage_Pylo10",
                  "percentage_Pylo5",
                  "percentage_Pylo3") ~ "Phylogenetic distance",
    
    Option %in% c("random_mean_percentage_Random10",
                  "random_mean_percentage_Random5",
                  "random_mean_percentage_Random3") ~ "Random"
  )) %>%
  mutate("Plant Number" = case_when(
    Option %in% c("percentage_Abun10",
                  "percentage_Pylo10",
                  "random_mean_percentage_Random10") ~ "Top10",
    
    Option %in% c("percentage_Abun5",
                  "percentage_FlwShape5",
                  "percentage_Pylo5",
                  "random_mean_percentage_Random5") ~ "Top5",
    
    Option %in% c("percentage_Abun3",
                  "percentage_FlwShape3",
                  "percentage_Pylo3",
                  "random_mean_percentage_Random3") ~ "Top3"
  ))

data_summary$`Plant Number` <- factor(data_summary$`Plant Number`, levels = c("Top10", "Top5", "Top3"))


option_info <- data_summary %>%
  dplyr::select(Option, `Focus Species`, `Plant Number`)

data_long_clean <- data_long_clean %>%
  left_join(option_info, by = "Option")
data_long_clean$`Plant Number` <- factor(data_long_clean$`Plant Number`, levels = c("Top10", "Top5", "Top3"))


focus_species_colors <- c(
  "Flower abundance" = "#1B9E77",
  "Flower abundance + shapes" = "#9E9AC8",
  "Phylogenetic distance" = "#E8A3D1",
  "Random" = "grey70"
)

# 设置 factor 顺序
data_summary$`Focus Species` <- factor(data_summary$`Focus Species`,
                                       levels = c(
                                         "Flower abundance",
                                         "Flower abundance + shapes",
                                         "Phylogenetic distance",
                                         "Random"
                                       ))

data_long_clean$`Focus Species` <- factor(data_long_clean$`Focus Species`,
                                          levels = c(
                                            "Flower abundance",
                                            "Flower abundance + shapes",
                                            "Phylogenetic distance",
                                            "Random"
                                          ))

set.seed(2025)

plot <- ggplot(data_summary, aes(x = Option, y = Mean_Percentage, color = `Focus Species`)) +
  geom_jitter(data = data_long_clean,
              aes(x = Option, y = Percentage, color = `Focus Species`),
              width = 0.2, size = 1, alpha = 0.5) +
  geom_crossbar(aes(ymin = Mean_Percentage - SD_Percentage, ymax = SD_Percentage + Mean_Percentage),
                position = position_dodge(width = 0.3), width = 0.5, color = "grey40") +
  geom_text(aes(label = paste("n =", n_non_missing)),
            position = position_dodge(width = 0.5), vjust = -1, size = 4, color = "black") +
  geom_text(aes(label = paste(round(Mean_Percentage, 2))),
            position = position_dodge(width = 0.5), vjust = 1.5, size = 4, color = "black") +
  facet_wrap(~ `Plant Number`, nrow = 1, scales = "free_x") +
  scale_color_manual(values = focus_species_colors) +
  labs(x = "Subsampling strategy", y = "Percent of pollinator richness captured", color = "Plants selected based on") + #"Percent of Unique Interaction Captured"
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", size = 13),
    panel.background = element_blank(),    
    panel.grid = element_blank(),           
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),  
    plot.margin = margin(10, 10, 10, 10),   
    strip.background = element_rect(fill = "white", color = NA), 
    strip.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal", 
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )+ 
  guides(color = guide_legend(override.aes = list(size = 5)))


plot
#############
# ------------------------------------------
# Test for significant differences between groups
# ------------------------------------------

library(tidyverse)
library(car)
library(ggsignif)


# extract groups automatically
get_group <- function(option_name){
  data_long_clean %>%
    filter(Option == option_name) %>%
    pull(Percentage)
}


group_abun10 <- get_group("percentage_Abun10")
group_abun5  <- get_group("percentage_Abun5")
group_abun3  <- get_group("percentage_Abun3")

group_shape5 <- get_group("percentage_FlwShape5")
group_shape3 <- get_group("percentage_FlwShape3")

group_pylo10 <- get_group("percentage_Pylo10")
group_pylo5  <- get_group("percentage_Pylo5")
group_pylo3  <- get_group("percentage_Pylo3")

group_random10 <- get_group("random_mean_percentage_Random10")
group_random5  <- get_group("random_mean_percentage_Random5")
group_random3  <- get_group("random_mean_percentage_Random3")


# ------------------------------------------
# Levene variance test
# ------------------------------------------

compare_variance <- function(g1, g2, group_names){
  
  test_data <- data.frame(
    Percentage = c(g1, g2),
    Group = factor(c(
      rep(group_names[1], length(g1)),
      rep(group_names[2], length(g2))
    ))
  )
  
  levene_result <- car::leveneTest(
    Percentage ~ Group,
    data=test_data
  )
  
  pval <- levene_result[1,"Pr(>F)"]
  
  cat(
    group_names[1],
    "vs",
    group_names[2],
    "Levene p =",
    pval,
    "\n"
  )
  
  return(pval > 0.05)
}



# ------------------------------------------
# Statistical test selection
# ------------------------------------------

compare_groups <- function(g1,g2,group_names){
  
  g1 <- na.omit(g1)
  g2 <- na.omit(g2)
  
  
  normal1 <- shapiro.test(g1)$p.value > 0.05
  normal2 <- shapiro.test(g2)$p.value > 0.05
  
  equal_var <- compare_variance(
    g1,
    g2,
    group_names
  )
  
  
  if(normal1 & normal2 & equal_var){
    
    result <- t.test(
      g1,
      g2,
      var.equal=TRUE
    )
    
    method <- "T-test"
    
  }else{
    
    result <- wilcox.test(
      g1,
      g2,
      exact=FALSE
    )
    
    method <- "Wilcoxon"
    
  }
  
  
  data.frame(
    Comparison=paste(group_names, collapse=" vs "),
    Method=method,
    P_Value=result$p.value,
    Significance=case_when(
      result$p.value <0.001 ~ "***",
      result$p.value <0.01  ~ "**",
      result$p.value <0.05  ~ "*",
      TRUE ~ "ns"
    )
  )
}



# ------------------------------------------
# Perform comparisons
# ------------------------------------------

test_results <- bind_rows(
  
  # flower shape vs abundance
  compare_groups(
    group_shape5,
    group_abun5,
    c("percentage_FlwShape5",
      "percentage_Abun5")
  ),
  
  compare_groups(
    group_shape3,
    group_abun3,
    c("percentage_FlwShape3",
      "percentage_Abun3")
  ),
  
  
  # phylogenetic vs abundance
  compare_groups(
    group_pylo10,
    group_abun10,
    c("percentage_Pylo10",
      "percentage_Abun10")
  ),
  
  compare_groups(
    group_pylo5,
    group_abun5,
    c("percentage_Pylo5",
      "percentage_Abun5")
  ),
  
  compare_groups(
    group_pylo3,
    group_abun3,
    c("percentage_Pylo3",
      "percentage_Abun3")
  ),
  
  
  # abundance vs random
  compare_groups(
    group_abun10,
    group_random10,
    c("percentage_Abun10",
      "random_mean_percentage_Random10")
  ),
  
  compare_groups(
    group_abun5,
    group_random5,
    c("percentage_Abun5",
      "random_mean_percentage_Random5")
  ),
  
  compare_groups(
    group_abun3,
    group_random3,
    c("percentage_Abun3",
      "random_mean_percentage_Random3")
  )
  
)


print(test_results)



# ------------------------------------------
# Prepare significance annotations
# ------------------------------------------

significance_results <- test_results %>%
  mutate(
    
    PlantNumber = case_when(
      grepl("10", Comparison) ~ "Top10",
      grepl("5", Comparison) ~ "Top5",
      grepl("3", Comparison) ~ "Top3"
    ),
    
    group1=sub(" vs .*","",Comparison),
    
    group2=sub(".* vs ","",Comparison),
    
    p_value=P_Value,
    
    significance=Significance,
    
    y_position=case_when(
      PlantNumber=="Top10" ~ 105,
      PlantNumber=="Top5" ~ 100,
      PlantNumber=="Top3" ~ 95
    )
    
  ) %>%
  select(
    group1,
    group2,
    p_value,
    significance,
    PlantNumber,
    y_position
  )


significance_results$PlantNumber <- factor(
  significance_results$PlantNumber,
  levels=c("Top10","Top5","Top3")
)


significance_results <- significance_results %>%
  group_by(PlantNumber) %>%
  mutate(
    y_position=y_position + row_number()*5
  ) %>%
  ungroup()


significance_results
# ------------------------------------------
# Add significance annotations to plot (800x450)
# ------------------------------------------
set.seed(2025)

for(i in seq_len(nrow(significance_results))) {
  res <- significance_results[i, ]
  filtered_data <- data_long_clean %>% filter(`Plant Number` == as.character(res$PlantNumber))
  
  if(nrow(filtered_data) > 0) {
    sig_color <- ifelse(res$significance == "ns", "darkgray", "black")
    
    plot <- plot + 
      geom_signif(
        comparisons = list(c(as.character(res$group1), as.character(res$group2))),
        annotations = res$significance,
        y_position = res$y_position,
        tip_length = 0.02,
        color = sig_color,
        textsize = 4,
        data = filtered_data,
        inherit.aes = FALSE,
        aes(x = Option, y = Percentage)
      )
  }
}



print(plot)
##图大小800x450


ggsave("/Chap1_TargetPlant_to_monitor/result_260723/main_phylo.png", plot, width = 9, height = 6, units = "in", dpi = 300)
ggsave("/Chap1_TargetPlant_to_monitor/result_260723/main.png", plot, width = 7, height = 5, units = "in", dpi = 300)
ggsave("/Chap1_TargetPlant_to_monitor/result_260526/main_unique_interaction.png", plot, width = 7, height = 5, units = "in", dpi = 300)


#############
#---z-score---
#############

z_results<-readRDS("z_results.rds")
sample_size <- z_results %>%
  filter(!is.na(Z_score)) %>%
  group_by(
    Plant_number,
    Strategy
  ) %>%
  summarise(
    n = n(),
    y = max(Z_score) + 0.3,
    .groups="drop"
  )

z_results %>%
  group_by(
    Strategy,
    Plant_number
  ) %>%
  summarise(
    p_value = wilcox.test(
      Z_score,
      mu = 0
    )$p.value
  )

z_sig <- z_results %>%
  filter(!is.na(Z_score)) %>%
  group_by(
    Strategy,
    Plant_number
  ) %>%
  summarise(
    p_value =
      wilcox.test(
        Z_score,
        mu = 0
      )$p.value,
    .groups="drop"
  ) %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

y_pos <- z_results %>%
  filter(!is.na(Z_score)) %>%
  group_by(
    Strategy,
    Plant_number
  ) %>%
  summarise(
    y=max(Z_score)+0.3,
    .groups="drop"
  )

z_sig <- left_join(
  z_sig,
  y_pos,
  by=c(
    "Strategy",
    "Plant_number"
  )
)

zplot<-ggplot(
  z_results %>% filter(!is.na(Z_score)),
  aes(
    x = Plant_number,
    y = Z_score,
    fill = Strategy
  )
) +
  geom_violin(
    alpha=0.4
  ) +
  geom_boxplot(
    width=0.15,
    position=position_dodge(0.9)
  ) +
  geom_jitter(
    aes(color=Strategy),
    position=position_jitterdodge(
      jitter.width=0.15,
      dodge.width=0.9
    ),
    alpha=0.2,
    size=0.8
  ) +
  geom_text(
    data = sample_size,
    aes(
      y = y,
      label = paste0("n=", n)
    ),
    position = position_dodge(width = 0.9),
    size = 3
  ) +
  geom_text(
    data=z_sig,
    aes(
      x=Plant_number,
      y=y+0.3,
      label=label,
      color=Strategy
    ),
    position=position_dodge(0.9),
    size=5,
    inherit.aes=FALSE
  )+
  scale_fill_manual(
    values=focus_species_colors
  ) +
  scale_color_manual(
    values=focus_species_colors
  ) +
  geom_hline(
    yintercept=0,
    linetype="dashed"
  ) +
  theme_classic()

zplot

ggsave("/Chap1_TargetPlant_to_monitor/result_260723/zscore_phylo.png", zplot, width = 9, height = 6, units = "in", dpi = 300)


friedman.test(
  Z_score ~ Strategy | Study_Network_id,
  data = top10_data
)

############################################ 
# ------------------------------------------
#  Generate a supplementary statistical table
# ------------------------------------------

# 定义辅助函数：根据两组数据生成表格行
make_row <- function(g1, g2, label1, label2){
  
  test <- wilcox.test(g1, g2)
  
  data.frame(
    Comparison = paste(label1, "vs", label2),
    n = sum(!is.na(g1)),  # 两组一样的话用一个即可
    Median1 = round(median(g1, na.rm = TRUE), 2),
    Median2 = round(median(g2, na.rm = TRUE), 2),
    P_raw = test$p.value
  )
}



# ------------------------------------------
#  within & between strategies
# ------------------------------------------
supp_table <- bind_rows(
  
  # Within strategy
  make_row(group_abundant_top10, group_abundant_top5, "Top 10 abundant", "Top 5 abundant"),
  make_row(group_abundant_top5, group_abundant_top3, "Top 5 abundant", "Top 3 abundant"),
  make_row(group_flwshape_top5, group_flwshape_top3, "Top 5 flower-shape", "Top 3 flower-shape"),
  make_row(group_random_10, group_random_5, "Random 10", "Random 5"),
  make_row(group_random_5, group_random_3, "Random 5", "Random 3"),
  
  # Between strategy
  make_row(group_abundant_top10, group_random_10, "Top 10 abundant", "Random 10"),
  make_row(group_abundant_top5, group_random_5, "Top 5 abundant", "Random 5"),
  make_row(group_abundant_top3, group_random_3, "Top 3 abundant", "Random 3"),
  make_row(group_flwshape_top5, group_random_5, "Top 5 flower-shape", "Random 5"),
  make_row(group_flwshape_top3, group_random_3, "Top 3 flower-shape", "Random 3"),
  make_row(group_abundant_top5, group_flwshape_top5, "Top 5 abundant", "Top 5 flower-shape"),
  make_row(group_abundant_top3, group_flwshape_top3, "Top 3 abundant", "Top 3 flower-shape")
)

supp_table <- supp_table %>%
  mutate(
    P_adjusted = p.adjust(P_raw, method = "BH"),
    
    # 格式化 P 值（论文写法）
    P_value = case_when(
      P_adjusted < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", P_adjusted)
    )
  ) %>%
  
  dplyr::select(Comparison, n, Median1, Median2, P_value)


#save
write.csv(
  supp_table,
  "./result_260526/Supplementary_Statistical_Table_AllComparisons.csv",
  row.names = FALSE
)


# ------------------------------------------
# generate word form
# ------------------------------------------
library(flextable)
library(officer)

ft <- flextable(supp_table)

ft <- theme_booktabs(ft)
ft <- autofit(ft)
ft <- fontsize(ft, size = 10, part = "all")
ft <- align(ft, align = "center", part = "all")
ft <- align(ft, j = "Comparison", align = "left", part = "all")
ft <- bold(ft, part = "header")

doc <- read_docx() %>%
  body_add_par(
    "Table S1. Pairwise comparisons of pollinator richness captured under different plant selection strategies. Values represent medians. P-values were adjusted using the Benjamini–Hochberg method.",
    style = "heading 2"
  ) %>%
  body_add_flextable(ft)

print(doc, target = "./result_260526/Supplementary_Table_Subsampling_AllComparisons.docx")


####海报绘图
# plot <- plot +
#   theme(
#     # 设置所有文字为白色
#     text = element_text(color = "white"),
#     axis.title = element_text(color = "white"),
#     axis.text = element_text(color = "white"),
#     strip.text = element_text(color = "white"),
#     legend.title = element_text(color = "white"),
#     legend.text = element_text(color = "white"),
#     
#     # 去掉背景使透明
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA),
#     panel.border = element_rect(color = "white", fill = NA, linewidth = 0.5),
#     
#     # 坐标网格线（如果想要也可用白色）
#     panel.grid.major = element_line(color = "grey30"),
#     panel.grid.minor = element_blank()
#   )
# 
# ggsave("result251105_published/main_poster.png", plot = plot,
#        bg = "transparent", width = 10, height = 5, dpi = 300)


############################################################

##----3.5 plot Relationship of the richness of pollinators of subsampled networks ---------------
#         with the richness of the full network for each subsampling strategy


library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
head(result_all)
result_all<-read.csv("result_all_published_phylo.csv",header=TRUE)


plot_rich_site<-result_all[,c("Study_Network_id","Abun_Top10_Rich","Abun_Top3_Rich","Abun_Top5_Rich",
                              "FlwShape_Top5_Rich","FlwShape_Top3_Rich","Total_Rich")]
head(plot_rich_site)

###################################################################

# ---------- 基础主题和颜色 ----------
my_colors <- c(
  "Top 10" = "#66C2A5", "Top 5" = "#4E79A7", "Top 3" = "#F28E2B",
  "flw 5" = "#4E79A7", "flw 3" = "#F28E2B",
  "Random 10" = "#66C2A5", "Random 5" = "#4E79A7", "Random 3" = "#F28E2B"
)

base_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.position = "none"
  )

# ---------- 函数：计算每个分组回归显著性 ----------
get_linetype <- function(df, x_col, y_col, group_col){
  
  df %>%
    group_by(across(all_of(group_col))) %>%
    summarise(
      lm_res = list(
        lm(
          as.formula(paste(y_col, "~", x_col)),
          data = cur_data()
        )
      ),
      .groups = "drop"
    ) %>%
    mutate(
      p = sapply(lm_res, \(x) summary(x)$coefficients[2,4]),
      linetype = ifelse(p < 0.05, "solid", "dashed")
    ) %>%
    dplyr::select(all_of(group_col), linetype)
}




# ---------- reshape Abundance, Flower-shape, Random ----------
rich_abun_long <- result_all %>%
  dplyr::select(Study_Network_id, Total_Rich, Abun_Top10_Rich, Abun_Top5_Rich, Abun_Top3_Rich) %>%
  pivot_longer(cols = c(Abun_Top10_Rich, Abun_Top5_Rich, Abun_Top3_Rich),
               names_to = "Top_Plant_Level", values_to = "Richness_Value") %>%
  mutate(Top_Plant_Level = factor(Top_Plant_Level, 
                                  levels = c("Abun_Top10_Rich","Abun_Top5_Rich","Abun_Top3_Rich"),
                                  labels = c("Top 10","Top 5","Top 3")))

rich_shape_long <- result_all %>%
  dplyr::select(Study_Network_id, Total_Rich, FlwShape_Top5_Rich, FlwShape_Top3_Rich) %>%
  pivot_longer(cols = c(FlwShape_Top5_Rich, FlwShape_Top3_Rich),
               names_to = "flw_Plant_Level", values_to = "Richness_Value") %>%
  mutate(Top_Plant_Level = factor(flw_Plant_Level,
                                  levels = c("FlwShape_Top5_Rich","FlwShape_Top3_Rich"),
                                  labels = c("flw 5","flw 3")))

rich_random_long <- result_all %>%
  dplyr::select(Study_Network_id, Total_Rich, Random_10_Rich, Random_5_Rich, Random_3_Rich) %>%
  pivot_longer(cols = c(Random_10_Rich, Random_5_Rich, Random_3_Rich),
               names_to = "Random_Level", values_to = "Richness_Value") %>%
  mutate(Top_Plant_Level = factor(Random_Level,
                                  levels = c("Random_10_Rich","Random_5_Rich","Random_3_Rich"),
                                  labels = c("Random 10","Random 5","Random 3")))

# ---------- 计算每组 linetype ----------
linetype_abun <- get_linetype(
  rich_abun_long,
  x_col = "Total_Rich",
  y_col = "Richness_Value",
  group_col = "Top_Plant_Level"
)

linetype_shape <- get_linetype(
  rich_shape_long,
  x_col = "Total_Rich",
  y_col = "Richness_Value",
  group_col = "Top_Plant_Level"
)

linetype_random <- get_linetype(
  rich_random_long,
  x_col = "Total_Rich",
  y_col = "Richness_Value",
  group_col = "Top_Plant_Level"
)



# ---------- 绘图函数，增加 stat_cor ----------
plot_rich <- function(df, linetype_df, ylab,
                      cor_label_x = "left",
                      cor_label_y = NULL) {
  
  n_group <- length(unique(df$Top_Plant_Level))
  
  if (is.null(cor_label_y)) {
    cor_label_y <- seq(0.98, 0.98 - 0.02*(n_group-1), length.out = n_group)
  }
  
  ggplot(df, aes(x = Total_Rich, y = Richness_Value, color = Top_Plant_Level)) +
    geom_point(shape = 1, size = 2,
               position = position_jitter(width = 0.03, height = 0.5)) +
    geom_smooth(
      method = "lm", se = FALSE,
      aes(linetype = Top_Plant_Level),
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    scale_linetype_manual(values = setNames(linetype_df$linetype, linetype_df[[1]])) +
    stat_cor(
      aes(group = Top_Plant_Level, color = Top_Plant_Level),
      method = "pearson", size = 4,
      label.x.npc = cor_label_x,
      label.y.npc = cor_label_y
    ) +
    scale_color_manual(values = my_colors) +
    labs(x = "Pollinator richness in full network", y = ylab) +
    base_theme+
    theme(
      plot.margin = margin(t = 0, r = 10, b = 0, l = 0)
    )
}


# ---------- 绘制三张图 ----------
library(cowplot)
library(ggplot2)

p_abun <- plot_rich(rich_abun_long, linetype_abun, ylab = "") +
  ggtitle("A. Abundance-based") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.text = element_text(size = 11)
  )+
  scale_y_continuous(limits = c(0, 150))


p_shape <- plot_rich(rich_shape_long, linetype_shape, ylab = "") +
  ggtitle("B. Abundance + flower shape") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.text = element_text(size = 11)
  )+
  scale_y_continuous(limits = c(0, 150))

p_random <- plot_rich(rich_random_long, linetype_random, ylab = "", cor_label_x = 0.33) +
  ggtitle("C. Random") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.text = element_text(size = 11)
  )+
  scale_y_continuous(limits = c(0, 150))

# ---------- 单独生成 legend ----------
legend_plot <- p_abun +
  labs(colour = "Number of plant species subsampled") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

legend <- cowplot::get_legend(legend_plot)


# ---------- 拼接三张图 + legend ----------
plots_row <- plot_grid(
  p_abun, p_shape, p_random,
  ncol = 3,
  align = "h",    # 水平对齐
  axis = "tb",    # 上下轴对齐
  rel_widths = c(1, 1, 1) # 三张图宽度相等
)

final_plot <- plot_grid(
  plots_row,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

# ---------- 添加大 Y 轴标题 ----------
final_plot <- ggdraw(final_plot) +
  draw_label(
    "Pollinator richness in subsampled network",
    x = 0, y = 0.5,          # 左侧居中纵向
    angle = 90,              # 竖直
    fontface = "bold",
    size = 12,
    hjust = 0.41,             # 完全水平居中
    vjust = 1.5              # 完全纵向居中
  )

final_plot


ggsave("/Chap1_TargetPlant_to_monitor/result_260526/rich_3.png", final_plot, width = 9.8, height = 5, units = "in", dpi = 300)
##==============================================================================

# --- 3. Network Metrics ------------------------------------------

##==============================================================================

#---3.1 calculate NODF for each networks ------------

#####
library(bipartite)
library(vegan)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)

#  Nestedness function
calculate_nestedness <- function(df) {
  
  dt <- as.data.table(df)
  
  interaction_matrix <- dcast(
    dt,
    Plant_accepted_name ~ Pollinator_accepted_name,
    value.var = "Interaction_addup",
    fill = 0
  )
  
  mat <- as.matrix(interaction_matrix[, -1, with = FALSE])
  mat[is.na(mat)] <- 0
  
  # 去掉空行列
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NA_real_)
  
  tryCatch(
    bipartite::nested(mat, method = "NODF"),
    error = function(e) NA_real_
  )
}

################
# ID skeleton

all_ids <- data.frame(
  Study_Network_id = sort(unique(result_orig$Study_Network_id))
)

#############
#  nestedness function

nested_df_from_data <- function(data, colname, id_skeleton) {
  
  df <- data %>%
    filter(
      !is.na(Pollinator_accepted_name),
      !is.na(Plant_accepted_name),
      !is.na(Interaction_addup)
    ) %>%
    dplyr::select(
      Plant_accepted_name,
      Pollinator_accepted_name,
      Interaction_addup,
      Study_Network_id
    ) %>%
    group_by(Study_Network_id) %>%
    summarise(
      value = calculate_nestedness(cur_data()),
      .groups = "drop"
    )
  
  # align skeleton（关键步骤）
  df_aligned <- id_skeleton %>%
    left_join(df, by = "Study_Network_id") %>%
    rename(!!colname := value)
  
  return(df_aligned)
}


datasets <- list(
  orig = result_orig,
  top10 = result10,
  top5 = result5,
  top3 = result3,
  op2_top5 = abun_species_top5_merge,
  op2_top3 = abun_species_top3_merge
)



nested_list <- imap(
  datasets,
  ~ nested_df_from_data(.x, .y, all_ids)
)

all_nestedness_df <- reduce(nested_list, left_join, by = "Study_Network_id") %>%
  rename(
    Nestedness_10 = top10,
    Nestedness_5 = top5,
    Nestedness_3 = top3,
    Nestedness_op2_top5 = op2_top5,
    Nestedness_op2_top3 = op2_top3,
    Nestedness_orig = orig
  )%>%
  mutate(across(
    starts_with("Nestedness"),
    ~ as.numeric(unname(.))
  ))

############################################################
## check if the ID match

id_check <- sapply(nested_list, function(df)
  identical(df$Study_Network_id, all_ids$Study_Network_id)
)

if(!all(id_check)){
  stop("❌ ID mismatch detected! Check alignment.")
} else {
  cat("✅ All datasets perfectly aligned with identical Study_Network_id\n")
}

# ✔ row number
cat("Total networks:", nrow(all_nestedness_df), "\n")

# ✔ repeat ID ?
if(any(duplicated(all_nestedness_df$Study_Network_id))){
  stop("❌ Duplicate Study_Network_id found!")
} else {
  cat("✅ No duplicated IDs\n")
}

###################
# check  NA

na_summary <- colSums(is.na(all_nestedness_df))
print(na_summary)

###################
#  Spot check

unique(all_nestedness_df$Study_Network_id)
net1<-result3 %>%
  filter(Study_Network_id == "6_Marini_UNIPD03 CD20_2020") %>%
  calculate_nestedness()
net1
net_check<-all_nestedness_df %>%
  filter(Study_Network_id == "6_Marini_UNIPD03 CD20_2020")
net_check
#抽个检查发现没问题

##################
# =========================
# Pearson correlations
# =========================

nodf_cols <- c("Nestedness_orig", "Nestedness_10", "Nestedness_5", "Nestedness_3",
               "Nestedness_op2_top5", "Nestedness_op2_top3")



cor_results <- lapply(nodf_cols[-1], function(x) {
  
  x_vec <- as.numeric(all_nestedness_df[[x]])
  y_vec <- as.numeric(all_nestedness_df$Nestedness_orig)
  
  # exclude NA
  valid_idx <- !is.na(x_vec) & !is.na(y_vec)
  
  x_clean <- x_vec[valid_idx]
  y_clean <- y_vec[valid_idx]
  
  if(length(x_clean) < 3 || length(unique(x_clean)) < 2){
    return(c(r = NA, p_value = NA, n = length(x_clean)))
  }
  
  test <- cor.test(x_clean, y_clean, method = "pearson")
  
  c(r = unname(test$estimate), p_value = test$p.value, n = length(x_clean))
})


# dataframe
cor_results_df <- as.data.frame(do.call(rbind, cor_results)) %>%
  mutate(
    Comparison = nodf_cols[-1],
    Metric = "Nestedness",
    r = as.numeric(r),
    p_value = as.numeric(p_value)
  )

rownames(cor_results_df) <- nodf_cols[-1]
# 
# print("Pearson correlation with full network:")
print(cor_results_df)

cor(all_nestedness_df$Nestedness_op2_top5,
    all_nestedness_df$Nestedness_orig,
    use = "complete.obs")

#############################################################

########  Visualization: 

# Noted: “R” and “p” indicate the Pearson correlation coefficient and its corresponding p-value. 
#          Lines represent linear model fits for visualization purposes. 

# Define color palette for different top-N plant levels
my_colors <- c(
  "Top 10" = "#66C2A5",  # Mint green-blue
  "Top 5"  = "#4E79A7",  # Calm blue
  "Top 3"  = "#F28E2B",  # Soft orange
  "flw 5"  = "#4E79A7",  # Same as Top 5
  "flw 3"  = "#F28E2B"   # Same as Top 3
)
# Convert nestedness data to long format for Top N plants
long_df <- all_nestedness_df %>%
  pivot_longer(
    cols = c(Nestedness_10, Nestedness_5, Nestedness_3),
    names_to = "Top_Plant_Level",
    values_to = "Nestedness_Value"
  ) %>%
  filter(!is.na(Nestedness_Value)) %>%
  mutate(
    Top_Plant_Level = factor(
      Top_Plant_Level,
      levels = c("Nestedness_10", "Nestedness_5", "Nestedness_3"),
      labels = c("Top 10", "Top 5", "Top 3")
    )
  )


# Plot nestedness comparisons for Top N plants vs full network
reg_significance <- long_df %>%
  group_by(Top_Plant_Level) %>%
  summarise(
    lm_res = list(lm(Nestedness_Value ~ Nestedness_orig, data = cur_data_all())),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = summary(lm_res)$coefficients[2,4],
    signif_line = ifelse(p_value < 0.05, "solid", "dashed")
  ) %>%
  dplyr::select(Top_Plant_Level, signif_line)

set.seed(2025)
p1 <- ggplot(long_df, aes(y = Nestedness_Value, x = Nestedness_orig, color = Top_Plant_Level)) +
  geom_point(shape = 1, size = 2, position = position_jitter(width = 0.03, height = 0.5)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    aes(linetype = Top_Plant_Level),  # 根据显著性映射线型
    linewidth = 0.8
  ) +
  scale_linetype_manual(values = setNames(reg_significance$signif_line, reg_significance$Top_Plant_Level)) +
  scale_y_continuous(limits = c(0, 60)) +
  stat_cor(aes(color = Top_Plant_Level), method = "pearson", 
           size = 4, 
           label.x.npc = "left",
           label.y.npc = c(0.98, 0.96,0.92)) +
  scale_color_manual(values = my_colors) +
  labs(
    y = "NODF in Subnetwork of Top-N Plants",
    x = "NODF of the Full Network",
    color = "Top N Level",
    linetype = "Regression significance"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )

# Reshape data for flower abundance top 5 and 3
long_df2 <- all_nestedness_df %>%
  pivot_longer(
    cols = c(Nestedness_op2_top5, Nestedness_op2_top3),
    names_to = "flw_Plant_Level",
    values_to = "Nestedness_Value"
  ) %>%
  filter(!is.na(Nestedness_Value)) %>%
  mutate(
    Top_Plant_Level = factor(
      flw_Plant_Level,
      levels = c("Nestedness_op2_top5", "Nestedness_op2_top3"),
      labels = c("flw 5", "flw 3")
    )
  )

reg_significance2 <- long_df2 %>%
  group_by(Top_Plant_Level) %>%
  summarise(
    lm_res = list(lm(Nestedness_Value ~ Nestedness_orig, data = cur_data_all())),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = summary(lm_res)$coefficients[2,4],
    signif_line = ifelse(p_value < 0.05, "solid", "dashed")
  ) %>%
  dplyr::select(Top_Plant_Level, signif_line)


# Plot nestedness for flower abundance based subnetworks
p2 <- ggplot(long_df2, aes(y = Nestedness_Value, x = Nestedness_orig, color = Top_Plant_Level)) +
  geom_point(shape = 1, size = 2, position = position_jitter(width = 0.03, height = 0.5)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    aes(linetype = Top_Plant_Level),
    linewidth = 0.8
  ) +
  scale_y_continuous(limits = c(0, 60)) +
  scale_linetype_manual(values = setNames(reg_significance2$signif_line, reg_significance2$Top_Plant_Level)) +
  stat_cor(aes(color = Top_Plant_Level), method = "pearson",
           size = 4, 
           label.x.npc = "left",
           label.y.npc = c(0.98, 0.96)) +
  scale_color_manual(values = my_colors) +
  labs(
    y = "NODF in Subnetwork of Flw-N Plants",
    x = "NODF of the Full Network",
    color = "Top N Level",
    linetype = "Regression significance"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )


# Display plots side by side
library(cowplot)
plot_grid(p1, p2, ncol = 2)

#########################################

#################################################################################

###-----3.2 Calculate connectance for each network -----------


library(dplyr)
library(purrr)
library(data.table)

# -----------------------------
# global ID reference
# -----------------------------
all_ids <- sort(unique(data_merge$Study_Network_id))

# -----------------------------
# 1️⃣ connectance function
# -----------------------------
calculate_connectance <- function(df) {
  
  if (is.null(df) || nrow(df) == 0) return(NA_real_)
  
  dt <- as.data.table(df)
  
  mat_df <- dcast(
    dt,
    Plant_accepted_name ~ Pollinator_accepted_name,
    value.var = "Interaction_addup",
    fill = 0
  )
  
  mat <- as.matrix(mat_df[, -1, with = FALSE])
  mat[is.na(mat)] <- 0
  
  # remove empty rows/cols
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NA_real_)
  
  sum(mat > 0) / (nrow(mat) * ncol(mat))
}

# -----------------------------

connectance_df_from_data <- function(data, colname, all_ids) {
  
  out <- data %>%
    dplyr::filter(
      !is.na(Pollinator_accepted_name),
      !is.na(Plant_accepted_name),
      !is.na(Interaction_addup)
    ) %>%
    dplyr::select(
      Plant_accepted_name,
      Pollinator_accepted_name,
      Interaction_addup,
      Study_Network_id
    ) %>%
    dplyr::group_by(Study_Network_id) %>%
    dplyr::group_modify(~{
      dplyr::tibble(
        value = calculate_connectance(.x)
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(!!colname := value)
  
  # align ID
  out <- tibble(Study_Network_id = all_ids) %>%
    left_join(out, by = "Study_Network_id")
  
  return(out)
}

# -----------------------------
# datasets
# -----------------------------
datasets <- list(
  orig = data_merge,
  top10 = result10,
  top5 = result5,
  top3 = result3,
  op2_top5 = abun_species_top5_merge,
  op2_top3 = abun_species_top3_merge
)

connect_list <- purrr::imap(
  datasets,
  ~ connectance_df_from_data(.x, .y, all_ids)
)

# -----------------------------
# strict merge
# -----------------------------
all_connectance_df <- purrr::reduce(
  connect_list,
  dplyr::full_join,
  by = "Study_Network_id"
)
all_connectance_df <- all_connectance_df %>%
  dplyr::rename(
    Connectance_orig = orig,
    Connectance_10 = top10,
    Connectance_5 = top5,
    Connectance_3 = top3,
    Connectance_op2_top5 = op2_top5,
    Connectance_op2_top3 = op2_top3
  )
# -----------------------------
#  ID strict validation
# -----------------------------

# ✔ 1. row number
if (nrow(all_connectance_df) != length(all_ids)) {
  stop("❌ Row number mismatch!")
}

# ✔ 2. ID
if (!setequal(all_connectance_df$Study_Network_id, all_ids)) {
  stop("❌ ID set mismatch!")
}

# ✔ 3. order
all_connectance_df <- all_connectance_df %>%
  dplyr::arrange(match(Study_Network_id, all_ids))

cat("✅ All connectance datasets aligned correctly\n")


all_connectance_df <- all_connectance_df %>%
  dplyr::mutate(
    dplyr::across(-Study_Network_id, ~ as.numeric(.))
  )


#⃣ quick sanity check

head(all_connectance_df)
summary(all_connectance_df)

##################
check_connectance <- function(id, data, table, col){
  
  raw <- data %>%
    filter(Study_Network_id == id) %>%
    calculate_connectance()
  
  stored <- table %>%
    filter(Study_Network_id == id) %>%
    pull({{col}})
  
  data.frame(
    id = id,
    raw = raw,
    stored = stored,
    diff = raw - stored
  )
}

unique(all_connectance_df$Study_Network_id)
check_connectance(
  "38_Maurer_R1",
  result10,
  all_connectance_df,
  Connectance_10
)
#----- Pearson correlation---------------

library(ggplot2)
library(dplyr)

na_counts <- sapply(all_connectance_df, function(x) sum(is.na(x)))
print("NA counts per column:")
print(na_counts)

connect_cols <- c("Connectance_orig", "Connectance_10", "Connectance_5", "Connectance_3",
                  "Connectance_op2_top5", "Connectance_op2_top3")

par(mfrow = c(2,3))  
for(col in connect_cols) {
  
  x <- all_connectance_df[[col]]
  
  x <- as.numeric(unlist(x))  
  
  x <- x[is.finite(x)]
  
  hist(x,
       breaks = 30,
       main = paste("Distribution of", col),
       xlab = "Connectance",
       col = "skyblue")
}


#  子网络 vs 全网络散点图
for(col in connect_cols[-1]) {
  plot(all_connectance_df[[col]], all_connectance_df$Connectance_orig,
       xlab = paste(col, "Connectance"),
       ylab = "Full Network Connectance",
       main = paste("Comparison:", col, "vs Full Network"),
       pch = 19, col = "steelblue")
  abline(0,1, col="red", lwd=2)
}

#  Pearson 
cor_results2 <- lapply(connect_cols[-1], function(x) {
  
  x_vec <- as.numeric(unlist(all_connectance_df[[x]]))
  y_vec <- as.numeric(all_connectance_df$Connectance_orig)
  
  valid_idx <- is.finite(x_vec) & is.finite(y_vec)
  
  x_clean <- x_vec[valid_idx]
  y_clean <- y_vec[valid_idx]
  
  if(length(x_clean) < 3 || length(unique(x_clean)) < 2){
    return(c(r = NA, p_value = NA))
  }
  
  test <- cor.test(x_clean, y_clean, method = "pearson")
  
  c(
    r = unname(test$estimate),
    p_value = test$p.value,
    n = length(x_clean)
  )
})


# dataframe
cor_results_df2 <- as.data.frame(do.call(rbind, cor_results2)) %>%
  mutate(
    Comparison = connect_cols[-1],
    Metric = "Connectance",
    r = as.numeric(r),
    p_value = as.numeric(p_value)
  )

rownames(cor_results_df2) <- connect_cols[-1]

print("Pearson correlation with full network:")
print(cor_results_df2)


library(tibble)
# merge
cor_results_all <- bind_rows(
  mutate(cor_results_df, Metric = "Nestedness"),
  mutate(cor_results_df2, Metric = "Connectance")
)

# 转 data.frame 并把行名变成 Subsampling 列
cor_results_all <- rownames_to_column(as.data.frame(cor_results_all), var = "Subsampling")
rownames(cor_results_all) <- NULL

write.csv(cor_results_all, "./result_260526/network_correlation_results.csv", row.names = FALSE)





##-----------  Visualization -------------

library(ggplot2)
library(ggpubr)

# Convert to long format for top N plants
long_connect_df1 <- all_connectance_df %>%
  pivot_longer(
    cols = c(Connectance_10, Connectance_5, Connectance_3),
    names_to = "Top_Plant_Level",
    values_to = "Connectance_Value"
  ) %>%
  filter(!is.na(Connectance_Value), Connectance_Value != 0) %>%
  mutate(
    Top_Plant_Level = factor(
      Top_Plant_Level,
      levels = c("Connectance_10", "Connectance_5", "Connectance_3"),
      labels = c("Top 10", "Top 5", "Top 3")
    )
  )

# Define custom colors for plots
my_colors <- c(
  "Top 10" = "#66C2A5",
  "Top 5"  = "#4E79A7",
  "Top 3"  = "#F28E2B",
  "flw 5"  = "#4E79A7",
  "flw 3"  = "#F28E2B"
)

reg_significance <- long_connect_df1 %>%
  group_by(Top_Plant_Level) %>%
  summarise(
    lm_res = list(lm(Connectance_Value ~ Connectance_orig, data = cur_data_all())),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = summary(lm_res)$coefficients[2, 4],
    signif_line = ifelse(p_value < 0.05, "solid", "dashed")
  ) %>%
  dplyr::select(Top_Plant_Level, signif_line)

set.seed(2025)
p3 <- ggplot(long_connect_df1, aes(y = Connectance_Value, x = Connectance_orig, color = Top_Plant_Level)) +
  geom_point(shape = 1, size = 2, position = position_jitter(width = 0.005, height = 0.005)) +
  # 根据显著性调整线型
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    aes(linetype = Top_Plant_Level),
    show.legend = TRUE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_linetype_manual(
    values = setNames(reg_significance$signif_line, reg_significance$Top_Plant_Level)
  ) +
  stat_cor(aes(color = Top_Plant_Level), method = "pearson", 
           size = 4,
           label.x.npc = "left",
           label.y.npc = c(0.98, 0.96,0.94)) +
  scale_color_manual(values = my_colors) +
  labs(
    y = "Connectance in Subnetwork of Top-N Plants",
    x = "Connectance of the Full Network",
    color = "Top N Level",
    linetype = "Significance"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )

print(p3)

# Convert to long format for flower abundance top N plants
long_connect_df2 <- all_connectance_df %>%
  pivot_longer(
    cols = c(Connectance_op2_top5, Connectance_op2_top3),
    names_to = "Flw_Plant_Level",
    values_to = "Connectance_Value"
  ) %>%
  filter(!is.na(Connectance_Value), Connectance_Value != 0) %>%
  mutate(
    Flw_Plant_Level = factor(
      Flw_Plant_Level,
      levels = c("Connectance_op2_top5", "Connectance_op2_top3"),
      labels = c("flw 5", "flw 3")
    )
  )
head(all_connectance_df)
reg_significance_flw <- long_connect_df2 %>%
  group_by(Flw_Plant_Level) %>%
  summarise(
    lm_res = list(lm(Connectance_Value ~ Connectance_orig, data = cur_data_all())),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = summary(lm_res)$coefficients[2, 4],
    signif_line = ifelse(p_value < 0.05, "solid", "dashed")
  ) %>%
  dplyr::select(Flw_Plant_Level, signif_line)

set.seed(2025)
p4 <- ggplot(long_connect_df2, aes( x = Connectance_orig,y = Connectance_Value, color = Flw_Plant_Level)) +
  geom_point(shape = 1, size = 2, position = position_jitter(width = 0.005, height = 0.005)) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    aes(linetype = Flw_Plant_Level),
    show.legend = TRUE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_linetype_manual(
    values = setNames(reg_significance_flw$signif_line, reg_significance_flw$Flw_Plant_Level)
  ) +
  stat_cor(
    aes(color = Flw_Plant_Level),
    method = "pearson",
    size = 4,
    label.x.npc = "left",
    label.y.npc = c(0.98, 0.96)) +
  scale_color_manual(values = my_colors) +
  labs(
    y = "Connectance in Subnetwork of Flw-N Plants",
    x = "Connectance of the Full Network",
    color = "Top N Level",
    linetype = "Significance"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )

print(p4) ###img size 800*500


##===============================

#---------3.3 Network metrics for random selection--------
#note: random process will take a few minutes

#================================
library(dplyr)
library(tidyr)
library(pbapply)
library(ggplot2)


set.seed(2025)

plant_pool <- data_count_scaled %>%
  distinct(Study_Network_id, Plant_species, Flower_count_scaled, .keep_all = TRUE)

plant_pool_split <- split(plant_pool, plant_pool$Study_Network_id)

network_list <- unique(plant_pool$Study_Network_id)

# sampling function
get_random_interactions <- function(n_sp){
  
  sampled_plants <- do.call(rbind, lapply(plant_pool_split, function(df){
    
    n_select <- min(n_sp, nrow(df))
    df[sample(nrow(df), n_select), ]
    
  }))
  
  data_interact %>%
    inner_join(
      sampled_plants,
      by = c(
        "Study_Network_id",
        "Plant_original_name" = "Plant_species"
      )
    )
}


calc_nodf <- function(df){
  
  dt <- as.data.table(df)
  
  mat_df <- dcast(
    dt,
    Plant_accepted_name ~ Pollinator_accepted_name,
    value.var = "Interaction_addup",
    fill = 0
  )
  
  mat <- as.matrix(mat_df[, -1, with = FALSE])
  mat[is.na(mat)] <- 0
  
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NA_real_)
  
  tryCatch(
    bipartite::nested(mat, method = "NODF"),
    error = function(e) NA_real_
  )
}

calc_connectance <- function(df){
  
  dt <- as.data.table(df)
  
  mat_df <- dcast(
    dt,
    Plant_accepted_name ~ Pollinator_accepted_name,
    value.var = "Interaction_addup",
    fill = 0
  )
  
  mat <- as.matrix(mat_df[, -1, with = FALSE])
  mat[is.na(mat)] <- 0
  
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NA_real_)
  
  sum(mat > 0) / (nrow(mat) * ncol(mat))
}

library(pbapply)

run_random_metrics <- function(n_sp, n_iter = 1000){
  
  pboptions(type = "timer")
  
  res <- pblapply(1:n_iter, function(i){
    
    random_df <- get_random_interactions(n_sp)
    
    random_df %>%
      group_by(Study_Network_id) %>%
      summarise(
        NODF = calc_nodf(cur_data()),
        Connectance = calc_connectance(cur_data()),
        .groups = "drop"
      )
    
  })
  
  bind_rows(res)
}

# run random 
# 
# random_10_all <- run_random_metrics(10, 1000)
# random_5_all  <- run_random_metrics(5, 1000)
# random_3_all  <- run_random_metrics(3, 1000)
# 
# random_10 <- random_10_all %>%
#   group_by(Study_Network_id) %>%
#   summarise(
#     Random_NODF = mean(NODF, na.rm = TRUE),
#     Random_Connectance = mean(Connectance, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(Method = "random 10")
# 
# random_5 <- random_5_all %>%
#   group_by(Study_Network_id) %>%
#   summarise(
#     Random_NODF = mean(NODF, na.rm = TRUE),
#     Random_Connectance = mean(Connectance, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(Method = "random 5")
# 
# random_3 <- random_3_all %>%
#   group_by(Study_Network_id) %>%
#   summarise(
#     Random_NODF = mean(NODF, na.rm = TRUE),
#     Random_Connectance = mean(Connectance, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(Method = "random 3")
# random_all <- bind_rows(random_10, random_5, random_3)
# 
# saveRDS(random_all,"random_all_metrics.rds")

library("dplyr")
random_all <- readRDS("random_all_metrics.rds") %>%
  dplyr::mutate(Method = trimws(Method)) %>%
  dplyr::mutate(Method = dplyr::case_when(
    Method %in% c("random10", "Random10", "random 10") ~ "Top 10",
    Method %in% c("random5", "Random5", "random 5")   ~ "Top 5",
    Method %in% c("random3", "Random3", "random 3")   ~ "Top 3",
    TRUE ~ Method
  ))

random_all$Method <- factor(random_all$Method,
                            levels = c("Top 10", "Top 5", "Top 3"))


random_nodf_df <- random_all %>%
  dplyr::select(Study_Network_id, Method, Random_NODF)

random_conn_df <- random_all %>%
  dplyr::select(Study_Network_id, Method, Random_Connectance)


nodf_random_df <- all_nestedness_df %>%
  dplyr::select(Study_Network_id, Nestedness_orig) %>%
  left_join(
    pivot_wider(
      random_nodf_df,
      names_from = Method,
      values_from = Random_NODF,
      values_fn = mean
    ),
    by = "Study_Network_id"
  ) %>%
  dplyr::rename(
    `Top 10` = `Top 10`,
    `Top 5`  = `Top 5`,
    `Top 3`  = `Top 3`
  )

conn_random_df <- all_connectance_df %>%
  dplyr::select(Study_Network_id, Connectance_orig) %>%
  left_join(
    pivot_wider(
      random_conn_df,
      names_from = Method,
      values_from = Random_Connectance,
      values_fn = mean   
    ),
    by = "Study_Network_id"
  )%>%
  dplyr::rename(
    `Top 10` = `Top 10`,
    `Top 5`  = `Top 5`,
    `Top 3`  = `Top 3`
  )

# -----------------------------
#long format
# -----------------------------
random_long_nodf <- nodf_random_df %>%
  pivot_longer(
    cols = c(`Top 10`, `Top 5`, `Top 3`),
    names_to = "Method",
    values_to = "Value"
  ) %>%
  mutate(Method = factor(Method,
                         levels = c("Top 10","Top 5","Top 3")))

random_long_conn <- conn_random_df %>%
  pivot_longer(
    cols = c(`Top 10`, `Top 5`, `Top 3`),
    names_to = "Method",
    values_to = "Value"
  ) %>%
  mutate(Method = factor(Method,
                         levels = c("Top 10","Top 5","Top 3")))

# -----------------------------
# NODF plot
# -----------------------------
p_random_nodf <- ggplot(random_long_nodf,
                        aes(x = Nestedness_orig,
                            y = Value,
                            color = Method)) +
  
  geom_point(shape = 1, size = 2,
             position = position_jitter(width = 0.02, height = 0.2)) +
  
  geom_smooth(method = "lm", se = FALSE,
              aes(linetype = Method),
              linewidth = 0.8) +
  
  stat_cor(aes(color = Method),
           method = "pearson",
           size = 4,
           label.x.npc = "left") +
  scale_y_continuous(limits = c(0, 60)) +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = c(
    "Top 10" = "solid",
    "Top 5"  = "solid",
    "Top 3"  = "solid"
  )) +
  
  labs(
    x = "NODF of the Full Network",
    y = "NODFin Subnetworks of Random-N Plants",
    color = "Method",
    linetype = "Method"
  )+
  theme_classic(base_size=12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )


# -----------------------------
# Connectance plot
# -----------------------------
p_random_conn <- ggplot(random_long_conn,
                        aes(x = Connectance_orig,
                            y = Value,
                            color = Method)) +
  
  geom_point(shape = 1, size = 2,
             position = position_jitter(width = 0.005, height = 0.005)) +
  
  geom_smooth(method = "lm", se = FALSE,
              aes(linetype = Method),
              linewidth = 0.8) +
  
  stat_cor(aes(color = Method),
           method = "pearson",
           size = 4,
           label.x.npc = "left") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = c(
    "Top 10" = "solid",
    "Top 5"  = "solid",
    "Top 3"  = "solid"
  )) +
  
  labs(
    x = "Connectance of the Full Network",
    y = "Connectance in Subnetworks of Random-N Plants",
    color = "Method",
    linetype = "Method"
  )+
  theme_classic(base_size=12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )


p_random_nodf
p_random_conn

##################
#----------plot all---------
#############
library(cowplot)

p1.1 <- p1 + theme(legend.position = "none")
p2.1 <- p2 + theme(legend.position = "none" 
                   #,axis.title.y = element_blank()
)
p3.1 <- p3 + theme(legend.position = "none")
p4.1 <- p4 + theme(legend.position = "none")

p_random_nodf.1<-p_random_nodf+ theme(legend.position = "none")
p_random_conn.1<-p_random_conn + theme(legend.position = "none")


top_row <- plot_grid(
  p1.1, p2.1, p_random_nodf.1,
  ncol = 3,
  align = "hv",
  axis = "tblr"
)

bottom_row <- plot_grid(
  p3.1, p4.1, p_random_conn.1,
  ncol = 3,
  align = "hv",
  axis = "tblr"
)
legend_plot <- p_random_conn +
  theme(
    legend.position = "bottom"
  ) +
  guides(
    color = guide_legend(nrow = 1),
    linetype = guide_legend(nrow = 1)
  )

legend <- cowplot::get_legend(legend_plot)

combined_plot <- plot_grid(
  top_row,
  bottom_row,
  legend,
  ncol = 1,
  rel_heights = c(1, 1, 0.12)
)

combined_plot#size:980*870
ggsave("./result_260526/network.matrix.png", combined_plot, width = 15, 
       height = 9.5, units = "in", dpi = 300)


########################################################################################

saveRDS(result10,"data/processd/result10.rds")  
saveRDS(data_merge,"data/processd/data_merge.rds") 
##########################################################################


#Shannon test
# #plot
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# 
# result_all<-read.csv("result_all_published.csv",header=TRUE)
# plotshannon<-result_all[,c("Study_Network_id","shannon_index_abun5","shannon_index_abun3","shannon_index_flw5","shannon_index_flw3")]
# 
# plotshannon<-result_all[,c("Study_Network_id","shannon_index_abun5","shannon_index_abun3","shannon_index_flw5","shannon_index_flw3")]%>% 
#   filter(Study_Network_id %in% Op1better$Study_Network_id)
# 
# 
# shannon_long <- pivot_longer(plotshannon, -Study_Network_id, names_to = "Option", values_to = "Shannon_index")
# shannon_long_clean <- shannon_long[complete.cases(shannon_long), ]
# print(shannon_long_clean)
# 
# data_summary <- shannon_long_clean %>%
#   group_by(Option) %>%
#   summarise(
#     Mean_Shannon_index = mean(Shannon_index, na.rm = TRUE),
#     SD_Shannon_index = sd(Shannon_index, na.rm = TRUE)
#   )
# 
# print(data_summary)
# 
# 
# custom_order <- c("shannon_index_abun5","shannon_index_flw5","shannon_index_abun3","shannon_index_flw3")
# 
# shannon_long_clean$Option <- factor(shannon_long_clean$Option, levels = custom_order)
# 
# non_missing_count <- shannon_long_clean %>%
#   group_by(Option) %>%
#   summarise(
#     n_non_missing = sum(!is.na(Shannon_index))
#   )
# 
# data_summary <- left_join(data_summary, non_missing_count, by = "Option") %>%
#   mutate("Focus Species" = case_when(
#     Option %in% c("shannon_index_abun5", "shannon_index_abun3") ~ "Flower abundance",
#     Option %in% c("shannon_index_flw5", "shannon_index_flw3") ~ "Flower shapes"
#   )) %>%
#   mutate("Plant Number" = case_when(
#     Option %in% c("shannon_index_abun5", "shannon_index_flw5") ~ "Top5",
#     TRUE ~ "Top3"
#   ))
# 
# data_summary$`Plant Number` <- factor(data_summary$`Plant Number`, levels = c( "Top5", "Top3"))
# 
# 
# option_info <- data_summary %>%
#   dplyr::select(Option, `Focus Species`, `Plant Number`)
# 
# shannon_long_clean <- shannon_long_clean %>%
#   left_join(option_info, by = "Option")
# shannon_long_clean$`Plant Number` <- factor(shannon_long_clean$`Plant Number`, levels = c( "Top5", "Top3"))
# 
# group_colors <- c(
#                   "shannon_index_abun5"  = "#1b9e77",  
#                   "shannon_index_abun3"  = "#1b9e77",  
#                   
#                   "shannon_index_flw5"  = "#9E9AC8",
#                   "shannon_index_flw3"  = "#9E9AC8"
# )
# 
# 
# focus_species_colors <- c(
#   "Flower abundance" = "#1B9E77",
#   "Flower shapes" = "#9E9AC8"
# )
# 
# # 设置 factor 顺序
# data_summary$`Focus Species` <- factor(data_summary$`Focus Species`,
#                                        levels = c("Flower abundance", "Flower shapes"))
# 
# shannon_long_clean$`Focus Species` <- factor(shannon_long_clean$`Focus Species`,
#                                           levels = c("Flower abundance", "Flower shapes"))
# 
# set.seed(2025)
# 
# plot <- ggplot(data_summary, aes(x = Option, y = Mean_Shannon_index, color = `Focus Species`)) +
#   geom_jitter(data = shannon_long_clean,
#               aes(x = Option, y = Shannon_index, color = `Focus Species`),
#               width = 0.2, size = 1, alpha = 0.5) +
#   geom_crossbar(aes(ymin = Mean_Shannon_index - SD_Shannon_index, ymax = SD_Shannon_index + Mean_Shannon_index),
#                 position = position_dodge(width = 0.3), width = 0.5, color = "grey40") +
#   geom_text(aes(label = paste("n =", n_non_missing)),
#             position = position_dodge(width = 0.5), vjust = -1, size = 4, color = "black") +
#   geom_text(aes(label = paste(round(Mean_Shannon_index, 2))),
#             position = position_dodge(width = 0.5), vjust = 1.5, size = 4, color = "black") +
#   facet_wrap(~ `Plant Number`, nrow = 1, scales = "free_x") +
#   scale_color_manual(values = focus_species_colors) +
#   labs(x = "Subsampling Option", y = "Shannon index Captured", color = "Plants selected based on") +
#   theme_minimal(base_size = 12) +
#   theme(
#     panel.background = element_blank(),    
#     panel.grid = element_blank(),           
#     panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),  
#     plot.margin = margin(10, 10, 10, 10),   
#     strip.background = element_rect(fill = "white", color = NA), 
#     strip.text = element_text(size = 12, face = "bold"),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     legend.position = "right",
#     legend.title = element_text(size = 11),
#     legend.text = element_text(size = 10)
#   )+ 
#   guides(color = guide_legend(override.aes = list(size = 5)))
# 
# 
# plot
# 
# ggsave("result251105_published/shannon.png", plot, width = 7.5, height = 4.5, units = "in", dpi = 300)
# 
# #############
# op2_data <- data_merge %>% filter(Study_Network_id %in% Op2better$Study_Network_id) %>%
#   left_join(traits %>% dplyr::select(Plant_accepted_name, flw_shape_revised), by = "Plant_accepted_name") %>%
#   filter(!is.na(flw_shape_revised), !is.na(Pollinator_accepted_name))
# 
# op1_data <- data_merge %>% filter(Study_Network_id %in% Op1better$Study_Network_id) %>%
#   left_join(traits %>% dplyr::select(Plant_accepted_name, flw_shape_revised), by = "Plant_accepted_name") %>%
#   filter(!is.na(flw_shape_revised), !is.na(Pollinator_accepted_name))
# 
# 
# ############ 哪个研究的site有超过50种同时开花？
# metadata<-readRDS("Interaction_data_published.rds")#The EuPPollNet interaction data
# meta_count<-readRDS("Flower_counts_published.rds")#The EuPPollNet flower data
# colnames(meta_count)
# a<-meta_count%>%filter(Study_id=="20_Hoiss")
# 
# over50 <- meta_count %>%
#   ungroup() %>%
#   filter(!is.na(Plant_species)) %>%
#   filter(Flower_count > 0) %>%
#   group_by(Study_id, Network_id, Day, Month, Year, Site_id) %>%
#   summarise(n_species = n_distinct(Plant_species), .groups = "drop") %>%
#   filter(n_species >= 50)


