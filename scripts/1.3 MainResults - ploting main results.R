
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

