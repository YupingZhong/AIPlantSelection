###############
#sampling completeness

#############################
# -----------------------------
# 完整版：Pollinator Sampling Completeness
# -----------------------------
library(dplyr)
library(iNEXT)
library(ggplot2)
data_merge <- readRDS("data/processed/data_merge.rds")
# -----------------------------
# 1. 准备每个网络的 pollinator abundance
# -----------------------------
pollinator_abundance_list <- data_merge %>%
  filter(!is.na(Pollinator_accepted_name)) %>%
  group_by(Study_Network_id, Pollinator_accepted_name) %>%
  summarise(abundance = sum(Interaction_addup), .groups = "drop") %>%
  group_by(Study_Network_id) %>%
  summarise(abundance_vec = list(abundance), .groups = "drop")

# -----------------------------
# 2. 计算每个网络的 sample coverage
# -----------------------------
library(dplyr)
library(iNEXT)

coverage_results <- lapply(1:nrow(pollinator_abundance_list), function(i){
  
  net_id <- pollinator_abundance_list$Study_Network_id[i]
  abund  <- pollinator_abundance_list$abundance_vec[[i]]
  
  # 防御式检查
  if(length(abund) < 2 || sum(abund) < 2){
    return(data.frame(
      Study_Network_id = net_id,
      coverage_est = NA
    ))
  }
  
  info <- DataInfo(abund, datatype = "abundance")
  
  data.frame(
    Study_Network_id = net_id,
    coverage_est = info$SC
  )
})

coverage_results <- do.call(rbind, coverage_results)

summary(coverage_results$coverage_est)

# -----------------------------
# 3. 保存 CSV
# -----------------------------
write.csv(coverage_results, "pollinator_sampling_completeness.csv", row.names = FALSE)

#################################################################

sampling_c <- ggplot(coverage_results, aes(x = coverage_est)) +
  geom_histogram(
    bins = 20,
    fill = "#A1D8E8",
    colour = "white",
    linewidth = 0.3
  ) +
  labs(
    x = "Pollinator sampling completeness",
    y = "Number of networks"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(2, "pt"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )
  
sampling_c
ggsave("/Chap1_TargetPlant_to_monitor/result_260723/sampling_complete.png", sampling_c, width = 5, height = 3, units = "in", dpi = 300)
