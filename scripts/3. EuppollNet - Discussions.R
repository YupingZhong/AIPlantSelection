
#-------------Disscussions--------------------------------------------------
library("dplyr")
getwd()
setwd("E:/Chap1_TargetPlant_to_monitor")
####################################################
# data prepare
data_count_scaled<-readRDS("data/raw/data_count_scaled_published.rds")#全部互作数据
data_interact<-readRDS("data/raw/data_interact_published.rds")
data_merge<- readRDS("data/processed/data_merge")

unique(data_interact$Study_Network_id)
unique(data_count_scaled_species$Study_Network_id)
unique(data_merge$Study_Network_id)


###===============================================================================

# ------------1. Pollinator Species constantly Missed by Abundant 10 subsampling strategy

###=====================================================================================
library(stringr)
library(tidyr)
library(dplyr)
result10<-readRDS("/data/processed/selected_plant_result_abun10.rds")
####################### order level
order_pollinator_10 <- result10 %>%
  ungroup() %>%
  filter(!is.na(Pollinator_accepted_name)) %>%
  distinct(Study_Network_id, Pollinator_order, Pollinator_family) %>%  # 去重
  group_by(Pollinator_order, Pollinator_family) %>%
  summarise(Captured = n(), .groups = "drop")

#  Total networks where Apidae were observed
order_pollinator_total <- data_interact %>%
  ungroup() %>%
  filter(!is.na(Pollinator_accepted_name)) %>%
  distinct(Study_Network_id, Pollinator_order, Pollinator_family) %>%  # 去重
  group_by(Pollinator_order, Pollinator_family) %>%
  summarise(Networks = n(), .groups = "drop")

# Join captured data with total and calculate missed counts and missed rate (%)
order_Pollinator_missed <- order_pollinator_total %>%
  left_join(order_pollinator_10, by = c("Pollinator_order", "Pollinator_family")) %>%
  replace_na(list(Captured = 0)) %>%
  mutate(
    Missed = Networks - Captured,
    Missed_rate = (Missed / Networks) * 100   # Corrected formula: missed / total, not captured / total
  ) %>%
  filter(!Missed_rate == 0,!is.na(Pollinator_family)) #只看

order_Pollinator_missed %>%
  filter(Networks < Captured)

# Export results for record
write.csv(order_Pollinator_missed, "/Chap1_TargetPlant_to_monitor/result_260526/order_Pollinator_missed.csv")

###################### plot
library(ggplot2)
library(dplyr)
library(tidyr)

# ⭐ 只保留四个主要目
target_orders <- c("Diptera", "Lepidoptera", "Hymenoptera", "Coleoptera")

order_plot_df <- order_Pollinator_missed %>%
  filter(Pollinator_order %in% target_orders) %>%   # ⭐ 筛选
  mutate(Captured = Networks - Missed) %>%
  dplyr::select(Pollinator_order, Pollinator_family, Captured, Missed) %>%
  pivot_longer(cols = c(Captured, Missed),
               names_to = "Status",
               values_to = "Count") %>%
  mutate(Status = factor(Status, levels = c("Missed", "Captured")))

order_plot_df <- order_plot_df %>%
  mutate(Pollinator_order = factor(Pollinator_order,
                                   levels = c("Diptera", "Lepidoptera",
                                              "Hymenoptera", "Coleoptera")))

# 排序
family_order <- order_Pollinator_missed %>%
  filter(Pollinator_order %in% target_orders) %>%   # ⭐ 同步筛选
  arrange(desc(Networks)) %>%
  pull(Pollinator_family)

order_plot_df <- order_plot_df %>%
  mutate(Pollinator_family = factor(Pollinator_family, levels = family_order))

# ⭐ 准备miss rate标签数据
label_df <- order_Pollinator_missed %>%
  filter(Pollinator_order %in% target_orders) %>%
  mutate(
    Miss_rate_pct = round(Missed / Networks * 100, 1),
    Pollinator_family = factor(Pollinator_family, levels = family_order)
  ) %>%
  dplyr::select(Pollinator_order, Pollinator_family, Networks, Miss_rate_pct)


# ⭐ 作图
p <- ggplot(order_plot_df, 
            aes(x = Pollinator_family, y = Count, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = c(
      Captured = "#1B9E77",
      Missed   = "#D95F02"
    ),
    labels = c(
      Captured = "Captured networks",
      Missed   = "Missed networks"
    )
  ) +
  facet_wrap(~ factor(Pollinator_order,
             levels = c("Diptera", "Lepidoptera", 
                        "Hymenoptera", "Coleoptera")),
             nrow   = 2,
             scales = "free_x") +
  labs(
    x    = NULL,
    y    = "Number of networks",
    fill = NULL
  ) +
  # geom_text(data = label_df,
  #           aes(x = Pollinator_family, 
  #               y = Networks + 5,  # 在柱子顶部上方
  #               label = paste0(Miss_rate_pct, "%"),
  #               fill = NULL),
  #           vjust = 0,
  #           size = 2.5,
  #           fontface = "bold",
  #           color = "#D95F02") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

p

ggsave("./result_260526/order_Pollinator_missed_bar.png",
       p, width = 9, height = 8,dpi = 300)


######################## for Apidae only
# Summarize captured counts for Apidae pollinators in top 10 selected networks
#  Captured counts in top 10 networks
length(unique(data_interact$Study_Network_id))
length(unique(result10$Study_Network_id))

pollinator_10 <- result10 %>%
  ungroup() %>%
  filter(!is.na(Pollinator_accepted_name),
         Pollinator_family == "Apidae") %>%
  distinct(Study_Network_id, Pollinator_genus, Pollinator_accepted_name) %>%  # 去重
  group_by(Pollinator_genus, Pollinator_accepted_name) %>%
  summarise(Captured = n(), .groups = "drop")

#  Total networks where Apidae were observed
pollinator_total <- data_interact %>%
  ungroup() %>%
  filter(!is.na(Pollinator_accepted_name),
         Pollinator_family == "Apidae") %>%
  distinct(Study_Network_id, Pollinator_genus, Pollinator_accepted_name) %>%  # 去重
  group_by(Pollinator_genus, Pollinator_accepted_name) %>%
  summarise(Networks = n(), .groups = "drop")

# Join captured data with total and calculate missed counts and missed rate (%)
Pollinator_missed <- pollinator_total %>%
  left_join(pollinator_10, by = c("Pollinator_accepted_name", "Pollinator_genus")) %>%
  replace_na(list(Captured = 0)) %>%
  mutate(
    Missed = Networks - Captured,
    Missed_rate = (Missed / Networks) * 100   # Corrected formula: missed / total, not captured / total
  ) %>%
  # Exclude species only identified to genus level (single word names)
  filter(str_count(Pollinator_accepted_name, "\\S+") != 1)%>%
  filter(!Missed_rate == 0) #只看

Pollinator_missed %>%
  filter(Networks < Captured)

# Export results for record
write.csv(Pollinator_missed, "/Chap1_TargetPlant_to_monitor/result_260526/Pollinator_missed.csv")

##output
library(flextable)
library(officer)
library(dplyr)

Pollinator_missed_table <- Pollinator_missed %>%
  filter(Networks > 5 )%>% #总共出现5次及以下的不纳入统计
  mutate(
    `Missed rate numeric` = Missed_rate
  ) %>%
  arrange(Pollinator_genus, desc(`Missed rate numeric`)) %>%
  group_by(Pollinator_genus) %>%   # ✅ 新增：用于 genus 分组显示
  mutate(
    `NO.` = row_number(),
    `Pollinator genus display` = ifelse(row_number() == 1, Pollinator_genus, "")
  ) %>%
  ungroup() %>%
  mutate(
    `Missed rate (%)` = sprintf("%.2f", Missed_rate)
  ) %>%
  rename(
    `Pollinator genus` = `Pollinator genus display`,   # ✅ 用 display 替换
    `Pollinator species` = Pollinator_accepted_name,
    `Networks total` = Networks,
    `Networks missed` = Missed
  ) %>%
  dplyr::select(`NO.`, `Pollinator genus`, `Pollinator species`,
         `Networks total`, `Networks missed`, `Missed rate (%)`)


# 创建 flextable（全部加 flextable:: 防冲突）
ft <- flextable::flextable(Pollinator_missed_table) %>%
  flextable::theme_booktabs() %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::font(fontname = "Arial", part = "all") %>%   # ✅ 关键修复点
  flextable::align(align = "center", part = "all") %>%
  flextable::set_header_labels(
    `NO.` = "NO.",
    `Pollinator genus` = "Pollinator\ngenus",
    `Pollinator species` = "Pollinator\nspecies",
    `Networks total` = "Networks\ntotal",   # ✅ 修正列名
    `Networks missed` = "Networks\nmissed",
    `Missed rate (%)` = "Missed rate\n(%)"
  )

#  species 斜体（必须单独写，不能放在 pipe 里）
ft <- flextable::compose(
  ft,
  j = "Pollinator species",
  value = flextable::as_paragraph(
    flextable::as_i(.data$`Pollinator species`)   # ✅ 正确写法
  )
)


ft <- flextable::autofit(ft)


doc <- read_docx()
doc <- body_add_par(doc, "Table S4. Missed Apidae species across networks.", style = "heading 2")
doc <- body_add_flextable(doc, value = ft)


print(doc, target = "./result_260526/Apidae_Network_Counts.docx")


# #####################
# ## Visualization setup
# library(rotl)
# library(ape)
# library(stringr)
# library(ggtree)
# library(ggtreeExtra)
# library(ggplot2)
# library(viridis)
# library(dplyr)
# # If planning to add multiple fill scales, load ggnewscale (not added here yet)
# # library(ggnewscale)
# 
# # 1. Clean and standardize species names to "Genus species" format with capitalization
# clean_species_name <- function(name) {
#   name <- str_squish(name)
#   words <- str_split(name, " ", simplify = TRUE)
#   if (ncol(words) < 2) return(NA)
#   paste(str_to_title(words[1]), tolower(words[2]))
# }
# 
# Pollinator_missed <- Pollinator_missed %>%
#   mutate(Clean_name = sapply(Pollinator_accepted_name, clean_species_name)) %>%
#   filter(!is.na(Clean_name))
# 
# # 2. Match species names to Open Tree Taxonomy (OTL) and filter valid matches
# matched_names_clean <- tnrs_match_names(unique(Pollinator_missed$Clean_name))
# 
# valid_matches <- matched_names_clean %>%
#   filter(!flags %in% c("unplaced", "sibling_higher", "incertae_sedis")) %>%
#   filter(!is.na(ott_id))
# 
# ott_ids_valid <- valid_matches$ott_id
# 
# # 3. Get phylogenetic subtree for matched species
# tree <- tol_induced_subtree(ott_ids = ott_ids_valid)
# 
# # 4. Join tree data and pollinator data, prepare labels for plotting
# Pollinator_captured_clean <- Pollinator_missed %>%
#   inner_join(valid_matches %>% select(unique_name, ott_id),
#              by = c("Clean_name" = "unique_name")) %>%
#   mutate(label = gsub(" ", "_", Pollinator_accepted_name))  # Replace spaces with underscores for tip labels
# 
# # 5. Clean tree tip labels by removing OTT IDs to match data labels
# tree$tip.label <- str_replace(tree$tip.label, "_ott[0-9]+$", "")
# 
# # 6. Keep only species present in both tree and data
# common_labels <- intersect(tree$tip.label, Pollinator_captured_clean$label)
# 
# tree_filtered <- drop.tip(tree, setdiff(tree$tip.label, common_labels))
# Pollinator_captured_filtered <- Pollinator_captured_clean %>% filter(label %in% common_labels)
# 
# # 7. Plot circular phylogenetic tree with species labels
# p <- ggtree(tree_filtered, layout = "circular") +
#   geom_tiplab2(aes(label = label), size = 3.2, offset = 0.1, align = TRUE)
# 
# # Add outer ring for Missed rate (%)
# p2 <- p + 
#   geom_fruit(
#     data = Pollinator_captured_filtered,
#     geom = geom_col,
#     mapping = aes(x = Missed_rate, y = label, fill = Missed_rate),
#     orientation = "y",
#     offset = 1.0,   # Adjust offset to avoid overlap with tip labels
#     width = 0.6
#   ) +
#   scale_fill_gradient(
#     name = "Missed rate (%)",
#     low = "#9E9AC8",  # or chose green color "#1B9E77"
#     high = "darkred"
#   )
# 
# # Print the plot with Missed rate visualization
# print(p2) # 1100*1000
# 
# ggsave("/Chap1_TargetPlant_to_monitor/result_260125_with_cap/Bee_missed.png",
#        p2, width = 11, height = 10, units = "in", dpi = 300)
# 
# ##################################################
# library(ggplot2)
# library(dplyr)
# library(forcats)
# 
# # 按 Missed_rate 排序
# Pollinator_missed <- Pollinator_missed %>%
#   mutate(Pollinator_species = Pollinator_accepted_name) %>%
#   arrange(Missed_rate) %>%
#   mutate(Pollinator_species = fct_reorder(Pollinator_species, Missed_rate))
# 
# # 绘图
# p_missed <- ggplot(Pollinator_missed, aes(x = Pollinator_species, y = Missed_rate, fill = Missed_rate)) +
#   geom_col() +
#   coord_flip() +
#   scale_y_continuous(labels = scales::percent_format(scale = 1)) +
#   scale_fill_gradient(low = "#9E9AC8", high = "darkred") +  # 渐变色
#   labs(
#     x = "Pollinator species of Apidae",
#     y = "Missed rate (%)",
#     fill = "Missed rate (%)",
#     title = "Missed rate of Apidae pollinators in top 10 plants"
#   ) +
#   theme_classic() +
#   theme(
#     axis.text.y = element_text(size = 8),
#     axis.text.x = element_text(size = 10),
#     legend.position = "right"
#   )
# 
# print(p_missed)
# 
# # 保存图片
# ggsave("Apidae_missed_rate_gradient.png", p_missed, width = 6, height = 8, dpi = 300)
# 
# 
# # 保存图片
# ggsave("Apidae_missed_rate.png", p_missed, width = 6, height = 8, dpi = 300)
# 


###=================================================================================

###-------------2. How many pollinator to identify to capture most of the interaction in Europe
##(Implication for classifier training)

#======================================================================================

library(ggplot2)
library(dplyr)
library(purrr)
#read data
metadata<-readRDS("Interaction_data_published.rds")%>%
  mutate(Study_Network_id = paste(Study_id, Network_id, sep = "_"))
meta_count<-readRDS("Flower_counts_published.rds")%>%
  mutate(Study_Network_id = paste(Study_id, Network_id, sep = "_"))
#首先删除那些pollinator只鉴定到属和大类的记录
Poll_sp<-unique(metadata$Pollinator_accepted_name)#2668 POLLINATORS 
Plant_sp<-unique(metadata$Plant_accepted_name)#1543 POLLINATORS

#check the crop site
flowercout_check<-meta_count %>%
  group_by(Study_Network_id) %>%
  mutate(sum_abun = sum(Flower_count),percent_abun = Flower_count/sum_abun)%>%
  ungroup()%>%
  filter(percent_abun>0.5 )

flowercoutsite<-flowercout_check%>%distinct(Study_Network_id,Plant_species)
flowercout_check%>%distinct(Study_Network_id)
flowercoutsite_list<-unique(flowercout_check$Plant_species)

#list argriculture plants with over 50% ralative abundance:
#this step is to exclude crop dominant site

crop_species <- c(
  "Medicago sativa",
  "Helianthus annuus",
  "Brassica napus",
  "Vicia faba",
  "Solanum tuberosum",
  "Asparagus officinalis",
  "Brassica oleracea"
)


# 1️⃣ site level filter
site_info <- metadata %>%
  group_by(Study_Network_id) %>%
  summarise(
    site_visit_time = sum(Interaction, na.rm = TRUE),
    plant_sp = n_distinct(Plant_accepted_name)
  ) %>%
  filter(plant_sp < 10)

#  abundance level filter（只保留 crop）
crop_dom <- meta_count %>%
  filter(Plant_species %in% crop_species) %>%
  group_by(Study_Network_id) %>%
  mutate(
    sum_abun = sum(Flower_count, na.rm = TRUE),
    percent_abun = Flower_count / sum_abun
  ) %>%
  filter(percent_abun > 0.5) %>%
  ungroup()

# 找交集 Study_Network_id
crop_dom_site <- intersect(
  site_info$Study_Network_id,
  unique(crop_dom$Study_Network_id)
)
crop_dom_site


dplyr::n_distinct(metadata$Study_Network_id)

##========================================

#how many pollinator to be identified

##=============================================
# Filter pollinator species
Interaction_data_pol_sp <- metadata %>%
  filter(!is.na(Pollinator_accepted_name)) %>%
  filter(Pollinator_rank == "SPECIES") # only species level

# Filter plant species
Interaction_data_plant_sp <- metadata %>%
  filter(!is.na(Plant_accepted_name)) %>%
  filter(Plant_rank == "SPECIES")    

# Unique pollinator and plant species
length(unique(Interaction_data_pol_sp$Pollinator_accepted_name)) 
# Result: 2,223 pollinators

length(unique(Interaction_data_plant_sp$Plant_accepted_name)) 
# Result: 1,411 plants

sum(metadata$Interaction)
sum(Interaction_data_pol_sp$Interaction)

#======================================================================================
# 🔴 方案A：平均相对丰度法（推荐）- 主要图表
#======================================================================================

pollinator_importance_A <- Interaction_data_pol_sp %>%
  filter(!Study_Network_id %in% crop_dom_site) %>%
  group_by(Study_Network_id) %>%
  mutate(network_total_inter = sum(Interaction, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Pollinator_accepted_name, Study_Network_id) %>%
  summarize(
    relative_abundance = sum(Interaction, na.rm = TRUE) / 
      first(network_total_inter),
    .groups = "drop"
  ) %>%
  # 补全所有组合（包括0交互）
  complete(Pollinator_accepted_name, Study_Network_id,
           fill = list(relative_abundance = 0)) %>%
  group_by(Pollinator_accepted_name) %>%
  summarize(
    Mean_relative_abundance = mean(relative_abundance, na.rm = TRUE),
    n_networks = n_distinct(Study_Network_id),
    n_networks_present = sum(relative_abundance > 0),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_relative_abundance)) %>%
  mutate(
    Rank = row_number(),
    Cumulative_mean_abundance = cumsum(Mean_relative_abundance),
    Cumulative_pct = (Cumulative_mean_abundance / sum(Mean_relative_abundance)) * 100
  )

# 找到95%的截点
target_pct_A <- 95
target_row_A <- pollinator_importance_A %>%
  filter(Cumulative_pct >= target_pct_A) %>%
  slice_min(Rank)

target_x_A <- target_row_A$Rank
target_y_A <- target_row_A$Cumulative_pct
target_species_A <- target_row_A$Pollinator_accepted_name

# 📊 绘制方案A的主图
plot_A <- ggplot(pollinator_importance_A, aes(x = Rank, y = Cumulative_pct)) +
  geom_point(size = 2, shape = 1, color = "#FFB300", alpha = 0.8) +
  geom_line(size = 0.8, color = "#FFB300", alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, linetype = "solid",
              fill = "#3498DB", color = "#3498DB", alpha = 0.15, size = 0.6) +
  
  # 参考线
  geom_vline(xintercept = target_x_A, linetype = "dashed", 
             color = "#66C2A5", size = 1, alpha = 0.7) +
  geom_hline(yintercept = target_y_A, linetype = "dashed", 
             color = "#66C2A5", size = 1, alpha = 0.7) +
  
  # 注释
  annotate("text",
           x = target_x_A, y = target_y_A,
           label = paste0(target_pct_A, "% captured\n", 
                          target_x_A, " species needed\n(",
                          sprintf("%.1f%%", (target_x_A/nrow(pollinator_importance_A))*100), 
                          " of total)"),
           color = "#66C2A5", fontface = "bold",
           hjust = -0.1, vjust = 1.5, size = 4.5,
           bbox = list(boxcolour = "white", alpha = 0.9)) +
  
  labs(
    #title = "Scheme A: Mean Relative Abundance Across Networks",
    #subtitle = "Crop-dominated sites excluded; Zero interactions included",
    x = "Species Rank (sorted by mean relative interactions)",
    y = "Cumulative Relative Interactions (%)"
  ) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 10)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
    axis.title = element_text(face = "bold"),
    #panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

print(plot_A)


cat("\n========== 🔴 方案A 结果 ==========\n")
cat(sprintf("总物种数: %d\n", nrow(pollinator_importance_A)))
cat(sprintf("捕获 %d%% 平均相对丰度需要: %d 个物种 (%0.1f%%)\n", 
            target_pct_A, target_x_A, (target_x_A/nrow(pollinator_importance_A))*100))
cat(sprintf("该物种为: %s\n\n", target_species_A))
cat("前15个最重要物种：\n")
print(pollinator_importance_A %>% 
        slice(1:15) %>%
        select(Rank, Pollinator_accepted_name, Mean_relative_abundance, 
               n_networks_present, Cumulative_pct))


#======================================================================================
# 🟢 分目展示：用方案A（推荐）
#======================================================================================

cat("\n\n========== 按目（Order）分组分析 - 方案A ==========\n\n")

pollinator_by_order_A <- Interaction_data_pol_sp %>%
  filter(!Study_Network_id %in% crop_dom_site) %>%
  filter(!is.na(Pollinator_order)) %>%
  group_by(Study_Network_id) %>%
  mutate(network_total = sum(Interaction, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Pollinator_order, Pollinator_accepted_name, Study_Network_id) %>%
  summarize(
    species_inter = sum(Interaction, na.rm = TRUE),
    network_total = first(network_total),
    .groups = "drop"
  ) %>%
  mutate(relative_abundance = species_inter / network_total) %>%
  group_by(Pollinator_order, Pollinator_accepted_name) %>%
  summarize(
    Mean_relative_abundance = mean(relative_abundance, na.rm = TRUE),
    n_networks_present = n_distinct(Study_Network_id),
    .groups = "drop"
  ) %>%
  group_by(Pollinator_order) %>%
  arrange(desc(Mean_relative_abundance), .by_group = TRUE) %>%
  mutate(
    Rank = row_number(),
    Cumulative_pct = (cumsum(Mean_relative_abundance) / sum(Mean_relative_abundance)) * 100
  ) %>%
  ungroup()

# 找每个目的95%截点
target_by_order_A <- pollinator_by_order_A %>%
  group_by(Pollinator_order) %>%
  filter(Cumulative_pct >= 95) %>%
  slice_min(Rank) %>%
  ungroup()

print(target_by_order_A %>%
        dplyr::select(Pollinator_order, Rank, Pollinator_accepted_name, 
               Mean_relative_abundance, Cumulative_pct))

# 选择主要的目（至少10个物种且有95%数据）
orders_to_plot <- target_by_order_A %>%
  inner_join(
    pollinator_by_order_A %>%
      group_by(Pollinator_order) %>%
      summarize(n_sp = n(), .groups = "drop"),
    by = "Pollinator_order"
  ) %>%
  filter(n_sp >= 10) %>%
  pull(Pollinator_order)

# 只保留数据量足够的目（至少10个物种）
# 手动指定目的顺序，替换原来的 orders_to_plot
orders_to_plot <- c("Diptera", "Lepidoptera", "Hymenoptera", "Coleoptera", "Hemiptera")

# 验证这5个目在数据里都存在
setdiff(orders_to_plot, unique(pollinator_by_order_A$Pollinator_order))
# 如果输出为 character(0) 说明全部存在，否则检查拼写


cat("绘制的目：", paste(orders_to_plot, collapse = ", "), "\n\n")

# 绘制函数
plot_order_accum_A <- function(ord) {
  df <- pollinator_by_order_A %>% 
    filter(Pollinator_order == ord) %>%
    mutate(Pollinator_accepted_name = forcats::fct_reorder(Pollinator_accepted_name, -Rank))
  
  tgt <- target_by_order_A %>% filter(Pollinator_order == ord)
  
  n_sp_total <- n_distinct(df$Pollinator_accepted_name)
  
  ggplot(df, aes(x = Rank, y = Cumulative_pct)) +
    geom_point(size = 1.5, shape = 1, color = "#FFB300", alpha = 0.8) +
    geom_line(size = 0.8, color = "#FFB300", alpha = 0.7) +
    geom_vline(xintercept = tgt$Rank, linetype = "dashed",
               color = "#27AE60", size = 1) +
    annotate("text",
             x = tgt$Rank,
             y = max(df$Cumulative_pct) * 0.85,
             label = paste0(tgt$Rank, " spp.\n(95%)"),
             color = "#27AE60", hjust = -0.1, fontface = "bold",
             size = 5) +
    labs(title = sprintf("%s (n=%d)", ord, n_sp_total), 
         x = "Species rank", 
         y = "Cumulative %") +
    scale_y_continuous(limits = c(0, 105)) +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 20),
      panel.grid = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(color = "black", size = 16),
      plot.margin = margin(10, 20, 10, 10)
    )
}

# 生成图列表
plots_list <- lapply(orders_to_plot, plot_order_accum_A)
supp_fig_A <- plot_grid(
  plotlist = plots_list,
  ncol = 2,
  labels = letters[seq_along(orders_to_plot)],
  label_size = 28
)

print(supp_fig_A)

# 保存图表
ggsave("result_260526/pollinator_accumulation_scheme_A_main.png",
       plot_A, width = 6, height = 4, dpi = 300)

ggsave("pollinator_accumulation_comparison.png",
       comparison_plot, width = 14, height = 6, dpi = 300)

ggsave("result_260526/pollinator_accumulation_by_order_A.png",
       supp_fig_A,
       width = 8,
       height = ceiling(length(orders_to_plot) / 2) * 3.5,
       dpi = 300)

cat("\n✅ 图表已保存")


#################### approach B: frequency
# Step 1: 对每个传粉昆虫统计总交互数 & 它的交互植物物种列表
pollinator_inter_times <- Interaction_data_pol_sp %>%
  filter(!Study_Network_id %in% crop_dom_site)%>% #exclude crop dominant site
  group_by(Pollinator_accepted_name) %>%
  summarize(
    Total_Interaction = sum(Interaction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Interaction)) %>%
  mutate(
    Rank = row_number(),  # 按交互次数排序
    Cumulative_interaction = cumsum(Total_Interaction)  # 累积交互次数
  )

# Step 2: 找到累积达到90%交互次数的点
target_y <- max(pollinator_inter_times$Cumulative_interaction) * 0.95
target_row <- pollinator_inter_times %>%
  filter(Cumulative_interaction >= target_y) %>%
  slice_min(Rank)

target_x <- target_row$Rank
target_y_value <- target_row$Cumulative_interaction

# Step 3: 绘图
accumulation_curve <- ggplot(pollinator_inter_times, aes(x = Rank, y = Cumulative_interaction)) +
  geom_point(size = 2, shape = 1, color = "#FFB300", alpha = 0.9) +
  geom_line(size = 0.7, color = "#FFB300") +
  geom_smooth(method = "loess", se = TRUE, linetype = "solid",
              fill = "#66C2A5", color = "#66C2A5", alpha = 0.2) +
  geom_vline(xintercept = target_x, linetype = "dashed", color = "#66C2A5", size = 1) +
  geom_hline(yintercept = target_y_value, linetype = "dashed", color = "#66C2A5") +
  geom_text(aes(x = target_x, y = target_y_value),
            label = paste0("95% of interactions captured\n", target_x, " species should be identified"),
            color = "#66C2A5", fontface = "bold", hjust = -0.1, vjust = 1.5, size = 5) +
  labs(
    x = "Number of Pollinator Species\n(Ranked by Interaction Frequency)",
    y = "Cumulative Interactions"
  ) +
  theme_bw(base_size=12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )
# Step 4: 展示图形
print(accumulation_curve)#600*370


ggsave("/Chap1_TargetPlant_to_monitor/result_260423/identify_sp.png",
       accumulation_curve, width = 5, height = 4, units = "in", dpi = 300)

############################ 分目展示
# ⭐ 需要先在 Interaction_data_pol_sp 里保留 Pollinator_order 信息
pollinator_inter_by_order <- Interaction_data_pol_sp %>%
  filter(!Study_Network_id %in% crop_dom_site) %>%
  filter(!is.na(Pollinator_order)) %>%
  group_by(Pollinator_order, Pollinator_accepted_name) %>%
  summarize(Total_Interaction = sum(Interaction, na.rm = TRUE), .groups = "drop") %>%
  group_by(Pollinator_order) %>%
  arrange(desc(Total_Interaction), .by_group = TRUE) %>%
  mutate(
    Rank = row_number(),
    Cumulative_interaction = cumsum(Total_Interaction),
    Pct = Cumulative_interaction / max(Cumulative_interaction)
  ) %>%
  ungroup()

# ⭐ 找每个目的95%截点
target_by_order <- pollinator_inter_by_order %>%
  group_by(Pollinator_order) %>%
  filter(Pct >= 0.95) %>%
  slice_min(Rank) %>%
  ungroup()

# 只保留数据量足够的目（至少10个物种）
# 手动指定目的顺序，替换原来的 orders_to_plot
orders_to_plot <- c("Diptera", "Lepidoptera", "Hymenoptera", "Coleoptera", "Hemiptera")

# 验证这5个目在数据里都存在
setdiff(orders_to_plot, unique(pollinator_inter_by_order$Pollinator_order))
# 如果输出为 character(0) 说明全部存在，否则检查拼写

plot_order_accum <- function(ord) {
  df <- pollinator_inter_by_order %>% filter(Pollinator_order == ord)
  tgt <- target_by_order %>% filter(Pollinator_order == ord)
  
  ggplot(df, aes(x = Rank, y = Cumulative_interaction)) +
    geom_line(color = "#FFB300", linewidth = 0.7) +
    geom_point(size = 1.5, shape = 1, color = "#FFB300", alpha = 0.9) +
    geom_vline(xintercept = tgt$Rank, linetype = "dashed",
               color = "#66C2A5", linewidth = 0.8) +
    annotate("text",
             x = tgt$Rank,
             y = max(df$Cumulative_interaction) * 0.85,  # ⭐ 从0.5→0.85，上移
             label = paste0(tgt$Rank, " spp.\n(95%)"),
             color = "#66C2A5", hjust = -0.1,
             size = 7) +                                 # ⭐ 从5→7
    labs(title = ord, x = "Species rank", y = "Cumulative interactions") +
    theme_bw(base_size = 22) +                           # ⭐ 从18→22
    theme(
      plot.title  = element_text(face = "bold.italic", size = 22),  # ⭐ 从14→22
      panel.grid  = element_blank(),
      axis.title  = element_text(size = 20),             # ⭐ 从13→20
      axis.text   = element_text(color = "black", size = 18),  # ⭐ 从11→18
      plot.margin = margin(10, 20, 10, 10)               # ⭐ 右边留更多空间给注释文字
    )
}

# 重新生成图列表
plots_list <- lapply(orders_to_plot, plot_order_accum)
supp_fig <- plot_grid(
  plotlist   = plots_list,
  ncol       = 2,
  labels     = c("a", "b", "c", "d", "e"),
  label_size = 30
)

ggsave("/Chap1_TargetPlant_to_monitor/result_260526/supp_accumulation_by_order.png",
       supp_fig,
       width  = 18,
       height = ceiling(5 / 2) * 6,   # ⭐ 写死5，ceil(5/2)=3行，高度=18
       dpi    = 300)

# ggsave("/Chap1_TargetPlant_to_monitor/result_260526/supp_accumulation_by_order.pdf",
#        supp_fig,
#        width  = 10,
#        height = ceiling(length(orders_to_plot) / 2) * 4)


###########################
#哪些种/属最需要被鉴定
# Step 1: 取出达到95%累积交互的前N种
library(dplyr)
library(flextable)
library(officer)

# ⭐ Step 1: 提取前N种并加上分类信息
pollinator_top_with_genus <- pollinator_inter_times %>%
  filter(Rank <= target_x) %>%
  mutate(Genus = sub(" .*", "", Pollinator_accepted_name)) %>%
  distinct(Pollinator_accepted_name, .keep_all = TRUE) %>%
  # 从原始数据里拿 Order 和 Family 信息
  left_join(
    Interaction_data_pol_sp %>%
      distinct(Pollinator_accepted_name, Pollinator_order, Pollinator_family),
    by = "Pollinator_accepted_name"
  ) %>%
  # ⭐ 按目、科、交互次数排序
  arrange(Pollinator_order, Pollinator_family, desc(Total_Interaction)) %>%
  select(
    Pollinator_order,
    Pollinator_family,
    Pollinator_accepted_name,
    Total_Interaction
  ) %>%
  rename(
    Order   = Pollinator_order,
    Family  = Pollinator_family,
    Species = Pollinator_accepted_name,
    `Total interactions` = Total_Interaction
  )

# ⭐ Step 2: 合并重复的 Order 和 Family 单元格（视觉分组）
# 找每个 Order 对应的行号
order_groups <- pollinator_top_with_genus %>%
  mutate(row_id = row_number()) %>%
  group_by(Order) %>%
  summarise(rows = list(row_id), .groups = "drop") %>%
  mutate(shade = row_number() %% 2 == 0)  # 奇偶交替

ft_table <- flextable(pollinator_top_with_genus) %>%
  merge_v(j = c("Order", "Family")) %>%
  theme_booktabs() %>%
  fontsize(size = 11, part = "all") %>%
  bold(part = "header") %>%
  italic(j = "Species", part = "body") %>%
  bold(j = "Order", part = "body") %>%
  valign(j = c("Order", "Family"), valign = "top", part = "body") %>%  # ⭐ 顶部对齐
  align(j = "Total interactions", align = "right", part = "all") %>%
  align(j = c("Order","Family","Species"), align = "left", part = "all") %>%
  width(j = "Order",              width = 1.3) %>%
  width(j = "Family",             width = 1.5) %>%
  width(j = "Species",            width = 3.0) %>%
  width(j = "Total interactions", width = 1.5) %>%
  # ⭐ 按 Order 交替背景色
  {
    ft <- .
    for (i in seq_len(nrow(order_groups))) {
      rows <- order_groups$rows[[i]]
      bg_color <- if (order_groups$shade[i]) "#F0F0F0" else "white"
      ft <- bg(ft, i = rows, bg = bg_color, part = "body")
    }
    ft
  } %>%
  set_caption(caption = paste0(
    "Table S1. The ", target_x,
    " pollinator species accounting for 95% of all recorded interactions, ",
    "organised by order and family."
  )) %>%
  fix_border_issues()

# ⭐ Step 3: 导出 Word
doc <- read_docx() %>%
  body_add_par("Table S1", style = "heading 1") %>%
  body_add_flextable(ft_table)

print(doc, target = "/Chap1_TargetPlant_to_monitor/result_260526/top_pollinator_species_table.docx")
message("✅ Word 表格已导出")