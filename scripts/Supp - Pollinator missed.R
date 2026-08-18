#####################
##Supp  pollinator Taxon/species constantly missed if only monitor top 10 abundant plants


library("dplyr")

####################################################
# data prepare
data_count_scaled<-readRDS("data/processed/data_count_scaled_published.rds")#全部互作数据
data_interact<-readRDS("data/processed/data_interact_published.rds")
data_merge<- readRDS("data/processed/data_merge.rds")

unique(data_interact$Study_Network_id)
unique(data_count_scaled_species$Study_Network_id)
unique(data_merge$Study_Network_id)


###===============================================================================

# ------------1. Pollinator Species constantly Missed by Abundant 10 subsampling strategy

###=====================================================================================
library(stringr)
library(tidyr)
library(dplyr)
result10<-readRDS("data/processed/selected_plant_result_abun10.rds")
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
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      Missed = "#E76F51",
      Captured = "#a1d8e8"
    ),
    breaks = c("Missed", "Captured"),
    labels = c(
      Missed = "Missed networks",
      Captured = "Captured networks"
    )
  ) +
  facet_wrap(
    ~ Pollinator_order,
    nrow = 2,
    scales = "free_x"
  ) +
  labs(
    x = NULL,
    y = "Number of networks",
    fill = NULL
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(2, "pt"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(6, 8, 6, 6)
  )

p

ggsave(
  "/Chap1_TargetPlant_to_monitor/result_260723/order_Pollinator_missed_bar.png",
  p, width = 9, height = 8, dpi = 300, bg = "white"
)



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


