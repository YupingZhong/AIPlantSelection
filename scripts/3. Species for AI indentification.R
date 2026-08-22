
###=================================================================================

### How many pollinator to identify to capture most of the interaction in Europe
##(Implication for classifier training)



library("dplyr")

####################################################
# data prepare
data_count_scaled<-readRDS("data/processed/data_count_scaled_published.rds")#全部互作数据
data_interact<-readRDS("data/processed/data_interact_published.rds")
data_merge<- readRDS("data/processed/data_merge.rds")

#======================================================================================

library(ggplot2)
library(dplyr)
library(purrr)
library(ggrepel)
#read data
metadata<-readRDS("data/raw/Interaction_data_published.rds")%>%
  mutate(Study_Network_id = paste(Study_id, Network_id, sep = "_"))
meta_count<-readRDS("data/raw/Flower_counts_published.rds")%>%
  mutate(Study_Network_id = paste(Study_id, Network_id, sep = "_"))
#首先删除那些pollinator只鉴定到属和大类的记录
Poll_sp<-unique(metadata$Pollinator_accepted_name)#2668 POLLINATORS 
Plant_sp<-unique(metadata$Plant_accepted_name)#1543 POLLINATORS

# #check the crop site
# flowercout_check<-meta_count %>%
#   group_by(Study_Network_id) %>%
#   mutate(sum_abun = sum(Flower_count),percent_abun = Flower_count/sum_abun)%>%
#   ungroup()%>%
#   filter(percent_abun>0.5 )
# 
# flowercoutsite<-flowercout_check%>%distinct(Study_Network_id,Plant_species)
# flowercout_check%>%distinct(Study_Network_id)
# flowercoutsite_list<-unique(flowercout_check$Plant_species)
# 
# #list argriculture plants with over 50% ralative abundance:
# #this step is to exclude crop dominant site
# 
# crop_species <- c(
#   "Medicago sativa",
#   "Helianthus annuus",
#   "Brassica napus",
#   "Vicia faba",
#   "Solanum tuberosum",
#   "Asparagus officinalis",
#   "Brassica oleracea"
# )
# 
# 
# # 1️⃣ site level filter
# site_info <- metadata %>%
#   group_by(Study_Network_id) %>%
#   summarise(
#     site_visit_time = sum(Interaction, na.rm = TRUE),
#     plant_sp = n_distinct(Plant_accepted_name)
#   ) %>%
#   filter(plant_sp < 10)
# 
# #  abundance level filter（只保留 crop）
# crop_dom <- meta_count %>%
#   filter(Plant_species %in% crop_species) %>%
#   group_by(Study_Network_id) %>%
#   mutate(
#     sum_abun = sum(Flower_count, na.rm = TRUE),
#     percent_abun = Flower_count / sum_abun
#   ) %>%
#   filter(percent_abun > 0.5) %>%
#   ungroup()
# 
# # 找交集 Study_Network_id
# crop_dom_site <- intersect(
#   site_info$Study_Network_id,
#   unique(crop_dom$Study_Network_id)
# )
# crop_dom_site
# 
# 
# dplyr::n_distinct(metadata$Study_Network_id)

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
  # filter(!Study_Network_id %in% crop_dom_site) %>%
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

top5_pollinators <- pollinator_importance_A %>%
  slice_head(n = 5)

# 数据中的传粉者物种总数
n_species_total <- nrow(pollinator_importance_A)

richness_colors <- c("Pollinator richness" = "#F1C453")
nestedness_color <- "#5A8BD4"

common_ylim <- c(0, 105)
common_breaks <- seq(0, 100, by = 20)
# ==============================================================================
# Top 5 pollinators for annotation
# ==============================================================================

top5_pollinators <- pollinator_importance_A %>%
  slice_head(n = 5) %>%
  mutate(
    label = paste0(
      "italic('",
      Pollinator_accepted_name,
      "')"
    )
  )
# ==============================================================================
# PANEL A
# Number of pollinators needed to capture most interactions
# ==============================================================================

p_capture <- ggplot(
  
  pollinator_importance_A,
  
  aes(
    x = Rank,
    y = Cumulative_pct
  )
  
) +
  
  # ----------------------------------------------------------------------------
# Cumulative interaction curve
# ----------------------------------------------------------------------------

geom_line(
  
  linewidth = 0.9,
  
  color =
    richness_colors["Pollinator richness"],
  
  alpha = 0.85
  
) +
  
  geom_point(
    
    shape = 16,
    
    size = 2.0,
    
    alpha = 0.65,
    
    color =
      richness_colors["Pollinator richness"]
    
  ) +
  
  # ----------------------------------------------------------------------------
# Highlight Top 5 species
# ----------------------------------------------------------------------------

geom_point(
  
  data =
    top5_pollinators,
  
  shape = 16,
  
  size = 3.0,
  
  color =
    nestedness_color
  
) +
  
  # ----------------------------------------------------------------------------
# Top 5 species labels
# ----------------------------------------------------------------------------

geom_text_repel(
  
  data = top5_pollinators,
  
  aes(
    label = label
  ),
  
  parse = TRUE,
  
  size = 3.8,
  
  color = "black",
  
  box.padding = 0.4,
  
  point.padding = 0.3,
  
  force = 1,
  
  min.segment.length = 0,
  
  segment.color = "grey50",
  
  direction = "y",
  
  seed = 123
  
)+
  
  # ----------------------------------------------------------------------------
# 95% reference lines
# ----------------------------------------------------------------------------

geom_vline(
  
  xintercept =
    target_x_A,
  
  linetype =
    "dashed",
  
  linewidth =
    0.8,
  
  color =
    nestedness_color,
  
  alpha =
    0.75
  
) +
  
  geom_hline(
    
    yintercept =
      target_pct_A,
    
    linetype =
      "dashed",
    
    linewidth =
      0.8,
    
    color =
      nestedness_color,
    
    alpha =
      0.75
    
  ) +
  
  # ----------------------------------------------------------------------------
# 95% annotation
# ----------------------------------------------------------------------------

annotate(
  
  "text",
  
  x =
    target_x_A +
    0.03 * max(pollinator_importance_A$Rank),
  
  y =
    85,
  
  label =
    paste0(
      target_pct_A,
      "% captured\n",
      target_x_A,
      " species (",
      sprintf(
        "%.1f%%",
        target_x_A /
          nrow(pollinator_importance_A) *
          100
      ),
      " of total)"
    ),
  
  size =
    4.0,
  
  fontface =
    "plain",
  
  color =
    nestedness_color,
  
  hjust =
    0,
  
  vjust =
    0.5
  
) +
  
  # ----------------------------------------------------------------------------
# X-axis
# ----------------------------------------------------------------------------

scale_x_continuous(
  
  expand =
    expansion(
      mult = c(
        0.01,
        0.02
      )
    )
  
) +
  
  # ----------------------------------------------------------------------------
# Y-axis
# ----------------------------------------------------------------------------

scale_y_continuous(
  
  limits =
    common_ylim,
  
  breaks =
    common_breaks,
  
  expand =
    expansion(
      mult = c(
        0,
        0.02
      )
    )
  
) +
  
  # ----------------------------------------------------------------------------
# Labels
# ----------------------------------------------------------------------------

labs(
  
  x =
    "Pollinator species identified",
  
  y =
    "Cumulative relative interactions (%)"
  
) +
  
  # ----------------------------------------------------------------------------
# Theme
# ----------------------------------------------------------------------------

theme_classic(
  
  base_size = 12
  
) +
  
  theme(
    
    axis.title.x =
      element_text(
        face = "plain",
        size = 13
      ),
    
    axis.title.y =
      element_text(
        face = "plain",
        size = 13
      ),
    
    axis.text.x =
      element_text(
        color = "black",
        size = 11
      ),
    
    axis.text.y =
      element_text(
        color = "black",
        size = 11
      ),
    
    legend.position =
      "none",
    
    plot.margin =
      margin(
        5,
        5,
        5,
        5
      )
    
  )

print(p_capture)


ggsave(
  "/Chap1_TargetPlant_to_monitor/result_260723/Pollinator_importance_curve.png",
  p_capture,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

#======================================================================================
# 🟢 分目展示：用方案A（推荐）
#======================================================================================

cat("\n\n========== 按目（Order）分组分析 - 方案A ==========\n\n")

pollinator_by_order_A <- Interaction_data_pol_sp %>%
  # filter(!Study_Network_id %in% crop_dom_site) %>%
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

library(cowplot)

curve_color <- "#F1C453"       # 金黄色：累计曲线
highlight_color <- "#5A8BD4"   # 蓝色：95% 阈值与文字

plot_order_accum_A <- function(ord, show_y_title = TRUE) {
  
  df <- pollinator_by_order_A %>%
    filter(Pollinator_order == ord)
  
  tgt <- target_by_order_A %>%
    filter(Pollinator_order == ord)
  
  n_sp_total <- n_distinct(df$Pollinator_accepted_name)
  
  ggplot(df, aes(x = Rank, y = Cumulative_pct)) +
    geom_point(
      size = 1.5, shape = 1,
      colour = curve_color, alpha = 0.8
    ) +
    geom_line(
      linewidth = 0.8,
      colour = curve_color, alpha = 0.7
    ) +
    geom_vline(
      xintercept = tgt$Rank,
      linetype = "dashed",
      colour = highlight_color,
      linewidth = 0.8,
      alpha = 0.75
    ) +
    geom_hline(
      yintercept = 95,
      linetype = "dashed",
      colour = highlight_color,
      linewidth = 0.8,
      alpha = 0.75
    ) +
    annotate(
      "text",
      x = tgt$Rank,
      y = 85,
      label = paste0(tgt$Rank, " spp.\n(95%)"),
      colour = highlight_color,
      hjust = -0.1,
      fontface = "plain",
      size = 4.5
    ) +
    labs(
      title = sprintf("%s (n = %d)", ord, n_sp_total),
      x = "Pollinator species identified",
      y = if (show_y_title) {
        "Cumulative relative interactions (%)"
      } else {
        NULL
      }
    ) +
    scale_y_continuous(
      limits = c(0, 105),
      breaks = seq(0, 100, by = 20)
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 13),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(colour = "black", size = 11),
      plot.margin = margin(5, 12, 5, 5)
    )
}

plots_list <- lapply(seq_along(orders_to_plot), function(i) {
  plot_order_accum_A(
    ord = orders_to_plot[i],
    show_y_title = i %% 3 == 1
  )
})

supp_fig_A <- plot_grid(
  plotlist = plots_list,
  ncol = 3,
  nrow = 2,
  labels = letters[seq_along(orders_to_plot)],
  label_size = 16
)

print(supp_fig_A)

ggsave(
  "/Chap1_TargetPlant_to_monitor/result_260723/pollinator_accumulation_by_order_A.png",
  supp_fig_A,
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

###########################
#哪些种/属最需要被鉴定
# Step 1: 取出达到95%累积交互的前N种
library(dplyr)
library(flextable)
library(officer)

# Step 1：提取达到 95% 累积平均相对互作量所需的物种
pollinator_top_with_genus <- pollinator_importance_A %>%
  filter(Rank <= target_x_A) %>%
  left_join(
    Interaction_data_pol_sp %>%
      # filter(!Study_Network_id %in% crop_dom_site) %>%
      distinct(
        Pollinator_accepted_name,
        Pollinator_order,
        Pollinator_family
      ),
    by = "Pollinator_accepted_name"
  ) %>%
  mutate(
    `Mean relative interactions (%)` =
      round(Mean_relative_abundance * 100, 3)
  ) %>%
  arrange(
    Pollinator_order,
    Pollinator_family,
    desc(Mean_relative_abundance)
  ) %>%
  dplyr::select(
    Pollinator_order,
    Pollinator_family,
    Pollinator_accepted_name,
    `Mean relative interactions (%)`,
    n_networks_present
  ) %>%
  rename(
    Order = Pollinator_order,
    Family = Pollinator_family,
    Species = Pollinator_accepted_name,
    `Networks present` = n_networks_present
  )

# Step 2：确定各 Order 的行，用于交替底色
order_groups <- pollinator_top_with_genus %>%
  mutate(row_id = row_number()) %>%
  group_by(Order) %>%
  summarise(rows = list(row_id), .groups = "drop") %>%
  mutate(shade = row_number() %% 2 == 0)

# Step 3：制作表格
ft_table <- flextable(pollinator_top_with_genus) %>%
  merge_v(j = c("Order", "Family")) %>%
  theme_booktabs() %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  italic(j = "Species", part = "body") %>%
  bold(j = "Order", part = "body") %>%
  valign(
    j = c("Order", "Family"),
    valign = "top",
    part = "body"
  ) %>%
  align(
    j = c("Mean relative interactions (%)", "Networks present"),
    align = "right",
    part = "all"
  ) %>%
  align(
    j = c("Order", "Family", "Species"),
    align = "left",
    part = "all"
  ) %>%
  set_header_labels(
    Order = "Order",
    Family = "Family",
    Species = "Species",
    `Mean relative interactions (%)` =
      "Mean relative\ninteractions (%)",
    `Networks present` = "Networks\npresent"
  ) %>%
  width(j = "Order", width = 1.15) %>%
  width(j = "Family", width = 1.35) %>%
  width(j = "Species", width = 2.5) %>%
  width(j = "Mean relative interactions (%)", width = 1.35) %>%
  width(j = "Networks present", width = 1.05) %>%
  {
    ft <- .
    for (i in seq_len(nrow(order_groups))) {
      ft <- bg(
        ft,
        i = order_groups$rows[[i]],
        bg = if (order_groups$shade[i]) "#F0F0F0" else "white",
        part = "body"
      )
    }
    ft
  } %>%
  set_caption(
    caption = paste0(
      "Table S1. The ", target_x_A,
      " pollinator species accounting for 95% of the mean relative interactions across networks, ",
      "organised by order and family."
    )
  ) %>%
  fix_border_issues() %>%
  autofit() %>%
  fit_to_width(max_width = 6.5)

# Step 4：导出 Word
doc <- read_docx() %>%
  body_add_flextable(ft_table)

print(
  doc,
  target = "/Chap1_TargetPlant_to_monitor/result_260723/top_pollinator_species_table.docx"
)


write.csv(
  pollinator_top_with_genus,
  "/Chap1_TargetPlant_to_monitor/result_260723/Supplementary_Dataset_S1_top_pollinator_species.csv",
  row.names = FALSE
)


###########short version
table_S1_summary <- pollinator_top_with_genus %>%
  group_by(Order, Family) %>%
  summarise(
    `Species required` = n(),
    `Mean relative interactions (%)` =
      round(sum(`Mean relative interactions (%)`), 2),
    .groups = "drop"
  ) %>%
  group_by(Order) %>%
  mutate(
    `Species required (%)` = round(
      `Species required` / sum(`Species required`) * 100,
      1
    )
  ) %>%
  ungroup() %>%
  arrange(Order, desc(`Species required`))

ft_summary <- flextable(table_S1_summary) %>%
  theme_booktabs() %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  merge_v(j = "Order") %>%
  valign(j = "Order", valign = "top", part = "body") %>%
  autofit() %>%
  set_caption(
    caption = paste0(
      "Table S1. Taxonomic composition of the ", target_x_A,
      " pollinator species required to account for 95% of mean relative interactions."
    )
  )
# Step 4：导出 Word
doc_ft_summary <- read_docx() %>%
  body_add_flextable(ft_summary)

print(
  doc_ft_summary,
  target = "/Chap1_TargetPlant_to_monitor/result_260723/top_pollinator_species_table.docx"
)



# #################### approach B: frequency
# # Step 1: 对每个传粉昆虫统计总交互数 & 它的交互植物物种列表
# pollinator_inter_times <- Interaction_data_pol_sp %>%
#   filter(!Study_Network_id %in% crop_dom_site)%>% #exclude crop dominant site
#   group_by(Pollinator_accepted_name) %>%
#   summarize(
#     Total_Interaction = sum(Interaction, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(Total_Interaction)) %>%
#   mutate(
#     Rank = row_number(),  # 按交互次数排序
#     Cumulative_interaction = cumsum(Total_Interaction)  # 累积交互次数
#   )
# 
# # Step 2: 找到累积达到90%交互次数的点
# target_y <- max(pollinator_inter_times$Cumulative_interaction) * 0.95
# target_row <- pollinator_inter_times %>%
#   filter(Cumulative_interaction >= target_y) %>%
#   slice_min(Rank)
# 
# target_x <- target_row$Rank
# target_y_value <- target_row$Cumulative_interaction
# 
# # Step 3: 绘图
# accumulation_curve <- ggplot(pollinator_inter_times, aes(x = Rank, y = Cumulative_interaction)) +
#   geom_point(size = 2, shape = 1, color = "#FFB300", alpha = 0.9) +
#   geom_line(size = 0.7, color = "#FFB300") +
#   geom_smooth(method = "loess", se = TRUE, linetype = "solid",
#               fill = "#66C2A5", color = "#66C2A5", alpha = 0.2) +
#   geom_vline(xintercept = target_x, linetype = "dashed", color = "#66C2A5", size = 1) +
#   geom_hline(yintercept = target_y_value, linetype = "dashed", color = "#66C2A5") +
#   geom_text(aes(x = target_x, y = target_y_value),
#             label = paste0("95% of interactions captured\n", target_x, " species should be identified"),
#             color = "#66C2A5", fontface = "bold", hjust = -0.1, vjust = 1.5, size = 5) +
#   labs(
#     x = "Number of Pollinator Species\n(Ranked by Interaction Frequency)",
#     y = "Cumulative Interactions"
#   ) +
#   theme_bw(base_size=12) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
#     axis.title = element_text(face = "bold", size = 15),
#     axis.text = element_text(color = "black", size = 12),
#     panel.background = element_blank(),
#     panel.grid = element_blank()
#   )
# # Step 4: 展示图形
# print(accumulation_curve)#600*370
# 
# 
# ggsave("/Chap1_TargetPlant_to_monitor/result_260423/identify_sp.png",
#        accumulation_curve, width = 5, height = 4, units = "in", dpi = 300)

# ############################ interaction frequency 分目展示
# # ⭐ 需要先在 Interaction_data_pol_sp 里保留 Pollinator_order 信息
# pollinator_inter_by_order <- Interaction_data_pol_sp %>%
#   filter(!Study_Network_id %in% crop_dom_site) %>%
#   filter(!is.na(Pollinator_order)) %>%
#   group_by(Pollinator_order, Pollinator_accepted_name) %>%
#   summarize(Total_Interaction = sum(Interaction, na.rm = TRUE), .groups = "drop") %>%
#   group_by(Pollinator_order) %>%
#   arrange(desc(Total_Interaction), .by_group = TRUE) %>%
#   mutate(
#     Rank = row_number(),
#     Cumulative_interaction = cumsum(Total_Interaction),
#     Pct = Cumulative_interaction / max(Cumulative_interaction)
#   ) %>%
#   ungroup()
# 
# # ⭐ 找每个目的95%截点
# target_by_order <- pollinator_inter_by_order %>%
#   group_by(Pollinator_order) %>%
#   filter(Pct >= 0.95) %>%
#   slice_min(Rank) %>%
#   ungroup()
# 
# # 只保留数据量足够的目（至少10个物种）
# # 手动指定目的顺序，替换原来的 orders_to_plot
# orders_to_plot <- c("Diptera", "Lepidoptera", "Hymenoptera", "Coleoptera", "Hemiptera")
# 
# # 验证这5个目在数据里都存在
# setdiff(orders_to_plot, unique(pollinator_inter_by_order$Pollinator_order))
# # 如果输出为 character(0) 说明全部存在，否则检查拼写
# 
# library(cowplot)
# 
# curve_color <- "#F1C453"       # 金黄色：累计曲线
# highlight_color <- "#5A8BD4"   # 蓝色：95% 阈值与文字
# 
# plot_order_accum <- function(ord) {
#   
#   df <- pollinator_inter_by_order %>%
#     filter(Pollinator_order == ord)
#   
#   tgt <- target_by_order %>%
#     filter(Pollinator_order == ord)
#   
#   n_sp_total <- n_distinct(df$Pollinator_accepted_name)
#   
#   ggplot(df, aes(x = Rank, y = Cumulative_interaction)) +
#     geom_point(
#       size = 1.5, shape = 1,
#       colour = curve_color, alpha = 0.8
#     ) +
#     geom_line(
#       linewidth = 0.8,
#       colour = curve_color, alpha = 0.7
#     ) +
#     geom_vline(
#       xintercept = tgt$Rank,
#       linetype = "dashed",
#       colour = highlight_color,
#       linewidth = 0.8,
#       alpha = 0.75
#     ) +
#     geom_hline(
#       yintercept = tgt$Cumulative_interaction,
#       linetype = "dashed",
#       colour = highlight_color,
#       linewidth = 0.8,
#       alpha = 0.75
#     ) +
#     annotate(
#       "text",
#       x = tgt$Rank,
#       y = max(df$Cumulative_interaction) * 0.85,
#       label = paste0(tgt$Rank, " spp.\n(95%)"),
#       colour = highlight_color,
#       hjust = -0.1,
#       size = 4.5
#     ) +
#     labs(
#       title = sprintf("%s (n = %d)", ord, n_sp_total),
#       x = "Pollinator species identified",
#       y = "Cumulative interactions"
#     ) +
#     scale_x_continuous(
#       expand = expansion(mult = c(0.02, 0.15))
#     ) +
#     theme_classic(base_size = 12) +
#     theme(
#       plot.title = element_text(face = "bold.italic", size = 13),
#       axis.title = element_text(size = 12),
#       axis.text = element_text(colour = "black", size = 11),
#       plot.margin = margin(5, 12, 5, 5)
#     )
# }
# 
# plots_list <- lapply(orders_to_plot, plot_order_accum)
# 
# supp_fig <- plot_grid(
#   plotlist = plots_list,
#   ncol = 3,
#   nrow = 2,
#   labels = letters[seq_along(orders_to_plot)],
#   label_size = 16
# )
# 
# print(supp_fig)
# 
# ggsave(
#   "/Chap1_TargetPlant_to_monitor/result_260723/supp_accumulation_by_order.png",
#   supp_fig,
#   width = 15,
#   height = 8,
#   units = "in",
#   dpi = 600,
#   bg = "white"
# )
# 

# ggsave("/Chap1_TargetPlant_to_monitor/result_260526/supp_accumulation_by_order.pdf",
#        supp_fig,
#        width  = 10,
#        height = ceiling(length(orders_to_plot) / 2) * 4)

