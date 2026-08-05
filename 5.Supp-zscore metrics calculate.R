library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(bipartite)
library(ggplot2)
library(cowplot)
library(ggpubr)

# ---------------------------
# 0️⃣ 数据集合
# ---------------------------
datasets <- list(
  orig = data_merge,
  top10 = result10,
  top5 = result5,
  top3 = result3,
  op2_top5 = abun_species_top5_merge,
  op2_top3 = abun_species_top3_merge
)

# ---------------------------
# 1️⃣ 构建 network matrix
# ---------------------------
make_web_list <- function(data) {
  data %>%
    filter(!is.na(Pollinator_accepted_name)) %>%
    dplyr::select(Plant_accepted_name, Pollinator_accepted_name,
           Interaction_addup, Study_Network_id) %>%
    group_by(Study_Network_id) %>%
    group_split() %>%
    setNames(
      data %>%
        filter(!is.na(Pollinator_accepted_name)) %>%
        distinct(Study_Network_id) %>%
        pull(Study_Network_id)
    ) %>%
    lapply(function(df) {
      mat <- dcast(
        as.data.table(df),
        Plant_accepted_name ~ Pollinator_accepted_name,
        value.var = "Interaction_addup",
        fill = 0
      )
      mat <- as.matrix(mat[, -1])
      mat[is.na(mat)] <- 0
      mat <- mat[rowSums(mat) > 0, , drop = FALSE]
      mat <- mat[, colSums(mat) > 0, drop = FALSE]
      if (nrow(mat) < 2 || ncol(mat) < 2) return(NULL)
      return(mat)
    }) %>% purrr::compact()
}

# ---------------------------
# 2️⃣ observed metrics
# ---------------------------
calc_network_index <- function(webs, index_name) {
  lapply(webs, function(mat) bipartite::networklevel(mat, index = index_name))
}

# ---------------------------
# 3️⃣ null models
# ---------------------------
generate_nulls <- function(webs, method = "r2dtable", N = 500) {
  lapply(webs, function(mat) bipartite::nullmodel(mat, method = method, N = N))
}

# ---------------------------
# 4️⃣ null metric
# ---------------------------
calc_null_index <- function(nulls, index_name) {
  lapply(nulls, function(null_list) {
    sapply(null_list, function(mat) bipartite::networklevel(mat, index = index_name))
  })
}

# ---------------------------
# 5️⃣ z-score calculation
# ---------------------------
calc_zscore <- function(obs, nulls) {
  mapply(function(o, n) {
    n_vec <- unlist(n)
    if (length(n_vec) < 2 || sd(n_vec, na.rm = TRUE) == 0 || is.na(o)) return(NA)
    (o - mean(n_vec, na.rm = TRUE)) / sd(n_vec, na.rm = TRUE)
  }, obs, nulls, SIMPLIFY = FALSE)
}

# ---------------------------
# 6️⃣ z-score wrapper for all datasets
# ---------------------------
calc_zscore_all <- function(datasets, index_name) {
  results <- lapply(datasets, function(data) {
    webs <- make_web_list(data)
    if (length(webs) == 0) return(NULL)
    
    obs <- calc_network_index(webs, index_name)
    nulls <- generate_nulls(webs, N = 500)          # ✅ 设置 500 次
    null_index <- calc_null_index(nulls, index_name)
    z <- calc_zscore(obs, null_index)
    
    z_vec <- sapply(seq_along(webs), function(i) if (is.null(z[[i]]) || length(z[[i]]) == 0) NA else z[[i]])
    data.frame(Study_Network_id = names(webs), z_value = z_vec)
  })
  df <- reduce(results, left_join, by = "Study_Network_id")
  colnames(df) <- c("Study_Network_id", "z_orig", "z_top10", "z_top5", "z_top3", "z_op2_top5", "z_op2_top3")
  return(df)
}

# ---------------------------
# 7️⃣ 计算 nestedness & connectance
# ---------------------------
all_z_nested <- calc_zscore_all(datasets, index_name = "nestedness")
all_z_connectance <- calc_zscore_all(datasets, index_name = "connectance")

unique(abun_species_top3_merge$Study_Network_id)
saveRDS(all_z_nested,"all_z_nested.rds")
saveRDS(all_z_connectance,"all_z_connectance.rds")

# ---------------------------
# 8️⃣ 作图准备
# ---------------------------
my_colors <- c(
  "Top 10" = "#66C2A5", "Top 5" = "#4E79A7", "Top 3" = "#F28E2B",
  "flw 5" = "#4E79A7", "flw 3" = "#F28E2B"
)

# nestedness TopN
long_nested <- all_z_nested %>%
  pivot_longer(cols = c(z_top10, z_top5, z_top3), names_to = "Top_Plant_Level", values_to = "z_value") %>%
  filter(!is.na(z_value)) %>%
  mutate(Top_Plant_Level = factor(Top_Plant_Level, levels = c("z_top10","z_top5","z_top3"), labels = c("Top 10","Top 5","Top 3")))

# nestedness flw
long_nested_flw <- all_z_nested %>%
  pivot_longer(cols = c(z_op2_top5, z_op2_top3), names_to = "flw_Plant_Level", values_to = "z_value") %>%
  filter(!is.na(z_value)) %>%
  mutate(Top_Plant_Level = factor(flw_Plant_Level, levels = c("z_op2_top5","z_op2_top3"), labels = c("flw 5","flw 3")))

# connectance TopN
long_connect <- all_z_connectance %>%
  pivot_longer(cols = c(z_top10, z_top5, z_top3), names_to = "Top_Plant_Level", values_to = "z_value") %>%
  filter(!is.na(z_value)) %>%
  mutate(Top_Plant_Level = factor(Top_Plant_Level, levels = c("z_top10","z_top5","z_top3"), labels = c("Top 10","Top 5","Top 3")))

# connectance flw
long_connect_flw <- all_z_connectance %>%
  pivot_longer(cols = c(z_op2_top5, z_op2_top3), names_to = "flw_Plant_Level", values_to = "z_value") %>%
  filter(!is.na(z_value)) %>%
  mutate(Top_Plant_Level = factor(flw_Plant_Level, levels = c("z_op2_top5","z_op2_top3"), labels = c("flw 5","flw 3")))

# ---------------------------
# 9️⃣ 绘图函数
# ---------------------------
plot_z <- function(long_df, y_label) {
  reg_sign <- long_df %>%
    group_by(Top_Plant_Level) %>%
    summarise(lm_res = list(lm(z_value ~ z_orig, data = cur_data_all())), .groups = "drop") %>%
    rowwise() %>%
    mutate(p_value = summary(lm_res)$coefficients[2,4],
           signif_line = ifelse(p_value < 0.05, "solid", "dashed"))
  
  ggplot(long_df, aes(x = z_orig, y = z_value, color = Top_Plant_Level)) +
    geom_point(shape = 1, size = 2, position = position_jitter(width=0.05,height=0.05)) +
    geom_smooth(method = "lm", se = FALSE, aes(linetype=Top_Plant_Level), linewidth=0.8) +
    scale_color_manual(values = my_colors) +
    scale_linetype_manual(values = setNames(reg_sign$signif_line, reg_sign$Top_Plant_Level)) +
    ggpubr::stat_cor(aes(color=Top_Plant_Level), method="pearson",
                     label.x.npc="left", label.y.npc=seq(0.98,0.90,length.out=nlevels(long_df$Top_Plant_Level))) +
    labs(x = "Z-scored full network", y = y_label, color="Top N Level", linetype="Regression significance") +
    theme_classic(base_size=12)
}

# nestedness
p1 <- plot_z(long_nested, "Z-scored nestedness (subnetwork)")
p2 <- plot_z(long_nested_flw, "Z-scored nestedness (subnetwork)")
# connectance
p3 <- plot_z(long_connect, "Z-scored connectance (subnetwork)")
p4 <- plot_z(long_connect_flw, "Z-scored connectance (subnetwork)")

# ---------------------------
# 10️⃣ 拼图
# ---------------------------
p1.1 <- p1 + theme(legend.position = "none")
p2.1 <- p2 + theme(legend.position = "none")
p3.1 <- p3 + theme(legend.position = "none")
p4.1 <- p4 + theme(legend.position = "none")

combined_plot <- plot_grid(
  plot_grid(p1.1, p3.1, ncol=1),
  plot_grid(p2.1, p4.1, ncol=1),
  ncol=2
)

combined_plot


combined_plot#size:980*870
ggsave("/Chap1_TargetPlant_to_monitor/result_260423/zscore-network.matrix.png", combined_plot, width = 9.8, 
       height = 9.1, units = "in", dpi = 300)

