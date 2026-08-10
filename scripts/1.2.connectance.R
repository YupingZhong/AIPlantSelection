
############################################################################
# 3. Network Metrics 
# calculation of connectance and nestedness for complete networks and sub-networks
###################################################################################

# ==========================================================
# Required input files
# ==========================================================

# - data/processed/data_count_scaled_published.rds
# - data/processed/data_interact_published.rds
# - data/processed/Subsampled_networks.rds

# ==========================================================
# Load libraries
# ==========================================================
library(bipartite)
library(vegan)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)

# ==========================================================
# Load data
# ==========================================================
data_count_scaled <- readRDS("data/processed/data_count_scaled_published.rds")
data_interact <- readRDS("data/processed/data_interact_published.rds")
datasets <- readRDS("data/processed/Subsampled_networks.rds" )

# ==========================================================
# Network IDs
# ==========================================================

result_orig <- data_interact%>%
  filter(!Interaction_addup == 0)%>%
  group_by(Study_Network_id) %>%
  ungroup()

# ==========================================================
# Function 1: calculate connectance for one network
# ==========================================================

all_ids <- sort(unique(data_interact$Study_Network_id))

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

# ==========================================================
# Function 2: calculate connectance for all networks
# ==========================================================

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


connect_list <- purrr::imap(
  datasets,
  ~ connectance_df_from_data(.x, .y, all_ids))

# ==========================================================
# Combine results
# ==========================================================

all_connectance_df <- purrr::reduce(
  connect_list,
  dplyr::full_join,
  by = "Study_Network_id"
)
all_connectance_df <- all_connectance_df %>%
  dplyr::rename(
    Connectance_orig = orig,
    Connectance_Abun10 = Abun10,
    Connectance_Abun5 = Abun5,
    Connectance_Abun3 = Abun3,
    Connectance_FlwShape5 = FlwShape5,
    Connectance_FlwShape3 = FlwShape3,
    Connectance_Pylo10 = Pylo10,
    Connectance_Pylo5 = Pylo5,
    Connectance_Pylo3 = Pylo3)

# ==========================================================
# Save results
# ==========================================================

saveRDS(all_connectance_df, "data/processed/all_connectance_df.rds")

# ==========================================================
# Quality checks
# ==========================================================
# Ensure the number of networks matches
stopifnot(
  nrow(all_connectance_df) == length(all_ids)
)

# Ensure all network IDs match
stopifnot(
  setequal(
    all_connectance_df$Study_Network_id,
    all_ids
  )
)

# Ensure there are no duplicated network IDs
stopifnot(
  !anyDuplicated(all_connectance_df$Study_Network_id)
)

# Arrange rows in the same order as the ID skeleton
all_connectance_df <- all_connectance_df %>%
  arrange(
    match(Study_Network_id, all_ids)
  )

# Ensure connectance columns are numeric
all_connectance_df <- all_connectance_df %>%
  mutate(
    across(
      starts_with("Connectance"),
      as.numeric
    )
  )

# Summarise missing values
colSums(is.na(all_connectance_df))

# ==========================================================
# Calculate Spearman correlations
# ==========================================================
# all connectance metrics
connect_cols <- c(
  "Connectance_orig",
  "Connectance_Abun10",
  "Connectance_Abun5",
  "Connectance_Abun3",
  "Connectance_FlwShape5",
  "Connectance_FlwShape3",
  "Connectance_Pylo10",
  "Connectance_Pylo5",
  "Connectance_Pylo3")

# Pearson correlation with full network

cor_results2 <- lapply(
  connect_cols[-1],
  function(x){
    
    x_vec <- as.numeric(all_connectance_df[[x]])
    y_vec <- as.numeric(all_connectance_df$Connectance_orig)
    
    valid_idx <- is.finite(x_vec) &
      is.finite(y_vec)
    
    x_clean <- x_vec[valid_idx]
    y_clean <- y_vec[valid_idx]
    
    
    if(
      length(x_clean) < 3 ||
      length(unique(x_clean)) < 2
    ){
      
      return(
        c(
          r = NA,
          p_value = NA,
          n = length(x_clean)
        )
      )
      
    }
    
    
    test <- cor.test(
      x_clean,
      y_clean,
      method = "pearson"
    )
    
    
    c(
      r = unname(test$estimate),
      p_value = test$p.value,
      n = length(x_clean)
    )
    
  }
)



# dataframe

cor_results_df2 <- as.data.frame(
  do.call(rbind, cor_results2)
) %>%
  mutate(
    Comparison = connect_cols[-1],
    Metric = "Connectance",
    r = as.numeric(r),
    p_value = as.numeric(p_value),
    n = as.numeric(n)
  )


rownames(cor_results_df2) <- connect_cols[-1]


print(
  "Pearson correlation with full network:"
)

print(cor_results_df2)


# =========================
# Save results
# =========================

saveRDS(cor_results_df2, "data/processed/cor_results_df2.rds")


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
# saveRDS(random_all,"data/processed/random_all_metrics.rds")



















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


