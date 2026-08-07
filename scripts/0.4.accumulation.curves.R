############################################################
# Optimising plant species selection for automated monitoring
# Networks overview
############################################################

# ==========================================================
# Load libraries
# ==========================================================

library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tidyr)
library(cowplot)

# ==========================================================
# Load data
# ==========================================================

data_count_scaled <- readRDS("data/processed/data_count_scaled_published.rds")#全部互作数据
data_interact <- readRDS("data/processed/data_interact_published.rds")

# ==========================================================
# Prepare data
# ==========================================================
# Merge data
data_merge <- merge(data_interact, data_count_scaled[,c("Flower_data_merger","Flower_count_scaled","Plant_species","Study_Network_id")], 
                    by = "Flower_data_merger",all = TRUE)%>%
  filter(!is.na(Flower_data_merger))%>%
  mutate(
    Study_Network_id = coalesce(Study_Network_id.x, Study_Network_id.y)) %>%
  select(-Study_Network_id.x, -Study_Network_id.y)%>%
  mutate(Interaction_addup = ifelse(is.na(Interaction_addup), 0, Interaction_addup))%>% # 替换 `Interaction_addup` 为 NA 的值为 0
  mutate(
    Plant_accepted_name = str_squish(str_replace_all(replace_na(Plant_accepted_name, ""), "×", ""))) 
# Rank plant species by abundance within each network
plant_rank_filt <- data_count_scaled %>%
  group_by(Study_Network_id, Plant_species) %>%
  summarise(
    Abundance = first(Flower_count_scaled),
    .groups = "drop") %>%
  group_by(Study_Network_id) %>%
  arrange(desc(Abundance), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()
# Total number of pollinator species per network
poll_pool <- data_interact %>%
  filter(!is.na(Pollinator_accepted_name)) %>%
  group_by(Study_Network_id) %>%
  summarise(
    total_poll = n_distinct(Pollinator_accepted_name),
    .groups = "drop")
# Total number of plant species per network
valid_nets <- plant_rank %>%
  group_by(Study_Network_id) %>%
  summarise(
    n_plants_total = n_distinct(Plant_species),
    .groups = "drop")
#List of poll species observed on each plant species per network
poll_map <- data_interact %>%
  group_by(Study_Network_id, Plant_species = Plant_original_name) %>%
  summarise(
    poll = list(unique(Pollinator_accepted_name)),
    .groups = "drop")

acc_curve <- plant_rank_filt %>%
  # Split dataset per network and calculate accumulation curve for each network
  group_by(Study_Network_id) %>%
  group_split() %>%
  map_dfr(function(net) {
    
    # Get network ID
    net_id <- unique(net$Study_Network_id)
    
    # Get plant species in order of abundance
    plants <- net %>%
      arrange(desc(Abundance)) %>%
      pull(Plant_species)
    
    # Get total number of pollinator species in the network
    total_poll <- poll_pool %>%
      filter(Study_Network_id == net_id) %>%
      pull(total_poll)
    
    # Empty vector to store pollinator species
    poll_seen <- character()
    
    # Empty list to store results
    results <- vector("list", length(plants))
    
    # Add plants sequentially and accumulate pollinator species
    for (k in seq_along(plants)) {
      
      new_poll <- poll_map %>%
        filter(
          Study_Network_id == net_id,
          Plant_species == plants[k]
        ) %>%
        pull(poll) %>%
        unlist() %>%
        unique()
      
      # Add newly encountered pollinators
      poll_seen <- union(poll_seen, new_poll)
      
      # Store result for this number of plants
      results[[k]] <- tibble(
        Study_Network_id = net_id,
        n_plants = k,
        pollinator_prop = length(poll_seen) / total_poll
      )
    }
    
    # Combine results for this network
    bind_rows(results)
  })

# Summary
summary_curve <- acc_curve %>%
  group_by(n_plants) %>%
  summarise(
    mean = mean(pollinator_prop, na.rm = TRUE),
    lower = quantile(pollinator_prop, 0.05, na.rm = TRUE),
    upper = quantile(pollinator_prop, 0.95, na.rm = TRUE),
    .groups = "drop")

# Plot accumulation curve
ggplot(summary_curve, aes(x = n_plants)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = mean), linewidth = 1.2) +
  labs(
    x = "Number of plant species included",
    y = "Proportion of pollinator richness captured") +
  theme_minimal()

# Plot individual accumulation curves for each network
ggplot(acc_curve,
       aes(n_plants, pollinator_prop)) +
  geom_line(
    aes(group = Study_Network_id),
    method = "loess",
    se = FALSE,
    span = 0.4,
    color = "grey55",
    alpha = 0.2,
    linewidth = 0.4) +
  theme_minimal()

# Check what percentage of sampling to get 80% of pollinator richness
# Minimum number of plants needed to reach 80% pollinator richness
threshold_80 <- acc_curve %>%
  filter(pollinator_prop >= 0.80) %>%
  group_by(Study_Network_id) %>%
  summarise(
    n_plants_80 = min(n_plants),
    .groups = "drop"
  ) %>%
  left_join(valid_nets, by = "Study_Network_id") %>%
  mutate(
    prop_plants = n_plants_80 / n_plants_total)

# Mean percentage of flowering species required
threshold_80 %>%
  summarise(
    mean_percent = mean(prop_plants) * 100,
    median_percent = median(prop_plants) * 100,
    sd_percent = sd(prop_plants) * 100)

# Random plant selection benchmark
# ==========================================================
# Create fast plant-pollinator lookup
# ==========================================================

poll_lookup <- poll_map %>%
  split(.$Study_Network_id) %>%
  map(function(x) {
    setNames(x$poll, x$Plant_species)
  })


# ==========================================================
# Random accumulation curve
# ==========================================================

random_K_curve <- function(net_id, n_perm = 200) {
  
  # Plant species available in this network
  net_plants <- plant_rank %>%
    filter(Study_Network_id == net_id) %>%
    pull(Plant_species) %>%
    unique()
  
  # Total pollinator richness
  total_poll <- poll_pool %>%
    filter(Study_Network_id == net_id) %>%
    pull(total_poll)
  
  # Plant-pollinator lookup for this network
  plant_polls <- poll_lookup[[net_id]]
  
  results <- vector(
    "list",
    n_perm
  )
  
  for (p in seq_len(n_perm)) {
    
    # Random order of plants
    sampled_order <- sample(net_plants)
    
    # Start with no pollinators
    poll_seen <- character()
    
    perm_results <- vector(
      "list",
      length(sampled_order)
    )
    
    for (k in seq_along(sampled_order)) {
      
      # Pollinators associated with the newly added plant
      new_poll <- plant_polls[[sampled_order[k]]]
      
      if (is.null(new_poll)) {
        new_poll <- character()
      }
      
      # Add only newly encountered pollinators
      poll_seen <- union(
        poll_seen,
        new_poll)
      
      perm_results[[k]] <- tibble(
        Study_Network_id = net_id,
        perm = p,
        n_plants = k,
        pollinator_prop = length(poll_seen) / total_poll)
    }
    
    results[[p]] <- bind_rows(
      perm_results)
  }
   bind_rows(results)}

# ==========================================================
# Run across networks
# ==========================================================

set.seed(123)

random_curve <- valid_nets$Study_Network_id %>%
  map_dfr(
    ~ random_K_curve(
      .x,
      n_perm = 100))
# Check proportion of pollinator richness captured for 10 plants
random_curve %>%
  filter(n_plants == 10) %>%
  summarise(
    mean_random = mean(pollinator_prop))

# Check proportion of pollinator richness captured for 10 plants
acc_curve %>%
  filter(n_plants == 10) %>%
  summarise(
    mean_abundance = mean(pollinator_prop))

# Summarise both curves for plotting
obs_summary <- acc_curve %>%
  group_by(n_plants) %>%
  summarise(
    mean = mean(pollinator_prop),
    lower = quantile(pollinator_prop, 0.05),
    upper = quantile(pollinator_prop, 0.95),
    .groups = "drop") %>%
  mutate(type = "abundance")

random_summary <- random_curve %>%
  group_by(n_plants) %>%
  summarise(
    mean = mean(pollinator_prop),
    lower = quantile(pollinator_prop, 0.05),
    upper = quantile(pollinator_prop, 0.95),
    .groups = "drop") %>%
  mutate(type = "random")

# Bind in a single datasets
plot_df <- bind_rows(obs_summary, random_summary)

# saveRDS(plot_df,"select_vs_radom_accum.rds")
# plot_df<-readRDS("select_vs_radom_accum.rds")

Poll_mean_20 <- ggplot(plot_df, aes(x = n_plants, y = mean, color = type, fill = type)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  scale_color_manual(values = c(
    abundance = "#1B9E77",
    random = "grey60"),
    labels = c(
      abundance = "Abundance-based",
      random = "Random")) +
  scale_fill_manual(values = c(
    abundance = "#1B9E77",
    random = "grey60"),
    labels = c(
      abundance = "Abundance-based",
      random = "Random")) +
  geom_vline(xintercept = c(3, 5, 10),
             linetype = "dashed",
             color = "black",
             alpha = 0.4) +
  geom_point(
    data = plot_df %>% filter(n_plants %in% c(3,5,10)),
    aes(color = type),
    size = 2)+
  annotate("text",
           x = c(3,5,10),
           y = Inf,
           label = c("3","5","10"),
           vjust = 1.5,
           size = 5) +
  labs(
    x = "Number of plant species sampled",
    y = "Percent of pollinator richness captured",
    color = "Plants selected",
    fill = "Plants selected")+
  scale_y_continuous(
    labels = function(x) x * 100,
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5))

Poll_mean_20
# Combine the three plots vertically

#ggsave("/Chap1_TargetPlant_to_monitor/result_260723/Poll_mean_20.png",Poll_mean_20,width = 9.5, height = 4.5
#       , units = "in", dpi = 300)

# histograms_combined <- plot_grid(histogram_1, histogram_2, ncol = 1, align = "v")
# 
# combined_plot <- plot_grid(histograms_combined, shape_barplot, ncol = 1, rel_heights = c(2, 1))

histograms_combined <- plot_grid(histogram_1, Poll_mean_20, ncol = 1, align = "v")
histograms_combined2 <- plot_grid(shape_barplot, shape_class, ncol = 1, align = "v")

histograms_combined3 <- plot_grid(histogram_1,shape_barplot,histogram_2,   ncol = 1, align = "v")

final_plot <- plot_grid(histograms_combined, histograms_combined2, ncol = 2, rel_widths = c(2.1, 2.1))

#ggsave("/Chap1_TargetPlant_to_monitor/result_260526/Fig1_four_penal.png", final_plot, width = 10.4, height = 6.7, units = "in", dpi = 300)

#ggsave("/Chap1_TargetPlant_to_monitor/result_260526/Distri_plant.png", histograms_combined3, width = 5, height = 6.5, units = "in", dpi = 300)
################3

#  spot check

################
#  
unique(acc_curve$Study_Network_id)

net_id <- "6_Marini_UNIPD03 CD26_2020"   # 改成你自己的ID

obs_net <- acc_curve %>%
  filter(Study_Network_id == net_id)

rand_net <- random_curve %>%
  filter(Study_Network_id == net_id)

obs_sum <- obs_net %>%
  group_by(n_plants) %>%
  summarise(
    mean = mean(pollinator_prop, na.rm = TRUE),
    lower = quantile(pollinator_prop, 0.05, na.rm = TRUE),
    upper = quantile(pollinator_prop, 0.95, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(type = "abundance")

rand_sum <- rand_net %>%
  group_by(n_plants) %>%
  summarise(
    mean = mean(pollinator_prop, na.rm = TRUE),
    lower = quantile(pollinator_prop, 0.05, na.rm = TRUE),
    upper = quantile(pollinator_prop, 0.95, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(type = "random")

plot_df <- bind_rows(obs_sum, rand_sum)

ggplot(plot_df, aes(x = n_plants, y = mean, color = type, fill = type)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    title = paste("Network:", net_id),
    x = "Number of plant species sampled",
    y = "Proportion of pollinator richness captured")

ggplot() +
  geom_line(data = rand_net,
            aes(n_plants, pollinator_prop, group = perm),
            alpha = 0.1, color = "grey60") +
  geom_line(data = obs_net,
            aes(n_plants, pollinator_prop),
            color = "red", linewidth = 1.2) +
  theme_minimal()

############## SES curve
random_stats <- random_curve %>%
  group_by(n_plants) %>%
  summarise(
    rand_mean = mean(pollinator_prop, na.rm = TRUE),
    rand_sd   = sd(pollinator_prop, na.rm = TRUE),
    .groups = "drop")

obs_stats <- acc_curve %>%
  group_by(n_plants) %>%
  summarise(
    obs_mean = mean(pollinator_prop, na.rm = TRUE),
    .groups = "drop")

ses_curve <- obs_stats %>%
  inner_join(random_stats, by = "n_plants") %>%
  mutate(
    SES = (obs_mean - rand_mean) / rand_sd)

sesplot <- ggplot(ses_curve, aes(x = n_plants, y = SES)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_line(linewidth = 1.2, color = "#2C7FB8") +
  geom_smooth(se = FALSE, color = "black", linewidth = 0.8) +
  theme_minimal(base_size = 13) +
  labs(
    x = "Number of plant species sampled (K)",
    y = "SES (Abundance-based vs Random)",
    title = "SES Accumulation Curve")

plot_combined <- plot_grid(accuplot, sesplot, ncol = 2, align = "v")

