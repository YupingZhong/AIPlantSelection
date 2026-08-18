############################################################
# Optimising plant species selection for automated monitoring
# Analysis by POLLINATOR ORDER
############################################################

# ==========================================================
# Load libraries and data
# ==========================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(vegan)
library(V.PhyloMaker2)
library(ape)

# Load the source data
data_count_scaled <- readRDS("data/processed/data_count_scaled_published.rds")
data_interact <- readRDS("data/processed/data_interact_published.rds")
traits <- read.csv("data/processed/merge.trait.csv", header = TRUE, fileEncoding = "UTF-8")

# ==========================================================
# Data structure verification
# ==========================================================
# Expected columns in data_interact:
# [31] "Study_Network_id"            
# [23] "Pollinator_order"            
# [33] "Flower_data_merger"         
# [36] "flw_shape_revised"           

colnames(data_interact)
unique(data_interact$Pollinator_order)

# Phylogenetic tree (pre-built from main script)
data("GBOTB.extended.TPL")

# ==========================================================
# MAIN FUNCTION: Calculate all metrics for one pollinator order
# ==========================================================

analyze_pollinator_order <- function(pollinator_order_name) {
  
  cat("\n=== Processing:", pollinator_order_name, "===\n")
  
  # =========================================================
  # Step 1: Filter data by pollinator order
  # =========================================================
  
  data_interact_filt <- data_interact %>%
    filter(Pollinator_order == pollinator_order_name)
  
  # Get networks that have this pollinator group
  networks_with_order <- unique(data_interact_filt$Study_Network_id)
  
  cat("Networks with", pollinator_order_name, ":", 
      length(networks_with_order), "\n")
  
  if (length(networks_with_order) == 0) {
    cat("WARNING: No networks for", pollinator_order_name, "\n")
    return(NULL)
  }
  
  # =========================================================
  # Step 2: Calculate total pollinator richness per network
  # =========================================================
  
  total_number <- data_interact_filt %>%
    group_by(Study_Network_id) %>%
    filter(!Interaction_addup == 0) %>%
    summarise(
      total_pollinator_count = n_distinct(
        Pollinator_accepted_name, na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  # =========================================================
  # Option 1: Abundance-based selection (Top 10, 5, 3)
  # =========================================================
  
  calc_abundance_strategy <- function(n_top, strategy_name) {
    
    top_species <- data_count_scaled %>%
      distinct(Study_Network_id, Plant_species, 
               Flower_count_scaled, .keep_all = FALSE) %>%
      group_by(Study_Network_id) %>%
      slice_max(order_by = Flower_count_scaled, n = n_top, 
                with_ties = FALSE) %>%
      ungroup() %>%
      group_by(Study_Network_id) %>%
      filter(n_distinct(Plant_species) > (n_top - 1)) %>%
      ungroup() %>%
      select(Study_Network_id, Plant_species)
    
    result <- left_join(
      top_species,
      data_interact_filt,
      by = c("Plant_species" = "Plant_original_name", 
             "Study_Network_id" = "Study_Network_id")
    ) %>%
      filter(!Interaction_addup == 0) %>%
      group_by(Study_Network_id, Pollinator_accepted_name) %>%
      summarise(freq = sum(Interaction_addup, na.rm = TRUE), 
                .groups = "drop") %>%
      group_by(Study_Network_id) %>%
      summarise(
        pollinator_count = n_distinct(Pollinator_accepted_name, na.rm = TRUE),
        .groups = "drop"
      )
    
    result %>%
      left_join(total_number, by = "Study_Network_id") %>%
      replace_na(list(pollinator_count = 0)) %>%
      mutate(
        percentage = (pollinator_count / total_pollinator_count) * 100,
        strategy = strategy_name,
        n_plants = n_top
      ) %>%
      select(Study_Network_id, pollinator_count, total_pollinator_count, 
             percentage, strategy, n_plants)
  }
  
  abun_10 <- calc_abundance_strategy(10, "Abundance")
  abun_5  <- calc_abundance_strategy(5, "Abundance")
  abun_3  <- calc_abundance_strategy(3, "Abundance")
  
  # =========================================================
  # Option 2: Flower shape-based selection (Top 5, 3)
  # =========================================================
  
  calc_shape_strategy <- function(n_top, strategy_name) {
    
    abun_species <- data_count_scaled %>%
      group_by(Study_Network_id, flw_shape_revised) %>%
      arrange(Study_Network_id, flw_shape_revised, 
              desc(Flower_count_scaled)) %>%
      top_n(1, wt = Flower_count_scaled) %>%
      distinct(Study_Network_id, Plant_species, Flower_count_scaled) %>%
      group_by(Study_Network_id) %>%
      mutate(rank = rank(desc(Flower_count_scaled), ties.method = "first")) %>%
      filter(rank < (n_top + 1)) %>%
      ungroup() %>%
      group_by(Study_Network_id) %>%
      filter(n_distinct(Plant_species) > (n_top - 1)) %>%
      ungroup() %>%
      select(Plant_species, Study_Network_id)
    
    result <- left_join(
      abun_species,
      data_interact_filt,
      by = c("Plant_species" = "Plant_original_name", 
             "Study_Network_id" = "Study_Network_id")
    ) %>%
      filter(Interaction_addup > 0) %>%
      group_by(Study_Network_id, Pollinator_accepted_name) %>%
      summarise(freq = sum(Interaction_addup, na.rm = TRUE), 
                .groups = "drop") %>%
      group_by(Study_Network_id) %>%
      summarise(
        pollinator_count = n_distinct(Pollinator_accepted_name, na.rm = TRUE),
        .groups = "drop"
      )
    
    result %>%
      left_join(total_number, by = "Study_Network_id") %>%
      replace_na(list(pollinator_count = 0)) %>%
      mutate(
        percentage = (pollinator_count / total_pollinator_count) * 100,
        strategy = strategy_name,
        n_plants = n_top
      ) %>%
      select(Study_Network_id, pollinator_count, total_pollinator_count, 
             percentage, strategy, n_plants)
  }
  
  shape_5 <- calc_shape_strategy(5, "Flower shape")
  shape_3 <- calc_shape_strategy(3, "Flower shape")
  
  # =========================================================
  # Option 3: Phylogenetic diversity selection (Top 10, 5, 3)
  # =========================================================
  
  calc_phylo_strategy <- function(n_top, strategy_name) {
    
    # Build phylogenetic tree
    plant_sp <- data_count_scaled %>%
      distinct(Plant_species) %>%
      na.omit() %>%
      filter(
        str_count(Plant_species, "\\S+") >= 2,
        !grepl("\\bsp\\.?\\b|cf\\.|aff\\.", Plant_species, 
               ignore.case = TRUE)
      ) %>%
      rename(species = Plant_species) %>%
      mutate(species = gsub(" ", "_", species)) %>%
      left_join(
        tips.info.TPL[, c("species", "genus", "family")],
        by = "species"
      )
    
    phylo_result <- phylo.maker(
      sp.list = plant_sp,
      tree = GBOTB.extended.TPL,
      scenarios = "S3"
    )
    phylo_tree <- phylo_result$scenario.3
    
    # Prepare data
    data_PD <- data_count_scaled %>%
      filter(
        !is.na(Plant_species),
        str_count(Plant_species, "\\S+") >= 2,
        !grepl("\\bsp\\.?\\b|cf\\.|aff\\.", Plant_species, 
               ignore.case = TRUE)
      ) %>%
      mutate(Plant_species_tree = gsub(" ", "_", Plant_species))
    
    data_PD_tree <- data_PD %>%
      filter(Plant_species_tree %in% phylo_tree$tip.label)
    
    # Selection function
    select_PD_species <- function(species_pool, phylo_tree, n_select){
      
      species_pool <- intersect(species_pool, phylo_tree$tip.label)
      
      if(length(species_pool) < n_select){
        return(NA)
      }
      
      selected <- sample(species_pool, 1)
      
      while(length(selected) < n_select){
        candidates <- setdiff(species_pool, selected)
        pd_gain <- sapply(candidates, function(x){
          tips <- c(selected, x)
          sub_tree <- drop.tip(phylo_tree, 
                               setdiff(phylo_tree$tip.label, tips))
          sum(sub_tree$edge.length)
        })
        selected <- c(selected, candidates[which.max(pd_gain)])
      }
      selected
    }
    
    # Select species
    PD_selected <- data_PD_tree %>%
      group_by(Study_Network_id) %>%
      summarise(
        Plant_species_tree = list(
          select_PD_species(Plant_species_tree, phylo_tree, n_top)
        ),
        .groups = "drop"
      ) %>%
      filter(!is.na(Plant_species_tree)) %>%
      unnest(Plant_species_tree) %>%
      mutate(Plant_species = gsub("_", " ", Plant_species_tree)) %>%
      select(Plant_species, Study_Network_id)
    
    # Calculate capture
    result <- left_join(
      PD_selected,
      data_interact_filt,
      by = c("Plant_species" = "Plant_original_name", 
             "Study_Network_id" = "Study_Network_id")
    ) %>%
      group_by(Study_Network_id) %>%
      summarise(
        pollinator_count = n_distinct(
          Pollinator_accepted_name[Interaction_addup > 0 & 
                                     !is.na(Pollinator_accepted_name)]
        ),
        .groups = "drop"
      ) %>%
      mutate(pollinator_count = replace_na(pollinator_count, 0))
    
    result %>%
      left_join(total_number, by = "Study_Network_id") %>%
      mutate(
        percentage = (pollinator_count / total_pollinator_count) * 100,
        strategy = strategy_name,
        n_plants = n_top
      ) %>%
      select(Study_Network_id, pollinator_count, total_pollinator_count, 
             percentage, strategy, n_plants)
  }
  
  phylo_10 <- calc_phylo_strategy(10, "Phylogenetic")
  phylo_5  <- calc_phylo_strategy(5, "Phylogenetic")
  phylo_3  <- calc_phylo_strategy(3, "Phylogenetic")
  
  # =========================================================
  # Option 4: Random selection (1000 replicates)
  # =========================================================
  
  set.seed(2025)
  
  plant_pool <- data_count_scaled %>%
    distinct(Study_Network_id, Plant_species, Flower_count_scaled, 
             .keep_all = TRUE)
  
  plant_pollinator <- data_interact_filt %>%
    select(Study_Network_id, Plant_original_name, Pollinator_accepted_name)
  
  get_pollinator_count <- function(n_sp){
    
    random_sp <- plant_pool %>%
      group_by(Study_Network_id) %>%
      group_modify(~ {
        df <- .x
        n_select <- min(n_sp, nrow(df))
        slice_sample(df, n = n_select)
      }) %>%
      ungroup()
    
    poll_sp <- plant_pollinator %>%
      inner_join(random_sp, 
                 by = c("Study_Network_id",
                        "Plant_original_name" = "Plant_species")) %>%
      group_by(Study_Network_id) %>%
      summarise(
        pollinator_count = n_distinct(Pollinator_accepted_name),
        .groups = "drop"
      )
    
    tibble(Study_Network_id = unique(plant_pool$Study_Network_id)) %>%
      left_join(poll_sp, by = "Study_Network_id") %>%
      mutate(pollinator_count = ifelse(is.na(pollinator_count), 0, 
                                       pollinator_count))
  }
  
  poll_sp_10 <- replicate(1000, get_pollinator_count(10), 
                          simplify = FALSE) %>% bind_rows() %>%
    group_by(Study_Network_id) %>%
    summarise(random_mean = mean(pollinator_count),
              random_sd = sd(pollinator_count), .groups = "drop")
  
  poll_sp_5 <- replicate(1000, get_pollinator_count(5), 
                         simplify = FALSE) %>% bind_rows() %>%
    group_by(Study_Network_id) %>%
    summarise(random_mean = mean(pollinator_count),
              random_sd = sd(pollinator_count), .groups = "drop")
  
  poll_sp_3 <- replicate(1000, get_pollinator_count(3), 
                         simplify = FALSE) %>% bind_rows() %>%
    group_by(Study_Network_id) %>%
    summarise(random_mean = mean(pollinator_count),
              random_sd = sd(pollinator_count), .groups = "drop")
  
  random_10 <- poll_sp_10 %>%
    left_join(total_number, by = "Study_Network_id") %>%
    mutate(
      percentage = (random_mean / total_pollinator_count) * 100,
      strategy = "Random",
      n_plants = 10
    ) %>%
    select(Study_Network_id, pollinator_count = random_mean, 
           total_pollinator_count, percentage, strategy, n_plants)
  
  random_5 <- poll_sp_5 %>%
    left_join(total_number, by = "Study_Network_id") %>%
    mutate(
      percentage = (random_mean / total_pollinator_count) * 100,
      strategy = "Random",
      n_plants = 5
    ) %>%
    select(Study_Network_id, pollinator_count = random_mean, 
           total_pollinator_count, percentage, strategy, n_plants)
  
  random_3 <- poll_sp_3 %>%
    left_join(total_number, by = "Study_Network_id") %>%
    mutate(
      percentage = (random_mean / total_pollinator_count) * 100,
      strategy = "Random",
      n_plants = 3
    ) %>%
    select(Study_Network_id, pollinator_count = random_mean, 
           total_pollinator_count, percentage, strategy, n_plants)
  
  # =========================================================
  # Combine all results
  # =========================================================
  
  result_all <- bind_rows(
    abun_10, abun_5, abun_3,
    shape_5, shape_3,
    phylo_10, phylo_5, phylo_3,
    random_10, random_5, random_3
  ) %>%
    mutate(pollinator_order = pollinator_order_name)
  
  return(result_all)
}


# ==========================================================
# Run analysis for all four pollinator orders
# ==========================================================

pollinator_orders <- c(
  "Hymenoptera",
  "Diptera",
  "Coleoptera",
  "Lepidoptera"
)

# Run for each order and combine
all_results <- map_df(
  pollinator_orders,
  analyze_pollinator_order
)


random_networks <- data_interact %>%
  filter(
    Interaction_addup > 0,
    !is.na(Pollinator_order)
  ) %>%
  distinct(
    pollinator_order = Pollinator_order,
    Study_Network_id
  )


result_all_corrected <- all_results %>%
  semi_join(
    networks_by_order,
    by = c(
      "pollinator_order",
      "Study_Network_id"
    )
  )

# ==========================================================
# Save results
# ==========================================================

write.csv(
  result_all_corrected,
  "data/processed/result_all_by_pollinator_order.csv",
  row.names = FALSE
)

saveRDS(
  result_all_corrected,
  "data/processed/result_all_by_pollinator_order.rds"
)

cat("\n✓ Analysis complete! Results saved.\n")
