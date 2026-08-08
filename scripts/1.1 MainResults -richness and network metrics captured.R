############################################################
# Optimising plant species selection for automated monitoring
# Main results
############################################################

# ==========================================================
# Load libraries
# ==========================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(giscoR)
library(patchwork)
library(ggstar)
library(scales)
library(viridis)
library(vegan)

# ==========================================================
# Load data
# ==========================================================

data_count_scaled <- readRDS("data/processed/data_count_scaled_published.rds")
data_interact <- readRDS("data/processed/data_interact_published.rds")
unique(data_count_scaled$Study_Network_id)
colnames(data_count_scaled)

unique(data_interact$Pollinator_rank)
# [1] "Hymenoptera"  "Diptera"      "Coleoptera"(26site)   "Lepidoptera"  
# less than 10 networks: 
# "Hemiptera"   "Thysanoptera"
# less than 6 networks:
# [6] "Neuroptera"   "Orthoptera"   "Mecoptera"    "Odonata"   "Dermaptera"   "Blattodea" 

#data_interact<-data_interact%>%filter(Pollinator_order == "Thysanoptera")

data_interact %>%
  group_by(Pollinator_order) %>%
  summarise(n_networks = n_distinct(Study_Network_id))

data_interact %>%
  distinct(Pollinator_order, Pollinator_accepted_name,Plant_accepted_name) %>%
  group_by(Pollinator_order) %>%
  summarise(unique_interactions = n())

length(unique(data_interact$Study_Network_id) )
length(unique(data_count_scaled$Study_Network_id) )

traits <- read.csv("data/processed/merge.trait.csv", header = TRUE, fileEncoding = "UTF-8")

# include plant with 0 visit into analysis
data_merge <- merge(data_interact, data_count_scaled[,c("Flower_data_merger","Flower_count_scaled","Plant_species","Study_Network_id")], 
                   by = "Flower_data_merger",all = TRUE)%>%
  filter(!is.na(Flower_data_merger))%>%
  mutate(
    Study_Network_id = coalesce(Study_Network_id.x, Study_Network_id.y)
  ) %>%
  dplyr::select(-Study_Network_id.x, -Study_Network_id.y)%>%
  mutate(Interaction_addup = ifelse(is.na(Interaction_addup), 0, Interaction_addup))%>% # 替换 `Interaction_addup` 为 NA 的值为 0
  mutate(
    Plant_accepted_name = str_squish(str_replace_all(replace_na(Plant_accepted_name, ""), "×", ""))
  ) 

# Checks
nrow(data_count_scaled)
nrow(data_interact)
nrow(data_merge)

data_count_scaled%>%filter(Flower_count_scaled == 0 )
unique(na.omit(data_merge$flw_shape_revised))

# flower shape check
shape_count <- data_count_scaled %>%
  group_by(Study_Network_id) %>%
  summarise(
    n_shape = n_distinct(flw_shape_revised),
    .groups = "drop"
  )

# Check how many networks contain 10, 5, and 3 flower shapes
shape_count %>%
  summarise(
    n_networks = n(),
    
    n_ge10 = sum(n_shape >= 10),
    prop_ge10 = mean(n_shape >= 10),
    
    n_ge5 = sum(n_shape >= 5),
    prop_ge5 = mean(n_shape >= 5),
    
    n_ge3 = sum(n_shape >= 3),
    prop_ge3 = mean(n_shape >= 3))

#====================================================
# Percentage of Pollinator Captured under Subsampling Strategies
#====================================================
### Plant identified to only genus level were excluded before we select the most abundant plant species.
### 在筛选前十时，必须先摆脱genus, 这一步在筛选的时候已经做了，将plant data都只保留属（interaction的没动）

# ----2.1 Option 1--------------------------------------
#Sort species within each site in descending order of abundance and select the top 10/5/3 species

colnames(data_count_scaled)

top_10_species <- data_count_scaled %>%
  distinct(Study_Network_id, Plant_species, Flower_count_scaled, .keep_all = FALSE) %>%
  group_by(Study_Network_id) %>%
  slice_max(order_by = Flower_count_scaled, n = 10, with_ties = FALSE) %>%  # Keep only the top 10 species per network (no ties allowed)
  ungroup()%>%
  group_by(Study_Network_id) %>%
  filter(n_distinct(Plant_species) > 9) %>%
  ungroup()%>%
  dplyr::select(Study_Network_id, Plant_species)

unique(top_10_species$Study_Network_id)


top_5_species <- data_count_scaled %>%
  distinct(Study_Network_id, Plant_species, Flower_count_scaled, .keep_all = TRUE) %>%
  group_by(Study_Network_id) %>%
  slice_max(order_by = Flower_count_scaled, n = 5, with_ties = FALSE) %>%  # Keep top 5 species per network (no ties allowed)
  ungroup()%>%
  group_by(Study_Network_id) %>%
  filter(n_distinct(Plant_species) > 4) %>%
  ungroup()%>%
  dplyr::select(Study_Network_id, Plant_species)

unique(top_5_species$Study_Network_id)
top_5_species%>%filter(is.na(Plant_species))

top_3_species <- data_count_scaled %>%
  distinct(Study_Network_id, Plant_species, Flower_count_scaled, .keep_all = FALSE) %>%
  group_by(Study_Network_id) %>%
  slice_max(order_by = Flower_count_scaled, n = 3, with_ties = FALSE) %>%  # Keep top 3 species per network (no ties allowed)
  ungroup()%>%
  group_by(Study_Network_id) %>%
  filter(n_distinct(Plant_species) > 2) %>%
  ungroup()%>%
  dplyr::select(Study_Network_id, Plant_species)

unique(top_3_species$Study_Network_id)
top_3_species%>%filter(is.na(Plant_species))
###### species number of the pollinators
result_orig<-  data_interact%>%
  filter(!Interaction_addup == 0)%>%
  group_by(Study_Network_id) %>%
  ungroup()

result10 <-  
  left_join(top_10_species, data_interact,
            by = c("Plant_species" = "Plant_original_name", 
                   "Study_Network_id" = "Study_Network_id")) %>%
  filter(!Interaction_addup == 0)%>%
  group_by(Study_Network_id) %>%
  ungroup()

#saveRDS(result10,"data/processed/selected_plant_result_abun10.rds")

sp_number10 <- result10 %>%
  ungroup%>%
  group_by(Study_Network_id) %>%
  filter(!Interaction_addup == 0)%>%
  summarise(pollinator_count = n_distinct(Pollinator_accepted_name,na.rm = TRUE))


total_number <- data_interact %>%
  ungroup%>%
  group_by(Study_Network_id) %>%
  filter(!Interaction_addup == 0)%>%
  summarise(total_pollinator_count = n_distinct(Pollinator_accepted_name,na.rm = TRUE))

# merge
merged_data10 <- merge(sp_number10,total_number,by.x = "Study_Network_id", by.y = "Study_Network_id", all = TRUE)%>%
  replace_na(list(pollinator_count = 0))


result5 <-  
  left_join(top_5_species, data_interact,
            by = c("Plant_species" = "Plant_original_name", 
                   "Study_Network_id" = "Study_Network_id")) %>%
  filter(!Interaction_addup == 0)%>%
  group_by(Study_Network_id) %>%
  ungroup()

sp_number5 <- result5 %>%
  group_by(Study_Network_id, Pollinator_accepted_name) %>%
  summarise(freq = sum(Interaction_addup, na.rm = TRUE), .groups = "drop") %>%
  group_by(Study_Network_id) %>%
  summarise(pollinator_count = n_distinct(Pollinator_accepted_name,na.rm = TRUE),
            shannon_index_abun5 = diversity(freq, index = "shannon"))
  

# merge
merged_data5 <- merge(sp_number5,total_number,by.x = "Study_Network_id", by.y = "Study_Network_id", all = TRUE)%>%
  replace_na(list(pollinator_count = 0))

result3 <-  
  left_join(top_3_species, data_interact,
            by = c("Plant_species" = "Plant_original_name", 
                   "Study_Network_id" = "Study_Network_id")) %>%
  filter(!Interaction_addup == 0)%>%
  group_by(Study_Network_id) %>%
  ungroup()

sp_number3 <- result3 %>%
  group_by(Study_Network_id, Pollinator_accepted_name) %>%
  summarise(freq = sum(Interaction_addup, na.rm = TRUE), .groups = "drop") %>%
  group_by(Study_Network_id) %>%
  summarise(pollinator_count = n_distinct(Pollinator_accepted_name,na.rm = TRUE),
            shannon_index_abun3 = diversity(freq, index = "shannon"))

# merge
merged_data3 <- merge(sp_number3,total_number,by.x = "Study_Network_id", by.y = "Study_Network_id", all = TRUE)%>%
  replace_na(list(pollinator_count = 0))

####### proportion of the pollinator species captured
percent_10 <- merged_data10 %>%
  mutate(percentage = (pollinator_count / total_pollinator_count) * 100)

percent_5 <- merged_data5 %>%
  mutate(percentage = (pollinator_count / total_pollinator_count) * 100)

percent_3 <- merged_data3 %>%
  mutate(percentage = (pollinator_count / total_pollinator_count) * 100)

##################################################################

# saveRDS(percent_10,"percent_10.rds")

########################################################################
#proportion of the unique interaction captured
unic_inter <- function(subset_data, full_data) {
  
  subset_n <- subset_data %>%
    ungroup() %>%
    distinct(
      Study_Network_id,
      Plant_accepted_name,
      Pollinator_accepted_name
    ) %>%
    count(Study_Network_id, name = "n_subset")
  
  full_n <- full_data %>%
    ungroup() %>%
    distinct(
      Study_Network_id,
      Plant_accepted_name,
      Pollinator_accepted_name
    ) %>%
    count(Study_Network_id, name = "n_total")
  
  subset_n %>%
    left_join(full_n, by = "Study_Network_id") %>%
    mutate(prop_interactions = n_subset / n_total)
}

cov_10 <- unic_inter(result10, data_interact)
cov_5  <- unic_inter(result5, data_interact)
cov_3  <- unic_inter(result3, data_interact)
# 
# colnames(cov_10)[colnames(cov_10) == "prop_interactions"] <- "Abundant_Top10"
# colnames(cov_5)[colnames(cov_5) == "prop_interactions"] <- "Abundant_Top5"
# colnames(cov_3)[colnames(cov_3) == "prop_interactions"] <- "Abundant_Top3"


################################################################################

###---2.1 Option 2--------------------------------------------


##only considered the interactions of the most abundant plant species for each flower shape type
####choose top 5/3 flower shape to make it comparable with op1
# Method: For each flower shape, select the species with the highest abundance, but only select 3 species in total.
# That is, if the network has more than 3 flower shape types, first select the most abundant species from each flower shape,
# then rank these species by abundance and choose the top 3 species.
# If a network has fewer than 3 flower shape types, it will be excluded.

##notice: 
# By default, the rank() function uses ties.method = "average",
# which means species with the same Flower_count_scaled will receive the same (average) rank.
# As a result, more than 3 species may be selected if multiple species are tied for the 3rd highest abundance.
# To ensure only 3 species are selected per network, even in the case of ties,
# consider using ties.method = "first" in the rank() function,
# or alternatively, use slice_head(n = 3) after arranging by abundance to strictly limit the output to the top 3.
data_merge_trait <- left_join(data_merge, traits, by = "Plant_accepted_name")

data_count_scaled_trait<-data_count_scaled%>%
  left_join(traits %>% rename(Plant_species = Plant_accepted_name), by = "Plant_species")

colnames(data_count_scaled)
nrow(data_count_scaled_trait)
nrow(data_count_scaled)

## overview of the flower shape

flw.no <- data_count_scaled%>%
  group_by(Study_Network_id) %>%
  summarise(n_flower_shapes = n_distinct(flw_shape_revised))

flw_networks <- n_distinct(flw.no$Study_Network_id)

pct_3 <- mean(flw.no$n_flower_shapes >= 3) * 100
pct_5 <- mean(flw.no$n_flower_shapes >= 5) * 100
pct_10 <- mean(flw.no$n_flower_shapes >= 10) * 100

pct_3
pct_5
pct_10

#Identify the most abundant species in each floral shape
abun_species_top5 <- data_count_scaled %>%
  group_by(Study_Network_id, flw_shape_revised) %>%
  arrange(Study_Network_id, flw_shape_revised, desc(Flower_count_scaled)) %>%
  top_n(1, wt = Flower_count_scaled) %>%
  distinct(Study_Network_id, Plant_species, Flower_count_scaled) %>%  # .keep_all = TRUE (keeps all columns)
  group_by(Study_Network_id) %>%
  mutate(rank = rank(desc(Flower_count_scaled),ties.method = "first")) %>% 
  filter(rank < 6) %>%
  ungroup() %>%
  group_by(Study_Network_id) %>%
  filter(n_distinct(Plant_species) > 4) %>%
  ungroup()%>%
  dplyr::select(Plant_species,Study_Network_id)


n_distinct(abun_species_top5$Study_Network_id)

abun_species_top3 <- data_count_scaled %>%
  group_by(Study_Network_id, flw_shape_revised) %>%
  arrange(Study_Network_id, flw_shape_revised, desc(Flower_count_scaled)) %>%
  top_n(1, wt = Flower_count_scaled) %>%
  distinct(Study_Network_id, Plant_species, Flower_count_scaled) %>%  # .keep_all = TRUE 保留所有列
  group_by(Study_Network_id) %>%
  mutate(rank = rank(desc(Flower_count_scaled),ties.method = "first")) %>% #只选择一个而非多个
  #mutate(rank = rank(desc(Flower_count_scaled,ties.method = "first"))) %>% # Here, if we set ties.method = "first" to ensure only 3 flower shapes are selected,
  # so even if two flower shapes are tied for 3rd place, only one will be chosen.
  filter(rank < 4) %>%
  ungroup() %>%
  group_by(Study_Network_id) %>%
  filter(n_distinct(Plant_species) > 2) %>%
  ungroup()%>%
  dplyr::select(Plant_species,Study_Network_id)

length(unique(abun_species_top3$Study_Network_id))
length(unique(abun_species_top5$Study_Network_id))

######

#pollinators that interact with these plant species at each site
abun_species_top3_merge <-  left_join(abun_species_top3, data_interact,
                            by = c("Plant_species" = "Plant_original_name", 
                                   "Study_Network_id" = "Study_Network_id")) %>%
  group_by(Study_Network_id) %>%
  ungroup()

op2_sp_number_top3<-abun_species_top3_merge%>%
  filter(Interaction_addup >0) %>%
  group_by(Study_Network_id, Pollinator_accepted_name) %>%
  summarise(freq = sum(Interaction_addup, na.rm = TRUE), .groups = "drop") %>%
  group_by(Study_Network_id) %>%
  summarise(pollinator_count = n_distinct(Pollinator_accepted_name,na.rm = TRUE),
            shannon_index_flw3 = diversity(freq, index = "shannon"))

abun_species_top5_merge <- left_join(abun_species_top5, data_interact,
                                     by = c("Plant_species" = "Plant_original_name", 
                                            "Study_Network_id" = "Study_Network_id")) %>%
  group_by(Study_Network_id) %>%
  ungroup()

op2_sp_number_top5<-abun_species_top5_merge%>%
  filter(Interaction_addup >0) %>%
  group_by(Study_Network_id, Pollinator_accepted_name) %>%
  summarise(freq = sum(Interaction_addup, na.rm = TRUE), .groups = "drop") %>%
  group_by(Study_Network_id) %>%
  summarise(pollinator_count = n_distinct(Pollinator_accepted_name,na.rm = TRUE),
            shannon_index_flw5 = diversity(freq, index = "shannon"))



#计算占比
#Calculate the percentage
op2_merged_data_top3 <- merge(op2_sp_number_top3,total_number,by.x = "Study_Network_id", by.y = "Study_Network_id", all.x = TRUE)
op2_percent_shape_top3 <- op2_merged_data_top3 %>%
  mutate(percentage = (pollinator_count / total_pollinator_count) * 100)
head(op2_percent_shape_top3)

op2_merged_data_top5 <- merge(op2_sp_number_top5,total_number,by.x = "Study_Network_id", by.y = "Study_Network_id", all.x = TRUE)
op2_percent_shape_top5 <- op2_merged_data_top5 %>%
  mutate(percentage = (pollinator_count / total_pollinator_count) * 100)
head(op2_percent_shape_top5)




#######################

#------unique interactions

######################
head(abun_species_top3_merge)
cov_floral_5  <- unic_inter(abun_species_top5_merge, data_interact)
cov_floral_3  <- unic_inter(abun_species_top3_merge, data_interact)

colnames(cov_floral_5)[colnames(cov_floral_5) == "prop_interactions"] <- "FlwShape_Top5"
colnames(cov_floral_3)[colnames(cov_floral_3) == "prop_interactions"] <- "FlwShape_Top3"



############################################################################
#----2.3 Phylo distance-----------------------

library(dplyr)
library(stringr)
library(tidyr)
library(V.PhyloMaker2)
library(ape)
#library("devtools")
#devtools::install_github("jinyizju/V.PhyloMaker2")

#------------------------------------
# 1. Prepare species list for phylogeny
#------------------------------------

plant_sp <- data_count_scaled %>%
  distinct(Plant_species) %>%
  na.omit() %>%
  filter(
    str_count(Plant_species, "\\S+") >= 2,
    !grepl("\\bsp\\.?\\b|cf\\.|aff\\.",
           Plant_species,
           ignore.case = TRUE)
  ) %>%
  rename(species = Plant_species) %>%
  mutate(
    species = gsub(" ", "_", species)
  ) %>%
  left_join(
    tips.info.TPL[, c("species", "genus", "family")],
    by = "species"
  )


#------------------------------------
# 2. Build phylogenetic tree
#------------------------------------

data("GBOTB.extended.TPL")
#data("nodes.info.TPL")


phylo_result <- phylo.maker(
  sp.list = plant_sp,
  tree = GBOTB.extended.TPL,
  scenarios = "S3"
)

phylo_tree <- phylo_result$scenario.3


# optional check
phylo_tree


#------------------------------------
# 3. Prepare PD dataset
#------------------------------------

data_PD <- data_count_scaled %>%
  filter(
    !is.na(Plant_species),
    str_count(Plant_species, "\\S+") >= 2,
    !grepl("\\bsp\\.?\\b|cf\\.|aff\\.",
           Plant_species,
           ignore.case = TRUE)
  ) %>%
  mutate(
    Plant_species_tree = gsub(" ", "_", Plant_species)
  )


# species that are present in tree

data_PD_tree <- data_PD %>%
  filter(
    Plant_species_tree %in% phylo_tree$tip.label
  )


#------------------------------------
# 4. Check phylogenetic matching
#------------------------------------

match_check <- data_PD %>%
  group_by(Study_Network_id) %>%
  summarise(
    n_total = n_distinct(Plant_species_tree),
    n_tree = sum(
      Plant_species_tree %in% phylo_tree$tip.label
    ),
    match_rate = n_tree / n_total,
    .groups="drop"
  )

summary(match_check$n_tree)
head(match_check)



#------------------------------------
# 5. Select phylogenetic diversity species
#------------------------------------

select_PD_species <- function(species_pool, phylo_tree, n_select){
  
  species_pool <- intersect(
    species_pool,
    phylo_tree$tip.label
  )
  
  if(length(species_pool) < n_select){
    return(NA)
  }
  
  selected <- sample(species_pool,1)
  
  while(length(selected) < n_select){
    
    candidates <- setdiff(
      species_pool,
      selected
    )
    
    pd_gain <- sapply(
      candidates,
      function(x){
        
        tips <- c(selected,x)
        
        sub_tree <- drop.tip(
          phylo_tree,
          setdiff(
            phylo_tree$tip.label,
            tips
          )
        )
        
        sum(sub_tree$edge.length)
      }
    )
    
    selected <- c(
      selected,
      candidates[which.max(pd_gain)]
    )
  }
  
  selected
}



#------------------------------------
# 6. PD top10, top5, top3
#------------------------------------

PD_selected <- lapply(
  c(10,5,3),
  function(n){
    
    data_PD_tree %>%
      group_by(Study_Network_id) %>%
      summarise(
        Plant_species_tree = list(
          select_PD_species(
            Plant_species_tree,
            phylo_tree,
            n
          )
        ),
        .groups="drop"
      ) %>%
      filter(!is.na(Plant_species_tree)) %>%
      unnest(Plant_species_tree) %>%
      mutate(
        Plant_species = gsub("_"," ",Plant_species_tree),
        PD_size = n
      )
  }
)


names(PD_selected) <- c(
  "PD_top10",
  "PD_top5",
  "PD_top3"
)



#------------------------------------
# 7. Calculate pollinator capture
#------------------------------------

calc_capture <- function(selected_data){
  
  left_join(
    selected_data,
    data_interact,
    by=c(
      "Plant_species"="Plant_original_name",
      "Study_Network_id"
    )
  ) %>%
    group_by(Study_Network_id) %>%
    summarise(
      pollinator_count =
        n_distinct(
          Pollinator_accepted_name[
            Interaction_addup > 0 &
              !is.na(Pollinator_accepted_name)
          ]
        ),
      .groups="drop"
    ) %>%
    mutate(
      pollinator_count =
        replace_na(pollinator_count,0)
    )
}



PD_sp_number10 <- calc_capture(PD_selected$PD_top10)
PD_sp_number5  <- calc_capture(PD_selected$PD_top5)
PD_sp_number3  <- calc_capture(PD_selected$PD_top3)



#------------------------------------
# 8. Calculate percentage
#------------------------------------

PD_percent10 <- PD_sp_number10 %>%
  left_join(total_number,
            by="Study_Network_id") %>%
  mutate(
    percentage =
      pollinator_count /
      total_pollinator_count * 100
  )


PD_percent5 <- PD_sp_number5 %>%
  left_join(total_number,
            by="Study_Network_id") %>%
  mutate(
    percentage =
      pollinator_count /
      total_pollinator_count * 100
  )


PD_percent3 <- PD_sp_number3 %>%
  left_join(total_number,
            by="Study_Network_id") %>%
  mutate(
    percentage =
      pollinator_count /
      total_pollinator_count * 100
  )

#------------------------------------
# 9. Create full interaction datasets for network metrics
#------------------------------------

result_PD10 <- left_join(
  PD_selected$PD_top10,
  data_interact,
  by=c(
    "Plant_species"="Plant_original_name",
    "Study_Network_id"
  )
)

result_PD5 <- left_join(
  PD_selected$PD_top5,
  data_interact,
  by=c(
    "Plant_species"="Plant_original_name",
    "Study_Network_id"
  )
)

result_PD3 <- left_join(
  PD_selected$PD_top3,
  data_interact,
  by=c(
    "Plant_species"="Plant_original_name",
    "Study_Network_id"
  )
)

######interaction coverage
PD_cov10 <- unic_inter(
  left_join(
    PD_selected$PD_top10,
    data_interact,
    by=c(
      "Plant_species"="Plant_original_name",
      "Study_Network_id"="Study_Network_id"
    )
  ),
  data_interact
)


PD_cov5 <- unic_inter(
  left_join(
    PD_selected$PD_top5,
    data_interact,
    by=c(
      "Plant_species"="Plant_original_name",
      "Study_Network_id"="Study_Network_id"
    )
  ),
  data_interact
)


PD_cov3 <- unic_inter(
  left_join(
    PD_selected$PD_top3,
    data_interact,
    by=c(
      "Plant_species"="Plant_original_name",
      "Study_Network_id"="Study_Network_id"
    )
  ),
  data_interact
)


colnames(PD_cov10)[4] <- "Phylo_Top10"
colnames(PD_cov5)[4]  <- "Phylo_Top5"
colnames(PD_cov3)[4]  <- "Phylo_Top3"


###################################################################

#----2.4 Random selection----------------------

library(dplyr)

set.seed(2025)

# Step 1 预处理数据（提高速度）
plant_pool <- data_count_scaled %>%
  distinct(Study_Network_id, Plant_species, Flower_count_scaled, .keep_all = TRUE)

plant_pollinator <- data_interact %>%
  dplyr::select(Study_Network_id,
         Plant_original_name,
         Pollinator_accepted_name)


# Step 2 定义函数：计算pollinator richness

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
  
  #  强制补全 309 个 network
  poll_sp_complete <- tibble(Study_Network_id = unique(plant_pool$Study_Network_id)) %>%
    left_join(poll_sp, by = "Study_Network_id") %>%
    mutate(pollinator_count = ifelse(is.na(pollinator_count), 0, pollinator_count))
  
  return(poll_sp_complete)
}



# Step 3 random sample for 1000 replicates

poll_sp_10_all <- replicate(1000, get_pollinator_count(10), simplify = FALSE) %>%
  bind_rows()

poll_sp_5_all <- replicate(1000, get_pollinator_count(5), simplify = FALSE) %>%
  bind_rows()

poll_sp_3_all <- replicate(1000, get_pollinator_count(3), simplify = FALSE) %>%
  bind_rows()

# Step 4 Average value

poll_sp_10 <- poll_sp_10_all %>%
  group_by(Study_Network_id) %>%
  summarise(
    random_mean = mean(pollinator_count),
    random_sd = sd(pollinator_count),
    .groups = "drop"
  )


poll_sp_5 <- poll_sp_5_all %>%
  group_by(Study_Network_id) %>%
  summarise(
    random_mean = mean(pollinator_count),
    random_sd = sd(pollinator_count),
    .groups = "drop"
  )


poll_sp_3 <- poll_sp_3_all %>%
  group_by(Study_Network_id) %>%
  summarise(
    random_mean = mean(pollinator_count),
    random_sd = sd(pollinator_count),
    .groups = "drop"
  )

head(poll_sp_10)

merged_random10 <- merge(poll_sp_10,total_number,by.x = "Study_Network_id", by.y = "Study_Network_id", all = TRUE)
merged_random5 <- merge(poll_sp_5,total_number,by.x = "Study_Network_id", by.y = "Study_Network_id", all = TRUE)
merged_random3 <- merge(poll_sp_3,total_number,by.x = "Study_Network_id", by.y = "Study_Network_id", all = TRUE)

####### proportion of the pollinator species captured
random_10 <- merged_random10 %>%
  mutate(
    random_mean_percentage =
      random_mean / total_pollinator_count * 100,
    
    random_sd_percentage =
      random_sd / total_pollinator_count * 100
  )

random_5 <- merged_random5 %>%
  mutate(
    random_mean_percentage =
      random_mean / total_pollinator_count * 100,
    
    random_sd_percentage =
      random_sd / total_pollinator_count * 100
  )

random_3 <- merged_random3 %>%
  mutate(
    random_mean_percentage =
      random_mean / total_pollinator_count * 100,
    
    random_sd_percentage =
      random_sd / total_pollinator_count * 100
  )
unique(random_3$Study_Network_id)

## unique interaction coverage
# main function (ONLY interactions)
#========================
get_interaction_metrics <- function(n_sp){
  
  # sample plants
  random_sp <- plant_pool %>%
    group_by(Study_Network_id) %>%
    group_modify(~ {
      
      df <- .x
      n_select <- min(n_sp, nrow(df))
      
      slice_sample(df, n = n_select)
    }) %>%
    ungroup()
  
  # subset interactions
  sampled_data <- plant_pollinator %>%
    inner_join(
      random_sp,
      by = c("Study_Network_id",
             "Plant_original_name" = "Plant_species")
    ) %>%
    distinct(
      Study_Network_id,
      Plant_original_name,
      Pollinator_accepted_name
    )
  
  # interaction richness per network
  subset_metrics <- sampled_data %>%
    group_by(Study_Network_id) %>%
    summarise(
      interaction_richness = n(),
      .groups = "drop"
    )
  
  # full network baseline
  full_metrics <- plant_pollinator %>%
    distinct(
      Study_Network_id,
      Plant_original_name,
      Pollinator_accepted_name
    ) %>%
    group_by(Study_Network_id) %>%
    summarise(
      total_interactions = n(),
      .groups = "drop"
    )
  
  # merge + coverage
  tibble(Study_Network_id = unique(plant_pool$Study_Network_id)) %>%
    left_join(subset_metrics, by = "Study_Network_id") %>%
    left_join(full_metrics, by = "Study_Network_id") %>%
    mutate(
      interaction_richness = ifelse(is.na(interaction_richness), 0, interaction_richness),
      interaction_coverage = interaction_richness / total_interactions
    )
}

#========================
#replicate function
#========================
run_rep <- function(n_sp, n_rep = 1000){
  
  replicate(n_rep, get_interaction_metrics(n_sp), simplify = FALSE) %>%
    bind_rows()
}

#========================
# Step 4: run scenarios
#========================

cov_random_10 <- run_rep(10)
cov_random_5  <- run_rep(5)
cov_random_3  <- run_rep(3)
colnames(cov_random_10)[colnames(cov_random_10) == "interaction_coverage"] <- "Random_10"
colnames(cov_random_5)[colnames(cov_random_5) == "interaction_coverage"] <- "Random_5"
colnames(cov_random_3)[colnames(cov_random_3) == "interaction_coverage"] <- "Random_3"

saveRDS(cov_random_10,"cov_random_10.rds")
saveRDS(cov_random_5,"cov_random_5.rds")
saveRDS(cov_random_3,"cov_random_3.rds")

cov_random_10<-readRDS("cov_random_10.rds")
cov_random_5<-readRDS("cov_random_5.rds")
cov_random_3<-readRDS("cov_random_3.rds")
#######
library(dplyr)
library(purrr)

# merge the results
percentage_datasets <- list(
  percent_5,
  percent_3,
  percent_10,
  op2_percent_shape_top5,
  op2_percent_shape_top3,
  PD_percent10,
  PD_percent5,
  PD_percent3,
  random_10,
  random_5,
  random_3
)
head(random_10)
names(percentage_datasets) <- c(
  "Abun5",
  "Abun3",
  "Abun10",
  "FlwShape5",
  "FlwShape3",
  "Pylo10",
  "Pylo5",
  "Pylo3",
  "Random10",
  "Random5",
  "Random3"
)
percentage_datasets <- lapply(percentage_datasets, function(df) {
  df %>%
    mutate(Study_Network_id = as.character(Study_Network_id))
})
percentage_datasets <- imap(percentage_datasets, function(df, nm) {
  
  df %>%
    rename_with(
      function(cols) paste0(cols, "_", nm),
      -Study_Network_id
    )
  
})

all_networks <- data.frame(
  Study_Network_id = as.character(unique(plant_pool$Study_Network_id))
)

result_all <- reduce(percentage_datasets, full_join, by = "Study_Network_id") %>%
  right_join(all_networks, by = "Study_Network_id")
head(result_all)
n_distinct(result_all$Study_Network_id)

subsampled_networks<-list(
  orig = result_orig,
  Abun10 = result10,
  Abun5 = result5,
  Abun3 = result3,
  FlwShape5 = abun_species_top5_merge,
  FlwShape3 = abun_species_top3_merge,
  Pylo10 = result_PD10,
  Pylo5 = result_PD5,
  Pylo3 = result_PD3
)
  
  
write.csv(result_all,"data/processed/result_all_published_PD.csv", row.names = TRUE)   
saveRDS(percentage_datasets,"data/processed/MainResults_richness.rds")
saveRDS(subsampled_networks,"data/processed/Subsampled_networks.rds")


############ unique coverage result merge

library(dplyr)
library(purrr)

unic_inter_datasets <- list(
  cov_10,
  cov_5,
  cov_3,
  cov_floral_5,
  cov_floral_3,
  cov_random_10,
  cov_random_5,
  cov_random_3
)

# 修正：用 unic_inter_datasets，不是 datasets
unic_inter_datasets <- lapply(unic_inter_datasets, function(df) {
  df %>%
    mutate(Study_Network_id = as.character(Study_Network_id)) %>%
    distinct(Study_Network_id, .keep_all = TRUE)
})

unic_inter_all_networks <- data.frame(
  Study_Network_id = as.character(unique(plant_pool$Study_Network_id))
)

# ✔关键修复：用 full_join，不用 bind_rows
unic_inter_result_all <- reduce(
  unic_inter_datasets,
  full_join,
  by = "Study_Network_id"
) %>%
  right_join(unic_inter_all_networks, by = "Study_Network_id")

n_distinct(unic_inter_result_all$Study_Network_id)

write.csv(
  unic_inter_result_all,
  "unic_inter_result_all_published.csv",
  row.names = FALSE
)

unique(unic_inter_result_all$Study_Network_id)


#############################
# ------------------------------------------------------------
# Z-score standardization relative to random benchmark
# ------------------------------------------------------------
# Formula:
# Z = (X_strategy - mean_random) / sd_random
#
# where:
# X_strategy:
#   The percentage of pollinator species captured by the selected
#   plant species using an informed strategy (e.g., flower abundance,
#   flower traits, or phylogenetic distance).
#
# mean_random:
#   The mean percentage of pollinator species captured across 1000
#   random plant selections within the same network.
#
# sd_random:
#   The standard deviation of the percentage captured across the
#   1000 random selections within the same network.


calc_zscore <- function(observed,
                        random,
                        strategy,
                        number){
  
  observed %>%
    left_join(
      random %>%
        select(
          Study_Network_id,
          random_mean_percentage,
          random_sd_percentage
        ),
      by="Study_Network_id"
    ) %>%
    mutate(
      percentage =
        pollinator_count /
        total_pollinator_count *100,
      
      Z_score =
        (percentage-random_mean_percentage)/
        random_sd_percentage,
      
      Strategy=strategy,
      Plant_number=number
    )
}

Abun_z10 <- calc_zscore(percent_10,random_10,"Flower abundance","Top10")
Abun_z5 <- calc_zscore(percent_5,random_5,"Flower abundance","Top5")
Abun_z3 <- calc_zscore(percent_3,random_3,"Flower abundance","Top3")
Shape_z5 <- calc_zscore(op2_percent_shape_top5,random_5,"Flower abundance + shapes","Top5")
Shape_z3 <- calc_zscore(op2_percent_shape_top3,random_3,"Flower abundance + shapes","Top3")
PD_z10 <- calc_zscore(PD_percent10,random_10,"Phylogenetic distance","Top10")
PD_z5 <- calc_zscore(PD_percent5,random_5,"Phylogenetic distance","Top5")
PD_z3 <- calc_zscore(PD_percent3,random_3,"Phylogenetic distance","Top3")

z_results <- bind_rows(
  Abun_z10,
  Abun_z5,
  Abun_z3,
  Shape_z5,
  Shape_z3,
  PD_z10,
  PD_z5,
  PD_z3
)

z_results <- z_results %>%
  mutate(
    Z_score = ifelse(
      random_sd_percentage == 0,
      NA,
      Z_score
    )
  )%>%
  select(
    Study_Network_id,
    Strategy,
    Plant_number,
    percentage,
    random_mean_percentage,
    random_sd_percentage,
    Z_score
  )%>%
  mutate(
    no_random_variation = random_sd_percentage == 0,
    Plant_number = factor(
      Plant_number,
      levels = c("Top10", "Top5", "Top3")
    )
  )

head(z_results)
saveRDS(z_results,"data/processed/z_results.rds")

