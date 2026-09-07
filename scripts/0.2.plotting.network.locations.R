############################################################
# Optimising plant species selection for automated monitoring
# Part 2: Network overview and network locations
############################################################

# ==========================================================
# Load libraries
# ==========================================================

library(dplyr)
library(ggplot2) #For plotting
library(giscoR) #For plotting
library(sf) #For handling coordinates and plotting
library(patchwork)#For plotting (binding plots)
library(ggstar) #For plotting (cool shapes)
library(scales) #For plotting (decimals on axes)
library(tidyr)
library(viridis)
library(stringr)
library(purrr)
library(cowplot)

# ==========================================================
# Load data
# ==========================================================

data_count_scaled <- readRDS("data/processed/data_count_scaled_published.rds")
data_interact <- readRDS("data/processed/data_interact_published.rds")

n_distinct(data_interact$Study_id)

# ==========================================================
# Prepare data
# ==========================================================

data_merge <- merge(data_interact, data_count_scaled[,c("Flower_data_merger","Flower_count_scaled","Plant_species","Study_Network_id")], 
                   by = "Flower_data_merger",all = TRUE)%>%
  filter(!is.na(Flower_data_merger))%>%
  mutate(
    Study_Network_id = coalesce(Study_Network_id.x, Study_Network_id.y)) %>%
  select(-Study_Network_id.x, -Study_Network_id.y)%>%
  mutate(Interaction_addup = ifelse(is.na(Interaction_addup), 0, Interaction_addup))%>% # replace `Interaction_addup` == NA with 0
  mutate(
    Plant_accepted_name = str_squish(str_replace_all(replace_na(Plant_accepted_name, ""), "×", ""))) 

data <- data_interact%>% 
  mutate(Network_id = paste0(Study_id, Network_id))  

#Fix dates
data = data %>%
  dplyr::mutate(Year = lubridate::year(Date), 
                Month = lubridate::month(Date), 
                Day = lubridate::day(Date))

#Prepare coordinates
for_plotting = data %>% 
  select(Study_id, Latitude, Longitude) %>% 
  distinct()
#Get countries and transfrom coordinates
countries <- gisco_get_countries(
  resolution = 20) %>%
  st_transform(3035) %>% 
  sf::st_as_sf(crs = 4326)


# ==========================================================
# Summarise networks
# ==========================================================
#Number of plant and poll species per network
#Pollinators------
#Check unique number of pollinator species
PollNumber_by_ID = data %>% 
  select(Pollinator_rank, Pollinator_accepted_name, Study_id) %>% 
  filter(Pollinator_rank == "SPECIES") %>% 
  select(Pollinator_accepted_name, Study_id) %>% 
  group_by(Study_id) %>% 
  summarise(Pollinator_counts = n_distinct(Pollinator_accepted_name))
#Merge now with spp counts by network dataset
poll_data = left_join(PollNumber_by_ID, for_plotting, by = "Study_id")
#Create duplicate cord columns to operate on them
poll_data$coords <- data.frame(
  Longitude = poll_data$Longitude, Latitude = poll_data$Latitude)
#Set coordinates
poll_data$coords_sf = st_as_sf(poll_data$coords, 
                               coords = c("Longitude", "Latitude"), crs = 4326)
#Plot map of Europe
#Very noisy with many coordinates per point,
#Let's select a unique coordinate by study
poll_data1 =  poll_data %>% 
  distinct(Study_id, .keep_all = TRUE)

#-----------------------#
#Plants-----
#-----------------------#
#Check unique number of plant species
PlantNumber_by_ID = data %>% 
  select(Plant_rank, Plant_accepted_name, Study_id) %>% 
  filter(Plant_rank == "SPECIES") %>% 
  select(Plant_accepted_name, Study_id) %>% 
  group_by(Study_id) %>% 
  summarise(Plant_counts = n_distinct(Plant_accepted_name))
#Merge now with spp counts by network dataset
plant_data = left_join(PlantNumber_by_ID, for_plotting, by = "Study_id")
#Create duplicate cord columns to operate on them
plant_data$coords <- data.frame(
  Longitude = plant_data$Longitude, Latitude = plant_data$Latitude)
#Set coordinates
plant_data$coords_sf = st_as_sf(plant_data$coords, 
                                coords = c("Longitude", "Latitude"), crs = 4326)
#Plot map of Europe
#Very noisy with many coordinates per point,
#Let's select a unique coordinate by study
plant_data1 =  plant_data %>% 
  distinct(Study_id, .keep_all = TRUE)

#Add vars
plant_data1$Group = "Plant"
poll_data1$Group = "Pollinator"
#Rename cols
plant_data2 = plant_data1 %>% 
  rename(Counts = Plant_counts)
poll_data2 = poll_data1 %>% 
  rename(Counts = Pollinator_counts)
#Bind rows
all_data = bind_rows(plant_data2, poll_data2)

# ==========================================================
# Plot network map
# ==========================================================
#Try to plot both together (same plot)
map = ggplot(countries) +
  geom_sf(fill = "floralwhite", color = "black", alpha = 1, linewidth = 0.1) +
  geom_point(data = all_data, alpha = 1, aes(x = Longitude, y = Latitude,
                                             size = Counts, color = Group), shape = 1, stroke = 1) +
  coord_sf(crs = st_crs(3035), default_crs = st_crs(4326),
           xlim = c(-16, 39), ylim = c(37, 70),
           expand = FALSE) +
  scale_size_continuous(name = "Species counts", range = c(1, 12)) +
  scale_color_manual(name = "Taxonomic group",
                     values = c("Plant" = "#009E73", "Pollinator" = "#D55E00")) +
  guides(
    size = guide_legend(nrow = 1),
    color = guide_legend(override.aes = list(shape = 1,size = 6, stroke = 1))  
  ) +
  scale_x_continuous(breaks = seq(0, 20, by = 10)) +
  scale_y_continuous(breaks = seq(40, 60, by = 10)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "aliceblue"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.box.spacing = unit(0.1, "cm"),  
    legend.spacing.y = unit(0, "cm"),   
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, -5, -5, -5), 
    legend.title.align = 0.5,
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(vjust = -5, hjust = 0.057, size = 14))

map
#ggsave("/Chap1_TargetPlant_to_monitor/result_260526/map_1.png", map, width = 6, height = 6, units = "in", dpi = 600)

