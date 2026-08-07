############################################################
# Optimising plant species selection for automated monitoring
# Networks overview
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

# Number of plant species per network
plant_diversity <- data_count_scaled %>% 
  filter(Flower_count != 0) %>%
  group_by(Study_Network_id) %>%
  summarise(plant_sp_number = n_distinct(Plant_species, na.rm = TRUE))
# Number of pollinator species per network
pollinator_diversity <- data_interact %>%
  group_by(Study_Network_id) %>%
  summarise(pollinator_sp_number = n_distinct(Pollinator_accepted_name, na.rm = TRUE))
# Number of flower shapes per network
trait_diversity <- data_count_scaled %>%
  filter(Flower_count_scaled != 0, !is.na(Plant_species)) %>%
  group_by(Study_Network_id) %>%
  summarise(shape_number = n_distinct(flw_shape_revised, na.rm = TRUE))
# Check representation of families per flower shape
dominant_family <- data_merge %>%
  filter(Flower_count_scaled != 0, !is.na(Plant_species)) %>%
  filter(!is.na(Flower_count_scaled),!is.na(flw_shape_revised), !is.na(Plant_family)) %>%
  group_by(flw_shape_revised, Plant_family) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(flw_shape_revised) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()
# Check freq of shapes per network
shape_freq <- trait_diversity %>%
  group_by(shape_number) %>%
  summarise(freq = n())

# Check total number of shapes in the dataset
shape_class_freq <- data_count_scaled %>%
  filter(Flower_count_scaled != 0,
         !is.na(Plant_species),
         !is.na(flw_shape_revised)) %>%
  distinct(Study_Network_id, flw_shape_revised) %>%
  count(flw_shape_revised, name = "network_number") %>%
  arrange(desc(network_number))

# -----------------------------
#Explore pollinator frequencies
#Ensure same axis limits
x_max <- max(
  max(plant_diversity$plant_sp_number, na.rm = TRUE),
  max(pollinator_diversity$pollinator_sp_number, na.rm = TRUE))

breaks_all <- seq(0, x_max, by = 10)
# Plant distribution
histogram_1 <- ggplot(plant_diversity, aes(x = plant_sp_number, fill = "#009E73")) +
  geom_histogram(binwidth = 2, color = "black", boundary = 0) +
  labs(x = "Flower species", y = "Number of Networks") +
  scale_fill_identity() +
  scale_x_continuous(breaks = breaks_all, limits = c(0, x_max), expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 12, face = "bold"))
histogram_1
# -----------------------------
# Pollinator distribution
histogram_2 <- ggplot(pollinator_diversity, aes(x = pollinator_sp_number, fill = "#D55E00")) +
  geom_histogram(binwidth = 2, color = "black", boundary = 0) +
  labs(x = "Pollinator species in interaction data", y = "Number of Networks") +
  scale_fill_identity() +
  scale_x_continuous(breaks = breaks_all, limits = c(0, x_max), expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 12, face = "bold"))
histogram_2
# -----------------------------
# Flower shape distribution
shape_freq$shape_number <- as.factor(shape_freq$shape_number)
sum(as.numeric(as.character(shape_freq$shape_number)) * shape_freq$freq) /
  sum(shape_freq$freq)

shape_barplot <- ggplot(shape_freq, aes(x = shape_number, y = freq, fill = "#9E9AC8")) +
  geom_col(color = "black") +
  labs(
    x = "Number of flower shapes",
    y = "Number of Networks") +
  scale_fill_identity() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 12, face = "bold"))
shape_barplot
# -----------------------------
# 先从shape_class_freq中提取排序顺序 - 改为降序
shape_order <- shape_class_freq %>%
  arrange(network_number) %>%  # 改为降序
  pull(flw_shape_revised)

shape_class_freq <- shape_class_freq %>%
  filter(! flw_shape_revised %in% c("brush flowers","trap flowers"))#Trap flowers and brush flowers were excluded because each received fewer than four visits.

# ===== 重新绘制 shape_class =====

shape_class <- ggplot(
  shape_class_freq,
  aes(
    x = factor(flw_shape_revised, levels = shape_order),
    y = network_number,
    fill = "white")) +
  geom_col(color = "black") +
  labs(
    x = "Floral morphology",
    y = "Number of Networks") +
  scale_fill_identity() +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 13),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5))

shape_class
#ggsave("/Chap1_TargetPlant_to_monitor/result_260526/shape_class.png", shape_class, width = 4, 
#       height = 4, units = "in", dpi = 300)  

