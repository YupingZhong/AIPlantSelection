#Networks Overview Map
##1.Data overview
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

##################plant and pollinator species distribution
# read data
data_count_scaled<-readRDS("data/processed/data_count_scaled_published_0526.rds")#全部互作数据
data_interact<-readRDS("data/processed/data_interact_published_0526.rds")

n_distinct(data_interact$Study_id)

#0visit的,只出现在植物调查之中的植物也应该参与此分析
data_merge<- merge(data_interact, data_count_scaled[,c("Flower_data_merger","Flower_count_scaled","Plant_species","Study_Network_id")], 
                   by = "Flower_data_merger",all = TRUE)%>%
  filter(!is.na(Flower_data_merger))%>%
  mutate(
    Study_Network_id = coalesce(Study_Network_id.x, Study_Network_id.y)
  ) %>%
  select(-Study_Network_id.x, -Study_Network_id.y)%>%
  mutate(Interaction_addup = ifelse(is.na(Interaction_addup), 0, Interaction_addup))%>% # 替换 `Interaction_addup` 为 NA 的值为 0
  mutate(
    Plant_accepted_name = str_squish(str_replace_all(replace_na(Plant_accepted_name, ""), "×", ""))
  ) 


data<-data_interact%>% 
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

#Number of plant and poll species per network
#-----------------------#
#Pollinators------
#-----------------------#
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
    color = guide_legend(override.aes = list(shape = 1,size = 6, stroke = 1))  # 图例中颜色点变成空心圆圈
  ) +
  scale_x_continuous(breaks = seq(0, 20, by = 10)) +
  scale_y_continuous(breaks = seq(40, 60, by = 10)) +
  xlab("Longitude") +
  ylab("Latitude") +
  #ggtitle("(a)") +
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
    legend.box.spacing = unit(0.1, "cm"),  # 缩小行间距
    legend.spacing.y = unit(0, "cm"),     # 减小图例内部上下间距
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, -5, -5, -5), # 调整图例盒子边距防止重叠
    legend.title.align = 0.5,
    axis.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(vjust = -5, hjust = 0.057, size = 14)
  )

map
#ggsave("/Chap1_TargetPlant_to_monitor/result_260526/map_1.png", map, width = 6, height = 6, units = "in", dpi = 300)


##################################

plant_diversity <- data_count_scaled %>% 
  filter(Flower_count != 0) %>%
  group_by(Study_Network_id) %>%
  summarise(plant_sp_number = n_distinct(Plant_species, na.rm = TRUE))

pollinator_diversity <- data_interact %>%
  group_by(Study_Network_id) %>%
  summarise(pollinator_sp_number = n_distinct(Pollinator_accepted_name, na.rm = TRUE))

trait_diversity<-data_count_scaled%>%
  filter(Flower_count_scaled != 0, !is.na(Plant_species)) %>%
  group_by(Study_Network_id) %>%
  summarise(shape_number = n_distinct(flw_shape_revised, na.rm = TRUE))


dominant_family <- data_merge %>%
  filter(Flower_count_scaled != 0, !is.na(Plant_species)) %>%
  filter(!is.na(Flower_count_scaled),!is.na(flw_shape_revised), !is.na(Plant_family)) %>%
  group_by(flw_shape_revised, Plant_family) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(flw_shape_revised) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

shape_freq <- trait_diversity %>%
  group_by(shape_number) %>%
  summarise(freq = n())

shape_class_freq <- data_count_scaled %>%
  filter(Flower_count_scaled != 0,
         !is.na(Plant_species),
         !is.na(flw_shape_revised)) %>%
  distinct(Study_Network_id, flw_shape_revised) %>%
  count(flw_shape_revised, name = "network_number") %>%
  arrange(desc(network_number))

x_max <- max(
  max(plant_diversity$plant_sp_number, na.rm = TRUE),
  max(pollinator_diversity$pollinator_sp_number, na.rm = TRUE)
)

breaks_all <- seq(0, x_max, by = 10)

# -----------------------------
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
        strip.text = element_text(size = 12, face = "bold")
        )



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
        strip.text = element_text(size = 12, face = "bold")
  )

# -----------------------------
# Flower shape distribution
shape_freq$shape_number <- as.factor(shape_freq$shape_number)
sum(as.numeric(as.character(shape_freq$shape_number)) * shape_freq$freq) /
  sum(shape_freq$freq)

shape_barplot <- ggplot(shape_freq, aes(x = shape_number, y = freq, fill = "#9E9AC8")) +
  geom_col(color = "black") +
  labs(
    x = "Number of flower shapes",
    y = "Number of Networks"
  ) +
  scale_fill_identity() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(size = 12, face = "bold")
  )

# -----------------------------
# 先从shape_class_freq中提取排序顺序 - 改为降序
shape_order_final <- shape_class_freq %>%
  arrange(network_number) %>%  # 改为降序
  pull(flw_shape_revised)

shape_class_freq <- shape_class_freq %>%
  filter(! flw_shape_revised %in% c("brush flowers","trap flowers"))#Trap flowers and brush flowers were excluded because each received fewer than four visits.

# ===== 重新绘制 shape_class =====
shape_order

shape_class <- ggplot(
  shape_class_freq,
  aes(
    x = factor(flw_shape_revised, levels = shape_order),
    y = network_number,
    fill = "white"
  )
) +
  geom_col(color = "black") +
  labs(
    x = "Floral morphology",
    y = "Number of Networks"
  ) +
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
    plot.margin = margin(5, 15, 5, 5)
  )


shape_class

ggsave("/Chap1_TargetPlant_to_monitor/result_260526/shape_class.png", shape_class, width = 4, 
       height = 4, units = "in", dpi = 300)  


#--------------------------
# abundant_base pollinator richness accumulation curve VS random
data_merge %>%
  filter(Interaction_addup == 0)

a<-data_count_scaled%>%group_by(Study_Network_id)%>%n_distinct(Plant_species)

plant_rank <- data_count_scaled %>%
  group_by(Study_Network_id, Plant_species) %>%
  summarise(
    Abundance = first(Flower_count_scaled),
    .groups = "drop"
  ) %>%
  group_by(Study_Network_id) %>%
  arrange(desc(Abundance), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

poll_pool <- data_interact %>%
  filter(!is.na(Pollinator_accepted_name)) %>%
  group_by(Study_Network_id) %>%
  summarise(
    total_poll = n_distinct(Pollinator_accepted_name),
    .groups = "drop"
  )

library(dplyr)
library(purrr)
valid_nets <- plant_rank %>%
  group_by(Study_Network_id) %>%
  summarise(
    n_plants_total = n_distinct(Plant_species),
    .groups = "drop"
  ) #%>%
  #filter(n_plants_total >= 20)

plant_rank_filt <- plant_rank %>%
  filter(Study_Network_id %in% valid_nets$Study_Network_id)


poll_map <- data_interact %>%
  group_by(Study_Network_id, Plant_species = Plant_original_name) %>%
  summarise(
    poll = list(unique(Pollinator_accepted_name)),
    .groups = "drop"
  )

acc_curve <- plant_rank_filt %>%
  group_by(Study_Network_id) %>%
  group_split() %>%
  map_dfr(function(net){
    
    net_id <- unique(net$Study_Network_id)
    
    plants <- net %>%
      arrange(desc(Abundance)) %>%
      pull(Plant_species)
    
    total_poll <- poll_pool$total_poll[poll_pool$Study_Network_id == net_id][1]
    
    poll_seen <- character()
    
    map_dfr(seq_along(plants), function(k){
      
      new_poll <- poll_map %>%
        filter(
          Study_Network_id == net_id,
          Plant_species == plants[k]
        ) %>%
        pull(poll) %>%
        unlist() %>%
        unique()
      
      poll_seen <<- union(poll_seen, new_poll)
      
      tibble(
        Study_Network_id = net_id,
        n_plants = k,
        pollinator_prop = length(poll_seen) / total_poll
      )
    })
  })

unique(acc_curve$Study_Network_id)
library(dplyr)

Poll_mean<-ggplot(summary_curve, aes(x = n_plants)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(aes(y = mean), linewidth = 1.2) +
  labs(
    x = "Number of plant species included",
    y = "Proportion of pollinator richness captured"
  ) +
  theme_minimal()

Poll_mean

library(dplyr)
library(ggplot2)

ggplot(acc_curve,
       aes(n_plants, pollinator_prop)) +
  
  geom_line(
    aes(group = Study_Network_id),
    method = "loess",
    se = FALSE,
    span = 0.4,
    color = "grey55",
    alpha = 0.2,
    linewidth = 0.4
  ) +
  
  theme_minimal()


#random
valid_nets <- data_merge %>%
  group_by(Study_Network_id) %>%
  summarise(
    n_plants_total = n_distinct(Plant_species),
    .groups = "drop"
  ) 


plant_poll_list <- poll_map %>%
  group_by(Study_Network_id, Plant_species) %>%
  summarise(polls = list(unique(unlist(poll))), .groups = "drop")

random_K_curve <- function(net_id, n_perm = 200){
  
  net_plants <- plant_rank %>%
    filter(Study_Network_id == net_id) %>%
    pull(Plant_species) %>%
    unique()
  
  total_poll <- poll_pool$total_poll[
    poll_pool$Study_Network_id == net_id
  ][1]
  
  n_max <- length(net_plants)
  
  map_dfr(1:n_perm, function(p){
    
    map_dfr(1:n_max, function(k){
      
      sampled_plants <- sample(net_plants, k)
      
      poll_k <- poll_map %>%
        filter(
          Study_Network_id == net_id,
          Plant_species %in% sampled_plants   # ⭐这里必须一致
        ) %>%
        pull(poll) %>%
        unlist() %>%
        unique()
      
      tibble(
        Study_Network_id = net_id,
        perm = p,
        n_plants = k,
        pollinator_prop = length(poll_k) / total_poll
      )
    })
  })
}

random_curve <- valid_nets$Study_Network_id %>%
  map_dfr(~ random_K_curve(.x, n_perm = 100))

#####点抽查
random_curve %>%
  filter(n_plants == 10) %>%
  summarise(
    mean_random = mean(pollinator_prop)
  )
acc_curve %>%
  filter(n_plants == 10) %>%
  summarise(
    mean_abundance = mean(pollinator_prop)
  )

length(unique(plant_pool$Study_Network_id))
length(unique(valid_nets$Study_Network_id))

#################
obs_summary <- acc_curve %>%
  group_by(n_plants) %>%
  summarise(
    mean = mean(pollinator_prop),
    lower = quantile(pollinator_prop, 0.05),
    upper = quantile(pollinator_prop, 0.95),
    .groups = "drop"
  ) %>%
  mutate(type = "abundance")

random_summary <- random_curve %>%
  group_by(n_plants) %>%
  summarise(
    mean = mean(pollinator_prop),
    lower = quantile(pollinator_prop, 0.05),
    upper = quantile(pollinator_prop, 0.95),
    .groups = "drop"
  ) %>%
  mutate(type = "random")

plot_df <- bind_rows(obs_summary, random_summary)

# saveRDS(plot_df,"select_vs_radom_accum.rds")
# plot_df<-readRDS("select_vs_radom_accum.rds")

accuplot<-ggplot(plot_df, aes(x = n_plants, y = mean, color = type, fill = type)) +
  
  geom_line(linewidth = 1.2) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2,
              color = NA) +
  
  theme_minimal(base_size = 13) +
  
  labs(
    x = "Number of plant species sampled",
    y = "Proportion of pollinator richness captured",
    color = "",
    fill = ""
  )

Poll_mean_20<-ggplot(plot_df, aes(x = n_plants, y = mean, color = type, fill = type)) +
  
  # 主线
  geom_line(linewidth = 1.2) +
  
  # 置信区间
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  
  # ⭐改颜色
  scale_color_manual(values = c(
    abundance = "#1B9E77",
    random = "grey60"
  ),
  labels = c(
    abundance = "Abundance-based",
    random = "Random"
  )) +
  scale_fill_manual(values = c(
    abundance = "#1B9E77",
    random = "grey60"
  ),
  labels = c(
    abundance = "Abundance-based",
    random = "Random"
  )) +
  
  # ⭐标出 x = 3,5,10
  geom_vline(xintercept = c(3, 5, 10),
             linetype = "dashed",
             color = "black",
             alpha = 0.4) +
  
  # ⭐在这些点加“点”（分别画两条线的对应值）
  geom_point(
    data = plot_df %>% filter(n_plants %in% c(3,5,10)),
    aes(color = type),
    size = 2
  )+
  annotate("text",
           x = c(3,5,10),
           y = Inf,
           label = c("3","5","10"),
           vjust = 1.5,
           size = 5)+
  
  #theme_minimal(base_size = 13) +
  labs(
    x = "Number of plant species sampled",
    y = "Percent of pollinator richness captured",
    color = "Plants selected",
    fill = "Plants selected"
  )+
  scale_y_continuous(
    labels = function(x) x * 100,
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  )+
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(color = "black", size = 12),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(5, 15, 5, 5)
  )
Poll_mean_20
# Combine the three plots vertically
library(cowplot)

ggsave("/Chap1_TargetPlant_to_monitor/result_260723/Poll_mean_20.png",Poll_mean_20,width = 9.5, height = 4.5
       , units = "in", dpi = 300)

# histograms_combined <- plot_grid(histogram_1, histogram_2, ncol = 1, align = "v")
# 
# combined_plot <- plot_grid(histograms_combined, shape_barplot, ncol = 1, rel_heights = c(2, 1))

histograms_combined <- plot_grid(histogram_1, Poll_mean_20, ncol = 1, align = "v")
histograms_combined2<- plot_grid(shape_barplot, shape_class, ncol = 1, align = "v")

histograms_combined3<- plot_grid(histogram_1,shape_barplot,histogram_2,   ncol = 1, align = "v")

final_plot <- plot_grid(histograms_combined, histograms_combined2, ncol = 2, rel_widths = c(2.1, 2.1))

final_plot# size: 1040*670

ggsave("/Chap1_TargetPlant_to_monitor/result_260526/Fig1_four_penal.png", final_plot, width = 10.4, height = 6.7, units = "in", dpi = 300)

ggsave("/Chap1_TargetPlant_to_monitor/result_260526/Distri_plant.png", histograms_combined3, width = 5, height = 6.5, units = "in", dpi = 300)
################3

#  spot check

################
#""    
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
    .groups = "drop"
  ) %>%
  mutate(type = "abundance")

rand_sum <- rand_net %>%
  group_by(n_plants) %>%
  summarise(
    mean = mean(pollinator_prop, na.rm = TRUE),
    lower = quantile(pollinator_prop, 0.05, na.rm = TRUE),
    upper = quantile(pollinator_prop, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(type = "random")

plot_df <- bind_rows(obs_sum, rand_sum)

library(ggplot2)

ggplot(plot_df, aes(x = n_plants, y = mean, color = type, fill = type)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(
    title = paste("Network:", net_id),
    x = "Number of plant species sampled",
    y = "Proportion of pollinator richness captured"
  )


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
    .groups = "drop"
  )

obs_stats <- acc_curve %>%
  group_by(n_plants) %>%
  summarise(
    obs_mean = mean(pollinator_prop, na.rm = TRUE),
    .groups = "drop"
  )
ses_curve <- obs_stats %>%
  inner_join(random_stats, by = "n_plants") %>%
  mutate(
    SES = (obs_mean - rand_mean) / rand_sd
  )

library(ggplot2)

sesplot<-ggplot(ses_curve, aes(x = n_plants, y = SES)) +
  
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  
  geom_line(linewidth = 1.2, color = "#2C7FB8") +
  
  geom_smooth(se = FALSE, color = "black", linewidth = 0.8) +
  
  theme_minimal(base_size = 13) +
  
  labs(
    x = "Number of plant species sampled (K)",
    y = "SES (Abundance-based vs Random)",
    title = "SES Accumulation Curve"
  )

plot_combined <- plot_grid(accuplot, sesplot, ncol = 2, align = "v")
