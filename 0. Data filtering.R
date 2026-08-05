##The full script
##require data: Interaction_data.rds,Flower_counts.rds,test.merge.trait.csv

### METHODS
#### 1. data filtering
getwd()
library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(ggplot2)
#read data
metadata<-readRDS("Interaction_data_published.rds")%>%
  mutate(Study_Network_id = paste(Study_id, Network_id, sep = "_")) #The EuPPollNet interaction data

length(unique(metadata$Study_Network_id))  # 1630

meta_count<-readRDS("Flower_counts_published.rds")#The EuPPollNet flower data
colnames(metadata)

test<-meta_count%>%filter(Study_id == "20_Hoiss")
sum(metadata$Interaction) #623476
nrow(metadata) #623476
########
# Typo revision
library(dplyr)
library(stringr)

meta_count <- meta_count %>%
  mutate(
    # 去掉末尾的 sp./Sp./SP./1/?/spec.
    Plant_species = str_remove(Plant_species, "\\s+(sp\\.|Sp\\.|SP\\.|1|3|\\?|spec\\.|/|spp\\.?|spps\\.?|sp\\s*\\d+)$"),
    # 去掉中间的 " sp."
    Plant_species = str_replace(Plant_species, "\\s+sp\\.\\s+", " "),
    # 删除末尾的 sp / Sp / spp
    Plant_species = str_remove(Plant_species, "\\s+(sp|Sp|spp)$"),
    # 删除末尾的下划线
    Plant_species = str_remove(Plant_species, "_$")
  ) %>%
  # 去掉特定无效值
  filter(!is.na(Plant_species), !Plant_species %in% c("NO FLOWERS", "[Nothing] [Nothing]")) %>%
  mutate(
    Plant_species = case_when(
      Plant_species == "Ribes nigrum.andega"       ~ "Ribes nigrum",
      Plant_species == "Ribes nigrum.blackdown"   ~ "Ribes nigrum",
      Plant_species == "Capsella bursa.pastoris"  ~ "Capsella bursa-pastoris",
      Plant_species == "Weigela kosteriana.variegata" ~ "Weigela kosteriana",
      Plant_species == "Silene flos.cuculi"       ~ "Silene flos-cuculi",
      Plant_species == "Rorippa x.anceps"        ~ "Rorippa x anceps",
      Plant_species == "Arenaria serpyllifolia.aggr" ~ "Arenaria serpyllifolia",
      Plant_species == "Viola tricolor.aggr"     ~ "Viola tricolor",
      Plant_species == "Citrus x.limon"          ~ "Citrus x limon",
      Plant_species == "Berberis thunbergii.atropurpurea" ~ "Berberis thunbergii",
      Plant_species == "Oenothera parviflora.aggr" ~ "Oenothera parviflora",
      Plant_species == "Linum nervosum.perenne"  ~ "Linum nervosum",
      Plant_species == "Bituminaria bituminosa (L.) C.H.Stirt." ~ "Bituminaria bituminosa",
      Plant_species == "Vaccinium vitis_idaea" ~ "Vaccinium vitis-idaea",
      TRUE ~ Plant_species  # 其他行保持原样
    )
  )%>%
  # 创建新列只保留属名 + 种加词
  mutate(
    n_words = str_count(Plant_species, "\\S+"),  # 计算词数
    Plant_species = if_else(
      n_words >= 2,
      str_c(word(Plant_species, 1), word(Plant_species, 2), sep = " "),  # 前两个词
      Plant_species  # 只有一个词则保留
    )
  ) %>%
  select(-n_words)  %>%
  mutate(
    Plant_species = str_remove(Plant_species, "/$"),
    Plant_species = str_squish(Plant_species)  # 去掉尾部多余空格
  )
#########################################

#1.1 unify plant name and merge the interaction data and flower data
species_correct<-read.csv("species_checked_wof.csv")%>%
  select(WOF_name,original_name)%>%
  mutate(Plant_species = original_name)

meta_count <- meta_count %>%
  left_join(species_correct %>% select(Plant_species, WOF_name), by = "Plant_species") %>%
  mutate(Plant_species = coalesce(WOF_name, Plant_species)) %>%  # 如果有对应 WOF_name 就替换
  select(-WOF_name)  # 去掉临时列

##########
# Site sampled in multiple years should be seperated
########

year_flag <- meta_count %>%
  mutate(
    Study_Network_id_noyear = paste(Study_id, Network_id, sep = "_")
  ) %>%
  group_by(Study_Network_id_noyear) %>%
  summarise(
    has_year = any(!is.na(Year)),
    .groups = "drop"
  )


metadata <- metadata %>%
  mutate(
    Year = year(Date),
    Study_Network_id_noyear = paste(Study_id, Network_id, sep = "_")
  ) %>%
  left_join(year_flag, by = "Study_Network_id_noyear") %>%
  mutate(
    Study_Network_id = case_when(
      
      # ❗ 两边都有 year → 用 site-year
      has_year & !is.na(Year) ~
        paste(Study_id, Network_id, Year, sep = "_"),
      
      # ❗ flower data 没 year → 强制降级
      TRUE ~
        Study_Network_id_noyear
    )
  )

meta_count <- meta_count %>%
  mutate(
    Study_Network_id_noyear = paste(Study_id, Network_id, sep = "_")
  ) %>%
  left_join(year_flag, by = "Study_Network_id_noyear") %>%
  mutate(
    Study_Network_id = case_when(
      has_year ~ paste(Study_id, Network_id, Year, sep = "_"),
      TRUE ~ Study_Network_id_noyear
    )
  ) %>%
  filter(!is.na(Flower_count))


#check the record
unique(metadata$Study_id)#54
unique(metadata$Study_Network_id)#1740
unique(meta_count$Study_Network_id)#1740
b<-meta_count%>%
  filter(!is.na(Flower_count))%>%
  distinct(Study_Network_id)
study_network_ids<-b$Study_Network_id

#site_id selected
#select all the sites with plant data
data_interact <- metadata[metadata$Study_Network_id %in% study_network_ids, ] #data_interact是筛选出的，含有plant_data的networks

data_count <- meta_count%>%
  filter(Flower_count!=0)%>%
  filter(!is.na(Flower_count))

unique(data_count$Study_Network_id)

study_network_id_count <- unique(data_count$Study_Network_id)
study_network_id_interact <- unique(data_interact$Study_Network_id)

# Match two datasets
count_not_in_interact <- setdiff(study_network_id_count, study_network_id_interact)
interact_not_in_count <- setdiff(study_network_id_interact, study_network_id_count)

count_not_in_interact #20
interact_not_in_count #0

unique(metadata$Study_Network_id)
unique(meta_count$Study_Network_id)


# sampling_check_site <- meta_count %>%
#   distinct(Study_Network_id, Year, Month) %>%
#   count(Study_Network_id, Year, name = "n_months_sampled") %>%
#   filter(n_months_sampled > 1) %>%
#   arrange(desc(n_months_sampled))
# 
# sampling_check_site
# 
# multi_month_detail <- meta_count %>%
#   distinct(Study_Network_id, Year, Month) %>%
#   group_by(Study_Network_id, Year) %>%
#   summarise(
#     months = paste(sort(unique(Month)), collapse = ", "),
#     n_months = n(),
#     .groups = "drop"
#   ) %>%
#   filter(n_months > 1) %>%
#   arrange(desc(n_months))
# 
# multi_month_detail
# 
# year_detail <- meta_count %>%
#   distinct(Study_Network_id, Year) %>%
#   group_by(Study_Network_id) %>%
#   summarise(
#     years = paste(sort(unique(Year)), collapse = ", "),
#     n_years = n(),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(n_years))
# 
# year_detail
# 
# 
# year_detail_in <- metadata %>%
#   select(Study_Network_id, Year) %>%
#   distinct() %>%
#   group_by(Study_Network_id) %>%
#   summarise(
#     years = paste(sort(unique(Year)), collapse = ", "),
#     n_years = n_distinct(Year),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(n_years))
# 
# 
# year_detail_in
# # 
# colnames(metadata)
# colnames(meta_count)
# multi_month_detail_in <- metadata %>%
#   distinct(Study_Network_id, Year, Date) %>%
#   group_by(Study_Network_id, Year) %>%
#   summarise(
#     months = paste(sort(unique(Date)), collapse = ", "),
#     n_months = n(),
#     .groups = "drop"
#   ) %>%
#   filter(n_months > 1) %>%
#   arrange(desc(n_months))
# 
# multi_month_detail_in
# 
# 
# test<-metadata%>%filter(Study_Network_id == "22_Kallnik_h1")


#data filtering and scaling 
#calculate average plant abundance for each plant species in the plant survey data. 
data_count_scaled <- meta_count%>%
  #filter(!Study_Network_id%in%count_not_in_interact) %>%
  group_by(Plant_species,Study_Network_id, Year, Month, Day) %>%
  mutate(Flower_count = sum(Flower_count, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Plant_species, Study_Network_id) %>% # Year seperate
  filter (!Flower_count == 0) %>%
  mutate(Flower_count_scaled = mean(Flower_count, na.rm = TRUE))%>%
  distinct(Flower_count_scaled, Plant_species, Study_Network_id, .keep_all = TRUE)%>%
  group_by(Flower_data_merger) %>%
  filter(n() == 1) %>%
  ungroup()


length(unique(data_count_scaled$Study_Network_id))#1055
unique(data_interact$Study_id)#36


#Removed networks with no separate plant survey. 
data_count_scaled<-data_count_scaled%>%
  filter(!Study_Network_id%in%count_not_in_interact)
length(unique(data_count_scaled$Study_Network_id))#1035
unique(data_interact$Study_id)#36

length(unique(data_count_scaled$Study_Network_id))#36

test<-data_count_scaled%>%filter(Study_id == "20_Hoiss")
#Now we are left with 36 studies with 1035 networks.


##########################################

#=========================screening==============================================
#######
# 检查剩余网络数
length(unique(data_count_scaled$Study_Network_id))# 1035

# remove all wind flowers from interaction data
traits<-read.csv("test.merge.trait.csv", header = TRUE, fileEncoding = "UTF-8")
# 
data_interact_trait<-left_join(data_interact, traits, by = "Plant_accepted_name")
data_count_scaled_trait<-left_join(data_count_scaled, traits, by = c("Plant_species" = "Plant_accepted_name"))

data_interact<-data_interact_trait%>%
  filter(!flw_shape_revised == "wind flowers")
  # group_by(Plant_accepted_name,Study_Network_id, Pollinator_accepted_name)%>%
  # mutate(Interaction_addup=sum(Interaction))%>%
  # ungroup()%>%
  # distinct(Interaction_addup, Plant_accepted_name,Study_Network_id, Pollinator_accepted_name, .keep_all = TRUE)

data_count_scaled<-data_count_scaled_trait%>%
  filter(!flw_shape_revised == "wind flowers")

length(unique(data_interact$Study_Network_id))
length(unique(data_count_scaled$Study_Network_id))# 1035
length(unique(data_count_scaled$Study_id))# 36


#1.2 checked if:
#•	less than 75% plant names in the interaction data are not in the plant list.  
#•	Plant species involved in the most interactions in the network are not in the plant list.
###(1)•	Remove -more than 25% plant names in the interaction data are not in the plant list
##To check how many records doesn't merge
colnames(data_interact)
colnames(data_count_scaled)
length(unique(data_count_scaled$Study_Network_id))# 1030
library(dplyr)
library(stringr)
plant_inter<-data_interact%>%select(Study_Network_id,Plant_original_name)%>%distinct()
plant_flower<-data_count_scaled%>%select(Study_Network_id,Plant_species)%>%distinct()

# 清理 plant_inter
plant_inter <- plant_inter %>%
  mutate(
    Plant_original_name = str_trim(Plant_original_name), # 去掉前后空格
    Plant_original_name = str_to_lower(Plant_original_name) # 全小写
  )

# 清理 plant_flower
plant_flower <- data_count_scaled %>%
  select(Study_Network_id, Plant_species, Flower_count) %>%
  distinct() %>%
  na.omit() %>%
  mutate(
    Plant_species = str_trim(Plant_species),
    Plant_species = str_to_lower(Plant_species)
  )

# 再 join
plant_merge <- plant_inter %>%
  left_join(plant_flower, 
            by = c("Study_Network_id", "Plant_original_name" = "Plant_species"))

# 检查缺失情况
NA1_proportion <- plant_merge %>%
  group_by(Study_Network_id) %>%
  summarise(
    total_plants = n(),
    missing_plants = sum(is.na(Flower_count)),
    proportion_missing = missing_plants / total_plants
  )

# 筛掉 >25% 缺失的网络
dropNetwork <- NA1_proportion %>%
  filter(proportion_missing > 0.25) %>%
  pull(Study_Network_id)

dropNetwork #452

# 更新数据
data_interact <- data_interact %>%
  filter(!Study_Network_id %in% dropNetwork)

data_count_scaled <- data_count_scaled %>%
  filter(!Study_Network_id %in% dropNetwork)


# 检查剩余网络数
length(unique(data_count_scaled$Study_Network_id))#583
length(unique(data_count_scaled$Study_id))#31
length(unique(data_interact$Study_Network_id))#581

#1.3  looked for studies that contain few interactions overall (data poor studies).  
###(1)•	remove <30 interactions
library(ggplot2)
library(dplyr)
library(cowplot)  

plant_visit <- data_interact %>%
  group_by(Study_Network_id) %>%
  summarise(plant_visit_time = sum(Interaction, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(flag = plant_visit_time <= 30)
  
main_data <- plant_visit %>% filter(plant_visit_time <= 500)
tail_data <- plant_visit %>% filter(plant_visit_time > 500)


# 主图
p_main <- ggplot(main_data, aes(x = plant_visit_time, fill = flag)) +
  geom_histogram(binwidth = 10, color = "black") +
  
  scale_fill_manual(values = c(
    "TRUE" = "orange",
    "FALSE" = "#66C2A5"
  )) +
  
  labs(
    x = "Total interactions per network (≤500)",
    y = "Number of networks",
    fill = "Low interaction (≤30)"
  ) +
  
  theme_classic() +   # ✅ 先放主题
  
  theme(
    legend.position = "bottom"
  ) +
  
  guides(
    fill = guide_legend(override.aes = list(size = 5))
  )

# inset 小图
p_tail <- ggplot(tail_data, aes(x = plant_visit_time)) +
  geom_histogram(binwidth = 100, fill = "grey50", color = "black") +
  theme_minimal(base_size = 8) +
  labs(x = ">500", y = NULL)

# 拼 inset（放右上角）
histogram_1 <- ggdraw() +
  draw_plot(p_main) +
  draw_plot(p_tail, x = 0.55, y = 0.5, width = 0.4, height = 0.4)



print(histogram_1)#750*350

ggsave("/Chap1_TargetPlant_to_monitor/result_260526/hist1.png", histogram_1, width = 7.5, height = 3.5, units = "in", dpi = 300)


## filtering
data_interact<-data_interact%>%
  left_join(
    plant_visit,by="Study_Network_id")%>%
  filter(plant_visit_time>=30)

data_count_scaled<-data_count_scaled%>%
  left_join(plant_visit,by="Study_Network_id")%>%
  filter(plant_visit_time>=30)

length(unique(data_interact$Study_Network_id))#left with 459
unique(data_interact$Study_id)#  studies 31

#drop 124 Networks (<30 interactions)
#left with  459 Networks

###(2)Removed studies with fewer than 10 plant species and fewer than 10 pollinator species.

#figures
data_count_scaled_species <- data_count_scaled %>%
  filter(!is.na(Plant_species)) %>%       # 排除 NA
  filter(str_detect(Plant_species, " "))  # 保留含空格的名字（双名）

plant_diversity<-data_count_scaled_species%>%
  filter(!is.na(Plant_species)) %>%
  group_by(Study_Network_id) %>%
  summarise(plant_sp_number = n_distinct(Plant_species,na.rm = TRUE)) %>%
  mutate(flag = plant_sp_number <= 10)
nrow(plant_diversity)# 542 461

saveRDS(plant_diversity,"plant_diversity_461networks.rds")


breaks_2 <- seq(0, max(plant_diversity$plant_sp_number), by = 10)
unique(plant_diversity$Study_Network_id)#left with 350 Networks

histogram_2 <- ggplot(
  plant_diversity,
  aes(x = plant_sp_number, fill = flag)
) +
  
  geom_histogram(
    binwidth = 2,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "TRUE" = "orange",
      "FALSE" = "#66C2A5"
    )
  ) +
  
  guides(
    fill = guide_legend(override.aes = list(size = 5))
  ) +
  
  labs(
    x = "Number of plant species",
    y = "Number of networks",
    fill = "Low richness (≤10)"
  ) +
  
  theme_classic() +   # ✅ 和 p_main 一致
  
  theme(
    legend.position = "bottom"   # ✅ 一致
  )


print(histogram_2)#750*350

ggsave("/Chap1_TargetPlant_to_monitor/result_260526/hist2.png", histogram_2, width = 7.5, height = 3.5, units = "in", dpi = 300)


few_sp <- data_interact %>%
  filter(!is.na(Plant_accepted_name)) %>%       # 排除 NA
  filter(str_detect(Plant_accepted_name, " ")) %>%  # 保留含空格的名字（双名）
  group_by(Study_Network_id) %>%
  summarize(
    plant_sp_number = n_distinct(Plant_accepted_name, na.rm = TRUE),
    pollinator_sp_number = n_distinct(Pollinator_accepted_name, na.rm = TRUE)
  ) %>%
  filter(plant_sp_number<10| pollinator_sp_number < 10)

data_interact <-data_interact%>%
  filter(!Study_Network_id%in%few_sp$Study_Network_id)
data_count_scaled <-data_count_scaled%>%
  filter(!Study_Network_id%in%few_sp$Study_Network_id)%>%
  filter(!is.na(Plant_species)& Flower_count != 0) %>%       # 排除 NA
  filter(str_detect(Plant_species, " "))  # 保留含空格的名字（双名）, 去掉所有只鉴定到属的记录


# select_network<-unique(data_interact$Study_Network_id)
length(unique(data_interact$Study_id))
unique(data_interact$Study_Network_id)#332个网络 270，30个研究 27

length(unique(data_count_scaled$Plant_species))# 1396 1029


###(2)Removed studies with fewer than 10 plant species in plant survey data.
few_sp_in_plant_survey <- data_count_scaled %>%
  group_by(Study_Network_id) %>%
  summarize(
    plant_sp_number = n_distinct(Plant_species, na.rm = TRUE)
  ) %>%
  filter(plant_sp_number<10)

few_sp_in_plant_survey ## 4 record


data_interact <-data_interact%>%
  filter(!Study_Network_id%in%few_sp_in_plant_survey$Study_Network_id)
data_count_scaled <-data_count_scaled%>%
  filter(!Study_Network_id%in%few_sp_in_plant_survey$Study_Network_id)%>%
  filter(!is.na(Plant_species)& Flower_count != 0) %>%       # 排除 NA
  filter(str_detect(Plant_species, " "))  # 保留含空格的名字（双名）, 去掉所有只鉴定到属的记录


# select_network<-unique(data_interact$Study_Network_id)
unique(data_interact$Study_id)
unique(data_interact$Study_Network_id)

unique(data_count_scaled$Plant_species)

#### summarize
select_network<-unique(data_interact$Study_Network_id)
unique(data_interact$Study_id)
unique(data_interact$Study_Network_id)#最终选择328个网络，30个研究
unique(data_count_scaled$Study_Network_id)
unique(data_count_scaled$Plant_species)# 1396
unique(data_interact$Plant_accepted_name)# 1028

data_interact<-data_interact%>%
  group_by(Plant_accepted_name,Study_Network_id, Pollinator_accepted_name)%>%
  mutate(Interaction_addup=sum(Interaction))%>%
  ungroup()%>%
  distinct(Interaction_addup, Plant_accepted_name,Study_Network_id, Pollinator_accepted_name, .keep_all = TRUE)

final_study_id<-unique(data_interact$Study_id)

##################################

#-----------Check plant sampling methods

###################################

colnames(data_count_scaled)
Plant_method<-data_count_scaled%>%
  select("Study_id","Network_id",  "Units","Comments" )%>%
  distinct()

uniq<-Plant_method%>%
  select("Study_id","Units" )%>%
  distinct()
unique(uniq$Study_id)

library(dplyr)
library(flextable)
library(officer)

plant_unit <- read.csv("plant_sampling_unit.csv") %>%
  select("Study_id","Flower.sampling.methods")%>%
  filter (Study_id %in% final_study_id)
setdiff(final_study_id, plant_unit$`Study ID`)
names(plant_unit) <- c("Study ID", "Flower sampling methods")

ft <- flextable(plant_unit) %>%
  theme_booktabs() %>%                     # 三线表
  fontsize(size = 11, part = "all") %>%    # 字号
  font(fontname = "Arial", part = "all") %>%
  
  # 期刊风格列宽
  width(j = "Study ID", width = 1.5) %>%
  width(j = "Flower sampling methods", width = 7) %>%
  
  # 对齐
  align(j = "Study ID", align = "center") %>%
  align(j = "Flower sampling methods", align = "left")

doc <- read_docx() %>%
  body_add_par("Table 1. Flower sampling methods per study", style = "heading 2") %>%
  body_add_flextable(ft)

print(doc, target = "./result_260526/plant_sampling_unit.docx")

#save
 saveRDS(data_count_scaled,"data_count_scaled_published_0526.rds")
 saveRDS(data_interact,"data_interact_published_0526.rds")
