# ==========================================================
# Table S4: Flower shape categories based on Kugler's methodology
# ==========================================================

library(dplyr)
library(tidyr)
library(flextable)
library(officer)

# Load data
data_count_scaled <- readRDS("data/processed/data_count_scaled_published.rds")
data_interact <- readRDS("data/processed/data_interact_published.rds")
traits <- read.csv("data/processed/merge.trait.csv", header = TRUE, fileEncoding = "UTF-8")

# ==========================================================
# 1. Merge data to get Plant_family information
# ==========================================================

# Extract Plant_family from data_interact (unique pairs)
plant_family_map <- data_interact %>%
  select(Plant_original_name, Plant_family) %>%
  distinct() %>%
  rename(Plant_species = Plant_original_name)

# Merge data_count_scaled with family information (keep all rows from data_count_scaled)
data_with_family <- data_count_scaled %>%
  left_join(
    plant_family_map,
    by = c("Plant_species" = "Plant_species")
  ) %>%
  # Replace NA family with "Unknown"
  mutate(Plant_family = ifelse(is.na(Plant_family), "Unknown", Plant_family))

# ==========================================================
# 2. Verification: Check which flower shapes have family info
# ==========================================================

cat("\nFlower shape categories and coverage:\n")
print(
  data_with_family %>%
    filter(!is.na(flw_shape_revised)) %>%
    group_by(flw_shape_revised) %>%
    summarise(
      total_records = n(),
      with_family = sum(Plant_family != "Unknown"),
      without_family = sum(Plant_family == "Unknown"),
      .groups = "drop"
    ) %>%
    arrange(desc(total_records))
)

# ==========================================================
# 3. Create Table S5 - Simplified version
# ==========================================================

# ==========================================================
# 1. Get plant info from data_count_scaled
#    (includes families and all surveyed plants)
# ==========================================================

plant_info <- data_count_scaled %>%
  filter(!is.na(flw_shape_revised)) %>%
  left_join(
    data_interact %>%
      select(Plant_accepted_name, Plant_family) %>%
      distinct(),
    by = c("Plant_species" = "Plant_accepted_name")
  ) %>%
  mutate(Plant_family = ifelse(is.na(Plant_family), "Unknown", Plant_family))


plant_species_count <- plant_info %>%
  group_by(flw_shape_revised) %>%
  summarise(
    n_plant_species = n_distinct(Plant_species),
    
    # Extract top 3 families (from all surveyed plants, not just those with interactions)
    main_families = {
      families <- sort(table(Plant_family[Plant_family != "Unknown"]), 
                       decreasing = TRUE)
      if(length(families) >= 3) {
        paste(names(families[1:3]), collapse = ", ")
      } else if(length(families) > 0) {
        paste(names(families), collapse = ", ")
      } else {
        "Unknown"
      }
    },
    
    .groups = "drop"
  )


# ==========================================================
# 2. Get interaction network info from data_interact
# ==========================================================

interaction_counts <- data_interact %>%
  filter(Interaction_addup > 0) %>%  # 只统计有互作的记录
  filter(!is.na(flw_shape_revised)) %>%
  group_by(flw_shape_revised) %>%
  summarise(
    n_networks = n_distinct(Study_Network_id),
    
    # Extract top 3 most abundant species (by interactions)
    example_species = paste(
      names(sort(table(Plant_accepted_name), decreasing = TRUE)[1:3]),
      collapse = "; "
    ),
    
    n_interactions = n(),
    .groups = "drop"
  )


# ==========================================================
# 3. Merge both information sources
# ==========================================================

table_s5 <- plant_species_count %>%
  left_join(
    interaction_counts,
    by = "flw_shape_revised"
  ) %>%
  
  # Rank by number of plant species
  arrange(desc(n_plant_species)) %>%
  
  # Select and rename columns
  select(
    `Flower Shape Category` = flw_shape_revised,
    `Plant Species (n)` = n_plant_species,
    `Networks (n)` = n_networks,
    `Main Families (top 3)` = main_families,
    `Example Species` = example_species
  )

print(table_s5)

# ==========================================================
# 2. Export as CSV
# ==========================================================

write.csv(
  table_s5,
  "data/processed/Table_S5_Flower_Shape_Categories.csv",
  row.names = FALSE
)

# ==========================================================
# 3. Create formatted Word table
# ==========================================================

ft_table_s5 <- flextable(table_s5) %>%
  
  theme_booktabs() %>%
  
  # 字体大小 - 补充表格通常用9-10pt
  fontsize(size = 9, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  
  bold(part = "header") %>%
  
  # 对齐方式
  align(j = 1, align = "left", part = "all") %>%
  align(j = 2:5, align = "center", part = "all") %>%
  
  # 列宽设置 - 根据内容合理设置，不强行占满页面
  # 总宽度约 7.5 英寸（标准Word页面 8.5" - 两侧边距 1"）
  width(j = 1, width = 2.2) %>%      # Flower Shape Category
  width(j = 2, width = 1.2) %>%      # Plant Species (n)
  width(j = 3, width = 1.0) %>%      # Networks (n)
  width(j = 4, width = 2.0) %>%      # Main Families (top 3)
  width(j = 5, width = 1.8)          # Example Species


# Create Word document
doc_table_s5 <- read_docx() %>%
  
  body_add_par(
    "Table S5",
    style = "heading 2"
  ) %>%
  
  body_add_par(
    paste0(
      "Flower shape categories following Kugler's methodology ",
      "(Kugler 1955, 1970) and BiolFlor database classification. ",
      "The table shows flower shape categories ranked by frequency, ",
      "along with the number of plant species, networks where each category ",
      "occurs, the three most common plant families, and example species."
    ),
    style = "Normal"
  ) %>%
  
  body_add_par("", style = "Normal") %>%
  
  body_add_flextable(ft_table_s5)

# Save Word document
print(
  doc_table_s5,
  target = "data/processed/Table_S5_Flower_Shape_Categories.docx"
)

cat("\n✓ Table S5 created and saved!\n")
cat("Files saved:\n")
cat("  - Table_S5_Flower_Shape_Categories.csv\n")
cat("  - Table_S5_Flower_Shape_Categories.docx\n")