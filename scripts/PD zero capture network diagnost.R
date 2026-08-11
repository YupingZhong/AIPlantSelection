############################################################
# Function: diagnose PD failures
############################################################

make_PD_diagnostic <- function(
    PD_percent,
    PD_selected_df,
    top_n
){
  
  # --------------------------------------------------
  # Networks with PD = 0
  # --------------------------------------------------
  
  PD_fail <- PD_percent %>%
    filter(
      percentage == 0 |
        is.na(percentage)
    ) %>%
    select(
      Study_Network_id,
      PD_percentage = percentage,
      PD_pollinator_count = pollinator_count
    )
  
  # --------------------------------------------------
  # PD plants
  # --------------------------------------------------
  
  PD_species <- PD_selected_df %>%
    group_by(
      Study_Network_id
    ) %>%
    summarise(
      PD_plants =
        paste(
          Plant_species,
          collapse = "; "
        ),
      .groups = "drop"
    )
  
  # --------------------------------------------------
  # Abundance plants
  # --------------------------------------------------
  
  Abun_species <- data_count_scaled %>%
    group_by(
      Study_Network_id
    ) %>%
    arrange(
      desc(Flower_count_scaled),
      .by_group = TRUE
    ) %>%
    slice_head(
      n = top_n
    ) %>%
    summarise(
      Abun_plants =
        paste(
          Plant_species,
          collapse = "; "
        ),
      .groups = "drop"
    )
  
  # --------------------------------------------------
  # Abundance pollinator capture
  # --------------------------------------------------
  
  Abun_capture <- data_count_scaled %>%
    group_by(
      Study_Network_id
    ) %>%
    arrange(
      desc(Flower_count_scaled),
      .by_group = TRUE
    ) %>%
    slice_head(
      n = top_n
    ) %>%
    select(
      Study_Network_id,
      Plant_species
    ) %>%
    left_join(
      data_interact,
      by = c(
        "Study_Network_id",
        "Plant_species" =
          "Plant_original_name"
      )
    ) %>%
    filter(
      Interaction_addup > 0
    ) %>%
    group_by(
      Study_Network_id
    ) %>%
    summarise(
      Abun_pollinator_count =
        n_distinct(
          Pollinator_accepted_name
        ),
      .groups = "drop"
    )
  
  # --------------------------------------------------
  # Abundance percentage column
  # --------------------------------------------------
  
  abun_col <- paste0(
    "percentage_Abun",
    top_n
  )
  
  Abun_percent <- result_all %>%
    select(
      Study_Network_id,
      !!sym(abun_col)
    )
  
  names(Abun_percent)[2] <-
    "Abun_percentage"
  
  # --------------------------------------------------
  # Total pollinator richness
  # --------------------------------------------------
  
  Total_pollinator <- total_number %>%
    select(
      Study_Network_id,
      total_pollinator_count
    )
  
  # --------------------------------------------------
  # Network size
  # --------------------------------------------------
  
  Network_size <- data.frame(
    Study_Network_id = unique(data_count_scaled$Study_Network_id)
  ) %>%
    
    left_join(
      
      data_count_scaled %>%
        group_by(
          Study_Network_id
        ) %>%
        summarise(
          n_plants_count =
            n_distinct(
              Plant_species
            ),
          .groups = "drop"
        ),
      
      by = "Study_Network_id"
    ) %>%
    
    left_join(
      
      data_interact %>%
        group_by(
          Study_Network_id
        ) %>%
        summarise(
          n_plants_interact =
            n_distinct(
              Plant_accepted_name
            ),
          .groups = "drop"
        ),
      
      by = "Study_Network_id"
    )
  
  # --------------------------------------------------
  # Final table
  # --------------------------------------------------
  
  PD_fail %>%
    left_join(
      Total_pollinator,
      by = "Study_Network_id"
    ) %>%
    left_join(
      PD_species,
      by = "Study_Network_id"
    ) %>%
    left_join(
      Abun_species,
      by = "Study_Network_id"
    ) %>%
    left_join(
      Abun_capture,
      by = "Study_Network_id"
    ) %>%
    left_join(
      Abun_percent,
      by = "Study_Network_id"
    ) %>%
    arrange(
      desc(Abun_percentage)
    )%>%
    left_join(
      Network_size,
      by = "Study_Network_id"
    )
}

############################################################
# Generate diagnostics
############################################################

PD10_diagnostic <- make_PD_diagnostic(
  PD_percent = PD_percent10,
  PD_selected_df = PD_selected$PD_top10,
  top_n = 10
)

PD5_diagnostic <- make_PD_diagnostic(
  PD_percent = PD_percent5,
  PD_selected_df = PD_selected$PD_top5,
  top_n = 5
)

PD3_diagnostic <- make_PD_diagnostic(
  PD_percent = PD_percent3,
  PD_selected_df = PD_selected$PD_top3,
  top_n = 3
)

############################################################
# View
############################################################

View(PD10_diagnostic)
View(PD5_diagnostic)
View(PD3_diagnostic)
