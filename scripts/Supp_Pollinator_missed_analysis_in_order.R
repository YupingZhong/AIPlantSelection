###############################################################################
# Supplementary analysis: pollinator species missed when monitoring only the
# 10 most abundant plant species
#
# Primary unit of analysis: pollinator species x network
#
# Required input objects/files:
#   1. data/processed/data_interact_published.rds
#      Complete interaction data.
#   2. data/processed/selected_plant_result_abun10.rds
#      Interactions retained when only the 10 most abundant plant species are
#      monitored.
#
# Main outputs:
#   - species-level missed-network summary
#   - family-level species-by-network summary
#   - order-level species-by-network summary
#   - bee-species summary
#   - Supplementary Fig. S6 (counts) and S7 (proportional miss rates)
#   - optional legacy family-by-network disappearance summary
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

# -----------------------------------------------------------------------------
# 0. User settings
# -----------------------------------------------------------------------------

project_dir <- "D:/Chap1_TargetPlant_to_monitor/Project_PlantSelection"
interaction_file <- file.path(project_dir, "data/processed/data_interact_published.rds")
top10_file <- file.path(project_dir, "data/processed/selected_plant_result_abun10.rds")
output_dir <- "/Chap1_TargetPlant_to_monitor/result_260723"

target_orders <- c(
  "Diptera",
  "Lepidoptera",
  "Hymenoptera",
  "Coleoptera"
)

# Families generally treated as bees. Check this list against the taxonomy and
# scope used in the manuscript before publication.
bee_families <- c(
  "Andrenidae",
  "Apidae",
  "Colletidae",
  "Halictidae",
  "Megachilidae",
  "Melittidae",
  "Stenotritidae"
)

# Primary species-level analysis excludes records identified only to genus or a
# higher rank. Set FALSE if single-word accepted names are valid species IDs in
# the underlying dataset.
require_binomial_species_names <- TRUE

# These thresholds are used only for a descriptive sensitivity table. They do
# not affect the main results or figures.
minimum_networks_for_consistent_miss <- 5L
minimum_miss_rate_for_consistent_miss <- 50

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Read and validate data
# -----------------------------------------------------------------------------

if (!file.exists(interaction_file)) {
  stop("Complete interaction file not found: ", interaction_file)
}

if (!file.exists(top10_file)) {
  stop("Top-10 result file not found: ", top10_file)
}

data_interact <- readRDS(interaction_file)
result10 <- readRDS(top10_file)

required_total_columns <- c(
  "Study_Network_id",
  "Pollinator_order",
  "Pollinator_family",
  "Pollinator_accepted_name"
)

required_top10_columns <- c(
  "Study_Network_id",
  "Pollinator_accepted_name"
)

missing_total_columns <- setdiff(required_total_columns, names(data_interact))
missing_top10_columns <- setdiff(required_top10_columns, names(result10))

if (length(missing_total_columns) > 0L) {
  stop(
    "Columns missing from data_interact: ",
    paste(missing_total_columns, collapse = ", ")
  )
}

if (length(missing_top10_columns) > 0L) {
  stop(
    "Columns missing from result10: ",
    paste(missing_top10_columns, collapse = ", ")
  )
}

# Clean only whitespace. Do not silently alter capitalization or taxonomy.
data_interact_clean <- data_interact %>%
  mutate(
    Study_Network_id = str_squish(as.character(Study_Network_id)),
    Pollinator_order = str_squish(as.character(Pollinator_order)),
    Pollinator_family = str_squish(as.character(Pollinator_family)),
    Pollinator_accepted_name = str_squish(
      as.character(Pollinator_accepted_name)
    )
  )

result10_clean <- result10 %>%
  mutate(
    Study_Network_id = str_squish(as.character(Study_Network_id)),
    Pollinator_accepted_name = str_squish(
      as.character(Pollinator_accepted_name)
    )
  )

# Check whether one accepted name has conflicting order/family assignments.
taxonomy_conflicts <- data_interact_clean %>%
  filter(
    !is.na(Pollinator_accepted_name),
    Pollinator_accepted_name != ""
  ) %>%
  distinct(
    Pollinator_accepted_name,
    Pollinator_order,
    Pollinator_family
  ) %>%
  count(Pollinator_accepted_name, name = "Taxonomic_assignments") %>%
  filter(Taxonomic_assignments > 1L)

if (nrow(taxonomy_conflicts) > 0L) {
  write.csv(
    taxonomy_conflicts,
    file.path(output_dir, "CHECK_taxonomy_conflicts.csv"),
    row.names = FALSE
  )
  stop(
    "Some accepted species names have multiple order/family assignments. ",
    "Resolve CHECK_taxonomy_conflicts.csv before rerunning the analysis."
  )
}

# -----------------------------------------------------------------------------
# 2. Construct the complete species x network observation universe
# -----------------------------------------------------------------------------

species_network_total <- data_interact_clean %>%
  filter(
    !is.na(Study_Network_id),
    Study_Network_id != "",
    !is.na(Pollinator_order),
    Pollinator_order != "",
    !is.na(Pollinator_family),
    Pollinator_family != "",
    !is.na(Pollinator_accepted_name),
    Pollinator_accepted_name != ""
  ) %>%
  distinct(
    Study_Network_id,
    Pollinator_order,
    Pollinator_family,
    Pollinator_accepted_name
  )

if (require_binomial_species_names) {
  species_network_total <- species_network_total %>%
    filter(str_count(Pollinator_accepted_name, "\\S+") >= 2L)
}

# A species x network observation is captured if that same accepted species is
# present in result10 for the same network. Taxonomy is taken from the complete
# dataset to avoid mismatches caused by metadata differences between files.
species_network_captured_keys <- result10_clean %>%
  filter(
    !is.na(Study_Network_id),
    Study_Network_id != "",
    !is.na(Pollinator_accepted_name),
    Pollinator_accepted_name != ""
  ) %>%
  distinct(Study_Network_id, Pollinator_accepted_name)

# Diagnostic: top-10 keys that do not match the cleaned complete-data universe.
# A nonzero result can indicate inconsistent names/network IDs or deliberate
# exclusions such as genus-only identifications.
unmatched_top10_keys <- species_network_captured_keys %>%
  anti_join(
    species_network_total %>%
      distinct(Study_Network_id, Pollinator_accepted_name),
    by = c("Study_Network_id", "Pollinator_accepted_name")
  )

write.csv(
  unmatched_top10_keys,
  file.path(output_dir, "CHECK_unmatched_top10_species_network_keys.csv"),
  row.names = FALSE
)

species_network_status <- species_network_total %>%
  left_join(
    species_network_captured_keys %>% mutate(Captured = TRUE),
    by = c("Study_Network_id", "Pollinator_accepted_name")
  ) %>%
  mutate(
    Captured = replace_na(Captured, FALSE),
    Status = if_else(Captured, "Captured", "Missed")
  )

# Internal checks: there must be exactly one row per species x network, and no
# captured count can exceed its complete-data denominator.
duplicate_species_network_keys <- species_network_status %>%
  count(
    Study_Network_id,
    Pollinator_accepted_name,
    name = "Rows_per_key"
  ) %>%
  filter(Rows_per_key > 1L)

if (nrow(duplicate_species_network_keys) > 0L) {
  stop("Duplicate species x network keys remain after data preparation.")
}

# -----------------------------------------------------------------------------
# 3. Species-level results
# -----------------------------------------------------------------------------

species_missed <- species_network_status %>%
  group_by(
    Pollinator_order,
    Pollinator_family,
    Pollinator_accepted_name
  ) %>%
  summarise(
    Networks_total = n(),
    Networks_captured = sum(Captured),
    Networks_missed = sum(!Captured),
    Missed_rate = 100 * Networks_missed / Networks_total,
    .groups = "drop"
  ) %>%
  arrange(
    Pollinator_order,
    Pollinator_family,
    desc(Missed_rate),
    desc(Networks_total),
    Pollinator_accepted_name
  )

if (any(species_missed$Networks_captured > species_missed$Networks_total)) {
  stop("Validation failed: a species has more captured than total networks.")
}

bee_species_missed <- species_missed %>%
  filter(Pollinator_family %in% bee_families)

consistently_missed_species <- species_missed %>%
  filter(
    Networks_total >= minimum_networks_for_consistent_miss,
    Missed_rate >= minimum_miss_rate_for_consistent_miss
  )

# -----------------------------------------------------------------------------
# 4. Family-level summaries based on species x network records
# -----------------------------------------------------------------------------

family_missed <- species_network_status %>%
  group_by(Pollinator_order, Pollinator_family) %>%
  summarise(
    Species = n_distinct(Pollinator_accepted_name),
    Networks = n_distinct(Study_Network_id),
    Species_network_records = n(),
    Captured = sum(Status == "Captured"),
    Missed = sum(Status == "Missed"),
    Missed_rate = 100 * Missed / Species_network_records,
    .groups = "drop"
  ) %>%
  arrange(Pollinator_order, desc(Missed_rate), desc(Species_network_records))

# A complementary, species-weighted family summary. Unlike the occurrence-
# weighted rate above, each species contributes equally regardless of how many
# networks it occurred in. This is useful as a sensitivity analysis.
family_species_weighted <- species_missed %>%
  group_by(Pollinator_order, Pollinator_family) %>%
  summarise(
    Species = n(),
    Mean_species_miss_rate = mean(Missed_rate),
    Median_species_miss_rate = median(Missed_rate),
    Species_ever_missed = sum(Networks_missed > 0L),
    Proportion_species_ever_missed = 100 * Species_ever_missed / Species,
    .groups = "drop"
  ) %>%
  arrange(Pollinator_order, desc(Mean_species_miss_rate))

# -----------------------------------------------------------------------------
# 5. Order-level summaries based on species x network records
# -----------------------------------------------------------------------------

order_missed <- species_network_status %>%
  group_by(Pollinator_order) %>%
  summarise(
    Families = n_distinct(Pollinator_family),
    Species = n_distinct(Pollinator_accepted_name),
    Networks = n_distinct(Study_Network_id),
    Species_network_records = n(),
    Captured = sum(Status == "Captured"),
    Missed = sum(Status == "Missed"),
    Missed_rate = 100 * Missed / Species_network_records,
    .groups = "drop"
  ) %>%
  arrange(desc(Missed_rate), desc(Species_network_records))

order_species_weighted <- species_missed %>%
  group_by(Pollinator_order) %>%
  summarise(
    Species = n(),
    Mean_species_miss_rate = mean(Missed_rate),
    Median_species_miss_rate = median(Missed_rate),
    Species_ever_missed = sum(Networks_missed > 0L),
    Proportion_species_ever_missed = 100 * Species_ever_missed / Species,
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_species_miss_rate))


# -----------------------------------------------------------------------------
# Export tables
# -----------------------------------------------------------------------------

write.csv(
  species_network_status,
  file.path(output_dir, "species_network_status.csv"),
  row.names = FALSE
)

write.csv(
  species_missed,
  file.path(output_dir, "species_missed_summary.csv"),
  row.names = FALSE
)

write.csv(
  bee_species_missed,
  file.path(output_dir, "bee_species_missed_summary.csv"),
  row.names = FALSE
)

write.csv(
  consistently_missed_species,
  file.path(output_dir, "species_consistently_missed_sensitivity.csv"),
  row.names = FALSE
)

write.csv(
  family_missed,
  file.path(output_dir, "family_species_network_summary.csv"),
  row.names = FALSE
)

write.csv(
  family_species_weighted,
  file.path(output_dir, "family_species_weighted_summary.csv"),
  row.names = FALSE
)

write.csv(
  order_missed,
  file.path(output_dir, "order_species_network_summary.csv"),
  row.names = FALSE
)

write.csv(
  order_species_weighted,
  file.path(output_dir, "order_species_weighted_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 8-9. Two supplementary figures
# -----------------------------------------------------------------------------
# Layout: Diptera | Hymenoptera; Lepidoptera | Coleoptera.
# Both figures use the same family ordering; no families are removed.
layout_orders <- c("Diptera", "Hymenoptera", "Lepidoptera", "Coleoptera")
family_for_plot <- family_missed %>%
  filter(Pollinator_order %in% layout_orders)
if (nrow(family_for_plot) == 0L) stop("No target orders available for plotting.")

figure_width_cm <- 18
family_spacing_mm <- 3.4
base_font_pt <- 9
family_font_pt <- 9
size_limit <- max(family_for_plot$Species_network_records)
size_breaks <- c(1, 10, 100, 1000)
size_breaks <- size_breaks[size_breaks <= size_limit]

base_theme <- theme_classic(base_size = base_font_pt, base_family = "sans") +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.y = element_text(size = family_font_pt),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 9, margin = margin(t = 5)),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 5)),
    plot.margin = margin(5, 7, 5, 3),
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  )

count_panels <- rate_panels <- vector("list", 4L)
family_counts <- integer(4L)
for (i in seq_along(layout_orders)) {
  order_name <- layout_orders[i]
  df <- family_for_plot %>%
    filter(Pollinator_order == order_name) %>%
    arrange(desc(Species_network_records), Pollinator_family)
  family_counts[i] <- nrow(df)
  panel_title <- paste0("(", letters[i], ") ", order_name)
  if (!nrow(df)) {
    empty <- ggplot() + annotate("text", x = 0, y = 0, label = "No eligible records", size = 3) +
      labs(title = panel_title) + theme_void() +
      theme(plot.title = element_text(size = 11, face = "bold"))
    count_panels[[i]] <- rate_panels[[i]] <- empty
    next
  }
  df <- df %>% mutate(
    Pollinator_family = factor(Pollinator_family, levels = rev(Pollinator_family))
  )
  long_df <- df %>%
    select(Pollinator_family, Captured, Missed) %>%
    pivot_longer(c(Captured, Missed), names_to = "Status", values_to = "Count") %>%
    mutate(Status = factor(Status, levels = c("Captured", "Missed")))

  count_panels[[i]] <- ggplot(long_df, aes(Count, Pollinator_family, fill = Status)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = c(Captured = "#A1D8E8", Missed = "#E76F51"),
                      breaks = c("Captured", "Missed"), drop = FALSE) +
    scale_x_continuous(limits = c(0, max(df$Species_network_records) * 1.04),
                       breaks = scales::breaks_pretty(n = 3),
                       labels = scales::label_number(accuracy = 1, big.mark = ","),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(drop = FALSE, expand = expansion(add = 0.6)) +
    labs(title = panel_title, x = "Species-by-network records", y = NULL, fill = NULL) +
    base_theme

  rate_panels[[i]] <- ggplot(df, aes(Missed_rate, Pollinator_family)) +
    geom_point(aes(size = Species_network_records), colour = "#D95F45", alpha = 0.85) +
    scale_size_continuous(range = c(1.3, 6), limits = c(1, size_limit),
                          breaks = size_breaks, name = "Species-by-network records") +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       labels = scales::label_percent(scale = 1),
                       expand = expansion(add = 4)) +
    scale_y_discrete(drop = FALSE, expand = expansion(add = 0.6)) +
    labs(title = panel_title, x = "Missed records (%)", y = NULL) + base_theme
}

# Each grid row has enough height for its larger family list.
row_heights_mm <- c(max(family_counts[1:2]), max(family_counts[3:4])) *
  family_spacing_mm + 22

extract_legend <- function(p) {
  g <- ggplotGrob(p + theme(legend.position = "bottom"))
  ids <- grep("^guide-box", g$layout$name)
  ids <- ids[!vapply(g$grobs[ids], inherits, logical(1), "zeroGrob")]
  if (!length(ids)) return(grid::nullGrob())
  g$grobs[[ids[1]]]
}

make_2x2 <- function(panels, legend = NULL) {
  grobs <- lapply(panels, function(p) ggplotGrob(p + theme(legend.position = "none")))
  # Align panel edges within columns and across rows.
  for (ids in list(c(1, 3), c(2, 4))) {
    widths <- do.call(grid::unit.pmax, lapply(grobs[ids], function(g) g$widths))
    for (i in ids) grobs[[i]]$widths <- widths
  }
  for (ids in list(c(1, 2), c(3, 4))) {
    heights <- do.call(grid::unit.pmax, lapply(grobs[ids], function(g) g$heights))
    for (i in ids) grobs[[i]]$heights <- heights
  }
  heights_mm <- c(row_heights_mm[1], 4, row_heights_mm[2])
  if (!is.null(legend)) heights_mm <- c(heights_mm, 10)
  out <- gtable::gtable(
    widths = grid::unit.c(grid::unit(1, "null"), grid::unit(5, "mm"), grid::unit(1, "null")),
    heights = grid::unit(heights_mm, "mm")
  )
  for (i in 1:4) {
    out <- gtable::gtable_add_grob(out, grobs[[i]],
      t = if (i <= 2) 1 else 3, l = if (i %% 2 == 1) 1 else 3, clip = "off")
  }
  if (!is.null(legend)) out <- gtable::gtable_add_grob(out, legend, t = 4, l = 1, r = 3)
  out
}

legend_index <- which(family_counts > 0)[1]
p_count <- make_2x2(count_panels, extract_legend(count_panels[[legend_index]]))
p_rate <- make_2x2(rate_panels, extract_legend(rate_panels[[legend_index]]))

save_figure <- function(plot, stem, extra_mm) {
  height_cm <- (sum(row_heights_mm) + 4 + extra_mm) / 10
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot,
         width = figure_width_cm, height = height_cm, units = "cm",
         device = grDevices::cairo_pdf, family = "Arial", bg = "white")
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot,
         width = figure_width_cm, height = height_cm, units = "cm", dpi = 600, bg = "white")
  ggsave(file.path(output_dir, paste0(stem, "_preview.png")), plot,
         width = figure_width_cm, height = height_cm, units = "cm", dpi = 150, bg = "white")
}
save_figure(p_count, "Figure_S6_counts_2x2_v2", 10)
save_figure(p_rate, "Figure_S6.2_miss_rates_2x2_v2", 10)

writeLines(c(
  "Figure S6. Capture and omission of pollinator species-by-network records when monitoring the ten most abundant flowering plant species. Panels show (a) Diptera, (b) Hymenoptera, (c) Lepidoptera and (d) Coleoptera. Each record represents one pollinator species observed in one network. Bars show captured (blue) and missed (orange) records, aggregated by family. Families are ordered by total record count within each order. Count axes differ among panels.",
  "",
  "Figure S6.2. Percentage of pollinator species-by-network records missed when monitoring the ten most abundant flowering plant species. Panels and family ordering match Figure S6. For each family, the percentage is calculated as missed records divided by all species-by-network records, multiplied by 100; it is not the percentage of networks in which the entire family is absent. All panels share a 0-100% scale. Point size indicates total species-by-network record count, using a common size scale across all panels; corresponding captured and missed counts are shown in Figure S6."
), file.path(output_dir, "Figure_S6_S7_captions.txt"))
# Preview in RStudio:
# grid::grid.newpage(); grid::grid.draw(p_count)
# grid::grid.newpage(); grid::grid.draw(p_rate)
