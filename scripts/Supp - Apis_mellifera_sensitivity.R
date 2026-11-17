# R1.6: sensitivity of the classifier-priority analysis to Apis mellifera.
# Run from the original project root:
# source("path/to/Apis_mellifera_sensitivity.R")
# Or: Rscript Apis_mellifera_sensitivity.R input.rds output_directory
# Only base R is required for calculations; ggplot2 is optional for PNG plots.

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[1] else "data/raw/Interaction_data_published.rds"
output_dir <- if (length(args) >= 2L) args[2] else "results/Apis_mellifera_sensitivity"
target_pct <- 95
# Do not automatically exclude Bombus terrestris: species identity alone does
# not distinguish managed from wild individuals.
excluded_species <- "Apis mellifera"

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
metadata <- readRDS(input_file)
required <- c("Study_id", "Network_id", "Pollinator_accepted_name",
              "Pollinator_rank", "Interaction")
if (!all(required %in% names(metadata))) {
  stop("Missing columns: ", paste(setdiff(required, names(metadata)), collapse = ", "))
}
if (!is.numeric(metadata$Interaction)) stop("Interaction must be numeric.")
if (any(!is.finite(metadata$Interaction)) || any(metadata$Interaction < 0)) {
  stop("Interaction contains missing, nonfinite or negative values; resolve before analysis.")
}
if (anyNA(metadata$Study_id) || anyNA(metadata$Network_id)) stop("Missing network identifiers.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
metadata$Study_Network_id <- paste(metadata$Study_id, metadata$Network_id, sep = "_")
# Same species-rank filter as the supplied script. Positive interactions define
# species presence; zero-only records cannot establish occurrence.
d <- metadata[!is.na(metadata$Pollinator_accepted_name) &
                !is.na(metadata$Pollinator_rank) &
                metadata$Pollinator_rank == "SPECIES", , drop = FALSE]
d$Pollinator_accepted_name <- trimws(gsub("[[:space:]]+", " ", d$Pollinator_accepted_name))
d <- d[nzchar(d$Pollinator_accepted_name) & d$Interaction > 0, , drop = FALSE]
if (!nrow(d)) stop("No positive species-level interaction records.")
if (!any(d$Pollinator_accepted_name %in% excluded_species)) {
  stop("Apis mellifera not found under its exact accepted name. Check taxonomic names.")
}

write_out <- function(x, name) {
  write.csv(x, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE, na = "")
}

# Zero-filling is implicit: sum(network-relative interactions) / N networks
# equals averaging after completing every species x network combination with 0.
analyse <- function(x, scenario) {
  totals <- aggregate(x$Interaction, list(Network = x$Study_Network_id), sum)
  names(totals)[2] <- "Total_interactions"
  counts <- aggregate(x$Interaction,
                      list(Species = x$Pollinator_accepted_name,
                           Network = x$Study_Network_id), sum)
  names(counts)[3] <- "Interactions"
  counts$Relative <- counts$Interactions / totals$Total_interactions[
    match(counts$Network, totals$Network)]
  ranks <- aggregate(counts$Relative, list(Species = counts$Species), sum)
  names(ranks)[2] <- "Mean_relative_interactions"
  ranks$Mean_relative_interactions <- ranks$Mean_relative_interactions / nrow(totals)
  present <- table(counts$Species)
  ranks$Networks_present <- as.integer(present[ranks$Species])
  # Alphabetical tie-break makes rankings reproducible.
  ranks <- ranks[order(-ranks$Mean_relative_interactions, ranks$Species), ]
  rownames(ranks) <- NULL
  ranks$Rank <- seq_len(nrow(ranks))
  ranks$Cumulative_pct <- 100 * cumsum(ranks$Mean_relative_interactions) /
    sum(ranks$Mean_relative_interactions)
  n95 <- which(ranks$Cumulative_pct >= target_pct)[1]
  ranks$Priority95 <- ranks$Rank <= n95
  ranks$Scenario <- scenario
  stopifnot(abs(sum(ranks$Mean_relative_interactions) - 1) < 1e-8)
  list(ranking = ranks, network_totals = totals,
       summary = data.frame(Scenario = scenario, Networks = nrow(totals),
                            Species = nrow(ranks), Interactions = sum(x$Interaction),
                            Species_for_95pct = n95,
                            Percent_species_for_95pct = 100 * n95 / nrow(ranks),
                            Coverage_at_cutoff = ranks$Cumulative_pct[n95]))
}

remaining <- d[!d$Pollinator_accepted_name %in% excluded_species, , drop = FALSE]
if (!nrow(remaining)) stop("No interactions remain after exclusion.")
all_networks <- sort(unique(d$Study_Network_id))
retained_networks <- sort(unique(remaining$Study_Network_id))
lost_networks <- setdiff(all_networks, retained_networks)
write_out(data.frame(Network = lost_networks), "networks_empty_after_exclusion")

# Primary comparison: same paired network set in both scenarios. Networks with
# only honeybee interactions have undefined non-honeybee relative proportions;
# they are reported and excluded from BOTH sides of this paired comparison.
paired_full <- d[d$Study_Network_id %in% retained_networks, , drop = FALSE]
full <- analyse(paired_full, "Including Apis mellifera")
without <- analyse(remaining, "Excluding Apis mellifera")
original <- analyse(d, "Original full dataset")
write_out(original$ranking, "original_full_ranking")
write_out(original$summary, "original_full_summary")
write_out(full$ranking, "paired_including_Apis_ranking")
write_out(without$ranking, "paired_excluding_Apis_ranking")
summary_table <- rbind(full$summary, without$summary)
write_out(summary_table, "paired_scenario_summary")

# Honeybee contribution, both pooled-count and equal-network weighted.
apis_network <- aggregate(d$Interaction[d$Pollinator_accepted_name %in% excluded_species],
                          list(Network = d$Study_Network_id[d$Pollinator_accepted_name %in% excluded_species]), sum)
names(apis_network)[2] <- "Apis_interactions"
network_audit <- merge(original$network_totals, apis_network, by = "Network", all.x = TRUE)
network_audit$Apis_interactions[is.na(network_audit$Apis_interactions)] <- 0
network_audit$Apis_pct <- 100 * network_audit$Apis_interactions / network_audit$Total_interactions
network_audit$Retained_in_paired_analysis <- network_audit$Network %in% retained_networks
write_out(network_audit, "Apis_contribution_by_network")
write_out(data.frame(Networks = length(all_networks),
                     Networks_with_Apis = sum(network_audit$Apis_interactions > 0),
                     Networks_empty_after_exclusion = length(lost_networks),
                     Pooled_Apis_interaction_pct = 100 * sum(network_audit$Apis_interactions) /
                       sum(network_audit$Total_interactions),
                     Mean_network_Apis_interaction_pct = mean(network_audit$Apis_pct)),
          "Apis_contribution_summary")

comparison <- merge(full$ranking[, c("Species", "Rank", "Mean_relative_interactions", "Priority95")],
                    without$ranking[, c("Species", "Rank", "Mean_relative_interactions", "Priority95")],
                    by = "Species", all = TRUE, suffixes = c("_including", "_excluding"))
comparison$Rank_change <- comparison$Rank_excluding - comparison$Rank_including
write_out(comparison, "species_rank_and_priority_comparison")
# Compare non-honeybee lists only; distinguish reranking from simply dropping Apis.
a <- full$ranking$Species[full$ranking$Priority95 & !full$ranking$Species %in% excluded_species]
b <- without$ranking$Species[without$ranking$Priority95]
shared <- intersect(a, b)
common <- comparison[complete.cases(comparison[, c("Rank_including", "Rank_excluding")]), ]
write_out(data.frame(Non_Apis_priority_including = length(a),
                     Priority_excluding = length(b), Shared_priority = length(shared),
                     Jaccard = length(shared) / length(union(a, b)),
                     Species_entering = length(setdiff(b, a)), Species_leaving = length(setdiff(a, b)),
                     Spearman_rank_correlation = if (nrow(common) > 1L)
                       cor(common$Rank_including, common$Rank_excluding, method = "spearman") else NA_real_,
                     Original_non_Apis_priority_coverage_after_exclusion_pct =
                       100 * sum(without$ranking$Mean_relative_interactions[without$ranking$Species %in% a])),
          "priority_list_comparison_summary")
write_out(data.frame(Species = setdiff(b, a)), "species_entering_priority_list")
write_out(data.frame(Species = setdiff(a, b)), "species_leaving_priority_list")

curves <- rbind(full$ranking, without$ranking)
curves$Percent_species <- 100 * curves$Rank /
  summary_table$Species[match(curves$Scenario, summary_table$Scenario)]
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  plot_curve <- function(x, title) {
    ggplot(curves, aes(x = .data[[x]], y = Cumulative_pct, colour = Scenario)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = target_pct, linetype = "dashed", colour = "grey50") +
      scale_colour_manual(values = c("Including Apis mellifera" = "#D49A20",
                                     "Excluding Apis mellifera" = "#3876B5")) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
      labs(x = title, y = "Cumulative mean relative interactions (%)", colour = NULL) +
      theme_classic(base_size = 12) + theme(legend.position = "bottom")
  }
  ggsave(file.path(output_dir, "Figure5_Apis_sensitivity.png"),
         plot_curve("Rank", "Pollinator species identified"), width = 7, height = 5, dpi = 300)
  ggsave(file.path(output_dir, "Figure5_Apis_sensitivity_percent_species.png"),
         plot_curve("Percent_species", "Pollinator species identified (% of scenario total)"),
         width = 7, height = 5, dpi = 300)
} else message("ggplot2 unavailable: CSV results saved; PNG plots skipped.")

notes <- c(
  "R1.6 sensitivity analysis: Figure 5 / classifier-priority list only.",
  "Within each scenario, recompute network denominators, then average relative interactions equally across networks, including absence zeros.",
  "The cumulative curve represents mean within-network interaction proportions, NOT the proportion of pooled raw interaction counts.",
  "Paired analysis excludes from BOTH scenarios networks with no remaining non-Apis interactions; see audit CSV.",
  "Original full-dataset ranking is also exported to check reproduction of the manuscript baseline.",
  "Names are whitespace-normalized and matched exactly; synonyms need prior taxonomic harmonization.",
  "No other managed species are automatically excluded, because management cannot be inferred from species identity.",
  "Original by-order code averages only networks where a species occurs; it differs from the main zero-filled calculation and is not reproduced here.",
  "This script does not rerun plant-selection, richness or network-structure analyses from other manuscript sections.",
  "Threshold-list membership can depend on tied ranks; this script uses alphabetical tie-breaking.")
writeLines(notes, file.path(output_dir, "analysis_notes.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
print(summary_table)
message("Results saved to: ", normalizePath(output_dir))
