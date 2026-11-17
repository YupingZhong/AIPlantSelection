# Independent audit: never overwrites project input files or existing results.
# Run in RStudio with Source. Requires dplyr, stringr and tidyr.
suppressPackageStartupMessages({library(dplyr); library(stringr)})
project_dir <- 'D:/Chap1_TargetPlant_to_monitor/Project_PlantSelection'
output_dir <- Sys.getenv('TAXONOMY_OUTPUT_DIR', unset = file.path(
  project_dir, 'results_uncertain_identification_sensitivity',
  format(Sys.time(), '%Y%m%d_%H%M%S')))
if (dir.exists(output_dir)) stop('Output directory already exists; rerun later.')
dir.create(output_dir, recursive=TRUE)

read_input <- function(name) readRDS(file.path(project_dir, 'data', name))
full <- read_input('processed/data_interact_published.rds')
top10 <- read_input('processed/selected_plant_result_abun10.rds')
raw <- read_input('raw/Interaction_data_published.rds')

flag_taxonomy <- function(d) {
  required <- c('Pollinator_accepted_name','Pollinator_original_name',
                'Pollinator_rank','Pollinator_unsure_id','Pollinator_uncertainty_type')
  stopifnot(all(required %in% names(d)))
  d %>% mutate(
    accepted = str_squish(as.character(Pollinator_accepted_name)),
    rank_species = !is.na(Pollinator_rank) & toupper(Pollinator_rank)=='SPECIES',
    name_species = !is.na(accepted) & str_count(accepted,'\\S+')>=2L,
    uncertainty_present = !is.na(Pollinator_uncertainty_type) &
      str_trim(as.character(Pollinator_uncertainty_type))!='',
    uncertain = (!is.na(Pollinator_unsure_id) &
      tolower(str_trim(as.character(Pollinator_unsure_id)))=='yes') | uncertainty_present,
    # Conservative operational filter, not independent taxonomic validation.
    strict_species = rank_species & name_species & !uncertain &
      !is.na(Pollinator_unsure_id) &
      tolower(str_trim(as.character(Pollinator_unsure_id)))=='no')
}
full <- flag_taxonomy(full)
raw <- flag_taxonomy(raw)
top10 <- flag_taxonomy(top10)
save_table <- function(x, name) write.csv(as.data.frame(x),
  file.path(output_dir,paste0(name,'.csv')),row.names=FALSE,na='')

# STEP 1: describe resolution BEFORE the extra species-only filter.
# Processed rows are aggregated interaction rows, NOT individual insects.
resolution_table <- function(d) d %>% group_by(Pollinator_order) %>% summarise(
  Interaction_rows=n(), Species_rank_rows=sum(rank_species),
  Coarser_or_unknown_rank_rows=sum(!rank_species),
  Species_rank_percent=100*mean(rank_species),
  Species_rank_uncertain_rows=sum(rank_species & uncertain),
  Strict_species_rows=sum(strict_species), Strict_species_percent=100*mean(strict_species),
  .groups='drop')
processed_resolution <- resolution_table(full)
save_table(processed_resolution,'01_processed_resolution_by_order')
save_table(resolution_table(raw),'02_raw_database_resolution_by_order')
save_table(full %>% count(Pollinator_rank,Pollinator_unsure_id,
                         Pollinator_uncertainty_type,name='Interaction_rows'),
           '03_rank_and_uncertainty')
save_table(full %>% filter(rank_species != name_species | (name_species & uncertain)) %>%
  distinct(Pollinator_original_name,Pollinator_accepted_name,Pollinator_rank,
           Pollinator_unsure_id,Pollinator_uncertainty_type),'04_names_for_review')

# STEP 2: compare the EXISTING miss-rate name filter with strict screening.
# Keep the saved top-10 plant selection fixed; no reselection of plants.
valid_metadata <- function(d) d %>% filter(
  !is.na(Study_Network_id),Study_Network_id!='',
  !is.na(Pollinator_order),Pollinator_order!='',
  !is.na(Pollinator_family),Pollinator_family!='',
  !is.na(accepted),accepted!='')
get_captured_keys <- function(method, positive_only=FALSE) {
  t <- top10
  if (method %in% c('Species_rank_and_name','Original_name_filter')) {
    t <- if (method=='Original_name_filter') filter(t,name_species) else
      filter(t,rank_species & name_species)
  } else if(method=='Strict_species') {
    t <- filter(t,strict_species)
  }
  if(positive_only) t <- filter(t,!is.na(Interaction_addup),Interaction_addup>0)
  t %>% transmute(Study_Network_id,accepted) %>% distinct()
}
miss_status <- function(d,method) {
  captured_keys <- get_captured_keys(method)
  universe <- valid_metadata(d) %>% distinct(Study_Network_id,Pollinator_order,
                                            Pollinator_family,accepted)
  conflicts <- universe %>% count(Study_Network_id,accepted) %>% filter(n>1)
  if(nrow(conflicts)) stop('Conflicting order/family assignments in miss universe.')
  universe %>% left_join(captured_keys %>% mutate(Captured=TRUE),
                         by=c('Study_Network_id','accepted')) %>%
    mutate(Captured=coalesce(Captured,FALSE), Method=method)
}
statuses <- bind_rows(
  miss_status(filter(full,name_species),'Original_name_filter'),
  miss_status(filter(full,rank_species & name_species),'Species_rank_and_name'),
  miss_status(filter(full,strict_species),'Strict_species'))
order_miss <- statuses %>% group_by(Method,Pollinator_order) %>% summarise(
  Taxa=n_distinct(accepted), Networks=n_distinct(Study_Network_id),
  Taxon_network_records=n(), Missed=sum(!Captured),
  Miss_rate_percent=100*mean(!Captured),.groups='drop')
save_table(order_miss,'05_miss_rates_by_order')
save_table(statuses %>% group_by(Method,Pollinator_order,Pollinator_family,accepted) %>%
  summarise(Networks=n(),Missed_networks=sum(!Captured),
            Miss_rate_percent=100*mean(!Captured),.groups='drop'),
  '06_species_miss_rates')

# STEP 3: top-10 richness coverage, with the SAME plant sample and positive links.
# Networks with no eligible pollinator taxa are retained with NA coverage, not zero.
network_ids <- full %>% distinct(Study_Network_id)
richness <- function(d,method) {
  captured_keys <- get_captured_keys(method,positive_only=TRUE)
  universe <- d %>% filter(!is.na(accepted),accepted!='',
                           !is.na(Interaction_addup),Interaction_addup>0) %>%
    distinct(Study_Network_id,accepted) %>%
    left_join(captured_keys %>% mutate(Captured=TRUE),
              by=c('Study_Network_id','accepted')) %>%
    mutate(Captured=coalesce(Captured,FALSE)) %>% group_by(Study_Network_id) %>%
    summarise(Total_taxa=n(),Captured_taxa=sum(Captured),.groups='drop')
  network_ids %>% left_join(universe,by='Study_Network_id') %>% mutate(
    Total_taxa=coalesce(Total_taxa,0L),Captured_taxa=coalesce(Captured_taxa,0L),
    Coverage_percent=if_else(Total_taxa>0,100*Captured_taxa/Total_taxa,NA_real_),
    Method=method)
}
richness_results <- bind_rows(richness(full,'Original_all_accepted_taxa'),
  richness(filter(full,rank_species & name_species),'Species_rank_and_name'),
  richness(filter(full,strict_species),'Strict_species'))
save_table(richness_results,'07_top10_richness_by_network')
richness_summary <- richness_results %>% group_by(Method) %>% summarise(
  Networks_with_eligible_taxa=sum(Total_taxa>0),Mean_total_taxa=mean(Total_taxa),
  Mean_coverage_percent=mean(Coverage_percent,na.rm=TRUE),
  Median_coverage_percent=median(Coverage_percent,na.rm=TRUE),.groups='drop')
save_table(richness_summary,'08_top10_richness_summary')

# Additional diagnostics: strict filtering acts on RECORDS in both full and
# selected datasets. A species may remain if it also has unflagged records.
species_base <- filter(full,rank_species & name_species)
species_strict <- filter(full,strict_species)
excluded_records <- filter(full,rank_species & !strict_species)
save_table(excluded_records,'11_species_rank_records_excluded_in_strict')
save_table(filter(full,rank_species & uncertain) %>%
  group_by(Pollinator_order,Pollinator_original_name,accepted,Pollinator_uncertainty_type) %>%
  summarise(Interaction_rows=n(),Networks=n_distinct(Study_Network_id),
    Interactions=sum(Interaction_addup,na.rm=TRUE),.groups='drop'),
  '12_uncertain_names_and_types')
base_keys <- species_base %>% distinct(Study_Network_id,accepted)
strict_keys <- species_strict %>% distinct(Study_Network_id,accepted)
save_table(anti_join(base_keys,strict_keys,by=c('Study_Network_id','accepted')),
  '13_species_network_occurrences_lost')
save_table(anti_join(distinct(species_base,accepted),distinct(species_strict,accepted),
  by='accepted'),'14_species_no_longer_represented')
comparison <- richness_results %>% select(Study_Network_id,Method,Total_taxa,
  Captured_taxa,Coverage_percent) %>% tidyr::pivot_wider(names_from=Method,
  values_from=c(Total_taxa,Captured_taxa,Coverage_percent)) %>%
  mutate(Delta_coverage_percentage_points=Coverage_percent_Strict_species-
    Coverage_percent_Species_rank_and_name)
save_table(comparison,'15_paired_network_richness_changes')
order_comparison <- order_miss %>% filter(Method!='Original_name_filter') %>%
  select(Method,Pollinator_order,Taxon_network_records,Missed,Miss_rate_percent) %>%
  tidyr::pivot_wider(names_from=Method,
    values_from=c(Taxon_network_records,Missed,Miss_rate_percent)) %>%
  mutate(Delta_miss_percentage_points=Miss_rate_percent_Strict_species-
    Miss_rate_percent_Species_rank_and_name)
save_table(order_comparison,'16_paired_order_miss_changes')

# STEP 4: AI list sensitivity. Mirrors script 3's mean within-network relative
# interaction contribution; study/network IDs intentionally do not split years.
# The denominator is recomputed within each eligible dataset, as in script 3.
priority <- function(d,method) {
  a <- d %>% mutate(Study_Network_id=paste(Study_id,Network_id,sep='_')) %>%
    filter(!is.na(accepted),accepted!='') %>% group_by(Study_Network_id,accepted) %>%
    summarise(Visits=sum(Interaction,na.rm=TRUE),.groups='drop') %>%
    group_by(Study_Network_id) %>% mutate(Total=sum(Visits)) %>% ungroup()
  if(any(a$Total<=0)) stop('Nonpositive network total in priority analysis.')
  nn <- n_distinct(a$Study_Network_id)
  a %>% mutate(Relative=Visits/Total) %>% group_by(accepted) %>%
    summarise(Mean_relative=sum(Relative)/nn,.groups='drop') %>%
    arrange(desc(Mean_relative),accepted) %>% mutate(Rank=row_number(),
      Cumulative_percent=100*cumsum(Mean_relative)/sum(Mean_relative),Method=method)
}
priority_results <- bind_rows(priority(filter(raw,rank_species),'Original_species_rank'),
                              priority(filter(raw,strict_species),'Strict_species'))
priority_summary <- priority_results %>% group_by(Method) %>% summarise(
  Total_taxa=n(),Taxa_to_95_percent=min(Rank[Cumulative_percent>=95]),.groups='drop')
save_table(priority_results,'09_AI_priority_rankings')
save_table(priority_summary,'10_AI_priority_summary')
saveRDS(list(order_miss=order_miss,richness=richness_results,
             priority=priority_results),file.path(output_dir,'audit_results.rds'))
writeLines(c('This audit does not change or overwrite project inputs.',
  'Strict means species-rank, multiword accepted name, explicit No uncertainty flag, and no uncertainty type.',
  'Strict filtering is not independent verification of identification accuracy.',
  'Both full and selected top-10 records are filtered separately; uncertain captured records cannot count as strict captures.',
  'The top-10 plant sample is held fixed. Other sampling strategies and structural metrics are NOT rerun.',
  'Accepted morphospecies labels may have been collapsed upstream; inspect original names.',
  'Priority-list calculations match the raw-data, mean-relative-contribution approach in script 3.',
  'A 95% cutoff refers to cumulative mean within-network relative contribution, not necessarily pooled visits.'),
  file.path(output_dir,'README.txt'))
cat('\nSTEP 1: resolution\n');print(processed_resolution)
cat('\nSTEP 2: miss rates\n');print(order_miss %>% filter(Pollinator_order %in%
  c('Diptera','Hymenoptera','Coleoptera','Lepidoptera')),n=20)
cat('\nSTEP 3: top-10 richness\n');print(richness_summary)
cat('\nSTEP 4: priority list\n');print(priority_summary)
cat('\nResults saved to:\n',output_dir,'\n')

# Compare list membership at each method's own 95% cutoff.
priority_membership <- priority_results %>% group_by(Method) %>%
  mutate(In_priority_list=Rank<=min(Rank[Cumulative_percent>=95])) %>% ungroup() %>%
  dplyr::select(accepted,Method,Rank,Mean_relative,In_priority_list) %>%
  tidyr::pivot_wider(names_from=Method,
    values_from=c(Rank,Mean_relative,In_priority_list)) %>%
  mutate(In_priority_list_Original_species_rank=coalesce(In_priority_list_Original_species_rank,FALSE),
         In_priority_list_Strict_species=coalesce(In_priority_list_Strict_species,FALSE))
save_table(priority_membership,'17_AI_priority_membership_comparison')
input_paths <- file.path(project_dir,'data',c('processed/data_interact_published.rds',
 'processed/selected_plant_result_abun10.rds','raw/Interaction_data_published.rds'))
save_table(data.frame(File=input_paths,MD5=unname(tools::md5sum(input_paths)),
 Bytes=file.info(input_paths)$size),'18_input_file_fingerprints')
writeLines(capture.output(sessionInfo()),file.path(output_dir,'R_session_info.txt'))
cat('\nSpecies-rank rows:',sum(full$rank_species),
    '| flagged uncertain:',sum(full$rank_species & full$uncertain),
    '| strict retained:',sum(full$strict_species),'\n')
print(order_comparison,width=Inf)
cat('\nAll checks completed. No project inputs overwritten.\n')

