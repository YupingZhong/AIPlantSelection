# AIPlantSelection

R code supporting the manuscript **Optimizing plant species selection for automated monitoring of plant–pollinator interactions**.

This project evaluates plant-selection strategies for pollinator monitoring using EuPPollNet interaction and floral-abundance data. Analyses compare abundance-based, flower-shape abundance-based, phylogenetically informed, and random selection; assess pollinator richness and network structure; identify rare but highly attractive plants and potentially missed pollinators; and rank pollinator species for AI classifier development.

## Data

Source interaction data are available from EuPPollNet: https://doi.org/10.5281/zenodo.15183272.

The plant-selection analyses use a filtered dataset of 267 networks from 27 studies in 13 European countries. The continental AI-priority analysis uses the broader interaction database, rather than only this filtered subset. See the manuscript and scripts for filtering and aggregation rules.

The data-preparation script expects:

```text
data/raw/Interaction_data_published.rds
data/raw/Flower_counts_published.rds
data/raw/species_checked_wof.csv
data/raw/plant_sampling_unit.csv
data/processed/merge.trait.csv
```

The auxiliary CSV files and trait table are project-specific inputs and are not necessarily contained in the source database deposit. They must be available at the paths above before running the workflow. Downloading the source database alone does not necessarily supply all project inputs.

## Software

The analyses were conducted in R 4.4.0. RStudio can be used by opening `Project_PlantSelection.Rproj`.

Package requirements are declared by `library()` calls in the individual scripts. Major analysis packages include `bipartite`, `vegan`, `lme4`, `lmerTest`, `emmeans`, `betapart`, `permute`, `ape`, and `V.PhyloMaker2`. Data processing and plotting additionally use packages such as `dplyr`, `tidyr`, `stringr`, `purrr`, `lubridate`, `ggplot2`, and `cowplot`. Consult each script for its complete requirements; phylogenetic and trait-extraction packages may require separate installation instructions.

To record the environment used for a run:

```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

## Project structure

```text
Project_PlantSelection.Rproj
scripts/          # Data processing, analyses, and figure generation
data/raw/         # Source and auxiliary input files
data/processed/   # Derived datasets, traits, and cached analysis objects
results/          # Analysis outputs and figures
```

Some scripts use additional result directories specified within the script.

## Running the analyses

Work from a local copy of the repository. Set the working directory to the project root, or open the RStudio project. Relative input paths are resolved from that directory.

Start with data preparation and the plant-subsampling calculations:

```r
source("scripts/0.1.data.preparation.R")
source("scripts/1.1 MainResults -richness captured by subsamplings(zero_capture_fixed).R")
```

Then run the modules needed for the desired results. The following is a guide to their roles, not a single unattended master pipeline:

| Module | Scripts |
| --- | --- |
| Geographic coverage and overview | `0.2 network_locations_bioregion.R`; `0.3. networks.overview.R` |
| Network metrics | `1.2.connectance.R`; `1.2.nestedness.R`; `1.2. H2.R`; `1.2.merge.and.random.metrics.R` |
| Richness comparisons and Holm-adjusted figure | `1.3 MainResults - Richness figure output_Holm.R` |
| Progressive and proportional plant sampling | Scripts beginning `1.4` and `1.5` |
| Rare-attractive plants and sampling effectiveness | Scripts beginning `2.1` and `2.2` |
| Flower-shape composition and dispersion analyses | `2.3 FollowingAnalysis - Flower_shape_analysis_WITH_betadisper.R` |
| Continental AI-priority ranking | `3. Species for AI indentification.R` |
| Miss rates and taxonomic summaries | `Supp_Pollinator_missed_analysis_in_order_revised.R` |
| Identification-uncertainty sensitivity | `taxonomy_check.R` |
| Honeybee-exclusion sensitivity | `Supp - Apis_mellifera_sensitivity.R` |
| Other supplementary analyses | Remaining scripts beginning `Supp` |

Read each script's input section before execution. Several scripts depend on processed objects created by earlier modules. Some scripts may contain machine-specific paths or cache switches; update these to the local project location and rerun cached analyses when inputs change. Running scripts may replace derived files, so retain a separate copy of inputs and results before rerunning.

Random subsampling and permutation analyses are stochastic. Preserve any seeds set in the scripts and record the R/package environment when reproducing results.

## Interpretation

- Floral abundance is used for ranking plant species within networks; raw abundance values are not compared across studies using different measurement units.
- Species-level rank and identification uncertainty are separate database fields. Consult the taxonomic-sensitivity script for its additional screening rules.
- The AI-priority threshold concerns cumulative mean within-network relative interaction contributions, not a percentage of the total pollinator species pool or necessarily pooled interaction counts.

## Citation and code archive

GitHub repository: https://github.com/YupingZhong/AIPlantSelection

For the archived code DOI, use the Zenodo record associated with the corresponding GitHub release once that release has been archived. The source-data DOI above identifies the input database, not this code archive.

## Contact

For questions or reproducibility issues, open an issue in the GitHub repository and include the script name, error message, and output of `sessionInfo()`.
