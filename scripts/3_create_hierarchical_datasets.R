# create hierarchical datasets

######################################################################
# Script: create_hierarchical_datasets.R
# Author: Emilio Skarwan (emilio.skarwan@scilifelab.se)
# Date: 2026-1-16
# Purpose: Takes sample level info and creates the coarser hierarchical datasets.
# Input:  ./local/integrated_data/filtered_samples_nTPM.tsv,
#         ./local/integrated_data/filtered_samples_var.tsv,
#         ./local/integrated_data/filtered_samples_obs.tsv
# Output:
#         ./local/integrated_data/samples_obs.tsv
#         ./local/integrated_data/samples_var.tsv,
#         ./local/integrated_data/samples_nTPM.tsv,
#
######################################################################

devtools::load_all("/Users/emilioskarwan/Documents/SciLifeDrive/AnnDatR")
library(dplyr)
library(tibble)
library(readr)
library(tidyr)


#library(patchwork)


data_loc <- "./local/integrated_data/"

# load in data

# Open Normalized pseudobulk cluster data
adata <- AnnDatR$new(
  prefix_name = "filtered_samples",
  layer = "nTPM",
  var_names = "ensembl_id",
  file_dir = data_loc
)


obs_tissue_name <- adata$obs %>% select(tissue_name,
                                        #region_tissue_name,
                                        consensus_tissue_name,
                                        #tissue_group,
                                        organ_name, organ_color,
                                        #region_color
)  %>% distinct()

X_tissue <- adata$X %>%
  pivot_longer(
    cols = -ensembl_id,
    names_to = "ID",
    values_to = "nTPM"
  ) %>%
  left_join(adata$obs %>% select(ID, tissue_name)) %>%
  group_by(ensembl_id, tissue_name) %>%
  summarise(nTPM = mean(nTPM)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = tissue_name,
    values_from = nTPM,
    values_fill = 0
  ) %>% select(ensembl_id,obs_tissue_name$tissue_name, everything())


obs_tissue_name %>% write_tsv(file.path(data_loc, "tissue_obs.tsv"))
X_tissue %>% write_tsv(file.path(data_loc, "tissue_nTPM.tsv"))
adata$var %>% write_tsv(file.path(data_loc, "tissue_var.tsv"))




obs_consensus_name <- adata$obs %>% select(consensus_tissue_name,
                                           organ_name, organ_color)  %>% distinct()

X_consensus <- X_tissue %>%
  pivot_longer(
    cols = -ensembl_id,
    names_to = "tissue_name",
    values_to = "nTPM"
  ) %>%
  left_join(adata$obs %>% select(tissue_name, consensus_tissue_name)) %>%
  group_by(ensembl_id, consensus_tissue_name) %>%
  summarise(nTPM = max(nTPM)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = consensus_tissue_name,
    values_from = nTPM,
    values_fill = 0
  ) %>% select(ensembl_id,obs_tissue_name$consensus_tissue_name, everything())


obs_consensus_name %>% write_tsv(file.path(data_loc, "consensus_obs.tsv"))
X_consensus %>% write_tsv(file.path(data_loc, "consensus_nTPM.tsv"))
adata$var %>% write_tsv(file.path(data_loc, "consensus_var.tsv"))




