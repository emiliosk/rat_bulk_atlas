# TMM Normalisation

######################################################################
# Script: TMM_normalisation.R
# Author: Emilio Skarwan (emilio.skarwan@scilifelab.se)
# Date: 2026-1-13
# Purpose: Takes TPM normalized data of protein coding genes previously pooled
#         sample data (output of 1_processing_kallisto_output.R), filters out
#         samples that did not pass QC and outputs the TMM factors, and
#         TMM normalized data.
# Input:  ./local/integrated_data/allsample_gene_tpm.tsv,
#         ./local/integrated_data/allsample_gene_var.tsv,
#         ./local/integrated_data/allsample_gene_obs.tsv
# Output: ./local/other_data/TMM_factors_by_ID.csv,
#         ./local/integrated_data/filtered_samples_nTPM.tsv,
#         ./local/integrated_data/filtered_samples_pTPM.tsv,
#         ./local/integrated_data/filtered_samples_obs.tsv,
#         ./local/integrated_data/filtered_samples_var.tsv,
#       figures in ./figures/TMM_normalisation/
#
######################################################################

# ------------------------------------------------------------------------------
# Load Libraries
# ------------------------------------------------------------------------------
devtools::load_all("/Users/emilioskarwan/Documents/SciLifeDrive/AnnDatR")
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(ggplot2)


# ------------------------------------------------------------------------------
# Constants and Configuration
# ------------------------------------------------------------------------------
data_loc <- "./local/integrated_data/"
save_dir <- "./figures/TMM_normalisation/"
ensure_dir(save_dir)

organ_order <- c(
  "Brain",
  "Eye",
  "Endocrine tissues",
  "Respiratory system",
  "Proximal digestive tract",
  "Gastrointestinal tract",
  "Liver & gallbladder",
  "Pancreas",
  "Kidney & urinary bladder",
  "Male reproductive system",
  "Breast & female reproductive system",
  "Muscle & vascular tissue",
  "Connective & soft tissue",
  "Skin",
  "Bone marrow & immune system"
)
# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------
# Open Normalized pseudobulk cluster data
adata <- AnnDatR$new(
  prefix_name = "allsamples_gene",
  layer = "pTPM",
  var_names = "ensembl_id",
  file_dir = data_loc
)
other_meta <- read_tsv("metadata/rat_sampleid_tissue_organ.tsv")
# ------------------------------------------------------------------------------
# Filter out samples which failed QC
# ------------------------------------------------------------------------------

#ignore samples that failed QC starting with "fail" under QC_260114
#samples_to_ignore <- adata$obs %>%  filter(grepl('fail', QC_260114)) %>% pull(ID)
samples_to_ignore <- adata$obs %>% filter(keep_260126 == FALSE) %>% pull(ID)
adata$obs$keep <- !(adata$obs$ID %in% samples_to_ignore)
adata$filter_obs('keep', TRUE)
# Export filtered data
write_tsv(adata$var,
          file.path(data_loc, "filtered_samples_var.tsv")
)
write_tsv(adata$obs, file.path(data_loc, "filtered_samples_obs.tsv"))
write_tsv(adata$X,
          file.path(data_loc, "filtered_samples_pTPM.tsv")
)

X_df <- adata$X %>%
  column_to_rownames(var = adata$var_names_col)

# ------------------------------------------------------------------------------
# Compute TMM Factors
# ------------------------------------------------------------------------------
# Compute TMM factors
tmm_factors <- calc_tmm_normfactors(
  as.matrix(X_df),
  method = "TMM", # checks out
  refColumn = "median", # checks out
  # Trim parameters:
  logratioTrim = 0.3,
  sumTrim = 0.3,
  # Weighting should be done only if count data:
  doWeighting = FALSE # checks out, because weighting variance model is designed on raw counts, we are running on normalised counts
) %>%
  enframe("ID", "tmm_factors")


# ------------------------------------------------------------------------------
# Apply Normalization
# ------------------------------------------------------------------------------
# Apply TMM factors to normalise expression data
expression_data_long <- adata$X %>%
  pivot_longer(cols = -1, names_to = "ID", values_to = "pTPM") %>%
  left_join(tmm_factors, by = "ID") %>%
  mutate(nTPM = pTPM / tmm_factors)


expression_data <- expression_data_long %>%
  select(adata$var_names_col, ID, nTPM) %>%
  pivot_wider(names_from = ID, values_from = nTPM)

#re order columns in expression data
expression_data <- expression_data %>%
  select(adata$var_names_col, adata$obs[["ID"]] %>% unique(), everything())

# ------------------------------------------------------------------------------
# Export Data
# ------------------------------------------------------------------------------
expression_data_long %>%
  select(ID, tmm_factors) %>%
  distinct() %>%
  write_csv("./local/other_data/TMM_factors_by_ID.csv")

expression_data %>%
  write_tsv("./local/integrated_data/filtered_samples_nTPM.tsv")


# ------------------------------------------------------------------------------
# Visualisation
# ------------------------------------------------------------------------------
adata$obs <- adata$obs %>%
  mutate(organ_name = factor(organ_name, levels = organ_order))

tissue_order <- adata$obs %>%
  arrange(
    organ_name,
    #consensus_tissue_name,
    tissue_rank,
    #tissue_group,
   # region_tissue_name,
    tissue_name
  ) %>%
  pull(tissue_name) %>%
  unique()

ID_order <- adata$obs %>%
  arrange(
    organ_name,
    #consensus_tissue_name,
    tissue_rank,
    #tissue_group,
  #  region_tissue_name,
    tissue_name,
    ID
  ) %>%
  select(ID) %>%
  pull(ID)


plot_pal <- get_pal(adata$obs, 'organ_name', 'organ_color')


expression_data_long %>%
  select(ensembl_id, ID, pTPM, nTPM) %>%
  pivot_longer(
    cols = c(-ensembl_id, -ID),
    names_to = "norm_method",
    values_to = "value"
  ) %>%
  mutate(norm_method = factor(norm_method, levels = c("pTPM", "nTPM"))) %>%
  left_join(adata$obs %>% select(ID, tissue_name, organ_name), by = "ID") %>%
  mutate(ID = factor(ID, levels = rev(ID_order))) %>%
  ggplot(aes(x = value, y = ID, fill = organ_name)) +
  geom_boxplot(draw_quantiles = 0.5, outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_x_log10() +
  theme(
    legend.position = "none"
  ) +
  facet_wrap(~norm_method, ncol = 2) +
  scale_fill_manual(values = plot_pal) +
  ggtitle("TMM factors by normalisation method")

ggsave(
  paste0(save_dir, "TMM_factors_by_ID_by_normmethod.pdf"),
  width = 6,
  height = 28
)


adata$obs %>%
  left_join(
    expression_data_long %>% select(ID, tmm_factors) %>% distinct(),
    by = "ID"
  ) %>%
  mutate(tissue_name = factor(tissue_name, levels = rev(tissue_order))) %>%
  ggplot(aes(x = tmm_factors, y = tissue_name, fill = organ_name)) +
  geom_dotplot(dotsize = 0.5, stackdir = 'center', binaxis = 'x') +
  theme(
    legend.position = "none"
  ) +
  scale_fill_manual(values = plot_pal) +
  ggtitle(paste0('TMM factors of all samples (n = ', dim(adata$obs)[1], ')'))

ggsave(
  paste0(save_dir, "TMM_factors_by_tissue_dotplot.pdf"),
  width = 6,
  height = 20
)
