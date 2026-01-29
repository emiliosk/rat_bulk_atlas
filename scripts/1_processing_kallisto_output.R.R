# Processing kallisto output

######################################################################
# Script: processing_kallisto_output.R
# Author: Emilio Skarwan (emilio.skarwan@scilifelab.se)
# Date: 2026-1-13
# Purpose: Pools kallisto output files into gene and transcript level.
#
# Input:  kallisto output files under ./local/kallisto_out/
#         sample metadata ./metadata/mouse_sample_meta_20251212.tsv
#         metadata ensemble transcript ids ./metadata/ensembl_gene_transcript_mouse.tsv
# Output: ./local/integrated_data/allsamples_gene_pTPM.tsv,
#         ./local/integrated_data/allsamples_gene_var.tsv,
#         ./local/integrated_data/allsamples_gene_obs.tsv
#         ./local/integrated_data/sample_transcript_TPM.tsv,
#         ./local/integrated_data/sample_transcript_var.tsv,
#         ./local/integrated_data/sample_transcript_obs.tsv
######################################################################
# ------------------------------------------------------------------------------
# Load Libraries
# ------------------------------------------------------------------------------
library(tidyr)
library(readr)
library(dplyr)
library(purrr)
library(jsonlite)


# ------------------------------------------------------------------------------
# Constants and Configuration
# ------------------------------------------------------------------------------
SPECIES <- "rnorvegicus"
VERSION <- "109"
TRANSCRIPTOME_SIZE <- 	2647915728
COVERAGE_BASES <- 200
data_dir <- "./local/kallisto_out/"
count_export_dir <- "./local/integrated_data/"

sample_metadata_path <- "./metadata/rat_sampleid_tissue_organ_updated_colors.tsv"
pcg_path <- "./metadata/ensembl_gene_transcript_rat.tsv"
biotypes_to_include <- c(
  "protein_coding",
  "protein_coding_LoF",
  "IG_V_gene",
  "IG_D_gene",
  "IG_C_gene",
  "IG_J_gene",
  "TR_V_gene",
  "TR_J_gene",
  "TR_D_gene",
  "TR_C_gene"
)



# ------------------------------------------------------------------------------
# Load Metadata
# ------------------------------------------------------------------------------
# Read in Metadata files
sample_metadata <- read_tsv(sample_metadata_path) %>%
  rename(sample_id = scilifelab_id,
         tissue_name = from_sample_tissue,
         organ_rank = organ_order,
         tissue_rank = tissue_order)
sample_codes <- sample_metadata$sample_id

sample_metadata <- sample_metadata %>%
  mutate(ID = paste0(gsub(" ", "-", tissue_name), "_", sample_id)) %>%
  select(ID, everything())

sample_codes <- setNames(sample_codes, sample_codes)

# run_infos <- sample_codes %>%
#   map(\(sample) fromJSON(file.path(data_dir, sample, "run_info.json")))
# run_infos <- tibble(sample_id = names(run_infos), repo = run_infos) %>%
#   tidyr::unnest_wider(repo) # %>% select(-`k-mer length`)
# run_infos <- run_infos %>%
#   mutate(coverage = (n_pseudoaligned * COVERAGE_BASES) / TRANSCRIPTOME_SIZE)

# sample_metadata <- sample_metadata %>% left_join(run_infos, by = "sample_id")

hpa_pcgs <- read_delim(
  pcg_path,
  delim = " ",
  col_names = c("transcript_id", "gene_id", "gene_type", "transcript_type")
)
hpa_pcgs <- hpa_pcgs %>%
  filter(transcript_type %in% biotypes_to_include) %>%
  filter(gene_type %in% biotypes_to_include)


# ------------------------------------------------------------------------------
# Load Expression Data
# ------------------------------------------------------------------------------
# Function to read and process a single file
#  <- function(id) {
#   file_path <- file.path(data_dir, id, "abundance.tsv")
#   read_tsv(file_path, show_col_types = FALSE) %>%
#     mutate(sample_id = id) %>%
#     rename(ensembl_transcript_id_version = target_id)
# }
# # Read all samples using map and bind_rows
# sample_counts <- map_df(sample_codes, read_sample)
sample_counts <- read_tsv("./local/raw_out/hpa_rat_all_samples_all_transcripts_tpm_109.tsv") %>%
  rename(sample_id = scilifelab_id,
         ensembl_transcript_id = enst_id)
# split sample_comment into tissue_name and sex, where there is a '|'
sample_counts <- sample_counts %>% mutate(
  tissue_name = sub("\\|.*", "", sample_comment),
  sex = sub(".*\\|", "", sample_comment)) %>%
  select(-sample_comment)


# Extract `adata.var`
adata.var <- sample_counts %>% select(1) %>% distinct()
sample_sex_map <- sample_counts %>% select(sample_id, sex) %>%
  distinct()
# Pivot to wide format
adata.transcript_tpm <- sample_counts %>%
  select(ensembl_transcript_id, sample_id, tpm) %>%
  left_join(sample_metadata %>% select(ID, sample_id)) %>%
  select(-sample_id) %>%
  pivot_wider(names_from = ID, values_from = tpm)


# ------------------------------------------------------------------------------
# Fetch Ensembl Annotations
# ------------------------------------------------------------------------------

ensembl <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = paste0(SPECIES, "_gene_ensembl"),
  version = VERSION
)
gene_metadata <- biomaRt::getBM(
  attributes = c(
    "ensembl_transcript_id_version",
    "ensembl_transcript_id",
    "ensembl_gene_id",
    "external_gene_name",
    "transcript_length",
    "chromosome_name"
  ),
  mart = ensembl
)
gene_metadata <- gene_metadata %>%
  mutate(pcg = ensembl_transcript_id %in% hpa_pcgs$transcript_id)

adata.var <- adata.var %>%
  left_join(gene_metadata %>% select(-ensembl_transcript_id_version) %>% distinct(),
            by = "ensembl_transcript_id")

get_hs_homologs <- function(mart_conn) {
  biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "hsapiens_homolog_ensembl_gene",
      "hsapiens_homolog_associated_gene_name",
      "hsapiens_homolog_orthology_type"
    ),
    mart = mart_conn
  ) %>%
    mutate(
      gene_symbol = if_else(
        external_gene_name == "",
        ensembl_gene_id,
        external_gene_name
      )
    )
}


hs_homologs <- get_hs_homologs(ensembl)


# ------------------------------------------------------------------------------
# Process Gene-Level Data
# ------------------------------------------------------------------------------
adata.obs <- sample_metadata %>%
  left_join(
    sample_counts %>%
      group_by(sample_id) %>%
      summarise(total_counts_raw = sum(est_count)),
    by = "sample_id"
  )

adata.obs <- adata.obs %>% left_join(sample_sex_map)

adata.gene_ptpm <- adata.transcript_tpm %>%
  left_join(
    adata.var %>%
      select(ensembl_transcript_id, ensembl_gene_id, pcg)
  ) %>%
  #1: filter protein coding transcripts
  filter(pcg) %>%
  select(-pcg) %>%
  pivot_longer(
    cols = c(-ensembl_transcript_id, -ensembl_gene_id),
    names_to = "sample",
    values_to = "tpm"
  ) %>%
  group_by(sample, ensembl_gene_id) %>%
  #2: summarise the total transcripts by summing total in a gene
  summarise(tpm = sum(tpm)) %>%
  ungroup() %>%
  group_by(sample) %>%
  #3: normalize to create pTPM
  mutate(ptpm = tpm * (1e6 / sum(tpm))) %>%
  ungroup() %>%
  select(-tpm) %>%
  pivot_wider(names_from = sample, values_from = ptpm)

adata.obs <- adata.obs %>%
  left_join(
    adata.gene_ptpm %>%
      mutate_if(is.numeric, \(x) x >= 1) %>%
      pivot_longer(cols = -1, names_to = "ID", values_to = "detected") %>%
      group_by(ID) %>%
      summarise("genes_detected" = sum(detected)),
    by = "ID"
  ) %>%
  select(ID, sample_id, everything())

adata.obs <- adata.obs %>% rename(consensus_tissue_name = tissue)

adata.var <- adata.var %>%
  rename(gene_symbol = external_gene_name) %>%
  mutate(
    gene_symbol = case_when(
      gene_symbol == "" ~ ensembl_gene_id,
      .default = gene_symbol
    )
  )


# ------------------------------------------------------------------------------
# Export Data
# ------------------------------------------------------------------------------

# export transcript level data
write_tsv(
  adata.var %>% arrange(ensembl_transcript_id),
  file.path(count_export_dir, "sample_transcript_var.tsv")
)
write_tsv(adata.obs, file.path(count_export_dir, "sample_transcript_obs.tsv"))
write_tsv(
  adata.transcript_tpm %>%
    select(ensembl_transcript_id, adata.obs$ID, everything()) %>%
    arrange(ensembl_transcript_id),
  file.path(count_export_dir, "sample_transcript_TPM.tsv")
)


# export gene level data
adata.var <- adata.var %>%
  filter(pcg) %>%
  select(ensembl_id = ensembl_gene_id, gene_symbol, chromosome_name) %>%
  distinct()

adata.var <- adata.var %>%
  left_join(
    hs_homologs %>%
      group_by(ensembl_id = ensembl_gene_id) %>%
      summarise(
        hsapiens_ensembl_id = paste(
          hsapiens_homolog_ensembl_gene,
          collapse = ", "
        ),
        hsapiens_gene_symbol = paste(
          hsapiens_homolog_associated_gene_name,
          collapse = ", "
        ),
        hsapiens_homolog_orthology_type = paste(
          unique(hsapiens_homolog_orthology_type),
          collapse = ", "
        )
      )
  )

write_tsv(
  adata.var %>% arrange(ensembl_id),
  file.path(count_export_dir, "allsamples_gene_var.tsv")
)
write_tsv(adata.obs, file.path(count_export_dir, "allsamples_gene_obs.tsv"))
write_tsv(
  adata.gene_ptpm %>%
    select(ensembl_id = ensembl_gene_id, adata.obs$ID, everything()) %>%
    arrange(ensembl_id),
  file.path(count_export_dir, "allsamples_gene_pTPM.tsv")
)

