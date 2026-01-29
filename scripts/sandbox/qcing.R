# qc

devtools::load_all("/Users/emilioskarwan/Documents/SciLifeDrive/AnnDatR")
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(ggplot2)
library(purrr)
library(ggridges)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(ggplotify)

data_loc <- "./local/integrated_data/"


organ_order <- c(
  'Brain',
  'Eye',
  'Endocrine tissues',
  'Respiratory system',
  'Proximal digestive tract',
  'Gastrointestinal tract',
  'Liver & gallbladder',
  'Pancreas',
  'Kidney & urinary bladder',
  'Male reproductive system',
  'Breast & female reproductive system',
  'Muscle & vascular tissue',
  'Connective & soft tissue',
  'Skin',
  'Bone marrow & immune system'
)


#

# load in data

# Open Normalized pseudobulk cluster data
adata <- AnnDatR$new(
  prefix_name = "filtered_samples",
  layer = "nTPM",
  var_names = "ensembl_id",
  file_dir = data_loc
)
adata_tpm <- AnnDatR$new(
  prefix_name = "filtered_samples",
  layer = "pTPM",
  var_names = "ensembl_id",
  file_dir = data_loc
)
# this just to make the checkup easier
#adata$obs <- adata$obs %>% mutate(ID = paste0(ID, '_',as.character(as.numeric(include))))
#adata$X <- adata$X %>% rename_with(~ paste0(.x, '_',as.character(as.numeric(adata$obs$include))), -ensembl_id)
#adata_tpm$obs <- adata_tpm$obs %>% mutate(ID = paste0(ID, '_',as.character(as.numeric(include))))
#adata_tpm$X <- adata_tpm$X %>% rename_with(~ paste0(.x, '_',as.character(as.numeric(adata_tpm$obs$include))), -ensembl_id)

dardel_meta <- read_tsv("./metadata/rat_dardel_filename_metadata_20260126.tsv")

rat_map <- dardel_meta %>% select(sample_id = scilifeid, id = sampleid) %>% unique() %>%
  mutate(rat = stringr::str_extract(id, "[^.]+$")) %>% select(-id)


adata$obs <- adata$obs %>% left_join(rat_map)

adata.gene_tpm <- adata_tpm$X

adata$obs <- adata$obs %>%
  mutate(organ_name = factor(organ_name, levels = organ_order))


tissue_order <- adata$obs %>%
  arrange(
    organ_name,
    # consensus_tissue_name,
    tissue_rank,
    #  tissue_group,
    #  region_tissue_name,
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
    # region_tissue_name,
    tissue_name,
    ID
  ) %>%
  select(ID) %>%
  pull(ID)

mt_genes <- adata$var %>%
  filter(chromosome_name == 'MT') %>%
  pull(ensembl_id) %>%
  unique()
adata$obs$total_mt_counts <- adata$X %>%
  filter(ensembl_id %in% mt_genes) %>%
  column_to_rownames('ensembl_id') %>%
  colSums()
adata$obs$total_counts <- adata$X %>%
  column_to_rownames('ensembl_id') %>%
  colSums()
adata$obs$pct_counts_mt <-
  adata$obs$total_mt_counts /  adata$obs$total_counts *  100

adata$var$total_counts <- adata$X %>% column_to_rownames('ensembl_id') %>%
  rowSums()

#### maybe don'¡t need this
y_genes <- adata$var %>%
  filter(chromosome_name == 'Y') %>%
  pull(ensembl_id) %>%
  unique()
adata$obs$total_y_counts <- adata$X %>%
  filter(ensembl_id %in% y_genes) %>%
  column_to_rownames('ensembl_id') %>%
  colSums()
adata$obs$pct_counts_y <-
  adata$obs$total_y_counts /  adata$obs$total_counts *  100

# calculate dim reductions
set_PCA(adata, nPcs = 180, pca_method = 'svd', scale_by = 'sample') #

set_UMAP(adata, pc_lim = kaisers_PCA_rule(adata$uns$pca), n_epochs = 1000)

sample_correlation <- adata$X  %>%
  column_to_rownames('ensembl_id') %>%
  cor(method = 'spearman')


library(purrr)
library(broom)
long_data <- adata$X %>%
  mutate(across(where(is.numeric), ~log2(.x + 1))) %>%
  pivot_longer(-ensembl_id, names_to = "ID", values_to = "expression") %>%
  left_join(select(adata$obs, ID, rat), by = "ID")


rat_variance_res <- long_data %>%
  group_by(ensembl_id) %>%
  do(glance(lm(expression ~ rat, data = .))) %>% # Uses broom to get R-squared
  ungroup() %>%
  arrange(desc(adj.r.squared)) # Higher R-squared = gene varies mostly due to rat
rat_variance_res <- rat_variance_res %>% left_join(adata$var %>% select(ensembl_id, gene_symbol, chromosome_name))


# Save dim reduction plots
save_dir <- "./figures/dim_reduction/sample/"
ensure_dir(save_dir)

plot_PCA(
  adata,
  color_by = 'organ_name',
  color_code = 'organ_color',
  alpha = 0.7,
  leg_ncol = 2
)
ggsave(file.path(save_dir, 'PCA.pdf'), width = 6, height = 8)


plot_UMAP(
  adata,
  color_by = 'organ_name',
  color_code = 'organ_color',
  alpha = 0.7,
  leg_ncol = 2
)
ggsave(file.path(save_dir, 'UMAP.pdf'), width = 6, height = 8)

#asign color code to based on entry under 'rat' column, based on MET color palette
adata$obs <- adata$obs %>%
  mutate(
    rat_color = case_when(
      rat == '1' ~ '#1f77b4',
      rat == '2' ~ '#17becf',
      rat == '3' ~ '#2ca02c',
      rat == 'M3' ~ '#9467bd',

      TRUE ~ '#000000'
    )
  )

plot_UMAP(
  adata,
  color_by = 'rat',
  color_code = 'rat_color',
  alpha = 0.7,
  leg_ncol = 2
)
ggsave(file.path(save_dir, 'rat_UMAP.pdf'), width = 6, height = 8)

p <- plot_UMAP_plotly(
  adata,
  color_by = 'organ_name',
  color_code = 'organ_color',
  hover_text_columns = c('sample_id', 'tissue_name'),
  size = 10,
  alpha = 1
)
htmlwidgets::saveWidget(
  p,
  file.path(save_dir, 'UMAP.html'),
  selfcontained = TRUE
)


# umap highligh not super needed atm, since QC sheets generated

# save_dir <- "./figures/dim_reduction/sample/by_tissue/UMAP_highlight"
# ensure_dir(save_dir)
#
# tissues <- unique(adata$obs$tissue_name)
#
# map(tissues, \(sel_tissue) {
#   plot_UMAP.highlight(adata, 'tissue_name', sel_tissue)
#   ggsave(file.path(save_dir, paste0( sel_tissue, '.pdf')), width = 6, height = 8)
#   })
#
# plot_UMAP.highlight(adata, 'tissue_name', sel_tissue)

get_corr <- function(wide_data, var_names = NULL) {
  if (is.null(var_names)) {
    var_names = names(adata.gene_tpm)[1]
  }
  sample_correlation <- wide_data %>%
    column_to_rownames(var_names) %>%
    cor(method = 'spearman')
  return(sample_correlation)
}

sample_to_tissue_correlation_plot <- function(
    sample_correlation,
    obs,
    element_of_interest,
    column_name
) {
  obs_name <- names(obs)[1]
  n_samples <- obs %>%
    filter(!!sym(column_name) == element_of_interest) %>%
    nrow()

  samples_of_interest <- obs %>%
    filter(!!sym(column_name) == element_of_interest) %>%
    pull(ID)

  #samples_of_interest.corr <-
  #  sample_correlation %>% as.data.frame() %>%
  #  rownames_to_column(obs_name) %>%
  #  select(1, samples_of_interest) %>%
  #  filter(!!sym(obs_name) %in% samples_of_interest) %>%
  #  column_to_rownames(obs_name) %>%
  #  pheatmap(treeheight_row = 0, treeheight_cols = 25, cluster_rows = n_samples > 1, cluster_cols = n_samples > 1)

  pl <- sample_correlation %>%
    as_tibble(rownames = 'ID') %>%
    select(1, samples_of_interest) %>%
    filter(!!sym(obs_name) %in% samples_of_interest) %>%
    gather(sample_id_target, corr, -1) %>%
    ggplot(aes(x = ID, y = sample_id_target, fill = corr)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(name = "Spearman rho") +
    theme_minimal()
  return(pl)
}
top_correlating_samples_to_tissue_plot <-
  function(sample_correlation, obs, column_name, element_of_interest) {
    obs_name <- names(obs)[1]
    n_samples <- obs %>%
      filter(!!sym(column_name) == element_of_interest) %>%
      nrow()

    samples_of_interest <- obs %>%
      filter(!!sym(column_name) == element_of_interest) %>%
      pull(ID)

    top_correlating_samples <-
      sample_correlation %>%
      as_tibble(rownames = obs_name) %>%
      select(c(obs_name, samples_of_interest)) %>%
      gather(sample_counter, roh, -1) %>%
      group_by(sample_counter) %>%
      arrange(desc(roh)) %>%
      slice_head(n = 10) %>%
      ungroup() %>%
      pull(ID) %>%
      unique()

    top_correlating_samples.corr <-
      sample_correlation %>%
      as.data.frame() %>%
      rownames_to_column(obs_name) %>%
      select(1, top_correlating_samples) %>%
      filter(!!sym(obs_name) %in% samples_of_interest) %>%
      column_to_rownames(obs_name) %>%
      pheatmap(
        treeheight_row = 0,
        treeheight_cols = 25,
        cluster_rows = n_samples > 1
      )

    return(top_correlating_samples.corr)
  }
obs_highlighted_scatter_plot <- function(
    obs,
    column_name,
    element_of_interest
) {

  obs_name <- names(obs)[1]
  plot_data <- obs %>%
    mutate(
      Color = case_when(
        .[[column_name]] == element_of_interest ~ "#ff0000",
        TRUE ~ "#a9a9a9"
      )
    )
  pl <- ggplot(data = plot_data) +
    geom_point(
      aes(pct_counts_mt, genes_detected, colour = Color),
      alpha = 0.6
    ) +
    scale_color_identity() +
    geom_text_repel(
      data = plot_data %>% filter(!!sym(column_name) == element_of_interest),
      aes(pct_counts_mt, genes_detected, label = !!sym(obs_name))
    ) +
    theme_classic()
  return(pl)
}

tpm_ridge_plot <- function(wide_data, obs, element_of_interest, column_name) {
  obs_name <- names(obs)[1]

  samples_of_interest <- obs %>%
    filter(!!sym(column_name) == element_of_interest) %>%
    pull(ID)

  pl <- wide_data %>%
    select(1, samples_of_interest) %>%
    gather(!!sym(obs_name), tpm, -1) %>%
    mutate(log1p_tpm = log1p(tpm)) %>%
    ggplot(aes(x = log1p_tpm, y = ID)) +
    stat_density_ridges(quantile_lines = TRUE, quantiles = 4)
  return(pl)
}


## QC plots


save_dir <- "./figures/QC/"
ensure_dir(save_dir)

plot_pal <- get_pal(adata$obs, 'organ_name', 'organ_color')

adata$obs %>%
  mutate(ID = factor(ID, levels = rev(ID_order))) %>%
  ggplot(aes(x = pct_counts_mt, y = ID, fill = organ_name)) +
  geom_col() +
  scale_fill_manual(values = plot_pal) +
  guides(fill = 'none') +
  ggtitle('Percentage of mitochondrial counts per sample')

ggsave(file.path(save_dir, 'mt_content.pdf'), width = 6, height = 40)

adata$var %>% filter(chromosome_name == 'X') %>% arrange(-total_counts)

adata$X %>% filter(ensembl_id == 'ENSRNOG00000060617') %>%
  gather(ID, counts, -ensembl_id) %>%
  left_join(adata$obs, by = 'ID') %>%
  mutate(ID = factor(ID, levels = ID_order)) %>%
  ggplot(aes(x = ID, y = counts, fill = organ_name)) +
  geom_col() +
  scale_fill_manual(values = plot_pal) +
  theme_minimal() +
  facet_grid(. ~ rat, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('Counts of UTY (Y gene) per sample')
ggsave(file.path(save_dir, 'ygene_exp_by_sex.pdf'), width = 40, height = 7)

adata$X %>% filter(ensembl_id == 'ENSRNOG00000055939') %>%
  gather(ID, counts, -ensembl_id) %>%
  left_join(adata$obs, by = 'ID') %>%
  mutate(ID = factor(ID, levels = ID_order)) %>%
  ggplot(aes(x = ID, y = counts, fill = organ_name)) +
  geom_col() +
  scale_fill_manual(values = plot_pal) +
  theme_minimal() +
  facet_grid(. ~ rat, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('Counts of Zcchc18 (X gene) per sample')
ggsave(file.path(save_dir, 'xgene_exp.pdf'), width = 40, height = 7)


sample_correlation %>%
  pheatmap::pheatmap() %>%
  as.ggplot() +
  ggtitle('Sample to sample Spearman correlation heatmap')

ggsave(file.path(save_dir, 'corr_pheatmap.pdf'), height = 28, width = 28)


gene_to_plot <- 'ENSRNOG00000063068'
gene_name <- adata$var %>% filter(ensembl_id == gene_to_plot) %>% pull(gene_symbol)
adata$X %>% filter(ensembl_id == gene_to_plot) %>%
  gather(ID, counts, -ensembl_id) %>%
  left_join(adata$obs, by = 'ID') %>%
  mutate(ID = factor(ID, levels = ID_order)) %>%
  ggplot(aes(x = ID, y = counts, fill = organ_name)) +
  geom_col() +
  scale_fill_manual(values = plot_pal) +
  theme_minimal() +
  facet_grid(. ~ rat, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(paste0('Counts of ', gene_name, '  per sample'))
ggsave(file.path(save_dir, paste0(gene_name, '_gene_exp.pdf')), width = 40, height = 7)




tissues <- unique(adata$obs$tissue_name)

for (selected_tissue in tissues) {
  print(selected_tissue)

  pl1 <- plot_PCA.highlight(
    adata,
    obs_column = 'tissue_name',
    element_to_highlight = selected_tissue,
    plot_names = TRUE,
    fixed_coords = TRUE,
    pca_id = 'pca',
    leg_ncol = 3
  )
  pl2 <- plot_UMAP.highlight(adata, 'tissue_name', selected_tissue)

  pl3 <- sample_to_tissue_correlation_plot(
    sample_correlation,
    adata$obs,
    selected_tissue,
    'tissue_name'
  )

  pl4 <- top_correlating_samples_to_tissue_plot(
    sample_correlation = sample_correlation,
    obs = adata$obs,
    column_name = 'tissue_name',
    element_of_interest = selected_tissue
  ) %>%
    as.ggplot()

  pl5 <- tpm_ridge_plot(
    adata.gene_tpm,
    adata_tpm$obs,
    element_of_interest = selected_tissue,
    column_name = 'tissue_name'
  )
  pl6 <- obs_highlighted_scatter_plot(adata$obs, 'tissue_name', selected_tissue)

  (pl1 + pl2) / (pl3 + pl4) / (pl5 + pl6)

  ggsave(
    paste0('./figures/QC/by_tissue/', selected_tissue, '.pdf'),
    height = 15,
    width = 15
  )
}



get_hs_homologs <- function(version = '109', dataset = 'rnorvegicus_gene_ensembl' ) {
  ensembl <- biomaRt::useEnsembl(biomart = 'genes',
                                 dataset = dataset,
                                 version = version)
  hs_homologs <- biomaRt::getBM(attributes=c( 'ensembl_gene_id', 'external_gene_name', 'hsapiens_homolog_ensembl_gene', 'hsapiens_homolog_associated_gene_name','hsapiens_homolog_orthology_type'),
                                mart=ensembl)

  hs_homologs <- hs_homologs |> mutate(gene_symbol = case_when(
    external_gene_name == '' ~ ensembl_gene_id,
    .default = external_gene_name))
  return(hs_homologs)
}


htm <- function(gene_name, hs_homologs = NULL){
  if (is.null(hs_homologs) && !exists("hs_homologs", envir = .GlobalEnv)) {
    hs_homologs <- get_hs_homologs(version = '109', dataset = 'rnorvegicus_gene_ensembl')
    # Optionally assign it to the global environment so it's available for future calls
    assign("hs_homologs", hs_homologs, envir = .GlobalEnv)
  } else if (is.null(hs_homologs)) {
    # If it's not passed as an argument, use the global hs_homologs
    hs_homologs <- get("hs_homologs", envir = .GlobalEnv)
  }

  homologs <- hs_homologs |> dplyr::filter((hsapiens_homolog_ensembl_gene == !!gene_name) |
                                             (hsapiens_homolog_associated_gene_name == !!gene_name) ) |> pull(gene_symbol)

  return(homologs)
}


html <- function(gene_list, var = NULL){
  homolog_list <- purrr::map(gene_list, htm)
  #homolog_list <- purrr::imap(homolog_list,  ~ set_names(.x, rep(.y, length(.x)))) %>% flatten()
  homolog_list <- unlist(homolog_list)
  #homolog_list <- unlist(purrr::map(gene_list, htm), recursive = FALSE, use.names = TRUE)
  if (!is.null(var)){
    homolog_list <- unlist(map(homolog_list, \(x)  if (x %in% var$gene_symbol) x ))
  }

  return(homolog_list)
}




expression_heatmap <- function(wide_data,
                               obs,
                               var,
                               gene_list,
                               sample_id_column,
                               group_by_column,
                               group_by_order,
                               color_by,
                               color_code
) {
  require(ComplexHeatmap)

  gene_list.ensemble <-
    var |> filter(gene_symbol %in% !!gene_list) |> pull(ensembl_id)

  matrix_data <-
    wide_data |> left_join(var) |>
    filter(ensembl_id %in% gene_list.ensemble) |>
    select(-ensembl_id) |>
    column_to_rownames('gene_symbol') |>
    log1p() |>
    t() |>
    scale() |>
    as.data.frame()

  obs <- obs |> mutate(sample_id = !!sym(sample_id_column),
                       group_by_column = !!sym(group_by_column),
                       color_by = !!sym(color_by),
                       color_code = !!sym(color_code))

  obs <- obs |> mutate(group_by_column = factor(group_by_column, levels = group_by_order))


  matrix_data <- matrix_data[obs |> pull(sample_id), gene_list]

  # obs |> select(tissue_name,region_tissue_name, organ_name, organ_color) |> distinct()  |>  arrange(organ_name, region_tissue_name, tissue_name) |> pull(organ_color)


  heatmap <- Heatmap(
    matrix_data,
    name = "Expression",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "top",
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    col = circlize::colorRamp2(
      c(-2, -1, 0, 1, 2, 3),
      c("white", "white", "#f0ed6e", "#f16f1f", "#300b5a", "#000000")
    ),
    heatmap_legend_param = list(title = "z-score(log1p(TPM))",
                                legend_direction = "vertical"),
    column_split = data.frame(factor(names(gene_list), levels = names(gene_list) |> unique())),
    row_split = data.frame(obs$group_by_column),
    border = TRUE,
    row_title_rot = 0,
    column_title_rot = 90,
    row_title_gp = gpar(
      fontface = 'bold'
      #fill = obs |>
      #  select(group_by_column,region_tissue_name, organ_name, organ_color) |> distinct()  |>
      #  #arrange(organ_name, region_tissue_name, group_by_column) |>
      #  pull(organ_color)
    ),
    left_annotation = rowAnnotation(
      group =  obs |> select(sample_id, color_by) |> pull(sample_id),
      col = list(group = setNames(obs$color_code, obs$sample_id)),
      show_legend = FALSE
    )
  )

  return(heatmap)
}



save_dir <- "./figures/QC/expression_heatmaps/"
ensure_dir(save_dir)

gene_list <- c(
  'adipocytes' = 'ADIPOQ',
  'adipocytes' = 'FABP4',
  'adipocytes' = 'LIPE',
  'adipocytes' = 'PLIN1',
  'mesothelial cells' = 'MSLN',
  'mesothelial cells' = 'ITLN1',
  'mesothelial cells' = 'UPK3B',
  'mesothelial cells' = 'WT1'

)

gene_list <- html(gene_list, adata$var)

heatmap <- expression_heatmap(wide_data = adata$X,
                              obs = adata$obs,
                              var = adata$var %>% select(ensembl_id, gene_symbol),
                              gene_list = gene_list,
                              sample_id_column = 'ID',
                              group_by_column = 'tissue_name',
                              group_by_order = adata$obs |> select(organ_name, tissue_rank, tissue_name) |> arrange(organ_name, tissue_rank, tissue_name) |> distinct() |>  pull(tissue_name),
                              color_by = 'organ_name',
                              color_code = 'organ_color'
)


pdf(file.path(save_dir, 'adipos_mesothelial_markers.pdf'),  height = 50, width = 12)
heatmap
dev.off()


gene_list <- c(
  'adipocytes' = 'ADIPOQ',
  'adipocytes' = 'FABP4',
  'adipocytes' = 'LIPE',
  'adipocytes' = 'PLIN1',
  'mesothelial cells' = 'MSLN',
  'mesothelial cells' = 'ITLN1',
  'mesothelial cells' = 'UPK3B',
  'mesothelial cells' = 'WT1',

  'Vascular ECs' = 'PECAM1',
  'Vascular ECs' = 'VWF',
  'Venous ECs' = 'ACKR1',
  'Venous ECs' = 'POSTN',
  'Arterial ECs' = 'HEY1',
  'Arterial ECs' = 'HEY2',
  'Arterial ECs' = 'SEMA3G',
  'Lymphatic EC' = ' LYVE1',
  'Lymphatic EC' = 'PROX1',
  'VSmooth muscle' = 'ACTA2',
  'VSmooth muscle' = 'CNN1',
  'VSmooth muscle' = 'MYH11',
  'Pericytes' = 'PDGFRB',
  'Pericytes' = 'RGS5',
  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'COL1A1',
  'tear gland' = 'LYZ',
  #'tear gland' = 'LCN1',
  'tear gland' = 'LTF',
  'tear gland' = 'LACRT',
  'tear gland' = 'AQP5',
  'tear gland' = 'SLC12A2',
  'tear gland' = 'NKCC1',
  'tear gland' = 'EPCAM',
  'tear gland' = 'KRT18',
  'tear gland' = 'KRT19',
  'tear gland' = 'CHRM3',

  'Squamous epithelium' = 'KRT14',
  'Squamous epithelium' = 'IVL',
  'Squamous epithelium' = 'KRT1',
  'Squamous epithelium' = 'KRT10',

  'Mucosa ligning' = 'KRT13',
  'Mucosa ligning' = 'KRT4',

  'Mesothelial cells' = 'MSLN',
  'Mesothelial cells' = 'ITLN1',



  'Smooth muscle' = 'ACTA2',
  'Smooth muscle' = 'ACTG2',
  'Smooth muscle' = 'DES',
  'Smooth muscle' = 'MYH11',
  'Heart' = 'HAS1',
  'Heart' = 'RYR2',
  'Heart' = 'LEPR',
  'Synovial' = 'POSTN',
  'Synovial' = 'PRG4',
  'Hair keratins' = 'KRT71',
  'Hair keratins' = 'KRT25',
  'Hair keratins' = 'KRT27',
  'Hair keratins' = 'LGR5',
  'Hair keratins' = 'LHCGR',
  'No hair keratins' = 'KRT9'

)

gene_list <- html(gene_list, adata$var)

heatmap <- expression_heatmap(wide_data = adata$X,
                              obs = adata$obs,
                              var = adata$var %>% select(ensembl_id, gene_symbol),
                              gene_list = gene_list,
                              sample_id_column = 'ID',
                              group_by_column = 'tissue_name',
                              group_by_order = adata$obs |> select(organ_name, tissue_name) |> arrange(organ_name, tissue_name) |> distinct() |>  pull(tissue_name),
                              color_by = 'organ_name',
                              color_code = 'organ_color'
)


pdf(file.path(save_dir, 'teargland_markers.pdf'),  height = 40, width = 12)
heatmap
dev.off()



markers_ESj <- read_delim("~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Documents/SciLifeDrive/macaque_bulk/macaque_bulk_git/experiment_plots_before_20250129/checking_markers/markers_to_check_20240912_ESj.csv", delim = ';')
colnames(markers_ESj) <- c( "order","organ_system", "cell_type", "marker", "source", "x1","x2"    )
markers_ESj <- markers_ESj |> filter(marker != 'KRT6A')
hs_markers <- setNames( markers_ESj |> select(cell_type, marker) |> distinct() |> drop_na() |> pull(marker) ,
                        markers_ESj |> select(cell_type, marker) |> distinct() |> drop_na() |> pull(cell_type))

gene_list <- html(hs_markers, adata$var)

heatmap <- expression_heatmap(wide_data = adata$X,
                              obs = adata$obs,
                              var = adata$var %>% select(ensembl_id, gene_symbol),
                              gene_list = gene_list,
                              sample_id_column = 'ID',
                              group_by_column = 'tissue_name',
                              group_by_order = adata$obs |> select(organ_name, tissue_rank, tissue_name) |> arrange(organ_name, tissue_rank, tissue_name) |> distinct() |>  pull(tissue_name),
                              color_by = 'organ_name',
                              color_code = 'organ_color'
)


pdf(file.path(save_dir, 'esj_markers.pdf'),  height = 50, width = 60)
heatmap
dev.off()




gene_list <- c(

  'Glutamtergic Neurons' = 'SLC17A6',
  'Glutamtergic Neurons' = 'SLC17A7',
  'Glutamtergic Neurons' = 'VGLUT1',
  'Glutamtergic Neurons' = 'SLC17A8',
  'GABA Neurons' = 'GAD1',
  'GABA Neurons' = 'GAD2',
  'GABA Neurons' = 'SLC32A1',
  'GABA Neurons' = 'SLC6A1',
  'Glycine Neurons' = 'SLC6A9',
  'Glycine Neurons' = 'SLC32A1',
  'Catecholamine Neurons' = 'TH',
  'Catecholamine Neurons' = 'DDC',
  'Dopamine Neurons' = 'SLC6A3',
  'Noradrenergic Neurons' = 'DBH',
  'Noradrenergic Neurons' = 'SLC6A2',
  'Adrenergic Neurons' = 'PNMT',
  'Serotonergic Neurons' = 'TPH2',
  'Serotonergic Neurons' = 'SLC6A4',
  'Cholinergic Neurons' = 'SLC5A7',
  'Cholinergic Neurons' = 'CHAT',
  'Cholinergic Neurons' = 'ACHE',
  'NPY Neurons' = 'NPY',
  'CRH Neurons' = 'CRH',
  'CRH Neurons' = 'CRHR1',
  'CRH Neurons' = 'CRHR2',
  'OXT Neurons' = 'OXTR',
  'OXT Neurons' = 'OXT',
  'OXT Neurons' = 'NEUROD6',
  'NMDA receptors' = 'GRIN1',
  'NMDA receptors' = 'GRIN2',
  'AMPA receptors' = 'GRIA1',
  'AMPA receptors' = 'GRIA2',
  'LAMP5' = 'LAMP5',
  'Interneurons' = 'SST',
  'Interneurons' = 'VIP',
  'Interneurons' = 'PVALB',
  'CALB1' = 'CALB1',
  'PV' = 'PV',
  'CR' = 'CR',
  'Substance P' = 'TAC1',
  #'CALM' = 'CALM1',
 # 'CALM' = 'CALM2',
 # 'CALM' = 'CALM3',
  'Myeline' = 'MBP',
  'Myeline' = 'PLP1',
  'Myeline' = 'MAG',
  'axon proteins' = 'NEFH',
  'axon proteins' = 'NEFM',
  'axon proteins' = 'NEFL',
  'axon proteins' = 'DCLK1',
  'axon proteins' = 'AP2M1',
  'axon guidance' = 'SEMA3A',
  'axon guidance' = 'SEMA4D',
  'axon guidance' = 'EPHA4',
  'axon guidance' = 'EPHB2',
  'axon guidance' = 'NTN1',
  'axon guidance' = 'NTN3',
  'axon guidance' = 'SLIT2',
  'axon guidance' = 'RHOA',
  'axon guidance' = 'RAC1',


  # amugdala is hard

  #basal ganglia
  'Basal ganglia - Striatum' = 'GPR88',
  'Basal ganglia - Striatum' = 'ADORA2A',
  'Basal ganglia - Striatum' =  'NPPB'   ,
  'Basal ganglia - Striatum' =    'DACT2'  ,

  #  'Basal ganglia - GP' =  'DRD1'   ,
  #  'Basal ganglia - GP' =  'DRD2'   ,
  #  'Basal ganglia - GP' =  'PENK'   ,

  #'Basal ganglia - Striatum' =  'SBSN'    ,
  #'Basal ganglia - Striatum' =  'TRBV25-1'    ,
  'Septum' = 'PRDM12', #hypothalamus too
  #'Basal ganglia - Septum' = 'RTN4RL2',
  #'Basal ganglia - Septum' = 'NEUROD6',
  'Cerebellum' = 'PRR35',
  'Cerebellum' = 'CDH15',
  'Cerebellum' = 'GABRA6',
  'Cerebellum' = 'BARHL2',
  'Cerebellum' = 'PMIS2',
  'Cerebellum' = 'CBLN3',
  'Cerebellum' = 'NRK',
  'Cerebellum' = 'EFCAB9',
  'Cerebellum' = 'KASH5',
  'Cortex - claustrum' = 'NR4A2',

  #'Hippoc' = 'SPATA31C2', #hpa enriched
  #'Hippoc' = 'SPATA31C1', #hpa enriched
  #'Hippoc' = 'SPATA31E1', #hpa enriched
  #'Hippoc' = 'C1QTNF9B',#hpa enriched
  #'Hippoc' = 'NEUROG3',
  #'Hippoc' = 'CHRNA1',
  #'Hippoc' = 'NTF3',
  #'Hippoc' = 'TNNT2',
  #'Hippoc' = 'ARC',
  #'Hippoc' = 'FNDC1',
  #'Hippoc' = 'SCGN',#stronger in hypoth
  #'Subiculum' = 'TEX101',
  #'Subiculum' = 'ZFHX3',
  #'Subiculum' = 'SKOR2',
  #'Subiculum' = 'IGHJ3',
  #'Subiculum' = 'ENSG00000288681',
  #'Subiculum' = 'FEZF2',
  #'entorhi' = 'ROS1',
  #'entorhi' = 'CPLX3',
  #'entorhi' = 'GUCA1C'
  'Hypothalamus' = 'GHRH',
  'Hypothalamus' = 'NR5A1',
  'Hypothalamus' = 'FEZF1',
  'Hypothalamus' = 'SIX6',
  'Hypothalamus' = 'AGRP',
  'Hypothalamus' = 'NPVF',
  'Hypothalamus' = 'GAL',
  'Hypothalamus' = 'MC3R',
  'Hypothalamus' = 'GZMB',
  'Hypothalamus' = 'NPY',
  'Hypothalamus' = 'POMC',
  'Hypothalamus' = 'GSX1',
  'Hypothalamus p' = 'AVP',
  'Hypothalamus p' = 'OXT',

  'Midbrain' = 'EN1',
  #'Midbrain' = 'TFF3',
  #'Midbrain' = 'CLEC3A',
  'Midbrain' = 'SLC6A4',
  'Midbrain' = 'TPH2',
  #'Midbrain' = 'FEV',

  'Midbrain - sn' = 'SLC6A3',
  'Midbrain - sn' = 'PITX3',
  'Midbrain - sn' = 'NTSR1',
  'Midbrain - sn' = 'SDC1',

  'Midbrain colliculus' = 'TFAP2D',
  'Midbrain colliculus' = 'TFAP2B',
  'Midbrain colliculus' = 'POU4F2',


  'Olfactory' = 'SCGN',
  'Olfactory' = 'DCX',
  'Olfactory' = 'S100A5',
  'Medulla' = 'GCG',
  'Medulla' = 'PRLH',
  'Pons' = 'SLC6A2',
  'Pons' = 'DBH',
  'Pons' = 'RLN3',
  'Pons' = 'HOXC4',
  #'Pons' = 'PHOX2A',
  #'Pons' = 'C1QL4',
  #'Pons' = 'EXD1',
  'TH Domaine' = 'TH',
  'spinal cord' = 'HOXC10',
  'spinal cord' = 'HOXD9',
  'spinal cord' = 'MNX1',
  'spinal cord' = 'HOXA9',


  'Thalamus' = 'MYF6',
  'Thalamus' = 'RGS16',
  'Thalamus' = 'TAFA4',


  'White matter' = 'CLDN11',
  'White matter' = 'CNP',
  'HPCA' = 'HPCA'

)


gene_list <- html(gene_list, adata$var)

heatmap <- expression_heatmap(wide_data = adata$X,
                              obs = adata$obs,
                              var = adata$var %>% select(ensembl_id, gene_symbol),
                              gene_list = gene_list,
                              sample_id_column = 'ID',
                              group_by_column = 'tissue_name',
                              group_by_order = adata$obs |> select(organ_name, tissue_rank, tissue_name) |> arrange(organ_name, tissue_rank, tissue_name) |> distinct() |>  pull(tissue_name),
                              color_by = 'organ_name',
                              color_code = 'organ_color'
)


pdf(file.path(save_dir, 'brain_markers.pdf'),  height = 50, width = 50)
heatmap
dev.off()




gene_list <- c(
  'Amygdala' = 'OXTR',
  'BG/Septum' = 'TRPC4',
  'Cortex' = 'CNTN3',
  'Cortex' = 'CUX1',
  'Cortex' = 'CUX2',
  'Cortex' = 'LAMP5',
  'Cerebellum' = 'EOMES',
  'Cerebellum' = 'GABRA6',
  'Cerebellum' = 'PCP2',
  'Hippocampus' = 'NEUROD6',
  'Hypothalamus' = 'AVP',
  'Hypothalamus' = 'FEZF1',
  'Hypothalamus' = 'OXT',
  'Midbrain' = 'EN1',
  'Midbrain' = 'POU4F2',
  'Olfactory bulb' = 'SCGN',
  'Olfactory bulb' = 'S100A5',
  'Medulla' = 'PRLH',
  'Pons' = 'DBH',
  'Pons' = 'SLC6A2',
  'TH DOpa'  = 'TH',
  'spinal cord' = 'HOXC10',
  'spinal cord' = 'HOXD9',
  'Thalamus' = 'RGS16',
  'Thalamus' = 'TAFA4',
  'CornealFibro' = 'KERA',
  'Melanocytes' = 'MLANA',
  'Melanocytes' = 'PAX3',
  'ratment Epithelia' = 'RBP1',
  'RodPhotoreceptors' = 'RHO',
  'RetinalAmacrine' = 'C1QL2',
  'RetinalHorizontal' = 'ONECUT1'

)


gene_list <- html(gene_list, adata$var)

heatmap <- expression_heatmap(wide_data = adata$X,
                              obs = adata$obs,
                              var = adata$var %>% select(ensembl_id, gene_symbol),
                              gene_list = gene_list,
                              sample_id_column = 'ID',
                              group_by_column = 'tissue_name',
                              group_by_order =adata$obs |> select(organ_name, tissue_rank, tissue_name) |> arrange(organ_name, tissue_rank, tissue_name) |> distinct() |>  pull(tissue_name),
                              color_by = 'organ_name',
                              color_code = 'organ_color'
)


pdf(file.path(save_dir, 'brain_markers_2.pdf'),  height = 50, width = 15)
heatmap
dev.off()




gene_list <- c(
  'Adrenal cortex cells' = 'MC2R',
  'Adrenal cortex cells' =  'STAR',
  'Adrenal medulla cells' = 'CHGB',
  'Cortico/Gonado/Lacto/Somatotrophs' = 'POMC',
  'Cortico/Gonado/Lacto/Somatotrophs' = 'FSHB',
  'Cortico/Gonado/Lacto/Somatotrophs' = 'DIO2',
  'Thyroid gland' = 'TG',
  'Thyroid gland' = 'TPO',
  'Alveolar cells' = 'UPK3B',
  'Alveolar cells' = 'SFTPC',
  'Alveolar cells' = 'NAPSA',

  'Ciliated cells' = 'FOXJ1',
  'Ciliated cells' = 'RSPH1',
  'Mucus' = 'MUC5B',
  'Mucus' = 'BPIFA1',
  'Skeletal muscle' = 'CKM',
  'Skeletal muscle' = 'MB',
  'Skeletal muscle' = 'TNNI1',
  'Skeletal muscle' = 'MYH7',

  'Basal cells' = 'KRT5',
  'Basal cells' = 'TP63',
  'Basal cells' = 'KRT14',
  'Basal cells' = 'KRT15',


  'Olfactory epithelium' = 'OMP',
  'Olfactory epithelium' = 'ADCY3',
  'Olfactory epithelium' = 'KRT18',

  'Olfacotry receptors' = 'OR4X2',
  'Olfacotry receptors' = 'OR2D2',
  'Olfacotry receptors' = 'OR10H1',

  'Squamous epithelium' = 'KRT14',
  'Squamous epithelium' = 'IVL',
  'Squamous epithelium' = 'KRT1',
  'Squamous epithelium' = 'KRT10',

  'Mucosa ligning' = 'KRT13',
  'Mucosa ligning' = 'KRT4',

  'Mesothelial cells' = 'MSLN',
  'Mesothelial cells' = 'ITLN1',



  'Smooth muscle' = 'ACTA2',
  'Smooth muscle' = 'ACTG2',
  'Smooth muscle' = 'DES',
  'Smooth muscle' = 'MYH11',
  'Heart' = 'HAS1',
  'Heart' = 'RYR2',
  'Heart' = 'LEPR',
  'Synovial' = 'POSTN',
  'Synovial' = 'PRG4',
  'Hair keratins' = 'KRT71',
  'Hair keratins' = 'KRT25',
  'Hair keratins' = 'KRT27',
  'Hair keratins' = 'LGR5',
  'Hair keratins' = 'LHCGR',
  'No hair keratins' = 'KRT9',
  'Enterocytes' = 'ANPEP',
  'Enterocytes' = 'FABP2',
  'Enterocytes' = 'RBP2',
  'Colonocytes' = 'CA2',
  'Colonocytes' = 'FABP1',
  'Colonocytes' = 'SLC26A2',

  'Pancreas' = 'PRSS1',
  'Pancreas' = 'PRSS2',
  'Pancreas' = 'SST',
  'Pancreas' = 'INS',
  'Adipocytes' = 'ADIPOQ',
  'Adipocytes' = 'LIPE',

  'Urothelial cells' = 'KRT7',
  'Urothelial cells' = 'UPK1A',
  'Urothelial cells' = 'UPK3A',

  'Epydimis pcells' = 'DEFB19',
  'Epydimis pcells' = 'DEFB113',
  'Spermatocytes' = 'SPATA3',
  # 'Spermatocytes' = 'SYCP3',
  'Spermatocytes' = 'MEIOB',
  'Prostate' = 'ratR',
  'Prostate' = 'SCGB3A1',


  'breast' = 'LALBA',
  'breast' = 'OXTR',

  'ooocytes' = 'DAZL',
  'ooocytes' = 'FIGLA',
  'ooocytes' = 'ZP4',

  'Vascular ECs' = 'PECAM1',
  'Vascular ECs' = 'VWF',
  'Venous ECs' = 'ACKR1',
  'Venous ECs' = 'POSTN',
  'Arterial ECs' = 'HEY1',
  'Arterial ECs' = 'HEY2',
  'Arterial ECs' = 'SEMA3G',
  'Lymphatic EC' = ' LYVE1',
  'Lymphatic EC' = 'PROX1',
  'VSmooth muscle' = 'ACTA2',
  'VSmooth muscle' = 'CNN1',
  'VSmooth muscle' = 'MYH11',
  'Pericytes' = 'PDGFRB',
  'Pericytes' = 'RGS5',
  'Fibroblasts' = 'DCN',
  'Fibroblasts' = 'COL1A1',

  'T-cells' = 'CD3D',

  'T-cells' = 'CD4',

  # 'neurotrphil progenitors' = 'CAMP',
  'neurotrphil progenitors' = 'CEACAM8',
  'Monocyte progenitors' = 'CSF1R',
  'Monocyte progenitors' = 'IRF8',
  'Monocytes'= 'S100A8',
  'Monocytes'= 'FCN1',

  'Lacrimal acinar cells' = 'BIPFA2',
  'Lacrimal acinar cells' = 'LTF',
  'Lacrimal acinar cells' = 'PIP',
  'Lacrimal acinar cells' = 'AMY1',
  'Cartilage' = 'COL2A1',
  'Cartilage' = 'ACAN',
  'Cartilage' = 'HAPLN1'


)


gene_list <- html(gene_list, adata$var)

#append to gene list
#gene_list <- c( gene_list,'Seminal vescile' =  'Svs4', 'Seminal vescile' =  'Svs5', 'Seminal vescile' =  'Svs6')
heatmap <- expression_heatmap(wide_data = adata$X,
                              obs = adata$obs,
                              var = adata$var %>% select(ensembl_id, gene_symbol),
                              gene_list = gene_list,
                              sample_id_column = 'ID',
                              group_by_column = 'tissue_name',
                              group_by_order = adata$obs |> select(organ_name, tissue_rank, tissue_name) |> arrange(organ_name, tissue_rank, tissue_name) |> distinct() |>  pull(tissue_name),
                              color_by = 'organ_name',
                              color_code = 'organ_color'
)


pdf(file.path(save_dir, 'endocrine_markers.pdf'),  height = 50, width = 30)
heatmap
dev.off()

