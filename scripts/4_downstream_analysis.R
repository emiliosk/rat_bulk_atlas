# overview categories

devtools::load_all("/Users/emilioskarwan/Documents/SciLifeDrive/AnnDatR")
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(ggplot2)

#library(patchwork)


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
# load in data

# Open Normalized pseudobulk cluster data
adata <- AnnDatR$new(
  prefix_name = "tissue",
  layer = "nTPM",
  var_names = "ensembl_id",
  file_dir = data_loc
)
adata$obs  <- adata$obs %>% mutate(organ_name = factor(organ_name, levels = organ_order))
adata$obs %>% arrange(organ_name)
tissue_order <- adata$obs %>%
  arrange(
    organ_name,
    consensus_tissue_name,
    #tissue_group,
    # tissue_rank,
    #  region_tissue_name,
    tissue_name
  ) %>%
  pull(tissue_name) %>%
  unique()

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


set_PCA(adata, nPcs = 100, pca_method = 'svd', scale_by = 'sample') #
set_UMAP(adata, pc_lim = kaisers_PCA_rule(adata$uns$pca), n_epochs = 1000)
set_dendrogram(adata)

hpa_gene_classification(adata, enr_fold = 4, max_group_n = 5, det_lim = 1, inplace = T)




#----- Figures -----

save_dir <- "./figures/dim_reduction/tissue/"
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


p <- plot_UMAP_plotly(
  adata,
  color_by = 'organ_name',
  color_code = 'organ_color',
  hover_text_columns = c(),
  size = 10,
  alpha = 1
)
htmlwidgets::saveWidget(
  p,
  file.path(save_dir, 'UMAP.html'),
  selfcontained = TRUE
)

plot_cor_pheatmap(adata, color_by = 'organ_name', color_code = 'organ_color') %>% ggplotify::as.ggplot()
ggsave(file.path(save_dir, 'correlation.pdf'), width = 16, height = 12)

plot_dendrogram(adata, color_code = 'organ_color', x_expansion = 0.5, y_expansion = 0.01 )
ggsave(file.path(save_dir, 'dendrogram.pdf'), width = 10, height = 12)




#----- Category overviews -----
plot_specificity_barplot(adata$var, title = 'Rat: Speceificity category per tissue')
ggsave( file.path(save_dir, 'speceficity_barplot.pdf'), height = 11, width = 5)

plot_specificity_tau(adata$var,  title = 'Rat: Tau across speceificity category')
ggsave( file.path(save_dir, 'tauspeceficity_violinplot.pdf'), height = 5, width = 6)

plot_distribution_barplot(adata$var, title = 'Rat: Distribution category per tissue')
ggsave( file.path(save_dir, 'distribution_barplot.pdf'), height = 11, width = 5)

plot_alluvial_categories(adata$var)
ggsave( file.path(save_dir, 'alluvial categories.pdf'), height = 4, width = 7)





plot_data <- adata$X %>% column_to_rownames('ensembl_id')
colnames(plot_data) <- stringr::str_to_sentence(colnames(plot_data))
tissue_dendro_ward <- hclust4RNAseq_ward(plot_data)

circular_dendrogram_retinastyle_2(
  clust = tissue_dendro_ward,
  color_mapping = adata$obs %>% select(tissue_name, organ_color) %>%
    mutate(tissue_name = stringr::str_to_sentence(tissue_name)),
  label_col = "tissue_name",
  color_col = "organ_color",
  scale_expansion = c(0.6, 0.6),
  text_size = 2.4,
  width_range = c(0.5, 6),
  arc_strength = 0.4,
  default_color = "gray80")

ggsave( file.path(save_dir, 'retinogramm.pdf'), height = 7, width = 7)


adata$obs


# Consensus level -------------------------------------------------------------



adata <- AnnDatR$new(
  prefix_name = "consensus",
  layer = "nTPM",
  var_names = "ensembl_id",
  file_dir = data_loc
)
adata$obs  <- adata$obs %>% mutate(organ_name = factor(organ_name, levels = organ_order))



set_PCA(adata, nPcs = 80, pca_method = 'svd', scale_by = 'sample') #
set_UMAP(adata, pc_lim = kaisers_PCA_rule(adata$uns$pca), n_epochs = 1000)
set_dendrogram(adata)

hpa_gene_classification(adata, enr_fold = 4, max_group_n = 5, det_lim = 1, inplace = T)




save_dir <- "./figures/dim_reduction/consensus/"
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


p <- plot_UMAP_plotly(
  adata,
  color_by = 'organ_name',
  color_code = 'organ_color',
  hover_text_columns = c(),
  size = 10,
  alpha = 1
)
htmlwidgets::saveWidget(
  p,
  file.path(save_dir, 'UMAP.html'),
  selfcontained = TRUE
)

plot_cor_pheatmap(adata, color_by = 'organ_name', color_code = 'organ_color') %>% ggplotify::as.ggplot()
ggsave(file.path(save_dir, 'correlation.pdf'), width = 16, height = 12)

AnnDatR:::plot_dendrogram(adata, color_code = 'organ_color', x_expansion = 0.5, y_expansion = 0.01 )
ggsave(file.path(save_dir, 'dendrogram.pdf'), width = 10, height = 12)




#----- Category overviews -----
plot_specificity_barplot(adata$var, title = 'Rat: Speceificity category per tissue')
ggsave( file.path(save_dir, 'speceficity_barplot.pdf'), height = 11, width = 5)

plot_specificity_tau(adata$var,  title = 'Rat: Tau across speceificity category')
ggsave( file.path(save_dir, 'tauspeceficity_violinplot.pdf'), height = 5, width = 6)

plot_distribution_barplot(adata$var, title = 'Rat: Distribution category per tissue')
ggsave( file.path(save_dir, 'distribution_barplot.pdf'), height = 11, width = 5)

plot_alluvial_categories(adata$var)
ggsave( file.path(save_dir, 'alluvial categories.pdf'), height = 4, width = 7)



plot_2_column_alluvial(
  adata$var %>% mutate(spec_category = stringr::str_to_sentence(spec_category),
                       dist_category = stringr::str_to_sentence(dist_category)),
  id_column = 'ensembl_id',
  columns_to_plot  = c('spec_category', 'dist_category'),
  column_labels = c('Specificity', 'Distribution'),
  order_values = rev(names(gene_category_pal)),
  color_pallets = gene_category_pal,
  plot_title = 'Comparison: MS vs RNA specificity catoegires',
  legend_title = NULL,
  alluvial_width = 0.4,
  show_labels = TRUE
)
ggsave( file.path(save_dir, 'alluvial categories_2.pdf'), height = 6, width = 6)




plot_specificity_barplot(adata$var, title = 'Rat: Speceificity category per tissue') +
  coord_flip() +
  ggtitle('Rat: Speceificity category per tissue') +
  theme(
    axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)
  )


plot_specificity_barplot_wlabel(adata$var, title = 'Rat: Speceificity category per tissue') +
  theme(
    panel.background = ggplot2::element_blank(),   # Removes the inner plot background
    plot.background = ggplot2::element_blank(),    # Removes the outer background
    panel.grid.major = ggplot2::element_blank(),   # Removes major gridlines
    panel.grid.minor = ggplot2::element_blank(),   # Removes minor gridlines
  )
ggsave( file.path(save_dir, 'speceficity_barplot.pdf'), height = 7, width = 6)


plot_specificity_barplot(adata$var, title = 'Rat: Speceificity category per tissue') +
  ggplot2::coord_cartesian() +
  coord_flip() +
  ggtitle('Rat: Speceificity category per tissue') +
  theme(
    axis.text.x = ggplot2::element_text(angle = 60, vjust = 1, hjust = 1),
    panel.background = ggplot2::element_blank(),   # Removes the inner plot background
    plot.background = ggplot2::element_blank(),    # Removes the outer background
    panel.grid.major = ggplot2::element_blank(),   # Removes major gridlines
    panel.grid.minor = ggplot2::element_blank(),   # Removes minor gridlines
  )
ggsave( file.path(save_dir, 'speceficity_barplot_r.pdf'), height = 6, width = 8)
ggsave( file.path(save_dir, 'speceficity_barplot_r.png'), height = 6, width = 8)


specificity_data <- plot_specificity_barplot(adata$var, title = 'Rat: Speceificity category per tissue', return_data = T)
specificity_data %>% filter(enriched_tissues == 'olfactory epithelium')



# Filtered samples level -------------------------------------------------------------

adata <- AnnDatR$new(
  prefix_name = "filtered_samples",
  layer = "nTPM",
  var_names = "ensembl_id",
  file_dir = data_loc
)
#adata$obs  <- adata$obs %>% rename(organ_rank = organ_order, tissue_rank = tissue_order)
#adata$obs  <- adata$obs %>% mutate(organ_name = factor(organ_name, levels = organ_order))



set_PCA(adata, nPcs = 170, pca_method = 'svd', scale_by = 'sample') #
set_UMAP(adata, pc_lim = kaisers_PCA_rule(adata$uns$pca), n_epochs = 1000)
set_cor(adata)





save_dir <- "./figures/dim_reduction/filtered_samples/"
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


p <- plot_UMAP_plotly(
  adata,
  color_by = 'organ_name',
  color_code = 'organ_color',
  hover_text_columns = c(),
  size = 10,
  alpha = 1
)
htmlwidgets::saveWidget(
  p,
  file.path(save_dir, 'UMAP.html'),
  selfcontained = TRUE
)

plot_cor_pheatmap(adata, color_by = 'organ_name', color_code = 'organ_color') %>% ggplotify::as.ggplot()
ggsave(file.path(save_dir, 'correlation.pdf'), width = 35, height = 30)




save_directory = save_dir

file_name_sspearman = file.path(save_directory, 'spearman_sample.csv')
if(file.exists(file_name_sspearman) ){
  sample_tmm_spearman <- read_csv(file_name_sspearman)
} else {
  adata$obsm$correlation <-  adata$X %>%
    column_to_rownames(adata$var_names_col) %>%
    cor(method = "spearman", use = "pairwise.complete.obs") %>%
    as.data.frame() %>%
    as_tibble(rownames = "sample_name")
  write_csv(as.data.frame(adata$obsm$correlation) %>% as_tibble(),file_name_sspearman)
  sample_tmm_spearman <- read_csv(file_name_sspearman)
}
metadata = adata$obs
#metadata = metadata %>% rename(ID = sample_id) %>% select(-id)


correlation_to_different_organs <- tibble(sample_name = c(), correlation = c())
for (sample in metadata$ID){
  sample_organ <- metadata %>% filter(ID == sample) %>% .$organ_name %>% unique()
  different_organ_samples <- metadata %>% filter(organ_name != sample_organ) %>% .$ID
  sample_correlation <- sample_tmm_spearman %>% filter(sample_name %in% different_organ_samples) %>% select("sample_name", sample) %>% rename("correlation" = sample)
  correlation_to_different_organs <- rbind(correlation_to_different_organs, sample_correlation)
}

correlation_to_same_organ <- tibble(sample_name = c(), correlation = c())
for (sample in metadata$ID){
  sample_organ <- metadata %>% filter(ID == sample) %>% .$organ_name %>% unique()
  same_organ_samples <- metadata %>% filter(organ_name == sample_organ) %>% filter(ID != sample) %>% .$ID
  sample_correlation <- sample_tmm_spearman %>% filter(sample_name %in% same_organ_samples) %>% select("sample_name", sample) %>% rename("correlation" = sample)
  correlation_to_same_organ <- rbind(correlation_to_same_organ, sample_correlation)
}

correlation_to_same_tissue <- tibble(sample_name = c(), correlation = c())
for (sample in metadata$ID){
  sample_tissue <- metadata %>% filter(ID == sample) %>% .$tissue_name %>% unique()
  same_tissue_samples <- metadata %>% filter(tissue_name == sample_tissue) %>% filter(ID != sample) %>% .$ID
  sample_correlation <- sample_tmm_spearman %>% filter(sample_name %in% same_tissue_samples) %>% select("sample_name", sample) %>% rename("correlation" = sample)
  correlation_to_same_tissue <- rbind(correlation_to_same_tissue, sample_correlation)
}

correlation_to_same_tissue_by_tissue <- tibble(sample_name = c(), correlation = c(), tissue_name = c())
for (sample in metadata$ID){
  sample_tissue <- metadata %>% filter(ID == sample) %>% .$tissue_name %>% unique()
  same_tissue_samples <- metadata %>% filter(tissue_name == sample_tissue) %>% .$ID
  sample_correlation <- sample_tmm_spearman %>% filter(sample_name %in% same_tissue_samples) %>% select("sample_name", sample) %>% filter(sample_name != sample) %>%  left_join(metadata %>% select(ID, tissue_name), by= c("sample_name" = "ID")) %>% rename("correlation" = sample)
  correlation_to_same_tissue_by_tissue <- rbind(correlation_to_same_tissue_by_tissue, sample_correlation)
}

correlation_to_same_tissue_by_tissue <- correlation_to_same_tissue_by_tissue %>%
  mutate(tissue_name = stringr::str_to_sentence(tissue_name)) %>%
  group_by(tissue_name) %>%
  mutate(min = min(correlation)) %>%
  ungroup() %>%
  arrange(min)


# correlation_to_same_ctissue_by_ctissue <- tibble(
#   sample_name = c(),
#   correlation = c(),
#   consensus_tissue_name = c()
# )
#
# for (sample in metadata$ID) {
#   sample_tissue <- metadata %>% filter(ID == sample) %>% .$consensus_tissue_name %>% unique()
#   same_tissue_samples <- metadata %>% filter(consensus_tissue_name == sample_tissue) %>% .$ID
#   sample_correlation <- sample_tmm_spearman %>% filter(sample_name %in% same_tissue_samples) %>%
#     select("sample_name", sample) %>%
#     filter(sample_name != sample) %>%
#     left_join(metadata %>% select(ID, consensus_tissue_name),                                                                                                                                                                                 by = c("sample_name" = "ID")) %>% rename("correlation" = sample)
#   correlation_to_same_ctissue_by_ctissue <- rbind(correlation_to_same_ctissue_by_ctissue, sample_correlation)
# }
#
# correlation_to_same_ctissue_by_ctissue <- correlation_to_same_ctissue_by_ctissue %>%
#   mutate(consensus_tissue_name = stringr::str_to_sentence(consensus_tissue_name)) %>%
#   group_by(consensus_tissue_name) %>%
#   mutate(min = min(correlation)) %>%
#   ungroup() %>%
#   arrange(min)
#



p1 <- correlation_to_different_organs %>%
  ggplot(aes(correlation)) +
  geom_histogram(bins = 100)  + xlim(0.5,1)+
  theme_classic() +
  theme(panel.background = element_rect("gray90"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y.left = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
  ggtitle("Different organs")

p2 <- correlation_to_same_organ %>%
  ggplot(aes(correlation)) +
  geom_histogram(bins = 100) + xlim(0.5,1) +
  theme_classic() +
  theme(panel.background = element_rect("gray90"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
  ggtitle("Same organs")

p3 <- correlation_to_same_tissue %>%
  ggplot(aes(correlation)) +
  geom_histogram(bins = 100) + xlim(0.5,1) +
  theme_classic() +
  theme(panel.background = element_rect("gray90"),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
  ggtitle("Same tissue")

p4 <- correlation_to_same_tissue_by_tissue %>%
  mutate(tissue_name = factor(tissue_name, correlation_to_same_tissue_by_tissue$tissue_name %>%
                                unique())) %>%
  ggplot(aes(x = correlation, y = tissue_name)) +
  geom_boxplot(outlier.size=0.3) + xlim(0.5,1) +
  theme_bw() +
  xlab ("Spearman's roh") +
  theme(axis.title.y = element_blank())
library(patchwork)
p1 / p2 / p3 / p4 + plot_layout(heights = c(1, 1, 1, 15))
ggsave( file.path(save_directory, 'spearman_sample.pdf'), height = 14, width = 4)

p3 <- correlation_to_same_tissue %>%
  ggplot(aes(correlation)) +
  geom_histogram(bins = 100) + xlim(0.5,1) +
  theme_classic() +
  theme(panel.background = element_rect("gray90"),
        axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        #axis.title.x = element_blank()
  ) +
  xlab ("Spearman's roh") +
  scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
  ggtitle("Same tissue")

p1 / p2 / p3
ggsave( file.path(save_directory, 'spearman_sample_s.pdf'), height = 4, width = 4)

correlation_to_same_tissue_by_tissue %>%
  mutate(tissue_name = factor(tissue_name, correlation_to_same_tissue_by_tissue$tissue_name %>%
                                unique() %>% rev())) %>%
  ggplot(aes(y = correlation, x = tissue_name, fill = tissue_name)) +
  geom_boxplot(outlier.size=0.3) + ylim(0.8,1) +
  theme_bw() +
  scale_fill_manual(values = adata$obs %>%
                      mutate(tissue_name = stringr::str_to_sentence(tissue_name)) %>%
                      select(tissue_name, organ_color) %>% unique() %>% deframe())+
  ylab ("Spearman's roh") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none")

ggsave( file.path(save_directory, 'spearman_sample_vertical.pdf'), height = 4, width = 10)


#
# correlation_to_same_ctissue_by_ctissue %>%
#   mutate(consensus_tissue_name = factor(consensus_tissue_name, correlation_to_same_ctissue_by_ctissue$consensus_tissue_name %>%
#                                           unique())) %>%
#   ggplot(aes(x = correlation, y = consensus_tissue_name, fill = consensus_tissue_name)) +
#   geom_boxplot(outlier.size=0.3) + xlim(0.5,1) +
#   theme_bw() +
#   scale_fill_manual(values = adata$obs %>%
#                       mutate(consensus_tissue_name = stringr::str_to_sentence(consensus_tissue_name)) %>%
#                       select(consensus_tissue_name, organ_color) %>% unique() %>% deframe())+
#   xlab ("Spearman's roh") +
#   theme(axis.title.y = element_blank(),
#         legend.position = "none")
# ggsave( file.path(save_directory, 'spearman_sample_consensus.pdf'), height = 6, width = 5)
#


