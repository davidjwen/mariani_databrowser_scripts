library(monocle3)
library(tidyverse)
library(patchwork)

setwd("~/databrowser/appdata")
rm(list=ls())


# quick wrapper to read in data
create_cds = function(path){
  load_mm_data(
    mat_path = paste0(path, "matrix.mtx"),
    feature_anno_path = paste0(path, "features.tsv"),
    cell_anno_path = paste0(path, "barcodes.tsv")
  )
}

# load in the data from the data directory
CONTROL_DIR = "../data/max4dpr/control/filtered_feature_bc_matrix/"
control.cds = create_cds(CONTROL_DIR)

KO_DIR = "../data/max4dpr/sox9ko/filtered_feature_bc_matrix/"
ko.cds = create_cds(KO_DIR)

# combine the two data sets to be processed together
comb.cds = combine_cds(list(control.cds, ko.cds))
# note: need a gene_short_name column added in the data in the gene_metadata
# just needs to be renamed

preprocess_data = function(cds){
  names(cds@rowRanges@elementMetadata@listData)[1] = "gene_short_name"
  cds = preprocess_cds(cds, num_dim = 25)
  cds = reduce_dimension(cds, max_components = 2)
  cds = cluster_cells(cds, resolution=1e-5)
  
  return(cds)
}

comb.cds = preprocess_cds(comb.cds)


# grab the cells in cluster 2 that are missing in the Sox9 KO data set
in_polygon = function(x, y, xmin, xmax, ymin, ymax){
  slope = (ymin-ymax)/(xmax-xmin)
  return(between(x, xmin, xmax) & between(y, ymin, ymax) & (y > (slope*x - slope*xmin + ymax)))
}

colData(comb.cds)$sample.class = factor(ifelse(colData(comb.cds)$sample == 1, "control", "sox9ko"))

in_polygon = Vectorize(in_polygon, vectorize.args = c("x","y"))

UMAPdims = comb.cds@int_colData$reducedDims$UMAP
UMAPdims = data.frame(UMAP_1 = UMAPdims[, 1], UMAP_2 = UMAPdims[, 2])
UMAPdims$dataset = colData(comb.cds)$sample.class

sox9ko_population_class = factor(ifelse(
  in_polygon(UMAPdims$UMAP_1, UMAPdims$UMAP_2, 0, 5.0, 4, 7) & colData(comb.cds)$sample.class == "control", 
  "sox9missing", "other"))

colData(comb.cds)$sox9ko_cells = sox9ko_population_class

p1 = plot_cells(comb.cds, color_cells_by = "sample.class") + 
  theme(legend.position="right")
p2 = plot_cells(comb.cds, color_cells_by="sox9ko_cells") + 
  theme(legend.position="right")

p1/p2

top_missing_markers = top_markers(comb.cds, 
                                  group_cells_by = "sox9ko_cells",
                                  cores=10)

tmm <- top_missing_markers %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group)

tmm_ids <- unique(tmm %>% pull(gene_id))

# noticed that Klf4 missing
plot_cells(comb.cds, genes="Klf4")


subset_barcodes = as_tibble(colData(comb.cds)) %>%
  mutate(cell = rownames(colData(comb.cds)), 
                         cluster = clusters(comb.cds)) %>%
  filter(cluster == 2) %>%
  select(cell) %>%
  unlist()

subset.cds = comb.cds[, colnames(comb.cds) %in% subset_barcodes]
subset.cds = preprocess_cds(subset.cds, num_dim = 25)
subset.cds = reduce_dimension(subset.cds, max_components = 2)
subset.cds = cluster_cells(subset.cds, resolution=1e-4)

plot_cells(subset.cds, group_cells_by = "cluster") + 
  theme(legend.position="right")

p1 = plot_cells(subset.cds, color_cells_by = "sample.class") +
  facet_grid(rows = vars(sample.class)) + 
  theme(legend.position = "right")

p2 = plot_cells(subset.cds, genes = c("Myc")) + 
  facet_grid(rows = vars(sample.class))

p1+p2

top_subset_markers = top_markers(subset.cds, 
                                  group_cells_by = "cluster",
                                  cores=10)

tsm <- top_subset_markers %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  arrange(as.numeric(cell_group), pseudo_R2, specificity)

top_missing_markers = top_markers(subset.cds, 
                                  group_cells_by = "sox9ko_cells",
                                  cores=10)

tmm <- top_missing_markers %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  arrange(cell_group, -pseudo_R2)

tmm_ids <- unique(tmm %>% pull(gene_id))

control.cds = preprocess_data(control.cds)
plot_cells(comb.cds, genes=c("Ppbp", "Tyrobp", "Ccl6", "Mmp9", "Atp6v0d2", "Mmp13", "Sat1", "Acp5", "Ctsk"))
