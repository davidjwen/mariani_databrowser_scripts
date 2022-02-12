# library(Seurat)
# library(patchwork)
library(monocle3)
library(tidyverse)
# library(SingleR)

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
names(comb.cds@rowRanges@elementMetadata@listData)[1] = "gene_short_name"
comb.cds = preprocess_cds(comb.cds, num_dim = 25)
comb.cds = reduce_dimension(comb.cds, max_components = 2)
comb.cds = cluster_cells(comb.cds, resolution=1e-5)

# rename the sample numbers
colData(comb.cds)$sample.class = factor(ifelse(colData(comb.cds)$sample == 1, "control", "sox9ko"))


# see what cell types are each cluster
top_cluster_markers = top_markers(comb.cds, 
                                  group_cells_by = "cluster",
                                  cores=10)

tcm <- top_cluster_markers %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  arrange(as.numeric(cell_group), pseudo_R2, specificity)

# manual cell type annotation
colData(comb.cds)$assigned.cell.type = clusters(comb.cds) %>%
  recode(
    "1" = "Osteoblasts and Connective Tissue Cells",
    "2" = "Connective Tissue Cells",
    "3" = "Macrophages",
    "4" = "Osteoclasts",
    "5" = "Muscle Cells 1",
    "6" = "Erythrocytes",
    "7" = "Connective Tissue Cells 2", 
    "8" = "Myeloid Cells",
    "9" = "T-cells",
    "10" = "Endothelial Cells",
    "11" = "Muscle/Connective Tissue Cells",
    "12" = "White Blood Cells",
    "13" = "Keratinocytes",
    "14" = "Muscle Cells 2"
  )

plot_cells(comb.cds, color_cells_by  = "assigned.cell.type") + 
  theme(legend.position="right")

UMAPdims = comb.cds@int_colData$reducedDims$UMAP
UMAPdims = data.frame(UMAP_1 = UMAPdims[, 1], UMAP_2 = UMAPdims[, 2])
UMAPdims$dataset = colData(comb.cds)$sample.class


######################
# WRITE OUT DATA SET #
######################

umap_data = as_tibble(cbind(colData(comb.cds), UMAPdims)) %>% 
  select(!sample & !dataset & !n.umi & !Size_Factor) %>%
  mutate(cell = rownames(colData(comb.cds)),
         cluster = paste0("(", clusters(comb.cds), ") ", assigned.cell.type),
         cluster_num = clusters(comb.cds))

# rename the rownames
# they are the same
exprs_data = assays(comb.cds)$counts
rownames(exprs_data) = comb.cds@rowRanges@elementMetadata@listData$gene_short_name

markers = tcm %>%
  arrange(as.numeric(cell_group), (-pseudo_R2)) %>%
  select(gene_short_name, cell_group) 

write.csv(umap_data, "~/databrowser_deployment/appdata/monocle_clusters.csv", row.names=F)
write.csv(markers, "~/databrowser_deployment/appdata/monocle_markers.csv", row.names=F)
saveRDS(exprs_data, "~/databrowser_deployment/appdata/monocle_rna_counts.RDS")

# the RDS object is smaller
# Matrix::writeMM(exprs_data, "~/databrowser_deployment/appdata/monocle_rna_counts.mm")
