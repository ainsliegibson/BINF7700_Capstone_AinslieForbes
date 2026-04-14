library(Seurat)
library(ggplot2)
library(dplyr)

path_to_f<-"~/capstone_dir/axo_data/seurat_files/all_umap.rds"

obj <- readRDS(path_to_f)
obj

#Name of clusters
Idents(obj) <- "seurat_clusters"

#Cells of interst based on UMAP from Edha UMAP
cells_to_subset <- WhichCells(obj, idents = c(1, 3, 4, 5, 6, 9, 10, 11, 16))
sub_obj <- subset(obj, cells = cells_to_subset)

cells_to_subset_strict <- WhichCells(obj, idents = c(1, 3, 4, 10, 11))
sub_obj_st <- subset(obj, cells = cells_to_subset_strict)

# use integrated assay for reclustering
DefaultAssay(sub_obj) <- "integrated"


sub_obj <- RunPCA(sub_obj, verbose = FALSE)
sub_obj_st <- RunPCA(sub_obj_st, verbose = FALSE)

v_dim<-VizDimLoadings(sub_obj, dims = 1:2, reduction = "pca")
d_heat<-DimHeatmap(sub_obj, dims = 1:18, cells = 500, balanced = TRUE)
e_plot<-ElbowPlot(sub_obj)

v_dim_st<-VizDimLoadings(sub_obj_st, dims = 1:2, reduction = "pca")
d_heat_st<-DimHeatmap(sub_obj_st, dims = 1:18, cells = 500, balanced = TRUE)
e_plot_st<-ElbowPlot(sub_obj_st)

#Based on results from the heatmap and elbow map it was decided to run clustering with 12 dimentions
sub_obj <- FindNeighbors(sub_obj, dims = 1:12)
sub_obj <- FindClusters(sub_obj, resolution = 0.4)  
sub_obj <- RunUMAP(sub_obj, dims = 1:12)

saveRDS(sub_obj, "subset_reclustered_relaxed.rds")

sub_obj_st <- FindNeighbors(sub_obj_st, dims = 1:8)
sub_obj_st <- FindClusters(sub_obj_st, resolution = 0.4)
sub_obj_st <- RunUMAP(sub_obj_st, dims = 1:8)

saveRDS(sub_obj_st, "subset_reclustered_stricter.rds")

umap_plot<-DimPlot(sub_obj, label = TRUE, repel = TRUE) 
umap_plot_st<-DimPlot(sub_obj_st, label = TRUE, repel = TRUE)

#output plots

d_heat
ggsave("elbow_plot_rec_relax.png", e_plot)
ggsave("umap_plot_rec_relax.png", umap_plot)

d_heat_st
ggsave("elbow_plot_rec_stricter.png", e_plot_st)
ggsave("umap_plot_rec_stricter.png", umap_plot_st)
