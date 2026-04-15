BiocManager::install(version = "3.22", ask = FALSE)

library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(celldex)

#pick one of the subset files
#path_f<-"/data/subset_reclustered_stricter.rds"
path_f<-"/data/subset_reclustered_relaxed.rds"

sub_obj<-readRDS(path_f)

library(SingleCellExperiment)

mat <- LayerData(sub_obj, assay = "SCT", layer = "data")

sce <- SingleCellExperiment(mat)

colData(sce)$orig.ident  <- sub_obj$orig.ident
colData(sce)$cluster <- sub_obj$seurat_clusters

ref.set <- celldex::HumanPrimaryCellAtlasData()

pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.fine, assay.type.test=1)

lbls.keep <- table(pred.cnts$labels)>10
sub_obj$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')
DimPlot(sub_obj, reduction='umap', group.by='SingleR.labels')

table_of_values<-table(Cluster = sub_obj$seurat_clusters, Label = sub_obj$SingleR.labels)
table_of_values
