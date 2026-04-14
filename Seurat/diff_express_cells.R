library(Seurat)
library(ggplot2)
library(dplyr)

path_to_f<-"/home/forbes.ai/capstone_dir/Seurat/subset_reclustered_relaxed.rds"
path_to_f_st<-"/home/forbes.ai/capstone_dir/Seurat/subset_reclustered_stricter.rds"
sub_obj <- readRDS(path_to_f)
sub_obj_st <- readRDS(path_to_f_st)

#annotation
DefaultAssay(sub_obj) <- "SCT"
sub_obj <- PrepSCTFindMarkers(sub_obj)
Idents(sub_obj) <- "seurat_clusters"

DefaultAssay(sub_obj_st) <- "SCT"
sub_obj_st <- PrepSCTFindMarkers(sub_obj_st)
Idents(sub_obj_st) <- "seurat_clusters"

##Used AI to help with crashing errors
library(future)
plan(sequential)

# Restrict features to variable genes for speed/memory (highly recommended)
feat_test <- VariableFeatures(sub_obj)
if (length(feat_test) == 0) {
  # fall back: compute var features on SCT if missing
  sub_obj <- FindVariableFeatures(sub_obj, assay = "SCT", nfeatures = 3000)
  feat_test <- VariableFeatures(sub_obj)
}

feat_test_st <- VariableFeatures(sub_obj_st)
if (length(feat_test_st) == 0) {
  # fall back: compute var features on SCT if missing
  sub_obj_st <- FindVariableFeatures(sub_obj_st, assay = "SCT", nfeatures = 3000)
  feat_test_st <- VariableFeatures(sub_obj_st)
}
###Stop of AI Use

#Find differentially expressed features
sub_obj.markers <- FindAllMarkers(
   sub_obj, 
   only.pos = TRUE,
   logfc.threshold = 0.25,
   test.use = "wilcox_limma" 
)

sub_obj_st.markers <- FindAllMarkers(
   sub_obj_st,
   only.pos = TRUE,
   logfc.threshold = 0.25,
   test.use = "wilcox_limma"
)

#Added this line because program was having trouble finding clusters column name, 
#Sometimes it can use a name
cluster_col <- intersect(c("cluster", "group", "ident"), colnames(sub_obj.markers))[1]
cluster_col_st <- intersect(c("cluster", "group", "ident"), colnames(sub_obj_st.markers))[1]

#Groups the cells by cluster, 
#Filters so only genes with sig diff between clusters are kept 
#Keeps only top 10 genes per cluster
top_10<-sub_obj.markers %>%
    group_by(.data[[cluster_col]]) %>%
    dplyr::filter(avg_log2FC > 1)  %>%
    slice_head(n = 10) %>%
    ungroup()

top_10_st<-sub_obj_st.markers %>%
    group_by(.data[[cluster_col_st]]) %>%
    dplyr::filter(avg_log2FC > 1)  %>%
    slice_head(n = 10) %>%
    ungroup()

saveRDS(top_10, "diff_exp_cluster_top_10_relaxed.rds")
saveRDS(sub_obj, "diff_exp_reclustered_relaxed.rds")
saveRDS(top_10_st, "diff_exp_cluster_top_10_stricter.rds")
saveRDS(sub_obj_st, "diff_exp_reclustered_stricter.rds")

tg_cluster <- top_10 %>%
    group_by(.data[[cluster_col]]) %>%
    slice_head(n = 1) %>% 
    pull(gene)

tg_cluster_st <- top_10_st %>%
    group_by(.data[[cluster_col_st]]) %>%
    slice_head(n = 1) %>%
    pull(gene)
#Write genses to csv file
write.csv(top_10, "cluster_top_markers_relaxed.csv", row.names = FALSE)
write.csv(top_10_st, "cluster_top_markers_stricter.csv", row.names = FALSE)

#Graphs for visualization
heat_map <- DoHeatmap(sub_obj, features = top_10$gene) + 
  theme(
    axis.text.y = element_text(size = 5), 
    axis.text.x = element_text(size = 5)  
  )

heat_map

feat_plot<-FeaturePlot(sub_obj, features = tg_cluster)
feat_plot

met_genes <- c("CDH1", "CDH2", "EPCAM", "FN1", "LOX", "MMP2") 
#("OCLN", "PRSS8", "VCAN", "VIM", "ZEB1")

met_cluster <- c("ADAMTS9", "CXCL11", "TNFAIP6", "QRFPR", "ITGA2",
"AP1S2","ASPSCR1","BIN1","CALU","CAPRIN2","CD6","CDH11","CDH2","CDH3","CDH1",
"CDK1","CEP170","CHN1","CHST3","COL4A1","COL4A2","COL6A3","COL7A1","CREB3L1",
"CXCL2","DAB2","DDR1","DNAJB4","DOCK4","DSC2","EML1","EMP3","EPSTI1","FBN1",
"FBN2","FN1","FLOT1","GBP6","GNG11","HAS2","HMOX1","ICAM1","IFIH1","IGFBP7",
"IL4R","ITGB4","JAG2","KLF4","KRT17","LAPTM4B","LMBR1L","LOX","LRRC32","LTBP2",
"MAP1B","MAP7","MMP2","MYO1B","NID2","NLGN2","NOTCH4","NUP62","P2RX5","PCOLCE",
"PCOLCE2","PDIA4","PKP2","PLIN2","PLLP","PML","PMP22","PRR5L","PTGFR","PTPN6",
"RAB11FIP4","RAB20","RBPMS2","RHOB","RRM2B","SDC1","SERPINB1","SGK3","SLC39A3",
"SLCO2A1","SORL1","SPARC","ST14","TAGLN","TIMP2","TNC","TNFSF10","TRIM29",
"TSPAN5","TUBA1A","UBAC2","VEGFC","VIM","WNT5A","ZC3H12A","ZNF697"
) #second line is start of high and cluster 2


feat_plot2<-FeaturePlot(sub_obj, features = met_genes)
feat_plot3<-FeaturePlot(sub_obj, features = met_cluster)

feat_plot2

feat_plot3

epithelial_markers<-c("CDH1", "DSP", "TJP1")
mesenchymal_markers<-c("VIM", "CDH2", "FOXC2", "SNAI1", "SNAI2", "TWIST1", "GSC", "FN1", "MMP2")
# "MMP3", "MMP9", "SOX10")


feat_plot_ep<-FeaturePlot(sub_obj, features = epithelial_markers)

feat_plot_ep

feat_plot_me<-FeaturePlot(sub_obj, features = mesenchymal_markers)

feat_plot_me


heat_map <- DoHeatmap(sub_obj_st, features = top_10_st$gene) +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5)
  )

heat_map

feat_plot2<-FeaturePlot(sub_obj_st, features = met_genes)
feat_plot3<-FeaturePlot(sub_obj_st, features = met_cluster)

feat_plot2

feat_plot3


feat_plot_ep<-FeaturePlot(sub_obj_st, features = epithelial_markers)

feat_plot_ep

feat_plot_me<-FeaturePlot(sub_obj_st, features = mesenchymal_markers)

feat_plot_me
