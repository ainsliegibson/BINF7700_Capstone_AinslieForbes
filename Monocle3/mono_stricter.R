library(Seurat)
library(monocle3)
library(ggplot2)

path_to_f<-"/home/forbes.ai/capstone_dir/Seurat/diff_exp_reclustered_stricter.rds"
sub_obj <- readRDS(path_to_f)

expr <- sub_obj[["RNA"]]$counts
cell_metadata <- sub_obj@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(expr),
  row.names = rownames(expr)
)
cds <- new_cell_data_set(
  expr,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds) #added

cds <- align_cds(cds, alignment_group ="orig.ident", residual_model_formula_str = "~ percent.mt + percent.ribo")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds) #added

seurat_cluster_plot<-plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_clusters",
           graph_label_size=1.5, label_cell_groups =FALSE)
seurat_cluster_plot + ggtitle("seurat_cluster_plot")
#I want to look at epithelial and MSC clusters closer


part_cluster_plot<-plot_cells(cds, color_cells_by = "partition",
           graph_label_size=1.5, label_cell_groups =FALSE)
part_cluster_plot + ggtitle("part_cluster_plot")
#partiton 3 & 2 + 1, 2 & 5

cds_p23 <- cds[, partitions(cds) %in% c("2", "3")]
cds_p23 <-preprocess_cds(cds_p23, num_dim = 50)
plot_pc_variance_explained(cds_p23)
cds_p23 <-reduce_dimension(cds_p23)
cds_p23 <- cluster_cells(cds_p23, partition_qval = 1)

cds_p123 <- cds[, partitions(cds) %in% c("1", "2", "3")]
cds_p123 <-preprocess_cds(cds_p123, num_dim = 50)
plot_pc_variance_explained(cds_p123)
cds_p123 <-reduce_dimension(cds_p123)
cds_p123 <- cluster_cells(cds_p123, partition_qval = 1)


part_cluster_plot_2<-plot_cells(cds_p23, color_cells_by = "partition",
           graph_label_size=1.5, label_cell_groups =FALSE)
part_cluster_plot_2 + ggtitle("part_cluster_plot_23")


part_cluster_plot_3<-plot_cells(cds_p123, color_cells_by = "partition",
           graph_label_size=1.5, label_cell_groups =FALSE)
part_cluster_plot_3 + ggtitle("part_cluster_plot_123")

cds <- learn_graph(cds)
cds_p23 <- learn_graph(cds_p23)

mono_cluster_plot<-plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, label_cell_groups =FALSE)
mono_cluster_plot + ggtitle("mono_cluster_plot")

seur_cluster_plot_2<-plot_cells(cds_p23,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, label_cell_groups =FALSE)
seur_cluster_plot_2 + ggtitle("seurat_cluster_plot_2")

time_cluster_plot<-plot_cells(cds,
           color_cells_by = "orig.ident",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
time_cluster_plot + ggtitle("time_cluster_plot")

#root 0 for stricter
root_cells <- colnames(cds)[cds$seurat_clusters == "0"]
cds <- order_cells(cds, root_cells = root_cells)

root_graph_plot<-plot_cells(cds,graph_label_size=1.5, label_cell_groups =FALSE)
root_graph_plot+ ggtitle("root_graph_plot")

puesdotime_cluster_plot<-plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) 
puesdotime_cluster_plot + ggtitle("puesdotime_cluster_plot")


root_cells_2 <- colnames(cds_p23)[cds_p23$seurat_clusters == "0"]
cds_p23 <- order_cells(cds_p23, root_cells = root_cells_2)

root_graph_plot_2<-plot_cells(cds_p23,graph_label_size=1.5, label_cell_groups =TRUE)
root_graph_plot_2+ ggtitle("root_graph_plot_2")

puesdotime_cluster_plot_2<-plot_cells(cds_p23,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
puesdotime_cluster_plot_2 + ggtitle("puesdotime_cluster_plot_2")

q()
