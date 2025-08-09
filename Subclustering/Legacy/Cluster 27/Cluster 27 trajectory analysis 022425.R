{
  library(parallel)
  library(clusterProfiler)
  library(blme)
  library(Seurat)
  library(tidyverse)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(SeuratObject)
  library(Signac)
  library(glmGamPoi)
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)
  library(emmeans)
  library(CytoTRACE)
  library(ggrepel)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(Polychrome)
  library(scCustomize)
  
  P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
  swatch(P40)
  names(P40) <- NULL
  
  mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
  prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
  define_degs_prop<- readRDS('Functions/define_degs_prop.rds')
  mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
  prop_deg_function.rds<- readRDS('Functions/DEG_functions/prop_deg_function.rds')
  define_behavior_degs<- readRDS('Functions/define_behavior_degs')
  clown_go<- readRDS('Functions/clown_go')
  define_degs<- readRDS('Functions/define_degs')
  
}

obj <- readRDS('C:/Users/Gabe/Desktop/RNA Object.rds')


library(monocle3)

### Trajectory
subset_obj <- obj[,!obj$harmony.wnn_res0.4_clusters %in% c(4, 29, 14, 18, 26, 22,2, 30, 15)] #remove nonneuronal clusters

DimPlot(subset_obj, group.by = 'harmony.wnn_res0.4_clusters', label = T)

### First, I want to see what it considers basal
subset_obj$cell_type <- ifelse(subset_obj$harmony.wnn_res0.4_clusters==27, 'immature_27',NA)

##Learn graph
gene_meta_data <- data.frame(row.names = rownames(subset_obj@assays$RNA$data),
                             gene_short_name=
                               rownames(subset_obj@assays$RNA$data)
)

cds <- new_cell_data_set(subset_obj@assays$RNA$data,
                         cell_metadata = subset_obj@meta.data,
                         gene_metadata =gene_meta_data)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "orig.ident") 

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
get_earliest_principal_node <- function(cds=cds, var=c(unique(subset_obj$harmony.wnn_res0.4_clusters))){
  cell_ids <- which(colData(cds)[, "harmony.wnn_res0.4_clusters"] %in% var)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))


plot_cells(cds, color_cells_by = 'pseudotime')
plot_cells(cds, color_cells_by = 'pseudotime', show_trajectory_graph = F)
plot_cells(cds, color_cells_by= 'cell_type', show_trajectory_graph = F)
plot_cells(cds, color_cells_by= 'harmony.wnn_res0.4_clusters', show_trajectory_graph = F)
plot_cells(cds, color_cells_by= 'harmony.wnn_res0.4_clusters', show_trajectory_graph = T)

### unbiased it thinks 7 is basal and doesnt lead anywhere
###forcing 27 to be basal
get_earliest_principal_node <- function(cds=cds, var=c(27)){
  cell_ids <- which(colData(cds)[, "harmony.wnn_res0.4_clusters"] %in% var)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds2 <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
plot_cells(cds2, color_cells_by = 'pseudotime')
plot_cells(cds2, color_cells_by = 'pseudotime', show_trajectory_graph = F)
plot_cells(cds2, color_cells_by= 'cell_type', show_trajectory_graph = F)
plot_cells(cds2, color_cells_by= 'harmony.wnn_res0.4_clusters', show_trajectory_graph = F)
plot_cells(cds2, color_cells_by= 'harmony.wnn_res0.4_clusters', show_trajectory_graph = T)

#alright it isnt able to put it on a trajectory so that sucks...
