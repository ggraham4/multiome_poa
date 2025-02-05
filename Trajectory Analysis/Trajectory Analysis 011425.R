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
  library(CytoTRACE)
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
  library(BiocManager)
  library(monocle3)
  library(SingleCellExperiment)

}
#install BiocManager dependencies
#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                       'terra', 'ggrastr'))

#install monocle3
#library("devtools")
#devtools::install_github('cole-trapnell-lab/monocle3')

#load in object
obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

###subset to only radial glia and neurons
subset_obj <- subset(obj, 
                     #oligos
                     harmony.wnn_res0.4_clusters!=4&
                     #fibros
                     harmony.wnn_res0.4_clusters!=29&
                     #microglia
                    harmony.wnn_res0.4_clusters!=14&
                    #opcs
                    harmony.wnn_res0.4_clusters!=18&
                    #Leukocytes
                    harmony.wnn_res0.4_clusters!=26&
                    #ependymal
                    harmony.wnn_res0.4_clusters!=22
                    
                    )
DotPlot(object = subset_obj, 
        group.by = "harmony.wnn_res0.4_clusters", 
        features = c('gad2','LOC111584103','slc17a6b','elavl3'),
        cols = c("#D2B4DE", "#8E44AD", "#6C3483")
) + 
  coord_flip()
###Idk why elavl3 is so low in 28 but I think it is still neuronal

gene_meta_data <- data.frame(row.names = rownames(subset_obj@assays$RNA$counts),
                                                  gene_short_name=
                                                    rownames(subset_obj@assays$RNA$counts)
)

gene_meta_data$cell_type <- ifelse(gene_meta_data$gene_short_name == 'LOC111577263', 'radial glia', NA)
                                                  
cds <- new_cell_data_set(subset_obj@assays$RNA$counts,
                         cell_metadata = subset_obj@meta.data,
                         gene_metadata =gene_meta_data)
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "batch")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)
                    
## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
cds <- order_cells(cds, reduction = 'UMAP')

plot_cells(cds)
                    
plot_cells(cds, genes=c("LOC111577263"), show_trajectory_graph= F)
                  

plot_cells(cds, genes=c("hmx3a",'hmx2','ar'), show_trajectory_graph= F)

###find marker genes
marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=7)

top_specific_markers <- marker_test_res %>%
                            filter(fraction_expressing >= 0.10) %>%
                            group_by(cell_group) %>%
                            top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
top_specific_marker_ids <- append(top_specific_marker_ids, c('ar', 'hmx2','hmx3a'))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)

FeaturePlot(obj, 'gpm6bb')

FindMarkers(obj, 19)

FeaturePlot(obj, 'galr1a')

FeaturePlot(obj, 'kcnh1b')

plot_cells(cds, color_cells_by="pseudotime", show_trajectory_graph= T)
plot_cells(cds, color_cells_by="pseudotime", show_trajectory_graph= F)

FindMarkers(obj, 27)

DotPlot(object = subset_obj, 
        group.by = "harmony.wnn_res0.4_clusters", 
        features = c('kcnh1b')
) + 
  coord_flip()

cds_sub <- choose_graph_segments(cds)

plot_cells(cds_sub, color_cells_by="pseudotime", show_trajectory_graph= T)

cds2 <- preprocess_cds(cds, num_dim = 50)
cds2 <- align_cds(cds2, alignment_group = "harmony.wnn_res0.4_clusters")
cds2 <- reduce_dimension(cds2)

plot_cells(cds2, label_groups_by_cluster=T,  color_cells_by = "harmony.wnn_res0.4_clusters")

clust_19_genes <- c("ar-1",
                    "hmx2",
                    "hmx3a",
                    "esr2b",
                    "ush2a")

plot_cells(cds2,
           genes=clust_19_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

cds2 <- cluster_cells(cds2)
plot_cells(cds2, color_cells_by = "partition")

cds2<- learn_graph(cds2)
plot_cells(cds2,
           color_cells_by = "harmony.wnn_res0.4_clusters",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

cds2 <- order_cells(cds2)
plot_cells(cds2, color_cells_by = 'pseudotime', show_trajectory_graph = T)


# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds=cds2, var=c(0,1,19,15,17,23,25,27,28)){
  cell_ids <- which(colData(cds)[, "harmony.wnn_res0.4_clusters"] %in% var)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds2 <- order_cells(cds2, root_pr_nodes=get_earliest_principal_node(cds2))

plot_cells(cds2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

cds3 <- order_cells(cds2, root_pr_nodes=get_earliest_principal_node(cds2, var = c(0:31)))

plot_cells(cds3,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

cds4 <- cds3
cds4 <- order_cells(cds4)
plot_cells(cds4,
           color_cells_by = 'pseudotime')
plot_cells(cds4, genes = 'LOC111577263', show_trajectory_graph = F)


