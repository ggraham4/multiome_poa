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

subset_obj <- FindSubCluster(subset_obj, 2, graph.name ='harmony.wsnn')

gene_meta_data <- data.frame(row.names = rownames(subset_obj@assays$RNA$counts),
                                                  gene_short_name=
                                                    rownames(subset_obj@assays$RNA$counts)
)

cds <- new_cell_data_set(subset_obj@assays$RNA$counts,
                         cell_metadata = subset_obj@meta.data,
                         gene_metadata =gene_meta_data)
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100, method = 'PCA')

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "orig.ident")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds, reduction_method = 'UMAP')

#i get very different results depending on which reduction I use I dont like this

## Step 4: Cluster the cells
cds <- cluster_cells(cds)
                    
## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
get_earliest_principal_node <- function(cds=cds, var=c('2_9')){
  cell_ids <- which(colData(cds)[, "sub.cluster"] == var)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  return(root_pr_nodes)
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

cds <- order_cells(cds)

plot_cells(cds, 
           color_cells_by = 'harmony.wnn_res0.4_clusters',
           label_branch_points  = F,
           label_leaves = F,
           label_groups_by_cluster = T,
            group_label_size = 5
           )
plot_cells(cds, color_cells_by = 'pseudotime')
#why does it think 5 is basal
###why cant i replicate anything
