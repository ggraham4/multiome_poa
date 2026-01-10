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
  library(SeuratWrappers)
  }

obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds") 
# need this obj because it has the pca

#find rgc subclusters
Idents(obj) = 'res0.8_50nn_40PC_45LSI'
obj <- FindSubCluster(obj, 1, 'harmony.wsnn')
DimPlot(obj, group.by = 'sub.cluster', reduction = 'harmony_wnn.umap')

###subset to only radial glia and neurons
subset_obj <- subset(obj, 
                     #oligos
                     res0.8_50nn_40PC_45LSI!=2&
                     #microglia
                    res0.8_50nn_40PC_45LSI!=11&
                    #opcs
                    res0.8_50nn_40PC_45LSI!=13&
                    #dividing glia
                    res0.8_50nn_40PC_45LSI!=26&
                    #leuko
                    res0.8_50nn_40PC_45LSI!=20&
                    #ependymal
                    res0.8_50nn_40PC_45LSI!=15
                    
                    )


# prepare meta data
gene_meta_data <- data.frame(row.names = rownames(subset_obj@assays$RNA$counts),
                                                  gene_short_name=
                                                    rownames(subset_obj@assays$RNA$counts)
)
cell_metadata <- subset_obj@meta.data


#convert to cds
cds <- as.cell_data_set(subset_obj)
#add in reductions
subset_obj = FindVariableFeatures(subset_obj)

subset_obj = RunPCA(subset_obj, dim = 50, verbose = TRUE, assay = "RNA",               
 features = VariableFeatures(object = subset_obj), reduction.name = "pca_monocle",
 reduction.key = "pca_monocle_")

ElbowPlot(subset_obj, reduction = 'pca_monocle')

reducedDims(cds)$PCA <- Embeddings(subset_obj, "pca_monocle")
reducedDims(cds)$UMAP <- Embeddings(subset_obj, "harmony_wnn.umap")

# Process the dataset with monocle3 functions
cds <- preprocess_cds(cds, num_dim = 5) # based on elbow
cds <- cluster_cells(cds)

# Set up cluster information
clusters <- subset_obj$sub.cluster
cell_types <-subset_obj$sub.cluster

# Ensure names are set correctly
names(clusters) <- colnames(subset_obj)
names(cell_types) <- colnames(subset_obj)

# Assign clusters - make sure we're using the monocle3 structure correctly
colData(cds)$seurat_clusters <- clusters
colData(cds)$cell_type <- cell_types

# Properly assign partitions using colData instead of the older syntax
colData(cds)$partition <- as.factor(cell_types)

# Learn the trajectory graph with optimized parameters
set.seed(0)  # Set seed for reproducibility
cds <- learn_graph(cds, 
                 use_partition = FALSE,
                 close_loop = FALSE,
                 learn_graph_control = list(
                   minimal_branch_len = 15,
                   prune_graph = TRUE,
                   geodesic_distance_ratio = 0.333
                 ))

get_earliest_principal_node <- function(cds=cds, var=c('1_1')){
  cell_ids <- which(colData(cds)[, "sub.cluster"] == var)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  return(root_pr_nodes)
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds), reduction_method = 'UMAP')

plot_cells(cds, 
           color_cells_by = 'res0.8_50nn_40PC_45LSI',
           label_branch_points  = F,
           label_leaves = F,
           label_groups_by_cluster = T,
            group_label_size = 5
           )
plot_cells(cds, color_cells_by = 'pseudotime')

#plot_cells(cds, color_cells_by = 'pseudotime', reduction_method = 'PCA')
subset_obj@meta.data$Status = factor(subset_obj@meta.data$Status , levels = c('NRM',
                                                                              'M',
                                                                              'D',
                                                                              'E',
                                                                              'NF',
                                                                              "F"))
subset_obj$pseudotime = pseudotime(cds)

plot_dat =subset_obj@meta.data%>%
  group_by(individual, sub.cluster, Status)%>%
  summarize(mean_pseudotime = mean(pseudotime))

library(forcats)
ggplot(plot_dat, aes(x = fct_reorder(sub.cluster,  mean_pseudotime), y = (mean_pseudotime)))+
  geom_boxplot()+
  geom_point()

ggplot(subset(plot_dat, Status %in% c('M','D','F')), aes(x = fct_reorder(sub.cluster,  mean_pseudotime), y = (mean_pseudotime), color= Status))+
  geom_boxplot()

ggplot(subset(plot_dat, sub.cluster ==6), aes(x = Status, y = (mean_pseudotime), color= Status))+
  geom_boxplot()+
  geom_point()

ggplot(subset(plot_dat, sub.cluster ==5), aes(x = Status, y = (mean_pseudotime), color= Status))+
  geom_boxplot()+
  geom_point()# could be something

ggplot(subset(plot_dat, sub.cluster ==24), aes(x = Status, y = (mean_pseudotime), color= Status))+
  geom_boxplot()+
  geom_point()# surely

ggplot(subset(plot_dat, sub.cluster ==19), aes(x = Status, y = (mean_pseudotime), color= Status))+
  geom_boxplot()+
  geom_point() # definitely

ggplot(subset_obj@meta.data, aes(x = fct_reorder(sub.cluster,  pseudotime), y = (pseudotime)))+
  geom_boxplot()+
  geom_point()

ggplot(subset(subset_obj@meta.data, Status %in% c('M','D','E','F') ), aes(x = fct_reorder(sub.cluster,  pseudotime), y = (pseudotime),color = Status))+
  geom_boxplot()

ggplot(subset(subset_obj@meta.data,res0.8_50nn_40PC_45LSI==9 & Status %in% c('M','D','E','F') ), aes(x = Status, y = (pseudotime),color = Status))+
  geom_boxplot()+
  geom_point()

ggplot(subset(subset_obj@meta.data,res0.8_50nn_40PC_45LSI==5 & Status %in% c('M','D','E','F') ), aes(x = Status, y = (pseudotime),color = Status))+
  geom_boxplot()+
  geom_point()


