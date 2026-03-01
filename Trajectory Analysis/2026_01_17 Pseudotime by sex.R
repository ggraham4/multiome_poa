#> sex differences in pseudotime
#> here, I want to see if 
#> 1) The sexes differ in cell density along the pseudotime (do they get hung up)
#> 2) do they have differential expression along the pseudotime
#> 3) Are there sex specific branches in the pseudotime

{
  library(emmeans)
  library(forcats)
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

# read in data
obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds") 
# need this obj because it has the pca

#find rgc subclusters - this tells me the source node (1_1 is the most immature I know)
Idents(obj) = 'res0.8_50nn_40PC_45LSI'
obj <- FindSubCluster(obj, 1, 'harmony.wsnn')
DimPlot(obj, group.by = 'sub.cluster', reduction = 'harmony_wnn.umap')

# factorize status to have plots properly ordered
obj@meta.data$Status = factor(obj@meta.data$Status , levels = c('NRM',
                                                                              'M',
                                                                              'D',
                                                                              'E',
                                                                              'NF',
                                                                              "F"))


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


### Learn the trajectory ####
# prepare meta data
gene_meta_data <- data.frame(row.names = rownames(subset_obj@assays$RNA$counts),
                                                  gene_short_name=
                                                    rownames(subset_obj@assays$RNA$counts)
)
cell_metadata <- subset_obj@meta.data


#convert to cds
cds <- as.cell_data_set(subset_obj)
#add in reductions, here I am going to run a new PCA just to be sure I am getting good PCs
#find var features and run PCA
subset_obj = FindVariableFeatures(subset_obj)
subset_obj = RunPCA(subset_obj, dim = 50, verbose = TRUE, assay = "RNA",               
 features = VariableFeatures(object = subset_obj), reduction.name = "pca_monocle",
 reduction.key = "pca_monocle_")


ElbowPlot(subset_obj, reduction = 'pca_monocle')

#add these reductions into the cds object
reducedDims(cds)$PCA <- Embeddings(subset_obj, "pca_monocle")
reducedDims(cds)$UMAP <- Embeddings(subset_obj, "harmony_wnn.umap")

# Process the dataset with monocle3 functions
cds <- preprocess_cds(cds, num_dim = 5) # based on elbow, I think 5 is a good number
# cluster
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

# here, I define a function to set 1_1 as the earliest principal node,
#since my other work has shown this to be a stem cell or at least the most
# immature

get_earliest_principal_node <- function(cds=cds, var=c('1_1')){ 
  cell_ids <- which(colData(cds)[, "sub.cluster"] == var) # pick cell barcodes in 1_1
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ]) # find closest vertex
  
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))] # return nodes
  
  return(root_pr_nodes)
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds), reduction_method = 'UMAP')

# plot to make sure its the same
plot_cells(cds, 
           color_cells_by = 'res0.8_50nn_40PC_45LSI',
           label_branch_points  = F,
           label_leaves = F,
           label_groups_by_cluster = T,
            group_label_size = 5
           )
plot_cells(cds, color_cells_by = 'pseudotime')
# indeed, 1_1 is the node, and it looks like differentiating cells transit through
# 24 and 0. I strongly believe 22 is also some kind of immature cell, but I only have
# my speculation to support that for now

## add to seurat object
subset_obj$pseudotime = pseudotime(cds)%>%as.numeric()

#### ANALYSIS 1: Difference in Sex Density Across Pseudotime ####
# extract meta data
meta <- subset_obj@meta.data %>%
  subset(Status != 'NRM')%>%
  dplyr::select(pseudotime, Status, individual) %>% 
  filter(is.finite(pseudotime))

n_bins <- 50  # Resolution of analysis
max_time <- max(meta$pseudotime)
min_time <- min(meta$pseudotime)

window_width <- (max_time - min_time) / 20 
step_size <- (max_time - min_time) / n_bins

# 2. Define the Rolling Window parameters
n_bins <- 50  # Resolution of your analysis
max_time <- max(meta$pseudotime)
min_time <- min(meta$pseudotime)
# Create overlapping windows (smoother than rigid bins)
window_width <- (max_time - min_time) / 20 
step_size <- (max_time - min_time) / n_bins

# Define center points for the windows
eval_points <- seq(min_time, max_time, by = step_size)

# 3. The "Rolling" Calculation
# For every individual, what % of its total cells fall in Window X?

total_cells = meta%>%
  group_by(individual, Status)%>%
  summarize(n_cells = n())

rolling_calc = function(window){
  # so at each eval point + / - window width/2, calculate proportion
  
  #define boundaries
  window_start = window - (window_width/2)
  window_end = window + (window_width/2)
  
  #extract cells that meet those boundaries
  meta_window = subset(meta, pseudotime >= window_start & pseudotime <= window_end)
  
  #calculate cells in window
  cells_in_window = meta_window%>%
    group_by(individual)%>%
    summarize(cells_in = n())
    
  #calculate proportion, handle 0s
  joint = cells_in_window%>%
    right_join(total_cells, by = 'individual')%>%
    mutate(cells_in = ifelse(is.na(cells_in), 0, cells_in)) %>% 
    mutate(proportion_cells = cells_in/n_cells)
  joint$window = window
  
  return(joint)
}

rolling_proportion = lapply(eval_points, rolling_calc)
rolling_proportion_bound = do.call(rbind, rolling_proportion)

ggplot(rolling_proportion_bound, aes(x = window, y = proportion_cells))+
  geom_point(aes(color = Status))+
  geom_smooth(aes(color = Status), se = F)

ggplot(rolling_proportion_bound, aes(x = window, y = proportion_cells))+
  geom_smooth(aes(color = Status), se = F)
#they all seem VERY similar, but lets do statistics anyway

rolling_stats_out = data.frame()
for(i in eval_points){
  subset_eval= subset(rolling_proportion_bound, window ==i & Status %in% c('M','D','F'))
  
  model_mat = cbind(subset_eval$cells_in, subset_eval$n_cells-subset_eval$cells_in)
  model = glm(model_mat~Status, data = subset_eval, family = 'binomial')
  av_ = anova(model, test = 'Chisq')
  
  pair = pairs(emmeans(model, 'Status'), adjust = 'none')%>%as.data.frame()
  
  newd = data.frame(window = i,
                    av_p = av_$`Pr(>Chi)`[2],
                    m_d_p = pair$p.value[pair$contrast == 'M - D'],
                    m_f_p = pair$p.value[pair$contrast == 'M - F'],
                    d_f_p = pair$p.value[pair$contrast == 'D - F'],
                    m_d_estimate = pair$estimate[pair$contrast == 'M - D'],
                    m_f_estimate = pair$estimate[pair$contrast == 'M - F'],
                    d_f_estimate = pair$estimate[pair$contrast == 'D - F']
                    )
  rolling_stats_out= rbind(rolling_stats_out, newd)
  

}
# p value correction
rolling_stats_out$av_q = p.adjust(rolling_stats_out$av_p, 'fdr', nrow(rolling_stats_out))
rolling_stats_out$signif = ifelse(rolling_stats_out$av_q <0.05, '*',NA)

# plot significant windows
ggplot(rolling_proportion_bound, aes(x = window, y = proportion_cells))+
  geom_point(aes(color = Status))+
  geom_smooth(aes(color = Status), se = F)+
  geom_vline(data = rolling_stats_out, aes(xintercept =window[signif =='*']))


# assign windows to clusters
window_to_cluster = function(window){
  
  #define boundaries
  window_start = window - (window_width/2)
  window_end = window + (window_width/2)
  
  #extract cells that meet those boundaries
  meta_window = subset(meta, pseudotime >= window_start & pseudotime <= window_end)
  meta_cells = rownames(meta_window)
  
  clusters = table(subset_obj@meta.data$sub.cluster[rownames(subset_obj@meta.data)%in% meta_cells])%>%
    as.data.frame()
  sorted_clusters <- clusters[order(clusters$Freq, decreasing = T), ]

    top_cluster = sorted_clusters$Var1[1]
    top_cluster_prop = sorted_clusters$Freq[1]/sum(sorted_clusters$Freq)
  
    second_cluster = sorted_clusters$Var1[2]
    second_cluster_prop = sorted_clusters$Freq[2]/sum(sorted_clusters$Freq)

newd = data.frame(window = window,
                  top_cluster = top_cluster,
                  top_cluster_prop = top_cluster_prop,
                  second_cluster = second_cluster,
                  second_cluster_prop=second_cluster_prop)

return(newd)
  
}

# plot windows
plot_window = function(wind){
  plot_data = subset(rolling_proportion_bound, window == wind)
  
  p = ggplot(plot_data, aes(x = Status, y = proportion_cells))+
    geom_boxplot()+
    geom_point()
  return(p)
  
}

#get significant windows
significant_windows =rolling_stats_out$window[rolling_stats_out$av_q<0.05]

sig_windows_labeled = lapply(significant_windows, window_to_cluster)
sig_windows_labeled_bound = do.call(rbind,sig_windows_labeled )

rolling_stats_out_labeled = rolling_stats_out%>%
  left_join(sig_windows_labeled_bound, by = 'window')

plot_window(significant_windows[1])
plot_window(significant_windows[6])

# differences in cluster 5
plot_window(significant_windows[4])# fewer in females
plot_window(significant_windows[5]) # fewer in females

# cluster 4
plot_window(significant_windows[8]) #fewer in doms 
plot_window(significant_windows[9]) # fewer in doms

# 1_0 (potentially gliogenic)
plot_window(significant_windows[2]) # more in doms, fewer in females
plot_window(significant_windows[3]) #same as above

# I also want to say that what we are saying here is that there is that there 
# are differences at this developmental timepoint

hist(meta$pseudotime)
# so we see differences at early pseudotimepoints (4,5,6,7, and late ish (21, 24, 25))

plot_cells(cds, 
           color_cells_by = 'res0.8_50nn_40PC_45LSI',
           label_branch_points  = F,
           label_leaves = F,
           label_groups_by_cluster = T,
            group_label_size = 5
           )+
plot_cells(cds, color_cells_by = 'pseudotime')
DimPlot(subset_obj, reduction = 'harmony_wnn.umap', label = T)

DimPlot(subset(subset_obj,sub.cluster =='5'), reduction = 'harmony_wnn.umap', label = T)
# 5 is an interesting cluster that kind of spans across 0
DimPlot(subset(subset_obj,sub.cluster =='5' & Status %in% c('M','D','F')), reduction = 'harmony_wnn.umap', label = F, group.by = 'Status')
FeaturePlot(subset(subset_obj,sub.cluster =='5' & Status %in% c('M','D','F')), reduction = 'harmony_wnn.umap', label = F, feature = 'pseudotime')


DimPlot(subset(subset_obj,sub.cluster%in%c('5','0')), reduction = 'harmony_wnn.umap', label = T, alpha = 0.1)

# OK well a lot to chew on, but for now I guess its time to move on

# I think based on the UMAP, slingshot is a better way to do this



