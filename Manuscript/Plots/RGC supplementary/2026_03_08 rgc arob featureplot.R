library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)
mecp = readRDS('Functions/mean_expression_cluster_plot.rds')
mecd = readRDS('Functions/mean_expression_cluster_data.rds')

#read in obj
obj  = readRDS("A:/optimal_clustering_05_06_2025/nemo.orig_harmony.integration_all_testd_clusters.rds")

# subcluster
sub_1 = FindSubCluster(obj,1, graph.name='harmony.wsnn')
Idents(sub_1) <- 'sub.cluster'
sub_1 = subset(sub_1, res0.8_50nn_40PC_45LSI ==1)
sub_1$Status = factor(sub_1$Status, levels = c('NRM','M',"D",'E','NF','F'))

# plot by aroB
arob = FeaturePlot(sub_1, 'LOC111577263', reduction = 'harmony_wnn.umap', pt.size = 0.1, alpha = 1)+
  theme_void()+
  xlim(-5,5)+
  ylim(1,10)+
  theme(title = element_blank())

arob

ggsave(plot = arob,
      file = "UMAP_rgc arob featureplot.tiff",
     device = "tiff",
    units = "in",
   width =3.16 ,
  height = 2.08,      
  path = "Manuscript/Plots/")

