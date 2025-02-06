#Trying out cellchat
setwd('~/Desktop/multiome_poa')
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
  library(hdWGCNA)
  library(WGCNA)
  library(CellChat)


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

#read in seurat object
obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

#OK 2 things, first, coltans code uses ortholog gene names in mus like
#"Gpc1", meaning I need to create an object with changed gene names

#for some reason I need to do this because it doesnt tolerate 0s in cluster labels even if they are characters ??
obj$cluster_labels <-as.character(obj$harmony.wnn_res0.4_clusters)
obj$cluster_labels <- gsub("\\b0\\b", "o", obj$cluster_labels)
obj$cluster_labels <- paste0('Cluster ',obj$cluster_labels)
                             

cellchat <- createCellChat(object = obj, 
                           #meta = meta, 
                           group.by = "cluster_labels")

#need to get code on how to convert to mmusculus

