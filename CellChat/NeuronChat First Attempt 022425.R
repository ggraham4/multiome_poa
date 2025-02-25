### Running NeuronChat
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
obj <- readRDS('C:/Users/Gabe/Desktop/RNA object mouse names.rds')
#obj <- subset(obj, harmony.wnn_res0.4_clusters != 29)
#for some reason I need to do this because it doesnt tolerate 0s in cluster labels even if they are characters ??
obj$cluster_labels <-as.character(obj$harmony.wnn_res0.4_clusters)
obj$cluster_labels <- gsub("\\b0\\b", "o", obj$cluster_labels)
obj$cluster_labels <- paste0(obj$cluster_labels)


library(NeuronChat)

normalized_count_mtx <- obj@assays$RNA$data

# creat NeuronChat object 
x <- createNeuronChat(normalized_count_mtx,DB='mouse',group.by = obj$harmony.wnn_res0.4_clusters) # use DB='human' for human data

# calculation of communication networks  
x <- run_NeuronChat(x,M=100) #really im not even allowed to subset my fucking object fuck this pacakge

# aggregating networks over interaction pairs
net_aggregated_x <- net_aggregation(x@net,method = 'weight')

# visualization
netVisual_circle_neuron(net_aggregated_x)

group = as.character(obj$harmony.wnn_res0.4_clusters)
group = ifelse(group %in% c('3','6','7','8','9','10','11','12','13','15','16','17','23','24','27'),
                      'Glutamatergic',
                      ifelse(group %in% c('1','5','20','25','31'),
                             'GABAergic',
                             ifelse(group %in% c('1','19'), 
                                    'Mixed',
                                    'Non-Neuronal')
                      )
)
names(group) = as.character(obj$harmony.wnn_res0.4_clusters)




#I dont understand what to do with this grouping argument
###--- THis grouping variable is fucking ridiculous, it wants me to do the labeling by hand??? That defeats the whole point of 
#fucking using neuronchat ####---
unique(names(x@net))
netVisual_chord_neuron(x,interaction_use='Tac2_Tacr3',group=,lab.cex = 1)

####

g1 <- rankNet_Neuron(x,slot.name = "net",measure = c("weight"),mode='single',font.size = 5) 
g2 <- rankNet_Neuron(x,slot.name = "net",measure = c("count"),mode='single',font.size = 5)
g1+g2

#error
x<- identifyCommunicationPatterns_Neuron(x, slot.name = "net", pattern = c("outgoing"), k=4,height = 18)
x<- identifyCommunicationPatterns_Neuron(x, slot.name = "net", pattern = c("incoming"), k=4,height = 18)
#error

library(ggalluvial)
netAnalysis_river_Neuron(x,slot.name = "net", pattern = c("outgoing"),font.size = 2.5,cutoff.1 = 0.5,cutoff.2=0.5)
#this package sucks so fucking bad
netAnalysis_river_Neuron(x,slot.name = "net", pattern = c("incoming"),font.size = 2.5,cutoff.1 = 0.5,cutoff.2=0.5)

x <- computeNetSimilarity_Neuron(x,type='functional')
x  <- CellChat::netEmbedding(x, slot.name = "net_analysis", type = "functional")
#> Manifold learning of the signaling networks for a single dataset
x <- CellChat::netClustering(x, type='functional',slot.name = "net_analysis",k=5)
#> Classification learning of the signaling networks for a single dataset
netVisual_embedding_Neuron(x, type = "functional", label.size = 5,pathway.remove.show = F)

netVisual_embeddingZoomIn_Neuron(x, type = "functional", nCol = 2,label.size = 3)

