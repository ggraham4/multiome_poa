library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
  clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)



colors = c("red", "#006400", "blue",'#000000', 'purple','gray','brown','orange')

# DEG enrichment? 
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
# cluster 6
all_degs_6_go = clown_go(degs$gene[degs$cluster==6])%>%dotplot()+
  labs(title = '6_POA_2 All DEGs')

degs_6_go_df=clown_go(degs$gene[degs$cluster==6])
d = data.frame(degs_6_go_df@result$Count,
degs_6_go_df@result$Description)
  
ggsave(plot = all_degs_6_go,
       file = "all_degs_6_go.svg",
         device = 'svg',
         units = "in",
         width = 5,
         height = 2.7,
         path = "Manuscript/Plots/Manuscript v1.2.1/6_deg_go")

