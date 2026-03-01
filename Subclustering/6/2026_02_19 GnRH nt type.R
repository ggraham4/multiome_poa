library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
  clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)



obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_6 = FindSubCluster(obj, 6, graph.name='harmony.wsnn')
Idents(sub_6) <- 'sub.cluster'
sub_6 = subset(sub_6, final_clusters ==6)
sub_6$sub.cluster = factor(sub_6$sub.cluster, levels = c(paste0('6_',0:3)))
sub_6$Status = factor(sub_6$Status, levels = c('NRM','M',"D",'E','NF','F'))

### split into gnrh + and gnrh -####
sub_6$gnrh =ifelse(sub_6@assays$RNA$data['LOC111571064',], 'GnRH', 'Not')

Idents(sub_6) <- 'gnrh'

DotPlot(sub_6, 
        features = c('elavl3',
                      'LOC111588076',
                    'gad2',# vgat
                     'LOC111584103',
                     'slc17a6b', # gad1b
                     'slc17a7a'
                     ), scale = F)+
  coord_flip()


