library(Seurat)
library(ggplot2)
library(tidyverse)
library(scCustomize)

multiome <- readRDS("/Users/ggraham/Desktop/optimal_clustering_rna_only.rds")

vis = readRDS("~/Desktop/Visium/vis_dissection.rds")
Idents(vis) = 'seurat_clusters.projected0.4'

marks_vis = FindAllMarkers(vis, only.pos = T)

marks_vis_top = marks_vis%>%
  filter(!str_detect(gene, "LOC"))%>%
  group_by(cluster)%>%
  slice_min(order_by = p_val_adj, n = 2)

marks_mul = FindAllMarkers(multiome, only.pos = T, group.by = 'final_clusters')

marks_mul_top <- marks_mul %>%
  filter(!str_detect(gene, "LOC")) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup()

features = c(
  's100b', # ependymoglia 
    'gfap', #rgc
  'crocc2', # ependymal cell,
    'top2a', #DG
  'mpz', # oligo,
  'cspg4', #opc
  'ptprc', # immune
  'p2ry12', # immune
    'elavl3', # neuron

  
  #neurotransmitters
    'slc17a7a',
  'slc17a6b',
  'gad1a',
  'gad2',
  
  #regions
    'sp9', # subpallial
  'six3a',
  'lhx9',
  
  'tshz1', # dm
  'eomesb',
  'prkcbb', # dm-3
  
  'sox5', # dd
  'rftn1a',
  
  'ptk2bb', # dln # not sequenced
    'eomesa',
  'npy8br',
  
  'prkcq',
  'LOC111578282',
  
  'c1ql3b', #dl-g,
  'foxo1a',
  'reln',
  
  'calb2a',
  
  'foxa2', # not sure
  'prdx1',
  'nxph1',
  'slc6a1b',
  
  'hmx2',
  'hmx3a',
  'hmx3b', 
  'gal',# preoptic area

  'npy',
  'sst1.1' # vsst

)

DotPlot(multiome, features, dot.min = 0.1)+
  coord_flip()

Clustered_DotPlot(multiome, features, k = 7)+
  coord_flip()

Clustered_DotPlot(multiome, marks_vis_top$gene,
                    row_label_size = 5)+
  coord_flip()

Clustered_DotPlot(multiome, marks_mul_top$gene)+
  coord_flip()

marks_ranked <- marks_mul %>%
  filter(!str_detect(gene, "LOC")) %>%
  subset(pct.1 > 0.1)%>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val_adj, .by_group = TRUE) %>%
  mutate(rank_in_cluster = row_number()) %>%
  slice_head(n = 3) %>%
  ungroup()

DotPlot(multiome, marks_ranked$gene, dot.min = 0.1
        )+
  coord_flip()

Clustered_DotPlot(multiome, marks_ranked$gene,
                    cluster_ident = F,
                    cluster_feature = F,
        )+
  coord_flip()


write.csv( data.frame( gene = marks_ranked$gene,
                       cluster = marks_ranked$cluster), 'markers.csv')

markers2 =read.csv('markers.csv')

Clustered_DotPlot(multiome, markers2$gene,
                    cluster_ident = F,
                    cluster_feature = F,
        )+
  coord_flip()


