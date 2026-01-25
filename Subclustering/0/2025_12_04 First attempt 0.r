library(Seurat)
library(CytoTRACE)
library(lme4)
library(tidyverse)
library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")

a <- listAttributes(ensembl)
biomart_basic <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'entrezgene_description'))


namer = function(gene){
  gene_name = biomart_basic$entrezgene_description[biomart_basic$entrezgene_accession==gene]
  return(gene_name[1])
  
}


obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

marks_0 = FindMarkers(obj, 0)
marks_0$name = sapply(FUN=namer, X=rownames(marks_0))

findsub_0 = FindSubCluster(obj, 0, graph.name = 'harmony.wsnn')%>%
  subset(final_clusters ==0)

DimPlot(findsub_0, reduction = 'harmony_wnn.umap', group = 'sub.cluster')

findsub_0$cyto = CytoTRACE(findsub_0@assays$RNA$data%>%as.matrix())$CytoTRACE

findsub_0@meta.data$Status = factor(findsub_0@meta.data$Status, levels = c('NRM','M','D','E','NF','F'))
dat = findsub_0@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(mean_cyto = mean(cyto),
            se_cyto =sd(cyto)/sqrt(n()))

ggplot(dat, aes(x = sub.cluster, y = mean_cyto, color = Status))+
  geom_boxplot()

allcell_findsub_0 = findsub_0@meta.data%>%
  group_by(individual, Status)%>%
  summarize(total_cells = n())

subcell_findsub_0 = findsub_0@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(cells = n())

prop_cells_findsub0 = subcell_findsub_0%>%
  right_join(allcell_findsub_0, by = 'individual')

prop_cells_findsub0%>%glimpse()

ggplot(prop_cells_findsub0, aes(x = sub.cluster, y = cells/total_cells, color = Status.x))+
  geom_boxplot()

findsub_0_markers = FindAllMarkers(findsub_0, group.by = 'sub.cluster', only.pos = T)

top10_markers <- findsub_0_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 3)

DotPlot(findsub_0, features = unique(top10_markers$gene), group.by = 'sub.cluster') + 
  coord_flip()

DimPlot(obj)

DimPlot(findsub_0, group.by ='sub.cluster')

very_strong_markers_findsub_0 = subset(findsub_0_markers, pct.1>0.9)

# Are the syb



