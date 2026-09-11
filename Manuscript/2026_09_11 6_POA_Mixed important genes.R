library(Seurat)

obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
obj = subset(obj, final_clusters==6)

obj= FindSubCluster(obj, '6', graph.name = 'harmony.wsnn')
Idents(obj) ='sub.cluster'
DotPlot(obj, c(
  'drd3',
  'nmbr',
  'npy7r',
  'tacr3a',
  'cckb',
  'ar',
  'esr2b',
  'pgr'
))+
  coord_flip()
›