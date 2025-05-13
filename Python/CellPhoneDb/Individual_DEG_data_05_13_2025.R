library(Seurat)
library(SeuratObject)
library(Matrix)

options(Seurat.object.assay.version = "v3")

human_named = readRDS( 'A:/optimal_clustering_05_06_2025/RNA_object_human_names.rds')

human_named$res0.8_50nn_40PC_45LSI <- as.character(human_named$res0.8_50nn_40PC_45LSI)
Idents(human_named) = human_named$res0.8_50nn_40PC_45LSI


for(i in unique(human_named$individual)){
  
  if(i == 'GH'){next}
  print(i)
  subset_obj = subset(human_named, individual == i)
  
  subset_obj@meta.data$Cell = rownames(subset_obj@meta.data)
  df = subset_obj@meta.data[, c('Cell', 'res0.8_50nn_40PC_45LSI')]
  
  meta_path = paste0('A:/cellphonedb_05_12_2025/',i,'/human_named_meta.tsv')
  
  write.table(df, file =meta_path, sep = '\t', quote = F, row.names = F)
  
  
  DEGs <- FindAllMarkers(subset_obj, 
                         verbose = T, 
                         only.pos = T, 
                         random.seed = 1, 
                         logfc.threshold = 0.2, 
                         min.pct = 0.1, 
                         return.thresh = 0.05)
  
  fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC  > 0.1)
  
  # 1st column = cluster; 2nd column = gene 
  fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 
  
  file_path = paste0('A:/cellphonedb_05_12_2025/',i,'/DEGs.tsv')
  write.table(fDEGs, file =file_path , sep = '\t', quote = F, row.names = F)
  
}
