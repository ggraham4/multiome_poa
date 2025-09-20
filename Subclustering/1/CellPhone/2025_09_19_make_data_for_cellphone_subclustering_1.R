#make data 6 subclusters for cellphone

library(Seurat)
library(SeuratObject)
library(Matrix)

options(Seurat.object.assay.version = "v3")

obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

#subcluster
obj_subclustered_1 <- FindSubCluster(obj, 1, 'harmony.wsnn')

human_named = readRDS( 'A:/optimal_clustering_05_06_2025/RNA_object_human_names.rds')

human_named$sub.cluster <- as.character(obj_subclustered_1$sub.cluster)
Idents(human_named) = human_named$sub.cluster

# Save normalised counts - NOT scaled!
writeMM(human_named@assays$RNA$data, file = 'A:/cellphone_1_09_19_2025/whole_matrix.mtx')
# save gene and cell names
write(x = rownames(human_named@assays$RNA$data), file = "A:/cellphone_1_09_19_2025/whole_features.tsv")
write(x = colnames(human_named@assays$RNA$data), file = "A:/cellphone_1_09_19_2025/whole_barcodes.tsv")


table(human_named@meta.data$sub.cluster)

human_named@meta.data$Cell = rownames(human_named@meta.data)
df = human_named@meta.data[, c('Cell', 'sub.cluster')]
write.table(df, file ='A:/cellphone_1_09_19_2025/whole_human_named_meta.tsv', sep = '\t', quote = F, row.names = F)

#library('BiocManager')
#BiocManager::install('limma', force =T)
library(limma)

DEGs <- FindAllMarkers(human_named, 
                       group.by = 'sub.cluster',
                       verbose = T, 
                       only.pos = T, 
                       random.seed = 1, 
                       logfc.threshold = 0.2, 
                       min.pct = 0.1, 
                       return.thresh = 0.05)

fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC  > 0.1)

# 1st column = cluster; 2nd column = gene 
fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 

write.table(fDEGs, file ='A:/cellphone_1_09_19_2025/whole_human_named_DEGs.tsv', sep = '\t', quote = F, row.names = F)

for(i in unique(human_named$individual)){
  
  if(i == 'GH'){next}
  print(i)
  dir.create(paste0('A:/cellphone_1_09_19_2025/',i,'/'))
  
  subset_obj = subset(human_named, individual == i)
  
  subset_obj@meta.data$Cell = rownames(subset_obj@meta.data)
  df = subset_obj@meta.data[, c('Cell', 'sub.cluster')]
  meta_path = paste0('A:/cellphone_1_09_19_2025/',i,'/whole_human_named_meta.tsv')
  write.table(df, file =meta_path, sep = '\t', quote = F, row.names = F)
  
  
  DEGs <- FindAllMarkers(subset_obj, 
                         group.by = 'sub.cluster',
                         verbose = T, 
                         only.pos = T, 
                         random.seed = 1, 
                         logfc.threshold = 0.2, 
                         min.pct = 0.1, 
                         return.thresh = 0.05)
  
  fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC  > 0.1)
  
  # 1st column = cluster; 2nd column = gene 
  fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 
  
  file_path = paste0('A:/cellphone_1_09_19_2025/',i,'/whole_DEGs.tsv')
  write.table(fDEGs, file =file_path , sep = '\t', quote = F, row.names = F)
  
}

