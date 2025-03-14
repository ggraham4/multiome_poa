library(Seurat)
library(SeuratObject)
library(Matrix)

options(Seurat.object.assay.version = "v3")

human_named = readRDS("C:/Users/Gabe/Desktop/RNA_object_human_names.rds")

human_named$harmony.wnn_res0.4_clusters <- as.character(human_named$harmony.wnn_res0.4_clusters)
Idents(human_named) = human_named$harmony.wnn_res0.4_clusters

# Save normalised counts - NOT scaled!
writeMM(human_named@assays$RNA$data, file = 'A:/CellPhoneDB 030225/matrix.mtx')
# save gene and cell names
write(x = rownames(human_named@assays$RNA$data), file = "A:/CellPhoneDB 030225/features.tsv")
write(x = colnames(human_named@assays$RNA$data), file = "A:/CellPhoneDB 030225/barcodes.tsv")


table(human_named@meta.data$harmony.wnn_res0.4_clusters)

human_named@meta.data$Cell = rownames(human_named@meta.data)
df = human_named@meta.data[, c('Cell', 'harmony.wnn_res0.4_clusters')]
write.table(df, file ='A:/CellPhoneDB 030225/human_named_meta.tsv', sep = '\t', quote = F, row.names = F)

library('BiocManager')
BiocManager::install('limma', force =T)
library(limma)

DEGs <- FindAllMarkers(human_named, 
                                               verbose = T, 
                                               only.pos = T, 
                                               random.seed = 1, 
                                               logfc.threshold = 0.2, 
                                               min.pct = 0.1, 
                                               return.thresh = 0.05)

fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC  > 0.1)

# 1st column = cluster; 2nd column = gene 
fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 

write.table(fDEGs, file ='A:/CellPhoneDB 030225/human_named_DEGs.tsv', sep = '\t', quote = F, row.names = F)


                       