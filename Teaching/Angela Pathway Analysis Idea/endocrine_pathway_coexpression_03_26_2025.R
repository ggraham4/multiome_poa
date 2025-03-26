### steroid pathway via coexpression
library(biomaRt)
library(Seurat)
library(tidyverse)

obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

obj$endocrine <- ifelse(obj@assays$RNA$data['hsd17b14',] >0, 'endocrine', NA)
obj$endocrine <- ifelse(obj@assays$RNA$data['hsd3b1',] >0, 'endocrine', obj$endocrine)
obj$endocrine <- ifelse(obj@assays$RNA$data['hsd17b1',] >0, 'endocrine', obj$endocrine)
obj$endocrine <- ifelse(obj@assays$RNA$data['LOC111577263',] >0, 'endocrine', obj$endocrine)
obj$endocrine <- ifelse(is.na(obj$endocrine), 'no_endocrine', obj$endocrine)

DimPlot(obj, split.by = 'endocrine')

endocrine_data_cluster_level <- table(obj@meta.data$endocrine, obj@meta.data$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  pivot_wider(names_from = 'Var1', values_from = 'Freq')%>%
  mutate(prop_pos = endocrine/(endocrine+no_endocrine))
#overwhelmingly glial

### FindMarkers
#Find genes differentially expressed between IEG and IEG- clusters
endocrine_marker_all <-data.frame()
for(i in unique(obj$harmony.wnn_res0.4_clusters)){
  print(i)
  temp_obj <- subset(obj, harmony.wnn_res0.4_clusters==i)
  Idents(temp_obj) = 'endocrine'
  endocrine_markers <- FindAllMarkers(temp_obj, 
                                assay = 'RNA', 
                                group.by = 'endocrine',
                                logfc.threshold = 0,
                                min.pct = 1/204
  ) #204 cells in smallest cluster
  endocrine_markers$harmony.wnn_res0.4_clusters =i
  endocrine_marker_all <- rbind(endocrine_marker_all, endocrine_markers)
}
#extract significant endocrine like genes
endocrine_marker_all_signif <- endocrine_marker_all[endocrine_marker_all$p_val_adj<0.05 &endocrine_marker_all$cluster=='endocrine' ,]

#find genes significant in over half of clusters
endocrine_like_genes_data <- table(endocrine_marker_all_signif$gene,
                             endocrine_marker_all_signif$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  group_by(Var1)%>% #cluster
  summarize(n_clusters = sum(Freq))%>%
  subset(n_clusters > length(unique(obj$harmony.wnn_res0.4_clusters))/4 #reducing stringency to 1/4
  ) 

#endocrine like genes
endocrine_like_genes <- endocrine_like_genes_data$Var1%>%droplevels()
# ok well this includes EAAT1 so that might not be good
#also laminin
