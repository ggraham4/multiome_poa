## find good markers for dm

library(Seurat)
library(tidyverse)
library(dplyr)

old_obj = readRDS("A:/Anemonefish POA Legacy R Objects/Old object.rds")

DimPlot(old_obj, label = T, group.by = 'clusters49')

markers <- FindMarkers(old_obj, 2, group.by = 'clusters49')

# 50% foreground
markers_50 = subset(markers, pct.1 > 0.5)

#3.3 fold enrichment 
markers_50$fold_enrichment = markers_50$pct.1/markers_50$pct.2

markers_3.3 = subset(markers_50, fold_enrichment>3.3)

#p value less than 0.01
markers_pval = subset(markers_3.3, p_val_adj <0.01)

library(biomaRt)

#biomart
clown_mart <- useEnsembl(biomart = 'genes',
                         dataset = 'apercula_gene_ensembl')
att = listAttributes(clown_mart)

clown_data <- getBM(mart = clown_mart, #call the human mart
                    attributes = c('external_gene_name', #entrezgene_acession is the NIH gene name
                                   'transcript_length'  ))

#select only the shortest transcript
clown_data_shortest <- clown_data %>%
  group_by(external_gene_name) %>%
  slice_min(order_by = transcript_length, n = 1, with_ties = FALSE) %>%
  ungroup()

#select genes with transcript length greater than 500
markers_3.3$gene = rownames(markers_3.3)
markers_3.3_length = markers_3.3%>%
  right_join(clown_data_shortest, by = join_by('gene' == 'external_gene_name'))%>%
  subset(!is.na(fold_enrichment))%>%
  subset(transcript_length>500)

mean_normalized_expression_clust <- function(gene){
  # Get the indices of the cells in the specified cluster
    cells_in_cluster <- which(old_obj@meta.data$clusters49 == 2)
  
  # Get the column index for the gene
  gene_index <- which(rownames(old_obj@assays$RNA$data) == gene)
  
  # Pull the values directly (gene = row, cells = columns)
  gene_counts <- counts_mat[gene_index, cells_in_cluster]
  
  # Return the max
  return(mean(gene_counts))
  
  
}

counts_mat = as.matrix(old_obj@assays$RNA$data)

markers_3.3_length$mean_norm <- mapply(mean_normalized_expression_clust,
                                       gene = markers_3.3_length$gene)
# expression is fine

markers_3.3_trimmed = markers_3.3_length[1:4,]

#rite.csv(markers_3.3_trimmed, "MERFISH/markers_dm.csv")
