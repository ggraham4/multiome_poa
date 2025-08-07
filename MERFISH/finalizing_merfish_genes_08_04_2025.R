## all markers
library(Seurat)
library(tidyverse)
library(biomaRt)

#data
obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

obj_subclustered_6 <- FindSubCluster(obj, 6, 'harmony.wsnn')

#define functions
`%notin%` = Negate(`%in%`)
clown_mart <- useEnsembl(biomart = 'genes',
                         dataset = 'aocellaris_gene_ensembl')
att = listAttributes(clown_mart)

clown_data <- getBM(mart = clown_mart, #call the human mart
                    attributes = c('entrezgene_accession', #entrezgene_acession is the NIH gene name
                                   'transcript_length'  ))

return_transcript_length = function(gene){
  return(min(clown_data$transcript_length[clown_data$entrezgene_accession==gene]))
}

#get mean normalized expression in a cluster
mean_normalized_expression_clust <- function(gene, cluster){
  message(cluster)
  # Get the indices of the cells in the specified cluster
  cells_in_cluster <- which(obj_subclustered_6@meta.data$sub.cluster == cluster)
  
  # Get the column index for the gene
  gene_index <- which(rownames(obj_subclustered_6@assays$RNA$data) == gene)
  
  # Pull the values directly (gene = row, cells = columns)
  gene_counts <- counts_mat[gene_index, cells_in_cluster]
  
  # Return the max
  return(mean(gene_counts))
}

fold_enrichment = function(gene, cluster){
  print(gene)
  cells_in_cluster <- which(obj_subclustered_6@meta.data$sub.cluster == cluster)
  cells_notin_cluster <- which(obj_subclustered_6@meta.data$sub.cluster != cluster)
  
  gene_index <- which(rownames(obj_subclustered_6@assays$RNA$data) == gene)
  
  pct.1 = mean(obj_subclustered_6@assays$RNA$data[gene_index,cells_in_cluster ])
  pct.2 = mean(obj_subclustered_6@assays$RNA$data[gene_index,cells_notin_cluster ])
  
  fold_enrichment = pct.1/pct.2
  return(fold_enrichment)
}

mean_expression = function(gene, cluster){
  print(gene)
  cells_in_cluster <- which(obj_subclustered_6@meta.data$sub.cluster == cluster)

  gene_index <- which(rownames(obj_subclustered_6@assays$RNA$data) == gene)
  
  pct.1 = mean(obj_subclustered_6@assays$RNA$data[gene_index,cells_in_cluster ])

  return(pct.1)
}

identify_strongest_cluster = function(gene){
  print(gene)
  gene_index = which(rownames(obj_subclustered_6@assays$RNA$data)==gene)
  
  newd = data.frame(clusters = obj_subclustered_6$sub.cluster,
                    expression = obj_subclustered_6@assays$RNA$data[gene,])
  
  newd_grouped = newd%>%
    group_by(clusters)%>%
    summarize(mean_exp = mean(expression))
  
  return(newd_grouped$clusters[newd_grouped$mean_exp == max(newd_grouped$mean_exp)])
  
}

cells_in_cluster_expressing = function(gene, cluster){
  gene_index = which(rownames(obj_subclustered_6@assays$RNA$data)==gene)
  cells_in_cluster <- which(obj_subclustered_6@meta.data$sub.cluster == cluster)
  
  expression_in_cluster=sum(obj_subclustered_6@assays$RNA$data[gene_index, cells_in_cluster]>0)
  pct.1 = expression_in_cluster/length(cells_in_cluster)
  return(pct.1)
  
}

cells_notin_cluster_expressing = function(gene, cluster){
  print(gene)
  gene_index = which(rownames(obj_subclustered_6@assays$RNA$data)==gene)
  cells_notin_cluster <- which(obj_subclustered_6@meta.data$sub.cluster != cluster)
  
  expression_notin_cluster=sum(obj_subclustered_6@assays$RNA$data[gene_index, cells_notin_cluster]>0)
  pct.1 = expression_notin_cluster/length(cells_notin_cluster)
  return(pct.1)
  
}
dm_markers = read_csv("MERFISH/markers_dm.csv")
dm_markers$cluster ='Dm'

multiome_markers = read_csv("MERFISH/markers_multiome_obj_only.csv")

other_multiome_markers <- read_csv("A:/merfish_mutual_information_07_31_2025/markers_multiome_all_possible_08_04_2025.csv")

all_markers <- append(c(dm_markers$gene), multiome_markers$`0`)


good_markers_6_0 = FindMarkers(obj_subclustered_6, '6_0', group.by = 'sub.cluster')
markers_6_0vs_6_1 = FindMarkers(obj_subclustered_6, ident.1='6_0', ident.2 = '6_1', group.by = 'sub.cluster')
markers_6_0vs_6_1_clean = subset(markers_6_0vs_6_1, pct.1 > pct.2 & p_val_adj <0.001)

markers_6_0vs_6_2 = FindMarkers(obj_subclustered_6, ident.1='6_0', ident.2 = '6_2', group.by = 'sub.cluster')
markers_6_0vs_6_2_clean = subset(markers_6_0vs_6_2, pct.1 > pct.2 & p_val_adj <0.001)

markers_6_0vs_6_3 = FindMarkers(obj_subclustered_6, ident.1='6_0', ident.2 = '6_3', group.by = 'sub.cluster')
markers_6_0vs_6_3_clean = subset(markers_6_0vs_6_3, pct.1 > pct.2 & p_val_adj <0.001)

unique_6_0_markers <- intersect(rownames(markers_6_0vs_6_1_clean), rownames(markers_6_0vs_6_2_clean))
unique_6_0_markers_2 <- intersect(rownames(markers_6_0vs_6_3_clean), rownames(markers_6_0vs_6_2_clean))

unique_6_0_markers_full <-  intersect(unique_6_0_markers_2, rownames(markers_6_0vs_6_3_clean))

markers_6_0vs_6_3_clean%>%head(10)
markers_6_0vs_6_2_clean%>%head(10)
markers_6_0vs_6_1_clean%>%head(10)
## ok so I will add zfhx4, zfhx3b, LOC11157418, ar, esr2b

## I alsp want to add some from 6 at large
good_markers_6 = FindMarkers(obj, 6)

#hmx3a, hmx3b, hmx2, ar is already listed, pgr, gal, esr2a, eomesb, tacr3a

good_markers_6_0$gene = rownames(good_markers_6_0)
good_markers_6_0$transcript_length = lapply(X=good_markers_6_0$gene, FUN = return_transcript_length)
good_markers_6_0_filtered <- head(good_markers_6_0, 100)
good_markers_6_0_filtered$fold_enrichment = mapply(FUN = fold_enrichment,
                                          gene = good_markers_6_0_filtered$gene,
                                          cluster = '6_0')
good_markers_6_0_filtered$mean_expression = mapply(FUN = mean_expression,
                                                   gene = good_markers_6_0_filtered$gene,
                                                   cluster = '6_0')

good_markers_6_1 = FindMarkers(obj_subclustered_6,'6_1', group.by ='sub.cluster')
good_markers_7 = FindMarkers(obj_subclustered_6,'7', group.by ='sub.cluster')
good_markers_7[rownames(good_markers_7)%notin% all_markers,]

good_markers_3 = FindMarkers(obj_subclustered_6,'3', group.by ='sub.cluster')
good_markers_9 = FindMarkers(obj_subclustered_6,'9', group.by ='sub.cluster')

final_genes_list <- unique(append(all_markers, c('zfhx4',
                                          'zfhx3b',
                                          'hmx3a',
                                          'hmx3b',
                                          'ar',
                                          'esr2a',
                                          'esr2b',
                                          'cckbrb',
                                          'tacr3a',
                                          'galr1a',
                                          'pgr',
                                          'nfixa',
                                          'znf536',
                                          'slc32a1',
                                          'nrxn3a',
                                          'chrdl2',
                                          'wnk1b')))

final_df = data.frame(gene = final_genes_list)
final_df$transcript_length = lapply(X=final_df$gene, FUN = return_transcript_length)
final_df$strongest_cluster = lapply(X=final_df$gene, FUN = identify_strongest_cluster)
final_df$mean_exp_strongest_clust =  mapply( FUN = mean_expression, gene = final_df$gene,cluster = final_df$strongest_cluster)
final_df$expression_enrichment =  mapply( FUN = fold_enrichment, gene = final_df$gene,cluster = final_df$strongest_cluster)
final_df$pct_in_cluster_expressing =  mapply( FUN = cells_in_cluster_expressing, gene = final_df$gene,cluster = final_df$strongest_cluster)
final_df$pct_notin_cluster_expressing =  mapply( FUN = cells_notin_cluster_expressing, gene = final_df$gene,cluster = final_df$strongest_cluster)
final_df$pct_fold_enrichment = final_df$pct_in_cluster_expressing/final_df$pct_notin_cluster_expressing
final_df$transcript_length <- unlist(final_df$transcript_length)

for(i in colnames(final_df)){
  final_df[[i]] = unlist(final_df[[i]])
}

#write.csv(final_df, 'C:/Users/Gabe/Desktop/multiome_poa/MERFISH/final_gene_selection_08_04_2025.csv')

table(final_df$strongest_cluster)
