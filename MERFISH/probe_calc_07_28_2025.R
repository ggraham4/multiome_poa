library(Seurat)
library(tidyverse)
library(dplyr)

obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

obj_subclustered_6 <- FindSubCluster(obj, 6, 'harmony.wsnn')

markers <- FindAllMarkers(obj_subclustered_6, group.by = 'sub.cluster')

#write.csv(markers, "A:/merfish_mutual_information_07_31_2025/all_markers_subclustered_multiome_08_05_2025.csv")

## i am going to make the counts matrix for all marker genes cause the whole
#> matrix is way too big to work with
marker_indicies = which(rownames(obj@assays$RNA$data) %in% unique(markers$gene[markers$p_val_adj<0.0001]))

markers_counts <- obj@assays$RNA$data[marker_indicies,]
markers_counts_matrix <- as.matrix(markers_counts)
#write.csv(markers_counts_matrix, "A:/merfish_mutual_information_07_31_2025/marker_genes.csv")

marker_genes <- as.data.frame(rownames(markers_counts_matrix))
#write.csv(marker_genes, "A:/merfish_mutual_information_07_31_2025/marker_genes_names.csv")
marker_genes <- read_csv("A:/merfish_mutual_information_07_31_2025/marker_genes_names.csv")

### back to regularly scheduled programming

# 50% foreground
markers_50 = subset(markers, pct.1 > 0.5)

#3.3 fold enrichment 
markers_50$fold_enrichment = markers_50$pct.1/markers_50$pct.2

markers_3.3 = subset(markers_50, fold_enrichment>3.3)

#p value less than 0.01
markers_pval = subset(markers_3.3, p_val_adj <0.01)

#clusters
unique(markers_pval$cluster)
length(unique(markers_pval$cluster))

unique(markers$cluster)
length(unique(markers$cluster))

#cluster 9 is missing from stringent markers, I will add those back in manually
setdiff(unique(markers$cluster), unique(markers_pval$cluster))

# im going to try with 2.8 to see if that does better
markers_2.8 = subset(markers_50, fold_enrichment>2.8)
# better for 6 subclusters, still cant find 9

#bring in rna transcript lenght (>500) 
library(biomaRt)

#biomart
clown_mart <- useEnsembl(biomart = 'genes',
                         dataset = 'aocellaris_gene_ensembl')
att = listAttributes(clown_mart)

clown_data <- getBM(mart = clown_mart, #call the human mart
                    attributes = c('entrezgene_accession', #entrezgene_acession is the NIH gene name
                                   'transcript_length'  ))

#select only the shortest transcript
clown_data_shortest <- clown_data %>%
  group_by(entrezgene_accession) %>%
  slice_min(order_by = transcript_length, n = 1, with_ties = FALSE) %>%
  ungroup()

#select genes with transcript length greater than 500
markers_2.8_length = markers_2.8%>%
  right_join(clown_data_shortest, by = join_by('gene' == 'entrezgene_accession'))%>%
  subset(!is.na(fold_enrichment))%>%
  subset(transcript_length>500)

##and highest counts <3000
#only run the counts mat once
counts_mat <- obj_subclustered_6@assays$RNA$counts

#get maximal counts
max_count_func <- function(gene, cluster) {
  # Get the indices of the cells in the specified cluster
  cells_in_cluster <- which(obj_subclustered_6@meta.data$sub.cluster == cluster)
  
  # Get the column index for the gene
  gene_index <- which(rownames(obj_subclustered_6@assays$RNA$counts) == gene)
  
  # Pull the values directly (gene = row, cells = columns)
  gene_counts <- counts_mat[gene_index, cells_in_cluster]
  
  # Return the max
  return(max(gene_counts))
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

max_count_func('esr2b','6_1')

#markers_2.8_length$max_count <- mapply(max_count_func,
#                                       gene = markers_2.8_length$gene,
#                                       cluster = markers_2.8_length$cluster)
markers_2.8_length$mean_norm <- mapply(mean_normalized_expression_clust,
                                       gene = markers_2.8_length$gene,
                                       cluster = markers_2.8_length$cluster)

### find distribution of gene expression, then exclude any degs that are in the top 20%
dist<-rowMeans(obj_subclustered_6@assays$RNA$data)

hist(dist)

bottom_95_percent = quantile(dist, probs = c(.95))
## too low, maybe just markers ig

bottom_90_markers = quantile(markers_2.8_length$mean_norm, probs = c(0.9))

markers_2.8_sorted = markers_2.8_length %>%
  filter(mean_norm <= bottom_90_markers) %>%
  filter(!(gene %in% gene[duplicated(gene) | duplicated(gene, fromLast = TRUE)])) %>%
  group_by(cluster) 

length(unique(markers_2.8_sorted$cluster))
unique(markers_2.8_sorted$cluster)

## have now lost two clusters... I hate this

## anyway, pick top 4 genes per cluster
markers_top4 = markers_2.8_sorted%>%
  filter(!(gene %in% gene[duplicated(gene) | duplicated(gene, fromLast = TRUE)])) %>%
  group_by(cluster) %>%
  arrange(desc(fold_enrichment))%>%
  slice_max(order_by = fold_enrichment, n = 4)

#write.csv(markers_top4,"A:/merfish_mutual_information_07_31_2025/markers_top4_08_04_2025.csv" )
all_possible = markers_2.8_sorted%>%
  filter(!(gene %in% gene[duplicated(gene) | duplicated(gene, fromLast = TRUE)])) %>%
  group_by(cluster) %>%
  arrange(desc(fold_enrichment))

#write.csv(all_possible,"A:/merfish_mutual_information_07_31_2025/markers_multiome_all_possible_08_04_2025.csv" )

#write.csv(gene_df, "A:/merfish_mutual_information_07_31_2025/binarized_expression_subclusters.csv")
#write.csv(as.data.frame(clusters), "A:/merfish_mutual_information_07_31_2025/subclusters.csv")


### write whole expression matrix
#write.csv(as.matrix(obj_subclustered_6@assays$RNA$data),"A:/merfish_mutual_information_07_31_2025/normalized_expression.csv" )

