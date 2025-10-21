{
  library(ggsignif)
  library(patchwork)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(CytoTRACE)
  library(BiocManager)
  library(lme4)
  library(enrichR)
  `%notin%` = Negate(`%in%`)
  library(car)
  library(emmeans)
  library(glmGamPoi)
  library(scran)
  library(parallel)
    library(parallel)

  library(Seurat)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(emmeans)
  library(openxlsx)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
  library(biomaRt)
  library(hdWGCNA)
  
  `%notin%` <- Negate(`%in%`)
  
  geneNamer = function(gene){
  names = read.csv('Reference/genes updated.csv')
  
  name = names$NIH_description[names$NIH_accession==gene][1]
  
  if(is.na(name)){name = gene}
  return(name)
}
add_transcript_length = function(gene){
  length = biomart_basic$transcript_length[biomart_basic$entrezgene_accession==gene]%>%min()
  if(length == Inf){length = 0}
  return(length)
}
add_max_expression = function(gene){
  print(gene)
  return(max(counts_matrix[gene,]))
  
}

## whip out biomart
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")

att = listAttributes(ensembl)

biomart_basic <-getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'transcript_length'))
}
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
counts_matrix = obj@assays$RNA$counts%>%as.matrix()

### start -----
#bring in safe genes
preScreenedGenes = read.csv("/Users/ggraham/MERFISH/2025_10_06_sparse_mat_rownames.csv")
preScreenedGenes$transcript_length = sapply(preScreenedGenes$x, add_transcript_length)
preScreenedGenes$max_expression = sapply(preScreenedGenes$x, add_max_expression)

preScreenedGenesGood = subset(preScreenedGenes, max_expression < 138 & transcript_length>= 1440)

Idents(obj) = 'final_clusters'
marks_6_0 = FindMarkers(obj, ident.1=6, ident.2 =c(0,3),features = preScreenedGenesGood$x)

#candidate genes
candidate_gene = read.csv("/Users/ggraham/Desktop/multiome_poa/MERFISH/2025_10_14 NS_MI_genes.csv")
candidate_genes = data.frame(gene =candidate_gene$gene%>%unique())

candidate_genes <- candidate_genes%>%
  distinct(gene)

mer_genes =c(as.list(candidate_genes), as.list(head(rownames(marks_6_0),51)))%>%unlist()
merfish_genes = data.frame(gene = unique(mer_genes))

# test find transfer labels -- thank you claude
library(Seurat)
library(dplyr)
library(ggplot2)

# MERFISH Gene Panel Validation
# Tests if selected genes can accurately map held-out cells back to original clusters

# Parameters
holdout_fraction <- 0.2  # Fraction of cells to hold out for testing
n_neighbors <- 30        # Number of neighbors for label transfer
seed <- 42

set.seed(seed)

# Ensure merfish_genes exist in the object
merfish_genes_filtered <- merfish_genes$gene[merfish_genes$gene %in% rownames(obj)]
cat(sprintf("Using %d/%d MERFISH genes found in object\n", 
            length(merfish_genes_filtered), length(merfish_genes$gene)))

# Split cells into reference and query sets
all_cells <- colnames(obj)
n_test <- round(length(all_cells) * holdout_fraction)
test_cells <- sample(all_cells, n_test)
ref_cells <- setdiff(all_cells, test_cells)

cat(sprintf("\nReference cells: %d\nTest cells: %d\n", 
            length(ref_cells), length(test_cells)))

# Create reference and query objects
ref_obj <- subset(obj, cells = ref_cells)
query_obj <- subset(obj, cells = test_cells)

# Store original cluster labels from query
true_labels <- query_obj$final_clusters
names(true_labels) <- colnames(query_obj)

# Subset to MERFISH genes only
ref_subset <- subset(ref_obj
                     #, features = merfish_genes_filtered
                     )
query_subset <- subset(query_obj
                      # , features = merfish_genes_filtered
                       )

# Normalize and scale data
ref_subset <- NormalizeData(ref_subset, verbose = FALSE)
ref_subset <- ScaleData(ref_subset, verbose = FALSE)
query_subset <- NormalizeData(query_subset, verbose = FALSE)

# Find anchors and transfer labels
cat("\nFinding transfer anchors...\n")
anchors <- FindTransferAnchors(
  reference = ref_subset,
  query = query_subset,
  dims = 1:min(30, length(merfish_genes_filtered) - 1),
  features = merfish_genes$gene,
  verbose = FALSE
)

cat("Transferring labels...\n")
predictions <- TransferData(
  anchorset = anchors,
  refdata = ref_subset$final_clusters,
  dims = 1:min(30, length(merfish_genes_filtered) - 1),
  k.weight = n_neighbors,
  verbose = FALSE
)

# Add predictions to query object
query_obj$predicted_clusters <- predictions$predicted.id
query_obj$prediction_score <- predictions$prediction.score.max

# Calculate accuracy metrics
accuracy <- mean(query_obj$predicted_clusters == true_labels)
cat(sprintf("\n=== RESULTS ===\n"))
cat(sprintf("Overall Accuracy: %.2f%%\n", accuracy * 100))

# Per-cluster accuracy
cluster_accuracy <- query_obj@meta.data %>%
  group_by(final_clusters) %>%
  summarise(
    n_cells = n(),
    correct = sum(predicted_clusters == final_clusters),
    accuracy = mean(predicted_clusters == final_clusters),
    mean_score = mean(prediction_score)
  ) %>%
  arrange(desc(accuracy))

print(cluster_accuracy)

# Confusion matrix
conf_matrix <- table(True = true_labels, Predicted = query_obj$predicted_clusters)
cat("\nConfusion Matrix:\n")
print(conf_matrix)

# Calculate mean prediction score
mean_score <- mean(query_obj$prediction_score)
cat(sprintf("\nMean Prediction Score: %.3f\n", mean_score))

# Visualize results
p1 <- ggplot(query_obj@meta.data, aes(x = final_clusters, fill = predicted_clusters)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "Prediction Distribution by True Cluster",
       x = "True Cluster", y = "Proportion", fill = "Predicted") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(query_obj@meta.data, aes(x = prediction_score)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Prediction Scores",
       x = "Prediction Score", y = "Count") +
  geom_vline(xintercept = mean_score, linetype = "dashed", color = "red")

print(p1)
#print(p2)

# Summary assessment
cat("\n=== ASSESSMENT ===\n")
if (accuracy > 0.9 & mean_score > 0.7) {
  cat("✓ EXCELLENT: Gene panel is highly suitable for MERFISH mapping\n")
} else if (accuracy > 0.8 & mean_score > 0.6) {
  cat("✓ GOOD: Gene panel should work well for MERFISH mapping\n")
} else if (accuracy > 0.7 & mean_score > 0.5) {
  cat("⚠ MODERATE: Gene panel may need refinement\n")
} else {
  cat("✗ POOR: Consider adding more informative genes\n")
}

# Identify problematic clusters
low_acc_clusters <- cluster_accuracy %>%
  filter(accuracy < 0.7) %>%
  pull(final_clusters)

if (length(low_acc_clusters) > 0) {
  cat(sprintf("\nClusters with <70%% accuracy: %s\n", 
              paste(low_acc_clusters, collapse = ", ")))
}

length(merfish_genes$gene)
cluster_accuracy[cluster_accuracy$final_clusters=='6',] 


# which set of genes gives the best accuracy
# ========================
# Iterative gene removal loop
# ========================

results <- data.frame(
  step = integer(),
  removed_gene = character(),
  n_genes = integer(),
  cluster6_accuracy = numeric(),
  stringsAsFactors = FALSE
)

# Copy original gene list so we don't overwrite it
current_genes <- merfish_genes$gene

# Function to run label transfer and get cluster 6 accuracy
get_cluster6_accuracy <- function(gene_list) {
  anchors <- FindTransferAnchors(
    reference = ref_subset,
    query = query_subset,
    dims = 1:min(30, length(gene_list) - 1),
    features = gene_list,
    verbose = FALSE
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = ref_subset$final_clusters,
    dims = 1:min(30, length(gene_list) - 1),
    k.weight = n_neighbors,
    verbose = FALSE
  )
  
  query_obj$predicted_clusters <- predictions$predicted.id
  query_obj$prediction_score <- predictions$prediction.score.max
  
  cluster_acc <- query_obj@meta.data %>%
    group_by(final_clusters) %>%
    summarise(
      accuracy = mean(predicted_clusters == final_clusters)
    )
  
  acc6 <- cluster_acc$accuracy[cluster_acc$final_clusters == 6]
  return(acc6)
}

# Initial accuracy for cluster 6
initial_acc6 <- get_cluster6_accuracy(current_genes)
cat(sprintf("Initial cluster 6 accuracy: %.2f%%\n", initial_acc6*100))

# Loop to remove genes one by one from the end
step <- 1
acc6 <- initial_acc6

while (length(current_genes) > 1 && acc6 >= 0.7) {
  # Remove the last gene
  removed_gene <- tail(current_genes, 1)
  current_genes <- head(current_genes, -1)
  
  # Get new accuracy
  acc6 <- get_cluster6_accuracy(current_genes)
  
  # Save results
  results <- rbind(
    results,
    data.frame(
      step = step,
      removed_gene = removed_gene,
      n_genes = length(current_genes),
      cluster6_accuracy = acc6,
      stringsAsFactors = FALSE
    )
  )
  
  cat(sprintf("Step %d | Removed: %s | Genes left: %d | Acc6: %.2f%%\n",
              step, removed_gene, length(current_genes), acc6*100))
  
  # Break condition
  if (acc6 < 0.7) {
    cat("\n⚠ Stopping: Cluster 6 accuracy dropped below 70%\n")
    break
  }
  
  step <- step + 1
}

# Show results
print(results)
#write.csv(results, "merfish_gene_removal_cluster6_accuracy.csv", row.names = FALSE)
# so the set with scn4ba last which is equivalent to 

optimal_genes = c(as.list(candidate_genes), as.list(head(rownames(marks_6_0),47)))%>%unlist()%>%unique()

# final test ----
merfish_genes = data.frame(gene = unique(c(optimal_genes,
                           'tshz1',
                           'ccdc80',
                           'ntng1a',
                           'LOC111577263'
                           ))
                           )

# test find transfer labels -- thank you claude
library(Seurat)
library(dplyr)
library(ggplot2)

# MERFISH Gene Panel Validation
# Tests if selected genes can accurately map held-out cells back to original clusters

# Parameters
holdout_fraction <- 0.2  # Fraction of cells to hold out for testing
n_neighbors <- 30        # Number of neighbors for label transfer
seed <- 42

set.seed(seed)

# Ensure merfish_genes exist in the object
merfish_genes_filtered <- merfish_genes$gene[merfish_genes$gene %in% rownames(obj)]
cat(sprintf("Using %d/%d MERFISH genes found in object\n", 
            length(merfish_genes_filtered), length(merfish_genes$gene)))

# Split cells into reference and query sets
all_cells <- colnames(obj)
n_test <- round(length(all_cells) * holdout_fraction)
test_cells <- sample(all_cells, n_test)
ref_cells <- setdiff(all_cells, test_cells)

cat(sprintf("\nReference cells: %d\nTest cells: %d\n", 
            length(ref_cells), length(test_cells)))

# Create reference and query objects
ref_obj <- subset(obj, cells = ref_cells)
query_obj <- subset(obj, cells = test_cells)

# Store original cluster labels from query
true_labels <- query_obj$final_clusters
names(true_labels) <- colnames(query_obj)

# Subset to MERFISH genes only
ref_subset <- subset(ref_obj
                     #, features = merfish_genes_filtered
                     )
query_subset <- subset(query_obj
                      # , features = merfish_genes_filtered
                       )

# Normalize and scale data
ref_subset <- NormalizeData(ref_subset, verbose = FALSE)
ref_subset <- ScaleData(ref_subset, verbose = FALSE)
query_subset <- NormalizeData(query_subset, verbose = FALSE)

# Find anchors and transfer labels
cat("\nFinding transfer anchors...\n")
anchors <- FindTransferAnchors(
  reference = ref_subset,
  query = query_subset,
  dims = 1:min(30, length(merfish_genes_filtered) - 1),
  features = merfish_genes$gene,
  verbose = FALSE
)

cat("Transferring labels...\n")
predictions <- TransferData(
  anchorset = anchors,
  refdata = ref_subset$final_clusters,
  dims = 1:min(30, length(merfish_genes_filtered) - 1),
  k.weight = n_neighbors,
  verbose = FALSE
)

# Add predictions to query object
query_obj$predicted_clusters <- predictions$predicted.id
query_obj$prediction_score <- predictions$prediction.score.max

# Calculate accuracy metrics
accuracy <- mean(query_obj$predicted_clusters == true_labels)
cat(sprintf("\n=== RESULTS ===\n"))
cat(sprintf("Overall Accuracy: %.2f%%\n", accuracy * 100))

# Per-cluster accuracy
cluster_accuracy <- query_obj@meta.data %>%
  group_by(final_clusters) %>%
  summarise(
    n_cells = n(),
    correct = sum(predicted_clusters == final_clusters),
    accuracy = mean(predicted_clusters == final_clusters),
    mean_score = mean(prediction_score)
  ) %>%
  arrange(desc(accuracy))

print(cluster_accuracy)

# Confusion matrix
conf_matrix <- table(True = true_labels, Predicted = query_obj$predicted_clusters)
cat("\nConfusion Matrix:\n")
print(conf_matrix)

# Calculate mean prediction score
mean_score <- mean(query_obj$prediction_score)
cat(sprintf("\nMean Prediction Score: %.3f\n", mean_score))

# Visualize results
p1 <- ggplot(query_obj@meta.data, aes(x = final_clusters, fill = predicted_clusters)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "Prediction Distribution by True Cluster",
       x = "True Cluster", y = "Proportion", fill = "Predicted") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(query_obj@meta.data, aes(x = prediction_score)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Prediction Scores",
       x = "Prediction Score", y = "Count") +
  geom_vline(xintercept = mean_score, linetype = "dashed", color = "red")

print(p1)
#print(p2)

# Summary assessment
cat("\n=== ASSESSMENT ===\n")
if (accuracy > 0.9 & mean_score > 0.7) {
  cat("✓ EXCELLENT: Gene panel is highly suitable for MERFISH mapping\n")
} else if (accuracy > 0.8 & mean_score > 0.6) {
  cat("✓ GOOD: Gene panel should work well for MERFISH mapping\n")
} else if (accuracy > 0.7 & mean_score > 0.5) {
  cat("⚠ MODERATE: Gene panel may need refinement\n")
} else {
  cat("✗ POOR: Consider adding more informative genes\n")
}

# Identify problematic clusters
low_acc_clusters <- cluster_accuracy %>%
  filter(accuracy < 0.7) %>%
  pull(final_clusters)

if (length(low_acc_clusters) > 0) {
  cat(sprintf("\nClusters with <70%% accuracy: %s\n", 
              paste(low_acc_clusters, collapse = ", ")))
}

length(merfish_genes$gene)
cluster_accuracy[cluster_accuracy$final_clusters=='6',] 

### output final gene list ----
final_gene_list = data.frame(gene = merfish_genes$gene)
final_gene_list$transcript_length = sapply(final_gene_list$gene, add_transcript_length)
final_gene_list$max_expression = sapply(final_gene_list$gene, add_max_expression)

write.csv(final_gene_list, 'MERFISH/2025_10_14 Final Genes List.csv')
