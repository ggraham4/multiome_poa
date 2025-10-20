library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)

# Load data
obj <- readRDS("~/Desktop/optimal_clustering_rna_only.rds")

# Calculate max expression for all genes
counts_matrix <- obj@assays$RNA$counts %>% as.matrix()
max_expression <- apply(counts_matrix, 1, max)

# Filter genes by max expression threshold
genes_passing_threshold <- names(max_expression[max_expression <= 137])

cat(sprintf("Total genes: %d\n", length(max_expression)))
cat(sprintf("Genes with max expression <= 137: %d\n", length(genes_passing_threshold)))

# Find markers for cluster 6 vs 0 and 3
all_marks_6_vs_0_3 <- FindMarkers(obj, 
                                  ident.1 = '6', 
                                  ident.2 = c('0', '3'),
                                  group.by = 'final_clusters',
                                  logfc.threshold = 0.25)

# Filter markers to only include genes passing expression threshold
all_marks_6_vs_0_3 <- all_marks_6_vs_0_3[rownames(all_marks_6_vs_0_3) %in% genes_passing_threshold, ]

cat(sprintf("Markers for 6 vs 0+3 (after expression filter): %d\n", nrow(all_marks_6_vs_0_3)))

# Load base candidate genes
candidate_gene <- read.xlsx("MERFISH/candidate_genes.xlsx")
base_candidate_genes <- data.frame(gene = candidate_gene$gene %>% unique())
base_candidate_genes$gene[base_candidate_genes$gene == 'fgf12a'] <- 'hmx2'
base_candidate_genes$gene[base_candidate_genes$gene == 'grb10b'] <- 'hmx3b'
base_candidate_genes <- base_candidate_genes %>% distinct(gene)

# Filter base genes to only include those passing expression threshold
base_candidate_genes <- base_candidate_genes[base_candidate_genes$gene %in% genes_passing_threshold, ]

# Keep only first 97 genes (leaving room for 23 test genes)
base_genes <- base_candidate_genes[1:min(97, nrow(base_candidate_genes))]

cat(sprintf("Base genes (after expression filter): %d\n", length(base_genes)))

# Get available markers for testing (filter out genes already in base)
available_markers <- all_marks_6_vs_0_3[!(rownames(all_marks_6_vs_0_3) %in% base_genes), ]
available_markers <- available_markers[available_markers$p_val_adj < 0.05, ]

cat(sprintf("Base genes: %d\n", length(base_genes)))
cat(sprintf("Available markers for testing: %d\n", nrow(available_markers)))

# Function to test gene panel
test_gene_panel <- function(test_genes, base_genes, obj, seed = 42) {
  
  # Combine base + test genes
  full_gene_list <- c(base_genes, test_genes)
  
  # Filter to genes in object
  genes_filtered <- full_gene_list[full_gene_list %in% rownames(obj)]
  
  if(length(genes_filtered) < 50) {
    return(list(accuracy_cluster6 = 0, 
                overall_accuracy = 0,
                n_genes = length(genes_filtered)))
  }
  
  set.seed(seed)
  
  # Split data
  holdout_fraction <- 0.2
  all_cells <- colnames(obj)
  n_test <- round(length(all_cells) * holdout_fraction)
  test_cells <- sample(all_cells, n_test)
  ref_cells <- setdiff(all_cells, test_cells)
  
  # Create reference and query
  ref_obj <- subset(obj, cells = ref_cells)
  query_obj <- subset(obj, cells = test_cells)
  
  true_labels <- query_obj$final_clusters
  
  # Subset to test genes
  ref_subset <- subset(ref_obj, features = genes_filtered)
  query_subset <- subset(query_obj, features = genes_filtered)
  
  # Normalize
  ref_subset <- NormalizeData(ref_subset, verbose = FALSE)
  ref_subset <- ScaleData(ref_subset, verbose = FALSE)
  query_subset <- NormalizeData(query_subset, verbose = FALSE)
  
  # Find anchors
  tryCatch({
    anchors <- FindTransferAnchors(
      reference = ref_subset,
      query = query_subset,
      dims = 1:min(30, length(genes_filtered) - 1),
      verbose = FALSE
    )
    
    predictions <- TransferData(
      anchorset = anchors,
      refdata = ref_subset$final_clusters,
      dims = 1:min(30, length(genes_filtered) - 1),
      k.weight = 30,
      verbose = FALSE
    )
    
    query_obj$predicted_clusters <- predictions$predicted.id
    
    # Calculate cluster 6 specific accuracy
    cluster_6_cells <- true_labels == '6'
    cluster_6_accuracy <- mean(query_obj$predicted_clusters[cluster_6_cells] == '6')
    
    # Overall accuracy
    overall_accuracy <- mean(query_obj$predicted_clusters == true_labels)
    
    return(list(
      accuracy_cluster6 = cluster_6_accuracy,
      overall_accuracy = overall_accuracy,
      n_genes = length(genes_filtered),
      n_cluster6_cells = sum(cluster_6_cells)
    ))
    
  }, error = function(e) {
    return(list(accuracy_cluster6 = 0, 
                overall_accuracy = 0,
                n_genes = length(genes_filtered),
                error = as.character(e)))
  })
}

# Optimization loop
n_iterations <- 100  # Number of random gene sets to test
n_genes_to_test <- 23

results <- data.frame(
  iteration = integer(),
  accuracy_cluster6 = numeric(),
  overall_accuracy = numeric(),
  n_genes = integer(),
  genes = character(),
  stringsAsFactors = FALSE
)

cat("\n=== Starting Optimization ===\n")
cat(sprintf("Testing %d random combinations of %d genes\n\n", 
            n_iterations, n_genes_to_test))

set.seed(123)  # For reproducibility of gene selection

for(i in 1:n_iterations) {
  
  # Randomly sample 23 genes from available markers
  if(nrow(available_markers) >= n_genes_to_test) {
    sampled_genes <- sample(rownames(available_markers), n_genes_to_test)
  } else {
    sampled_genes <- rownames(available_markers)
  }
  
  # Test this combination
  result <- test_gene_panel(sampled_genes, base_genes, obj)
  
  # Store results
  results <- rbind(results, data.frame(
    iteration = i,
    accuracy_cluster6 = result$accuracy_cluster6,
    overall_accuracy = result$overall_accuracy,
    n_genes = result$n_genes,
    genes = paste(sampled_genes, collapse = ","),
    stringsAsFactors = FALSE
  ))
  
  # Progress update
  if(i %% 10 == 0) {
    cat(sprintf("Completed %d/%d iterations. Best cluster 6 accuracy so far: %.2f%%\n", 
                i, n_iterations, max(results$accuracy_cluster6) * 100))
  }
}

# Find best result
best_result <- results[which.max(results$accuracy_cluster6), ]
best_genes <- strsplit(best_result$genes, ",")[[1]]

cat("\n=== OPTIMIZATION RESULTS ===\n")
cat(sprintf("Best Cluster 6 Accuracy: %.2f%%\n", best_result$accuracy_cluster6 * 100))
cat(sprintf("Overall Accuracy: %.2f%%\n", best_result$overall_accuracy * 100))
cat(sprintf("Iteration: %d\n\n", best_result$iteration))

cat("Best 23 genes:\n")
print(best_genes)

# Show top 10 results
cat("\n=== Top 10 Gene Combinations ===\n")
top_results <- results %>% 
  arrange(desc(accuracy_cluster6)) %>% 
  head(10) %>%
  select(iteration, accuracy_cluster6, overall_accuracy)
print(top_results)

# Plot results
p1 <- ggplot(results, aes(x = accuracy_cluster6 * 100)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = best_result$accuracy_cluster6 * 100, 
             color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Distribution of Cluster 6 Accuracy Across Gene Combinations",
       x = "Cluster 6 Accuracy (%)",
       y = "Count") +
  annotate("text", 
           x = best_result$accuracy_cluster6 * 100, 
           y = Inf, 
           label = sprintf("Best: %.1f%%", best_result$accuracy_cluster6 * 100),
           vjust = 2, color = "red")

p2 <- ggplot(results, aes(x = accuracy_cluster6, y = overall_accuracy)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_point(data = best_result, aes(x = accuracy_cluster6, y = overall_accuracy),
             color = "red", size = 3) +
  theme_minimal() +
  labs(title = "Cluster 6 Accuracy vs Overall Accuracy",
       x = "Cluster 6 Accuracy",
       y = "Overall Accuracy")

print(p1)
print(p2)

# Create final optimized gene list
optimized_genes <- data.frame(gene = c(base_genes, best_genes))

# Save results
write.csv(results, "MERFISH/gene_optimization_results.csv", row.names = FALSE)
write.csv(optimized_genes, "MERFISH/optimized_gene_panel.csv", row.names = FALSE)

cat("\n=== FILES SAVED ===\n")
cat("All results: MERFISH/gene_optimization_results.csv\n")
cat("Optimized gene panel: MERFISH/optimized_gene_panel.csv\n")

# Final validation with best genes
cat("\n=== FINAL VALIDATION ===\n")
final_result <- test_gene_panel(best_genes, base_genes, obj, seed = 999)
cat(sprintf("Cluster 6 Accuracy (validation): %.2f%%\n", 
            final_result$accuracy_cluster6 * 100))
cat(sprintf("Overall Accuracy (validation): %.2f%%\n", 
            final_result$overall_accuracy * 100))