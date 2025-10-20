library(Seurat)
library(dplyr)
library(ggplot2)

# ============ PARAMETERS ============
holdout_fraction <- 0.2
n_neighbors <- 30
target_accuracy <- 0.75
max_expression_threshold <- 137
max_iterations <- 50   # safety stop

set.seed(42)

# ============ INITIAL SETUP ============
all_cells <- colnames(obj)
n_test <- round(length(all_cells) * holdout_fraction)
test_cells <- sample(all_cells, n_test)
ref_cells <- setdiff(all_cells, test_cells)

ref_obj <- subset(obj, cells = ref_cells)
query_obj <- subset(obj, cells = test_cells)
true_labels <- query_obj$final_clusters
names(true_labels) <- colnames(query_obj)

# Ensure genes exist in object
merfish_genes_filtered <- merfish_genes$gene[merfish_genes$gene %in% rownames(obj@assays$RNA$data)]
cat(sprintf("Using %d/%d MERFISH genes found in object\n", 
            length(merfish_genes_filtered), length(merfish_genes$gene)))

# Pre-normalize reference (static)
ref_subset <- subset(ref_obj)
ref_subset <- NormalizeData(ref_subset, verbose = FALSE)
ref_subset <- ScaleData(ref_subset, verbose = FALSE)

# ============ HELPER: accuracy calculation ============
cluster_accuracy_function <- function(query_obj) {
  query_obj@meta.data %>%
    group_by(final_clusters) %>%
    summarise(
      n_cells = n(),
      correct = sum(predicted_clusters == final_clusters),
      accuracy = mean(predicted_clusters == final_clusters),
      mean_score = mean(prediction_score)
    )
}

# ============ HELPER: classification run ============
run_label_transfer <- function(ref_subset, query_obj, gene_panel) {
  query_subset <- subset(query_obj)
  query_subset <- NormalizeData(query_subset, verbose = FALSE)
  
  dims_to_use <- min(30, length(gene_panel) - 1)
  
  anchors <- FindTransferAnchors(
    reference = ref_subset,
    query = query_subset,
    dims = 1:dims_to_use,
    features = gene_panel,
    verbose = FALSE
  )
  
  predictions <- TransferData(
    anchorset = anchors,
    refdata = ref_subset$final_clusters,
    dims = 1:dims_to_use,
    k.weight = n_neighbors,
    verbose = FALSE
  )
  
  query_subset$predicted_clusters <- predictions$predicted.id
  query_subset$prediction_score <- predictions$prediction.score.max
  return(query_subset)
}

# ============ MAIN LOOP ============
iteration <- 0
replacement_genes <- c()
used_genes <- c()

current_genes <- merfish_genes_filtered
gene_panel_length <- length(current_genes)

# Number of top markers to test each iteration
top_n_candidates <- 20

repeat {
  iteration <- iteration + 1
  cat("\n--- Iteration", iteration, "---\n")
  
  # Run current panel
  query_result <- run_label_transfer(ref_subset, query_obj, current_genes)
  
  # Assign tmp_ident for FindMarkers
  misclassified_6 <- rownames(query_result@meta.data[
    query_result$final_clusters == 6 & query_result$predicted_clusters != 6
  ,])
  cells_0 <- rownames(query_result@meta.data[
    query_result$final_clusters == 0
    , ])
  
  query_result$tmp_ident <- "other"
  query_result$tmp_ident[misclassified_6] <- "misclassified6"
  query_result$tmp_ident[cells_0] <- "cells_0"
  Idents(query_result) <- "tmp_ident"
  
  # Compute accuracy
  acc_table <- cluster_accuracy_function(query_result)
  cluster6_acc <- acc_table$accuracy[acc_table$final_clusters == 6]
  cluster6_acc <- ifelse(length(cluster6_acc) == 0, 0, cluster6_acc)
  cat(sprintf("Cluster 6 accuracy: %.2f%%\n", cluster6_acc * 100))
  
  if (cluster6_acc >= target_accuracy | iteration >= max_iterations) break
  
  # Find top N candidate genes for misclassified cells
  candidate_markers <- FindMarkers(
    query_result,
    ident.1 = "misclassified6",
    ident.2 = "cells_0",
    group.by = "tmp_ident"
  )
  candidate_markers <- candidate_markers[order(-candidate_markers$avg_log2FC), , drop=FALSE]
  candidate_genes <- setdiff(rownames(candidate_markers), c(current_genes, used_genes))
  candidate_genes <- candidate_genes[1:min(top_n_candidates, length(candidate_genes))]
  
  if (length(candidate_genes) == 0) {
    cat("⚠ No new candidate genes available.\n")
    break
  }
  
  # Test each candidate by simulating replacement
  best_gene <- NA
  best_acc <- cluster6_acc
  row_to_replace <- gene_panel_length - iteration + 1
  
  for (gene in candidate_genes) {
    test_panel <- current_genes
    test_panel[row_to_replace] <- gene
    test_result <- run_label_transfer(ref_subset, query_obj, test_panel)
    test_acc <- cluster_accuracy_function(test_result)$accuracy[
      cluster_accuracy_function(test_result)$final_clusters == 6
    ]
    if (test_acc > best_acc) {
      best_acc <- test_acc
      best_gene <- gene
    }
  }
  
  if (!is.na(best_gene)) {
    cat("Replacing gene at row", row_to_replace,
        current_genes[row_to_replace], "with", best_gene,
        "→ new cluster 6 accuracy:", best_acc*100, "%\n")
    current_genes[row_to_replace] <- best_gene
    replacement_genes <- c(replacement_genes, best_gene)
    used_genes <- c(used_genes, best_gene)
  } else {
    cat("⚠ No candidate improved accuracy this iteration. Skipping.\n")
  }
}

# ============ REPORT ============
final_acc_table <- cluster_accuracy_function(query_result)

cat("\n===== FINAL REPORT =====\n")
cat(sprintf("Final Cluster 6 Accuracy: %.2f%%\n", 
            final_acc_table$accuracy[final_acc_table$final_clusters == 6] * 100))
cat(sprintf("Total Replacements: %d\n", length(replacement_genes)))
cat("Genes added/replaced:\n")
print(replacement_genes)

replacement_df <- data.frame(
  iteration = seq_along(replacement_genes),
  replaced_gene = replacement_genes
)
print(replacement_df)
