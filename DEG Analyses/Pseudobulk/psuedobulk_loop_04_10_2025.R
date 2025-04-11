#edgeR function 

library(edgeR)
library(Seurat)
library(tidyverse)

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

results_list <- list()
for(h in unique(obj$harmony.wnn_res0.4_clusters)){
  message(h)
  obj_subset <- subset(obj, harmony.wnn_res0.4_clusters == h & Status %in% c('M',"D","F"))
  
  counts <- obj_subset@assays$RNA$counts
  metadata <- obj_subset@meta.data
  
  individuals <- unique(metadata$individual)
  individual_counts <- matrix(0, nrow = nrow(counts), ncol = length(individuals))
  rownames(individual_counts) <- rownames(counts)
  colnames(individual_counts) <- individuals

  for(i in seq_along(individuals)) {
  ind <- individuals[i]
  cells_idx <- which(metadata$individual == ind)
  if(length(cells_idx) > 0) {
    individual_counts[, i] <- rowSums(counts[, cells_idx, drop = FALSE])
  }
  }
  
  individual_metadata <- data.frame(
  individual = individuals,
  Status = sapply(individuals, function(ind) {
    # Get status for this individual
    idx <- which(metadata$individual == ind)[1]
    metadata$Status[idx]
  })
)

  dge <- DGEList(counts = individual_counts)
  dge$samples$Status <- factor(individual_metadata$Status)
  
  # Apply TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")
  
  design <- model.matrix(~factor(dge$samples$Status))
  colnames(design) <- c("Intercept", "StatusF", "StatusM")
  
  # Estimate dispersion
  dge <- estimateDisp(dge, design = design)
  
  # Fit GLM
  fit <- glmQLFit(dge, design)

  # Overall test for any differences
  qlf <- glmQLFTest(fit)
  results <- topTags(qlf, n = nrow(dge), sort.by = "PValue")
  
  results_list[[paste0('cluster_',h)]]$results = results
  

}

#ok it works but its ridiculous, it doesnt think brain aromatase is significant and like that's out most robust result???