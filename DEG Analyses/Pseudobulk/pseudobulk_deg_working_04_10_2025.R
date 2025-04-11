library(edgeR)
library(Seurat)
library(tidyverse)

### Read in object
obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')
obj_19 <- subset(obj, harmony.wnn_res0.4_clusters==19 & Status %in% c('M','D','F'))

# Extract count matrix from Seurat object
counts <- obj_19@assays$RNA$counts
metadata <- obj_19@meta.data

# Aggregate counts by individual
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

# Create individual-level metadata
individual_metadata <- data.frame(
  individual = individuals,
  Status = sapply(individuals, function(ind) {
    # Get status for this individual
    idx <- which(metadata$individual == ind)[1]
    metadata$Status[idx]
  })
)

# Check number of samples per status
table(individual_metadata$Status)

# Create DGEList object
dge <- DGEList(counts = individual_counts)
dge$samples$Status <- factor(individual_metadata$Status)

# Apply TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# For overall test across all groups
design <- model.matrix(~factor(dge$samples$Status))
colnames(design) <- c("Intercept", "StatusF", "StatusM")

# Estimate dispersion
dge <- estimateDisp(dge, design = design)

# Fit GLM
fit <- glmQLFit(dge, design)

# Overall test for any differences
qlf <- glmQLFTest(fit)
results <- topTags(qlf, n = nrow(dge), sort.by = "PValue")
head(results$table)

# Now perform pairwise comparisons
# Create design matrix using treatment contrasts with M as reference
status_levels <- levels(factor(dge$samples$Status))

# For pairwise comparisons
# 1. M vs D
con_M_vs_D <- makeContrasts(StatusD, levels=design)
qlf_M_vs_D <- glmQLFTest(fit, contrast=con_M_vs_D)
results_M_vs_D <- topTags(qlf_M_vs_D, n=nrow(dge), sort.by="PValue")
cat("\nM vs D comparison:\n")
head(results_M_vs_D$table)
cat("Number of significant genes (FDR < 0.05):", sum(results_M_vs_D$table$FDR < 0.05), "\n")

# 2. M vs F
con_M_vs_F <- makeContrasts(StatusF, levels=design)
qlf_M_vs_F <- glmQLFTest(fit, contrast=con_M_vs_F)
results_M_vs_F <- topTags(qlf_M_vs_F, n=nrow(dge), sort.by="PValue")
cat("\nM vs F comparison:\n")
head(results_M_vs_F$table)
cat("Number of significant genes (FDR < 0.05):", sum(results_M_vs_F$table$FDR < 0.05), "\n")

# 3. D vs F
con_D_vs_F <- makeContrasts(StatusF-StatusD, levels=design)
qlf_D_vs_F <- glmQLFTest(fit, contrast=con_D_vs_F)
results_D_vs_F <- topTags(qlf_D_vs_F, n=nrow(dge), sort.by="PValue")
cat("\nD vs F comparison:\n")
head(results_D_vs_F$table)
cat("Number of significant genes (FDR < 0.05):", sum(results_D_vs_F$table$FDR < 0.05), "\n")