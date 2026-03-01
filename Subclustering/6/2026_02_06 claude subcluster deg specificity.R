

# 1. Calculate average expression and percent expressed per subcluster
avg_exp <- AverageExpression(sub_6, 
                              features = degs_6, 
                              group.by = 'sub.cluster',
                              return.seurat = FALSE)$RNA

pct_exp <- lapply(levels(sub_6$sub.cluster), function(cluster) {
  cells <- WhichCells(sub_6, idents = cluster)
  pct <- rowSums(GetAssayData(sub_6, layer = 'data')[degs_6, cells] > 0) / length(cells) * 100
  return(pct)
}) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()
colnames(pct_exp) <- levels(sub_6$sub.cluster)
colnames(avg_exp) <- levels(sub_6$sub.cluster)

# 2. Calculate specificity metrics for each gene
specificity_df <- data.frame(
  gene = degs_6,
  stringsAsFactors = FALSE
)

# For each subcluster, calculate how specific each gene is
for(cluster in levels(sub_6$sub.cluster)) {
  # Average expression in this cluster vs others
  exp_in <- avg_exp[, cluster]
  exp_out <- rowMeans(avg_exp[, colnames(avg_exp) != cluster, drop = FALSE])
  
  # Fold change (log2)
  fc <- log2((exp_in + 0.1) / (exp_out + 0.1))
  
  # Percent difference
  pct_in <- pct_exp[, cluster]
  pct_out <- rowMeans(pct_exp[, colnames(pct_exp) != cluster, drop = FALSE])
  pct_diff <- pct_in - pct_out
  
  # Combined specificity score (you can adjust the weighting)
  specificity <- fc * (pct_diff / 100)
  
  specificity_df[, paste0(cluster, "_fc")] <- fc
  specificity_df[, paste0(cluster, "_pct_diff")] <- pct_diff
  specificity_df[, paste0(cluster, "_specificity")] <- specificity
}

# 3. Assign each gene to its most specific subcluster
specificity_cols <- grep("_specificity$", colnames(specificity_df), value = TRUE)
specificity_df$max_cluster <- apply(specificity_df[, specificity_cols], 1, function(x) {
  cluster_names <- c('6_0', '6_1', '6_2', '6_3')
  cluster_names[which.max(x)]
})
specificity_df$max_specificity <- apply(specificity_df[, specificity_cols], 1, max)

# 4. Sort genes by cluster assignment and specificity
specificity_df <- specificity_df %>%
  arrange(factor(max_cluster, levels = c('6_0', '6_1', '6_2', '6_3')), 
          desc(max_specificity))

# 5. Create ordered gene list
ordered_genes <- specificity_df$gene

# 6. Enhanced DotPlot with ordered genes
degs_6_subclusters_ordered <- DotPlot(sub_6, 
                                       group.by = 'sub.cluster', 
                                       features = ordered_genes,
                                       dot.scale = 6) +
  coord_flip() +
  theme(axis.text.y = element_text(size = 8)) +
  labs(x = "Genes (ordered by subcluster specificity)",
       y = "Subcluster")

# 7. Create a heatmap showing specificity scores
library(pheatmap)
specificity_matrix <- specificity_df %>%
  dplyr::select(gene, ends_with("_specificity")) %>%
  column_to_rownames("gene") %>%
  as.matrix()
colnames(specificity_matrix) <- c('6_0', '6_1', '6_2', '6_3')
specificity_matrix <- specificity_matrix[ordered_genes, ]

pheatmap(specificity_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Gene Specificity Scores by Subcluster")

# 8. Create summary table of top markers per subcluster
top_markers_per_cluster <- specificity_df %>%
  group_by(max_cluster) %>%
  slice_max(order_by = max_specificity, n = 10) %>%
  select(gene, max_cluster, max_specificity, 
         contains("_fc"), contains("_pct_diff"))

print(top_markers_per_cluster)

# 9. Identify genes that are broadly expressed vs subcluster-specific
# Calculate coefficient of variation across clusters
specificity_df$cv_expression <- apply(avg_exp[specificity_df$gene, ], 1, function(x) {
  sd(x) / mean(x)
})

# High CV = more variable across subclusters (subcluster-specific)
# Low CV = more uniform across subclusters (broadly expressed in cluster 6)

specificity_df$expression_pattern <- ifelse(
  specificity_df$cv_expression < median(specificity_df$cv_expression),
  "Broadly expressed in cluster 6",
  "Subcluster-specific"
)

# 10. Visualize broad vs specific
ggplot(specificity_df, aes(x = cv_expression, y = max_specificity, 
                            color = max_cluster, label = gene)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = c('6_0' = 'red', '6_1' = 'blue', 
                                 '6_2' = 'green', '6_3' = 'purple')) +
  geom_vline(xintercept = median(specificity_df$cv_expression), 
             linetype = "dashed", alpha = 0.5) +
  labs(x = "Coefficient of Variation (high = subcluster-specific)",
       y = "Maximum Specificity Score",
       title = "DEG Expression Patterns") +
  theme_minimal()

# 11. Export results
#write.csv(specificity_df, "cluster6_DEG_specificity_analysis.csv", row.names = FALSE)

# 12. Alternative visualization: split DotPlot by expression pattern
degs_specific <- specificity_df %>% 
  filter(expression_pattern == "Subcluster-specific") %>% 
  pull(gene)

degs_broad <- specificity_df %>% 
  filter(expression_pattern == "Broadly expressed in cluster 6") %>% 
  pull(gene)

dotplot_specific <- DotPlot(sub_6, 
                            group.by = 'sub.cluster', 
                            features = degs_specific[1:min(50, length(degs_specific))]) +
  coord_flip() +
  ggtitle("Subcluster-Specific DEGs")

dotplot_broad <- DotPlot(sub_6, 
                         group.by = 'sub.cluster', 
                         features = degs_broad[1:min(50, length(degs_broad))]) +
  coord_flip() +
  ggtitle("Broadly Expressed DEGs in Cluster 6")

# Print summary statistics
cat("\nSummary of DEG patterns:\n")
table(specificity_df$expression_pattern)
cat("\nGenes per subcluster (by max specificity):\n")
table(specificity_df$max_cluster)

