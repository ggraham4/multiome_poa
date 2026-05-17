library(Seurat)
library(ggplot2)
library(pheatmap)
library(dplyr)

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
obj_6_only = subset(obj, final_clusters == 6)
degs = read.csv('~/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_6 = subset(degs, cluster == 6)



# --- Step 1: Get gene list from cluster 6 DEGs ---
genes_to_plot <- degs_6$gene
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(obj_6_only)]

# --- Step 2: Extract normalized expression and compute z-scores per status ---
expr_matrix <- GetAssayData(obj_6_only, assay = "RNA", layer = "data")[genes_to_plot, ]
status <- obj_6_only$Status

compute_group_means <- function(expr_matrix, status, group) {
  cells <- names(status)[status == group]
  rowMeans(expr_matrix[, cells, drop = FALSE])
}

mean_M <- compute_group_means(expr_matrix, status, "M")
mean_D <- compute_group_means(expr_matrix, status, "D")
mean_F <- compute_group_means(expr_matrix, status, "F")

mean_matrix <- cbind(M = mean_M, D = mean_D, F = mean_F)
zscore_matrix <- t(scale(t(mean_matrix)))
colnames(zscore_matrix) <- c("M", "I", "F")

# --- Step 3: Build annotation for gene labels (short_label) ---
label_map <- degs_6 %>%
  select(gene, short_label) %>%
  distinct() %>%
  filter(gene %in% rownames(zscore_matrix))

row_annotation <- data.frame(
  Label = label_map$short_label[match(rownames(zscore_matrix), label_map$gene)],
  row.names = rownames(zscore_matrix)
)

# --- Step 4: Define label colors using selected palette indices ---
full_palette <- c('#E8ECFB', '#D9CCE3', '#D1BBD7', '#CAACCB', '#BA8DB4',
                  '#AE76A3', '#AA6F9E', '#994F88', '#882E72', '#1965B0',
                  '#437DBF', '#5289C7', '#6195CF', '#7BAFDE', '#4EB265',
                  '#90C987', '#CAE0AB', '#F7F056', '#F7CB45', '#F6C141',
                  '#F4A736', '#F1932D', '#EE8026', '#E8601C', '#E65518',
                  '#DC050C', '#A5170E', '#72190E', '#42150A')

# Selected indices: 10, 14, 15, 17, 18, 26
selected_colors <- full_palette[c(10, 14, 15, 17, 18, 26)]

# Map unique labels to selected colors
unique_labels <- unique(na.omit(row_annotation$Label))
label_color_map <- setNames(selected_colors[seq_along(unique_labels)], unique_labels)

annotation_colors <- list(Label = label_color_map)

# --- Step 5: Plot pheatmap ---
color_palette <- colorRampPalette(c("#1965B0", "white", "#DC050C"))(100)

p = pheatmap(
  zscore_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  treeheight_col = 0,
  color = color_palette,
  breaks = seq(-2, 2, length.out = 101),
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 11,
  border_color = NA,
  main = "6_NPO DEGs — Mean Z-Score by Phase"
)

    ggsave(plot = p,
       file = paste0(i, 'deg_matrix_6.svg'),
       device = "svg",
       units = "in",
       width = 7.5,
       height = 3.8,
       path = paste0("Manuscript/Plots/Manuscript v1.2"))
