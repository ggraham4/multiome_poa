#### ATAC HEATMAP for CLUSTER 1 ####

#Read in stuff\
library(Seurat)
library(dplyr)
library(tidyverse)
library(pheatmap)

mecd_ATAC = function(object, gene, cluster, clustering = 'res0.8_50nn_40PC_45LSI'){
  library(stringr)    
  options(dplyr.summarise.inform = FALSE)
  
  counts <- t(object@assays$ATAC$data[,object@meta.data[[clustering]] == cluster]%>%as.matrix())
  Counts_of_interest <- as.data.frame(counts[,gene])
  Counts_of_interest$expression <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
  Counts_of_interest$Status <- object@meta.data$Status[object@meta.data[[clustering]] == cluster]
  
  results <-Counts_of_interest%>%
    group_by(individual, Status)%>%
    summarize(mean = mean(expression),
              se = sd(expression)/sqrt(n()))
  results$Sex <- results$Status
  
  
  results$Sex <- str_sub(results$individual, -1)
  results$Sex[results$individual == 'T17D'] = 'NF'
  results$Sex[results$individual == 'A12D'] = 'E'
  results$Sex[results$individual == 'T11D'] = 'E'
  results$Sex[results$individual == 'GH'] = 'NRM'
  return(results)
}


obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
obj_1_only = subset(obj, res0.8_50nn_40PC_45LSI == 1)
DefaultAssay(obj_1_only) =='ATAC'

dars = read.csv("Collaboration/all_clusters_DARs_peak_level_classified_with_support.csv")

dars_1 =subset(dars, cluster_id ==1)
dars_1

# --- Step 1: Get gene list from cluster 6 DARS ---
dars_to_plot <- dars_1$peak
dars_to_plot <- dars_to_plot[dars_to_plot %in% rownames(obj_1_only@assays$ATAC)]

# --- Step 2: Extract normalized expression and compute z-scores per status ---
expr_matrix <- GetAssayData(obj_1_only, assay = "ATAC", layer = "data")[dars_to_plot, ]
status <- obj_1_only$Status

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
label_map <- dars_1 %>%
  select(peak, classification, nearest_gene) %>%
  distinct() %>%
  filter(peak %in% rownames(zscore_matrix))

row_annotation <- data.frame(
  Label = label_map$classification[match(rownames(zscore_matrix), label_map$peak)],
  row.names = rownames(zscore_matrix)
)

# Map nearest_gene onto the current peak-based row order
gene_labels <- label_map$nearest_gene[match(rownames(zscore_matrix), label_map$peak)]

# --- Step 4: Define label colors using selected palette indices ---
full_palette <- c('#E8ECFB', '#D9CCE3', '#D1BBD7', '#CAACCB', '#BA8DB4',
                  '#AE76A3', '#AA6F9E', '#994F88', '#882E72', '#1965B0',
                  '#437DBF', '#5289C7', '#6195CF', '#7BAFDE', '#4EB265',
                  '#90C987', '#CAE0AB', '#F7F056', '#F7CB45', '#F6C141',
                  '#F4A736', '#F1932D', '#EE8026', '#E8601C', '#E65518',
                  '#DC050C', '#A5170E', '#72190E', '#42150A')

selected_colors <- full_palette[c(10,17, 29)]

unique_labels <- unique(na.omit(row_annotation$Label))
label_color_map <- setNames(selected_colors[seq_along(unique_labels)], unique_labels)

annotation_colors <- list(Label = label_color_map)

# --- Step 5: Plot pheatmap ---
color_palette <- colorRampPalette(c("#1965B0", "white", "#DC050C"))(100)

p = pheatmap(
  (zscore_matrix),
  cluster_rows = F,
  cluster_cols = FALSE,
  treeheight_row = 0,
  treeheight_col = 0,
  color = color_palette,
  breaks = seq(-2, 2, length.out = 101),
  show_rownames = TRUE,
  labels_row = gene_labels,
  show_colnames = TRUE,
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  fontsize_col = 11,
  border_color = NA,
  main = "1_RGC DARs — Mean Z-Score by Phase"
)

ggsave(plot = p,
       file = paste0('dar_matrix_1.svg'),
       device = "svg",
       units = "in",
       height = 2.75,
        width= 7,
       path = paste0("Manuscript/Plots/Manuscript v1.3"))

