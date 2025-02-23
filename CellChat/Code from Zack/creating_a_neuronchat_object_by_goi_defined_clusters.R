library(Seurat)

data <- readRDS("C:/Users/zvjohns/Desktop/projects/clownfish/rds/clownfish_OB_signal_qc_new_cluster_names.rds")


# Function to subset Seurat object based on gene expression
subset_seurat_by_gene <- function(seurat_obj, gene_name) {
  # Get expression data for the specific gene
  gene_expression <- seurat_obj@assays$RNA@data[gene_name, ]
  
  # Cells where expression of the gene is > 0
  cells_with_expression <- names(gene_expression[gene_expression > 0])
  
  # Subset the Seurat object based on these cells
  subset_obj <- subset(seurat_obj, cells = cells_with_expression)
  
  # Add a new column with the gene name
  subset_obj$gene_of_interest <- gene_name
  
  return(subset_obj)
}

# Given Seurat object and vector of genes
seurat_obj <- data # Replace this with your actual Seurat object
gene_csv <- read.csv("C:/Users/zvjohns/Desktop/projects/clownfish/csv/input/genes_interest.csv")
genes_interest <- subset(gene_csv, gene_csv$expressed==TRUE)
gene_vector <- as.vector(genes_interest$geneID) # Replace this with your actual list of genes

# Loop through the gene vector and create subsetted Seurat objects
list_of_seurat_objects <- lapply(gene_vector, function(gene) subset_seurat_by_gene(seurat_obj, gene))

# Merge all the subsetted Seurat objects into a single object
merged_seurat_obj <- list_of_seurat_objects[[1]]
for (i in 2:length(list_of_seurat_objects)) {
  merged_seurat_obj <- merge(merged_seurat_obj, y = list_of_seurat_objects[[i]])
}

# Now, 'merged_seurat_obj' contains the merged Seurat object

min(table(merged_seurat_obj$subsample_z, merged_seurat_obj$gene_of_interest))
min(table(merged_seurat_obj_filtered$subsample_z, merged_seurat_obj_filtered$gene_of_interest))

#saveRDS(merged_seurat_obj, "C:/Users/zvjohns/Desktop/projects/clownfish/rds/clownfish_OB_signal_qc_new_cluster_names_goi.rds")

## now subset only those clusters in which all individuals are represented

# Assuming you already have a Seurat object called 'seurat_obj'
# Replace 'seurat_obj' with the name of your actual Seurat object

# Check if the 'subsample_z' and 'gene_of_interest' columns exist in metadata

if ("subsample_z" %in% colnames(merged_seurat_obj@meta.data) && "gene_of_interest" %in% colnames(merged_seurat_obj@meta.data)) {
  
  # Get unique values in 'subsample_z' and 'gene_of_interest' columns
  unique_subsample_z <- unique(merged_seurat_obj$subsample_z)
  unique_gene_of_interest <- unique(merged_seurat_obj$gene_of_interest)
  
  # Initialize a vector to store gene names that meet the criteria
  selected_genes <- character(0)
  
  # Loop through unique values in 'gene_of_interest'
  for (gene_value in unique_gene_of_interest) {
    # Check if this gene value is detected in all unique 'subsample_z' values
    if (all(unique_subsample_z %in% merged_seurat_obj$subsample_z[merged_seurat_obj$gene_of_interest == gene_value])) {
      selected_genes <- unique(c(selected_genes, gene_value))
    }
  }
  
  # Subset the Seurat object based on selected genes
  merged_seurat_obj_filtered <- merged_seurat_obj[, merged_seurat_obj$gene_of_interest %in% selected_genes]
  
  # Print the number of selected genes
  cat("Number of selected genes:", length(selected_genes), "\n")
  
  # Optionally, you can assign the filtered Seurat object back to 'merged_seurat_obj'
  # merged_seurat_obj <- merged_seurat_obj_filtered
  
} else {
  cat("The 'subsample_z' and/or 'gene_of_interest' columns are not found in metadata.")
}



saveRDS(merged_seurat_obj_filtered, "C:/Users/zvjohns/Desktop/projects/clownfish/rds/clownfish_OB_signal_qc_new_cluster_names_goi_filtered_all_indivs.rds")

