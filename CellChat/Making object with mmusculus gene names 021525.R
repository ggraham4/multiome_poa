library(tidyverse)
library(readxl)
library(Seurat)
`%notin%` <- Negate(`%in%`)


###read in data
obj <-readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")


# first pull out rownames for fish
rownames <- rownames(obj)

#make a biomart of a ocellaris vs mmusculus gene names and vice versa ####

#start aocellaris biomart
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")
att <-listAttributes(ensembl)

#use ensemble_transcript id as a bridge between entrezgene_accession and mmusuculus homolog
biomart_gene_name <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c('ensembl_transcript_id',
                   'entrezgene_accession'))
biomart_mmusculus <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c('ensembl_transcript_id',
                   'mmusculus_homolog_associated_gene_name'))

#merge them
joined <- right_join(biomart_gene_name, biomart_mmusculus, by ='ensembl_transcript_id')[,2:3]
#remove duplicates
joined2 <- joined%>%dplyr::distinct(
  .keep_all = T)

#rename
translated <- joined2

#create translator function
mus_to_ocellaris_translator <- function(mus_gene){
  mus_gene_lower <- tolower(mus_gene)
  ocellaris_gene <- translated$entrezgene_accession[
    translated$mmusculus_homolog_associated_gene_name == mus_gene]
  return(ocellaris_gene)
}

ocellaris_to_mus_translator <- function(ocellaris_gene){
  mmusculus_gene <- translated$mmusculus_homolog_associated_gene_name[
    translated$entrezgene_accession == ocellaris_gene]
  return(mmusculus_gene)
}
ocellaris_to_mus_translator('haus2')
mus_to_ocellaris_translator('Haus2')
###ok they work

# loop over fish genes and replace with corresponding mouse gene
#chat gpt code cause this is rly annoying
# Get original row names (gene names)
gene_names <- rownames(obj@assays$RNA$data)

# Initialize vectors for new gene names and indices to remove
converted_names <- character(length(gene_names))
remove_indices <- c()

# Iterate over each gene name
for (i in seq_along(gene_names)) {
  print(paste0(i, " of ", length(gene_names)))  # Progress tracking
  
  converted_name <- ocellaris_to_mus_translator(gene_names[i])
  
  # If translation fails or returns empty, mark for removal
  if (is.null(converted_name) || length(converted_name) == 0 || all(converted_name == "")) {
    remove_indices <- c(remove_indices, i)
  } else {
    # If it returns a list, take the first item
    converted_names[i] <- if (is.list(converted_name)) converted_name[[1]] else converted_name
  }
}

# Remove invalid gene indices
filtered_data <- obj@assays$RNA$data[-remove_indices, ]

# Remove corresponding gene names
valid_names <- converted_names[-remove_indices]

# Ensure unique row names and remove any remaining empty or NA names
valid_names <- make.unique(valid_names)
valid_names <- valid_names[valid_names != "" & !is.na(valid_names)]  # Remove empty or NA names

# Ensure filtered_data and valid_names have the same length
filtered_data <- filtered_data[seq_along(valid_names), ]  # Trim matrix to match row names

# Assign unique rownames to filtered_data
rownames(filtered_data) <- valid_names

# Create a new Seurat object
new_obj <- CreateSeuratObject(counts = filtered_data, project = "renamed")
new_obj <- JoinLayers(new_obj)
new_obj@assays$RNA$data <- new_obj@assays$RNA$counts
# Final check
print(dim(new_obj@assays$RNA$counts))
print(head(rownames(new_obj@assays$RNA@counts)))

new_obj@meta.data <- obj@meta.data
new_obj@reductions <- obj@reductions

Idents(new_obj) <- 'harmony.wnn_res0.4_clusters'

DimPlot(new_obj)
FeaturePlot(new_obj, 'Cck' ) #HOLY FUCK IT WORKS THANK YOU JESUS
FeaturePlot(new_obj, 'Cckbr' ) 

### SAVE THE HOE WE CAN DO CELL CHAT TOMORROW!!!! 

saveRDS(new_obj, 'C:/Users/Gabe/Desktop/RNA object mouse names.rds')

