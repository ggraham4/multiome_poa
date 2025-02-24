### Trying to replicate zack's code
### neuronchat_by_individual_100423

library(CellChat)
library(patchwork)
library(stringr)


##gene info 
##I dont have this file but I think I know what it is so I can make that
#library(biomaRt)

##make biomart
#ocellaris_mart <- useEnsembl(biomart = 'genes',
#                             dataset = 'aocellaris_gene_ensembl')


##because ensembl id and gene name are on separate pages, we need to do this in two steps
#ocellaris_to_mmusculus1 <- getBM(mart = ocellaris_mart,
#                                attributes = c('ensembl_transcript_id',
#                                               'mmusculus_homolog_associated_gene_name')
#)

#ocellaris_to_mmusculus2 <- getBM(mart = ocellaris_mart,
#                                 attributes = c('ensembl_transcript_id',
#                                                'entrezgene_accession')
#)

#join them
#combined_data <- ocellaris_to_mmusculus1%>%
  #right_join(ocellaris_to_mmusculus2, by = 'ensembl_transcript_id')%>%
  #dplyr::select(-1)%>%
 # distinct()
#colnames(combined_data)= c('mmusculus_name', 'aocellaris_name')

##save
#write.csv(combined_data, 'C:/Users/Gabe/Desktop/multiome_poa/Reference/mmusculus_to_aocellaris.csv')

#only have to do all that once
combined_data <- read.csv('Reference/mmusculus_to_aocellaris.csv')

### read in object
obj <-readRDS("/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds")

clusters <- unique(obj$harmony.wnn_res0.4_clusters)

### im not really sure what zack is doing in his code with trials, but I think I can ignore it
obj$cell_barcode <- colnames(obj)
cell_ids <- obj$cell_barcode

data_input = data.frame(obj@assays$RNA$data)

##ok this i need to really understand
data_input$new_mouse_column = combined_data$mmusculus_name[match(rownames(data_input), combined_data$aocellaris_name)]
## what I dont understand is how hes going from rows to columns
rownames(data_input) ##wait ok this makes sense im just so not used to doing it this way, he made a new column with the mouse names so he can replace the rownames I assuem

##make meta data object to add later 
meta_data <- as.data.frame(obj@meta.data)

## removes rows with missing mouse gene names
data.input_NA_rm <-  data_input[complete.cases(data_input),]

### remove rows with 0 expression
n <- ncol(data.input_NA_rm)-1
data.input_NA_rm$rowsums <- rowSums(data.input_NA_rm[,1:n])

data.input_NA_rm_sort <- data.input_NA_rm[order(-data.input_NA_rm$rowsums),]
data.input_NA_rm_sort_no_zero <- data.input_NA_rm_sort[which(data.input_NA_rm_sort$rowsums!=0),]

### remove duplicate mouse columns
## i'm kind of concerned about this though cause like ccka and cckb will be listed as cck in the mouse gene
## so how do you know which one you are keeping
## that being said, it seems really annoying to fix that, my guess is youd have to sum by common gene name 
## and then join by cell or something, I'm not sure

##zacks method
#data.input_trimmed_zack <- data.input_NA_rm_sort_no_zero[!duplicated(data.input_NA_rm_sort_no_zero$new_mouse_column),]


data_matrix <- data.input_NA_rm_sort_no_zero[, -which(colnames(data.input_NA_rm_sort_no_zero) == "new_mouse_column")]
## this cannot be right why is it so fucking slow
aggregated_data <- rowsum(data_matrix, group = data.input_NA_rm_sort_no_zero$new_mouse_column, reorder = FALSE)
data.input_trimmed <- as.data.frame(aggregated_data)
rownames(data.input_trimmed) <- rownames(aggregated_data)

### lmfao they produce the exact same result
#unique(data.input_trimmed_zack[,-52519] == data.input_trimmed)

### make a seurat object and save
library(Seurat)
new_obj <- CreateSeuratObject(counts = data.input_trimmed[-1, -which(colnames(data.input_trimmed) == "rowsums")],
                              project = "renamed")
rownames(new_obj) <- rownames(data.input_trimmed[-1, -which(colnames(data.input_trimmed) == "rowsums")])
col_names <- colnames(data.input_trimmed[-1, -which(colnames(data.input_trimmed) == "rowsums")])
colnames(new_obj) <- sub("\\.(\\d+)$", "-\\1", col_names)

new_obj <- JoinLayers(new_obj)

new_obj@assays$RNA$data <- new_obj@assays$RNA$counts

new_obj@meta.data <- meta_data
new_obj@reductions <- obj@reductions
Idents(new_obj) <- 'harmony.wnn_res0.4_clusters'

DimPlot(new_obj)
FeaturePlot(new_obj, 'Cck') 
FeaturePlot(new_obj, 'Cckbr' ) 

#saveRDS(new_obj, '/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA object mouse names.rds')
marks <- FindAllMarkers(new_obj)

marks_19 <- subset(marks, cluster ==19)
#write.csv(marks_19, 'Reference/mmusculus_cluster_19_markers.csv')

FeaturePlot(new_obj, 'Kiss1r' ) 
FeaturePlot(new_obj, 'Tacr3' ) 
FeaturePlot(new_obj, 'Tac2' ) 

