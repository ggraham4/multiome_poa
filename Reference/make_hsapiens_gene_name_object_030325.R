### Trying to replicate zack's code
### neuronchat_by_individual_100423

library(patchwork)
library(stringr)


##gene info 
##I dont have this file but I think I know what it is so I can make that
library(biomaRt)

##make biomart
ocellaris_mart <- useEnsembl(biomart = 'genes',
                             dataset = 'aocellaris_gene_ensembl')


att = attributes(ocellaris_mart)
##because ensembl id and gene name are on separate pages, we need to do this in two steps
ocellaris_to_hsapiens1 <- getBM(mart = ocellaris_mart,
                                attributes = c('ensembl_transcript_id',
                                               'hsapiens_homolog_associated_gene_name'))

ocellaris_to_hsapiens2 <- getBM(mart = ocellaris_mart,
                                attributes = c('ensembl_transcript_id',
                                               'entrezgene_accession')
)

#join them
combined_data <- ocellaris_to_hsapiens1%>%
  right_join(ocellaris_to_hsapiens2, by = 'ensembl_transcript_id')%>%
  dplyr::select(-1)%>%
  distinct()
colnames(combined_data)= c('hsapiens_name', 'aocellaris_name')

##manually fix some problems 
##kisspeptin names are fucked 
combined_data$hsapiens_name <- ifelse(combined_data$aocellaris_name %in% c('kiss2', 'kiss1'), 'KISS1',combined_data$hsapiens_name )
combined_data$hsapiens_name <- ifelse(combined_data$aocellaris_name %in% c('kiss1ra', 'kiss1rb'), 'KISS1R',combined_data$hsapiens_name )
combined_data$hsapiens_name <- ifelse(combined_data$aocellaris_name %in% c('galr1a'), 'GALR1',combined_data$hsapiens_name )

##somehow it gets oxt and avp backward, maybe I need to think of a better approach
combined_data$hsapiens_name <- ifelse(combined_data$aocellaris_name %in% c('oxt'), 'OXT',combined_data$hsapiens_name )
combined_data$hsapiens_name <- ifelse(combined_data$aocellaris_name %in% c('avp'), 'AVP',combined_data$hsapiens_name )
combined_data$hsapiens_name <- ifelse(combined_data$aocellaris_name %in% c('avpr2b', 'avpr2l'), 'AVPR2',combined_data$hsapiens_name )

#it seems like most of my genes have good coverage re the cellphone db database so I think I'm otherwise ok

##save
write.csv(combined_data, 'C:/Users/Gabe/Desktop/multiome_poa/Reference/hsapiens_to_aocellaris.csv')

#only have to do all that once
combined_data <- read.csv('Reference/hsapiens_to_aocellaris.csv')

### read in object
obj <-readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

clusters <- unique(obj$harmony.wnn_res0.4_clusters)

### im not really sure what zack is doing in his code with trials, but I think I can ignore it
obj$cell_barcode <- colnames(obj)
cell_ids <- obj$cell_barcode

data_input = data.frame(obj@assays$RNA$data)

##ok this i need to really understand
data_input$new_human_column = combined_data$hsapiens[match(rownames(data_input), combined_data$aocellaris_name)]
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

### remove duplicate human columns


data_matrix <- data.input_NA_rm_sort_no_zero[, -which(colnames(data.input_NA_rm_sort_no_zero) == "new_human_column")]
## this cannot be right why is it so fucking slow
aggregated_data <- rowsum(data_matrix, group = data.input_NA_rm_sort_no_zero$new_human_column, reorder = FALSE)
data.input_trimmed <- as.data.frame(aggregated_data)
rownames(data.input_trimmed) <- rownames(aggregated_data)

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
FeaturePlot(new_obj, 'CCK') 
FeaturePlot(new_obj, 'CCKBR' ) 

saveRDS(new_obj, 'C:/Users/Gabe/Desktop/RNA_object_human_names.rds')
marks <- FindAllMarkers(new_obj)

marks_19 <- subset(marks, cluster ==19)
#write.csv(marks_19, 'Reference/mmusculus_cluster_19_markers.csv')

FeaturePlot(new_obj, 'Kiss1r' ) 
FeaturePlot(new_obj, 'Tacr3' ) 
FeaturePlot(new_obj, 'Tac2' ) 

