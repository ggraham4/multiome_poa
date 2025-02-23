library(NeuronChat)
library(patchwork)
library(stringr)

#gene_info <- read.table("E:/GT_laptop_transfer_062923/Z_drive/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
gene_info <- read.csv("C:/Users/zvjohns/Desktop/projects/clownfish/csv/input/percula_ortho_mmusculus_coltan_good.csv") 
#gene_info <- read.table("Z:/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 


combined <- readRDS("C:/Users/zvjohns/Desktop/projects/clownfish/rds/clownfish_OB_signal_qc_new_cluster_names_poa_neuropeptide_steroid_filtered_all_indivs.rds")
#combined <- merged_seurat_obj_filtered
#combined <- readRDS("clownfish_cluster2_neuropeptide_poa_cellchat.rds")
#combined <- readRDS("clownfish_cluster2_neuropeptide_poa_cellchat_alt_gnrh.rds")

`%notin%` <- Negate(`%in%`)

combined <- subset(combined, subset = clusters49 %notin% c(34,39,45,46,48))

#combined <- subset(combined, subset = clusters49 %notin% c(34,35, 39,45,46,48))
#combined$cellchat_real_labels <- str_c("cluster_",combined$new_values)
combined$cellchat_real_labels <- combined$gene_of_interest
#combined$cellchat_real_labels <- str_c("cluster_",combined$cellchat_label)

#combined <- readRDS("C:/Users/zjohnson37/Desktop/bb/bb_demux_012422.rds")
#cluster9_Glut <- subset(combined, subset = seuratclusters15 == 1)
#cluster9_Glut$cellchat_real_labels <- str_c("cluster_9_Glut")

#metadata <- data.frame(combined@meta.data)
#metadata$downsample <- 0

clusters <- unique(combined$cellchat_real_labels)
trials <- unique(combined$subsample_z)

colnames <- as.vector("trial_id")
for(i in 1:length(clusters)){
  for(j in 1:length(clusters)){
    colnames <- append(colnames, paste0(clusters[i],"_",clusters[j]))
  }
}

#colnames <- colnames[1:length(colnames)] #remove trial_id as first column because this one is across all trials pooled

out = data.frame(matrix(nrow=12,ncol=length(colnames)))

colnames[2:length(colnames)] <- sort(colnames[2:length(colnames)])

colnames(out) <- as.vector(colnames)




#for(z in 1:length(trials)){
for(z in 1:length(trials)){

set.seed(z)

paste0(print(trials[z]))

#setwd("/storage/coda1/p-js585/0/zjohnson37/cellchat/")

#gene_info <- read.table("gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 

#combined <- readRDS("bb_demux_012422.rds")

#z = 6
#for(z in 1:1000){

#gene_info <- read.csv("percula_ortho_mmusculus_coltan_good.csv") 

combined_trial <- subset(combined, subset = subsample_z == trials[z])
combined_trial$cell_barcode <- colnames(combined_trial)
#combined_trial <- subset(combined_trial, subset = seuratclusters15 %in% c(0,1,4))
cell_ids <- combined_trial$cell_barcode

data.input <- data.frame(combined_trial@assays$RNA@data)
data.input$new_mouse_column = gene_info$mmusculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
data.input_NA_rm <- data.input[complete.cases(data.input),]
#rownames(data.input_NA_rm) <- data.input_NA_rm$new_mouse_column



#combined_trial_trial$cell_barcode <- colnames(combined_trial_trial)
#combined_trial_trial <- subset(combined_trial_trial, subset = seuratclusters15 %in% c(0,1,4))
#cell_ids <- combined_trial_trial$cell_barcode


#data.input <- data.frame(combined_trial_trial@assays$RNA@data)
#data.input$clownfish_gene <- rownames(data.input)
#data.input$new_mouse_column = gene_info$musculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
#data.input_NA_rm <- data.input[complete.cases(data.input$new_mouse_column),]

n <- ncol(data.input_NA_rm)-1

data.input_NA_rm$rowsums <- rowSums(data.input_NA_rm[,1:n])
data.input_NA_rm_sort <- data.input_NA_rm[order(-data.input_NA_rm$rowsums),]
data.input_NA_rm_sort_no_zero <- data.input_NA_rm_sort[which(data.input_NA_rm_sort$rowsums!=0),]

#head(data.input_NA_rm_sort$rowsums)[1:20]
#length(unique(data.input_NA_rm_sort_no_zero$rowsums))

data.input_trimmed <- data.input_NA_rm_sort_no_zero[!duplicated(data.input_NA_rm_sort_no_zero$new_mouse_column),]
rownames(data.input_trimmed) <- data.input_trimmed$new_mouse_column
data.input_trimmed$new_mouse_column <- NULL
data.input_trimmed$rowsums <- NULL
#class(data.input_trimmed)
#nrow(data.input_trimmed)
#ncol(data.input_trimmed)
data.input.matrix <- as.matrix(data.input_trimmed)
#nrow(data.input.matrix)
#ncol(data.input.matrix)

### now create your cell chat object
meta_temp <- combined_trial@meta.data
#meta_temp$cellchat_real_labels

#cellchat <- createCellChat(object = data.input.matrix, meta = meta, group.by = "good_names")

cellchat <- createNeuronChat(object = data.input.matrix, meta = meta_temp, group.by = as.vector(combined_trial$cellchat_real_labels))
cellchat <- run_NeuronChat(cellchat,M=100)
net_aggregated_cellchat <- net_aggregation(cellchat@net,method = 'weight')

out[z,1] <- trials[z]
out[z,2:length(colnames)] <- as.vector(t(as.matrix(net_aggregated_cellchat)))



write.csv(out, "C:/Users/zvjohns/Desktop/projects/clownfish/cellchat_weights_real_poa_neuro_clusters_by_trial_id_connectivity_100423.csv")

}



######## NOW FOR PARENT CLUSTERS



#gene_info <- read.table("E:/GT_laptop_transfer_062923/Z_drive/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
gene_info <- read.csv("E:/GT_laptop_transfer_062923/Z_drive/clownfish/percula_ortho_mmusculus_coltan_good.csv") 
#gene_info <- read.table("Z:/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 


combined <- readRDS("E:/GT_laptop_transfer_062923/Z_drive/clownfish/clownfish_OB_signal_qc_new_cluster_names.rds")
#combined <- readRDS("clownfish_cluster2_neuropeptide_poa_cellchat.rds")
#combined <- readRDS("clownfish_cluster2_neuropeptide_poa_cellchat_alt_gnrh.rds")


`%notin%` <- Negate(`%in%`)

table(combined$clusters20_new_name, combined$subsample_z)

combined <- subset(combined, subset = clusters49 %notin% c(34,39,45,46,48))

#combined <- subset(combined, subset = clusters49 %notin% c(34,35, 39,45,46,48))
combined$cellchat_real_labels <- str_c("cluster_",combined$clusters20_new_name)
#combined$cellchat_real_labels <- str_c("cluster_",combined$cellchat_label)

#combined <- readRDS("C:/Users/zjohnson37/Desktop/bb/bb_demux_012422.rds")
#cluster9_Glut <- subset(combined, subset = seuratclusters15 == 1)
#cluster9_Glut$cellchat_real_labels <- str_c("cluster_9_Glut")

#metadata <- data.frame(combined@meta.data)
#metadata$downsample <- 0

clusters <- unique(combined$cellchat_real_labels)
trials <- unique(combined$subsample_z)

colnames <- as.vector("trial_id")
for(i in 1:length(clusters)){
  for(j in 1:length(clusters)){
    colnames <- append(colnames, paste0(clusters[i],"_",clusters[j]))
  }
}

#colnames <- colnames[1:length(colnames)] #remove trial_id as first column because this one is across all trials pooled

out = data.frame(matrix(nrow=12,ncol=length(colnames)))

colnames[2:length(colnames)] <- sort(colnames[2:length(colnames)])

colnames(out) <- as.vector(colnames)


for(z in 1:length(trials)){
  
  set.seed(z)
  
  paste0(print(trials[z]))
  
  #setwd("/storage/coda1/p-js585/0/zjohnson37/cellchat/")
  
  #gene_info <- read.table("gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
  
  #combined <- readRDS("bb_demux_012422.rds")
  
  #z = 6
  #for(z in 1:1000){
  
  #gene_info <- read.csv("percula_ortho_mmusculus_coltan_good.csv") 
  
  combined_trial <- subset(combined, subset = subsample_z == trials[z])
  combined_trial$cell_barcode <- colnames(combined_trial)
  #combined_trial <- subset(combined_trial, subset = seuratclusters15 %in% c(0,1,4))
  cell_ids <- combined_trial$cell_barcode
  
  data.input <- data.frame(combined_trial@assays$RNA@data)
  data.input$new_mouse_column = gene_info$mmusculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
  data.input_NA_rm <- data.input[complete.cases(data.input),]
  #rownames(data.input_NA_rm) <- data.input_NA_rm$new_mouse_column
  
  
  
  #combined_trial_trial$cell_barcode <- colnames(combined_trial_trial)
  #combined_trial_trial <- subset(combined_trial_trial, subset = seuratclusters15 %in% c(0,1,4))
  #cell_ids <- combined_trial_trial$cell_barcode
  
  
  #data.input <- data.frame(combined_trial_trial@assays$RNA@data)
  #data.input$clownfish_gene <- rownames(data.input)
  #data.input$new_mouse_column = gene_info$musculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
  #data.input_NA_rm <- data.input[complete.cases(data.input$new_mouse_column),]
  
  n <- ncol(data.input_NA_rm)-1
  
  data.input_NA_rm$rowsums <- rowSums(data.input_NA_rm[,1:n])
  data.input_NA_rm_sort <- data.input_NA_rm[order(-data.input_NA_rm$rowsums),]
  data.input_NA_rm_sort_no_zero <- data.input_NA_rm_sort[which(data.input_NA_rm_sort$rowsums!=0),]
  
  #head(data.input_NA_rm_sort$rowsums)[1:20]
  #length(unique(data.input_NA_rm_sort_no_zero$rowsums))
  
  data.input_trimmed <- data.input_NA_rm_sort_no_zero[!duplicated(data.input_NA_rm_sort_no_zero$new_mouse_column),]
  rownames(data.input_trimmed) <- data.input_trimmed$new_mouse_column
  data.input_trimmed$new_mouse_column <- NULL
  data.input_trimmed$rowsums <- NULL
  #class(data.input_trimmed)
  #nrow(data.input_trimmed)
  #ncol(data.input_trimmed)
  data.input.matrix <- as.matrix(data.input_trimmed)
  #nrow(data.input.matrix)
  #ncol(data.input.matrix)
  
  ### now create your cell chat object
  meta_temp <- combined_trial@meta.data
  #meta_temp$cellchat_real_labels
  
  #cellchat <- createCellChat(object = data.input.matrix, meta = meta, group.by = "good_names")
  
  cellchat <- createNeuronChat(object = data.input.matrix, meta = meta_temp, group.by = as.vector(combined_trial$cellchat_real_labels))
  cellchat <- run_NeuronChat(cellchat,M=100)
  net_aggregated_cellchat <- net_aggregation(cellchat@net,method = 'weight')
  
  out[z,1] <- trials[z]
  out[z,2:length(colnames)] <- as.vector(t(as.matrix(net_aggregated_cellchat)))
  
  
  
  write.csv(out, "C:/Users/zvjohns/Desktop/projects/clownfish/cellchat_weights_real_parent_clusters_by_trial_id_connectivity_100423.csv")
  
}










######## NOW FOR NEUROPEPTIDE CLUSTERS



#gene_info <- read.table("E:/GT_laptop_transfer_062923/Z_drive/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
gene_info <- read.csv("C:/Users/zvjohns/Desktop/projects/clownfish/csv/input/percula_ortho_mmusculus_coltan_good.csv") 
#gene_info <- read.table("Z:/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 

combined <- readRDS("C:/Users/zvjohns/Desktop/projects/clownfish/rds/clownfish_OB_signal_qc_new_cluster_names_poa_neuropeptide_steroid_filtered_all_indivs.rds")
#combined <- readRDS("clownfish_cluster2_neuropeptide_poa_cellchat.rds")
#combined <- readRDS("clownfish_cluster2_neuropeptide_poa_cellchat_alt_gnrh.rds")

combined$cellchat_label <- combined$gene_of_interest
table(combined$cellchat_label, combined$subsample_z)

`%notin%` <- Negate(`%in%`)

combined <- subset(combined, subset = clusters49 %notin% c(34,39,45,46,48))

#combined <- subset(combined, subset = clusters49 %notin% c(34,35, 39,45,46,48))
#combined$cellchat_real_labels <- str_c("cluster_",combined$clusters20_new_name)
#combined$cellchat_real_labels <- str_c("cluster_",combined$cellchat_label)

#combined <- readRDS("C:/Users/zjohnson37/Desktop/bb/bb_demux_012422.rds")
#cluster9_Glut <- subset(combined, subset = seuratclusters15 == 1)
#cluster9_Glut$cellchat_real_labels <- str_c("cluster_9_Glut")

#metadata <- data.frame(combined@meta.data)
#metadata$downsample <- 0

clusters <- unique(combined$cellchat_label)
trials <- unique(combined$subsample_z)

colnames <- as.vector("trial_id")
for(i in 1:length(clusters)){
  for(j in 1:length(clusters)){
    colnames <- append(colnames, paste0(clusters[i],"_",clusters[j]))
  }
}

#colnames <- colnames[1:length(colnames)] #remove trial_id as first column because this one is across all trials pooled

out = data.frame(matrix(nrow=12,ncol=length(colnames)))

colnames[2:length(colnames)] <- sort(colnames[2:length(colnames)])

colnames(out) <- as.vector(colnames)


for(z in 1:length(trials)){
  
  set.seed(z)
  
  paste0(print(trials[z]))
  
  #setwd("/storage/coda1/p-js585/0/zjohnson37/cellchat/")
  
  #gene_info <- read.table("gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
  
  #combined <- readRDS("bb_demux_012422.rds")
  
  #z = 6
  #for(z in 1:1000){
  
  #gene_info <- read.csv("percula_ortho_mmusculus_coltan_good.csv") 
  
  combined_trial <- subset(combined, subset = subsample_z == trials[z])
  combined_trial$cell_barcode <- colnames(combined_trial)
  #combined_trial <- subset(combined_trial, subset = seuratclusters15 %in% c(0,1,4))
  cell_ids <- combined_trial$cell_barcode
  
  
  data.input <- data.frame(as.matrix(combined_trial@assays$RNA@data))
  data.input$new_mouse_column = gene_info$mmusculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
  data.input_NA_rm <- data.input[complete.cases(data.input),]
  #rownames(data.input_NA_rm) <- data.input_NA_rm$new_mouse_column
  
  
  
  #combined_trial_trial$cell_barcode <- colnames(combined_trial_trial)
  #combined_trial_trial <- subset(combined_trial_trial, subset = seuratclusters15 %in% c(0,1,4))
  #cell_ids <- combined_trial_trial$cell_barcode
  
  
  #data.input <- data.frame(combined_trial_trial@assays$RNA@data)
  #data.input$clownfish_gene <- rownames(data.input)
  #data.input$new_mouse_column = gene_info$musculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
  #data.input_NA_rm <- data.input[complete.cases(data.input$new_mouse_column),]
  
  n <- ncol(data.input_NA_rm)-1
  
  data.input_NA_rm$rowsums <- rowSums(data.input_NA_rm[,1:n])
  data.input_NA_rm_sort <- data.input_NA_rm[order(-data.input_NA_rm$rowsums),]
  data.input_NA_rm_sort_no_zero <- data.input_NA_rm_sort[which(data.input_NA_rm_sort$rowsums!=0),]
  
  #head(data.input_NA_rm_sort$rowsums)[1:20]
  #length(unique(data.input_NA_rm_sort_no_zero$rowsums))
  
  data.input_trimmed <- data.input_NA_rm_sort_no_zero[!duplicated(data.input_NA_rm_sort_no_zero$new_mouse_column),]
  rownames(data.input_trimmed) <- data.input_trimmed$new_mouse_column
  data.input_trimmed$new_mouse_column <- NULL
  data.input_trimmed$rowsums <- NULL
  #class(data.input_trimmed)
  #nrow(data.input_trimmed)
  #ncol(data.input_trimmed)
  data.input.matrix <- as.matrix(data.input_trimmed)
  #nrow(data.input.matrix)
  #ncol(data.input.matrix)
  
  ### now create your cell chat object
  meta_temp <- combined_trial@meta.data
  #meta_temp$cellchat_real_labels
  
  #cellchat <- createCellChat(object = data.input.matrix, meta = meta, group.by = "good_names")
  
  cellchat <- createNeuronChat(object = data.input.matrix, meta = meta_temp, group.by = as.vector(combined_trial$cellchat_label))
  cellchat <- run_NeuronChat(cellchat,M=100)
  net_aggregated_cellchat <- net_aggregation(cellchat@net,method = 'weight')
  
  out[z,1] <- trials[z]
  out[z,2:length(colnames)] <- as.vector(t(as.matrix(net_aggregated_cellchat)))
  
  
  
  write.csv(out, "C:/Users/zvjohns/Desktop/projects/clownfish/cellchat_weights_real_poa_goi_parent_clusters_by_trial_id_connectivity_100423.csv")
  
}














######## NOW FOR WHOLE OBJECT INCLUDING NEUROPEPTIDE CLUSTERS



#gene_info <- read.table("E:/GT_laptop_transfer_062923/Z_drive/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
gene_info <- read.csv("C:/Users/zvjohns/Desktop/projects/clownfish/csv/input/percula_ortho_mmusculus_coltan_good.csv") 
#gene_info <- read.table("Z:/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 

combined <- readRDS("C:/Users/zvjohns/Desktop/projects/clownfish/rds/clownfish_OB_signal_qc_new_cluster_names_poa_neuropeptide_steroid_filtered_all_indivs.rds")
#combined <- readRDS("clownfish_cluster2_neuropeptide_poa_cellchat.rds")
#combined <- readRDS("clownfish_cluster2_neuropeptide_poa_cellchat_alt_gnrh.rds")

unique(combined$cellchat_label)
combined$cellchat_label <- combined$gene_of_interest
table(combined$cellchat_label, combined$subsample_z)

`%notin%` <- Negate(`%in%`)

combined <- subset(combined, subset = clusters49 %notin% c(34,39,45,46,48))

#combined <- subset(combined, subset = clusters49 %notin% c(34,35, 39,45,46,48))
#combined$cellchat_real_labels <- str_c("cluster_",combined$clusters20_new_name)
#combined$cellchat_real_labels <- str_c("cluster_",combined$cellchat_label)

#combined <- readRDS("C:/Users/zjohnson37/Desktop/bb/bb_demux_012422.rds")
#cluster9_Glut <- subset(combined, subset = seuratclusters15 == 1)
#cluster9_Glut$cellchat_real_labels <- str_c("cluster_9_Glut")

#metadata <- data.frame(combined@meta.data)
#metadata$downsample <- 0

clusters <- unique(combined$cellchat_label)
trials <- unique(combined$subsample_z)

colnames <- as.vector("trial_id")
for(i in 1:length(clusters)){
  for(j in 1:length(clusters)){
    colnames <- append(colnames, paste0(clusters[i],"_",clusters[j]))
  }
}

#colnames <- colnames[1:length(colnames)] #remove trial_id as first column because this one is across all trials pooled

out = data.frame(matrix(nrow=12,ncol=length(colnames)))

colnames[2:length(colnames)] <- sort(colnames[2:length(colnames)])

colnames(out) <- as.vector(colnames)


for(z in 1:length(trials)){
  
  set.seed(z)
  
  paste0(print(trials[z]))
  
  #setwd("/storage/coda1/p-js585/0/zjohnson37/cellchat/")
  
  #gene_info <- read.table("gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
  
  #combined <- readRDS("bb_demux_012422.rds")
  
  #z = 6
  #for(z in 1:1000){
  
  #gene_info <- read.csv("percula_ortho_mmusculus_coltan_good.csv") 
  
  combined_trial <- subset(combined, subset = subsample_z == trials[z])
  combined_trial$cell_barcode <- colnames(combined_trial)
  #combined_trial <- subset(combined_trial, subset = seuratclusters15 %in% c(0,1,4))
  cell_ids <- combined_trial$cell_barcode
  
  
  data.input <- data.frame(as.matrix(combined_trial@assays$RNA@data))
  data.input$new_mouse_column = gene_info$mmusculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
  data.input_NA_rm <- data.input[complete.cases(data.input),]
  #rownames(data.input_NA_rm) <- data.input_NA_rm$new_mouse_column
  
  
  
  #combined_trial_trial$cell_barcode <- colnames(combined_trial_trial)
  #combined_trial_trial <- subset(combined_trial_trial, subset = seuratclusters15 %in% c(0,1,4))
  #cell_ids <- combined_trial_trial$cell_barcode
  
  
  #data.input <- data.frame(combined_trial_trial@assays$RNA@data)
  #data.input$clownfish_gene <- rownames(data.input)
  #data.input$new_mouse_column = gene_info$musculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
  #data.input_NA_rm <- data.input[complete.cases(data.input$new_mouse_column),]
  
  n <- ncol(data.input_NA_rm)-1
  
  data.input_NA_rm$rowsums <- rowSums(data.input_NA_rm[,1:n])
  data.input_NA_rm_sort <- data.input_NA_rm[order(-data.input_NA_rm$rowsums),]
  data.input_NA_rm_sort_no_zero <- data.input_NA_rm_sort[which(data.input_NA_rm_sort$rowsums!=0),]
  
  #head(data.input_NA_rm_sort$rowsums)[1:20]
  #length(unique(data.input_NA_rm_sort_no_zero$rowsums))
  
  data.input_trimmed <- data.input_NA_rm_sort_no_zero[!duplicated(data.input_NA_rm_sort_no_zero$new_mouse_column),]
  rownames(data.input_trimmed) <- data.input_trimmed$new_mouse_column
  data.input_trimmed$new_mouse_column <- NULL
  data.input_trimmed$rowsums <- NULL
  #class(data.input_trimmed)
  #nrow(data.input_trimmed)
  #ncol(data.input_trimmed)
  data.input.matrix <- as.matrix(data.input_trimmed)
  #nrow(data.input.matrix)
  #ncol(data.input.matrix)
  
  ### now create your cell chat object
  meta_temp <- combined_trial@meta.data
  #meta_temp$cellchat_real_labels
  
  #cellchat <- createCellChat(object = data.input.matrix, meta = meta, group.by = "good_names")
  
  cellchat <- createNeuronChat(object = data.input.matrix, meta = meta_temp, group.by = as.vector(combined_trial$cellchat_label))
  cellchat <- run_NeuronChat(cellchat,M=100)
  net_aggregated_cellchat <- net_aggregation(cellchat@net,method = 'weight')
  
  out[z,1] <- trials[z]
  out[z,2:length(colnames)] <- as.vector(t(as.matrix(net_aggregated_cellchat)))
  
  
  
  write.csv(out, "C:/Users/zvjohns/Desktop/projects/clownfish/cellchat_weights_real_poa_goi_parent_clusters_by_trial_id_connectivity_100923.csv")
  
}






















######## NOW FOR GOI CLUSTERS



#gene_info <- read.table("E:/GT_laptop_transfer_062923/Z_drive/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
gene_info <- read.csv("C:/Users/zvjohns/Desktop/projects/clownfish/csv/input/percula_ortho_mmusculus_coltan_good.csv") 
#gene_info <- read.table("Z:/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 

combined <- readRDS("C:/Users/zvjohns/Desktop/projects/clownfish/rds/clownfish_OB_signal_qc_new_cluster_names_goi_filtered_all_indivs.rds")
#combined <- readRDS("clownfish_cluster2_neuropeptide_poa_cellchat.rds")
#combined <- readRDS("clownfish_cluster2_neuropeptide_poa_cellchat_alt_gnrh.rds")

combined$cellchat_label <- combined$gene_of_interest
table(combined$cellchat_label, combined$subsample_z)

`%notin%` <- Negate(`%in%`)

#combined <- subset(combined, subset = clusters49 %notin% c(34,39,45,46,48))

#combined <- subset(combined, subset = clusters49 %notin% c(34,35, 39,45,46,48))
#combined$cellchat_real_labels <- str_c("cluster_",combined$clusters20_new_name)
#combined$cellchat_real_labels <- str_c("cluster_",combined$cellchat_label)

#combined <- readRDS("C:/Users/zjohnson37/Desktop/bb/bb_demux_012422.rds")
#cluster9_Glut <- subset(combined, subset = seuratclusters15 == 1)
#cluster9_Glut$cellchat_real_labels <- str_c("cluster_9_Glut")

#metadata <- data.frame(combined@meta.data)
#metadata$downsample <- 0

clusters <- unique(combined$cellchat_label)
trials <- unique(combined$subsample_z)

colnames <- as.vector("trial_id")
for(i in 1:length(clusters)){
  for(j in 1:length(clusters)){
    colnames <- append(colnames, paste0(clusters[i],"_",clusters[j]))
  }
}

#colnames <- colnames[1:length(colnames)] #remove trial_id as first column because this one is across all trials pooled

out = data.frame(matrix(nrow=12,ncol=length(colnames)))

colnames[2:length(colnames)] <- sort(colnames[2:length(colnames)])

colnames(out) <- as.vector(colnames)


for(z in 1:length(trials)){
  
  set.seed(z)
  
  paste0(print(trials[z]))
  
  #setwd("/storage/coda1/p-js585/0/zjohnson37/cellchat/")
  
  #gene_info <- read.table("gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
  
  #combined <- readRDS("bb_demux_012422.rds")
  
  #z = 6
  #for(z in 1:1000){
  
  #gene_info <- read.csv("percula_ortho_mmusculus_coltan_good.csv") 
  
  combined_trial <- subset(combined, subset = subsample_z == trials[z])
  combined_trial$cell_barcode <- colnames(combined_trial)
  #combined_trial <- subset(combined_trial, subset = seuratclusters15 %in% c(0,1,4))
  cell_ids <- combined_trial$cell_barcode
  
  data.input <- data.frame(combined_trial@assays$RNA@data)
  data.input$new_mouse_column = gene_info$mmusculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
  data.input_NA_rm <- data.input[complete.cases(data.input),]
  #rownames(data.input_NA_rm) <- data.input_NA_rm$new_mouse_column
  
  
  
  #combined_trial_trial$cell_barcode <- colnames(combined_trial_trial)
  #combined_trial_trial <- subset(combined_trial_trial, subset = seuratclusters15 %in% c(0,1,4))
  #cell_ids <- combined_trial_trial$cell_barcode
  
  
  #data.input <- data.frame(combined_trial_trial@assays$RNA@data)
  #data.input$clownfish_gene <- rownames(data.input)
  #data.input$new_mouse_column = gene_info$musculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
  #data.input_NA_rm <- data.input[complete.cases(data.input$new_mouse_column),]
  
  n <- ncol(data.input_NA_rm)-1
  
  data.input_NA_rm$rowsums <- rowSums(data.input_NA_rm[,1:n])
  data.input_NA_rm_sort <- data.input_NA_rm[order(-data.input_NA_rm$rowsums),]
  data.input_NA_rm_sort_no_zero <- data.input_NA_rm_sort[which(data.input_NA_rm_sort$rowsums!=0),]
  
  #head(data.input_NA_rm_sort$rowsums)[1:20]
  #length(unique(data.input_NA_rm_sort_no_zero$rowsums))
  
  data.input_trimmed <- data.input_NA_rm_sort_no_zero[!duplicated(data.input_NA_rm_sort_no_zero$new_mouse_column),]
  rownames(data.input_trimmed) <- data.input_trimmed$new_mouse_column
  data.input_trimmed$new_mouse_column <- NULL
  data.input_trimmed$rowsums <- NULL
  #class(data.input_trimmed)
  #nrow(data.input_trimmed)
  #ncol(data.input_trimmed)
  data.input.matrix <- as.matrix(data.input_trimmed)
  #nrow(data.input.matrix)
  #ncol(data.input.matrix)
  
  ### now create your cell chat object
  meta_temp <- combined_trial@meta.data
  #meta_temp$cellchat_real_labels
  
  #cellchat <- createCellChat(object = data.input.matrix, meta = meta, group.by = "good_names")
  
  cellchat <- createNeuronChat(object = data.input.matrix, meta = meta_temp, group.by = as.vector(combined_trial$cellchat_label))
  cellchat <- run_NeuronChat(cellchat,M=100)
  net_aggregated_cellchat <- net_aggregation(cellchat@net,method = 'weight')
  
  out[z,1] <- trials[z]
  out[z,2:length(colnames)] <- as.vector(t(as.matrix(net_aggregated_cellchat)))
  
  
  
  write.csv(out, "C:/Users/zvjohns/Desktop/projects/clownfish/cellchat_weights_real_goi_clusters_by_trial_id_connectivity_100423.csv")
  
}




























colnames(net_aggregated_cellchat)




net_aggregated_cellchat

as.vector(t(as.matrix(net_aggregated_cellchat)))


par(mfrow=c(1,2))
netVisual_circle_neuron(net_aggregated_cellchat,vertex.label.cex = 1)
netVisual_chord_neuron(cellchat,method = 'weight',lab.cex = 1)
heatmap_aggregated(cellchat, method='weight')



g1 <- rankNet_Neuron(cellchat,slot.name = "net",measure = c("weight"),mode='single',font.size = 5) 
g2 <- rankNet_Neuron(cellchat,slot.name = "net",measure = c("count"),mode='single',font.size = 5) 
g1+g2

cellchat@net
cellchat <- identifyCommunicationPatterns_Neuron(cellchat, slot.name = "net", pattern = c("outgoing"), k=4,height = 12)
cellchat <- identifyCommunicationPatterns_Neuron(cellchat, slot.name = "net", pattern = c("incoming"), k=4,height = 12)

library(ggalluvial)

?netAnalysis_river_Neuron
netAnalysis_river_Neuron(cellchat,slot.name = "net", pattern = c("outgoing"),font.size = 2.5,cutoff.1 = 0.5,cutoff.2=0.5)
netAnalysis_river_Neuron(cellchat,slot.name = "net", pattern = c("incoming"),font.size = 2.5,cutoff.1 = 0.5,cutoff.2=0.5)

cellchat



cellchat@data.signaling
?net_aggregation
?netVisual_circle_neuron
?createNeuronChat
?heatmap_aggregated

nrow(cellchat@ligand.abundance)
nrow(cellchat@target.abundance)
nrow(cellchat@net_analysis)
cellchat@net_analysis$pattern$incoming
cellchat@fc
cellchat@info
cellchat@idents
cellchat@net$Vip_Vipr2
cellchat@

test <- data.frame(cellchat@fc)
test <- data.frame(cellchat@info)
rownames(test)





cellchat <- addMeta(cellchat, meta = meta_temp)
cellchat <- setIdent(cellchat, ident.use = "cellchat_real_labels") # set "labels" as default cell identity

levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB # simply use the default CellChatDBt
which(CellChatDB.use[["interaction"]]$ligand == "H2-BI") # 1887
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1887,]
which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps") #1900
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1900,]
cellchat@DB <- CellChatDB.use
CellChatDB$geneInfo

#showDatabaseCategory(CellChatDB)


cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
#meta$cellchat_real_labels = droplevels(meta$cellchat_real_labels, exclude = setdiff(levels(meta$cellchat_real_labels),unique(meta$cellchat_real_labels)))

cellchat <- computeCommunProb(cellchat, type =  "truncatedMean", trim = 0.1, population.size = FALSE)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#cellchat <- filterCommunication(cellchat, min.cells = 10)
groupSize <- as.numeric(table(cellchat@idents))



#out[z,1] <- trials[z]
head(cellchat@net_analysis$pattern$outgoing$pattern$
       out[1,(1:length(as.vector(t(as.matrix(cellchat@net$weight)))))] <- as.vector(t(as.matrix(cellchat@net$weight)))
     
     write.csv(out, "cellchat_weights_real_neuropeptide_connectivity_030623_anenomefish_gnrh_alt.csv")
     
     hist(as.vector(as.numeric(out[1,])), breaks = 18)
     
     t_mat <- data.frame(t(as.matrix(out)))
     
     t_mat_high <- subset(t_mat, t_mat$X1 > 0.5)
     
     
     #}