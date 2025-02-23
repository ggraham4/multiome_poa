library(NeuronChat)
library(patchwork)
library(stringr)
library(Seurat)

#### NOW ALL CHILD CLUSTERS EXCEPT 13d + 47 13d neuropeptide populations

## first how i made the object 

data <- readRDS("C:/Users/zvjohns/Desktop/projects/clownfish/rds/anenomefish_object_clean.rds")
data$gene_of_interest <- data$child_cluster

neuronchat_object <- data

length(unique(neuronchat_object$gene_of_interest))

### now run neuronchat ----

#gene_info <- read.table("E:/GT_laptop_transfer_062923/Z_drive/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
gene_info <- read.csv("C:/Users/zvjohns/Desktop/projects/clownfish/csv/input/percula_ortho_mmusculus_coltan_good.csv") 
#gene_info <- read.table("Z:/clownfish/gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 

data$neuronchat_label <- data$gene_of_interest
table(data$neuronchat_label, data$subject)

`%notin%` <- Negate(`%in%`)

#data <- subset(data, subset = clusters49 %notin% c(34,39,45,46,48))

#data <- subset(data, subset = clusters49 %notin% c(34,35, 39,45,46,48))
#data$neuronchat_real_labels <- str_c("cluster_",data$clusters20_new_name)
#data$neuronchat_real_labels <- str_c("cluster_",data$neuronchat_label)

#data <- readRDS("C:/Users/zjohnson37/Desktop/bb/bb_demux_012422.rds")
#cluster9_Glut <- subset(data, subset = seuratclusters15 == 1)
#cluster9_Glut$neuronchat_real_labels <- str_c("cluster_9_Glut")

#metadata <- data.frame(data@meta.data)
#metadata$downsample <- 0

clusters <- unique(data$neuronchat_label)

colnames <- as.vector("placeholder")
for(i in 1:length(clusters)){
  for(j in 1:length(clusters)){
    colnames <- append(colnames, paste0(clusters[i],"_",clusters[j]))
  }
}

#colnames <- colnames[1:length(colnames)] #remove trial_id as first column because this one is across all trials pooled

out = data.frame(matrix(nrow=length(colnames)-1,ncol=1))

colnames <- sort(colnames[2:length(colnames)])

colnames(out) <- "weight"
rownames(out) <- colnames

  #setwd("/storage/coda1/p-js585/0/zjohnson37/neuronchat/")
  
  #gene_info <- read.table("gene_info.txt", sep="\t", header = T, stringsAsFactors = F) 
  
  #data <- readRDS("bb_demux_012422.rds")
  
  #z = 6
  #for(z in 1:1000){
  
  #gene_info <- read.csv("percula_ortho_mmusculus_coltan_good.csv") 
  
  #data_trial <- subset(data, subset = subsample_z == trials[z])
  data$cell_barcode <- colnames(data)
  #data_trial <- subset(data_trial, subset = seuratclusters15 %in% c(0,1,4))
  cell_ids <- data$cell_barcode
  
  data.input <- data.frame(data@assays$RNA@data)
  data.input$new_mouse_column = gene_info$mmusculus_homolog_associated_gene_name[match(rownames(data.input), gene_info$external_gene_name)]
  data.input_NA_rm <- data.input[complete.cases(data.input),]
  #rownames(data.input_NA_rm) <- data.input_NA_rm$new_mouse_column
  
  
  
  #data_trial$cell_barcode <- colnames(data_trial)
  #data_trial <- subset(data_trial, subset = seuratclusters15 %in% c(0,1,4))
  #cell_ids <- data_trial$cell_barcode
  
  
  #data.input <- data.frame(data_trial@assays$RNA@data)
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
  meta_temp <- data@meta.data
  #meta_temp$neuronchat_real_labels
  
  #neuronchat <- createneuronchat(object = data.input.matrix, meta = meta, group.by = "good_names")
  
  neuronchat <- createNeuronChat(object = data.input.matrix, meta = meta_temp, group.by = as.vector(data$neuronchat_label))
  neuronchat <- run_NeuronChat(neuronchat,M=100)
  net_aggregated_neuronchat <- net_aggregation(neuronchat@net,method = 'weight')
  
  #head(net_aggregated_neuronchat)
  #net_aggregated_neuronchat["galn","2a"]
  
  library(reshape2)
  
  out <- melt(net_aggregated_neuronchat, varnames = c("sender", "receiver"), value.name = "weight")
  
  
  out[1:length(colnames),1] <- as.vector(t(as.matrix(net_aggregated_neuronchat)))
  
  split_function <- function(x) {
    split_string <- strsplit(x, "_")[[1]]
    return(split_string)
  }
  
  out$connection <- rownames(out)
  
  # Apply the function to each row
  out$sender <- sapply(out$connection, function(x) split_function(x)[1])
  out$receiver <- sapply(out$connection, function(x) split_function(x)[2])
  
  
  write.csv(out, "C:/Users/zvjohns/Desktop/projects/clownfish/neuronchat_weights_child_clusters_111723.csv")

  
  devtools::install_github("jokergoo/circlize")  
  
  library(circlize)
  
  #par(mfrow=c(1,2))
  # Visualization, circle plot, for the aggregated network
  netVisual_circle_neuron(net_aggregated_neuronchat,vertex.label.cex = 3, top=0.05)
  # Visualization, chordDiagram, for the aggregated network; also using neuronchat function netVisual_chord_cell_internal(net_aggregated_x, group = group,lab.cex=1)
  netVisual_chord_neuron(neuronchat,method = 'weight',lab.cex = 1, cut_off = 1)

  heatmap_aggregated(neuronchat, method='weight')
  
  ?netVisual_circle_neuron
  ?netVisual_chord_neuron
  weights_df <- read.csv("C:/Users/zvjohns/Desktop/projects/clownfish/neuronchat_weights_child_13d_neuropeptide_steroid_connectivity_111623.csv")

  hist(weights_df$weight, breaks = 200)  
  
  corr_values <- as.vector(net_aggregated_neuronchat)
  
  # Calculate the 95th percentile
  percentile_95 <- quantile(corr_values, 0.95)
  
  # Print the result
  print(paste("The 95th percentile value is:", percentile_95))
  