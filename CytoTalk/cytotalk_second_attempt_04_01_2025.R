library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(CytoTalk)

library(reticulate)
py_require('pcst_fast')
pcst_fast = import('pcst_fast')

#read in human named data
human_named =readRDS("C:/Users/Gabe/Desktop/RNA_object_human_names.rds")

#convert to SCE
human_named_sce =as.SingleCellExperiment(human_named)

#read into cytotalk
lst_scrna <- CytoTalk::from_single_cell_experiment(human_named_sce)

cells_27 = rownames(human_named@meta.data[human_named@meta.data$harmony.wnn_res0.4_clusters==27,])
cells_14 = rownames(human_named@meta.data[human_named@meta.data$harmony.wnn_res0.4_clusters==14,])


#run cytotalk, lets use 27 -> 14 as an example\
matrix_27 = as.matrix(human_named@assays$RNA$data[,human_named@meta.data$harmony.wnn_res0.4_clusters==27])
matrix_14 = as.matrix(human_named@assays$RNA$data[,human_named@meta.data$harmony.wnn_res0.4_clusters==14])

rownames(matrix_27) <- rownames(human_named@assays$RNA$data[,human_named@meta.data$harmony.wnn_res0.4_clusters==27])
rownames(matrix_14) <- rownames(human_named@assays$RNA$data[,human_named@meta.data$harmony.wnn_res0.4_clusters==14])

#write.csv(matrix_27, "A:/CytoTalk/scRNAseq-data/scRNAseq_27_sparse.csv")
#write.csv(matrix_14, "A:/CytoTalk/scRNAseq-data/scRNAseq_14_sparse.csv")

dir_in <- 'A:/CytoTalk/scRNAseq-data'
lst_scrna <- CytoTalk::read_matrix_folder(dir_in)
table(lst_scrna$cell_types)

type_a = '27_sparse'
type_b = '14_sparse'

results <- CytoTalk::run_cytotalk(lst_scrna, type_a, type_b,
                                  dir_out = 'A:/CytoTalk/27_14',
                                  cores = detectCores()-1)
# holy FUCK this thing is slow
#doesnt even work 

cytotalk_wrapper = function(seurat_object,
                            sender,
                            receiver){
  library(parallel)
  
  suppressWarnings(dir.create(paste0('A:/CytoTalk/',sender,'_',receiver)))

  sender_matrix = as.matrix(seurat_object@assays$RNA$data[,seurat_object@meta.data$harmony.wnn_res0.4_clusters==sender])
  rownames(sender_matrix) <- rownames(seurat_object@assays$RNA$data[,seurat_object@meta.data$harmony.wnn_res0.4_clusters==sender])
  write.csv(sender_matrix, paste0("A:/CytoTalk/scRNAseq-data/scRNAseq_",sender,".csv"))

  receiver_matrix <- as.matrix(seurat_object@assays$RNA$data[,seurat_object@meta.data$harmony.wnn_res0.4_clusters==receiver])
  rownames(receiver_matrix) <- rownames(seurat_object@assays$RNA$data[,seurat_object@meta.data$harmony.wnn_res0.4_clusters==receiver])
  write.csv(receiver_matrix, paste0("A:/CytoTalk/scRNAseq-data/scRNAseq_",receiver,".csv"))
  
  dir_in <- 'A:/CytoTalk/scRNAseq-data'
  lst_scrna <- CytoTalk::read_matrix_folder(dir_in)

  type_a = paste0(sender)
  type_b = paste0(receiver)
  
  results <- CytoTalk::run_cytotalk(lst_scrna, type_a, type_b,
                                    dir_out =paste0('A:/CytoTalk/',sender,'_',receiver),
                                    cores = detectCores()-1)
  
  
  print('Completed')
  
}

cytotalk_wrapper(human_named, 19, 5)

