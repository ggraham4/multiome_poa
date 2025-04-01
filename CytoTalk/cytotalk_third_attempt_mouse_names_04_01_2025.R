library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(CytoTalk)

library(reticulate)
py_require('pcst_fast')
pcst_fast = import('pcst_fast')

#read in human named data
mouse_named =readRDS("C:/Users/Gabe/Desktop/RNA object mouse names.rds")

cytotalk_wrapper = function(seurat_object,
                            sender,
                            receiver){
  library(parallel)
  library(CytoTalk)
  library(reticulate)
  py_require('pcst_fast')
  pcst_fast = import('pcst_fast')
  
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
  
  return(results)
}

results_19_5 <- cytotalk_wrapper(mouse_named, 19, 5)

results_19_5$pathways
results_27_14 <- cytotalk_wrapper(mouse_named, 27, 14)
results_27_14$pathways

