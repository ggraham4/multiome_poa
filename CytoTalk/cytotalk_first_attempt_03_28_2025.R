library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(CytoTalk)

#read in human named data
human_named =readRDS("C:/Users/Gabe/Desktop/RNA_object_human_names.rds")

#convert to SCE
human_named_sce =as.SingleCellExperiment(human_named)

#read into cytotalk
lst_scrna <- CytoTalk::from_single_cell_experiment(human_named_sce)

#run cytotalk, lets use 27 -> 14 as an example
type_a = '27'
type_b = '14'

results <- CytoTalk::run_cytotalk(lst_scrna, type_a, type_b)
# holy FUCK this thing is slow