# COLTAN G PARKER
# RHODES LAB UIUC


#==== Install CellChat from Source =============================================
# first make sure devtools is installed and loaded
#install.packages("devtools")

# need to install dependencies: ComplexHeatmap, BiocGenerics
#devtools::install_github("jokergoo/ComplexHeatmap")
#BiocManager::install("BiocGenerics")

# next try installing CellChat - will probably need to install Rtools first
#devtools::install_github("sqjin/CellChat")


#== LOAD LIBRARIES ==========================================================

#library(devtools)
library(tidyverse)
library(readxl)
library(Seurat)


#== 1. LOAD DATA ===============================================================
# load data and subset neurons
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/")
GT_data <- readRDS("GT_data_wClasses_2022-02-07.rds")
GT_neurons <- subset(GT_data, class_two == "Neurons")
rm(GT_data)





