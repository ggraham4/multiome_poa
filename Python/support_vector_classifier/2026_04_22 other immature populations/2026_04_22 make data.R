# libs
library(Seurat)
library(Matrix)
`%notin%` = Negate(`%in%`)

# read in object and subset to only neurons 
obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
Idents(obj) <- "res0.8_50nn_40PC_45LSI"
DimPlot(obj, reduction = 'harmony_wnn.umap')

obj_no_glia = subset(obj, 
                                (res0.8_50nn_40PC_45LSI %notin% c(
                                  2, # olig
                                  13, # opc
                                  20, # leuko
                                  26, # DG
                                  11, # mg
                                  15, # eg
                                  1 ))&# rgc
                                  Status != 'NRM')

DimPlot(obj_no_glia, reduction = 'harmony_wnn.umap')

obj_no_glia = FindVariableFeatures(obj_no_glia, layer = 'data')

variable_RNA = VariableFeatures(obj_no_glia)
RNA_data_matrix = obj_no_glia@assays$RNA$data[variable_RNA,]

obj_no_glia = FindVariableFeatures(obj_no_glia, layer = 'data', assay = 'ATAC')

variable_ATAC =VariableFeatures(obj_no_glia, assay = 'ATAC')

ATAC_data_matrix = obj_no_glia@assays$ATAC$data[variable_ATAC,]

together_variable_matrix = rbind(RNA_data_matrix, ATAC_data_matrix)

writeMM(together_variable_matrix, file="~/Desktop/SVC/2026_04_22 Variable RNA and ATAC features.mtx")
write.csv(obj_no_glia@meta.data, '~/Desktop/SVC/2026_04_22 Variable RNA and ATAC features metadata.csv')


##### make indivudal data for males, intermediates, and females ####

together_variable_matrix = readMM("~/Desktop/SVC/2026_04_22 Variable RNA and ATAC features.mtx")
meta.data = read.csv("~/Desktop/SVC/2026_04_22 Variable RNA and ATAC features metadata.csv")

meta.data.male.indicies = which(meta.data$Status=='M')
meta.data.female.indicies = which(meta.data$Status=='F')
meta.data.intermediate.indicies = which(meta.data$Status=='D')

meta.data.male = meta.data[meta.data.male.indicies,]
meta.data.female = meta.data[meta.data.female.indicies,]
meta.data.intermediate = meta.data[meta.data.intermediate.indicies,]

write.csv(meta.data.male,'/Users/ggraham/Desktop/SVC/2026_04_22 meta_data_male.csv')
write.csv(meta.data.female,'/Users/ggraham/Desktop/SVC/2026_04_22 meta_data_female.csv')
write.csv(meta.data.intermediate,'/Users/ggraham/Desktop/SVC/2026_04_22 meta_data_intermediate.csv')

male.matrix= together_variable_matrix[,meta.data.male.indicies]
female.matrix= together_variable_matrix[,meta.data.female.indicies]
intermediate.matrix= together_variable_matrix[,meta.data.intermediate.indicies]

writeMM(male.matrix, file="/Users/ggraham/Desktop/SVC/2026_04_22 Variable RNA and ATAC features male.mtx")
writeMM(female.matrix, file="/Users/ggraham/Desktop/SVC/2026_04_22 Variable RNA and ATAC features female.mtx")
writeMM(intermediate.matrix, file="/Users/ggraham/Desktop/SVC/2026_04_22 Variable RNA and ATAC features intermediate.mtx")


