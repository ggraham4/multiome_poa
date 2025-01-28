# COLTAN G PARKER
# RHODES LAB UIUC


#==== Install CellChat from Source =============================================
# first make sure devtools is installed and loaded
#install.packages("devtools")

# need to install dependencies: ComplexHeatmap, BiocGenerics
#devtools::install_github("jokergoo/ComplexHeatmap")
#BiocManager::install("BiocGenerics")
#BiocManager::install("Biobase")

# next try installing CellChat - will probably need to install Rtools first
#devtools::install_github("sqjin/CellChat")


#==== LOAD LIBRARIES ===========================================================
library(tidyverse)
library(readxl)
library(Seurat)
library(CellChat)
library(patchwork)


#== CellChat Tutorial ==========================================================
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

# load data
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/")
GT_data <- readRDS("GT_ortho_v06.rds")

# save object with clusters named labelparent and labelchild
#saveRDS(GT_data, "GT_ortho_v06.rds")

# try making CellChat object from Seurat object
cellchat <- createCellChat(object = GT_data, 
                           #meta = meta, 
                           group.by = "labelchild")
# verbose output indicates
# - using correct data slot
# - pulling metadata and groups correctly

# call up cellchat database
CellChatDB <- CellChatDB.mouse
dplyr::glimpse(CellChatDB$interaction)

# set this database as the database in the cellchat object
cellchat@DB <- CellChatDB

# need to run this subsetData command then preprocessing continues
cellchat <- subsetData(cellchat) # not sure what this does

# identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat) # this runs once for each cluster
cellchat <- identifyOverExpressedInteractions(cellchat) # this runs quickly

# compute communication probability - slow step, about ten minutes
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# the subsetCommunication() function could be useful at this point
# it is used to access inferred cell-cell communications of interest
# like subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
# which returns the inferred signaling from 1&2 to 4&5

# compute communication at a pathway level - much faster step
cellchat <- computeCommunProbPathway(cellchat)

# compute aggregated network 
cellchat <- aggregateNet(cellchat)


#==== CellChat Object V01 ======================================================
# save cellchat object at this point and test re-loading it
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/")

#saveRDS(cellchat, "GT_cellchat_v01.rds")
cellchat <- readRDS("GT_cellchat_v01.rds")

# first visualization of overall interaction space
groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# looking at the p-values there is no significant communications
cellchat@net$pval


#==== Re-do analysis with a more liberal "trim" value ==========================
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/")
GT_data <- readRDS("GT_ortho_v06.rds")

cellchat <- createCellChat(object = GT_data, 
                           #meta = meta, 
                           group.by = "labelchild")

# call up cellchat database
CellChatDB <- CellChatDB.mouse
# set this database as the database in the cellchat object
cellchat@DB <- CellChatDB
# need to run this subsetData command then preprocessing continues
cellchat <- subsetData(cellchat) # not sure what this does

# identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat) # this runs once for each cluster
cellchat <- identifyOverExpressedInteractions(cellchat) # this runs quickly

# now calculate with trim set to 0
# compute communication probability - slow step 
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "truncatedMean", trim = 0)

# compute communication at a pathway level
cellchat <- computeCommunProbPathway(cellchat)
# compute aggregated network 
cellchat <- aggregateNet(cellchat)

cellchat@net$pval


#==== Re-do analysis conservative trim value and parent clusters ===============

# try making CellChat object from Seurat object
cellchat <- createCellChat(object = GT_data, 
                           #meta = meta, 
                           group.by = "labelparent")
# call up cellchat database
CellChatDB <- CellChatDB.mouse
dplyr::glimpse(CellChatDB$interaction)

# set this database as the database in the cellchat object
cellchat@DB <- CellChatDB

# need to run this subsetData command then preprocessing continues
cellchat <- subsetData(cellchat) # not sure what this does

# identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat) # this runs once for each cluster
cellchat <- identifyOverExpressedInteractions(cellchat) # this runs quickly

# compute communication probability - slow step, about ten minutes
cellchat <- computeCommunProb(cellchat) # note not using raw

# compute communication at a pathway level
cellchat <- computeCommunProbPathway(cellchat)
# compute aggregated network 
cellchat <- aggregateNet(cellchat)


#==== ANALYSIS V01 =============================================================
# generate a cellchat object for males and females, then compare
# GT_ortho_v07 has all clusters labeled properly and unreliable removed

# load data
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/")
GT_data <- readRDS("GT_ortho_v07.rds")

# split to males and females
GT_data_F <- subset(GT_data, sex == "f")
GT_data_M <- subset(GT_data, sex == "m")

# for each sex, generate a CellChat object for parent and child clusters

# call up cellchat database
CellChatDB <- CellChatDB.mouse

# generate cellchat objects
# starting with female (F), then male (M), then all (A)

# decide on parameters for cellchat object generation
# - full ligand-receptor database
# - 1:1 orthologues with mouse database
# - no smoothing -- computeCommunProb(raw.use = TRUE)
# - use trimean default
# - do this for both parent and child clusters for both sexes
# - - remove unreliable cluster


# Female, parent cluster
chat_F_parent <- createCellChat(object = GT_data_F,
                                group.by = "labelparent")
chat_F_parent@DB <- CellChatDB
chat_F_parent <- subsetData(chat_F_parent)
chat_F_parent <- identifyOverExpressedGenes(chat_F_parent) # this runs once for each cluster
chat_F_parent <- identifyOverExpressedInteractions(chat_F_parent) # this runs quickly
chat_F_parent <- computeCommunProb(chat_F_parent)
chat_F_parent <- computeCommunProbPathway(chat_F_parent)
chat_F_parent <- aggregateNet(chat_F_parent)

# female, child clusters
chat_F_child <- createCellChat(object = GT_data_F,
                               group.by = "labelchild")
chat_F_child@DB <- CellChatDB
chat_F_child <- subsetData(chat_F_child)
chat_F_child <- identifyOverExpressedGenes(chat_F_child) # this runs once for each cluster
chat_F_child <- identifyOverExpressedInteractions(chat_F_child) # this runs quickly
chat_F_child <- computeCommunProb(chat_F_child)
chat_F_child <- computeCommunProbPathway(chat_F_child)
chat_F_child <- aggregateNet(chat_F_child)

# MALE, parent cluster
chat_M_parent <- createCellChat(object = GT_data_M,
                                group.by = "labelparent")
chat_M_parent@DB <- CellChatDB
chat_M_parent <- subsetData(chat_M_parent)
chat_M_parent <- identifyOverExpressedGenes(chat_M_parent) # this runs once for each cluster
chat_M_parent <- identifyOverExpressedInteractions(chat_M_parent) # this runs quickly
chat_M_parent <- computeCommunProb(chat_M_parent)
chat_M_parent <- computeCommunProbPathway(chat_M_parent)
chat_M_parent <- aggregateNet(chat_M_parent)

# MALE, child clusters
chat_M_child <- createCellChat(object = GT_data_M,
                               group.by = "labelchild")
chat_M_child@DB <- CellChatDB
chat_M_child <- subsetData(chat_M_child)
chat_M_child <- identifyOverExpressedGenes(chat_M_child) # this runs once for each cluster
chat_M_child <- identifyOverExpressedInteractions(chat_M_child) # this runs quickly
chat_M_child <- computeCommunProb(chat_M_child)
chat_M_child <- computeCommunProbPathway(chat_M_child)
chat_M_child <- aggregateNet(chat_M_child)

# ALL, parent cluster
chat_all_parent <- createCellChat(object = GT_data,
                                  group.by = "labelparent")
chat_all_parent@DB <- CellChatDB
chat_all_parent <- subsetData(chat_all_parent)
chat_all_parent <- identifyOverExpressedGenes(chat_all_parent) # this runs once for each cluster
chat_all_parent <- identifyOverExpressedInteractions(chat_all_parent) # this runs quickly
chat_all_parent <- computeCommunProb(chat_all_parent)
chat_all_parent <- computeCommunProbPathway(chat_all_parent)
chat_all_parent <- aggregateNet(chat_all_parent)

# ALL, child clusters
chat_all_child <- createCellChat(object = GT_data,
                               group.by = "labelchild")
chat_all_child@DB <- CellChatDB
chat_all_child <- subsetData(chat_all_child)
chat_all_child <- identifyOverExpressedGenes(chat_all_child) # this runs once for each cluster
chat_all_child <- identifyOverExpressedInteractions(chat_all_child) # this runs quickly
chat_all_child <- computeCommunProb(chat_all_child)
chat_all_child <- computeCommunProbPathway(chat_all_child)
chat_all_child <- aggregateNet(chat_all_child)


# save RDS objects for all these cellchats for posterity

saveRDS(chat_F_parent, "./results/cellchat/chat_f_parent_2022-03-18.rds")
saveRDS(chat_F_child, "./results/cellchat/chat_f_child_2022-03-18.rds")
saveRDS(chat_M_parent, "./results/cellchat/chat_m_parent_2022-03-18.rds")
saveRDS(chat_M_child, "./results/cellchat/chat_m_child_2022-03-18.rds")
saveRDS(chat_all_parent, "./results/cellchat/chat_all_parent_2022-03-18.rds")
saveRDS(chat_all_child, "./results/cellchat/chat_all_child_2022-03-18.rds")


# for saving output, change wd
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/results/cellchat/")

# steps to merge
chat_list_parent <- list(female = chat_F_parent, male = chat_M_parent)
chat_list_child <- list(female = chat_F_child, male = chat_M_child)

chat_merged_parent <- mergeCellChat(chat_list_parent, add.names = names(chat_list_parent))
chat_merged_child <- mergeCellChat(chat_list_child, add.names = names(chat_list_child))

# initial assessment of differences in connectivity
compareInteractions(chat_merged_parent, show.legend = F, group = c(1,2))
compareInteractions(chat_merged_parent, show.legend = F, group = c(1,2), measure = "weight")

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(chat_merged_parent, weight.scale = T)
netVisual_diffInteraction(chat_merged_parent, weight.scale = T, measure = "weight")

# initial assessment of differences in connectivity
compareInteractions(chat_merged_child, show.legend = F, group = c(1,2))
compareInteractions(chat_merged_child, show.legend = F, group = c(1,2), measure = "weight")

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(chat_merged_child, weight.scale = T)
netVisual_diffInteraction(chat_merged_child, weight.scale = T, measure = "weight")


# this seems off... huge difference between the sexes in interaction weight...
# possibly using the non-raw values will make a difference...


# realizing that if this is pulling log-normalized data that was generated
# on the reduced dataset, not full, then it's totally unrealiable

# need to make sure that I am working from a Seurat object that is:
# - the original fish data with all transformations applied
# - all possible mouse genes renamed

# as far as I can tell, it should work if we just rename the 1:1 genes
# but don't change anything else in the otherwise main seurat object


#==== V03 ====
# load data
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/")
GT_data <- readRDS("ocellaris2020_labeled_v02.rds")

# rename orthologous fish genes to mouse genes
# no need to remove other genes for cellchat purposes

# build Seurat object from scratch
# rename genes
# normalize and name cells based on updated object


#==== building seurat object from scratch ====
# set working directory
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/zack analysis/")

# load filtered feature matrices
m1.data <- Read10X(data.dir = "./filtered matrices/male1/")
m2.data <- Read10X(data.dir = "./filtered matrices/male2/")
f1.data <- Read10X(data.dir = "./filtered matrices/female1/")
f2.data <- Read10X(data.dir = "./filtered matrices/female2/")

setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/MOUSE FISH SYNTHESIS/")
orthologs <- read.csv("orthologous_table.csv")

ortho_fish <- orthologs$Orange.clownfish.gene.name


# CHANGE GENE NAMES TO CORRESPONDING MOUSE GENE NAMES
# f1.data
rownames <- rownames(f1.data)
for (i in rownames) {
  index <- as.numeric(match(i, rownames))
  if (sum(orthologs$Orange.clownfish.gene.name == i) != 0) {
    newname <- orthologs$Gene.name[orthologs$Orange.clownfish.gene.name == i]
    rownames[index] <- newname
  }
}
rownames(f1.data) <- rownames

# f2.data
rownames <- rownames(f2.data)
for (i in rownames) {
  index <- as.numeric(match(i, rownames))
  if (sum(orthologs$Orange.clownfish.gene.name == i) != 0) {
    newname <- orthologs$Gene.name[orthologs$Orange.clownfish.gene.name == i]
    rownames[index] <- newname
  }
}
rownames(f2.data) <- rownames

# m1.data
rownames <- rownames(m1.data)
for (i in rownames) {
  index <- as.numeric(match(i, rownames))
  if (sum(orthologs$Orange.clownfish.gene.name == i) != 0) {
    newname <- orthologs$Gene.name[orthologs$Orange.clownfish.gene.name == i]
    rownames[index] <- newname
  }
}
rownames(m1.data) <- rownames

# m2.data
rownames <- rownames(m2.data)
for (i in rownames) {
  index <- as.numeric(match(i, rownames))
  if (sum(orthologs$Orange.clownfish.gene.name == i) != 0) {
    newname <- orthologs$Gene.name[orthologs$Orange.clownfish.gene.name == i]
    rownames[index] <- newname
  }
}
rownames(m2.data) <- rownames

# add a step to add f_ or m_ before cell names
# f1.data
cellnames <- colnames(f1.data)
cellnames <- paste("f", cellnames, sep="_")
colnames(f1.data) <- cellnames

# f2.data
cellnames <- colnames(f2.data)
cellnames <- paste("f", cellnames, sep="_")
colnames(f2.data) <- cellnames

# m1.data
cellnames <- colnames(m1.data)
cellnames <- paste("m", cellnames, sep="_")
colnames(m1.data) <- cellnames

# m2.data
cellnames <- colnames(m2.data)
cellnames <- paste("m", cellnames, sep="_")
colnames(m2.data) <- cellnames


# convert to Seurat objects, merge
# then figure out how to only keep the cells that are present in the
# object that we have been working with from zack

#create objects
m1 <- CreateSeuratObject(counts = m1.data, project = "m1")
m2 <- CreateSeuratObject(counts = m2.data, project = "m2")
f1 <- CreateSeuratObject(counts = f1.data, project = "f1")
f2 <- CreateSeuratObject(counts = f2.data, project = "f2")

#name samples
m1$sample <- "m1"
m2$sample <- "m2"
f1$sample <- "f1"
f2$sample <- "f2"

#assign sex
m1$sex <- "male"
m2$sex <- "male"
f1$sex <- "female"
f2$sex <- "female"

combined <- merge(x=m1, y=c(m2,f1,f2), merge.data = TRUE)

# now need to keep only the cells in the real seurat object
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/")
GT_data <- readRDS("ocellaris2020_labeled_v02.rds")

GT_cells <- colnames(GT_data)

# then reduce "combined" to just cells included in GT_data
combined_cells <- colnames(combined)
combined_keepers <- combined[,colnames(combined) %in% GT_cells]


# since this = 0 then cell names are identical and in same order
sum(colnames(GT_data) != colnames(combined_keepers))


# normalize data
# add cluster IDs
# then move on to cellchat workflow

combined_keepers <- NormalizeData(combined_keepers, normalization.method = "LogNormalize")
combined_keepers@meta.data <- GT_data@meta.data


# save this object
saveRDS(combined_keepers, "ocellaris2020_seurat_for_cellchat_v01.rds")


# split to males and females
object_f <- subset(combined_keepers, sex == "f")
object_m <- subset(combined_keepers, sex == "m")


# generated cellchat objects as above

saveRDS(chat_F_parent, "./results/cellchat/chat_f_parent_2022-03-20.rds")
saveRDS(chat_F_child, "./results/cellchat/chat_f_child_2022-03-20.rds")
saveRDS(chat_M_parent, "./results/cellchat/chat_m_parent_2022-03-20.rds")
saveRDS(chat_M_child, "./results/cellchat/chat_m_child_2022-03-20.rds")
saveRDS(chat_all_parent, "./results/cellchat/chat_all_parent_2022-03-20.rds")
saveRDS(chat_all_child, "./results/cellchat/chat_all_child_2022-03-20.rds")


# for saving output, change wd
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/results/cellchat/")

# steps to merge
chat_list_parent <- list(female = chat_F_parent, male = chat_M_parent)
chat_list_child <- list(female = chat_F_child, male = chat_M_child)

chat_merged_parent <- mergeCellChat(chat_list_parent, add.names = names(chat_list_parent))
chat_merged_child <- mergeCellChat(chat_list_child, add.names = names(chat_list_child))

# initial assessment of differences in connectivity
compareInteractions(chat_merged_parent, show.legend = F, group = c(1,2))
compareInteractions(chat_merged_parent, show.legend = F, group = c(1,2), measure = "weight")

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(chat_merged_parent, weight.scale = T)
netVisual_diffInteraction(chat_merged_parent, weight.scale = T, measure = "weight")

# initial assessment of differences in connectivity
compareInteractions(chat_merged_child, show.legend = F, group = c(1,2))
compareInteractions(chat_merged_child, show.legend = F, group = c(1,2), measure = "weight")

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(chat_merged_child, weight.scale = T)
netVisual_diffInteraction(chat_merged_child, weight.scale = T, measure = "weight")




#==== CellChat for Thesis ====

# load in seurat object
# keep just a simplified list of clusters, per justin's rec
# start with just neuronal clusters that are sexually differentiated, plus key gene clusters

# sexually diff: male and female diff, high DEGs count
# GnRH cluster, Kiss cluster, Gal cluster

# Clusters to Keep
# 13 - GnRH
# 37 - Kiss
# 25 - Galanin/nonapeptides/steroid receptors
# 35 - male cluster
# 17, 20, 23, 31 - female clusters
# 45, 0, 11, 30, 26, 6, 7, 2 - high DEGs, some also male or female

# that is 16 clusters total, a good number to stop at I think

keep <- c("13", "37", "25", "35", "17", "20", "23", "31", "45", "0", "11", "30",
          "26", "6", "7", "2")


# step 1 - read in seurat object and reduce to just wanted clusters
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/")
gt_data <- 
  readRDS("ocellaris2020_seurat_for_cellchat_v01.rds") %>%
  subset(subset = clusters49 %in% keep)


# quick rename childclusters
clusteranalyses <- read_excel("cluster_analyses.xlsx")
meta <- gt_data@meta.data
list49 <- as.vector(meta$clusters49)
list49label <- c()
list49labelshort <- c()
clusteranalyses49 <- filter(clusteranalyses, clusterlevel == "child")
for (i in list49) {
  labl <- clusteranalyses49$label[clusteranalyses49$clusterID == i]
  list49label <- append(list49label, labl)
  
  labl_short <- clusteranalyses49$label_short[clusteranalyses49$clusterID == i]
  list49labelshort <- append(list49labelshort, labl_short)
}
meta$label_child <- list49label
meta$label_child_short <- list49labelshort
gt_data@meta.data <- meta


# step 2 - split to males and females
gt_data_f <- subset(gt_data, sex == "f")
gt_data_m <- subset(gt_data, sex == "m")

# overall more cells in females compared to males
# this makes sense, given many of the included clusters are overrep. in females


# step 3 - generate cellchat objects for each 
  # parameters
  # - use full ligand-receptor database (not limit to just secreted or just membrane-bound)
  # - genes used are 1:1 orthologues with mouse database
  # - allow smoothing in computeCommunProb() (account for dropouts)
  # - use trimean default


# call up cellchat database
CellChatDB <- CellChatDB.mouse

# create one chat object for each


# all (males and females lumped)
chat_all <- createCellChat(object = gt_data, group.by = "label_child_short")
chat_all@DB <- CellChatDB

chat_all <- 
  chat_all %>%
  subsetData() %>%
  identifyOverExpressedGenes() %>%
  identifyOverExpressedInteractions() %>%
  computeCommunProb() %>%
  computeCommunProbPathway() %>%
  aggregateNet()

saveRDS(chat_all, "./results/cellchat/chat-all-2022-06-03.rds")


# females
chat_f <- createCellChat(object = gt_data_f, group.by = "label_child_short")
chat_f@DB <- CellChatDB

chat_f <- 
  chat_f %>%
  subsetData() %>%
  identifyOverExpressedGenes() %>%
  identifyOverExpressedInteractions() %>%
  computeCommunProb() %>%
  computeCommunProbPathway() %>%
  aggregateNet()

saveRDS(chat_f, "./results/cellchat/chat-f-2022-06-03.rds")


# males
chat_m <- createCellChat(object = gt_data_m, group.by = "label_child_short")
chat_m@DB <- CellChatDB

chat_m <- 
  chat_m %>%
  subsetData() %>%
  identifyOverExpressedGenes() %>%
  identifyOverExpressedInteractions() %>%
  computeCommunProb() %>%
  computeCommunProbPathway() %>%
  aggregateNet()

saveRDS(chat_m, "./results/cellchat/chat_m-2022-06-03.rds")


# step 4 - merge male and female chat objects to do sex diff analysis

chat_list_fm <- list(male = chat_m, female = chat_f)

chat_merged_fm <- mergeCellChat(chat_list_fm, add.names = names(chat_list_fm))

saveRDS(chat_merged_fm, "./results/cellchat/chat_merged_fm-2022-06-03.rds")


# initial assessment of differences in connectivity
compareInteractions(chat_merged_fm, show.legend = F, group = c(1,2))
compareInteractions(chat_merged_fm, show.legend = F, group = c(1,2), measure = "weight")

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(chat_merged_fm, weight.scale = T)
netVisual_diffInteraction(chat_merged_fm, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(chat_merged_fm)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(chat_merged_fm, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2







