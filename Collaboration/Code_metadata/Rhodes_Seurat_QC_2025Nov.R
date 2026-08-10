

#### Load packages ####

library(tidyverse)
library(Seurat)
library(dittoSeq)
library(glmGamPoi)
library(patchwork)


# Change to appropriate working directory:

setwd("projects/rhodes/2025Nov-scRNASeq/results")


# NOTE: the codes are from the very beginning.


# BEGIN CODES ----

#### Read in data #### 
aggr_data <- Read10X(data.dir = "cellranger/Aggregated/outs/count/filtered_feature_bc_matrix/")

class(aggr_data)
# [1] "dgCMatrix"
# attr(,"package")
# [1] "Matrix"
# How many rows (genes) and columns (cells)?

dim(aggr_data)
# [1] 28004 25068


#### Create Seurat Object ####

# Seurat has their own S4 object designed
# to hold lots of different information
# needed for single cell experiments.
# We need to put our count matrix into
# a Seurat object. 

head(colnames(aggr_data)) # -1
tail(colnames(aggr_data)) # -4

# To view all of it
#View(as.data.frame(colnames(aggr_data)))

# The cell barcodes end in -1,-2,-3,-4 but which one is which sample?
# It is based on the order of the files in the aggregation.csv
# 10X helpfully puts a copy in the output:

sampOrder <- read.csv("cellranger/Aggregated/outs/aggregation.csv")
sampOrder

#     Genotype      Batch
# 1 JTS08_Female1     1
# 2 JTS08_Female2     1
# 3   JTS08_Male1     2
# 4   JTS08_male2     2


# Let's create our Seurat object. Our cell IDs
# have the barcode separated from the sample by
# a dash "-" and the sample IDs are in the second
# position. 


allcells <- CreateSeuratObject(counts = aggr_data, 
                               names.delim = "-",
                               names.field = 2, 
                               min.cells = 0)

allcells
# An object of class Seurat 
# 28004 features across 25068 samples within 1 assay 
# Active assay: RNA (28004 features, 0 variable features)
# 1 layer present: counts

levels(allcells$orig.ident)
#[1] "1" "2" "3" "4"

levels(allcells$orig.ident) <- sampOrder$sample_id



# How many cells does each sample have?

allcells$orig.ident %>% table()
# JTS08_Female1 JTS08_Female2   JTS08_Male1   JTS08_male2 
#     7751          6175          4441          6701 


Idents(allcells) <- "orig.ident"

# double check you didn't make an error:

Idents(allcells) %>% table()

# Get rid of unnecessary objects to save space: 

rm(aggr_data)
gc()  #to reclaim memory



#### Cell QC ####

# https://satijalab.org/stats/articles/pbmc3k_tutorial#qc-and-selecting-cells-for-further-analysis

# Let's look at the range of number of 
# UMIs (nCount_RNA) and number of
# genes detected (nFeature_RNA)

summary(allcells[["nCount_RNA"]])
#   nCount_RNA   
# Min.   :  500  
# 1st Qu.: 2435  
# Median : 4726  
# Mean   : 5659  
# 3rd Qu.: 7666  
# Max.   :54539  
summary(allcells[["nFeature_RNA"]])
# nFeature_RNA 
# Min.   : 270  
# 1st Qu.:1385  
# Median :2225  
# Mean   :2283  
# 3rd Qu.:3032  
# Max.   :9284  


#### Calculate MT percentages ####

# Need to get genes on MT chromosome; 
# Jenny found that for this clown fish genome, 
# they start with KEG47- and are the rownames:

rownames(allcells) %>% head()

MT.genes <- grep("^KEG47-", rownames(allcells), value = TRUE)
length(MT.genes)
# 37
MT.genes
# [1] "KEG47-t01" "KEG47-r02" "KEG47-t02" "KEG47-r01" "KEG47-t03" "KEG47-p13" "KEG47-t04" "KEG47-t05"
# [9] "KEG47-t06" "KEG47-p12" "KEG47-t07" "KEG47-t08" "KEG47-t09" "KEG47-t10" "KEG47-t11" "KEG47-p11"
# [17] "KEG47-t12" "KEG47-t13" "KEG47-p10" "KEG47-t14" "KEG47-p09" "KEG47-p08" "KEG47-p07" "KEG47-t15"
# [25] "KEG47-p06" "KEG47-t16" "KEG47-p05" "KEG47-p04" "KEG47-t17" "KEG47-t18" "KEG47-t19" "KEG47-p03"
# [33] "KEG47-p02" "KEG47-t20" "KEG47-p01" "KEG47-t21" "KEG47-t22"

#calculate percentage of UMIs in MT.genes

allcells[["percent.mt"]] <- PercentageFeatureSet(allcells, 
                                                 features = MT.genes)
head(allcells)

summary(allcells[["percent.mt"]])
# percent.mt     
# Min.   : 0.0000  
# 1st Qu.: 0.1659  
# Median : 0.3624  
# Mean   : 0.5890  
# 3rd Qu.: 0.7849  
# Max.   :16.0784



# Do violin plots of the values, which are a 
# mixture of dotplots and density plots

# x11(width = 12, height = 6) 
# quartz()
VlnPlot(allcells, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0)
ggsave("stats/QC_cell_prefilter.jpeg")

#x11(width = 12, height = 6) 
# quartz()
VlnPlot(allcells, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0, log = TRUE)
ggsave("stats/QC_cell_prefilter_log.jpeg")


# Read in cell barcodes to keep ----
# Code below assumes the files in individual_preprocess_outputs.zip
# are extracted and in a subdirectory named "preprocess"

# get file names

fn <- dir(path = "preprocess", pattern = "kept_barcodes.txt", full.names = TRUE)
fn # check to see whether we have all samples


# read in one file
keep.bc <- read.delim(fn[1])
head(keep.bc)

for (i in 2:length(fn)){
  temp <- read.delim(fn[i])
  keep.bc <- rbind(keep.bc, temp)
}

dim(keep.bc)
#20789     2
ncol(allcells)
#25068


sum(keep.bc$barcode %in% colnames(allcells))
#20789

# cut keep.bc down to only those in allcells

keep.bc <- keep.bc[keep.bc$barcode %in% colnames(allcells),]
dim(keep.bc)
#20789     2


# change keep.bc$orig.ident to factor with same levels as allcells$orig.ident

keep.bc$orig.ident <- factor(keep.bc$orig.ident, levels = levels(allcells$orig.ident))


# see how many currently, and how many would be kept

table(allcells$orig.ident)
#  JTS08_Female1 JTS08_Female2   JTS08_Male1   JTS08_male2 
#     7751          6175          4441          6701 

table(keep.bc$orig.ident)
# JTS08_Female1 JTS08_Female2   JTS08_Male1   JTS08_male2 
#     6871          4982          3595          5341 

# total % of cells to be kept
nrow(keep.bc) / ncol(allcells) *100
# 82.93043

# per-sample % kept

sort(table(keep.bc$orig.ident) / table(allcells$orig.ident) * 100)
#      JTS08_male2 JTS08_Female2   JTS08_Male1 JTS08_Female1 
#         79.70452      80.68016      80.95024      88.64663 

# a bit of range in the kept number of cells. 

# See where cells filtered out:

fn <- dir(path = "preprocess/", pattern = "ncells_filtered.txt", full.names = TRUE)

ncells <- read.delim(fn[1])

for (i in 2:length(fn)){
  temp <- read.delim(fn[i])
  ncells <- rbind(ncells, temp)
}

dim(ncells)
#32  5
head(ncells)
# These are the number of cells left after each filtering step
# Calculate the number lost in each

ncells2 <- data.frame(ncells[,1:2], 
                      QCfilt = ncells[,2] -  ncells[,3],
                      Doublet = ncells[,3] -  ncells[,4],
                      Kept = ncells[,4])

summary(ncells2$QCfilt / ncells2$ncell.beforeQC * 100)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  5.135  10.235  12.340  10.979  13.084  14.102


summary(ncells2$Doublet / ncells2$ncell.beforeQC * 100)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6.193   6.212   6.262   6.525   6.575   7.385

# This is a little lower than the expected 8%, but some
# doublets likely first filtered by high UMI counts


# plot all:

ncells3 <- tidyr::pivot_longer(ncells2, 
                               QCfilt:Kept,
                               names_to = "Fate", 
                               values_to = "Number")


# Convert variables to factors ----
ncells3$Fate <- factor(ncells3$Fate, levels = unique(ncells3$Fate))
ncells3$sample <- factor(ncells3$sample, levels = levels(allcells$orig.ident))

x11()
ggplot(ncells3, aes(x = sample, y = Number, fill = Fate)) +
  geom_bar(position="fill", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 8)) +
  ggtitle("Percentage of cells in filtering fates") +
  ylab("Percent of cells") +
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual("Fate",values = dittoSeq::dittoColors()[1:18])
ggsave("stats/Pct_Cells_filtered.jpeg")






# Filter out:

allcells <- subset(allcells, 
                   cells = keep.bc$barcode)

# add in separate celltype calls 

rownames(keep.bc) <- keep.bc$barcode


# Do violin plots of the values, which are a 
# mixture of dotplots and density plots

# NOTE for figures - most have two lines commented out before hand
# Not running either will open the plot in RStudio's small Plots tab
# x11() will open a separate plotting window on PCs
# quartz() will open a separate plotting window on Macs

# x11(width = 20, height = 6) 
# quartz()
VlnPlot(allcells, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0)
ggsave("stats/QC_cell_postfilter.jpeg")

#x11(width = 12, height = 6) 
# quartz()
VlnPlot(allcells, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0, log = TRUE)
ggsave("stats/QC_cell_postfilter_log.jpeg")



# Run pipeline
DefaultAssay(allcells) <- "RNA"
allcells <- NormalizeData(allcells) # should I run SCTransform?
allcells <- FindVariableFeatures(allcells)
allcells <- ScaleData(allcells, features = rownames(allcells))
allcells <- RunPCA(allcells)
allcells <- RunUMAP(allcells, dims = 1:30)
allcells <- FindNeighbors(allcells, dims = 1:30)
allcells <- FindClusters(allcells, resolution = 0.1)

table(allcells$seurat_clusters)


#### Visualize UMAP  ####

# x11()
# quartz()
# DimPlot(allcells, label = FALSE, reduction = "umap", group.by = "orig.ident") 
# ggsave("stats/UMAP_bySample_2025-11-13.jpeg")
# 
# 
# # Do the same plot, but split by sample:
# 
# 
# # x11(width = 18, height = 10)
# # quartz(width = 18, height = 10)
# DimPlot(allcells, label = FALSE, 
#         reduction = "umap", group.by = "orig.ident", 
#         split.by = "orig.ident", ncol = 2) + 
#   NoLegend()
# ggsave("stats/UMAP_bySample_split_2025-11-13.jpeg")
# 
# x11()
# DimPlot(allcells, label = FALSE, 
#         reduction = "umap", group.by = "seurat_clusters", 
#         ncol = 2) 

#by samples, no split
#x11()
# quartz()
dittoDimPlot(allcells, var = "orig.ident", reduction.use = "umap", 
             do.label = FALSE, labels.highlight = TRUE)
ggsave(paste0("stats/allcells_filtered_UMAP_samples.jpeg"),width = 12, height = 12)

#split by samples
#x11()
# quartz()
dittoDimPlot(allcells, var = "orig.ident", reduction.use = "umap", 
             do.label = FALSE, labels.highlight = TRUE, split.by = "orig.ident")
ggsave(paste0("stats/allcells_filtered_UMAP_samples_split.jpeg"),width = 12, height = 12)


#by clusters, no split
#x11()
# quartz()
dittoDimPlot(allcells, var = "seurat_clusters", reduction.use = "umap", 
             do.label = FALSE, labels.highlight = TRUE)
ggsave(paste0("stats/allcells_filtered_UMAP_clusters.jpeg"),width = 12, height = 12)

#split by clusters
#x11()
# quartz()
dittoDimPlot(allcells, var = "seurat_clusters", reduction.use = "umap", 
                  do.label = FALSE, labels.highlight = TRUE, split.by = "seurat_clusters")
ggsave(paste0("stats/allcells_filtered_UMAP_clusters_split.jpeg"),width = 12, height = 12)



# can save a copy of the Seurat object at this step before 
# checking whether integration needed or not

saveRDS(allcells, "stats/allcells_filtered_2025-11-13.Rds")
# to read back in:
#allcells <- readRDS("stats/allcells_filtered_2025-11-13.Rds")

