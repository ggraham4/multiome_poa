
#https://bioconductor.org/packages/release/bioc/vignettes/scds/inst/doc/scds.html

# Biocluster: module load R/4.4.0-IGB-gcc-8.2.0

print("starting R script; set working directory")

setwd("/home/groups/hpcbio/RNA-Seq/projects/rhodes/2025Nov-scRNASeq/results/")
#setwd("projects/rhodes/2025Nov-scRNASeq/results/")


# Create output directory
dir.create("preprocess")


# Load packages
suppressPackageStartupMessages({
library(scds)
library(scater)
library(rsvd)
library(Rtsne)
library(cowplot)
library(dittoSeq)
library(Seurat)
library(SeuratObject)
})

options(future.globals.maxSize = 8.0e10)


# Pull in arguments
args <- as.numeric(commandArgs(TRUE))

# Read in sample IDs and order

sampOrder <- read.csv("cellranger/Aggregated/outs/aggregation.csv")

# get sample ID based on input number
i <- sampOrder$sample_id[args[1]]

aggr_data <- Read10X(data.dir = paste0("cellranger/",i,"/outs/filtered_feature_bc_matrix/"))
# put same ending on cell barcodes
colnames(aggr_data) <- sub("1", args[1], colnames(aggr_data))

# create Seurat Object; keep only genes in 10 cells

allcells <- CreateSeuratObject(counts = aggr_data, 
                               names.delim = "-",
                               names.field = 2,
                               min.cells = 10)
allcells$orig.ident <- i

rm(aggr_data)
gc()  #to reclaim memory


#### Cell QC ####

#### Calculate MT percentages ####

MT.genes <- grep("^KEG47-", rownames(allcells), value = TRUE)
allcells[["percent.mt"]] <- PercentageFeatureSet(allcells, 
                                                 features = MT.genes)

# Find 3 MADs

temp1 <- scater::isOutlier(allcells$percent.mt, 
                           nmads = 3, type = "higher")
mt.threshold <- min(allcells$percent.mt[temp1])
mt.threshold
# 1.975851 - MAC

temp1 <- scater::isOutlier(allcells$nCount_RNA, 
                           nmads = 3, type = "higher")
UMI.threshold <- min(allcells$nCount_RNA[temp1])
UMI.threshold
# 16558 - MAC


# Make QC graphs

p.mt<- VlnPlot(allcells, features = c("percent.mt"),
               pt.size = 1, log = FALSE) + 
  geom_hline(yintercept = mt.threshold, 
             color = "red", size = 1.5) + NoLegend()

p.umi <- VlnPlot(allcells, features = c("nCount_RNA"),
                 pt.size = 1, log = FALSE) + 
  geom_hline(yintercept = UMI.threshold, 
             color = "red", size = 1.5) + NoLegend()

p.gene <- VlnPlot(allcells, features = c("nFeature_RNA"),
                  pt.size = 1, log = FALSE) + NoLegend()

ncell.beforeQC <- ncol(allcells)

allcells <- subset(allcells, 
                   subset = percent.mt < mt.threshold & 
                     nCount_RNA < UMI.threshold)

ncell.qcfilt <- ncol(allcells)

# Do doublet computation

sce <- as.SingleCellExperiment(allcells)
sce  <- cxds(sce,retRes = TRUE)
sce  <-  bcds(sce,retRes = TRUE)
sce  <-  cxds_bcds_hybrid(sce)

allcells$hybrid_score <- sce$hybrid_score


summary(allcells$hybrid_score)
# Mac
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000165 0.101497 0.259368 0.381782 0.542400 1.970795 

# get 3 MAD

temp1 <- scater::isOutlier(allcells$hybrid_score, 
                           nmads = 3, type = "higher")
db.threshold <- min(allcells$hybrid_score[temp1])
db.threshold
# 1.098414 - Mac


p.db<- VlnPlot(allcells, features = c("hybrid_score"),
               pt.size = 1, log = FALSE) + 
  geom_hline(yintercept = db.threshold, 
             color = "red", size = 1.5) + NoLegend()

combined_plot <- plot_grid(p.mt, p.umi, p.gene, p.db, ncol = 2, nrow = 2)

ggsave(paste0("preprocess/",i, "_QCplots_before.jpeg"), plot = combined_plot, 
       width = 12, height = 10, units = "in", dpi = 300)

allcells <- subset(allcells, 
                   subset = hybrid_score < db.threshold)


ncell.dbfilt <- ncol(allcells)

p.post <- VlnPlot(allcells, features = c("percent.mt", "nCount_RNA", "nFeature_RNA", "hybrid_score"),
                  pt.size = 1, log = FALSE, ncol = 2)
ggsave(paste0("preprocess/",i, "_QCplots_filt.jpeg"), plot = p.post, 
       width = 12, height = 10, units = "in", dpi = 300)


outNum <- data.frame(sample = i,
                     ncell.beforeQC, ncell.qcfilt, ncell.dbfilt,
                     mt.threshold, UMI.threshold, db.threshold)

write.table(outNum, file = paste0("preprocess/", i, "_ncells_filtered.txt"), row.names = FALSE, sep = "\t")


temp <- data.frame(barcode = colnames(allcells),
                   orig.ident = allcells$orig.ident)

# write out kept barcodes and number of cells filtered at each step

write.table(temp, file = paste0("preprocess/", i,"_kept_barcodes.txt"), row.names = FALSE, sep = "\t")


print("done")

