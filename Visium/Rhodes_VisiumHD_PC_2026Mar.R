

#set directory 

setwd("projects/rhodes/2026Mar-spatialRNASeq/")



#### Load packages ####

library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(dittoSeq) 
library(cowplot)
# library(SummarizedExperiment)
# library(muscat)
# library(SpatialExperiment)
# library(reshape2)


# 2 HD slides: H1-ZMF76DB (Run 1) and H1-9YFKHB7 (Run 2) has 4P10 and 6P17
# 4 capture areas total : 3P5 and 4P5 (Run 1), 4P10 and 6P17 (Run 2)
# Each capture area has 4 brains, barcodes for each gotten in Loupe Browser by Justin
# so 16 total samples; extra annotations for each sample given in the names

# Seurat vignette for HD with cell segmentation: https://satijalab.org/seurat/articles/visiumhd_analysis_cell_segmentations
# Seurat longer vignette for HD analysis: https://satijalab.org/seurat/articles/visiumhd_analysis_vignette


# Read in cell-sample info from Justin ----

brains_3P5 <- read.csv("results/spaceranger/Cells_IDs_per_sample/3P5_A1_brains.csv")
brains_4P5 <- read.csv("results/spaceranger/Cells_IDs_per_sample/4P5_D1_brains.csv")
brains_4P10 <- read.csv("results/spaceranger/Cells_IDs_per_sample/4P10_A1_brains.csv")
brains_6P17 <- read.csv("results/spaceranger/Cells_IDs_per_sample/6P17_D1_brains.csv")


table(brains_3P5$brains)
# 3P5_A1_cell_F_1_AtlasNA   3P5_A1_cell_F_2_AtlasNA 3P5_A1_cell_M_1_Atlas9_22 
#                    9999                     50831                     32443 
# 3P5_A1_cell_M_2_Atlas9 
#                  23707

table(brains_4P5$brains)
# 4P5_D1_cell_F_1_Atlas7 4P5_D1_cell_F_2_Atlas22  4P5_D1_cell_M_1_Atlas5 4P5_D1_cell_M_2_Atlas15 
#                  18171                   47640                   16482                   22411 

table(brains_4P10$brains)
# 4P10_A1_cell_F_1_Atlas13 4P10_A1_cell_F_2_Atlas24  4P10_A1_cell_M_1_Atlas6 4P10_A1_cell_M_2_Atlas20 
#                    17942                    56527                    23397                    3997

table(brains_6P17$brains)
# 6P17_D1_cell_F_1_Atlas10  6P17_D1_cell_M_1_Atlas8  6P17_D1_cell_M_2_Atlas9  6P17D1_cell_F_2_Atlas12 
#                    30012                    26140                    22669                    2477

# one of 6P17's is misformatted as 6P17D1

brains_6P17$brains <- sub("6P17D1", "6P17_D1", brains_6P17$brains)



# All have slice ID and capture area. F1, F2, M1 and M2 repeated per slice.
# Altas shows part of brain captured

# Compare with new files sent that have each brain region annotated:

brains_3P5b <- read.csv("results/spaceranger/Cells_IDs_per_sample/3P5_A1_brain_regions.csv")
brains_4P5b <- read.csv("results/spaceranger/Cells_IDs_per_sample/4P5_D1_brain_regions.csv")
brains_4P10b <- read.csv("results/spaceranger/Cells_IDs_per_sample/4P10_A1_brain_regions.csv")
brains_6P17b <- read.csv("results/spaceranger/Cells_IDs_per_sample/6P17_D1_brain_regions.csv")
names(brains_6P17b)[2] <- "brain.regions"

table(brains_3P5b$brain.regions)
# different Atlas areas for the same sample; no longer have AtlasNA so use these.

table(brains_4P5b$brain.regions)
table(brains_4P10b$brain.regions)
table(brains_6P17b$brain.regions)

sapply(strsplit(brains_3P5b$brain.regions, "_"), length) %>% table()
sapply(strsplit(brains_4P5b$brain.regions, "_"), length) %>% table()
sapply(strsplit(brains_4P10b$brain.regions, "_"), length) %>% table()
sapply(strsplit(brains_6P17b$brain.regions, "_"), length) %>% table()

# all have 7

nrow(brains_4P5)
#104704
nrow(brains_4P5b)
# 97971
sum(brains_4P5$Barcode %in% brains_4P5b$Barcode)
# 97918


# not the same number! Will need to use both

# separate out Important parts; one has an extra separator that looks to be a part of Atlas

brains_3P5 <- tidyr::separate_wider_delim(brains_3P5, cols = brains,
                                          delim = "_",
                                          names = c("Slice", "CA","Type","Sex","Individual","Atlas"),
                                          too_many = "merge",
                                          cols_remove = FALSE)
head(brains_3P5)


brains_3P5b <- tidyr::separate_wider_delim(brains_3P5b, cols = brain.regions,
                                          delim = "_",
                                          names = c("Slice", "CA","Type","Sex","Individual","Atlas","Region"),
                                          too_many = "merge",
                                          cols_remove = FALSE)
head(brains_3P5b)




# combine Sex, Individual, and Slice to make unique individual ID

brains_3P5 <- tidyr::unite(brains_3P5, col = "Sample", Sex, Individual, Slice, remove = FALSE)
head(brains_3P5)
brains_3P5b <- tidyr::unite(brains_3P5b, col = "Sample", Sex, Individual, Slice, remove = FALSE)
head(brains_3P5b)

# Start with all cells in the brain region one

use_3P5 <- brains_3P5b
names(use_3P5)[10] <- "fromCSV"

# modify the extra cells in the first one

temp <-  brains_3P5[! brains_3P5$Barcode %in% brains_3P5b$Barcode,]
temp$Region <- "none"
names(temp)[9] <- "fromCSV"
temp <- temp[,names(use_3P5)]
head(temp)

use_3P5 <- rbind(use_3P5, temp)
dim(use_3P5)
# 116996

table(use_3P5$Sample)
# F_1_3P5 F_2_3P5 M_1_3P5 M_2_3P5 
#    9999   50831   32459   23707


# do the same for the rest; 4P5 ----

# separate out Important parts; one has an extra separator that looks to be a part of Atlas

brains_4P5 <- tidyr::separate_wider_delim(brains_4P5, cols = brains,
                                          delim = "_",
                                          names = c("Slice", "CA","Type","Sex","Individual","Atlas"),
                                          too_many = "merge",
                                          cols_remove = FALSE)
brains_4P5b <- tidyr::separate_wider_delim(brains_4P5b, cols = brain.regions,
                                           delim = "_",
                                           names = c("Slice", "CA","Type","Sex","Individual","Atlas","Region"),
                                           too_many = "merge",
                                           cols_remove = FALSE)

# combine Sex, Individual, and Slice to make unique individual ID

brains_4P5 <- tidyr::unite(brains_4P5, col = "Sample", Sex, Individual, Slice, remove = FALSE)
brains_4P5b <- tidyr::unite(brains_4P5b, col = "Sample", Sex, Individual, Slice, remove = FALSE)

# Start with all cells in the brain region one

use_4P5 <- brains_4P5b
names(use_4P5)[10] <- "fromCSV"

# modify the extra cells in the first one

temp <-  brains_4P5[! brains_4P5$Barcode %in% brains_4P5b$Barcode,]
temp$Region <- "none"
names(temp)[9] <- "fromCSV"
temp <- temp[,names(use_4P5)]

use_4P5 <- rbind(use_4P5, temp)
dim(use_4P5)
# 104757

table(use_4P5$Sample)
# F_1_4P5 F_2_4P5 M_1_4P5 M_2_4P5 
#   18214   47650   16482   22411


# do the same for the rest; 4P10 ----

# separate out Important parts; one has an extra separator that looks to be a part of Atlas

brains_4P10 <- tidyr::separate_wider_delim(brains_4P10, cols = brains,
                                          delim = "_",
                                          names = c("Slice", "CA","Type","Sex","Individual","Atlas"),
                                          too_many = "merge",
                                          cols_remove = FALSE)
brains_4P10b <- tidyr::separate_wider_delim(brains_4P10b, cols = brain.regions,
                                           delim = "_",
                                           names = c("Slice", "CA","Type","Sex","Individual","Atlas","Region"),
                                           too_many = "merge",
                                           cols_remove = FALSE)

# combine Sex, Individual, and Slice to make unique individual ID

brains_4P10 <- tidyr::unite(brains_4P10, col = "Sample", Sex, Individual, Slice, remove = FALSE)
brains_4P10b <- tidyr::unite(brains_4P10b, col = "Sample", Sex, Individual, Slice, remove = FALSE)

# Start with all cells in the brain region one

use_4P10 <- brains_4P10b
names(use_4P10)[10] <- "fromCSV"

# modify the extra cells in the first one

temp <-  brains_4P10[! brains_4P10$Barcode %in% brains_4P10b$Barcode,]
temp$Region <- "none"
names(temp)[9] <- "fromCSV"
temp <- temp[,names(use_4P10)]

use_4P10 <- rbind(use_4P10, temp)
dim(use_4P10)
# 138116 

table(use_4P10$Sample)
# F_1_4P10 F_2_4P10 M_1_4P10 M_2_4P10 
#    17925    56696    23399    40096 


# do the same for the rest; 6P17 ----

# separate out Important parts; one has an extra separator that looks to be a part of Atlas

brains_6P17 <- tidyr::separate_wider_delim(brains_6P17, cols = brains,
                                          delim = "_",
                                          names = c("Slice", "CA","Type","Sex","Individual","Atlas"),
                                          too_many = "merge",
                                          cols_remove = FALSE)
brains_6P17b <- tidyr::separate_wider_delim(brains_6P17b, cols = brain.regions,
                                           delim = "_",
                                           names = c("Slice", "CA","Type","Sex","Individual","Atlas","Region"),
                                           too_many = "merge",
                                           cols_remove = FALSE)

# combine Sex, Individual, and Slice to make unique individual ID

brains_6P17 <- tidyr::unite(brains_6P17, col = "Sample", Sex, Individual, Slice, remove = FALSE)
brains_6P17b <- tidyr::unite(brains_6P17b, col = "Sample", Sex, Individual, Slice, remove = FALSE)

# Start with all cells in the brain region one

use_6P17 <- brains_6P17b
names(use_6P17)[10] <- "fromCSV"

# modify the extra cells in the first one

temp <-  brains_6P17[! brains_6P17$Barcode %in% brains_6P17b$Barcode,]
temp$Region <- "none"
names(temp)[9] <- "fromCSV"
temp <- temp[,names(use_6P17)]

use_6P17 <- rbind(use_6P17, temp)
dim(use_6P17)
# 103597

table(use_6P17$Sample)
# F_1_6P17 F_2_6P17 M_1_6P17 M_2_6P17 
#    30012    24776    26140    22669




# Read in 3P5----

s_3P5 <- Load10X_Spatial("results/spaceranger/3P5_A1/outs/", 
                          slice = "s_3P5", 
                         bin.size = c("polygons"),
                         image.name = "tissue_hires_image.png")
s_3P5
# An object of class Seurat 
# 28004 features across 117216 samples within 1 assay 
# Active assay: Spatial.Polygons (28004 features, 0 variable features)
# 1 layer present: counts
# 1 spatial field of view present: s_3P5.polygons

head(s_3P5)

# Compare to use.3P5

nrow(use_3P5)
#116996

sum(colnames(s_3P5) %in% use_3P5$Barcode)
# 116996

temp <- s_3P5@meta.data
temp <- tibble::rownames_to_column(temp, var = "Barcode")

temp <- left_join(temp, use_3P5)

head(temp)

temp$Sample <- tidyr::replace_na(temp$Sample, "NotBrain")
temp <- tibble::column_to_rownames(temp, "Barcode")

s_3P5@meta.data <- temp

s_3P5$orig.ident <- s_3P5$Sample

# check plot of NotBrain

x11()
SpatialDimPlot(s_3P5,
               group.by = "Sample",
               image.scale = "hires") +
  ggtitle("s_3P5") 
ggsave("results/stats/image_3P5_Samples.jpeg")

# Remove the NotBrain

s_3P5 <- subset(s_3P5, subset = Sample != "NotBrain")
s_3P5


# Read in 4P5----

s_4P5 <- Load10X_Spatial("results/spaceranger/4P5_D1/outs/", 
                         slice = "s_4P5", 
                         bin.size = c("polygons"),
                         image.name = "tissue_hires_image.png")
s_4P5
# An object of class Seurat 
# 28004 features across 104899 samples within 1 assay 
# Active assay: Spatial.Polygons (28004 features, 0 variable features)
# 1 layer present: counts
# 1 spatial field of view present: s_4P5.polygons

head(s_4P5)

# Compare to use.4P5

nrow(use_4P5)
#104757

sum(colnames(s_4P5) %in% use_4P5$Barcode)
# 104757

temp <- s_4P5@meta.data
temp <- tibble::rownames_to_column(temp, var = "Barcode")

temp <- left_join(temp, use_4P5)

head(temp)

temp$Sample <- tidyr::replace_na(temp$Sample, "NotBrain")
temp <- tibble::column_to_rownames(temp, "Barcode")

s_4P5@meta.data <- temp

s_4P5$orig.ident <- s_4P5$Sample

# check plot of NotBrain

x11()
SpatialDimPlot(s_4P5,
               group.by = "Sample",
               image.scale = "hires") +
  ggtitle("s_4P5") 
ggsave("results/stats/image_4P5_Samples.jpeg")

# Remove the NotBrain

s_4P5 <- subset(s_4P5, subset = Sample != "NotBrain")
s_4P5 



# Read in 4P10----

s_4P10 <- Load10X_Spatial("results/spaceranger/4P10_A1/outs/", 
                         slice = "s_4P10", 
                         bin.size = c("polygons"),
                         image.name = "tissue_hires_image.png")
#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
s_4P10
# An object of class Seurat 
# 28004 features across 138544 samples within 1 assay 
# Active assay: Spatial.Polygons (28004 features, 0 variable features)
# 1 layer present: counts
# 1 spatial field of view present: s_4P10.polygons

head(s_4P10)

# Compare to use.4P10

nrow(use_4P10)
#138116

sum(colnames(s_4P10) %in% use_4P10$Barcode)
# 138116

temp <- s_4P10@meta.data
temp <- tibble::rownames_to_column(temp, var = "Barcode")

temp <- left_join(temp, use_4P10)

head(temp)

temp$Sample <- tidyr::replace_na(temp$Sample, "NotBrain")
temp <- tibble::column_to_rownames(temp, "Barcode")

s_4P10@meta.data <- temp

s_4P10$orig.ident <- s_4P10$Sample

# check plot of NotBrain

x11()
SpatialDimPlot(s_4P10,
               group.by = "Sample",
               image.scale = "hires") +
  ggtitle("s_4P10") 
ggsave("results/stats/image_4P10_Samples.jpeg")

# Remove the NotBrain

s_4P10 <- subset(s_4P10, subset = Sample != "NotBrain")
s_4P10


# Read in 6P17----

s_6P17 <- Load10X_Spatial("results/spaceranger/6P17_D1/outs/", 
                         slice = "s_6P17", 
                         bin.size = c("polygons"),
                         image.name = "tissue_hires_image.png")
s_6P17
# An object of class Seurat 
# 28004 features across 133521 samples within 1 assay 
# Active assay: Spatial.Polygons (28004 features, 0 variable features)
# 1 layer present: counts
# 1 spatial field of view present: s_6P17.polygons

head(s_6P17)

# Compare to use.6P17

nrow(use_6P17)
#103597

sum(colnames(s_6P17) %in% use_6P17$Barcode)
# 103597

temp <- s_6P17@meta.data
temp <- tibble::rownames_to_column(temp, var = "Barcode")

temp <- left_join(temp, use_6P17)

head(temp)

temp$Sample <- tidyr::replace_na(temp$Sample, "NotBrain")
temp <- tibble::column_to_rownames(temp, "Barcode")

s_6P17@meta.data <- temp

s_6P17$orig.ident <- s_6P17$Sample

# check plot of NotBrain

x11()
SpatialDimPlot(s_6P17,
               group.by = "Sample",
               image.scale = "hires") +
  ggtitle("s_6P17") 
ggsave("results/stats/image_6P17_Samples.jpeg")

# Remove the NotBrain

s_6P17 <- subset(s_6P17, subset = Sample != "NotBrain")
s_6P17


# Check Atlas and Regions

x11()
SpatialDimPlot(s_3P5,
               group.by = "Atlas",
               image.scale = "lowres",
               cols = dittoColors()) +
  ggtitle("s_3P5") 
ggsave("results/stats/image_3P5_Atlas.jpeg")

x11()
SpatialDimPlot(s_3P5,
               group.by = "Region",
               image.scale = "hires",
               cols = dittoColors()) +
  ggtitle("s_3P5") 
ggsave("results/stats/image_3P5_Region.jpeg")



x11()
SpatialDimPlot(s_4P5,
               group.by = "Atlas",
               image.scale = "lowres",
               cols = dittoColors()) +
  ggtitle("s_4P5") 
ggsave("results/stats/image_4P5_Atlas.jpeg")

x11()
SpatialDimPlot(s_4P5,
               group.by = "Region",
               image.scale = "hires",
               cols = dittoColors()) +
  ggtitle("s_4P5") 
ggsave("results/stats/image_4P5_Region.jpeg")


x11()
SpatialDimPlot(s_4P10,
               group.by = "Atlas",
               image.scale = "lowres",
               cols = dittoColors()) +
  ggtitle("s_4P10") 
ggsave("results/stats/image_4P10_Atlas.jpeg")

x11()
SpatialDimPlot(s_4P10,
               group.by = "Region",
               image.scale = "hires",
               cols = dittoColors()) +
  ggtitle("s_4P10") 
ggsave("results/stats/image_4P10_Region.jpeg")


x11()
SpatialDimPlot(s_6P17,
               group.by = "Atlas",
               image.scale = "lowres",
               cols = dittoColors()) +
  ggtitle("s_6P17") 
ggsave("results/stats/image_6P17_Atlas.jpeg")

x11()
SpatialDimPlot(s_6P17,
               group.by = "Region",
               image.scale = "hires",
               cols = dittoColors()) +
  ggtitle("s_6P17") 
ggsave("results/stats/image_6P17_Region.jpeg")


table(paste(s_3P5$Sample, s_3P5$Atlas))





# Merge all samples together ----
# https://satijalab.org/seurat/articles/essential_commands#subsetting-and-merging

allcells <- merge(s_3P5,  y = c(s_4P5, s_4P10, s_6P17),
                  project = "Rhodes",
                  merge.data = FALSE)
# Warning: Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.

allcells
# An object of class Seurat 
# 28004 features across 463466 samples within 1 assay 
# Active assay: Spatial.Polygons (28004 features, 0 variable features)
# 4 layers present: counts.1, counts.2, counts.3, counts.4
# 4 spatial fields of view present: s_3P5.polygons s_4P5.polygons s_4P10.polygons s_6P17.polygons

# Temporarily join layers

allcells[["Spatial.Polygons"]] <- JoinLayers(allcells[["Spatial.Polygons"]])
# to undo:
# allcells[["Spatial.Polygons"]] <- split(allcells[["Spatial.Polygons"]], f = allcells$orig.ident) 


head(allcells[[]])
tail(allcells[[]])


#factor the sample idents

unique(allcells$Sample)
allcells$Sample <- factor(allcells$Sample, levels = paste0(c("F_1_", "F_2_","M_1_", "M_2_"), 
                                                             rep(c("3P5","4P5","4P10","6P17"), each = 4)))
allcells$orig.ident <- allcells$Sample
Idents(allcells) <- "Sample"


unique(allcells$Region) %>% sort()
# [1] "apoa"  "atn"   "ATn"   "cer"   "CM"    "dc"    "Dc"    "dd"    "dl"    "Dl"    "dlt"   "dm"   
# [13] "Dm"    "Dp"    "Gn"    "LFB"   "mpoa"  "NDILI" "NDILm" "nH"    "NLT"   "nMLF"  "none"  "NRL"  
# [25] "NRP"   "ON"    "OpT"   "OT"    "PGa"   "PGc"   "ppoa"  "PSm"   "Tla"   "TLa"   "vd"    "vl"   
# [37] "VMn"   "vs"    "vv"   

# some of these are the same except for capitalization. Asking Justin to clarify...




# Remove individual seurat objects to save space
ls()
#  [1] "allcells"     "brains_3P5"   "brains_3P5b"  "brains_4P10"  "brains_4P10b" "brains_4P5"  
# [7] "brains_4P5b"  "brains_6P17"  "brains_6P17b" "s_3P5"        "s_4P10"       "s_4P5"       
# [13] "s_6P17"       "temp"         "use_3P5"      "use_4P10"     "use_4P5"      "use_6P17"    
rm(list = ls()[2:18])
gc()


table(allcells$orig.ident)
#  F_1_3P5  F_2_3P5  M_1_3P5  M_2_3P5  F_1_4P5  F_2_4P5  M_1_4P5  M_2_4P5 F_1_4P10 F_2_4P10 M_1_4P10 
#     9999    50831    32459    23707    18214    47650    16482    22411    17925    56696    23399 
# M_2_4P10 F_1_6P17 F_2_6P17 M_1_6P17 M_2_6P17 
#    40096    30012    24776    26140    22669

#Any additional empty spots need to be filtered out? ----

# nCount_Spatial is the total number of molecules detected within a cell. 

x11()
VlnPlot(allcells, features = "nCount_Spatial.Polygons", pt.size = 0,
        cols = dittoColors())
#
# All samples are slightly bi-modal

x11()
VlnPlot(allcells, features = "nCount_Spatial.Polygons", log = TRUE, pt.size = 0,
        cols = dittoColors()) +
  geom_hline(yintercept = 25)

ggsave("results/stats/VlnPlot_totalUMIs_2026-05-11.jpeg")

#some very low for all


summary(allcells$nCount_Spatial.Polygons)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     1      61     113     153     196    4362


table(allcells$nCount_Spatial.Polygons < 25, allcells$orig.ident)
#        F_1_3P5 F_2_3P5 M_1_3P5 M_2_3P5 F_1_4P5 F_2_4P5 M_1_4P5 M_2_4P5 F_1_4P10 F_2_4P10 M_1_4P10
# FALSE    9884   49368   28489   22657   16977   45633   16176   21279    16781    52557    22973
# TRUE      115    1463    3970    1050    1237    2017     306    1132     1144     4139      426
# 
#       M_2_4P10 F_1_6P17 F_2_6P17 M_1_6P17 M_2_6P17
# FALSE    39485    27319    21445    23497    20470
# TRUE       611     2693     3331     2643     2199


# where are the spots?

allcells$below25 <- allcells$nCount_Spatial.Polygons < 25

x11(width = 14)
SpatialFeaturePlot(allcells, features = "nCount_Spatial.Polygons", ncol = 4, pt.size.factor = 1) 
# On the lower end of counts for most tissue

x11(width = 14)
SpatialDimPlot(allcells, group.by = "below25", ncol = 4, pt.size.factor = 1) +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave("results/stats/SpatialPlot_under25UMIs_2026-05-11.jpeg")
# Looks like its removing what I'd expect to be blank space



# nFeature_RNA is the number of genes detected in each cell. 

x11()
VlnPlot(allcells, features = "nFeature_Spatial.Polygons", pt.size = 0,
        cols = dittoColors())
#similar to nCount

x11()
VlnPlot(allcells, features = "nFeature_Spatial.Polygons", log = TRUE, pt.size = 0,
        cols = dittoColors()) +
  geom_hline(yintercept = 25)
# also some with very low genes detected

summary(allcells$nFeature_Spatial.Polygons)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1.0    58.0   104.0   136.2   177.0  2074.0

table(allcells$nFeature_Spatial.Polygons < 25, allcells$orig.ident)
#         F_1_3P5 F_2_3P5 M_1_3P5 M_2_3P5 F_1_4P5 F_2_4P5 M_1_4P5 M_2_4P5 F_1_4P10 F_2_4P10 M_1_4P10
# FALSE    9878   49245   28094   22539   16872   45494   16146   21142    16674    52245    22943
# TRUE      121    1586    4365    1168    1342    2156     336    1269     1251     4451      456
# 
# M_2_4P10 F_1_6P17 F_2_6P17 M_1_6P17 M_2_6P17
# FALSE    39387    27088    21007    22940    20183
# TRUE       709     2924     3769     3200     2486


x11(width = 14)
SpatialFeaturePlot(allcells, features = "nFeature_Spatial.Polygons", ncol = 4, pt.size.factor = 4) 



# Also check MT%. Read in gene info ----

geneInfo <- read.delim("data/references/GCF_022539595.1_ASM2253959v1_GeneInfo.txt")
dim(geneInfo)
# 27991     5
# slightly less than what came out of cell ranger!

nrow(allcells)
# 28004

head(geneInfo)
tail(geneInfo)
# need to replace _ with -

geneInfo$gene_id <- gsub("_", "-", geneInfo$gene_id)

sum(rownames(allcells) %in% geneInfo$gene_id)
# 27991

rownames(allcells)[!rownames(allcells) %in% geneInfo$gene_id]
# could look these up...

table(geneInfo$chr)
# named MT

MT.genes <- geneInfo$gene_id[geneInfo$chr %in% "MT"]
length(MT.genes)
# 37
sum(MT.genes %in% rownames(allcells))
# 37

#calculate percentage of UMIs in MT.genes

allcells[["percent.mt"]] <- PercentageFeatureSet(allcells, 
                                                 features = MT.genes)

summary(allcells[["percent.mt"]])
#    percent.mt     
# Min.   :  0.000  
# 1st Qu.:  1.258  
# Median :  2.996  
# Mean   :  4.925  
# 3rd Qu.:  6.167  
# Max.   :100.000

# x11(width = 10, height = 6) 
# quartz()
VlnPlot(allcells, 
        features = c("percent.mt"), 
        pt.size = 0,
        group.by = "orig.ident",
        cols = dittoColors())
ggsave("results/stats/VlnPlot_percentMT_2026-05-11.jpeg")


x11(width = 14)
SpatialFeaturePlot(allcells, features = "percent.mt", ncol = 4, pt.size.factor = 1) 
ggsave("results/stats/SpatialPlot_percentMT_2026-05-11.jpeg")
# definitely related to brain region

# based on the vlnplot, I would pick 25 but that might adverse affect certain regions
# Also the size of the points might be biasing my view.


allcells$MTbelow25 <- allcells$percent.mt < 25

x11(width = 14)
SpatialDimPlot(allcells, group.by = "MTbelow25", ncol = 4, pt.size.factor = 0.5) +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave("results/stats/SpatialPlot_MTunder25pct_2026-05-11.jpeg")
# Looks like its removing what I'd expect to be blank space


allcells$MTbelow50 <- allcells$percent.mt < 50
x11(width = 14)
SpatialDimPlot(allcells, group.by = "MTbelow50", ncol = 4, pt.size.factor = 0.5) +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave("results/stats/SpatialPlot_MTunder50pct_2026-05-11.jpeg")

# worst is 6P17, upper left. But this part was called "none" for the region


# Look at too high

x11()
VlnPlot(allcells, features = "nCount_Spatial.Polygons", pt.size = 0,
        cols = dittoColors())

allcells$above1000 <- allcells$nCount_Spatial.Polygons > 1000

x11(width = 14)
SpatialDimPlot(allcells, group.by = "above1000", ncol = 4, pt.size.factor = 0.5) +
  guides(fill = guide_legend(override.aes = list(size = 4)))
ggsave("results/stats/SpatialPlot_above1000_2026-05-11.jpeg")
# does appear to correlate with brain area - keep in.

# Final filtering ----

#remove cells/spots that do not have at least 25 UMIs
# or more than 50% MT
ncol(allcells)
# 463466
sum(allcells$nCount_Spatial.Polygons >= 25 &
      allcells$percent.mt < 50)
# 433999

allcells <- subset(allcells, subset = nCount_Spatial.Polygons >= 25 &
                     allcells$percent.mt < 50)
allcells
# An object of class Seurat 
# 28004 features across 433999 samples within 1 assay 
# Active assay: Spatial.Polygons (28004 features, 0 variable features)
# 4 spatial fields of view present: s_3P5.polygons s_4P5.polygons s_4P10.polygons s_6P17.polygons





x11()
SpatialFeaturePlot(allcells, features = "nCount_Spatial.Polygons", 
                   ncol = 4, pt.size.factor = 2,
                   ) 

x11()
SpatialFeaturePlot(allcells, features = "nFeature_Spatial.Polygons", 
                   ncol = 4, pt.size.factor = 4) 


# Check for NA, NaN, or infinite values

temp <-LayerData(allcells, assay = "Spatial.Polygons", layer = "counts")

sum(is.na(temp))
#0
sum(is.nan(temp))
#0
sum(is.infinite(temp))
#0

# Remove genes not present in at least 10 cells

table(rowSums(temp > 0) >= 10)
# FALSE  TRUE 
#  5641 22363

allcells <- allcells[rowSums(temp > 0) >= 10 ,]
allcells
# An object of class Seurat 
# 22363 features across 433999 samples within 1 assay 
# Active assay: Spatial.Polygons (22363 features, 0 variable features)
# 1 layer present: counts
# 4 spatial fields of view present: s_3P5.polygons s_4P5.polygons s_4P10.polygons s_6P17.polygons


# Fix the region abbreviations

allcells$Region <- factor(allcells$Region)
levels(allcells$Region)
#  [1] "apoa"  "atn"   "ATn"   "cer"   "CM"    "dc"    "Dc"    "dd"    "dl"    "Dl"   
# [11] "dlt"   "dm"    "Dm"    "Dp"    "Gn"    "LFB"   "mpoa"  "NDILI" "NDILm" "nH"   
# [21] "NLT"   "nMLF"  "none"  "NRL"   "NRP"   "ON"    "OpT"   "OT"    "PGa"   "PGc"  
# [31] "ppoa"  "PSm"   "Tla"   "TLa"   "vd"    "vl"    "VMn"   "vs"    "vv"

# Let's standarize to make the first letter capital

levels(allcells$Region) <- c("Apoa","Atn","Atn","Cer","Cm","Dc","Dc","Dd","Dl","Dl",
                             "Dlt","Dm","Dm","Dp","Gn","LFB","Mpoa","NDILI","NDILm","nH",
                             "NLT","nMLF","none","NRL","NRP","ON","OpT","OT","PGa","PGc",
                             "Ppoa","PSm","Tla","Tla","Vd","Vl","VMn","Vs","Vv")

nlevels(allcells$Region)
#34

# Build a named color vector for all 34 regions
regions <- levels(allcells$Region)  # or levels() if it's a factor
region_colors <- setNames(dittoColors()[seq_along(regions)], regions)
region_colors

# Switch "none" color with "OT"

region_colors["none"] <- "#8C8C8C"
region_colors["OT"] <- "#00F6B3"

# also make a dummy legend with all in them ----

# Dummy plot with all 34 regions — never displayed, just for the legend
dummy_df <- data.frame(x = seq_along(regions),
                       y = 1,
                       Region = factor(regions, levels = regions))

full_legend_region <- get_legend(
  ggplot(dummy_df, aes(x = x, y = y, fill = Region)) +
    geom_point(shape = 21, size = 3) +
    scale_fill_manual(values = region_colors) +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    theme(legend.position = "right")
)



p <- SpatialDimPlot(allcells, 
                    group.by = "Region", 
                    ncol = 4, 
                    pt.size.factor = 0.5,
                    cols = region_colors)
x11(width = 18)
plot_grid(p & theme(legend.position = "none"),
          full_legend_region,
          rel_widths = c(1, 0.12))
ggsave("results/stats/SpatialPlot_Region_2026-05-14.jpeg")



saveRDS(allcells, file = "results/stats/allcells_preInteg.rds")
#allcells <- readRDS("results/stats/allcells_preInteg.rds")



# See if integration needed ----
# https://satijalab.org/seurat/articles/integration_introduction#perform-integration-with-sctransform-normalized-datasets
# Steps for HD data: https://satijalab.org/seurat/articles/visiumhd_analysis_vignette#normalize-datasets





# Do standard steps for HD

allcells <- NormalizeData(allcells)
allcells <- FindVariableFeatures(allcells)
#allcells <- ScaleData(allcells), running this first uses a lot of memory


# Also according to the vignette, sketch analysis works better for rare populations:
# https://satijalab.org/seurat/articles/visiumhd_analysis_vignette#unsupervised-clustering

# optional memory access change, but will help speed things up ----
options(future.globals.maxSize = 24 * 1024^3)

allcells <- SketchData(allcells,
                       ncells = 50000,
                       method = "LeverageScore",
                       sketched.assay = "sketch",
                       features = VariableFeatures(allcells)
)

allcells
# An object of class Seurat 
# 44726 features across 433999 samples within 2 assays 
# Active assay: sketch (22363 features, 2000 variable features)
# 2 layers present: counts, data
# 1 other assay present: Spatial.Polygons
# 4 spatial fields of view present: s_3P5.polygons s_4P5.polygons s_4P10.polygons s_6P17.polygons


dim(allcells[["sketch"]]$counts)
#22363 50000

dim(allcells[["sketch"]]$data)
#22363 50000

DefaultAssay(allcells) <- "sketch"
allcells <- FindVariableFeatures(allcells)
allcells <- ScaleData(allcells, features = rownames(allcells))
allcells <- RunPCA(allcells, assay = "sketch", reduction.name = "pca.sketch")
allcells <- FindNeighbors(allcells, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
allcells <- FindClusters(allcells, resolution = 0.2, cluster.name = "seurat_cluster.sketched")

table(allcells$seurat_cluster.sketched)
#     0     1     2     3     4     5     6     7     8     9    10    11    12    13 
# 11216  7582  5327  4537  3819  3441  2839  2474  1633  1366  1340  1320  1267   659 
#    14    15    16    17 
#   381   348   246   205

# change to start at 1

levels(allcells$seurat_cluster.sketched)
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10"

# If already in numeric order, can change by:
levels(allcells$seurat_cluster.sketched) <- 1:nlevels(allcells$seurat_cluster.sketched)

# If not in numeric order (e.g., levels are 0, 1, 10, 11, 12, 2, etc.)
# do:
# allcells$seurat_cluster.sketched <- as.character(allcells$seurat_cluster.sketched) %>% as.numeric() %>% factor()
# levels(allcells$seurat_cluster.sketched) <- 1:nlevels(allcells$seurat_cluster.sketched)

# compare to UMAP
allcells <- RunUMAP(allcells, reduction = "pca.sketch", reduction.name = "umap.sketch", dims = 1:50, return.model = T)


DefaultAssay(allcells) <- "sketch"

x11()
DimPlot(allcells, reduction = "umap.sketch", label = T, 
        group.by = "seurat_cluster.sketched", cols = dittoColors()) + 
  ggtitle("Sketched clustering (50,000 cells)") +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave("results/stats/UMAP_sketch_res0.2_clusters_2026-05-14.jpeg")

x11()
DimPlot(allcells, reduction = "umap.sketch", label = T, 
        group.by = "seurat_cluster.sketched", 
        split.by = "seurat_cluster.sketched", cols = dittoColors(),
        ncol = 5) + 
  ggtitle("Sketched clustering (50,000 cells)") +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave("results/stats/UMAP_sketch_res0.2_clusters_split_2026-05-14.jpeg")



# check on uniqueness of top markers per cluster

Idents(allcells) <- allcells$seurat_cluster.sketched
allmarkers <- FindAllMarkers(allcells, only.pos = TRUE)


top10 <- allmarkers %>% group_by(cluster) %>%
  #  filter(avg_log2FC > 0) %>%  # don't need if only.pos = TRUE above
  slice_min(order_by =  p_val, n = 10, with_ties = FALSE) %>%
  ungroup()

x11()
DoHeatmap(allcells, features = top10$gene, 
          group.by = "seurat_cluster.sketched",
          assay = "sketch", slot = "scale.data", group.colors = dittoColors())
ggsave("results/stats/heatmap_sketch_res0.2_clusters.jpeg")
# 1 not super distinct 

# stop here and save

save.image("results/stats/temp.RData")


# check by sample:

x11()
DimPlot(allcells, reduction = "umap.sketch", label = F, 
        group.by = "Sample", cols = dittoColors()) + 
  ggtitle("Sketched clustering (50,000 cells)") +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave("results/stats/UMAP_sketch_samples_2026-05-14.jpeg")

x11(width=16, height = 16)
DimPlot(allcells, reduction = "umap.sketch", label = F, 
        group.by = "Sample", split.by = "Sample", cols = dittoColors(),
        ncol = 4) + 
  ggtitle("Sketched clustering (50,000 cells)") +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave("results/stats/UMAP_sketch_samples_split_2026-05-14.jpeg")
# this last one takes a very long time:
# Warning message:
# Removed 383999 rows containing missing values or values outside the scale range
# (`geom_scattermore()`).


# check by Region:

x11()
DimPlot(allcells, reduction = "umap.sketch", label = F, 
        group.by = "Region", cols = dittoColors()) + 
  ggtitle("Sketched clustering (50,000 cells)") +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave("results/stats/UMAP_sketch_region_2026-05-14.jpeg")
# a little separation by some regions

x11(width=18, height = 16)
DimPlot(allcells, reduction = "umap.sketch", label = F, 
        group.by = "Region", split.by = "Region", cols = dittoColors(),
        ncol = 6) + 
  ggtitle("Sketched clustering (50,000 cells)") +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave("results/stats/UMAP_sketch_region_split_2026-05-14.jpeg")





# now project to rest of cells ----


allcells <- ProjectData(
  object = allcells,
  assay = "Spatial.Polygons",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_clusters.projected = "seurat_cluster.sketched")
)

# Check the output

table(allcells$seurat_clusters.projected)
#      1     10     11     12     13     14     15     16     17     18      2      3 
# 140717   8061  10926  11063  10003   2435   1696   3149   1428    571  58525  27954 
# 4      5      6      7      8      9 
# 37847  23202  36097  29550  18602  12173
# Not in order either!!!

levels(allcells$seurat_clusters.projected)
# NULL

# not a factor?
class(allcells$seurat_clusters.projected)
# character

# [only works if already factor in order] If already in numeric order, can change by:
#levels(allcells$seurat_clusters.projected) <- 1:nlevels(allcells$seurat_clusters.projected)

# If not in numeric order (e.g., levels are 0, 1, 10, 11, 12, 2, etc.)
# do:
allcells$seurat_clusters.projected <- as.character(allcells$seurat_clusters.projected) %>% as.numeric() %>% factor()
levels(allcells$seurat_clusters.projected) <- 1:nlevels(allcells$seurat_clusters.projected)

table(allcells$seurat_cluster.sketched)
#     1     2     3     4     5     6     7     8     9    10    11    12    13    14 
# 23154  7784  5283  5117  2458  1226  1062   816   696   652   606   353   342   246 
# 15 
# 205
table(allcells$seurat_clusters.projected)
#    1     2     3     4     5     6     7     8     9    10    11    12    13    14 
# 11216  7582  5327  4537  3819  3441  2839  2474  1633  1366  1340  1320  1267   659 
# 15    16    17    18 
# 381   348   246   205

# Visualize to full projected dataset
DefaultAssay(allcells) <- "Spatial.Polygons"
Idents(allcells) <- "seurat_clusters.projected"

# Now can use dittoDimPlot for better colors
x11(width = 8, height = 7)
dittoDimPlot(allcells, var = "seurat_clusters.projected", reduction.use = "full.umap.sketch", 
             do.label = TRUE, labels.highlight = TRUE)
ggsave("results/stats/UMAP_projected_res0.2_clusters_2026-05-14.jpeg", width = 8, height = 8)

x11(width = 16, height = 15)
dittoDimPlot(allcells, var = "seurat_clusters.projected", reduction.use = "full.umap.sketch", 
             split.by = "seurat_clusters.projected", do.label = FALSE, labels.highlight = TRUE)
ggsave("results/stats/UMAP_projected_res0.2_clusters_split_2026-05-14.jpeg", width = 8, height = 8)


#See if integration is needed:
#x11()
dittoDimPlot(allcells, var = "orig.ident", reduction.use = "full.umap.sketch",
             do.label = FALSE, labels.highlight = FALSE)
#ggsave("results/stats/UMAP_projected_samples_2026-05-14.jpeg")

# looks fairly overlapping but hard to see with so many points

x11(width = 16, height = 15)
dittoDimPlot(allcells, var = "orig.ident", reduction.use = "full.umap.sketch",
             do.label = FALSE, labels.highlight = FALSE,
             split.by = "orig.ident", split.ncol = 4)
ggsave("results/stats/UMAP_projected_samples_split_2026-05-14.jpeg")




# Create down-sampled object to make visualization easier
DefaultAssay(allcells) <- "Spatial.Polygons"
Idents(allcells) <- "seurat_clusters.projected"
allcells_subset <- subset(allcells, cells = Cells(allcells[["Spatial.Polygons"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(allcells_subset) <- "Spatial.Polygons"
Idents(allcells_subset) <- "seurat_clusters.projected"
allcells_subset <- BuildClusterTree(allcells_subset, assay = "Spatial.Polygons", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(allcells_subset, assay = "Spatial.Polygons", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
#  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10proj

allcells_subset <- ScaleData(allcells_subset, assay = "Spatial.Polygons", features = top10proj$gene)
p <- DoHeatmap(allcells_subset, assay = "Spatial.Polygons", features = top10proj$gene, 
               size = 2.5, group.colors = dittoColors()) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
x11(width = 10, height = 8)
p
ggsave("results/stats/heatmap_projected_res0.2_clusters_1000sub_2026-05-14.jpeg", 
       plot = p, width = 10, height = 8)
# Still doesn't look great



# Also do clustering at 0.4 to get more ----

DefaultAssay(allcells) <- "sketch"
allcells <- FindClusters(allcells, resolution = 0.4, cluster.name = "seurat_cluster.sketched0.4")

table(allcells$seurat_cluster.sketched0.4)
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19 
#  7494 3976 3860 3206 3056 2798 2752 2697 2508 1643 1608 1371 1370 1310 1308 1242  989  964  844  741 
#   20   21   22   23   24   25   26   27   28   29 
#  680  661  642  443  411  386  353  246  236  205 

# change to start at 1

levels(allcells$seurat_cluster.sketched0.4)
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10"

# If already in numeric order, can change by:
levels(allcells$seurat_cluster.sketched0.4) <- 1:nlevels(allcells$seurat_cluster.sketched0.4)
Idents(allcells) <- allcells$seurat_cluster.sketched0.4


# now project 0.4 to rest of cells ----

allcells <- ProjectData(
  object = allcells,
  assay = "Spatial.Polygons",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_clusters.projected0.4 = "seurat_cluster.sketched0.4")
)

# Check the output

table(allcells$seurat_clusters.projected0.4)
# Not in order either!!!

levels(allcells$seurat_clusters.projected0.4)
# NULL

# not a factor?
class(allcells$seurat_clusters.projected0.4)
# character

# [only works if already factor in order] If already in numeric order, can change by:
#levels(allcells$seurat_clusters.projected) <- 1:nlevels(allcells$seurat_clusters.projected)

# If not in numeric order (e.g., levels are 0, 1, 10, 11, 12, 2, etc.)
# do:
allcells$seurat_clusters.projected0.4 <- as.character(allcells$seurat_clusters.projected0.4) %>% as.numeric() %>% factor()
levels(allcells$seurat_clusters.projected0.4) <- 1:nlevels(allcells$seurat_clusters.projected0.4)

table(allcells$seurat_clusters.projected0.4)
#     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16 
# 64262 62098 38291 26543 15598 32971 29109 15118 14442 12377 12410  8275  8230  8812 11308  9320 
# 17    18    19    20    21    22    23    24    25    26    27    28    29    30 
# 11946  7326  9076  3248  5413  5890  5491  3437  3832  1718  3237  1428  2217   576

table(allcells$seurat_clusters.projected, allcells$seurat_clusters.projected0.4)


# Visualize to full projected dataset
DefaultAssay(allcells) <- "Spatial.Polygons"
Idents(allcells) <- "seurat_clusters.projected0.4"

# Now can use dittoDimPlot for better colors
x11(width = 8, height = 7)
dittoDimPlot(allcells, var = "seurat_clusters.projected0.4", reduction.use = "full.umap.sketch", 
             do.label = TRUE, labels.highlight = TRUE)
ggsave("results/stats/UMAP_projected_res0.4_clusters_2026-05-14.jpeg", width = 8, height = 8)

x11(width = 16, height = 15)
dittoDimPlot(allcells, var = "seurat_clusters.projected0.4", reduction.use = "full.umap.sketch", 
             split.by = "seurat_clusters.projected0.4", do.label = FALSE, labels.highlight = TRUE)
ggsave("results/stats/UMAP_projected_res0.4_clusters_split_2026-05-14.jpeg", width = 8, height = 8)






# check QC metrics per cluster

x11(width = 15)
VlnPlot(allcells,features = "nCount_Spatial.Polygons", group.by = "seurat_clusters.projected", pt.size = 0, cols = dittoColors())

x11(width = 15)
VlnPlot(allcells,features = "nFeature_Spatial.Polygons", group.by = "seurat_clusters.projected", pt.size = 0, cols = dittoColors())

# cluster 10  a bit low for both. Where are these?


# Build a named color vector for all clusters
clusters <- levels(allcells$seurat_clusters.projected)  # or levels() if it's a factor
cluster_colors <- setNames(dittoColors()[seq_along(clusters)], clusters)
cluster_colors

# # Switch "none" color with "OT"
# 
# cluster_colors["none"] <- "#8C8C8C"
# cluster_colors["OT"] <- "#00F6B3"

# also make a dummy legend with all in them ----

# Dummy plot with all clusters — never displayed, just for the legend
dummy_cluster <- data.frame(x = seq_along(clusters),
                       y = 1,
                       Cluster = factor(clusters, levels = clusters))

full_legend_cluster <- get_legend(
  ggplot(dummy_cluster, aes(x = x, y = y, fill = Cluster)) +
    geom_point(shape = 21, size = 3) +
    scale_fill_manual(values = cluster_colors) +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    theme(legend.position = "right")
)



p <- SpatialDimPlot(allcells, 
                    group.by = "seurat_clusters.projected", 
                    ncol = 4, 
                    pt.size.factor = 0.5,
                    cols = cluster_colors)
x11(width = 18)
plot_grid(p & theme(legend.position = "none"),
          full_legend_cluster,
          rel_widths = c(1, 0.12))
ggsave("results/stats/SpatialPlot_res0.2_clusters_2026-05-14.jpeg")



# Do for 0.4 clusters ----

# Build a named color vector for all clusters
clusters <- levels(allcells$seurat_clusters.projected0.4)  # or levels() if it's a factor
cluster_colors <- setNames(dittoColors()[seq_along(clusters)], clusters)
cluster_colors

# # Switch "none" color with "OT"
# 
# cluster_colors["none"] <- "#8C8C8C"
# cluster_colors["OT"] <- "#00F6B3"

# also make a dummy legend with all in them ----

# Dummy plot with all clusters — never displayed, just for the legend
dummy_cluster <- data.frame(x = seq_along(clusters),
                            y = 1,
                            Cluster = factor(clusters, levels = clusters))

full_legend_cluster <- get_legend(
  ggplot(dummy_cluster, aes(x = x, y = y, fill = Cluster)) +
    geom_point(shape = 21, size = 3) +
    scale_fill_manual(values = cluster_colors) +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    theme(legend.position = "right")
)



p <- SpatialDimPlot(allcells, 
                    group.by = "seurat_clusters.projected0.4", 
                    ncol = 4, 
                    pt.size.factor = 0.5,
                    cols = cluster_colors)
x11(width = 18)
plot_grid(p & theme(legend.position = "none"),
          full_legend_cluster,
          rel_widths = c(1, 0.12))
ggsave("results/stats/SpatialPlot_res0.4_clusters_2026-05-14.jpeg")




# Write  .csv files of UMAP and clusters for Loupe ----
# NOTE: for HD will need separate files per slice and will need to
# take of bin endings that were added.

# the name given to the clusters in the dataframe below is what will be seen in Loupe

temp <- data.frame(Barcode = colnames(allcells),
                   Seurat_res0.2 = allcells$seurat_clusters.projected,
                   Seurat_res0.4 = allcells$seurat_clusters.projected0.4,
                   Sample = allcells$Sample,
                   Slice = allcells$Slice,
                   Sex = allcells$Sex,
                   Region = allcells$Region)
#change cluster numbers to "C01", etc.

temp$Seurat_res0.2 <- temp$Seurat_res0.2 %>% as.character() %>% as.numeric() %>% sprintf("c_%02d", .)
temp$Seurat_res0.4 <- temp$Seurat_res0.4 %>% as.character() %>% as.numeric() %>% sprintf("c_%02d", .)

head(temp)
tail(temp)

# remove the last 2 characters* of temp$Barcode
# * this only works if there are 9 or fewer samples

table(nchar(temp$Barcode))
# all 20

temp$Barcode <- substring(temp$Barcode, 1, 18)

# Do for all 4 Slices

ids <- unique(temp$Slice)

for(i in ids){
  temp2 <- temp[temp$Slice %in% i,]
  write.csv(temp2, file = paste0("results/stats/ForLoupe_Seurat_clusters_",i,"_2026_05_15.csv"),
            row.names = FALSE)
}


# Do the same for UMAP values

umap.vals <- Embeddings(object = allcells, reduction = "full.umap.sketch") %>% data.frame()
head(umap.vals)
umap.vals <- tibble::rownames_to_column(umap.vals, "Barcode")
umap.vals$Slice <- allcells$Slice

head(umap.vals)
tail(umap.vals)

umap.vals$Barcode <- substring(umap.vals$Barcode, 1, 18)

ids <- unique(umap.vals$Slice)

# for UMAP values, what you name the .csv file is what you see in Loupe

for(i in ids){
  temp2 <- umap.vals[umap.vals$Slice %in% i,]
  write.csv(temp2[,1:3], file = paste0("results/stats/ForLoupe_Seurat_UMAP_",i,"_2026_05_15.csv"),
            row.names = FALSE)
}





# Save object ----
saveRDS(allcells, file = "results/stats/allcells_norm_noInteg_2026-05-14.rds")
#allcells <- readRDS("results/stats/allcells_norm_noInteg_2026-05-14.rds")


# Jenny stopped here 2026-05-15 ----


# Codes to add additional cell annotation info ----

# Load in the Seurat object

allcells <- readRDS("results/stats/allcells_norm_noInteg_2026-05-14.rds")

# see cell metadata:
head(allcells)
tail(allcells)

# when the 4 capture areas were merged into one Seurat object,
# it automatically added _1, _2, _3 and _4 to cellid ends.
# The merge order was:
# 1 = s_3P, 
# 2 = s_4P5, 
# 3 = s_4P10, 
# 4 = s_6P17



# load .csv file output from Loupe
# for example, use previous ones

brains_3P5 <- read.csv("results/spaceranger/Cells_IDs_per_sample/3P5_A1_brains.csv")
brains_4P5 <- read.csv("results/spaceranger/Cells_IDs_per_sample/4P5_D1_brains.csv")
brains_4P10 <- read.csv("results/spaceranger/Cells_IDs_per_sample/4P10_A1_brains.csv")
brains_6P17 <- read.csv("results/spaceranger/Cells_IDs_per_sample/6P17_D1_brains.csv")


# add on the endings and then rbind together

head(brains_3P5)
# #            Barcode                  brains
# 1 cellid_000000009-1 3P5_A1_cell_F_1_AtlasNA
# 2 cellid_000000011-1 3P5_A1_cell_F_1_AtlasNA
# 3 cellid_000000012-1 3P5_A1_cell_F_1_AtlasNA
# 4 cellid_000000015-1 3P5_A1_cell_F_1_AtlasNA
# 5 cellid_000000016-1 3P5_A1_cell_F_1_AtlasNA
# 6 cellid_000000017-1 3P5_A1_cell_F_1_AtlasNA

head(brains_3P5)
head(brains_3P5)
head(brains_3P5)

# check/fix any differences in column names!


brains_3P5$Barcode <- paste0(brains_3P5$Barcode, "_1")
brains_4P5$Barcode <- paste0(brains_4P5$Barcode, "_2")
brains_4P10$Barcode <- paste0(brains_4P10$Barcode, "_3")
brains_6P17$Barcode <- paste0(brains_6P17$Barcode, "_4")

brains_all <- rbind(brains_3P5, brains_4P5, brains_4P10, brains_6P17)
rownames(brains_all) <- brains_all$Barcode


# If adding only 1 column, you can simply add it on:

# check current metadata column names:

names(allcells[[]])
#using one of these will replace the data in.
#using a new name will add it to the end

allcells$newname <- brains_all[colnames(allcells), "brains"]
head(allcells)












# codes below here not done but could be useful ----

# Check numbers/proportions of samples in each cluster
# using pipes:

NumCells <- table(allcells$seurat_clusters.projected, 
                  allcells$orig.ident) %>% 
  as.data.frame() 
head(NumCells)
colnames(NumCells) <- c("Cluster","Sample","NCells")
head(NumCells)

# get expected proportions per sample, which is the proportion of total cells

temp <- table(allcells$orig.ident) / ncol(allcells)
temp
#  s_3P5 Con_324_right  Con_341_left Con_341_right  Inf_335_left Inf_335_right  Inf_347_left Inf_347_right 
#     0.1135685     0.1397118     0.1321221     0.1448431     0.1179458     0.1171777     0.1178490     0.1167819


# Now use ggplot2 to make different bar plots:

x11(width = 10, height = 6)
# quartz(width = 10, height = 6)
ggplot(NumCells, aes(x = Cluster, 
                     y = NCells, 
                     fill = Sample)) +
  geom_bar(position="fill", stat="identity") +
  ylab("Proportion of cells in cluster") +
  geom_hline(yintercept = cumsum(rev(temp))) +
  scale_fill_manual(values = 1:8)
ggsave("results/stats/Cluster-Sample-Proportions_byCluster_2026-05-14.jpeg")
# clusters 1-4 almost perfect
# cluster 4 shifted but this is the yolk sac which isn't in all sections
# cluster 5 & 6 may be shifted towards controls - larger area?
# 10 & 11 are so small that props don't matter



#Do by sample

# don't add expected lines - they are too distracting
# temp <- table(allcells$seurat_clusters.projected) / ncol(allcells)

x11(width = 10, height = 6)
# quartz(width = 10, height = 6)
ggplot(NumCells, aes(x = Sample, 
                     y = NCells, 
                     fill = Cluster)) +
  geom_bar(position="fill", stat="identity") +
  ylab("Proportion of cells in cluster") +
#  geom_hline(yintercept = cumsum(rev(temp))) +
  scale_fill_manual(values = dittoSeq::dittoColors())
ggsave("results/stats/Cluster-Sample-Proportions_bySample_2026-05-14.jpeg")
# only diff in yellow obvious


# sum by treatment

NumCells <- table(allcells$seurat_clusters.projected, 
                  allcells$Treatment) %>% 
  as.data.frame() 
head(NumCells)
colnames(NumCells) <- c("Cluster","Treatment","NCells")
head(NumCells)

# get expected proportions per sample, which is the proportion of total cells

temp <- table(allcells$Treatment) / ncol(allcells)
temp
#    Control Influenza 
#  0.5302456 0.4697544


# Now use ggplot2 to make different bar plots:

x11(width = 10, height = 6)
# quartz(width = 10, height = 6)
ggplot(NumCells, aes(x = Cluster, 
                     y = NCells, 
                     fill = Treatment)) +
  geom_bar(position="fill", stat="identity") +
  ylab("Proportion of cells in cluster") +
  geom_hline(yintercept = cumsum(rev(temp))) +
  scale_fill_manual(values = c(2,4))
ggsave("results/stats/Cluster-Treatment-Proportions_byCluster_2026-05-14.jpeg")

# yes 5-9 slighly biased towards controls


#Do by sample

# don't add expected lines - they are too distracting
# temp <- table(allcells$seurat_clusters.projected) / ncol(allcells)

x11(width = 10, height = 6)
# quartz(width = 10, height = 6)
ggplot(NumCells, aes(x = Treatment, 
                     y = NCells, 
                     fill = Cluster)) +
  geom_bar(position="fill", stat="identity") +
  ylab("Proportion of cells in cluster") +
  #  geom_hline(yintercept = cumsum(rev(temp))) +
  scale_fill_manual(values = dittoSeq::dittoColors())
ggsave("results/stats/Cluster-Treatment-Proportions_byTreatment_2026-05-14.jpeg")
# can see combined effects of 5-9


# Go back and check on difference in QC metrics ----


x11()
VlnPlot(allcells, features = "nCount_Spatial.Polygons", group.by = "orig.ident", pt.size = 0,
        cols = 1:8)
# all control have more than all infected


x11()
VlnPlot(allcells, features = "nCount_Spatial.Polygons", group.by = "Treatment", pt.size = 0,
        cols = c(2,4))

# check per cluster, split by Treatment

x11(width = 12)
VlnPlot(allcells, features = "nCount_Spatial.Polygons", group.by = "seurat_clusters.projected", pt.size = 0,
        cols = dittoColors(), split.by = "Treatment")
ggsave("results/stats/VlnPlot_totalUMI_byCluster_splitTrt_2026-05-14.jpeg")




# PLOTS 

#Color images by clusters ----

# Make sure you are using the correct clusters
allcells$seurat_clusters <- allcells$seurat_clusters.projected %>% as.numeric() %>% factor()
Idents(allcells) <- "seurat_clusters"
table(allcells$seurat_clusters)
#      1      2      3      4      5      6      7      8      9     10     11 
# 488659 166403 217363  55877  35571  38163  35231  17496   8565    260     87

table(allcells$seurat_clusters, allcells$orig.ident)
# Some 0 once you hit cluster 10


#Hack to make cluster colors match UMAP plot;
#can also independently adjust pt.size.factor for each slice
#if want to call SpatialDimPlot(crop = TRUE) to maximize each tissue
plotstuff <- data.frame(image = unique(allcells$Slice), order = 1:8,
                        #place = c(1:7,9:11, 13:15),
                        pt.size.factor = 2.5)
plotstuff
#           image order pt.size.factor
# 1  s_3P5     1              1
# 2 Con_324_right     2              1
# 3  Con_341_left     3              1
# 4 Con_341_right     4              1
# 5  Inf_335_left     5              1
# 6 Inf_335_right     6              1
# 7  Inf_347_left     7              1
# 8 Inf_347_right     8              1

# # if would want to increase s_3P5 by a little 
plotstuff[1,3] <- 3

#get named vector of colors 
mycols <- dittoColors()[c(1:nlevels(allcells$seurat_clusters))]
names(mycols) <- levels(allcells$seurat_clusters)

#Do image plot with cluster colors:
plots.temp <- list()

for(i in 1:nrow(plotstuff)) {
  j <-as.character(plotstuff$image[i])
  temp <- levels(allcells$seurat_clusters[allcells$orig.ident == j, drop = TRUE]) %>% as.numeric()
  p <- SpatialDimPlot(allcells, #group.by = "seurat_clusters", 
                      cols = mycols[temp], shape = 22, crop = TRUE,
                      pt.size.factor = plotstuff$pt.size.factor[i],
                      images = paste0(j, ".008um"), label = TRUE, label.size = 3) + NoLegend()
  plots.temp[[i]] <- p + ggtitle(j)
}

x11(width = 28, height = 14)
wrap_plots( plots.temp, ncol = 4)
#NOTE: the size of the point in the saved graph below will
#be larger than what you see in R. Play around with the 
#pt.size.factor in plotstuff to get the desired size in
#the output .jpeg
ggsave("results/stats/Tissue_with_Clusters_withLabels_2026-05-14.jpeg")


# do without labels
for(i in 1:nrow(plotstuff)) {
  j <- as.character(plotstuff$image[i])
  temp <- levels(allcells$seurat_clusters[allcells$orig.ident == j, drop = TRUE]) %>% as.numeric()
  p <- SpatialDimPlot(allcells, #group.by = "seurat_clusters", 
                      cols = mycols[temp], shape = 22, crop = TRUE,
                      pt.size.factor = plotstuff$pt.size.factor[i],
                      images =paste0(j, ".008um"), label = FALSE, label.size = 2) + NoLegend()
  plots.temp[[i]] <- p + ggtitle(j)
}

x11(width = 28, height = 14)
wrap_plots( plots.temp, ncol = 4)
#NOTE: the size of the point in the saved graph below will
#be larger than what you see in R. Play around with the 
#pt.size.factor in plotstuff to get the desired size in
#the output .jpeg
ggsave("results/stats/Tissue_with_Clusters_noLabels_2026-05-14.jpeg")


#plot location of clusters, each cluster to one figure
# skip clusters 10-11

#for(k in 1:nlevels(allcells$seurat_clusters)) {
for(k in 1:9) {
    plots.temp <- list()
  for(i in 1:nrow(plotstuff)) {
    j <- as.character(plotstuff$image[i])
    cols.highlight = c(mycols[as.character(k)],"lightgray")
    names(cols.highlight) <- NULL
    p <- SpatialDimPlot(allcells, cells.highlight = CellsByIdentities(object = allcells, idents = k), 
                        cols.highlight = cols.highlight,
                        facet.highlight = TRUE, shape = 22, crop = TRUE,
                        pt.size.factor = plotstuff$pt.size.factor[i],
                        images = paste0(j, ".008um"), label = FALSE, label.size = 2) + NoLegend()
    plots.temp[[i]] <- p + ggtitle(j)
  }
  
    x11(width = 28, height = 14)
  print(wrap_plots( plots.temp, ncol = 4))
  ggsave(paste0("results/stats/cluster_",k,"_highlighted_2026-05-14.jpeg"))
  dev.off()
}



#plot location of clusters, each sample to one figure

for(j in plotstuff$image) {
  plots.temp <- list()
  temp <- levels(allcells$seurat_clusters[allcells$orig.ident == j, drop = TRUE]) %>% as.numeric()
#  for(i in 1:(nlevels(allcells$seurat_clusters))){
  for(i in 1:9) {
    
    if(i %in% temp){
      p <- SpatialDimPlot(allcells, cells.highlight = CellsByIdentities(object = allcells, idents = i), 
                          facet.highlight = TRUE, images =  paste0(j, ".008um"), cols.highlight = c("blue", "lightgrey"), 
                          pt.size.factor = plotstuff$pt.size.factor[plotstuff$image == j])
      plots.temp[[i ]] <- p
    } else {
      p <-  plot_spacer()
      plots.temp[[i ]] <- p
    }
  }
  x11(width = 14, height = 14)
  print(wrap_plots(plots.temp, ncol = 3))
  ggsave(paste0("results/stats/",j,"_highlighted_2026-05-14.jpeg"))
  dev.off()
}


# Within each cluster, do DE testing between treatments ----

# with only 2 groups, do pseudo-bulking 
# https://satijalab.org/seurat/articles/de_vignette#perform-de-analysis-after-pseudobulking

DefaultAssay(allcells) <- "Spatial.Polygons"

# remove cells of clusters 10 & 11

temp <- allcells[, allcells$seurat_clusters %in% 1:9]
temp

pseudocells <- AggregateExpression(temp, assays = "Spatial.Polygons", return.seurat = T, 
                                   group.by = c("Treatment", "orig.ident", "seurat_clusters"))
pseudocells
# An object of class Seurat 
# 18989 features across 72 samples within 1 assay 
# Active assay: Spatial.Polygons (18989 features, 0 variable features)
# 3 layers present: counts, data, scale.data


rm(temp)
gc()


tail(Cells(pseudocells))

head(pseudocells)

# Make combination of treatment and cluster

pseudocells$TrtClust <- paste(pseudocells$Treatment, pseudocells$seurat_clusters, sep = "_")
unique(pseudocells$TrtClust)

# factor into levels...
pseudocells$TrtClust <- factor(pseudocells$TrtClust, levels = paste(c("Control", "Influenza"), rep(1:9, each = 2), sep = "_"))
Idents(pseudocells) <- "TrtClust"

# Do for cluster 1 at first

i <- 1
Inf.vs.Con <- FindMarkers(pseudocells, 
                          ident.1 = paste0("Influenza_", i),
                          ident.2 = paste0("Control_", i), 
                          test.use = "DESeq2")
head(Inf.vs.Con)

Inf.vs.Con$cluster <- i
Inf.vs.Con <- tibble::rownames_to_column(Inf.vs.Con, var = "gene")

# check on 1 or 2 genes
allcells$TrtClust <- paste(allcells$Treatment, allcells$seurat_clusters, sep = "_")
allcells$TrtClust <- factor(allcells$TrtClust, levels = paste(c("Control", "Influenza"), rep(1:9, each = 2), sep = "_"))
Idents(allcells) <- "TrtClust"

x11()
VlnPlot(allcells, features = Inf.vs.Con$gene[1:2], 
        idents = c("Control_1", "Influenza_1"), group.by = "Treatment") 
# Somewhat compelling...

# do for rest

for(i in 2:9){
  temp <- FindMarkers(pseudocells, 
                         ident.1 = paste0("Influenza_", i),
                         ident.2 = paste0("Control_", i), 
                         test.use = "DESeq2")
  temp$cluster <- i
  temp <- tibble::rownames_to_column(temp, var = "gene")
  Inf.vs.Con <- rbind(Inf.vs.Con, temp)
}

Inf.vs.Con <- select(Inf.vs.Con, cluster, gene, avg_log2FC, p_val, p_val_adj, everything())
write.table(Inf.vs.Con, file = "results/stats/Inf.vs.Con_0.3res_clusters_2026-05-14.txt",
            row.names = FALSE, sep = "\t")
#Inf.vs.Con <- read.delim("results/stats/Inf.vs.Con_0.3res_clusters_2026-05-14.txt")

table(Inf.vs.Con$cluster)
#    1     2     3     4     5     6     7     8     9 
#18988 18973 18960 17720 18533 18783 18621 17837 15138

# only genes not tested are ones that have 0 for all samples


table(Inf.vs.Con$p_val_adj < 0.05,Inf.vs.Con$cluster )
#            1     2     3     4     5     6     7     8     9
# FALSE  17703 18835 18442 17714 18504 18595 18524 17787 15131
# TRUE    1284   136   514     5     7   188    97    50     7

# cluster 1 has the most, as would be expected from the first clustering 
# where they were different. 


# Also do regular testing to compare (and get proportions)
# set thresholds to zero to get all genes back

i <- 1
Inf.vs.Con.bad <- FindMarkers(allcells, paste0("Influenza_", i),
                              ident.2 = paste0("Control_", i), 
                              min.pct = 0, 
                              logfc.threshold = 0, 
                              assay = "Spatial.Polygons", slot = "data",
                              verbose = TRUE)

dim(Inf.vs.Con.bad)
#18989     5

head(Inf.vs.Con.bad)

Inf.vs.Con.bad$cluster <- i
Inf.vs.Con.bad <- tibble::rownames_to_column(Inf.vs.Con.bad, var = "gene")

for(i in 2:9){
  temp <- FindMarkers(allcells, 
                      ident.1 = paste0("Influenza_", i),
                      ident.2 = paste0("Control_", i), 
                      min.pct = 0, 
                      logfc.threshold = 0, 
                      assay = "Spatial.Polygons", slot = "data",
                      verbose = TRUE)
  temp$cluster <- i
  temp <- tibble::rownames_to_column(temp, var = "gene")
  Inf.vs.Con.bad <- rbind(Inf.vs.Con.bad, temp)
}

Inf.vs.Con.bad <- select(Inf.vs.Con.bad, cluster, gene, avg_log2FC, p_val, p_val_adj, everything())
write.table(Inf.vs.Con.bad, file = "results/stats/Inf.vs.Con_0.3res_clusters_perCell_2026-05-14.txt",
            row.names = FALSE, sep = "\t")
#Inf.vs.Con.bad <- read.delim("results/stats/Inf.vs.Con_0.3res_clusters_perCell_2026-05-14.txt")

table(Inf.vs.Con.bad$cluster)
#     1     2     3     4     5     6     7     8     9 
# 18989 18989 18989 18989 18989 18989 18989 18989 18989

# all tested with limits set to 0

tail(Inf.vs.Con.bad)


table(Inf.vs.Con.bad$p_val_adj < 0.05,Inf.vs.Con.bad$cluster )
#           1     2     3     4     5     6     7     8     9
# FALSE  7506  9650  9243 15853 14967 11657 13141 15152 18151
# TRUE  11483  9339  9746  3136  4022  7332  5848  3837   838

# of course, 0.05 is not the appropriate threshold

table(Inf.vs.Con.bad$p_val_adj < 1e-100,Inf.vs.Con.bad$cluster )
#           1     2     3     4     5     6     7     8     9
# FALSE 14151 17911 15825 18946 18977 18905 18924 18982 18989
# TRUE   4838  1078  3164    43    12    84    65     7     0


View(Inf.vs.Con.bad[Inf.vs.Con.bad$gene %in% Inf.vs.Con$gene[1:2],])
# In all, but low pct in both
# This might be an example of the effects lower transcripts in
# Inf, except this one has a positive FC, meaning more in Inf!!

# Put cell percentages into pseusobulk results

names(Inf.vs.Con.bad)[3:7] <- paste0("cell_", names(Inf.vs.Con.bad)[3:7])

Inf.vs.Con$clustgene <- paste(Inf.vs.Con$cluster, Inf.vs.Con$gene,sep = "_")
Inf.vs.Con.bad$clustgene <- paste(Inf.vs.Con.bad$cluster, Inf.vs.Con.bad$gene,sep = "_")

Inf.vs.Con <- left_join(Inf.vs.Con[,-c(6:7)], Inf.vs.Con.bad[,3:8])

head(Inf.vs.Con)

tail(Inf.vs.Con)

# also calcuate pct.diff

Inf.vs.Con$pct.diff <- Inf.vs.Con$cell_pct.1 - Inf.vs.Con$cell_pct.2

dim(Inf.vs.Con)
#163553     12

sum(is.na(Inf.vs.Con$p_val_adj))
#30
sum(is.na(Inf.vs.Con$cell_p_val_adj))
#0

View(Inf.vs.Con[is.na(Inf.vs.Con$p_val_adj),])

# these can be safely gotten rid of

Inf.vs.Con <- Inf.vs.Con[!is.na(Inf.vs.Con$p_val_adj),]


sum(Inf.vs.Con$p_val_adj > 0.5)
# 160255 

sum(Inf.vs.Con$cell_p_val_adj > 0.5)
#103201

sum(Inf.vs.Con$p_val_adj > 0.5 & Inf.vs.Con$cell_p_val_adj > 0.5)
#102629

# write out all results, then cut down

write.table(Inf.vs.Con, file = "results/stats/Inf.vs.Con_Perclusters_pseudoANDcell_2026-05-14.txt",
            row.names = FALSE, sep = "\t")


Inf.vs.Con <- Inf.vs.Con[!(Inf.vs.Con$p_val_adj > 0.5 & Inf.vs.Con$cell_p_val_adj > 0.5),]

dim(Inf.vs.Con)
#60894    12

# how many don't have the same direction?

table(sign(Inf.vs.Con$avg_log2FC) == sign(Inf.vs.Con$cell_avg_log2FC))
#FALSE  TRUE 
# 6697 54197

View(Inf.vs.Con[sign(Inf.vs.Con$avg_log2FC) != sign(Inf.vs.Con$cell_avg_log2FC),])
# generally these aren't genes that are significant


View(Inf.vs.Con)
# many of the most significant from the pseudobulk had very low
# pct in either group!


x11()
VlnPlot(allcells, features = "Pdia4", 
        idents = c("Control_3", "Influenza_3"), group.by = "Treatment") 
# This looks to be a case of fewer cells in Influenza expressing Pdia4,
# not a reduction in expression
x11()
VlnPlot(allcells, features = "Pdia4", pt.size = 0,
        idents = c("Control_3", "Influenza_3"), group.by = "Treatment") 
# Without points, you can now see the big difference

x11()
VlnPlot(allcells, features = "Serpine1", pt.size = 0,
        idents = c("Control_3", "Influenza_3"), group.by = "Treatment") 





# Seurat's pseudobulking doesn't seem to take into account how many
# cells there were. Does muscat?


library(SingleCellExperiment)
library(muscat)
library(DropletUtils)
library(scuttle)
library(scater)
library(scran)
library(dplyr)
library(purrr)
library(tidyr)
library(patchwork)

# muscat ----

temp <-LayerData(allcells, assay = "Spatial.Polygons", layer = "counts")


sce <- SingleCellExperiment(assays = list(counts = temp),
                           colData=DataFrame(allcells[[]]))

sce <- computeLibraryFactors(sce)

summary(sizeFactors(sce))

allcells$sizeFactors <- sizeFactors(sce)

x11()
VlnPlot(allcells, "sizeFactors", group.by = "orig.ident", pt.size = 0)
# looks very similar to total UMI counts - higher in controls


sce <- logNormCounts(sce)

temp <- logcounts(sce)

allcells$sceNormSum <- colSums(temp)
x11()
VlnPlot(allcells, "sceNormSum", group.by = "orig.ident", pt.size = 0)
#also still higher in Con

# remove clusters 10 and 11

sce <- sce[,sce$seurat_clusters %in% 1:9]
sce$seurat_clusters <- droplevels(sce$seurat_clusters)

levels(sce$seurat_clusters)


# prep for muscat by saying which factors are the clusters, treatment and sample ids
sce <- prepSCE(sce, 
                kid = "seurat_clusters", # cluster assignments
                gid = "Treatment",  # group IDs (ctrl/stim)
                sid = "orig.ident",   # sample IDs (ctrl/stim.1234)
                drop = FALSE)

# add c to beginning of clusters so they don't start with a number

levels(sce$cluster_id) <- paste0("C", 1:9)


nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# check to make sure enough bins

t(table(sce$cluster_id, sce$sample_id))
# all have over 500

# aggregate data
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per cluster
assayNames(pb)
#"C1" "C2" "C3" "C4" "C5" "C6" "C7" "C8" "C9"


# Pseudobulk-level MDS plot


pb_mds <- pbMDS(pb)
x11()
pb_mds + 
#  scale_shape_manual(values = c(17, 4)) +
  scale_color_manual(values = dittoColors())

# clusters separating; can see control-influenza split in several
# but oddly not so much in cluster 1


# Run pseudobulk analysis

res <- pbDS(pb, verbose = FALSE)
# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)


k1 <- tbl[[1]]
View(k1)
# same few top genes as seurat's pseudo bulk, although p-values are less significant

x11()
VlnPlot(allcells, features = "Arf2", pt.size = 0,
        idents = c("Control_1", "Influenza_1"), group.by = "Treatment") 
x11()
VlnPlot(allcells, features = "Arf2", pt.size = 0,
        idents = c("Control_1", "Influenza_1"), group.by = "orig.ident") 


x11()
VlnPlot(allcells, features = "Fabp3", pt.size = 0,
        idents = c("Control_1", "Influenza_1"), group.by = "Treatment") 


# Try muscat's differential detection analysis
# https://bioconductor.org/packages/release/bioc/vignettes/muscat/inst/doc/detection.html

# "a differential detection (DD) analysis aims to identify genes for which the average 
# fraction of cells in which the gene is detected changes between groups. "


pb_det <- aggregateData(sce,
                        assay="counts", fun="num.detected",
                        by=c("cluster_id", "sample_id"))
t(head(assay(pb_det)))


pb_mds_det <- pbMDS(pb_det)
x11()
pb_mds_det + 
  #  scale_shape_manual(values = c(17, 4)) +
  scale_color_manual(values = dittoColors())
# very similar to sum

res_DD <- pbDD(pb_det, min_cells=0, filter="none", verbose=FALSE)

tbl_DD <- res_DD$table[[1]]
# one data.frame per cluster
names(tbl_DD)

k1_DD <- tbl_DD[[1]]
View(k1_DD)


# Try both DS and DD together ----
# redo DS with min_cells = 0

res_DS <- pbDS(pb, min_cells=0, filter="none", verbose=FALSE)

res_both <- stagewise_DS_DD(res_DS, res_DD, verbose=FALSE)
head(res_both[[1]][[1]])

ps <- map_depth(res_both, 2, \(df) {
  data.frame(
    df[, c("gene", "cluster_id")],
    p_adj.stagewise=df$p_adj,
    p_adj.DS=df$res_DS$p_adj.loc,
    p_adj.DD=df$res_DD$p_adj.loc)
}) |> 
  lapply(do.call, what=rbind) |>
  do.call(what=rbind) |>
  data.frame(row.names=NULL)
head(ps)

tail(ps)



# for each approach & cluster, count number 
# of genes falling below 5% FDR threshold
x11(width = 17.5)
ns <- lapply(seq(0, 0.2, 0.005), \(th) {
  ps |>
    mutate(th=th) |>
    group_by(cluster_id, th) |>
    summarise(
      .groups="drop",
      across(starts_with("p_"), 
             \(.) sum(. < th, na.rm=TRUE)))
}) |> 
  do.call(what=rbind) |>
  pivot_longer(starts_with("p_"))
ggplot(ns, aes(th, value, col=name)) + 
  geom_line(linewidth=0.8, key_glyph="point") +
  geom_vline(xintercept=0.05, lty=2, linewidth=0.4) +
  guides(col=guide_legend(NULL, override.aes=list(size=3))) +
  labs(x="FDR threshold", y="number of significantly\ndifferential genes") +
  facet_wrap(~cluster_id, scales="free_y", nrow=2) + 
  theme_bw() + theme(
    panel.grid.minor=element_blank(),
    legend.key.size=unit(0.5, "lines"))


# subset adjuster p-values for cluster of interest
qs <- ps[grep("C1", ps$cluster_id), grep("p_", names(ps))]
# for each approach, extract genes at 5% FDR threshold
gs <- apply(qs, 2, \(.) ps$gene[. < 0.05])
# visualize set intersections between approaches
x11()
UpSetR::upset(UpSetR::fromList(gs), order.by="freq")
# vast majority in all 3

qs <- ps[grep("C3", ps$cluster_id), grep("p_", names(ps))]
# for each approach, extract genes at 5% FDR threshold
gs <- apply(qs, 2, \(.) ps$gene[. < 0.05])
# visualize set intersections between approaches
x11()
UpSetR::upset(UpSetR::fromList(gs), order.by="freq")
# only 11 only in stagewise - what are they?

sw <- grep("stagewise", names(gs))
setdiff(gs[[sw]], unlist(gs[-sw]))
#"Mlph"     "Gigyf2"   "Lrrc47"   "Snx16"    "Hes3"     "Dmrtb1"   "Slc25a28" "Ildr2"    "Vmn1r45"  "Tmem212"  "Dpp8"


qs <- ps[grep("C5", ps$cluster_id), grep("p_", names(ps))]
# for each approach, extract genes at 5% FDR threshold
gs <- apply(qs, 2, \(.) ps$gene[. < 0.05])
# visualize set intersections between approaches
x11()
UpSetR::upset(UpSetR::fromList(gs), order.by="freq")
# majority only in DD

# stopped here 2026-05-14 ----

saveRDS(allcells, file = "results/stats/allcells_norm_noInteg_2026-05-14.rds")
#allcells <- readRDS("results/stats/allcells_norm_noInteg_2026-05-14.rds")



# check on y genes ----

geneinfo <- read.delim("results/refdata-gex-mm10-2020-A_genes.txt")

head(geneinfo)

dim(geneinfo)
#32285  5     

sum(geneinfo$gene_name %in% rownames(allcells))
# 18998

sum(!rownames(allcells) %in% geneinfo$gene_name)
# 5

rownames(allcells)[!rownames(allcells) %in% geneinfo$gene_name]
# the .1 suggests two rows with the same gene name

grep("Pakap", rownames(allcells), value = TRUE)
# "Pakap"   "Pakap.1"

geneinfo[grep("Pakap", geneinfo$gene_name),]

sum(duplicated(geneinfo$gene_name))
# 40

table(table(geneinfo$gene_name))
#     1     2     3 
# 32207    36     2

# make these unique

geneinfo$gene_name_unique <- make.unique(geneinfo$gene_name)

sum(!rownames(allcells) %in% geneinfo$gene_name_unique)
# 0


geneinfo <- geneinfo[match(rownames(allcells), geneinfo$gene_name_unique), ]
dim(geneinfo)
#18989     6

table(geneinfo$seqnames)

Ygenes <- geneinfo$gene_name_unique[geneinfo$seqnames %in% "chrY"]
Ygenes

temp <-LayerData(allcells, assay = "Spatial.Polygons", layer = "counts")

allcells$sumYcounts <- colSums(temp[Ygenes,])

x11()
VlnPlot(allcells, features = "sumYcounts", group.by = "orig.ident", pt.size = 0)

summary(allcells$sumYcounts)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.00000 0.04259 0.00000 5.00000

table(allcells$orig.ident, allcells$sumYcounts > 0)
#                 FALSE   TRUE
# s_3P5  111313   9487
# Con_324_right 137883  10725
# Con_341_left  140507     28
# Con_341_right 143972  10094
# Inf_335_left  121549   3907
# Inf_335_right 120535   4104
# Inf_347_left  125339     14
# Inf_347_right 119703   4515

# maybe Con_341_left and Inf_347_left females and the rest males?


# very few y genes detected. Try X

Xgenes <- geneinfo$gene_name_unique[geneinfo$seqnames %in% "chrX"]
length(Xgenes)
# 700

allcells$sumXcounts <- colSums(temp[Xgenes,])
summary(allcells$sumXcounts)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.00    5.00   10.00   12.59   17.00  311.00

x11(width = 10)
VlnPlot(allcells, features = "sumXcounts", group.by = "orig.ident", pt.size = 0)

table(allcells$orig.ident, allcells$sumXcounts > 0)


# check pseudobulked data

head(pseudocells)
pseudocells$orig.ident2 <- substring(pseudocells$orig.ident, 1, nchar(pseudocells$orig.ident)-2)

x11(width = 20, height = 20)
VlnPlot(pseudocells, features = Ygenes, group.by = "orig.ident2")
ggsave("results/stats/Ygenes_expression_2026-07-14.jpeg")

# write out geneinfo again with new names

write.table(geneinfo, file = "results/stats/gene_info_filtered_2026-07-14.txt", row.names = FALSE, sep = "\t")




# package versions ----
sessionInfo()

# R version 4.5.0 (2026-04-11 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26100)
# 
# Matrix products: default
# LAPACK version 3.12.1
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/Chicago
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tidyr_1.3.1                 purrr_1.0.4                 scran_1.36.0                scater_1.36.0              
# [5] scuttle_1.18.0              DropletUtils_1.28.0         muscat_1.22.0               SingleCellExperiment_1.30.1
# [9] SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.0        
# [13] IRanges_2.42.0              S4Vectors_0.46.0            BiocGenerics_0.54.0         generics_0.1.4             
# [17] MatrixGenerics_1.20.0       matrixStats_1.5.0           dittoSeq_1.20.0             dplyr_1.1.4                
# [21] patchwork_1.3.1             ggplot2_3.5.2               Seurat_5.3.1.1000           SeuratObject_5.1.0         
# [25] sp_2.2-0                   
# 
# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.1-0     bitops_1.0-9              httr_1.4.7                RColorBrewer_1.1-3       
# [5] doParallel_1.0.17         numDeriv_2016.8-1.1       tools_4.5.0               sctransform_0.4.2        
# [9] backports_1.5.0           R6_2.6.1                  HDF5Array_1.36.0          lazyeval_0.2.2           
# [13] uwot_0.2.3                mgcv_1.9-1                rhdf5filters_1.20.0       GetoptLong_1.0.5         
# [17] withr_3.0.2               prettyunits_1.2.0         gridExtra_2.3             progressr_0.15.1         
# [21] cli_3.6.5                 Cairo_1.6-2               spatstat.explore_3.4-3    fastDummies_1.7.5        
# [25] labeling_0.4.3            mvtnorm_1.3-3             spatstat.data_3.1-6       blme_1.0-6               
# [29] ggridges_0.5.6            pbapply_1.7-2             R.utils_2.13.0            dichromat_2.0-0.1        
# [33] parallelly_1.45.0         limma_3.64.1              rstudioapi_0.17.1         shape_1.4.6.1            
# [37] gtools_3.9.5              ica_1.0-3                 spatstat.random_3.4-1     Matrix_1.7-3             
# [41] ggbeeswarm_0.7.2          abind_1.4-8               R.methodsS3_1.8.2         lifecycle_1.0.4          
# [45] edgeR_4.6.2               rhdf5_2.52.1              gplots_3.2.0              SparseArray_1.8.0        
# [49] Rtsne_0.17                grid_4.5.0                dqrng_0.4.1               promises_1.3.3           
# [53] crayon_1.5.3              miniUI_0.1.2              lattice_0.22-6            beachmat_2.24.0          
# [57] cowplot_1.1.3             metapod_1.16.0            pillar_1.11.0             ComplexHeatmap_2.24.1    
# [61] rjson_0.2.23              boot_1.3-31               corpcor_1.6.10            future.apply_1.20.0      
# [65] codetools_0.2-20          glue_1.8.0                spatstat.univar_3.1-3     data.table_1.17.6        
# [69] vctrs_0.6.5               png_0.1-8                 spam_2.11-1               Rdpack_2.6.4             
# [73] gtable_0.3.6              rbibutils_2.3             S4Arrays_1.8.1            mime_0.13                
# [77] reformulas_0.4.1          survival_3.8-3            pheatmap_1.0.13           iterators_1.0.14         
# [81] bluster_1.18.0            statmod_1.5.0             fitdistrplus_1.2-4        ROCR_1.0-11              
# [85] nlme_3.1-168              pbkrtest_0.5.4            EnvStats_3.1.0            progress_1.2.3           
# [89] RcppAnnoy_0.0.22          UpSetR_1.4.0              TMB_1.9.17                irlba_2.3.5.1            
# [93] vipor_0.4.7               KernSmooth_2.23-26        colorspace_2.1-1          ggrastr_1.0.2            
# [97] DESeq2_1.48.1             tidyselect_1.2.1          compiler_4.5.0            BiocNeighbors_2.2.0      
# [101] h5mread_1.0.1             DelayedArray_0.34.1       plotly_4.11.0             scales_1.4.0             
# [105] caTools_1.18.3            remaCor_0.0.18            lmtest_0.9-40             stringr_1.5.1            
# [109] digest_0.6.37             goftest_1.2-3             spatstat.utils_3.1-4      presto_1.0.0             
# [113] minqa_1.2.8               variancePartition_1.38.0  aod_1.3.3                 XVector_0.48.0           
# [117] RhpcBLASctl_0.23-42       htmltools_0.5.8.1         pkgconfig_2.0.3           lme4_1.1-37              
# [121] sparseMatrixStats_1.20.0  fastmap_1.2.0             rlang_1.1.6               GlobalOptions_0.1.2      
# [125] htmlwidgets_1.6.4         UCSC.utils_1.4.0          stageR_1.30.1             DelayedMatrixStats_1.30.0
# [129] shiny_1.11.1              farver_2.1.2              zoo_1.8-14                jsonlite_2.0.0           
# [133] BiocParallel_1.42.1       R.oo_1.27.1               BiocSingular_1.24.0       magrittr_2.0.3           
# [137] GenomeInfoDbData_1.2.14   dotCall64_1.2             Rhdf5lib_1.30.0           Rcpp_1.0.14              
# [141] viridis_0.6.5             reticulate_1.42.0         stringi_1.8.7             MASS_7.3-65              
# [145] plyr_1.8.9                parallel_4.5.0            listenv_0.9.1             ggrepel_0.9.6            
# [149] deldir_2.0-4              splines_4.5.0             tensor_1.5.1              hms_1.1.3                
# [153] circlize_0.4.16           locfit_1.5-9.12           igraph_2.1.4              spatstat.geom_3.4-1      
# [157] RcppHNSW_0.6.0            reshape2_1.4.4            ScaledMatrix_1.16.0       BiocManager_1.30.26      
# [161] nloptr_2.2.1              foreach_1.5.2             httpuv_1.6.16             RANN_2.6.2               
# [165] polyclip_1.10-7           future_1.58.0             clue_0.3-66               scattermore_1.2          
# [169] rsvd_1.0.5                broom_1.0.8               xtable_1.8-4              fANCOVA_0.6-1            
# [173] RSpectra_0.16-2           later_1.4.2               viridisLite_0.4.2         tibble_3.3.0             
# [177] lmerTest_3.1-3            glmmTMB_1.1.11            beeswarm_0.4.0            cluster_2.1.8.1          
# [181] globals_0.18.0