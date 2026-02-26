library(Seurat)
library(Signac)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(aricode)
library(RColorBrewer)
library(Biostrings)
library(GenomicRanges)
library(scales)
library(tidyverse)
library(scCustomize)
library(patchwork)
library(ggnewscale)
library(grid)
library(cowplot)
library(ggplotify)
library(ggrastr)
library(ggrepel)
library(CellChat)
library(reticulate)

mem.maxVSize(15 * 1024^3)
options(future.globals.maxSize = 15 * 1024^3)

setwd('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/hypo_pit')

homologs       <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/Mzebra_genome_files/Old_MZ_genome_files/mzebra_gene_homologs_annotation_info.rds')

source("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/scripts/netVisual_circle_mod.R")

# # # only run this once!!
# ---- install into the python that py_require() is forcing ----
# py <- "/storage/home/hcoda1/8/kleatherbury3/.cache/R/reticulate/uv/cache/archive-v0/F1PWlIMsxToeVI-RmP_yM/bin/python"
# 
# # Make sure pip exists
# system2(py, c("-m", "ensurepip", "--upgrade"), stdout = TRUE, stderr = TRUE)
# 
# # Upgrade pip tooling
# system2(py, c("-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"),
#         stdout = TRUE, stderr = TRUE)
# 
# # Force reinstall numpy==2.3.5 and the umap stack into THIS env
# system2(py, c("-m", "pip", "install", "--upgrade", "--force-reinstall", "numpy==2.3.5"),
#         stdout = TRUE, stderr = TRUE)
# 
# system2(py, c("-m", "pip", "install", "--upgrade", "--force-reinstall", "llvmlite", "numba", "umap-learn"),
#         stdout = TRUE, stderr = TRUE)
# 
# # OPTIONAL but recommended: remove the leftover "~umpy*" directories pip warned about
# site <- "/storage/home/hcoda1/8/kleatherbury3/.cache/R/reticulate/uv/cache/archive-v0/F1PWlIMsxToeVI-RmP_yM/lib/python3.12/site-packages"
# unlink(file.path(site, "~umpy"), recursive = TRUE, force = TRUE)
# unlink(file.path(site, "~umpy.libs"), recursive = TRUE, force = TRUE)
# 
# Sys.setenv(
#   RETICULATE_PYTHON = "/storage/home/hcoda1/8/kleatherbury3/.cache/R/reticulate/uv/cache/archive-v0/F1PWlIMsxToeVI-RmP_yM/bin/python",
#   PYTHONNOUSERSITE  = "1",
#   PYTHONPATH        = ""
# )

### ADD THIS RO .RPROFILE
## ---- Reticulate Python lock for CellChat / UMAP ----
library(reticulate)

# should point to F1PWlIM... python
py_config()

# must print numpy 2.3.5, and successfully import numba + umap
py_run_string("import numpy as np; import numba; import umap; print('numpy', np.__version__); print('numba', numba.__version__); print('umap', umap.__version__)")


############################################
## CellChat pipeline: male_control
## Retain only robust, cross-individual interactions
############################################
ptm = Sys.time()

cellchat_behave <- readRDS("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/hypo_pit/cellchat/Spawning_Female/child_clusters/cellchat_hypo_Spawning_Female_L9_Child_CellType_netCentrality.rds")
cellchat_hypo_male_behave <- readRDS("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/hypo_pit/cellchat/Spawning_Male/child_clusters/cellchat_hypo_Spawning_Male_L9_Child_CellType_netCentrality.rds")
cellchat_hypo_female_control <- readRDS("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/hypo_pit/cellchat/Non-Ovulating_Female/child_clusters/cellchat_hypo_Non-Ovulating_Female_L9_Child_CellType_netCentrality.rds")
cellchat_hypo_male_control <- readRDS("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/hypo_pit/cellchat/Control_Male/child_clusters/cellchat_hypo_Control_Male_L9_Child_CellType_netCentrality.rds")

common.groups <- intersect(levels(cellchat_hypo_male_behave@idents),levels(cellchat_behave@idents))

length(levels(cellchat_behave@idents))
length(levels(cellchat_hypo_male_behave@idents))

common.groups

object.list_group <- list(MB = cellchat_hypo_male_behave, FB = cellchat_behave, MC = cellchat_hypo_male_control, FC = cellchat_hypo_female_control)

object.list <- list()
for (nm in names(object.list_group)) {
  
  obj <- object.list_group[[nm]]
  
  lvls <- levels(obj@idents)
  group.existing <- rownames(obj@net$weight)
  print(length(lvls))
  print(length(group.exist
  
  obj@idents <- droplevels(factor(as.character(obj@idents)))
  
  obj <- netAnalysis_computeCentrality(obj, slot.name = "netP")
  object.list[[nm]] <- obj
  
}

saveRDS(object.list, "/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/hypo_pit/cellchat/HypoPit_CellChatObject_List_by_group.rds")

object.list <- readRDS("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/hypo_pit/cellchat/Spawning_Group/HypoPit_CellChatObject_List_by_group.rds")

common.groups <- intersect(levels(object.list[["MB"]]@idents),levels(object.list[["FB"]]@idents))

length(levels(object.list[["MB"]]@idents))
length(levels(object.list[["FB"]]@idents))

common.groups

object.list[["MB"]] <- liftCellChat(object.list[["MB"]], idents_ordered)
object.list[["MC"]] <- liftCellChat(object.list[["MC"]], idents_ordered)
object.list[["FB"]] <- liftCellChat(object.list[["FB"]], idents_ordered)
object.list[["FC"]] <- liftCellChat(object.list[["FC"]], idents_ordered)

cellchat_group <- mergeCellChat(object.list, add.names = names(object.list_group))
cellchat_group

saveRDS(cellchat_group, "/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/hypo_pit/cellchat/HypoPit_merged_CellChatObjects_by_group.rds")



# 
# object.list_male <- list(MB = cellchat_hypo_male_behave,
#                          MC = cellchat_hypo_male_control)
# cellchat_male <- mergeCellChat(object.list_male, add.names = names(object.list_male))
# cellchat_male

# object.list_female <- list(FB = cellchat_behave,
#                            FC = cellchat_hypo_female_control)
# cellchat_female <- mergeCellChat(object.list_female, add.names = names(object.list_female))
# cellchat_female
# 
# object.list_control <- list(MC = cellchat_hypo_male_control, FC = cellchat_hypo_female_control)
# cellchat_control <- mergeCellChat(object.list_control, add.names = names(object.list_control))
# cellchat_control


gc()


object.list_behave <- list(MB = object.list[["MB"]],
                           FB = object.list[["FB"]])
cellchat_behave <- mergeCellChat(object.list_behave, add.names = names(object.list_behave), cell.prefix = F)
cellchat_behave

gg1 <- compareInteractions(cellchat_behave, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_behave, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


idents_ordered <- c(
  "1.1.1-PV1/CAMK2G",
  "1.1.2-PV1/CAMK2G",
  "1.1.3-PV1/CACNB2",
  "1.2-PV1/SHOX2",
  "1.3-PV1/PCSK2",
  "2.1.1-LYPD1/FOXJ3",
  "2.1.2-LYPD1/RGS11",
  "2.2.1-LYPD1/PCSK2",
  "2.2.2-CCK",
  "2.3-PDYN/AgRP",
  "3.1.1-AgRP/NPY",
  "3.1.2-AgRP/NPY",
  "3.1.3-NPY",
  "3.2.1-SF1/BDNF",
  "3.2.2-PNOC/NPY",
  "3.3.1-GABA/TRH",
  "3.3.2-AgRP/NPY",
  "4-ZEB2/ST18",
  "5.1-CHRM4/FOXP4",
  "5.2-CHRM4/DGKZ",
  "6.1-FOXP4/PROX1",
  "6.2-FOXP4/PROX1",
  "6.3-FOXP4/FOSB",
  "7.1-ZFHX4/KCNH4",
  "7.2.1-ZFHX4/FOSB",
  "7.2.2-ZFHX4/ERBB4",
  "8.1-LMX1A",
  "8.2-LMX1A/vglut2",
  "8.3-LMX1A/vgat/vglut",
  "8.4-LMX1A/nos1",
  "9.1-Oligo/APOE",
  "9.1.1-Oligo/NFASC",
  "9.1.2-MF-Oligo-MPZ",
  "9.2-Oligo/CAMK2G",
  "9.3-Oligo-pnOls-PTPRD",
  "10.1-PNN",
  "10.2-PNN",
  "11.1.1-IRX1/CALN1",
  "11.1.2-IRX1/CALN1",
  "11.2-IRX1/PLA2G4A",
  "12.1-Magno-OXT",
  "12.2.1-Magno-AVP",
  "12.2.2-Magno-AVP",
  "13.1.1-Astro-Proto",
  "13.1.2-Astro-Trans",
  "13.2.1-a-Tany",
  "13.2.2-b-Tany",
  "13.3-Pituicyte-Pit",
  "14.1-MG",
  "14.2-MG",
  "15.1-CALB2/TSHZ2",
  "15.2-NXPH1/BCL11B",
  "16.1-SF1/PACAP",
  "16.2-NPY/CBLN1",
  "17.1.1-PACAP/ZNF536",
  "17.1.2-PACAP/ZNF536",
  "17.2-PRDM8/EBF1",
  "18-Astro-Fibro",
  "19.1-SATB1/PROX1",
  "19.2-ESR1",
  "19.3-LHX6/ARX",
  "20.1.1-GHRH",
  "20.1.2-FOSB",
  "20.2-KNDy",
  "20.3-FOXP1/ISL1",
  "20.4-LepR",
  "21.1.1-RORA/LEF1",
  "21.1.2-QKI/FAT2",
  "21.1.3-RORA/LEF1",
  "21.2-QKI/KCNC3",
  "21.3-GATA2",
  "22-TIDA",
  "23.1-Gonado-Pit",
  "23.2-Folliculo-Pit",
  "24-Thermo",
  "25-OPC",
  "26-FOXP4/RP3V",
  "27.1-Somato-Pit",
  "27.2.1-ptThyro-Pit",
  "27.2.2-pdThyro-Pit",
  "28-Lacto-Pit",
  "29.1-SHOX2/TCF7L2",
  "29.2-SHOX2/TCF7L2",
  "30-q-hmRG",
  "31.1-CP",
  "31.2-BMEC",
  "31.3-PVF",
  "32-ParvoPre-TRH",
  "33-Cortico-Pit",
  "34-CARTPT/CRHR1",
  "35-MEIS2",
  "36.1-Pit-SOX2/PROP1",
  "36.2-Melano-Pit",
  "37-TH",
  "38.1-NFOL",
  "38.2-NFOL",
  "39-Chol"
)

# Apply to CellChat object (safe ordering)
cellchat_behave@idents <- factor(cellchat_behave@idents, levels = idents_ordered)

all_nodes <- idents_ordered
all_nodes

color_key_child <- c(
  "1.1.1-PV1/CAMK2G"            = "#F2BEE1",
  "1.1.2-PV1/CAMK2G"            = "#E993CB",
  "1.1.3-PV1/CACNB2"            = "#ff9acd",
  "1.2-PV1/SHOX2"               = "#E275BD",
  "1.3-PV1/PCSK2"               = "#ffbfe2",
  "2.1.1-LYPD1/FOXJ3"           = "#DB1F83",
  "2.1.2-LYPD1/RGS11"           = "#E6589F",
  "2.2.1-LYPD1/PCSK2"           = "#EE88BA",
  "2.2.2-CCK"                   = "#F7BFDA",
  "2.3-PDYN/AgRP"               = "#F6A5C5",
  "3.1.1-AgRP/NPY"              = "#BBC7FF",
  "3.1.2-AgRP/NPY"              = "#91A7FF",
  "3.1.3-NPY"                   = "#5D7CFA",
  "3.2.1-SF1/BDNF"              = "#4C6EF5",
  "3.2.2-PNOC/NPY"              = "#364FC7",
  "3.3.1-GABA/TRH"              = "#4087FF",
  "3.3.2-AgRP/NPY"              = "#6CA3FF",
  "4-ZEB2/ST18"                 = "#BA2D22",
  "5.1-CHRM4/FOXP4"             = "#EDA0E4",
  "5.2-CHRM4/DGKZ"              = "#EC83E0",
  "6.1-FOXP4/PROX1"             = "#0d8e84",
  "6.2-FOXP4/PROX1"             = "#66c5bc",
  "6.3-FOXP4/FOSB"              = "#b2dfdb",
  "7.1-ZFHX4/KCNH4"             = "#DA81A8",
  "7.2.1-ZFHX4/FOSB"            = "#C94683",
  "7.2.2-ZFHX4/ERBB4"           = "#B80464",
  "8.1-LMX1A"                   = "#FABFB7",
  "8.2-LMX1A/vglut2"            = "salmon4",
  "8.3-LMX1A/vgat/vglut"        = "salmon1",
  "8.4-LMX1A/nos1"              = "salmon3",
  "9.1.1-Oligo/NFASC"           = "#87ab69",
  "9.1.2-MF-Oligo-MPZ"          = "#a3c585",
  "9.1-Oligo/APOE"              = "#4b6043",
  "9.2-Oligo/CAMK2G"            = "#658354",
  "9.3-Oligo-pnOls-PTPRD"       = "#c7ddb5",
  "10.1-PNN"                    = "#ED823E",
  "10.2-PNN"                    = "#F19C7C",
  "11.1.1-IRX1/CALN1"           = "#8b6c5c",
  "11.1.2-IRX1/CALN1"           = "#a08679",
  "11.2-IRX1/PLA2G4A"           = "#bca89f",
  "12.1-Magno-OXT"              = "#F51B1B",
  "12.2.1-Magno-AVP"            = "#F76060",
  "12.2.2-Magno-AVP"            = "#F76060",
  "13.1.1-Astro-Proto"          = "lightgoldenrod3",
  "13.1.2-Astro-Trans"          = "lightgoldenrod1",
  "13.2.1-a-Tany"               = "goldenrod3",
  "13.2.2-b-Tany"               = "#F2D67F",
  "13.3-Pituicyte-Pit"          = "#F2D61F",
  "14.1-MG"                     = "#20b3c4",
  "14.2-MG"                     = "cadetblue2",
  "15.1-CALB2/TSHZ2"            = "#E4928E",
  "15.2-NXPH1/BCL11B"           = "#FF9F9F",
  "16.1-SF1/PACAP"              = "#4444C4",
  "16.2-NPY/CBLN1"              = "#000080",
  "17.1.1-PACAP/ZNF536"         = "#B3544D",
  "17.1.2-PACAP/ZNF536"         = "#F7574A",
  "17.2-PRDM8/EBF1"             = "#E37E76",
  "18-Astro-Fibro"              = "#F2D67F",
  "19.1-SATB1/PROX1"            = "#ade5e2",
  "19.2-ESR1"                   = "#1763AB",
  "19.3-LHX6/ARX"               = "#1C7ED5",
  "20.1.1-GHRH"                 = "#A4D8FF",
  "20.1.2-FOSB"                 = "#339AF0",
  "20.2-KNDy"                   = "#1871C2",
  "20.3-FOXP1/ISL1"             = "#238BE6",
  "20.4-LepR"                   = "#74C0FC",
  "21.1.1-RORA/LEF1"            = "#C3DAFF",
  "21.1.2-QKI/FAT2"             = "#92b4fe",
  "21.1.3-RORA/LEF1"            = "#79B8F4",
  "21.2-QKI/KCNC3"              = "#A7C8FF",
  "21.3-GATA2"                  = "#729efd",
  "22-TIDA"                     = "steelblue3",
  "23.1-Gonado-Pit"             = "#9f87ca",
  "23.2-Folliculo-Pit"          = "#d9b9e2",
  "24-Thermo"                   = "lightsteelblue3",
  "25-OPC"                      = "#aedcae",
  "26-FOXP4/RP3V"               = "indianred",
  "27.1-Somato-Pit"             = "#8953a9",
  "27.2.1-ptThyro-Pit"          = "#bba8d9",
  "27.2.2-pdThyro-Pit"          = "#734f96",
  "28-Lacto-Pit"                = "#bb8fce",
  "29.1-SHOX2/TCF7L2"           = "#EB4F16",
  "29.2-SHOX2/TCF7L2"           = "#EB6E40",
  "30-q-hmRG"                   = "goldenrod1",
  "31.1-CP"                     = "#DD8C00",
  "31.2-BMEC"                   = "#EB9F62",
  "31.3-PVF"                    = "#D49200",
  "32-ParvoPre-TRH"             = "#E66055",
  "33-Cortico-Pit"              = "#7c6192",
  "34-CARTPT/CRHR1"             = "#D631C7",
  "35-MEIS2"                    = "#7865EF",
  "36.1-Pit-SOX2/PROP1"         = "#854abf",
  "36.2-Melano-Pit"             = "#c098f2",
  "37-TH"                       = "#4125F0",
  "38.1-NFOL"                   = "#ACCA82",
  "38.2-NFOL"                   = "#75975e",
  "39-Chol"                     = "#816FF0"
)

color_key_child_df <- data.frame(
  celltype = names(color_key_child),
  color    = unname(color_key_child),
  stringsAsFactors = FALSE
)

saveRDS(cellchat_behave, "/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/hypo_pit/cellchat/Control_Male/child_clusters/HypoPit_cellchatObject_merged_MF_SpawningGroup_L9_Child_level.rds")

groupSize <- table(factor(as.character(cellchat_behave@meta$L9_Child_CellType), levels = all_nodes))
groupSize <- as.numeric(groupSize)
groupSize

color.use.plot <- color_key_child[all_nodes]

pit_celltypes <- c(
  "13.3-Pituicyte-Pit",
  "23.1-Gonado-Pit",
  "23.2-Folliculo-Pit",
  "27.1-Somato-Pit",
  "27.2.1-ptThyro-Pit",
  "27.2.2-pdThyro-Pit",
  "28-Lacto-Pit",
  "33-Cortico-Pit",
  "36.1-Pit-SOX2/PROP1",
  "36.2-Melano-Pit"
)

outdir <- file.path("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/hypo_pit/cellchat/Spawning_Group")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

pdf(file.path(outdir, "SpawningGroup_netCellCellInteractions_L9_Child.pdf"), width = 16, height = 18)

par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction_mod(cellchat_behave, 
                          weight.scale = T, 
                          vertex.weight = groupSize,
                          vertex.size.max	= 5,
                          top = 0.05, 
                          # targets.use = pit_celltypes,
                          color.use = color.use.plot,
                          color.edge = c("royalblue1", "#B80464"),
                          angle_text = T,
                          label.offset.dist = 0.76,
                          label.edge = F
)
netVisual_diffInteraction_mod(cellchat_behave, 
                          weight.scale = T, 
                          vertex.weight = groupSize, 
                          vertex.size.max	= 5,
                          measure = "weight", 
                          top = 0.05,
                          # targets.use = pit_celltypes,
                          color.use = color.use.plot,
                          color.edge = c("royalblue1", "#B80464"),
                          angle_text = T,
                          label.offset.dist = 0.76,
                          label.edge = F)

par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction_mod(cellchat_behave, 
                              weight.scale = T, 
                              vertex.weight = groupSize,
                              vertex.size.max	= 5,
                              top = 0.01, 
                              targets.use = pit_celltypes,
                              color.use = color.use.plot,
                              color.edge = c("royalblue1", "#B80464"),
                              angle_text = T,
                              label.offset.dist = 0.76,
                              label.edge = F
)
netVisual_diffInteraction_mod(cellchat_behave, 
                              weight.scale = T, 
                              vertex.weight = groupSize, 
                              vertex.size.max	= 5,
                              measure = "weight", 
                              top = 0.01,
                              targets.use = pit_celltypes,
                              color.use = color.use.plot,
                              color.edge = c("royalblue1", "#B80464"),
                              angle_text = T,
                              label.offset.dist = 0.76,
                              label.edge = F)

gg1 <- netVisual_heatmap(
  cellchat_behave,
  color.use = color.use.plot,
  slot.name = 'net',
  measure = 'count',
  color.heatmap = c("#1763AB", "#B80464"),
  remove.isolate	= TRUE
)
print(gg1)

gg2 <- netVisual_heatmap(cellchat_behave,
                         color.use = color.use.plot,
                         slot.name = 'net',
                         color.heatmap = c("#1763AB", "#B80464"),
                         remove.isolate	= TRUE,
                         measure = "weight"
                         )
print(gg2)
print(gg1 + gg2)   # In the colorbar, red (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.

library(ComplexHeatmap)
library(circlize)

gg1 <- rankNet(cellchat_behave, 
               mode = "comparison", 
               measure = "weight", 
               sources.use = NULL, 
               targets.use = NULL, 
               stacked = T, 
               do.stat = TRUE,
               color.use = c("royalblue2", "hotpink3")
               )
gg1

gg2 <- rankNet(cellchat_behave, 
               mode = "comparison", 
               measure = "weight", 
               sources.use = NULL, 
               targets.use = NULL, 
               stacked = F, 
               do.stat = TRUE, 
               color.use = c("royalblue2", "hotpink3"))
gg2

gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:length(object.list)) {
  g <- netVisual_circle_mod(object.list[[i]]@net$count,
                            weight.scale = T,
                            label.edge= F,
                            vertex.label.cex = 0.7,
                            remove.isolate = TRUE,
                            edge.weight.max = weight.max[3],
                            edge.width.max = 12,
                            top = 0.075,
                            angle_text = TRUE,
                            arrow.size = 0.05,
                            title.name = paste0("Number of interactions - Top 10% : ", names(object.list)[i]))}


  print(g)
  
  par(mfrow = c(1,1), xpd=TRUE)
  for (i in 1:length(object.list)) {
    g <- netVisual_circle_mod(object.list[[i]]@net$weight,
                              weight.scale = T,
                              label.edge= F,
                              vertex.label.cex = 0.7,
                              remove.isolate = TRUE,
                              edge.weight.max = weight.max[3],
                              edge.width.max = 12,
                              top = 0.075,
                              angle_text = TRUE,
                              arrow.size = 0.05,
                              title.name = paste0("Interaction Weights/Strength - Top 5% : ", names(object.list)[i]))
    print(g)
    
  }


  pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[1], width = 35, height = 28)
  draw(ht1, ht_gap = unit(0.5, "cm"))
  
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[2], width = 35, height = 28)
  draw(ht2, ht_gap = unit(0.5, "cm"))
  
  ht3 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[3], width = 35, height = 28)
  draw(ht3, ht_gap = unit(0.5, "cm"))
  
  ht4 = netAnalysis_signalingRole_heatmap(object.list[[4]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[4], width = 35, height = 28)
  draw(ht4, ht_gap = unit(0.5, "cm"))
  
  dev.off()
  
  pdf(file.path(outdir, "BehaveGroup_OutgoingSignaling_Comparison_HeatMaps_L9_Child.pdf"), width = 31, height = 15)
  
  pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 35, height = 28)
  draw(ht1, ht_gap = unit(0.5, "cm"))
  
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2], width = 35, height = 28)
  draw(ht2, ht_gap = unit(0.5, "cm"))
  
  ht3 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[3], width = 35, height = 28)
  draw(ht3, ht_gap = unit(0.5, "cm"))
  
  ht4 = netAnalysis_signalingRole_heatmap(object.list[[4]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[4], width = 35, height = 28)
  draw(ht4, ht_gap = unit(0.5, "cm"))
  
  num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(object.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
  }
  patchwork::wrap_plots(plots = gg)
  
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  draw(ht3 + ht1, ht_gap = unit(0.5, "cm"))
  draw(ht4 + ht2, ht_gap = unit(0.5, "cm"))
  draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))
  
  dev.off()
  
  pdf(file.path(outdir, "FemaleGroup_IncomingSignaling_Comparison_HeatMaps_L9_Child.pdf"), width = 31, height = 15)
  
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  draw(ht3 + ht1, ht_gap = unit(0.5, "cm"))
  draw(ht4 + ht2, ht_gap = unit(0.5, "cm"))
  draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))
  
  dev.off()
  
  




p_list <- lapply(idents_ordered, function(ct) {
  tryCatch(
    netAnalysis_signalingChanges_scatter(cellchat_behave, idents.use = ct),
    error = function(e) NULL
  )
})
p_list <- Filter(Negate(is.null), p_list)


pdf(file.path(outdir, "SpawningGroup_PathwayComparisonAcrossSexes_perCellType.pdf"),
    width = 25, height = 80)

a <- wrap_plots(p_list, ncol = 4)  # 10 columns will be unreadable
print(a)

dev.off()


# netVisual_bubble(cellchat_behave, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)

pdf(file.path(outdir, paste0("SpawningGroup_PathwayComparisonAcrossSexes_DotPlot", ".pdf")),
    width = 12, height = 18)
gg1 <- netVisual_bubble(cellchat_behave, 
                        sources.use = pit_celltypes,
                        # targets.use = pit_celltypes,  
                        comparison = c(1, 2), 
                        max.dataset = 2, 
                        title.name = "Increased signaling in Spawning Female", 
                        angle.x = 45, 
                        remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat_behave, 
                        sources.use = pit_celltypes, 
                        # targets.use = pit_celltypes,  
                        comparison = c(1, 2), 
                        max.dataset = 1, 
                        title.name = "Decreased signaling in Spawning Female", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
print(gg1 + gg2)
dev.off()


pathways.show <- unique(intersect(cellchat_behave@netP[["MB"]][["pathways"]], cellchat_behave@netP[["FB"]][["pathways"]]))
pathways.show <- pathways.show[1:10]

pdf(file.path(outdir, paste0("SpawningGroup_PathwayComparisonAcrossSexes_HeatMap", ".pdf")),
    width = 30, height = 15)

for (pw in pathways.show) {
  
  ht <- list()
  
  object.list2 <- list(object.list[["MB"]], object.list[["FB"]])
  names(object.list2) <- c("MB", "FB")
  for (i in seq_along(object.list2)) {
    ht[[i]] <- netVisual_heatmap(
      object.list2[[i]],
      signaling = pw,
      color.heatmap = "Reds",
      color.use = color.use.plot,
      title.name = paste0(pw, " signaling - ", names(object.list2)[i])
    )
  }
  
  # combine into one figure (one page) — works for ComplexHeatmap objects
  ht_combined <- Reduce(`+`, ht)
  print(ht_combined)
  ComplexHeatmap::draw(ht_combined, merge_legends = TRUE)
}

dev.off()




cellchat_behave <- computeNetSimilarityPairwise(cellchat_behave, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat_behave <- netEmbedding(cellchat_behave, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat_behave <- netClustering(cellchat_behave, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_behave, type = "functional", label.size = 3.5)

# rankSimilarity(cellchat_behave, type = "functional")

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MB"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to cellchat_behave version < v2, cellchat_behave v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat_behave@var.features$LS.merged and cellchat_behave@var.features$LS.merged.info. 

cellchat_behave <- identifyOverExpressedGenes(cellchat_behave,  pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged cellchat_behave object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat_behave, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat_behave, net = net, datasets = "MB",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat_behave, net = net, datasets = "MB",ligand.logFC = -0.05, receptor.logFC = NULL)
netVisual_embeddingPairwiseZoomIn(cellchat_behave, type = "structural", nCol = 2)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat_behave)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_behave)

df <- findEnrichedSignaling(object.list_female[[2]], features = c("CCL19", "CXCL12"), idents = c("Inflam. FIB", "COL11A1+ FIB"), pattern ="outgoing")

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_behave, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list_female)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat_behave, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list_female)[2]))
#> Comparing communications on a merged object
gg1 + gg2


computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)