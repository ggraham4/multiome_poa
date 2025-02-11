{
  library(parallel)
  library(clusterProfiler)
  library(blme)
  library(Seurat)
  library(tidyverse)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(SeuratObject)
  library(Signac)
  library(glmGamPoi)
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)
  library(emmeans)
  library(CytoTRACE)
  library(ggrepel)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(Polychrome)
  library(scCustomize)

P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
names(P40) <- NULL

mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
define_degs_prop<- readRDS('Functions/define_degs_prop.rds')
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
prop_deg_function.rds<- readRDS('Functions/DEG_functions/prop_deg_function.rds')
define_behavior_degs<- readRDS('Functions/define_behavior_degs')
clown_go<- readRDS('Functions/clown_go')
define_degs<- readRDS('Functions/define_degs')

library(hdWGCNA)
library(WGCNA)

}

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')
radial_glia <- subset(obj, harmony.wnn_res0.4_clusters==2)
radial_glia <- FindSubCluster(radial_glia, cluster = 2, subcluster.name = 'sub', graph.name = 'harmony.wsnn')

#https://www.nature.com/articles/s41467-019-14198-8
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")
att <-listAttributes(ensembl)

#use ensemble_transcript id as a bridge between entrezgene_accession and mmusuculus homolog
biomart_gene_name <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c('ensembl_transcript_id',
                   'entrezgene_accession'))
biomart_mmusculus <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c('ensembl_transcript_id',
                   'mmusculus_homolog_associated_gene_name'))

#merge them
joined <- right_join(biomart_gene_name, biomart_mmusculus, by ='ensembl_transcript_id')[,2:3]
#remove duplicates
joined2 <- joined%>%dplyr::distinct(
                            .keep_all = T)

#rename
translated <- joined2

#create translator function
ocellaris_mus_translator <- function(mus_gene){
  mus_gene_lower <- tolower(mus_gene)
  ocellaris_gene <- translated$entrezgene_accession[
    translated$mmusculus_homolog_associated_gene_name == mus_gene]
  return(ocellaris_gene)
  }
ocellaris_mus_translator('Polr3a')



markers <- read_excel("Reference/mmusculus_astrocyte_markers_Batiuk_etal.xlsx",2)
genes <- markers$`Marker gene`[markers$`p-value (Wilcox test)`<0.05 & markers$`Fraction in the population of interest`>0.8]
genes_translated <- lapply(genes, ocellaris_mus_translator)
genes_translated <- unlist(genes_translated)
genes_translated<- unique(genes_translated)

DotPlot(radial_glia, group.by = 'sub',
       features = c('gfap',genes_translated))+
  coord_flip()

DotPlot(radial_glia, group.by = 'sub',
       features = 'LOC111577263')+
  coord_flip()

DotPlot(radial_glia, group.by = 'sub',
       features = c('LOC111580654',
                    'nr4a1',
                    'fbxo32',
                    'mideasa',
                    'mideasb',
                    'stox1',
                    'LOC111562412',
                    'alx1'
         ))+
  coord_flip()


