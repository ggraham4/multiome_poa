
#Cell type annotations 122224
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
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)

}
multiome_object <-readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")


#start aocellaris biomart
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
ocellaris_mus_translator('ND2')

#2) Based on hypomap poa marker genes
#read in hypomap genes
hypomap_allgenes <- read_excel('Cell Type Annotations/HypoMap/hypomap poa mapping.xlsx',
                               6)
#convert to ocellaris
hypomap_allgenes$ocellaris <- lapply(X=hypomap_allgenes$gene,FUN=ocellaris_mus_translator)

hypomap_allgenes <- hypomap_allgenes[hypomap_allgenes$p_val_adj<0.001 &
                                       hypomap_allgenes$specificity>5,]

#read in regions
hypomap_regions  <- read_excel('Cell Type Annotations/HypoMap/hypomap poa mapping.xlsx',
                               4)

### lpoa ####
hypomap_lPOA <- hypomap_regions$cluster_name[hypomap_regions$Region_summarized=='Lateral preoptic area']
hypomap_lPOA<- hypomap_lPOA[!is.na(hypomap_lPOA)]

lpoa_genes <- hypomap_allgenes$ocellaris[hypomap_allgenes$cluster_name %in%hypomap_lPOA]
lpoa_genes_final <- unlist(lpoa_genes)
lpoa_genes_final <- unique(lpoa_genes_final[lpoa_genes_final != ''])
module <- list()
module[[1]] <- c(lpoa_genes_final)

FeaturePlot(multiome_object,lpoa_genes_final[[1]]) 

multiome_object <- AddModuleScore(
  object = multiome_object,
  features = module,
  name = 'lPOA',
  
)
head(x = multiome_object[])

DotPlot(multiome_object, features = 'lPOA1')+
  coord_flip()


#### OK THIS NEEDS TO BE GENERALIZED #####
POA_regions <- c('Lateral preoptic area',
                 'Medial preoptic area',
                 'Periventricular hypothalamic nucleus, intermediate part',
                'Periventricular hypothalamic nucleus, posterior part',
                 'Paraventricular hypothalamic nucleus'
)

module <- list()
for(i in unique(hypomap_regions$Region_summarized)){
  b = 1
  #get list of names of the clusters in that region
  cluster_names <- hypomap_regions$cluster_name[hypomap_regions$Region_summarized==i]
  cluster_names<- cluster_names[!is.na(cluster_names)]
  
  
  region_genes <- hypomap_allgenes$ocellaris[hypomap_allgenes$cluster_name %in%cluster_names]
  region_genes_final <- unlist(region_genes)
  region_genes_final <- unique(region_genes_final[region_genes_final != ''])

  module[[paste0(i)]] <- c(region_genes_final) #this keeps overwriting
  
  b = b+1
  }
  

multiome_object <- AddModuleScore(
  object = multiome_object,
  features = module,
  name = 'POA_modules',
  
)
`%notin%` <- Negate(`%in%`)


 DotPlot(subset(multiome_object,harmony.wnn_res0.4_clusters %notin% c(2,4,29,30, 15, 22, 14, 26,18)), 
        features = c(paste0(rep('POA_modules',15),1:15)),
        dot.min = .5,
        col.min =  .5*2.5)+
  coord_flip()+
  scale_x_discrete(labels = names(module))+
   scale_color_gradientn(colors = c('white','red','blue'))
 
 
 #find top score for each cluster and each module

 DotPlot(subset(multiome_object,harmony.wnn_res0.4_clusters %notin% c(2,4,29,30, 15, 22, 14, 26,18)), 
         features = c(paste0(rep('POA_modules',15),1:15)),
         dot.min = 0,
         col.min =  0
 )+
   coord_flip()+
   scale_x_discrete(labels = names(module))+
   scale_color_gradientn(colors = c('white','red','blue'))
 

