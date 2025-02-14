
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
  #library(CytoTRACE)
 # SCTRransform_mean_plot <- readRDS("R/Gabe/SCTRransform_mean_plot.rds")
  #mac.neg.bin <- readRDS(file = 'R/Gabe/mac.neg.bin.rds')
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  #mean_cell <- readRDS('R/Gabe/mean_cell.rds')
  library(openxlsx)
  #clown_go <- readRDS('R/Gabe/clown_go.rds')
  
}
multiome_object <- readRDS('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/RNA Object.rds')


####OK I think I want to take a new approach to 
##the spatial mapping of poa clusters

#1) I will create a mouse to ocellaris gene translator (zack might have code for this)

#create a biomart of ocellaris to mmusculus orthologs

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

#2) Based on hypomap poa marker genes, for each subregion, I will extract the 
#marker genes, convert them to ocellaris, and see whicch cluster has the strongest expression

#read in hypomap genes
hypomap_allgenes <- read_excel('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/spatial mapping cichlid/hypomap poa mapping.xlsx',
                               6)
#convert to ocellaris
hypomap_allgenes$ocellaris <- lapply(X=hypomap_allgenes$gene,FUN=ocellaris_mus_translator)

#read in regions
hypomap_regions  <- read_excel('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/spatial mapping cichlid/hypomap poa mapping.xlsx',
                               4)
### lPOA ###
hypomap_lPOA <- hypomap_regions$cluster_name[hypomap_regions$Region_summarized=='Lateral preoptic area']
hypomap_lPOA<- hypomap_lPOA[!is.na(hypomap_lPOA)]

lpoa_genes <- hypomap_allgenes$ocellaris[hypomap_allgenes$cluster_name %in%hypomap_lPOA]
lpoa_genes_final <- unlist(lpoa_genes)
lpoa_genes_final <- unique(lpoa_genes_final[lpoa_genes_final != ''])


### by mean gene expression ###
expression <-  Seurat::AverageExpression(multiome_object, group.by = 'harmony.wnn_res0.4_clusters')[['RNA']]
colnames(expression)

lpoa_expression <- as.data.frame(expression)
lpoa_expression$gene <- rownames(lpoa_expression)
lpoa_expression <- lpoa_expression[lpoa_expression$gene%in%lpoa_genes_final,]
lpoa_expression_matrix <- as.matrix(lpoa_expression[,1:31])

cluster_means_lpoa <- colMeans(lpoa_expression_matrix)
sort(cluster_means_lpoa) ### 20, 30, 10 are teh most likely

### lets make this a function

likely_region_expression <- function(region){
  #extract region codes
  hypomap_region <- hypomap_regions$cluster_name[hypomap_regions$Region_summarized==region] 
  hypomap_region<- hypomap_region[!is.na(hypomap_region)]
  
  #select codes for that region
  region_genes <- hypomap_allgenes$ocellaris[hypomap_allgenes$cluster_name %in%hypomap_region&
                                               hypomap_allgenes$specificity>5]
  region_genes_final <- unlist(region_genes)
  region_genes_final <- unique(region_genes_final[region_genes_final != ''])
  
  #get average expression for that region
  region_expression <- as.data.frame(expression)
  region_expression$gene <- rownames(region_expression)
  region_expression <- region_expression[region_expression$gene%in%region_genes_final,]
  region_expression_matrix <- as.matrix(region_expression[,1:32])

  cluster_means_region <- colMeans(region_expression_matrix)
# sorted <- sort(cluster_means_region, decreasing = T) 
return(cluster_means_region)
}

likely_region_expression('Lateral preoptic area')
likely_region_expression('Medial preoptic area')
likely_region_expression('Periventricular hypothalamic nucleus, intermediate part')
likely_region_expression('Periventricular hypothalamic nucleus, posterior part')
likely_region_expression('Paraventricular hypothalamic nucleus')

#I want to make this a data frame
data <- data.frame()
for(i in unique(hypomap_regions$Region_summarized)){
  exp <- likely_region_expression(i)
  
  exp$region <-i
  
  data <- rbind(data, exp)
}

write.csv(data, '/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/spatial mapping cichlid/ranking regions based on hypomap expression.csv')
data <- na.omit(data)
                    
data_tall<- pivot_longer(data, cols=-33,
                         names_to = 'cluster')

data_tall$clust_fact <- as.numeric(str_sub(data_tall$cluster,start =2))

ggplot(data_tall, aes(x=clust_fact, y = region, fill = log(value)))+
  scale_x_continuous(breaks = c(0:31))+
  scale_fill_gradientn(colors = c('white','yellow','blue'))+
  geom_tile()

ggplot(subset(data_tall, region != 'Periventricular hypothalamic nucleus, intermediate part'), aes(x=cluster, y = region, fill = log(value)))+
  geom_tile()


heatmap(as.matrix(data[2:16,1:32]))

write.csv(data_tall,'/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/spatial mapping cichlid/mean hypomap marker expression.csv')


#### now I want to do the same with % of nuclei expressing the genes

library(purrr)

proportion_expressing <- function(object=multiome_object,
                         gene,
                         clustering='harmony.wnn_res0.4_clusters',
                         cluster
){
  
  Counts_of_interest <- as.data.frame(RNA_counts[,gene])
  Counts_of_interest[[gene]] <- Counts_of_interest[,1]
  
  temp_data <-Counts_of_interest

    if(nrow(temp_data)==0){next}
    temp_data <- temp_data[, c(gene)]    
    n_nuclei <- length(temp_data)

    n_expressing <- count(temp_data>0)

    proportion = n_expressing/n_nuclei
    data <- data.frame(cluster = cluster,
                       gene = gene,
                       proportion = proportion)
return(data)
}

ocellaris_genes <- unlist(unique(hypomap_allgenes$ocellaris))
  ocelaris_genes <- ocellaris_genes[ocellaris_genes!='']
  
prop_data <- data.frame()
for(i in 0:31){
  RNA_counts <- t(multiome_object@assays$RNA$counts[,multiome_object@meta.data[['harmony.wnn_res0.4_clusters']] == i])
  for(g in ocelaris_genes){
    print(g)
    # Use tryCatch to handle errors in proportion_expressing
    temp <- tryCatch(
      {
        proportion_expressing(gene = g, cluster = i)
      },
      error = function(e) {
        message(paste("Error in gene:", g, "cluster:", i, "-", e$message))
        return(NULL)
      }
    )
    # Only bind to prop_data if temp is not NULL
    if (!is.null(temp)) {
      prop_data <- rbind(prop_data, temp)
    }
  }
}
  write.csv(prop_data, '/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/spatial mapping cichlid/ranking regions based on hypomap prop nuclei expressing.csv')
proportion <- prop_data

likely_region_proportion <- function(region){
  #extract region codes
  hypomap_region <- hypomap_regions$cluster_name[hypomap_regions$Region_summarized==region] 
  hypomap_region<- hypomap_region[!is.na(hypomap_region)]
  
  #select codes for that region
  region_genes <- hypomap_allgenes$ocellaris[hypomap_allgenes$cluster_name %in%hypomap_region&
                                               hypomap_allgenes$specificity>5]
  region_genes_final <- unlist(region_genes)
  region_genes_final <- unique(region_genes_final[region_genes_final != ''])
  
  
  #get average expression for that region
  region_proportion <- proportion[proportion$gene%in%region_genes_final,]
  region_proportion_matrix <- pivot_wider(region_proportion, values_from = 'proportion',
                                          names_from = 'cluster',values_fn = {mean})

  cluster_means_region <- colMeans(as.matrix(region_proportion_matrix[2:32]))
  cluster_means_region<- as.data.frame(cluster_means_region)
  cluster_means_region$cluster <-colnames(region_proportion_matrix[2:32])
  
# sorted <- sort(cluster_means_region, decreasing = T) 
return(cluster_means_region)
}

####why is this not working
data_prop<- data.frame()
for(i in unique(hypomap_regions$Region_summarized)){
  print(i)
  if(is.na(i)){next}
  prop <- likely_region_proportion(i)
  prop <- as.data.frame(prop)
  prop$region <-i
  
  data_prop <- rbind(data_prop, prop)
}

ggplot(data_prop, aes(x=as.numeric(cluster), y = region, fill = sqrt(cluster_means_region)))+
  scale_x_continuous(breaks = c(0:31))+
  scale_fill_gradientn(colors = c('white','yellow','blue'))+
  geom_tile()


write.csv(data_prop,'/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/spatial mapping cichlid/mean hypomap marker proportions.csv')
#3 I will rank clusters based on likelihood of mPOA, aPOA, pPOA, PVN






