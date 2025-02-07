### I'm cutious about IT (oxt) and AVT (AVP)
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
  library(hdWGCNA)
  library(WGCNA)
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
}

#load in seurat obj
obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

DotPlot(obj, features=c('oxt','avp'), group.by = 'harmony.wnn_res0.4_clusters')+
  coord_flip()

# both expressed in cluster 0

FeaturePlot(obj, 'oxt')

FeaturePlot(obj, 'avp')

#both expressed in that little top thing

#read in proportion data from cluster 0

prop_0 <- read_csv("DEG Outputs/Prop DEGs 011525/prop_degs_cluster_0.csv")

prop_0_defined <- define_degs_prop(prop_0)
#neither differ in proportion

####See of any clusters differ by avp or oxt prop
together_data <- data.frame()
for(i in 0:31){
  data <- suppressMessages(define_degs_prop(
    read_csv(
          paste0(
            "DEG Outputs/Prop DEGs 011525/prop_degs_cluster_",i,".csv")  
    )
    )
  )
  data$cluster =i
together_data <- rbind(together_data,data)
  }

gene_query <- function(gene){
  sig_clusters <- together_data$cluster[together_data$gene == gene & 
                                          together_data$issignif == '*' &
                                          !is.na(together_data$issignif)]
  return(sig_clusters)
  
  }
gene_query('tacr3a')
gene_query('oxt')
gene_query('oxtr')
gene_query('galr1a')
gene_query('galr2a')
gene_query('cckb') #18 - OPCs
gene_query('npy')
gene_query('tac1')
gene_query('tac2')
gene_query('tac3')
gene_query('tac4')
gene_query('tacr1')
gene_query('kiss1')
gene_query('kiss1r')




