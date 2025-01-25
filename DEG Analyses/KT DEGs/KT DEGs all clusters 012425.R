#KTicular DEGs
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
  
  library(Polychrome)
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
define_kt_degs<- readRDS('Functions/define_kt_degs')

}

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

together_data <- data.frame()
for (i in c(0:31)) {
  print(i)
    data <- read.csv( 
              paste0('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/012425 KT DEGs/KT_results_cluster_', i, '.csv'))
       data <- define_kt_degs(data)
    assign(paste0('KT_degs_cluster_', i), data, envir = .GlobalEnv)
    together_data <- rbind(together_data, data)

}

together_data_summed <- together_data%>%
    subset(!is.na(class)&cluster !=30)%>%
  group_by(cluster, class)%>%
  summarize(class_count = n())

ggplot(together_data_summed, aes(x = cluster, y = class_count, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DEGs") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
    scale_x_continuous(breaks = c(0:31))+
  scale_y_continuous()+
  scale_fill_manual(values = P40)

### What are these DEGs ####
all_degs <- together_data$gene[together_data$cluster !=30 & together_data$cluster!=15 & together_data$issignif=='*'& !is.na(together_data$issignif)]
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")
biomart_basic <- getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'entrezgene_description'))

named <- biomart_basic[biomart_basic$entrezgene_accession %in%all_degs,]

###are there any repeated DEGs
length(all_degs)
length(unique(all_degs))
# no

together_data$class2 <- ifelse(together_data$m_d_q.value<0.05, '*','no')
together_data$gene[together_data$class2 == '*']

together_data_summed2 <- together_data%>%
    subset(class2 == '*'&cluster !=30)%>%
  group_by(cluster, class2)%>%
  summarize(class_count = n())

ggplot(together_data_summed2, aes(x = cluster, y = class_count, fill = class2)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DEGs") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
    scale_x_continuous(breaks = c(0:31))+
  scale_y_continuous()+
  scale_fill_manual(values = P40)

names2 <- together_data$gene[together_data$class2=='*' & !is.na(together_data$class2)&together_data$cluster !=30]
named2 <- biomart_basic[biomart_basic$entrezgene_accession %in%names2,]

### what about anova status
together_data$class3 <- ifelse(together_data$status_anova_q.value<0.05,'*','no')
status_anova_genes <- together_data$gene[together_data$class3 =='*'& together_data$cluster !=30]
#ok a lot

together_data_summed3 <- together_data%>%
    subset(class3 == '*'&cluster !=30)%>%
  group_by(cluster, class3)%>%
  summarize(class_count = n()) 






