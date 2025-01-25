#Testicular DEGs
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
define_testicular_degs<- readRDS('Functions/define_testicular_degs')

}

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

together_data <- data.frame()
for (i in c(0:31)) {
  print(i)
    data <- read.csv( 
              paste0('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/012325 Testicular DEGs/Test_results_cluster_', i, '.csv'))
       data <- define_testicular_degs(data)
    assign(paste0('Test_degs_cluster_', i), data, envir = .GlobalEnv)
    together_data <- rbind(together_data, data)

}

together_data_summed <- together_data%>%
    subset(!is.na(class))%>%
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

### lets plot some of them to make sure they arent gross
plot_gene_behavior <- function(object=obj, 
                               gene,
                               cluster){
  
  htr4_31 <- mean_expression_cluster_data(obj,
                                        gene,
                                        cluster)

behaviors <- data.frame(Percent_Testicular = obj@meta.data$Percent_Testicular,
                        individual = obj@meta.data$individual)%>%
  group_by(individual)%>%
  summarize(Percent_Testicular = mean(Percent_Testicular))

hrt4_31_test <- htr4_31%>%
  full_join(behaviors, by = 'individual')

return(ggplot(hrt4_31_test, aes(x = Percent_Testicular, y = mean, color = Sex))+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(title = paste0(gene,'_cluster_',cluster), y = paste0(gene)))
}

together_data$cluster[together_data$gene=='arg2' & together_data$issignif=='*']

plot_gene_behavior(obj,
                   'arg2',
                   29)

together_data$cluster[together_data$gene=='LOC111579916' & together_data$issignif=='*']

plot_gene_behavior(obj,
                   'LOC111579916',
                   5)

pos_cor_degs  <- together_data$gene[together_data$class!='Interaction' & together_data$issignif=='*'& !is.na(together_data$issignif)]
pos_cor_go <- clown_go(pos_cor_degs)
dotplot(pos_cor_go)

int_degs  <- together_data$gene[together_data$class=='Interaction' & together_data$issignif=='*'& !is.na(together_data$issignif)]
int_go <- clown_go(int_degs)
dotplot(int_go)

int_go$geneID

together_data$cluster[together_data$gene=='fmnl2a' & together_data$issignif=='*']

plot_gene_behavior(obj,
                   'fmnl2a',
                   5)

for(i in all_degs){
  cluster <- together_data$cluster[together_data$gene==i & together_data$issignif=='*'& !is.na(together_data$issignif)]

  print(plot_gene_behavior(obj,
                   i,
                   cluster))

}

cor_1 <- c('spaca6',
           'insm1b',
           'casp8ap2',
           'odr4',
           'tcta',
           'znf532')

cor_1_go <- clown_go(cor_1)

cor_2 <- all_degs[!(all_degs %in% cor_1)]
cor_2_go <- clown_go(cor_2)
dotplot(cor_2_go)#same thing 

## alright I dont rly think theres a lot here




