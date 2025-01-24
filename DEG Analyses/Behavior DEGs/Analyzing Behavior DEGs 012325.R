#behavior_deg_analysis
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

}

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

for (i in c(0:31)) {
    data <- read.csv( 
              paste0('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/012025 Behavior DEGs/behavior_results_cluster_', i, '.csv'))
       
    assign(paste0('beh_degs_cluster_', i), data, envir = .GlobalEnv)

}

together_data <-data.frame()
for(i in 0:31){
  
  data <- get(paste0('beh_degs_cluster_', i))
      data <- define_behavior_degs(data)
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

clust_7_genes <- together_data$gene[together_data$cluster ==7 & together_data$issignif=='*'& !is.na(together_data$issignif)]
clust_7_go <- clown_go(clust_7_genes)
dotplot(clust_7_go)

#what about all of them except 30
all_degs <- together_data$gene[together_data$cluster !=30 & together_data$cluster!=15 & together_data$issignif=='*'& !is.na(together_data$issignif)]
all_go <- clown_go(all_degs)
dotplot(all_go)

#nothin

clust_26_genes <- together_data$gene[together_data$cluster ==26 & together_data$issignif=='*'& !is.na(together_data$issignif)]
clust_26_go <- clown_go(clust_26_genes)
dotplot(clust_26_go)

ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")


biomart_basic <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'entrezgene_description'))

named <- biomart_basic[biomart_basic$entrezgene_accession %in%all_degs,]

#which one expresss serotonin receptor htr4
together_data$cluster[together_data$gene=='htr4' & together_data$issignif=='*']
#cluster 31 huh

htr4_31 <- mean_expression_cluster_data(multiome_obj,
                                        'htr4',
                                        31)

behaviors <- data.frame(behaviors = multiome_obj@meta.data$Behaviors_Day_2,
                        individual = multiome_obj@meta.data$individual)%>%
  group_by(individual)%>%
  summarize(behaviors = mean(behaviors))

hrt4_31_behaviors <- htr4_31%>%
  full_join(behaviors, by = 'individual')

ggplot(hrt4_31_behaviors, aes(x = behaviors, y = mean, color = Sex))+
  geom_point()+
  geom_smooth(method = 'lm')

plot_gene_behavior <- function(object, 
                               gene,
                               cluster){
  
  htr4_31 <- mean_expression_cluster_data(multiome_obj,
                                        gene,
                                        cluster)

behaviors <- data.frame(behaviors = multiome_obj@meta.data$Behaviors_Day_2,
                        individual = multiome_obj@meta.data$individual)%>%
  group_by(individual)%>%
  summarize(behaviors = mean(behaviors))

hrt4_31_behaviors <- htr4_31%>%
  full_join(behaviors, by = 'individual')

print(ggplot(hrt4_31_behaviors, aes(x = behaviors, y = mean, color = Sex))+
  geom_point()+
  geom_smooth(method = 'lm'))

  
}

together_data$cluster[together_data$gene=='hdac3' & together_data$issignif=='*']

plot_gene_behavior(multiome_obj,
                   'hdac3',
                   31)

together_data$cluster[together_data$gene=='saga' & together_data$issignif=='*']
plot_gene_behavior(multiome_obj,
                   'saga',
                   31)


together_data$cluster[together_data$gene=='tph1a' & together_data$issignif=='*']
plot_gene_behavior(multiome_obj,
                   'tph1a',
                   31)

together_data$cluster[together_data$gene=='neurod1' & together_data$issignif=='*']
plot_gene_behavior(multiome_obj,
                   'neurod1',
                   31)

clust_22_genes <- together_data$gene[together_data$cluster ==22 & together_data$issignif=='*'& !is.na(together_data$issignif)]
clust_22_go <- clown_go(clust_22_genes)
dotplot(clust_22_go)

clust_22_go$geneID[clust_22_go$Description == 'negative regulation of DNA-templated transcription'
]









