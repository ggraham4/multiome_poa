#IEG analysis second attempt
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
  library(CytoTRACE)
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
  library(glmnet)  
  `%notin%` <- Negate(`%in%`)
}

obj <-  readRDS("~/Desktop/snRNA-seq R Files 122524/RNA Object.rds")

obj$ieg <- ifelse(obj@assays$RNA$data['LOC111583367',] >0, 'ieg', NA)
obj$ieg <- ifelse(obj@assays$RNA$data['egr1',] >0, 'ieg', obj$ieg)
obj$ieg <- ifelse(obj@assays$RNA$data['npas4a',] >0, 'ieg', obj$ieg)
obj$ieg <- ifelse(is.na(obj$ieg), 'no_ieg', obj$ieg)

DimPlot(obj, split.by = 'ieg')

### filter out non-neuronal
neuronal_only <- subset(obj, harmony.wnn_res0.4_clusters %notin%
                          c(
                            2, #radial glia
                            4, #olig
                            14, #mg
                            18, #OPC
                            22, #ependymal
                            26, #leuko,
                            29, #fibro
                            15, #exclude
                            30 #OB
                          )
)
DimPlot(neuronal_only, split.by = 'ieg')
#### ok a few questions I can now ask

#1) Which clusters are enriched for proportion of cells expressing ieg
#>2) How do the sexes differ in these measures
#>3) Do the findmarkers things zack discusses in his methods
#>4) 

ieg_data_cluster_level <- table(neuronal_only@meta.data$ieg, neuronal_only@meta.data$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  pivot_wider(names_from = 'Var1', values_from = 'Freq')%>%
  mutate(prop_pos = ieg/no_ieg)

###CLUSTER 27 IS THE HIGHEST WAOW
### ITS NOT EVEN CLOSE

ggplot(ieg_data_cluster_level, aes(x = as.factor(Var2), y = prop_pos))+
  geom_bar(stat = 'identity')+
  geom_hline(yintercept = 0.5, linetype = 2)+
    geom_hline(yintercept =1, linetype = 2)+
  labs(x = 'Cluster', y = 'IEG + / IEG -')+
  theme_minimal()

### Lets do it with sex now

ieg_data_individual_level <- neuronal_only@meta.data%>%
  dplyr::select(c(individual, Status, harmony.wnn_res0.4_clusters, ieg))%>%
  group_by(Status, individual, harmony.wnn_res0.4_clusters)%>%
  summarize(ieg_pos = sum(ieg =='ieg'),
            ieg_neg = sum(ieg == 'no_ieg')
  )%>%
  mutate(prop_pos = ieg_pos/ieg_neg)%>%
  subset(Status %in% c('M','D','F'))%>%
  na.omit()%>%
  subset(prop_pos != Inf)

ieg_data_individual_level$Status <- factor(ieg_data_individual_level$Status, levels = c('M','D','F'))

ggplot(ieg_data_individual_level, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = prop_pos, shape = Status, color = Status))+
    geom_hline(yintercept = 0.5, linetype = 2)+
    geom_hline(yintercept =1, linetype = 2)+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'IEG + / IEG -')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))+
  scale_y_continuous(breaks = 0:7)


