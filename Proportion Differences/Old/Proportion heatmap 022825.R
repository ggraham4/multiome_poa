##making a proportion heatmap

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
  library(CellChat)
  
  
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
  `%notin%` <- Negate(`%in%`)
}

obj <- readRDS('C:/Users/Gabe/Desktop/RNA Object.rds')

cells_total <- obj@meta.data%>%
  group_by(individual)%>%
  summarize(cells = n())

cluster_populations <- obj@meta.data%>%
  group_by(individual, harmony.wnn_res0.4_clusters)%>%
  summarize(cells_in_cluster = n())

joined_pops <- cells_total%>%
  right_join(cluster_populations, by= 'individual')%>%
  mutate(prop_in_cluster = cells_in_cluster/cells)

corr.data = data.frame()
for(i in unique(joined_pops$harmony.wnn_res0.4_clusters)){
  print(i)
  cluster_1 = i
  temp_joined_data_1 <- subset(joined_pops, harmony.wnn_res0.4_clusters==cluster_1)
  cluster_1_props = temp_joined_data_1[,c(1, 5)]
  colnames(cluster_1_props) = c('individual', 'props_cluster_1')
  
  for(j in unique(joined_pops$harmony.wnn_res0.4_clusters)){
    print(j)
    cluster_2 = j
    temp_joined_data_2 <- subset(joined_pops, harmony.wnn_res0.4_clusters==cluster_2)
    cluster_2_props = temp_joined_data_2[,c(1, 5)]
    
    colnames(cluster_2_props) = c('individual', 'props_cluster_2')
    
    test_data <-cluster_1_props%>%
      right_join(cluster_2_props, by = 'individual')%>%
      na.omit()
    
    correlation <- cor.test(x = test_data$props_cluster_1, y = test_data$props_cluster_2)
    
    new_data <- data.frame('cluster_1' =cluster_1,
                           'cluster_2' = cluster_2,
                           correlation_coefficient = correlation$estimate,
                           correlation_p.value = correlation$p.value
      )
    corr.data <- rbind(corr.data, new_data)
    } 
  
}

corr.data$correlation_q.value<- p.adjust(corr.data$correlation_p.value, 'fdr', nrow(corr.data))

corr.data$correlation_coefficient = as.numeric(corr.data$correlation_coefficient)

corr_matrix <- corr.data%>%
  dplyr::select(c('cluster_1','cluster_2', 'correlation_coefficient'))%>%
  pivot_wider(names_from = 'cluster_2', values_from = 'correlation_coefficient')

corr_matrix <- sapply( corr_matrix, as.numeric )
rownames(corr_matrix) = as.factor(corr_matrix[,1])
corr_matrix = corr_matrix[,-1]
colnames(corr_matrix) = as.character(colnames(corr_matrix))

heatmap(corr_matrix[,-1])

library(viridis)
corr.data$issignif <- ifelse(corr.data$correlation_q.value <0.05, '*', NA)

ggplot(corr.data, aes(x= as.factor(as.numeric(cluster_1)), y = as.factor(as.numeric(cluster_2)), fill = correlation_coefficient))+
  geom_tile()+
  scale_fill_gradientn(colors = c('red',"white","blue"))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())+
  geom_text(aes(label = issignif), size =7, vjust = 0.75)
  

