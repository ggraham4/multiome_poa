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

obj <- readRDS('C:/Users/Gabe/Desktop/RNA Object.rds')

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
  mutate(prop_pos = ieg/(no_ieg+ieg))

###CLUSTER 27 IS THE HIGHEST WAOW

ggplot(ieg_data_cluster_level, aes(x = as.factor(Var2), y = prop_pos))+
  geom_bar(stat = 'identity')+
  geom_hline(yintercept = 0.5, linetype = 2)+
  labs(x = 'Cluster', y = 'IEG + / IEG -')+
  theme_minimal()

### Lets do it with sex now

ieg_data_individual_level <- neuronal_only@meta.data%>%
  dplyr::select(c(individual, Status, harmony.wnn_res0.4_clusters, ieg))%>%
  group_by(Status, individual, harmony.wnn_res0.4_clusters)%>%
  summarize(ieg_pos = sum(ieg =='ieg'),
            ieg_neg = sum(ieg == 'no_ieg')
  )%>%
  mutate(prop_pos = ieg_pos/(ieg_pos+ieg_neg))%>%
  subset(Status %in% c('M','D','F'))%>%
  na.omit()%>%
  subset(prop_pos != Inf)

ieg_data_individual_level$Status <- factor(ieg_data_individual_level$Status, levels = c('M','D','F'))

ggplot(ieg_data_individual_level, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = prop_pos, shape = Status, color = Status))+
  geom_hline(yintercept = 0.5, linetype = 2)+
  geom_hline(yintercept =1, linetype = 2)+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Proportion IEG+')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))+
  scale_y_continuous(breaks = c(0.5,0:7))

ggplot(ieg_data_individual_level, aes(x = Status, y = prop_pos, shape = Status, color = Status))+
  geom_hline(yintercept = 0.5, linetype = 2)+
  geom_hline(yintercept =1, linetype = 2)+
  geom_violin(alpha = 0, outlier.shape = NA)+
  geom_boxplot(size = 1)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Proportion IEG+')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))+
  scale_y_continuous(breaks = c(0.5,0:7))


### FindMarkers
#Find genes differentially expressed between IEG and IEG- clusters
ieg_marker_all <-data.frame()
for(i in unique(neuronal_only$harmony.wnn_res0.4_clusters)){
  print(i)
  temp_obj <- subset(neuronal_only, harmony.wnn_res0.4_clusters==i)
  Idents(temp_obj) = 'ieg'
  ieg_markers <- FindAllMarkers(temp_obj, 
                                assay = 'RNA', 
                                group.by = 'ieg',
                                logfc.threshold = 0,
                                min.pct = 1/204
  ) #204 cells in smallest cluster
  ieg_markers$harmony.wnn_res0.4_clusters =i
  ieg_marker_all <- rbind(ieg_marker_all, ieg_markers)
}
#extract significant ieg like genes
ieg_marker_all_signif <- ieg_marker_all[ieg_marker_all$p_val_adj<0.05 &ieg_marker_all$cluster=='ieg' ,]

#find genes significant in over half of clusters
ieg_like_genes_data <- table(ieg_marker_all_signif$gene,
                             ieg_marker_all_signif$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  group_by(Var1)%>% #cluster
  summarize(n_clusters = sum(Freq))%>%
  subset(n_clusters > length(unique(neuronal_only$harmony.wnn_res0.4_clusters))/2
  ) 

#ieg like genes
ieg_like_genes <- ieg_like_genes_data$Var1%>%droplevels()

neuronal_only$ieg2 <- NA
#Classify cells that express any of these iegs as ieg
for(i in ieg_like_genes){
  neuronal_only$ieg2 <-  ifelse(neuronal_only@assays$RNA$data[i,] >0, 'ieg', neuronal_only$ieg2)
  
}
neuronal_only$ieg2<-  ifelse(is.na(neuronal_only$ieg2), 'non_ieg', neuronal_only$ieg2)


ieg_data_cluster_level2 <- table(neuronal_only@meta.data$ieg2, neuronal_only@meta.data$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  pivot_wider(names_from = 'Var1', values_from = 'Freq')%>%
  mutate(prop_pos = ieg/(non_ieg+ieg))%>%
  na.omit()


ggplot(ieg_data_cluster_level2, aes(x = as.factor(Var2), y = prop_pos))+
  geom_bar(stat = 'identity')+
  geom_hline(yintercept = 0.5, linetype = 2)+
  geom_hline(yintercept =1, linetype = 2)+
  labs(x = 'Cluster', y = 'Prop IEG+')+
  theme_minimal()



ieg_data_individual_level2 <- neuronal_only@meta.data%>%
  dplyr::select(c(individual, Status, harmony.wnn_res0.4_clusters, ieg2))%>%
  group_by(Status, individual, harmony.wnn_res0.4_clusters)%>%
  summarize(ieg_pos = sum(ieg2 =='ieg'),
            ieg_neg = sum(ieg2 == 'non_ieg')
  )%>%
  mutate(prop_pos = ieg_pos/(ieg_neg+ieg_pos))%>%
  subset(Status %in% c('M','D','F'))%>%
  na.omit()%>%
  subset(prop_pos != Inf)

ieg_data_individual_level2$Status <- factor(ieg_data_individual_level2$Status, levels = c('M','D','F'))

ggplot(ieg_data_individual_level2, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = prop_pos, shape = Status, color = Status))+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Proportion IEG+')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))

ggplot(ieg_data_individual_level2, aes(x = Status, y = prop_pos, shape = Status, color = Status))+
  geom_hline(yintercept = 0.5, linetype = 2)+
  geom_hline(yintercept =1, linetype = 2)+
  geom_violin(alpha = 0, outlier.shape = NA)+
  geom_boxplot(size = 1)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Proportion IEG+')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))+
  scale_y_continuous(breaks = c(0.5,0:7))

RNA_data <- neuronal_only@assays$RNA$data
ieg_rows_rna <- RNA_data[which(rownames(RNA_data) %in% ieg_like_genes),]
ieg_scores <- colSums(ieg_rows_rna>0)

neuronal_only$ieg_scores = ieg_scores

ieg_scores <- neuronal_only@meta.data%>%
  group_by(individual, Status, harmony.wnn_res0.4_clusters)%>%
  summarize(mean_ieg_score = mean(ieg_scores))%>%
  subset(Status %in% c('M','D','F'))

ieg_scores$Status <- factor(ieg_scores$Status, levels = c('M','D','F'))


ieg_positive <- neuronal_only@meta.data%>%
  group_by(individual, Status, harmony.wnn_res0.4_clusters)%>%
  summarize(ieg_pos = sum(ieg_scores>0),
            ieg_neg = sum(ieg_scores ==0),
            prop_pos = ieg_pos/ieg_neg)%>%
  subset(Status %in% c('M','D','F'))%>%
  filter(prop_pos != Inf)

ieg_positive$Status <- factor(ieg_positive$Status, levels = c('M','D','F'))
ggplot(ieg_positive, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = ieg_pos/(ieg_pos+ieg_neg), shape = Status, color = Status))+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Prop IEG+')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))

### proportion of IEG+ cells
glmer_out <- list()
for(i in unique(neuronal_only$harmony.wnn_res0.4_clusters)){
  print(i)
  temp_ieg_positive <- subset(ieg_positive, harmony.wnn_res0.4_clusters==i)
  ieg_pos_matrix <- matrix(cbind(temp_ieg_positive$ieg_pos, temp_ieg_positive$ieg_neg), nrow(temp_ieg_positive),2)
  
  glmer_model <- glmer(ieg_pos_matrix~Status + (1|individual), data = temp_ieg_positive, family = binomial('logit'))
  
  glmer_out[[paste0('cluster_',i)]]$model <- glmer_model
  glmer_out[[paste0('cluster_',i)]]$av <- car::Anova(glmer_model, type ='III')
  glmer_out[[paste0('cluster_',i)]]$summary <- summary(glmer_model)
  
}

library(emmeans)
glmer_av_output <- data.frame()
for(i in unique(neuronal_only$harmony.wnn_res0.4_clusters)){
  print(i)
  av <- glmer_out[[paste0('cluster_',i)]]$av%>%as.data.frame()
  status_p <- av$`Pr(>Chisq)`[2]
  
  pairs <- pairs(emmeans(glmer_out[[paste0('cluster_',i)]]$model, 'Status'), adjust = 'tukey')%>%
    as.data.frame()
  
  m_d_p.value <- pairs$p.value[1]
  m_f_p.value <- pairs$p.value[2]
  d_f_p.value <- pairs$p.value[3]
  
  newd <- data.frame(cluster = i,
                     m_d_p.value,
                     m_f_p.value,
                     d_f_p.value,
                     status_p = status_p,
                     singular= isSingular(glmer_out[[paste0('cluster_',i)]]$model)
  )
  glmer_av_output <- rbind(glmer_av_output, newd)
}

library(ggsignif)
#plot clusters 10, 13, 21

clust_10 = subset(ieg_positive, harmony.wnn_res0.4_clusters==10)
clust_10$prop_pos = clust_10$ieg_pos/(clust_10$ieg_pos+clust_10$ieg_neg)

clust_10_plot <- ggplot(clust_10, aes(x = Status,  y = prop_pos))+
  geom_boxplot(outlier.shape = NA, aes(color = Status, fill = Status), alpha = 0.25)+
  geom_jitter(fill = NA, shape = 1, size = 2)+
  theme_classic()+
  labs(y = 'Proportion IEG+', x = 'Sex', fill = 'Sex',color ="Sex", title ="Cluster 10")+
  theme(legend.position = 'none', plot.title = element_text(hjust=0.5))+
  theme(legend.position = 'none')+
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+
  geom_signif(xmin = 1, xmax = 3, y_position = 0.7, annotation = '*', color = 'black', textsize = 7, tip_length =c(0,0))+
  ylim(0,1.2*0.7)
clust_10_plot

ggsave(plot = clust_10_plot,
       file = "prop_ieg_clust_10.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 2")



clust_13 = subset(ieg_positive, harmony.wnn_res0.4_clusters==13)
clust_13$prop_pos = clust_13$ieg_pos/(clust_13$ieg_pos+clust_13$ieg_neg)

clust_13_plot <- ggplot(clust_13, aes(x = Status,  y = prop_pos))+
  geom_boxplot(outlier.shape = NA, aes(color = Status, fill = Status), alpha = 0.25)+
  geom_jitter(fill = NA, shape = 1, size = 2)+
  theme_classic()+
  labs(y = 'Proportion IEG+', x = 'Sex', fill = 'Sex',color ="Sex", title ="Cluster 13")+
  theme(legend.position = 'none', plot.title = element_text(hjust=0.5))+
  theme(legend.position = 'none')+
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+
  geom_signif(xmin = 1, xmax = 3, y_position = 1.2, annotation = '*', color = 'black', textsize = 7, tip_length =c(0,0))+
  ylim(0,1.2*1.2)
clust_13_plot

ggsave(plot = clust_13_plot,
       file = "prop_ieg_clust_13.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 2")



clust_21 = subset(ieg_positive, harmony.wnn_res0.4_clusters==21)
clust_21$prop_pos = clust_21$ieg_pos/(clust_21$ieg_pos+clust_21$ieg_neg)

clust_21_plot <- ggplot(clust_21, aes(x = Status,  y = prop_pos))+
  geom_boxplot(outlier.shape = NA, aes(color = Status, fill = Status), alpha = 0.25)+
  geom_jitter(fill = NA, shape = 1, size = 2)+
  theme_classic()+
  labs(y = 'Proportion IEG+', x = 'Sex', fill = 'Sex',color ="Sex", title ="Cluster 21")+
  theme(legend.position = 'none', plot.title = element_text(hjust=0.5))+
  theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+
  geom_signif(xmin = 1, xmax = 1.9, y_position = .8, annotation = '*', color = 'black', textsize = 7, tip_length =c(0,0))+
  geom_signif(xmin = 2.1, xmax = 3, y_position = .8, annotation = '*', color = 'black', textsize = 7, tip_length =c(0,0))+
  ylim(0,1.2*.8)
clust_21_plot

ggsave(plot = clust_21_plot,
       file = "prop_ieg_clust_21.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 2")



