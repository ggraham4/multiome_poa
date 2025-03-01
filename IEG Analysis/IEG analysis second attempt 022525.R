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


###I think there is also a way to do this with a module score
iegs_list <- list()
iegs_list[[paste0('iegs')]] <- c(ieg_like_genes)
neuronal_only <- AddModuleScore(neuronal_only,features =  iegs_list, name = 'iegs')

DotPlot(neuronal_only, features ='iegs1')+
  coord_flip() #this doesnt seem right to me, it doesnt agree with the previous graph

###--- what the fuck is this beta binomial thing zack does ---####
###--- Ill figure that out later, for now I will do what coltan did ---####

##ok IEG score, the number of IEGs a nucleus expresses

RNA_data <- neuronal_only@assays$RNA$data
ieg_rows_rna <- RNA_data[which(rownames(RNA_data) %in% ieg_like_genes),]
ieg_scores <- colSums(ieg_rows_rna>0)

neuronal_only$ieg_scores = ieg_scores

FeaturePlot(neuronal_only, feature = 'ieg_scores')

##plot ieg_scores by sex
ieg_scores <- neuronal_only@meta.data%>%
  group_by(individual, Status, harmony.wnn_res0.4_clusters)%>%
  summarize(mean_ieg_score = mean(ieg_scores))%>%
  subset(Status %in% c('M','D','F'))

ieg_scores$Status <- factor(ieg_scores$Status, levels = c('M','D','F'))
ggplot(ieg_scores, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = mean_ieg_score, shape = Status, color = Status))+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Mean IEG Score')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))

ieg_positive <- neuronal_only@meta.data%>%
  group_by(individual, Status, harmony.wnn_res0.4_clusters)%>%
  summarize(ieg_pos = sum(ieg_scores>0),
            ieg_neg = sum(ieg_scores ==0),
            prop_pos = ieg_pos/ieg_neg)%>%
  subset(Status %in% c('M','D','F'))%>%
  filter(prop_pos != Inf)

ieg_positive$Status <- factor(ieg_positive$Status, levels = c('M','D','F'))
ggplot(ieg_positive, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = prop_pos, shape = Status, color = Status))+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Mean IEG Score')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))

###--- Statistics ---####
## based on coltans paper, negative binomial regression on ieg score
neg_binom_list <- list()
for(i in unique(neuronal_only$harmony.wnn_res0.4_clusters)){
  print(i)
neg_binom <- glmer.nb(ieg_scores ~ Status + log10GenesPerUMI +(1|individual), data = subset(neuronal_only@meta.data, Status %in% c('M','D','F')& harmony.wnn_res0.4_clusters==i))

neg_binom_list[[paste0(i)]]$summary <-summary(neg_binom)
neg_binom_list[[paste0(i)]]$av <-car::Anova(neg_binom, type = 'III')

}

av_p_values <- data.frame()
for(i in unique(neuronal_only$harmony.wnn_res0.4_clusters)){
  
 av <- neg_binom_list[[paste0(i)]]$av%>%as.data.frame()
 
 newd <- data.frame(cluster = i, 
                    Status_av_p = av$`Pr(>Chisq)`[2])
 av_p_values <- rbind(av_p_values, newd)
   
}
#nothing is going to be significant with a q value adjustment


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
  
  pairs <- pairs(emmeans(glmer_out[[paste0('cluster_',i)]]$model, 'Status'), adjust = 'none')%>%
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

glmer_av_output$m_d_q.value <- p.adjust(glmer_av_output$m_d_p.value, 'fdr', nrow(glmer_av_output))
glmer_av_output$m_f_q.value <- p.adjust(glmer_av_output$m_f_p.value, 'fdr', nrow(glmer_av_output))
glmer_av_output$d_f_q.value <- p.adjust(glmer_av_output$d_f_p.value, 'fdr', nrow(glmer_av_output))
glmer_av_output$status_q <- p.adjust(glmer_av_output$status_p, 'fdr', nrow(glmer_av_output))

glmer_av_output$issignif <- ifelse(glmer_av_output$m_d_q.value<0.05|
                                     glmer_av_output$m_f_q.value<0.05|
                                     glmer_av_output$d_f_q.value<0.05|
                                   glmer_av_output$status_q<0.05,
                                   '*',
                                   NA)
glmer_av_output$cluster <- as.numeric(glmer_av_output$cluster)
### nothing is significant

library(sjPlot)
plot_model(glmer_out[[paste0('cluster_','27')]]$model)


###trying out zack's beta binomial
library(PROreg)
for(i in unique(neuronal_only$harmony.wnn_res0.4_clusters)){
  temp_data <- subset(neuronal_only@meta.data, harmony.wnn_res0.4_clusters==i)%>%
    dplyr::select(individual, Status, ieg_scores)%>%
    subset(Status %in% c('M','D','F'))

  
model <- BBmm(fixed.formula = 
                temp_data$ieg_scores ~ temp_data$Status,
  random.formula =  ~temp_data$individual,
  data = temp_data,
  m = 21)  #max number of iegs


summary(model)

}
### the problem with this model is I cannot compute typeIII tests or pairwise comparisons


### what about glmmTMB
library(glmmTMB)
library(performance)
out_betabin <- data.frame()
for(i in unique(neuronal_only$harmony.wnn_res0.4_clusters)){
  print(i)
  temp_data <- subset(neuronal_only@meta.data, harmony.wnn_res0.4_clusters==i)%>%
    dplyr::select(individual, Status, ieg_scores)%>%
    subset(Status %in% c('M','D','F'))

  
model <- glmmTMB(cbind(ieg_scores,21-ieg_scores)  ~ Status + (1|individual),
  data = temp_data,
  family=betabinomial(link = "logit"))  

#hist(temp_data$ieg_scores)

summary(model)
pairs <- pairs(emmeans(model, 'Status'), adjust = 'none')%>%
  as.data.frame()

av <- car::Anova(model, type = 'III')

newd <- data.frame(cluster = i,
                   status_p.value = av$`Pr(>Chisq)`[2],
                   d_m_p.value = pairs$p.value[pairs$contrast=='D - M'],
                  f_m_p.value = pairs$p.value[pairs$contrast=='F - M'],
                   d_f_p.value = pairs$p.value[pairs$contrast=='D - F'],
                  singular =  check_singularity(model)

                  
)
out_betabin <- rbind(newd, out_betabin)


}
### ok this is definitely the way to go

out_betabin$status_q.value <- p.adjust(out_betabin$status_p.value, 'fdr', nrow(out_betabin))
out_betabin$d_m_q.value <- p.adjust(out_betabin$d_m_p.value, 'fdr', nrow(out_betabin))
out_betabin$f_m_q.value <- p.adjust(out_betabin$f_m_p.value, 'fdr', nrow(out_betabin))
out_betabin$d_f_q.value <- p.adjust(out_betabin$d_f_p.value, 'fdr', nrow(out_betabin))

out_betabin$issignif <- ifelse(out_betabin$status_q.value<0.05|
                                 out_betabin$d_m_q.value<0.05|
                                 out_betabin$f_m_q.value<0.05|
                                 out_betabin$d_f_q.value<0.05,
                               '*',
                               NA)

### nothing signif man 

out_betabin$lower_string_issignif <- ifelse(out_betabin$status_p.value<0.05|
                                 out_betabin$d_m_p.value<0.05|
                                 out_betabin$f_m_p.value<0.05|
                                 out_betabin$d_f_p.value<0.05,
                               '*',
                               NA)

out_betabin$issignif_0.1 <- ifelse(out_betabin$status_q.value<0.1|
                                 out_betabin$d_m_q.value<0.1|
                                 out_betabin$f_m_q.value<0.1|
                                 out_betabin$d_f_q.value<0.1,
                               '*',
                               NA)


