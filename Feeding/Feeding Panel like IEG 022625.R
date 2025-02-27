#> What I want to do is similar to zack's IEG panel, I want to look at genes
#>releted to feeding, and see how they differ between sexes and clusters
#>
#>Later, I will look at genes correlated with size change, and see how they 
#>overlap with the genes I find here
#>
#>I am going to be very conservative with this panel, using the review linked 
#>below, I will only pick genes that are positively correlated with feeding
#>
#>https://doi.org/10.1089/zeb.2006.3.13

## libraries
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

### mark nuclei as feeding+ or feeding-

neuronal_only$feeding <- ifelse(neuronal_only@assays$RNA$data['npy',] >0, 'feeding', NA)
neuronal_only$feeding <- ifelse(neuronal_only@assays$RNA$data['gal',] >0, 'feeding', neuronal_only$feeding)
neuronal_only$feeding <- ifelse(neuronal_only@assays$RNA$data['hcrt',] >0, 'feeding', neuronal_only$feeding)
neuronal_only$feeding <- ifelse(neuronal_only@assays$RNA$data['agrp',] >0, 'feeding', neuronal_only$feeding)

#i will subset to npy, galanin, agouti related protein, orexin (hcrt is orexin precursor)

neuronal_only$feeding <- ifelse(is.na(neuronal_only$feeding), 'no_feeding', neuronal_only$feeding)

DimPlot(neuronal_only, split.by = 'feeding')
### interesting

### before I find all the coexpressed genes, lets just look at these 4
fed_data_curated_cluster_level <- table(neuronal_only@meta.data$feeding, neuronal_only@meta.data$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  pivot_wider(names_from = 'Var1', values_from = 'Freq')%>%
  mutate(prop_pos = feeding/(feeding+no_feeding))

ggplot(fed_data_curated_cluster_level, aes(x = as.factor(Var2), y = prop_pos))+
  geom_bar(stat = 'identity')+
  geom_hline(yintercept = 0.5, linetype = 2)+
  labs(x = 'Cluster', y = 'Prop Feeding +')+
  theme_minimal()


feeding_data_3_individual_level <- neuronal_only@meta.data%>%
  dplyr::select(c(individual, Status, harmony.wnn_res0.4_clusters, feeding))%>%
  group_by(Status, individual, harmony.wnn_res0.4_clusters)%>%
  summarize(feed_pos = sum(feeding =='feeding'),
            feed_neg = sum(feeding == 'no_feeding')
  )%>%
  mutate(prop_pos = feed_pos/(feed_pos+feed_neg))%>%
  subset(Status %in% c('M','D','F'))%>%
  na.omit()%>%
  subset(prop_pos != Inf)

feeding_data_3_individual_level$Status <- factor(feeding_data_3_individual_level$Status, levels = c('M','D','F'))

ggplot(feeding_data_3_individual_level, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = prop_pos, shape = Status, color = Status))+
    geom_hline(yintercept = 0.5, linetype = 2)+
    geom_hline(yintercept =1, linetype = 2)+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Proportion Feed+')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))+
  scale_y_continuous(breaks = c(0.5,0:7))
#nothing screams significant to me other than the fact cluster 3 looks interesting


## ok I know cluster 20 is driven by it being npy+, finding coexpressed genes will be more interesting
feeding_marker_all <-data.frame()
for(i in unique(neuronal_only$harmony.wnn_res0.4_clusters)){
  print(i)
  temp_obj <- subset(neuronal_only, harmony.wnn_res0.4_clusters==i)
  feeding_markers <- FindAllMarkers(temp_obj, 
                              assay = 'RNA', 
                              group.by = 'feeding',
                              logfc.threshold = 0,
                              min.pct = 1/204) #204 cells in smallest cluster
  feeding_markers$harmony.wnn_res0.4_clusters =i
  feeding_marker_all <- rbind(feeding_marker_all, feeding_markers)
}

#extract significant feeding like genes
feeding_marker_all_signif <- feeding_marker_all[feeding_marker_all$p_val_adj<0.05 &feeding_marker_all$cluster=='feeding' ,]

#find genes significant in over half of clusters
feeding_like_genes_data <- table(feeding_marker_all_signif$gene,
      feeding_marker_all_signif$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  group_by(Var1)%>% #cluster
  summarize(n_clusters = sum(Freq))%>%
  subset(n_clusters > length(unique(neuronal_only$harmony.wnn_res0.4_clusters))/2
  ) 

### whelp I guess it ends here huh 
### maybe I need to revisit which genes I pick

#### What about genes negatively correlated with feeding, pomc and cart
FeaturePlot(neuronal_only, 'LOC111584656') ### pomc
FeaturePlot(neuronal_only, 'cart2') ### cart2
FeaturePlot(neuronal_only, 'cart3') ### cart3
FeaturePlot(neuronal_only, 'cart4') ### cart4
FeaturePlot(neuronal_only, 'LOC111573467') ### cart like
FeaturePlot(neuronal_only, 'cartl') ### cart like

neuronal_only$satiety <- ifelse(neuronal_only@assays$RNA$data['LOC111584656',] >0, 'satiety', NA)
neuronal_only$satiety <- ifelse(neuronal_only@assays$RNA$data['cart2',] >0, 'satiety', neuronal_only$satiety)
neuronal_only$satiety <- ifelse(neuronal_only@assays$RNA$data['cart3',] >0, 'satiety', neuronal_only$satiety)
neuronal_only$satiety <- ifelse(neuronal_only@assays$RNA$data['cart4',] >0, 'satiety', neuronal_only$satiety)
neuronal_only$satiety <- ifelse(neuronal_only@assays$RNA$data['LOC111573467',] >0, 'satiety', neuronal_only$satiety)
neuronal_only$satiety <- ifelse(neuronal_only@assays$RNA$data['cartl',] >0, 'satiety', neuronal_only$satiety)

neuronal_only$satiety <- ifelse(is.na(neuronal_only$satiety), 'no_satiety', neuronal_only$satiety)

DimPlot(neuronal_only, split.by='satiety')

### before I find all the coexpressed genes, lets just look at these 
sat_data_curated_cluster_level <- table(neuronal_only@meta.data$satiety, neuronal_only@meta.data$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  pivot_wider(names_from = 'Var1', values_from = 'Freq')%>%
  mutate(prop_pos = satiety/(satiety+no_satiety))

ggplot(sat_data_curated_cluster_level, aes(x = as.factor(Var2), y = prop_pos))+
  geom_bar(stat = 'identity')+
  geom_hline(yintercept = 0.5, linetype = 2)+
  labs(x = 'Cluster', y = 'Prop Satiety +')+
  theme_minimal()
#what in the fuck 31 lmfao

satiety_data_3_individual_level <- neuronal_only@meta.data%>%
  dplyr::select(c(individual, Status, harmony.wnn_res0.4_clusters, satiety))%>%
  group_by(Status, individual, harmony.wnn_res0.4_clusters)%>%
  summarize(sat_pos = sum(satiety =='satiety'),
            sat_neg = sum(satiety == 'no_satiety')
  )%>%
  mutate(prop_pos = sat_pos/(sat_pos+sat_neg))%>%
  subset(Status %in% c('M','D','F'))%>%
  na.omit()%>%
  subset(prop_pos != Inf)

satiety_data_3_individual_level$Status <- factor(satiety_data_3_individual_level$Status, levels = c('M','D','F'))

ggplot(satiety_data_3_individual_level, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = prop_pos, shape = Status, color = Status))+
    geom_hline(yintercept = 0.5, linetype = 2)+
    geom_hline(yintercept =1, linetype = 2)+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Proportion Feed+')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))+
  scale_y_continuous(breaks = c(0.5,0:7))
#looks like an increase in 19 and 17  but probably not significant

## anyway, its interesting 31 is a satiety associated cluster, does that mark somethign?

#coexpresed 
satiety_marker_all <-data.frame()
for(i in unique(neuronal_only$harmony.wnn_res0.4_clusters)){
  print(i)
  temp_obj <- subset(neuronal_only, harmony.wnn_res0.4_clusters==i)
  satiety_markers <- FindAllMarkers(temp_obj, 
                              assay = 'RNA', 
                              group.by = 'satiety',
                              logfc.threshold = 0,
                              min.pct = 1/204) #204 cells in smallest cluster
  satiety_markers$harmony.wnn_res0.4_clusters =i
  satiety_marker_all <- rbind(satiety_marker_all, satiety_markers)
}

#extract significant satiety like genes
satiety_marker_all_signif <- satiety_marker_all[satiety_marker_all$p_val_adj<0.05 &satiety_marker_all$cluster=='satiety' ,]

#find genes significant in over half of clusters
satiety_like_genes_data <- table(satiety_marker_all_signif$gene,
      satiety_marker_all_signif$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  group_by(Var1)%>% #cluster
  summarize(n_clusters = sum(Freq))%>%
  subset(n_clusters > length(unique(neuronal_only$harmony.wnn_res0.4_clusters))/4
  ) 

### nothin..

 table(satiety_marker_all_signif$gene,
      satiety_marker_all_signif$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  group_by(Var1)%>% #cluster
  summarize(n_clusters = sum(Freq))%>%
  subset(n_clusters > length(unique(neuronal_only$harmony.wnn_res0.4_clusters))/4
  ) 
 ##ok well vgf is here so thats something
 
 #im going to move forward with the 1/4 threshold just to see

 ##ok satiety score, the number of satietys a nucleus expresses
 satiety_like_genes <- satiety_like_genes_data$Var1

RNA_data <- neuronal_only@assays$RNA$data
satiety_rows_rna <- RNA_data[which(rownames(RNA_data) %in% satiety_like_genes),]
satiety_scores <- colSums(satiety_rows_rna>0)

neuronal_only$satiety_scores = satiety_scores

FeaturePlot(neuronal_only, feature = 'satiety_scores')

##plot satiety_scores by sex
satiety_scores <- neuronal_only@meta.data%>%
  group_by(individual, Status, harmony.wnn_res0.4_clusters)%>%
  summarize(mean_satiety_score = mean(satiety_scores))%>%
  subset(Status %in% c('M','D','F'))

satiety_scores$Status <- factor(satiety_scores$Status, levels = c('M','D','F'))
ggplot(satiety_scores, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = mean_satiety_score, shape = Status, color = Status))+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Mean satiety Score')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))
### im not convinced theres anything here


satiety_positive <- neuronal_only@meta.data%>%
  group_by(individual, Status, harmony.wnn_res0.4_clusters)%>%
  summarize(satiety_pos = sum(satiety_scores>0),
            satiety_neg = sum(satiety_scores ==0),
            prop_pos = satiety_pos/(satiety_pos+satiety_neg))%>%
  subset(Status %in% c('M','D','F'))%>%
  filter(prop_pos != Inf)

satiety_positive$Status <- factor(satiety_positive$Status, levels = c('M','D','F'))
ggplot(satiety_positive, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = prop_pos, shape = Status, color = Status))+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Proportion Nuclei Satiety+')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))
#one COULD argue that 27 differs in terms of proportion of nucleu that are satuety positive but I think thats really pushing it
# i feel like there might be something significant here


##### Statistcs ---######
#ok first I want to use a logistic to look at satiety+ vs satiety-  nuclei across clusters, 
#that to me is the most likely to be somethign

glmer_out <- list()
for(i in unique(neuronal_only$harmony.wnn_res0.4_clusters)){
  print(i)
temp_satiety_positive <- subset(satiety_positive, harmony.wnn_res0.4_clusters==i)
satiety_pos_matrix <- matrix(cbind(temp_satiety_positive$satiety_pos, temp_satiety_positive$satiety_neg), nrow(temp_satiety_positive),2)

glmer_model <- glmer(satiety_pos_matrix~Status + (1|individual), data = temp_satiety_positive, family = binomial('logit'))

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


glmer_av_output$issignif_0.1 <- ifelse(glmer_av_output$m_d_q.value<0.1|
                                     glmer_av_output$m_f_q.value<0.1|
                                     glmer_av_output$d_f_q.value<0.1|
                                   glmer_av_output$status_q<0.1,
                                   '*',
                                   NA)

### theres nothing here

glmer_av_output$issignif_p_0.5 <- ifelse(glmer_av_output$m_d_p.value<0.05|
                                     glmer_av_output$m_f_p.value<0.05|
                                     glmer_av_output$d_f_p.value<0.05|
                                   glmer_av_output$status_p<0.05,
                                   '*',
                                   NA)

#theres nothing here, but its interesting that 27 is p value significant, males are significantly 
#> lower than dominants which very weakly supports my idea, perhaps its time to look at genes
#> correlated with 
#> 
#> I just want to plot the significant ones to see if theres a trend

satiety_positive_subset <- subset(satiety_positive, harmony.wnn_res0.4_clusters%in%c(glmer_av_output$cluster[glmer_av_output$issignif_p_0.5=='*']))
ggplot(satiety_positive_subset, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = prop_pos, shape = Status, color = Status))+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Proportion Nuclei Satiety+')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))

## im not going to continue with this analysis, Im already reaching

### ok one more, collapsing across all clusters
satiety_pos_matrix <- matrix(cbind(satiety_positive$satiety_pos, satiety_positive$satiety_neg), nrow(satiety_positive),2)

glmer_model_all <- glmer(satiety_pos_matrix~Status + (1|individual), data = satiety_positive, family = binomial('logit'))
summary(glmer_model_all)
pairs(emmeans(glmer_model_all, 'Status'), adjust ='none')
#nope 





























