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

#### Summing Z-scores ####
#### Ok for each nucleus, lets find the z score of each ieg relative to the mean of that cluster IG

# z = (value- average)/ SD

#first, lets subset to only IEG+ cells
ieg_only_cells <- neuronal_only%>%
  subset(ieg == 'ieg')

## lets do test case with the whole umap
temp_obj = ieg_only_cells

#extract expression
temp_counts = temp_obj$RNA$data%>%
  as.data.frame()%>%
  t()%>%
  as.tibble()%>%
  dplyr::select_(.dots = ieg_like_genes)

#calculate means and SD for each gene
gene_measures = temp_counts%>%
  t()%>%
  as.data.frame()

gene_measures$SD = rowSds(gene_measures%>%as.matrix())
gene_measures$mean = rowMeans(gene_measures%>%as.matrix())

#ok I think now I make a new matrix and then do some matrix math

z_matrix = matrix(NA, 
                  nrow = nrow(gene_measures), 
                  ncol = ncol(gene_measures) - 2) 

for (i in 1:nrow(gene_measures)) {
  gene_data = as.numeric(gene_measures[i, 1:(ncol(gene_measures)-2)])
  gene_mean = gene_measures[i, "mean"]
  gene_sd = gene_measures[i, "SD"]
  
  # Calculate z-scores for this gene across all cells
  z_matrix[i, ] = (gene_data - gene_mean) / gene_sd
}

colnames(z_matrix) = colnames(temp_obj$RNA$data)
rownames(z_matrix) = colnames(temp_counts)

z_matrix = as.matrix(z_matrix)

hist(z_matrix[5,])

z_matrix_cells_as_rows = z_matrix%>%
  t()%>%
  as.data.frame()
  
z_matrix_cells_as_rows$sum_z_score = rowSums(z_matrix_cells_as_rows)

hist(z_matrix_cells_as_rows$sum_z_score)

temp_obj$ieg_z_score_sum = z_matrix_cells_as_rows$sum_z_score

FeaturePlot(temp_obj, 'ieg_z_score_sum')

data_for_plot = subset(temp_obj@meta.data, Status %in% c("M",'D',"F"))%>%
  group_by(harmony.wnn_res0.4_clusters, individual, Status)%>%
  summarize(mean_z_score_sum = mean(ieg_z_score_sum))
      
data_for_plot$Status = factor(data_for_plot$Status, levels = c('M','D','F'))                   
ggplot(data_for_plot, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = mean_z_score_sum, shape = Status, color = Status))+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Mean IEG Sum Z Score')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))

data_for_plot2 = subset(temp_obj@meta.data, Status %in% c("M",'D', 'E',"F"))%>%
  group_by(harmony.wnn_res0.4_clusters, Status)%>%
  summarize(mean_z_score_sum = mean(ieg_z_score_sum),
            se_z_score_sum = sd(ieg_z_score_sum)/sqrt(n()))
data_for_plot2$Status = factor(data_for_plot2$Status, levels = c('M','D', 'E','F'))                   

box_height =data_for_plot2%>%
  group_by(harmony.wnn_res0.4_clusters)%>%
  summarize(up = max(mean_z_score_sum)+max(se_z_score_sum),
            down = min(mean_z_score_sum)-max(se_z_score_sum))

ggplot(data_for_plot2, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = mean_z_score_sum))+
  geom_rect( aes(xmin = c(((1:92)/4)),
                 xmax = (c(1:92)/4)+0.5,
                 ymax = rep(box_height$up, each = 4), 
                 ymin = rep(box_height$down, each = 4)),
             fill = 'grey')+
  geom_hline(yintercept = 0, linetype =2)+
  geom_pointrange(aes(x = as.factor(harmony.wnn_res0.4_clusters),
                      y =mean_z_score_sum,
                      ymin = mean_z_score_sum -se_z_score_sum,
                      ymax = mean_z_score_sum +se_z_score_sum,
                     color = Status),
                  position = position_dodge(0.75))  +
  labs(x = 'Cluster', y = 'Mean IEG Sum Z Score')+
  theme_minimal()




