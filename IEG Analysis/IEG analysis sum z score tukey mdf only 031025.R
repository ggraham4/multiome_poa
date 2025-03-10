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

neuronal_only$ieg_z_score_sum = z_matrix_cells_as_rows$sum_z_score

FeaturePlot(neuronal_only, 'ieg_z_score_sum')

data_for_plot = subset(neuronal_only@meta.data, Status %in% c("M",'D', 'E','NF',"F"))%>%
  group_by(harmony.wnn_res0.4_clusters, individual, Status)%>%
  summarize(mean_z_score_sum = mean(ieg_z_score_sum))

data_for_plot$Status = factor(data_for_plot$Status, levels = c("M",'D', 'E','NF',"F"))  
data_for_plot$label = ifelse(data_for_plot$individual=='C13F', 'C13F',NA)
ggplot(data_for_plot, aes(x = as.factor(harmony.wnn_res0.4_clusters), y = mean_z_score_sum, shape = Status, color = Status))+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  #geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Mean IEG Sum Z Score')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3,4,5))+
  geom_text(aes(label =label))


ggplot(data_for_plot, aes(x = Status, y = mean_z_score_sum, shape = Status, color = Status))+
  geom_violin()+
  geom_boxplot(alpha = 0, outlier.shape = NA)

### Stats ####
library(emmeans)
### redoing z score calculation to be within cluster #####

#first, lets subset to only IEG+ cells
ieg_only_cells <- neuronal_only%>%
  subset(ieg == 'ieg')

sum_z_score_stats_cluster <- data.frame()
for(j in unique(ieg_only_cells$harmony.wnn_res0.4_clusters)){
  print(j)
  temp_obj = subset(ieg_only_cells, harmony.wnn_res0.4_clusters ==j)
  
  #extract expression
  temp_counts = temp_obj$RNA$data%>%
    as.data.frame()%>%
    t()%>%
    as.tibble()%>%
    dplyr::select(all_of(ieg_like_genes))
  
  #calculate means and SD for each gene
  gene_measures = temp_counts%>%
    na.omit()%>%
    t()%>%
    as.data.frame()
  
  #ok I think now I make a new matrix and then do some matrix math
  gene_measures$SD = rowSds(gene_measures%>%as.matrix())
  gene_measures$mean = rowMeans(gene_measures%>%as.matrix())
  
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
  
  z_matrix_cells_as_rows = z_matrix%>%
    t()%>%
    as.data.frame()
  
  z_matrix_cells_as_rows$sum_z_score = rowSums(z_matrix_cells_as_rows)
  
  temp_obj$ieg_z_score_sum = z_matrix_cells_as_rows$sum_z_score
  
  ### now do the stats ###
  data_for_stats = data.frame(
    sum_z_score = temp_obj$ieg_z_score_sum,
    individual = temp_obj$individual,
    Status = temp_obj$Status,
    cluster = temp_obj$harmony.wnn_res0.4_clusters)%>%
    subset(Status %in%c('M','D','F'))
  
  if(length(unique(data_for_stats$sum_z_score))==1){next}
  model = lmer(sum_z_score~Status + (1|individual), data = data_for_stats)
  av = car::Anova(model, type = 'III')%>% as.data.frame()
  pairs = pairs(emmeans(model, 'Status'), adjust = 'tukey')%>%as.data.frame()
  
  temp_data = data.frame(
    cluster = j,
    singular = isSingular(model),
    av_p.value = av$`Pr(>Chisq)`[2],
    d_f_p.value = pairs$p.value[pairs$contrast== 'D - F'],
    d_m_p.value = pairs$p.value[pairs$contrast== 'D - M'],
    f_m_p.value = pairs$p.value[pairs$contrast== 'F - M']
  )
  sum_z_score_stats_cluster <- rbind(sum_z_score_stats_cluster, temp_data)
  
}

sum_z_score_stats_cluster$issignif = ifelse(sum_z_score_stats_cluster$av_p.value<0.05, '*', NA)
view(sum_z_score_stats_cluster)

###OK THIS IS SOMETHING THIS IS NOT A DRILL, CLUSTER 27 DOES DIFFER BY SEX
# also 0 and 23 - POA clusters!!! 
#> but, pairwise comparisons not significant with tukey, what about without?
#> without only 23 has a pairwise difference. In this case, I am just supposed to look at the data

###plot 27
obj_27 = subset(ieg_only_cells, harmony.wnn_res0.4_clusters =='27')

#extract expression
temp_counts = obj_27$RNA$data%>%
  as.data.frame()%>%
  t()%>%
  as.tibble()%>%
  dplyr::select(all_of(ieg_like_genes))

#calculate means and SD for each gene
gene_measures = temp_counts%>%
  na.omit()%>%
  t()%>%
  as.data.frame()

#ok I think now I make a new matrix and then do some matrix math
gene_measures$SD = rowSds(gene_measures%>%as.matrix())
gene_measures$mean = rowMeans(gene_measures%>%as.matrix())

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

colnames(z_matrix) = colnames(obj_27$RNA$data)
rownames(z_matrix) = colnames(temp_counts)

z_matrix = as.matrix(z_matrix)

z_matrix_cells_as_rows = z_matrix%>%
  t()%>%
  as.data.frame()

z_matrix_cells_as_rows$sum_z_score = rowSums(z_matrix_cells_as_rows)

obj_27$ieg_z_score_sum = z_matrix_cells_as_rows$sum_z_score

### now do the stats ###
data_for_stats = data.frame(
  sum_z_score = obj_27$ieg_z_score_sum,
  individual = obj_27$individual,
  Status = obj_27$Status,
  cluster = obj_27$harmony.wnn_res0.4_clusters)%>%
  group_by(individual, Status)%>%
  summarize(mean_ieg_score = mean(sum_z_score),
            se_ieg_score = sd(sum_z_score)/sqrt(n())
  )%>%
  subset(Status != 'NRM')

data_for_stats$Status <- factor(data_for_stats$Status, levels = c('M','D','E','NF','F'))
ggplot(data_for_stats, aes(x =Status, y =mean_ieg_score, color = Status, shape = Status, group = Status))+
  geom_boxplot(fill = NA)+
  geom_pointrange(aes(x = Status,
                      y =mean_ieg_score,
                      ymin = mean_ieg_score- se_ieg_score,
                      ymax = mean_ieg_score+ se_ieg_score 
  ),position = position_jitterdodge(dodge.width = .7))

###plot 23
obj_23 = subset(ieg_only_cells, harmony.wnn_res0.4_clusters =='23')

#extract expression
temp_counts = obj_23$RNA$data%>%
  as.data.frame()%>%
  t()%>%
  as.tibble()%>%
  dplyr::select(all_of(ieg_like_genes))

#calculate means and SD for each gene
gene_measures = temp_counts%>%
  na.omit()%>%
  t()%>%
  as.data.frame()

#ok I think now I make a new matrix and then do some matrix math
gene_measures$SD = rowSds(gene_measures%>%as.matrix())
gene_measures$mean = rowMeans(gene_measures%>%as.matrix())

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

colnames(z_matrix) = colnames(obj_23$RNA$data)
rownames(z_matrix) = colnames(temp_counts)

z_matrix = as.matrix(z_matrix)

z_matrix_cells_as_rows = z_matrix%>%
  t()%>%
  as.data.frame()

z_matrix_cells_as_rows$sum_z_score = rowSums(z_matrix_cells_as_rows)

obj_23$ieg_z_score_sum = z_matrix_cells_as_rows$sum_z_score

### now do the stats ###
data_for_stats = data.frame(
  sum_z_score = obj_23$ieg_z_score_sum,
  individual = obj_23$individual,
  Status = obj_23$Status,
  cluster = obj_23$harmony.wnn_res0.4_clusters)%>%
  group_by(individual, Status)%>%
  summarize(mean_ieg_score = mean(sum_z_score),
            se_ieg_score = sd(sum_z_score)/sqrt(n())
  )%>%
  subset(Status != 'NRM')

data_for_stats$Status <- factor(data_for_stats$Status, levels = c('M','D','E','NF','F'))
ggplot(data_for_stats, aes(x =Status, y =mean_ieg_score, color = Status, shape = Status, group = Status))+
  geom_boxplot(fill = NA)+
  geom_pointrange(aes(x = Status,
                      y =mean_ieg_score,
                      ymin = mean_ieg_score- se_ieg_score,
                      ymax = mean_ieg_score+ se_ieg_score 
  ),position = position_jitterdodge(dodge.width = .7))

###plot 0
obj_0 = subset(ieg_only_cells, harmony.wnn_res0.4_clusters =='0')

#extract expression
temp_counts = obj_0$RNA$data%>%
  as.data.frame()%>%
  t()%>%
  as.tibble()%>%
  dplyr::select(all_of(ieg_like_genes))

#calculate means and SD for each gene
gene_measures = temp_counts%>%
  na.omit()%>%
  t()%>%
  as.data.frame()

#ok I think now I make a new matrix and then do some matrix math
gene_measures$SD = rowSds(gene_measures%>%as.matrix())
gene_measures$mean = rowMeans(gene_measures%>%as.matrix())

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

colnames(z_matrix) = colnames(obj_0$RNA$data)
rownames(z_matrix) = colnames(temp_counts)

z_matrix = as.matrix(z_matrix)

z_matrix_cells_as_rows = z_matrix%>%
  t()%>%
  as.data.frame()

z_matrix_cells_as_rows$sum_z_score = rowSums(z_matrix_cells_as_rows)

obj_0$ieg_z_score_sum = z_matrix_cells_as_rows$sum_z_score

### now do the stats ###
data_for_stats = data.frame(
  sum_z_score = obj_0$ieg_z_score_sum,
  individual = obj_0$individual,
  Status = obj_0$Status,
  cluster = obj_0$harmony.wnn_res0.4_clusters)%>%
  group_by(individual, Status)%>%
  summarize(mean_ieg_score = mean(sum_z_score),
            se_ieg_score = sd(sum_z_score)/sqrt(n())
  )%>%
  subset(Status != 'NRM')

data_for_stats$Status <- factor(data_for_stats$Status, levels = c('M','D','E','NF','F'))
ggplot(data_for_stats, aes(x =Status, y =mean_ieg_score, color = Status, shape = Status, group = Status))+
  geom_boxplot(fill = NA)+
  geom_pointrange(aes(x = Status,
                      y =mean_ieg_score,
                      ymin = mean_ieg_score- se_ieg_score,
                      ymax = mean_ieg_score+ se_ieg_score 
  ),position = position_jitterdodge(dodge.width = .7))




