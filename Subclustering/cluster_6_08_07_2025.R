### Subclustering Cluster 6 
#> GJG 08 07 2025
#> 
#> Cluster 6 is our main POA cluster of interest, so I want to know how it 
#> subclusters
#> 

#libraries
{
  library(ggridges)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(CytoTRACE)
  library(BiocManager)
  library(lme4)
  library(enrichR)
  `%notin%` = Negate(`%in%`)
  library(car)
  library(emmeans)
  
  clown_go = readRDS("Functions/clown_go2")  
  mecd = readRDS("Functions/mean_expression_cluster_data.rds")
}

#define functions
individual_mean_expression = function(object, gene, clustering= 'sub.cluster'){
  library(stringr)    
  options(dplyr.summarise.inform = FALSE)
  
  counts <- t(as.matrix(object@assays$RNA$data))
  Counts_of_interest <- as.data.frame(counts[,gene])
  Counts_of_interest$expression <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual
  Counts_of_interest$Status <- object@meta.data$Status
  Counts_of_interest$cluster <- object@meta.data[[clustering]]
    
  results <-Counts_of_interest%>%
    group_by(individual, Status, cluster)%>%
    summarize(mean = mean(expression),
              se = sd(expression)/sqrt(n()))
  results$Sex <- results$Status
  
  
  results$Sex <- str_sub(results$individual, -1)
  results$Sex[results$individual == 'T17D'] = 'NF'
  results$Sex[results$individual == 'A12D'] = 'E'
  results$Sex[results$individual == 'T11D'] = 'E'
  results$Sex[results$individual == 'GH'] = 'NRM'
  return(results)
  
  
}

plot_mean_expression <- function(data_frame){
  data_frame$Sex = factor(data_frame$Sex, levels = c('NRM','M', 'D','E','NF',"F"))
  plot <- ggplot(data_frame, aes(x = cluster, y = mean, color = Sex))+
    geom_point(position = position_dodge(0.75))+
    geom_boxplot(alpha = 0)+
    theme_minimal()+
    labs(y = 'Mean Normalized Expression', x ='Sub Cluster')
  
  return(plot)
  
  
}

#read in object
obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

#subcluster
obj_subclustered_6 <- FindSubCluster(obj, 6, 'harmony.wsnn')
#reduce size
obj_6_only <- subset(obj_subclustered_6, sub.cluster %in% c('6_0','6_1','6_2','6_3'))
rm(obj)

DimPlot(obj_6_only, reduction = 'harmony_wnn.umap', group.by = 'sub.cluster', label =T)
#interesting that 5_2 is so disconnected

### 1) Markers of subclusters ####
subcluster_markers <- FindAllMarkers(obj_6_only, group.by = 'sub.cluster')

#measure percent enrichment
subcluster_markers$fold_diff <- subcluster_markers$pct.1/subcluster_markers$pct.2

#find the strongest markers
strong_markers <- subset(subcluster_markers, (p_val_adj < 0.001) & (fold_diff>1) & (pct.1>pct.2))

# are there any repeated markers
repeats = table(strong_markers$gene)%>%sort()%>%as.data.frame()%>%subset(Freq>1)
#% interesting that htr2aa is here and egr1

#remove repeats
strong_specific_markers <- subset(subcluster_markers, gene %notin% repeats$Var1)

## use GO to analyze strong specific markers ##
#6_0 
clown_go(strong_specific_markers$gene[strong_specific_markers$cluster=='6_0'])%>%
  dotplot()
#% not really sure what to make of this

#6_1 
clown_go(strong_specific_markers$gene[strong_specific_markers$cluster=='6_1'])%>%
  dotplot()
#%neuropeptides

#6_2
clown_go(strong_specific_markers$gene[strong_specific_markers$cluster=='6_2'])%>%
  dotplot()
#% interesting that this one is glutamatergic while 6_0 is gabaergic, in my
#% previous analysis, I found that 6 was a mixed population, I should look at that
#% again

#6_3
clown_go(strong_specific_markers$gene[strong_specific_markers$cluster=='6_3'])%>%
  dotplot()
#% back to gaba

### 2) Hormone and neuropeptide receptors ####

#hormones
plot_mean_expression(individual_mean_expression(obj_6_only, 'esr1'))
plot_mean_expression(individual_mean_expression(obj_6_only, 'esr2a'))
plot_mean_expression(individual_mean_expression(obj_6_only, 'esr2b'))

plot_mean_expression(individual_mean_expression(obj_6_only, 'pgr'))

plot_mean_expression(individual_mean_expression(obj_6_only, 'ar'))

#neuropeptides
plot_mean_expression(individual_mean_expression(obj_6_only, 'galr1a'))

plot_mean_expression(individual_mean_expression(obj_6_only, 'cckar'))
plot_mean_expression(individual_mean_expression(obj_6_only, 'cckbra'))
plot_mean_expression(individual_mean_expression(obj_6_only, 'cckbrb'))
#interesting difference between brb and bra
plot_mean_expression(individual_mean_expression(obj_6_only, 'gnrhr1'))

plot_mean_expression(individual_mean_expression(obj_6_only, 'avpr2aa'))

plot_mean_expression(individual_mean_expression(obj_6_only, 'tacr1a'))
plot_mean_expression(individual_mean_expression(obj_6_only, 'tacr3a'))
plot_mean_expression(individual_mean_expression(obj_6_only, 'tacr3l'))

plot_mean_expression(individual_mean_expression(obj_6_only, 'kiss1ra'))
plot_mean_expression(individual_mean_expression(obj_6_only, 'kiss1rb'))

plot_mean_expression(individual_mean_expression(obj_6_only, 'npy8ar'))
plot_mean_expression(individual_mean_expression(obj_6_only, 'npy2rl'))
plot_mean_expression(individual_mean_expression(obj_6_only, 'npy'))

plot_mean_expression(individual_mean_expression(obj_6_only, 'avp'))
plot_mean_expression(individual_mean_expression(obj_6_only, 'oxt'))


#conclusions:
#> 6_0 has the majority of the steroid hormone receptors
#> neuropeptides are more interesting. The interesting ones are in 0 and 1, 
#> interesting that oxt and avp are in 6_2, those cells are known to be 
#> where?
#> 
#> Godwin et al 2012
#> distributed along the third ventricle in the ventral portion of the preoptic 
#> area in the anterior hypothalamus. The parvocellular AVT neuron population is
#>  located at the rostral and ventral portion of the preoptic area with 
#>  typically large numbers of smaller neurons closely apposed to the third 
#>  ventricle
#>  
#>  Are these differences due to differential cell maturity
#>  
### 2) CytoTRACE ####

cyto= CytoTRACE(as.matrix(obj_6_only@assays$RNA$data))

obj_6_only$cyto = cyto$CytoTRACE

#plot
cyto_data <- obj_6_only@meta.data%>%
  group_by(sub.cluster, Status, individual)%>%
  summarize(cyt = mean(cyto),
            se_cyt =sd(cyto)/sqrt(n()))

cyto_data$Status <-  factor(cyto_data$Status, levels = c('NRM','M', 'D','E','NF',"F"))
ggplot(cyto_data, aes(x = sub.cluster, y = cyt, color = Status))+
  geom_point(position = position_dodge(0.75))+
  geom_boxplot(alpha = 0)+
  theme_minimal()+
  labs(y = 'Mean Cyto', x ='Sub Cluster')
#% 1 = most immature - 6_2
#% Looks like there are sex differences in 6_2 and potentially 6_0

## statistical testing
#6_2
data_6_2 = subset(obj_6_only@meta.data, sub.cluster == '6_2')

model_6_2 = lmer(cyto~Status +(1|individual), data= data_6_2)
car::Anova(model_6_2, type = 3)
#indeed it is significant
pairs(emmeans(model_6_2, 'Status'), pairwise = T, adjust ='none')
#DE, DF, F NF (somehow?)
pairs(emmeans(model_6_2, 'Status'), pairwise = T) 
#DF
##not sure weather I should use adjustment or not, but either way the DF is the 
#relevant result

#my hunch is also 6_0
data_6_0 = subset(obj_6_only@meta.data, sub.cluster == '6_0')

model_6_0 = lmer(cyto~Status +(1|individual), data= data_6_0)
car::Anova(model_6_0, type = 3) #nothin

#61 also possible
data_6_1 = subset(obj_6_only@meta.data, sub.cluster == '6_1')

model_6_1 = lmer(cyto~Status +(1|individual), data= data_6_1)
car::Anova(model_6_1, type = 3) #nope, surprising though

#I wonder what the distribution of cytotrace scores looks like between subclusters
obj_6_only@meta.data$Status =factor(obj_6_only@meta.data$Status, levels = c('NRM','M','D','E','NF','F'))
ggplot(obj_6_only@meta.data, aes(x = cyto, y = Status, fill = Status))+
  geom_density_ridges()+
  theme_ridges()+
  theme(legend.position = "none")

obj_6_only@meta.data$Status =factor(obj_6_only@meta.data$Status, levels = c('NRM','M','D','NF','E','F'))
ggplot(obj_6_only@meta.data, aes(x = cyto, y = Status, fill = Status))+
  geom_density_ridges()+
  theme_ridges()+
  theme(legend.position = "none")

#put it back
obj_6_only@meta.data$Status =factor(obj_6_only@meta.data$Status, levels = c('NRM','M','D','E','NF','F'))

# what does a test ignoring subcluster look like
model_6 = lmer(cyto~Status +(1|individual), data= obj_6_only@meta.data)
car::Anova(model_6, 3)
pairs(emmeans(model_6, 'Status'))

#because it is bounded between 0 and 1, is it better to use a non parametric
#have to take the mean per individual because of pseudoreplication
krusk_data = obj_6_only@meta.data%>%
  group_by(individual, Status)%>%
  summarize(mean_cyto = mean(cyto))

krusk_cyto = kruskal.test(mean_cyto~Status, data = krusk_data)
# wait it has the exact same p value as the anova of the
#lmer with random effect ? wow

#> conclusions:
#> across clusters, dominants have the most immature cells, though that is only
#> significant in cluster 6_2

### 3) Proportion Differences ####
sub_cells = obj_6_only@meta.data%>%
  group_by(sub.cluster, Status, individual)%>%
  subset(Status %in% c('M','D','F'))%>%
  summarize(cells = n())

total_cells = obj_6_only@meta.data%>%
  group_by(Status, individual)%>%
  subset(Status %in% c('M','D','F'))%>%
  summarize(total_cells = n())

full_data <- sub_cells%>%
  right_join(total_cells, by = 'individual')

full_data$diff = full_data$total_cells-full_data$cells

glmer_output = data.frame()
for(i in 0 :3){
  subset_data = subset(full_data, sub.cluster  ==paste0('6_',i) )
  
  glmer_matrix = matrix(c(subset_data$cells,
                          subset_data$diff),
                        nrow(subset_data), 2)
  
  glmer_model = glmer(glmer_matrix~Status.x + (1|individual),
                      family = binomial('logit'),
                      data = subset_data)
  
  av= car::Anova(glmer_model, type = 3)
  
  newd = data.frame(cluster = i,
                    av = av$`Pr(>Chisq)`[2],
                    singular = isSingular(glmer_model))
  glmer_output = rbind(glmer_output, newd)
  
}
#6_3 but is singular, I think that is tolerable though

##plot
full_data$percent = full_data$cells/full_data$total_cells
full_data$Status.x = factor(full_data$Status.x, levels = c('NRM','M', 'D','E','NF',"F"))

ggplot(full_data, aes(x = sub.cluster, y = percent, color = Status.x))+
  geom_point(position = position_dodge(0.75))+
  geom_boxplot(alpha = 0)+
  theme_minimal()
# fascinating so females have more of the differentiated 6_3, sad I dont see 
# the complementary differences from the previous analysis though

#pairwise comparisons for 6_3

i = 3
  subset_data = subset(full_data, sub.cluster  ==paste0('6_',i) )
  
  glmer_matrix = matrix(c(subset_data$cells,
                          subset_data$diff),
                        nrow(subset_data), 2)
  
  glmer_model = glmer(glmer_matrix~Status.x + (1|individual),
                      family = binomial('logit'),
                      data = subset_data)
  
  av= car::Anova(glmer_model, type = 3)

  pairs(emmeans(glmer_model, 'Status.x')) # male female 

### 4) NT Expression ####
  
#split cells into glut, gaba , and mixed
obj_6_only$primary_neurotransmitter <- ifelse(((obj_6_only@assays$RNA$data['LOC111588076',]>0 
                                                | obj_6_only@assays$RNA$data['gad2',]>0) &
                                                 (obj_6_only@assays$RNA$data['LOC111584103',]==0 & 
                                                    obj_6_only@assays$RNA$data['slc17a6b',]==0 &  
                                                    obj_6_only@assays$RNA$data['slc17a7a',]==0)), 'GABA',NA)

obj_6_only$primary_neurotransmitter <- ifelse(((obj_6_only@assays$RNA$data['LOC111588076',]==0 
                                                & obj_6_only@assays$RNA$data['gad2',]==0) &
                                                 (obj_6_only@assays$RNA$data['LOC111584103',]>0 | 
                                                    obj_6_only@assays$RNA$data['slc17a6b',]>0 |  
                                                    obj_6_only@assays$RNA$data['slc17a7a',]>0)), 'GLUT',obj_6_only$primary_neurotransmitter)

obj_6_only$primary_neurotransmitter <- ifelse(((obj_6_only@assays$RNA$data['LOC111588076',]>0 
                                                 | obj_6_only@assays$RNA$data['gad2',]>0) &
                                                  (obj_6_only@assays$RNA$data['LOC111584103',]>0 | 
                                                     obj_6_only@assays$RNA$data['slc17a6b',]>0 |  
                                                     obj_6_only@assays$RNA$data['slc17a7a',]>0)), 'Mixed',obj_6_only$primary_neurotransmitter)

obj_6_only$primary_neurotransmitter <- ifelse(is.na(obj_6_only$primary_neurotransmitter), 'Neither', obj_6_only$primary_neurotransmitter)

table(obj_6_only$primary_neurotransmitter)

DimPlot(obj_6_only, label = F, reduction = 'harmony_wnn.umap', group.by ='primary_neurotransmitter')
DimPlot(obj_6_only, label = T, reduction = 'harmony_wnn.umap', group.by ='sub.cluster')

#plot
cluster_prim_type <- table(obj_6_only$primary_neurotransmitter, obj_6_only$sub.cluster)%>%as.data.frame()
cluster_prim_type$prop[cluster_prim_type$Var2 == '6_0'] = cluster_prim_type$Freq[cluster_prim_type$Var2 == '6_0']/ nrow(obj_6_only@meta.data[obj_6_only@meta.data$sub == '6_0',])
cluster_prim_type$prop[cluster_prim_type$Var2 == '6_1'] = cluster_prim_type$Freq[cluster_prim_type$Var2 == '6_1']/ nrow(obj_6_only@meta.data[obj_6_only@meta.data$sub == '6_1',])
cluster_prim_type$prop[cluster_prim_type$Var2 == '6_2'] = cluster_prim_type$Freq[cluster_prim_type$Var2 == '6_2']/ nrow(obj_6_only@meta.data[obj_6_only@meta.data$sub == '6_2',])
cluster_prim_type$prop[cluster_prim_type$Var2 == '6_3'] = cluster_prim_type$Freq[cluster_prim_type$Var2 == '6_3']/ nrow(obj_6_only@meta.data[obj_6_only@meta.data$sub == '6_3',])

ggplot(cluster_prim_type, aes(x = Var2, y = prop, group = interaction(Var2, Var1), fill = Var1))+
  geom_bar(stat='identity')+
  labs(x = 'Subcluster', y = 'Proportion')+
  theme_classic()


#compare it to cluster 0
obj_subclustered_6$primary_neurotransmitter <- ifelse(((obj_subclustered_6@assays$RNA$data['LOC111588076',]>0 
                                                | obj_subclustered_6@assays$RNA$data['gad2',]>0) &
                                                 (obj_subclustered_6@assays$RNA$data['LOC111584103',]==0 & 
                                                    obj_subclustered_6@assays$RNA$data['slc17a6b',]==0 &  
                                                    obj_subclustered_6@assays$RNA$data['slc17a7a',]==0)), 'GABA',NA)

obj_subclustered_6$primary_neurotransmitter <- ifelse(((obj_subclustered_6@assays$RNA$data['LOC111588076',]==0 
                                                & obj_subclustered_6@assays$RNA$data['gad2',]==0) &
                                                 (obj_subclustered_6@assays$RNA$data['LOC111584103',]>0 | 
                                                    obj_subclustered_6@assays$RNA$data['slc17a6b',]>0 |  
                                                    obj_subclustered_6@assays$RNA$data['slc17a7a',]>0)), 'GLUT',obj_subclustered_6$primary_neurotransmitter)

obj_subclustered_6$primary_neurotransmitter <- ifelse(((obj_subclustered_6@assays$RNA$data['LOC111588076',]>0 
                                                | obj_subclustered_6@assays$RNA$data['gad2',]>0) &
                                                 (obj_subclustered_6@assays$RNA$data['LOC111584103',]>0 | 
                                                    obj_subclustered_6@assays$RNA$data['slc17a6b',]>0 |  
                                                    obj_subclustered_6@assays$RNA$data['slc17a7a',]>0)), 'Mixed',obj_subclustered_6$primary_neurotransmitter)

cluster_prim_type_o <- table(obj_subclustered_6$primary_neurotransmitter, obj_subclustered_6$sub.cluster)%>%as.data.frame()
cluster_prim_type_o$prop = cluster_prim_type_o$Freq/ nrow(obj_subclustered_6@meta.data)

ggplot(cluster_prim_type_o, aes(x = Var2, y = prop, group = interaction(Var2, Var1), fill = Var1))+
  geom_bar(stat='identity', position = 'fill')+
  labs(x = 'Subcluster', y = 'Proportion')+
  theme_classic()

#Not sure if this is a really interesting finding or not, it seems like its not 
#super common to have both characteristics but not some exceptional finding either

# i want to compare the proportion of mixed cells in these subclusters 
# to other clusters, so I think fishers exact test for [mixed, unmixed]

cells_in_clusters <- obj_subclustered_6@meta.data%>%
  group_by()%>%
  summarize(mixed_cells = sum(primary_neurotransmitter=='Mixed', na.rm = T),
            n_cells = n())%>%
  mutate(non_mixed = n_cells-mixed_cells)%>%
  subset(sub.cluster %notin% c('11', # mg,
                               '1', #rg
                               '2', #og
                               '26',#dg
                               '13', #opc
                               '15', #epe
                               '20' #leuko
                               ))
mixed_test <- data.frame()
for(i in 0:3){
  cluster = paste0('6_',i)
fisher_mixed_matrix <- matrix(c(
  sum(cells_in_clusters$mixed_cells[cells_in_clusters$sub.cluster ==cluster]),
  sum(cells_in_clusters$non_mixed[cells_in_clusters$sub.cluster ==cluster]),
  sum(cells_in_clusters$mixed_cells[cells_in_clusters$sub.cluster !=cluster]),
  sum(cells_in_clusters$non_mixed[cells_in_clusters$sub.cluster !=cluster])),
  2,
  2,
  byrow = T)

fish = fisher.test(fisher_mixed_matrix)
newd = data.frame(cluster = cluster,
                  p = fish$p.value)
mixed_test = rbind(mixed_test, newd)

}
mixed_test$p_adj = p.adjust(mixed_test$p, 'fdr',4)
mixed_test$p_adj = round(mixed_test$p_adj, 10)

### 5) IEG expression ####

# i guess I can do zack version first
obj_subclustered_6$ieg <- ifelse(obj_subclustered_6@assays$RNA$data['LOC111583367',] >0, 'ieg', NA)
obj_subclustered_6$ieg <- ifelse(obj_subclustered_6@assays$RNA$data['egr1',] >0, 'ieg', obj_subclustered_6$ieg)
obj_subclustered_6$ieg <- ifelse(obj_subclustered_6@assays$RNA$data['npas4a',] >0, 'ieg', obj_subclustered_6$ieg)
obj_subclustered_6$ieg <- ifelse(is.na(obj_subclustered_6$ieg), 'no_ieg', obj_subclustered_6$ieg)

DimPlot(obj_subclustered_6, split.by = 'ieg')

### filter out non-neuronal
neuronal_only <- subset(obj_subclustered_6, sub.cluster %notin%
                          c('11', # mg,
                            '1', #rg
                            '2', #og
                            '26',#dg
                            '13', #opc
                            '15', #epe
                            '20' #leuko
                          ))

mixed_test <- data.frame()
DimPlot(neuronal_only, split.by = 'ieg')

ieg_data_cluster_level <- table(neuronal_only@meta.data$ieg, neuronal_only@meta.data$sub.cluster)%>%
  as.data.frame()%>%
  pivot_wider(names_from = 'Var1', values_from = 'Freq')%>%
  mutate(prop_pos = ieg/(no_ieg+ieg))

ggplot(ieg_data_cluster_level, aes(x = Var2, y = prop_pos))+
  geom_point()
# 6_1 and 6_2 off the charts as well as 22 which is a potential POA thing

## anyway I now need to do the findmarkers thing

#Find genes differentially expressed between IEG and IEG- clusters
ieg_marker_all <-data.frame()
for(i in unique(neuronal_only$sub.cluster)){
  print(i)
  temp_obj <- subset(neuronal_only, sub.cluster==i)
  Idents(temp_obj) = 'ieg'
  ieg_markers <- FindAllMarkers(temp_obj, 
                                assay = 'RNA', 
                                group.by = 'ieg',
                                logfc.threshold = 0,
                                min.pct = 1/127 #the smallest cluster
  ) #204 cells in smallest cluster
  ieg_markers$sub.cluster =i
  ieg_marker_all <- rbind(ieg_marker_all, ieg_markers)
}
#extract significant ieg like genes
ieg_marker_all_signif <- ieg_marker_all[ieg_marker_all$p_val_adj<0.05 &ieg_marker_all$cluster=='ieg' ,]

#find genes significant in over half of clusters
ieg_like_genes_data <- table(ieg_marker_all_signif$gene,
                             ieg_marker_all_signif$sub.cluster)%>%
  as.data.frame()%>%
  group_by(Var1)%>% #cluster
  summarize(n_clusters = sum(Freq))%>%
  subset(n_clusters > length(unique(neuronal_only$sub.cluster))/2
  ) 

#ieg like genes
ieg_like_genes <- ieg_like_genes_data$Var1%>%droplevels()

neuronal_only$ieg2 <- NA
#Classify cells that express any of these iegs as ieg
for(i in ieg_like_genes){
  neuronal_only$ieg2 <-  ifelse(neuronal_only@assays$RNA$data[i,] >0, 'ieg', neuronal_only$ieg2)
  
}
neuronal_only$ieg2<-  ifelse(is.na(neuronal_only$ieg2), 'non_ieg', neuronal_only$ieg2)

ieg_data_cluster_level2 <- table(neuronal_only@meta.data$ieg2, neuronal_only@meta.data$sub.cluster)%>%
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
  dplyr::select(c(individual, Status, sub.cluster, ieg2))%>%
  group_by(Status, individual, sub.cluster)%>%
  summarize(ieg_pos = sum(ieg2 =='ieg'),
            ieg_neg = sum(ieg2 == 'non_ieg')
  )%>%
  mutate(prop_pos = ieg_pos/(ieg_neg+ieg_pos))%>%
  subset(Status %in% c('M','D','F'))%>%
  na.omit()%>%
  subset(prop_pos != Inf)

ieg_data_individual_level2$Status <- factor(ieg_data_individual_level2$Status, levels = c('M','D','F'))

ggplot(ieg_data_individual_level2, aes(x = as.factor(sub.cluster), y = prop_pos, shape = Status, color = Status))+
  geom_boxplot(alpha = 0, outlier.shape = NA)+
  geom_point( position = position_dodge(1), size =1.25, color = 'black')+
  labs(x = 'Cluster', y = 'Proportion IEG+')+
  theme_minimal()+
  scale_shape_manual(values = c(1,2,3))


#How is IEG correlated to cyto?

### 6) DEGs ####
#> I did the deg analysis in another file in this folder, this is just the interesting part
deg_data <- read.csv("Subclustering/degs_6_defined_08_09_2025.csv")

#clown_go(deg_data$gene[deg_data$cluster=='6_0'])%>%dotplot()
#none somehow

clown_go(deg_data$gene[deg_data$cluster=='6_1'])%>%dotplot()

clown_go(deg_data$gene[deg_data$cluster=='6_2'])%>%dotplot()

clown_go(deg_data$gene[deg_data$cluster=='6_3'])%>%dotplot()

#im not convinced by most of these, I think just plot expression 

## plotting
plot_mean_expression(individual_mean_expression(obj_6_only, 'tacr3a'))+
  labs(title = 'tacr3a')
# holy jebus its upregulated in dominants, shocking its not significantly diff
#between males and dominants probably due to small cluster size

plot_mean_expression(individual_mean_expression(obj_6_only, 'pgr'))+
  labs(title = 'pgr')
# way upreguated in females, probaby not different between dominants and males

plot_mean_expression(individual_mean_expression(obj_6_only, 'LOC111562889'))+
  labs(title = 'rbfox1-homolog, LOC111562889')
  
