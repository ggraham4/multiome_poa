### GJG 2025 08 23
### I noticed dopamine receptor in cluster 6 was a deg, what's goin on w that

#libs
{
  library(Seurat)
  library(tidyverse)
  library(biomaRt)
  library(dplyr)
  library(lme4)
  library(car)
  library(emmeans)
  library(ggplot2)
  library(hdWGCNA)
  }

#functions
{
  `%notin%` <- Negate(`%in%`)
  ensembl <- useEnsembl(biomart = "genes", 
                        dataset = "aocellaris_gene_ensembl")
  a <- listAttributes(ensembl)
  biomart_basic <-
    getBM(
      mart = ensembl, #working mart 
      attributes = c("entrezgene_accession",
                     'entrezgene_description',
                     'external_gene_name'))
  mecp = readRDS("Functions/mean_expression_cluster_plot.rds")
  namer = function(gene){
    gene_name = biomart_basic$entrezgene_description[biomart_basic$entrezgene_accession==gene]
    return(gene_name[1])
    
  }
  namer2 = function(gene){
    gene_name = biomart_basic$external_gene_name[biomart_basic$entrezgene_accession==gene]
    return(gene_name[1])
    
  }
  gene_prop_data = function(object, gene){
    if(!require(tidyverse)){library(tidyverse)}
    
    #classify cells as positive or negative for the gene
    object$pos_neg = ifelse(object@assays$RNA$data[gene,]>0, TRUE, FALSE)
    
    #make a table
    model_data = object@meta.data%>%
      group_by(individual, Status)%>%
      subset(Status %in% c('M','D', 'F'))%>%
      summarize(n_pos = sum(pos_neg==TRUE),
                n_neg = sum(pos_neg==FALSE))
    
    return(model_data)
  }
  plot_module <-  function(module){
    #ur supposed to use ME not hME
    me_subset=MEs%>%
      group_by(individual, Status)%>%
      summarize(mean_module = mean(.data[[module]]),
                se_module = sd(.data[[module]])/sqrt(n()))
    
    me_subset$Status = factor(me_subset$Status, levels = c('NRM','M',"D",'E','NF','F'))
    
    plot <- ggplot(me_subset, aes(x = Status, y = mean_module, color = Status))+
      geom_boxplot(aes(group = Status, fill = Status), alpha = 0.25, outlier.shape = NA)+
      geom_pointrange(aes(x = Status, y = mean_module,
                          ymin = mean_module - se_module,
                          ymax = mean_module +se_module),
                      position = position_jitterdodge(3))+
      theme_classic()+
      labs(x  ='Status', y = module)+
      theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
    return(plot)
    
  }
  
  plot_module_cluster <-  function(module){
    #ur supposed to use ME not hME
    me_subset=MEs%>%
      group_by(individual, sub.cluster)%>%
      summarize(mean_module = mean(.data[[module]]),
                se_module = sd(.data[[module]])/sqrt(n()))
    
    plot <- ggplot(me_subset, aes(x = sub.cluster, y = mean_module, color = sub.cluster))+
      geom_boxplot(alpha = 0.25, outlier.shape = NA)+
      geom_pointrange(aes(x = sub.cluster, y = mean_module,
                          ymin = mean_module - se_module,
                          ymax = mean_module +se_module),
                      position = position_jitterdodge(1))+
      theme_classic()+
      labs(x  ='sub.cluster', y = module)+
      theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
    return(plot)
    
  }

  
}

#read in data
obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
#wgcna
obj_wgcna <- readRDS( 'A:/WGCNA_6/wgcna_seurat_obj_08_12_2025.rds')
hub_df = read.csv('A:/WGCNA_6/hubs_6_subclusters_08_12_2025.csv')
modules = read.csv('A:/WGCNA_6/modules_6_subclusters_08_12_2025.csv')
MEs <- GetMEs(obj_wgcna, harmonized=FALSE)

# how does drd3 expression change ####
# it increases as sex change goes on
mecp(obj, 'drd3',6, 'res0.8_50nn_40PC_45LSI')

#not a common cell type
FeaturePlot(obj, 'drd3', reduction = 'harmony_wnn.umap')

#even in cluster 6
Idents(obj) ='res0.8_50nn_40PC_45LSI'

FeaturePlot(subset(x=obj,idents = 6) , 'drd3', reduction = 'harmony_wnn.umap')

## what do drd3+ cluster 6 cells express? ####
obj_6 = subset(obj,res0.8_50nn_40PC_45LSI==6)

#mark cells as + or neg
obj_6$drd3 = ifelse(obj_6@assays$RNA$data['drd3',]>0, 'drd3', 'not')

#find markers
markers = FindAllMarkers(obj_6, group.by = 'drd3')
#not much, the second strongest marker is only in like 18% of them, I wonder
markers_drd3 = subset(markers, cluster == 'drd3' & p_val_adj <0.01)
# how it would be different if I did it without the 6 subcluster
#another time maybe

### lets name these genes
markers_drd3$gene_name= sapply(FUN= namer, X=markers_drd3$gene)
markers_drd3$gene_name2= sapply(FUN= namer2, X=markers_drd3$gene)

markers_neg = subset(markers, cluster == 'not' & p_val_adj <0.01)

### do they differ in proportion ####
# at a basic scale, yes
dop_table = table(obj_6$Status, obj_6$drd3)%>%as.data.frame.matrix()
dop_table$prop = dop_table$drd3/(dop_table$drd3+dop_table$not)

#what about in a glmer
drd3_prop_data = gene_prop_data(obj_6, 'drd3')

drd3_model = glmer(cbind(drd3_prop_data$n_pos,drd3_prop_data$n_neg)~Status + (1|individual), 
                   data = drd3_prop_data,
                   family = binomial())
drd3_prop_data$prop = drd3_prop_data$n_pos/(drd3_prop_data$n_pos+drd3_prop_data$n_neg)
drd3_prop_data$sum = (drd3_prop_data$n_pos+drd3_prop_data$n_neg)
drd3_prop_data$se = sqrt((drd3_prop_data$sum*drd3_prop_data$prop)*(1-drd3_prop_data$prop))

drd3_prop_data$Status = factor(drd3_prop_data$Status, levels = c('M','D','F'))

#plot
ggplot(drd3_prop_data, aes(x = Status, y = prop, color = Status))+
  geom_jitter(position = position_jitterdodge(0.5), size =3)+
  geom_boxplot(alpha = 0.25, aes(fill = Status))
        
#fold increase in females          
drd3_prop_data%>%
  group_by(Status)%>%
  summarise(n_pos = sum(n_pos),
            n_total = sum(n_pos+n_neg),
            prop = n_pos/n_total)
0.0458/0.0167

### is drd3 part of any interesting modules in the wgcna?? ####
#it is not

#what other interesting questions can I ask?
#which subclusters express 6 ####
obj_subclustered_6 <- FindSubCluster(obj, 6, 'harmony.wsnn')

plot = mecp(obj_subclustered_6, 'drd3',paste0('6_','0'),'sub.cluster')
for(i in 1:3){
plot = plot+mecp(obj_subclustered_6, 'drd3',paste0('6_',i),'sub.cluster')
}
plot

#what else is interesting about 6_0?

#it expresses ar
gene = 'ar'
plot = mecp(obj_subclustered_6, gene,paste0('6_','0'),'sub.cluster')
for(i in 1:3){
  plot = plot+mecp(obj_subclustered_6, gene,paste0('6_',i),'sub.cluster')
}
plot

#esr2b (though not as specifically as ar, also wait why does 6_3 have that pattern?)
gene = 'esr2b'
plot = mecp(obj_subclustered_6, gene,paste0('6_','0'),'sub.cluster')
for(i in 1:3){
  plot = plot+mecp(obj_subclustered_6, gene,paste0('6_',i),'sub.cluster')
}
plot

#pgr
gene = 'pgr'
plot = mecp(obj_subclustered_6, gene,paste0('6_','0'),'sub.cluster')
for(i in 1:3){
  plot = plot+mecp(obj_subclustered_6, gene,paste0('6_',i),'sub.cluster')
}
plot

#th
gene = 'th'
plot = mecp(obj_subclustered_6, gene,paste0('6_','0'),'sub.cluster')
for(i in 1:3){
  plot = plot+mecp(obj_subclustered_6, gene,paste0('6_',i),'sub.cluster')
}
plot

#gnrh1
gene = 'LOC111571064'
plot = mecp(obj_subclustered_6, gene,paste0('6_','0'),'sub.cluster')
for(i in 1:3){
  plot = plot+mecp(obj_subclustered_6, gene,paste0('6_',i),'sub.cluster')
}
plot


#it
gene = 'oxt'
plot = mecp(obj_subclustered_6, gene,paste0('6_','0'),'sub.cluster')
for(i in 1:3){
  plot = plot+mecp(obj_subclustered_6, gene,paste0('6_',i),'sub.cluster')
}
plot

#avt
gene = 'avp'
plot = mecp(obj_subclustered_6, gene,paste0('6_','0'),'sub.cluster')
for(i in 1:3){
  plot = plot+mecp(obj_subclustered_6, gene,paste0('6_',i),'sub.cluster')
}
plot

#cckbrb
gene = 'cckbrb'
plot = mecp(obj_subclustered_6, gene,paste0('6_','0'),'sub.cluster')
for(i in 1:3){
  plot = plot+mecp(obj_subclustered_6, gene,paste0('6_',i),'sub.cluster')
}
plot

#cckb
gene = 'cckb'
plot = mecp(obj_subclustered_6, gene,paste0('6_','0'),'sub.cluster')
for(i in 1:3){
  plot = plot+mecp(obj_subclustered_6, gene,paste0('6_',i),'sub.cluster')
}
plot

#th2
gene = 'th2'
plot = mecp(obj_subclustered_6, gene,paste0('6_','0'),'sub.cluster')
for(i in 1:3){
  plot = plot+mecp(obj_subclustered_6, gene,paste0('6_',i),'sub.cluster')
}
plot


### do DA neurons have different IEG signature to non DA? ####
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


# what about dopamine?
neuronal_only$drd3 <- ifelse(neuronal_only@assays$RNA$data['drd3',] >0, T, F)

ieg_data_cluster_level2 <- table(neuronal_only@meta.data$ieg2, neuronal_only@meta.data$drd3)%>%
  as.data.frame()%>%
  pivot_wider(names_from = 'Var1', values_from = 'Freq')%>%
  mutate(prop_pos = ieg/(non_ieg+ieg))%>%
  na.omit()


fisher_mat = cbind(ieg_data_cluster_level2$ieg, ieg_data_cluster_level2$non_ieg)
test =fisher.test(fisher_mat)
test
test$p.value

#glmer to be sure
glmer_dat_da = neuronal_only@meta.data%>%
  subset(sub.cluster %in%c('6_1','6_0','6_2','6_3'))%>%
  group_by(Status, individual,drd3)%>%
  summarize(n_pos = sum(ieg == 'ieg'),
            n_neg = sum(ieg =='no_ieg'))

da_ieg_model = glmer(cbind(glmer_dat_da$n_pos, glmer_dat_da$n_neg)~drd3*Status+(1|individual),
                     data = glmer_dat_da, family = binomial())
car::Anova(da_ieg_model, 3)

glmer_dat_da$Status = factor(glmer_dat_da$Status, levels = c('NRM',
                                                      'M',
                                                      'D',
                                                      'E',
                                                      'NF',
                                                      'F'))
ggplot(glmer_dat_da, aes(x = Status, y = n_pos/(n_neg+n_pos), shape = drd3, color = Status, fill = Status))+
  geom_boxplot(alpha = 0.2)+
  geom_point(position = position_jitterdodge(0.5), size =3)+
  labs(y = 'Proporton IEG+')


da_ieg_model2 = glmer(cbind(glmer_dat_da$n_pos, glmer_dat_da$n_neg)~drd3+(1|individual),
                     data = glmer_dat_da, family = binomial())
car::Anova(da_ieg_model2, 3)
