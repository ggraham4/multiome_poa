library(Seurat)
library(tidyverse)
library(lme4)
library(emmeans)
clown_go = readRDS('Functions/clown_go2')
library(ggplot2)
library(CytoTRACE)

# read in data
obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds") 
Idents(obj) = 'res0.8_50nn_40PC_45LSI'
# subset to 9 and only m, d, f
obj_9 = subset(obj, res0.8_50nn_40PC_45LSI==9 &
                 Status %in% c('M','D','F'))

## define a function to test if genes are differentially present
# use fixed logistic regression 

#meta data
meta_perm = obj_9@meta.data

#expression matrix
exp = obj_9@assays$RNA$data%>%as.matrix()
exp_bin = exp>0

genes = rownames(exp_bin[rowSums(exp_bin)>0,])

gene_presence = function(gene){
  print(gene)
  meta_temp = meta_perm
  meta_temp$expression = exp_bin[gene,]  
  
  meta_analysis = meta_temp%>%
    group_by(individual, Status)%>%
    summarize(successes = sum(expression),
              failures = n()-sum(expression))
  
  model_mat = cbind(meta_analysis$successes, meta_analysis$failures)
  
  model = glm(model_mat~ Status,
              family = binomial('logit'),
              data = meta_analysis)
  
  av_ = anova(model, test = 'Chisq')
  
  pairs = pairs(emmeans(model, 'Status'), adjust = 'none')%>%as.data.frame()
  
  newd = data.frame(
    gene = gene,
    av_p = av_$`Pr(>Chi)`[2],
    d_f_p = pairs$p.value[1],
    d_m_p = pairs$p.value[2],
    f_m_p = pairs$p.value[3],
    d_f_est = pairs$estimate[1],
    d_m_est = pairs$estimate[2],
    f_m_est = pairs$estimate[3]
  )  
  
}

run = lapply(X= genes, FUN =gene_presence)
run_bound = do.call(rbind, run)
run_bound$q = p.adjust(run_bound$av_p, 'fdr', nrow(run_bound))
run_bound$signif = ifelse(run_bound$q<0.05, '*', 'not')
# only one signif, mettl21a
run_bound$signif_0.1 = ifelse(run_bound$q<0.1, '*', 'not')

clown_go(run_bound$gene[run_bound$signif_0.1=='*'])%>%dotplot()

FeaturePlot(obj_9, 'mettl21a', split.by = 'Status')



# mettl21a
meta_temp = meta_perm
  meta_temp$expression = exp_bin['mettl21a',]  
  
  meta_analysis = meta_temp%>%
    group_by(individual, Status)%>%
    summarize(successes = sum(expression),
              failures = n()-sum(expression))
  meta_analysis$Status = factor(meta_analysis$Status, levels = c('M','D','F'))
  
  ggplot(meta_analysis, aes(x = Status, y = successes/(successes+failures)))+
    geom_boxplot()+
    geom_point()
  
  #Protein-lysine methyltransferase that selectively trimethylates residues in heat shock protein 70 (HSP70) family members.

    Idents(obj_9)='res0.8_50nn_40PC_45LSI'
  sub_9 = FindSubCluster(obj_9, 9, 'harmony.wsnn')
DimPlot(sub_9, group.by = 'sub.cluster', reduction = 'harmony_wnn.umap')

FeaturePlot(sub_9, 'elavl3', reductio = 'harmony_wnn.umap')
Idents(sub_9) = 'sub.cluster'
Marks = FindAllMarkers(sub_9, only.pos = T)
marks_9_2 = FindMarkers(sub_9, '9_2', only.pos = T)

### proportion differences? ####
# get metadata from sub_9 object
meta_sub = sub_9@meta.data

# get unique subclusters
subclusters = unique(meta_sub$sub.cluster)

subcluster_presence = function(subclust){
  print(subclust)
  
  meta_temp = meta_sub
  meta_temp$in_subcluster = meta_temp$sub.cluster == subclust
  
  meta_analysis = meta_temp %>%
    group_by(individual, Status) %>%
    summarize(successes = sum(in_subcluster),
              failures = n() - sum(in_subcluster))
  
  model_mat = cbind(meta_analysis$successes, meta_analysis$failures)
  
  model = glm(model_mat ~ Status,
              family = binomial('logit'),
              data = meta_analysis)
  
  av_ = anova(model, test = 'Chisq')
  
  pairs = pairs(emmeans(model, 'Status'), adjust = 'none') %>% as.data.frame()
  
  newd = data.frame(
    subcluster = subclust,
    av_p = av_$`Pr(>Chi)`[2],
    d_f_p = pairs$p.value[1],
    d_m_p = pairs$p.value[2],
    f_m_p = pairs$p.value[3],
    d_f_est = pairs$estimate[1],
    d_m_est = pairs$estimate[2],
    f_m_est = pairs$estimate[3]
  )
  
}

run_sub = lapply(X = subclusters, FUN = subcluster_presence)
run_sub_bound = do.call(rbind, run_sub)

# quick look
run_sub_bound %>% arrange(av_p)

# plot proportions for significant subclusters
signif_subs = run_sub_bound$subcluster[run_sub_bound$signif == '*']

for(subclust in signif_subs){
  meta_temp = meta_sub
  meta_temp$in_subcluster = meta_temp$sub.cluster == subclust
  
  meta_analysis = meta_temp %>%
    group_by(individual, Status) %>%
    summarize(successes = sum(in_subcluster),
              failures = n() - sum(in_subcluster))
  meta_analysis$Status = factor(meta_analysis$Status, levels = c('M','D','F'))
  
  p = ggplot(meta_analysis, aes(x = Status, y = successes/(successes+failures))) +
    geom_boxplot() +
    geom_point() +
    labs(title = paste('Subcluster', subclust), y = 'Proportion of cells')
  print(p)
}

# 9 2 transiently increases in dominants

clown_go(marks_9_2[marks_9_2$p_val_adj<0.05,]%>%rownames())%>%dotplot()

clown_go(Marks$gene[Marks$cluster=='9_1'& Marks$p_val_adj<0.05])%>%dotplot()
clown_go(Marks$gene[Marks$cluster=='9_2'& Marks$p_val_adj<0.05])%>%dotplot()
clown_go(Marks$gene[Marks$cluster=='9_0'& Marks$p_val_adj<0.05])%>%dotplot()
clown_go(Marks$gene[Marks$cluster=='9_3'& Marks$p_val_adj<0.05])%>%dotplot()

# cyto
# I hypothesize 9_2 is most immature
sub_9$Cyto = CytoTRACE(sub_9@assays$RNA$data%>%as.matrix())$CytoTRACE

VlnPlot(sub_9, 'Cyto')

DotPlot(sub_9, 'esr2b')

marks_9_0 = FindMarkers(sub_9, '9_0')
marks_9_1 = FindMarkers(sub_9, '9_1')
marks_9_3 = FindMarkers(sub_9, '9_3')


FeaturePlot(sub_9, 'rbfox3a', reduction = 'harmony_wnn.umap')
FeaturePlot(sub_9, 'sox4b', reduction = 'harmony_wnn.umap')
