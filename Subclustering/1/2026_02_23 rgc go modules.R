# generic module function where I input a GO term and it figures it out
library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)

go_module = function(term){
term2gene = readRDS("Function Scripts/Dependencies/Term2gene_clown_go2.rds")
term2name = readRDS('/Users/ggraham/Desktop/multiome_poa/Function Scripts/Dependencies/Term2name.rds')

go_terms = term2gene%>%
left_join(term2name, by = 'go_id')

term_genes = go_terms$aocellaris_name[ go_terms$go_id == term]
term_genes_in_obj = term_genes[term_genes %in% rownames(sub_1)]
print(paste0(length(term_genes_in_obj),' genes found'))

gene_pca_matrix = sub_1@assays$RNA$data[unique(term_genes_in_obj),]%>%t()%>%as.matrix()
gene_pca_matrix_no0 = gene_pca_matrix[,which(colSums(gene_pca_matrix)>0)]
pca = princomp(gene_pca_matrix_no0, cor = T)

  var_explained = pca$sdev^2
  max_var_pc = which.max(var_explained)
  print(paste0("PC", max_var_pc, " explains ", 
               round(var_explained[max_var_pc]/sum(var_explained)*100, 2), 
               "% of variance"))

if(mean(pca$loadings[,max_var_pc])>0){
  scores = pca$scores[,max_var_pc]
}else{
  scores=pca$scores[,max_var_pc]*-1
}

dupe = sub_1

dupe$module = scores

pca_ind =dupe@meta.data%>%
  group_by(Status, individual, sub.cluster)%>%
  summarize(mean_module = mean(module),
            se_module = sd(module)/sqrt(n()))

assign(paste0('module_', term), pca_ind, envir = .GlobalEnv)
p = ggplot(pca_ind, aes(x = Status, y= mean_module))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~sub.cluster, scales ='free')
return(p)
}


obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_1 = FindSubCluster(obj,1, graph.name='harmony.wsnn')
Idents(sub_1) <- 'sub.cluster'
sub_1 = subset(sub_1, final_clusters ==1)
sub_1$Status = factor(sub_1$Status, levels = c('NRM','M',"D",'E','NF','F'))


#go_module('GO:0031012')

#go_module('GO:0001508') #action potential

go_module('GO:0022008')# neurogenesis, 1_0, 1_3?, maybe 1_4

neurogenesis_module = data.frame()
for(i in paste0('1_', 0:5)){
  model_data = subset(`module_GO:0022008`, Status != 'NRM' & sub.cluster == i)
  model = aov(mean_module~Status, data =model_data)
 sum= summary(model)[[1]]
  newd = data.frame(cluster =i,
                    p = sum$`Pr(>F)`[1])
  neurogenesis_module = rbind(neurogenesis_module, newd)
  
}
# 1_0


go_module('GO:0007420') # brain development, 1_0

brain_dev_module = data.frame()
for(i in paste0('1_', 0:5)){
  model_data = subset(`module_GO:0007420`, Status != 'NRM' & sub.cluster == i)
  model = aov(mean_module~Status, data =model_data)
 sum= summary(model)[[1]]
  newd = data.frame(cluster =i,
                    p = sum$`Pr(>F)`[1])
  brain_dev_module = rbind(brain_dev_module, newd)
  
} # no

go_module('GO:0010975')# neuron projection development

neuron_proj_module = data.frame()
for(i in paste0('1_', 0:5)){
  model_data = subset(`module_GO:0010975`, Status != 'NRM' & sub.cluster == i)
  model = aov(mean_module~Status, data =model_data)
 sum= summary(model)[[1]]
  newd = data.frame(cluster =i,
                    p = sum$`Pr(>F)`[1])
  neuron_proj_module = rbind(neuron_proj_module, newd)
  
}# no

go_module('GO:0042551') #neuron maturation, 1_4

neuron_mat_module = data.frame()
for(i in paste0('1_', 0:5)){
  model_data = subset(`module_GO:0042551`, Status != 'NRM' & sub.cluster == i)
  model = aov(mean_module~Status, data =model_data)
 sum= summary(model)[[1]]
  newd = data.frame(cluster =i,
                    p = sum$`Pr(>F)`[1])
  neuron_mat_module = rbind(neuron_mat_module, newd)
  
} # no

go_module('GO:0016081') # synaptic vesicle docking, maybe

vesicule_module = data.frame()
for(i in paste0('1_', 0:5)){
  model_data = subset(`module_GO:0016081`, Status != 'NRM' & sub.cluster == i)
  model = aov(mean_module~Status, data =model_data)
 sum= summary(model)[[1]]
  newd = data.frame(cluster =i,
                    p = sum$`Pr(>F)`[1])
  vesicule_module = rbind(vesicule_module, newd)
  
} # no

#go_module('GO:0007218') #neuropeptide signaling

go_module('GO:0034056') #ERE binding, 1_1, 1_3?, 1_2

ere_module = data.frame()
for(i in paste0('1_', 0:5)){
  model_data = subset(`module_GO:0034056`, Status != 'NRM' & sub.cluster == i)
  model = aov(mean_module~Status, data =model_data)
 sum= summary(model)[[1]]
  newd = data.frame(cluster =i,
                    p = sum$`Pr(>F)`[1])
  ere_module = rbind(ere_module, newd)
  
} # 1_1, 1_3

#go_module('GO:0030521') # androgen receptor
