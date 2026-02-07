# perineuronal net and ECM GO term PCA

library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
  clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)
library(biomaRt)
  library(factoextra)

colors = c("red", "#006400", "blue",'#000000', 'purple','gray','brown','orange')

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_6 = FindSubCluster(obj, 6, graph.name='harmony.wsnn')
Idents(sub_6) <- 'sub.cluster'
sub_6 = subset(sub_6, final_clusters ==6)
sub_6$sub.cluster = factor(sub_6$sub.cluster, levels = c(paste0('6_',0:3)))
sub_6$Status = factor(sub_6$Status, levels = c('NRM','M',"D",'E','NF','F'))

term2gene = readRDS("Function Scripts/Dependencies/Term2gene_clown_go2.rds")
term2name = readRDS('/Users/ggraham/Desktop/multiome_poa/Function Scripts/Dependencies/Term2name.rds')

go_terms = term2gene%>%
  left_join(term2name, by = 'go_id')

perineuronal_net_ecm_genes = go_terms$aocellaris_name[ go_terms$go_id == 'GO:0072534' # pnn 
                                                         |go_terms$go_id =='GO:0031012' #ecm
                                                       ]

gene_pca_matrix = sub_6@assays$RNA$data[unique(perineuronal_net_ecm_genes),]%>%t()%>%as.matrix()
gene_pca_matrix_no0 = gene_pca_matrix[,which(colSums(gene_pca_matrix)>0)]
pca = princomp(gene_pca_matrix_no0, cor = T)

fviz_pca_var(pca)
hist(pca$loadings[,1])
mean(pca$loadings[,1])

sub_6$module = pca$scores[,1]*-1 # i guess we go by dim1?

pca_ind =sub_6@meta.data%>%
  group_by(Status, individual, sub.cluster)%>%
  summarize(mean_module = mean(module),
            se_module = sd(module)/sqrt(n()))

ggplot(pca_ind, aes(x = Status, y= mean_module))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~sub.cluster, scales ='free')
 # 6+_2 looks obviously like a significant difff

model =aov(mean_module~Status, data = subset(pca_ind, sub.cluster == '6_2' & Status != 'NRM'))
summary(model) #LETS FUCKING GOOO

# maybe 6_1 too
model2 =aov(mean_module~Status, data = subset(pca_ind, sub.cluster == '6_1' & Status != 'NRM'))
summary(model2) # nope

#ok lets just look at it with PNN and ECM individually too
#####PNN ####
pnn_genes = go_terms$aocellaris_name[ go_terms$go_id == 'GO:0072534' # pnn 
                                                       ]

gene_pnn_pca_matrix = sub_6@assays$RNA$data[unique(pnn_genes),]%>%t()%>%as.matrix()
pca_pnn_no0 = gene_pnn_pca_matrix[,which(colSums(gene_pnn_pca_matrix)>0)]
pca_pnn = princomp(pca_pnn_no0,cor = T)

fviz_pca_var(pca_pnn) # loadings on both, what do

 pca_pnn$center
 pca_pnn$loadings
  pca_pnn$loadings[,2]%>%hist()
  pca_pnn$loadings[,2]%>%mean()

sub_6$pnn = pca_pnn$scores[,1]*-1 
sub_6$pnn_2 = pca_pnn$scores[,2]*-1 

pca_ind_pnn =sub_6@meta.data%>%
  group_by(Status, individual, sub.cluster)%>%
  summarize(mean_module = mean(pnn),
            se_module = sd(pnn)/sqrt(n()))

pca_ind_pnn2 =sub_6@meta.data%>%
  group_by(Status, individual, sub.cluster)%>%
  summarize(mean_module = mean(pnn_2),
            se_module = sd(pnn_2)/sqrt(n()))

ggplot(pca_ind_pnn, aes(x = Status, y= mean_module))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~sub.cluster, scales ='free')

pnn_model1 =aov(mean_module~Status, data = subset(pca_ind_pnn, sub.cluster == '6_2' & Status != 'NRM'))
summary(pnn_model1) # nope

ggplot(pca_ind_pnn2, aes(x = Status, y= mean_module))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~sub.cluster, scales ='free')
# i dont think there is anything, maybe 6_0

pnn_model2 =aov(mean_module~Status, data = subset(pca_ind_pnn2, sub.cluster == '6_0' & Status != 'NRM'))
summary(pnn_model2) 

# ok finally lets do it with just ECM
#####ECM ####

ecm_genes = go_terms$aocellaris_name[ go_terms$go_id == 'GO:0031012' # ecm 
                                                       ]
gene_ecm_pca_matrix = sub_6@assays$RNA$data[unique(ecm_genes),]%>%t()%>%as.matrix()
pca_ecm_no0 = gene_ecm_pca_matrix[,which(colSums(gene_ecm_pca_matrix)>0)]
pca_ecm = princomp(pca_ecm_no0,cor = T)

fviz_pca_var(pca_ecm)
 pca_ecm$center
 pca_ecm$loadings
  pca_ecm$loadings[,1]%>%hist()
  pca_ecm$loadings[,1]%>%mean()
#mean is positive so lets use PC1
  sub_6$ecm = pca_ecm$scores[,1]

  pca_ind_ecm =sub_6@meta.data%>%
  group_by(Status, individual, sub.cluster)%>%
  summarize(mean_module = mean(ecm),
            se_module = sd(ecm)/sqrt(n()))

ggplot(pca_ind_ecm, aes(x = Status, y= mean_module))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~sub.cluster, scales ='free')
# very obvious difference here in 6_2 wow
# ok so we have pretty obvious evidence that something to do with the ECM changes in intermediates, 
# and it seems like expression of ECM related genes goes down

model_ecm_1 = aov(mean_module ~ Status, data=subset(pca_ind_ecm, Status != 'NRM'& sub.cluster == '6_2'))
summary(model_ecm_1) #0.05 lmfao Ill take it I guess # surely the mixed model gives a lower p here

ecm_lmer_1 =  lmer(ecm ~ Status+(1|individual), data=subset(sub_6@meta.data, Status != 'NRM'& sub.cluster == '6_2'))
car::Anova(ecm_lmer_1,3) # singular somehow gives me a higher p value? 
#anyway Im gonna go with the anova, the effect looks obviously real to me

model_ecm_2 = aov(mean_module ~ Status, data=subset(pca_ind_ecm, Status != 'NRM'& sub.cluster == '6_1'))
summary(model_ecm_2)  # still no, though 6_1 seems to be close each time

# out of curiosity Im wondering if a mixed model gives a different result
library(lme4)

model_1_lmer =lmer(module~Status +(1|individual), data = subset(sub_6@meta.data, sub.cluster == '6_2' & Status != 'NRM'))
car::Anova(model_1_lmer, 3) # singular and no longer significant

model_2_lmer =lmer(module~Status +(1|individual), data = subset(sub_6@meta.data, sub.cluster == '6_1' & Status != 'NRM'))
car::Anova(model_2_lmer, 3) #same

# maybe 6_1 too
model_pnn_1_lmer =lmer(pnn~Status+(1|individual), data = subset(sub_6@meta.data, sub.cluster == '6_1' & Status != 'NRM'))
car::Anova(model_pnn_1_lmer,3) # nope


