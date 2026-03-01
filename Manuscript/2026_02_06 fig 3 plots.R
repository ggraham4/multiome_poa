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
term_genes_in_obj = term_genes[term_genes %in% rownames(sub_6)]
print(paste0(length(term_genes_in_obj),' genes found'))

gene_pca_matrix = sub_6@assays$RNA$data[unique(term_genes_in_obj),]%>%t()%>%as.matrix()
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

dupe = sub_6

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


#colors = c("red", "#006400", "blue",'#000000', 'purple','gray','brown','orange')

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_6 = FindSubCluster(obj, 6, graph.name='harmony.wsnn')
Idents(sub_6) <- 'sub.cluster'
sub_6 = subset(sub_6, final_clusters ==6)
sub_6$sub.cluster = factor(sub_6$sub.cluster, levels = c(paste0('6_',0:3)))
sub_6$Status = factor(sub_6$Status, levels = c('NRM','M',"D",'E','NF','F'))

DimPlot(sub_6)
DotPlot(sub_6, group.by = 'sub.cluster',
        features = c('LOC111584656'), dot.min = .1)+
  coord_flip()

### lets look at NPO Markers
markers_6_subclusters = DotPlot(sub_6, group.by = 'sub.cluster',
        features = c('hmx2',
'hmx3a',
'hmx3b',
'lhx1a',
'esr2b',
'ar',
'pgr',
'LOC111562384', # ccka
'cckb',
'cckbra',

'tacr1a',
'tacr3a',


'LOC111583585', # sox2
'sox5',
'rorb',
'gal',
'npy',
'tac1',
'trhra',
'trhrb',
'prlh2',
'otpa',
'LOC111564688', #otpb
'fezf2',
'sim1a',

'sst1.1',
'avp',
'oxt',
'nts',
'trh',
'kiss1ra',
'sema6dl',
'six6a',
'roraa',
'cckbrb',
'LOC111575776',# penka
'LOC111571064', #gnrh1
'crhb'
        ),       dot.min = .05,
        scale=T,
        dot.scale = 3)+
  scale_x_discrete(labels = c('hmx2',
                      'hmx3a',
                      'hmx3b',
                      'lhx1a',
                      'esr2b',
                      'ar',
                      'pgr',
                      'ccka',
                      'cckb', # ccka
                      'cckbra',
                      'tacr1a',
'tacr3a',
                      'sox2',
                      'sox5',
                      'rorb',
                      'gal',
                      'npy',
'tac1',
                      'trhra',
                      'trhrb',
                      'prlh2',

                      'otpa',
                      'otpbl', #otpb
                      'fezf2',
                      'sim1a',

                      'sst1.1',
                      'avp',
                      'oxt',
                      'nts',
                      'trh',
'kiss1ra',

                     'sema6dl',
                      'six6a',
                      'roraa',
                     'cckbrb',
                      'penka',# penka
                      'gnrh1', #gnrh1
                      'crhb'
  ))+
  coord_flip()+
    labs(y = 'Subcluster')+
  theme(axis.title.y = element_blank())
markers_6_subclusters
#marks_1_2 =FindMarkers(sub_6, ident.1 = '6_1', ident.2 = '6_2')
  
ggsave(plot = markers_6_subclusters,
       file = "markers_6_subclusters.svg",
       device = "svg",
       units = "in",
       width = 3.5,
       height = 3.8,
       path = "Manuscript/Plots/Fig.3")

dimplot = DimPlot(sub_6)+
  theme_void()+
  theme(legend.position = "none")
dimplot
DimPlot(sub_6, label = T)

#ggsave(plot = dimplot,
#       file = "UMAP_sub6.tiff",
#       device = "tiff",
#       units = "in",
#       width = 4,
#       height = 4,
#       path = "Manuscript/Plots/Fig.3")

### DEG enrichgment
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_6 = degs$gene[degs$cluster==6]

degs_6_subclusters = DotPlot(sub_6, group.by = 'sub.cluster',
                             features = degs_6)+
  coord_flip()

degs_6_subclusters


### ECM ####

term2gene = readRDS("Function Scripts/Dependencies/Term2gene_clown_go2.rds")
term2name = readRDS('/Users/ggraham/Desktop/multiome_poa/Function Scripts/Dependencies/Term2name.rds')

go_terms = term2gene%>%
  left_join(term2name, by = 'go_id')

ecm_genes = go_terms$aocellaris_name[ go_terms$go_id == 'GO:0031012' # ecm 
                                                       ]
gene_ecm_pca_matrix = sub_6@assays$RNA$data[unique(ecm_genes),]%>%t()%>%as.matrix()
pca_ecm_no0 = gene_ecm_pca_matrix[,which(colSums(gene_ecm_pca_matrix)>0)]
pca_ecm = princomp(pca_ecm_no0,cor = T)

library(factoextra)
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

  dag =subset(pca_ind_ecm,sub.cluster =='6_2' & Status != 'NRM')
  dag$Status = as.character(dag$Status)
  dag$Phase = ifelse(dag$Status =='D', 'I', dag$Status)
  dag$Phase = ifelse(dag$Status =='E', 'LI', dag$Phase)
  
  dag$Phase = factor(dag$Phase, levels = c('M',
                                           "I",
                                           "LI",
                                           "NF",
                                           "F"))
  model_ecm_1 = aov(mean_module ~ Status, data=subset(pca_ind_ecm, Status != 'NRM'& sub.cluster == '6_2'))

  
  pair = pairs(emmeans(model_ecm_1, 'Status'), adjust = 'none')
textsize = 3

ecm_module_6_2 = ggplot(dag, aes(x = Phase, y= mean_module))+
      geom_boxplot(data = subset(dag, Phase != "NF"),
                 aes(x = Phase, y =mean_module,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
  labs(x = 'Phase', y = 'ECM PC1')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.9), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(1.1), annotation =c('p = 0.715'), color = "black",tip_length = c(0,0), textsize = textsize)+
    geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.9), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
    theme(legend.position = 'none')


#ggsave(plot = ecm_module_6_2,
#       file = "ecm_module_6_2.svg",
#       device = "svg",
#       units = "in",
#       width = 1.5,
#       height = 1.5,
#       path = "Manuscript/Plots/Fig.3")


go_module('GO:0034056') #ERE binding
ere_0 = aov(mean_module~Status, data = subset(`module_GO:0034056`,Status!= 'NRM' & sub.cluster == '6_0' ))
pair2 = pairs(emmeans(ere_0, 'Status'), adjust = 'none')

#ere_module_6_0

  dag2 =subset(`module_GO:0034056`,Status!= 'NRM' & sub.cluster == '6_0' )

  dag2$Status = as.character(dag2$Status)
  dag2$Phase = ifelse(dag2$Status =='D', 'I', dag2$Status)
  dag2$Phase = ifelse(dag2$Status =='E', 'LI', dag2$Phase)
  
  dag2$Phase = factor(dag2$Phase, levels = c('M',
                                           "I",
                                           "LI",
                                           "NF",
                                           "F"))
  
  ere_6_0= ggplot(dag2, aes(x = Phase, y= mean_module))+
      geom_boxplot(data = subset(dag2, Phase != "NF"),
                 aes(x = Phase, y =mean_module,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
  labs(x = 'Phase', y = 'ERE PC1')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(1.5), annotation =c('p = 0.194'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(1.8), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
    geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(1.5), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
    theme(legend.position = 'none')

ggsave(plot = ere_6_0,
       file = "ere_module_6_0.svg",
       device = "svg",
       units = "in",
       width = 1.5,
       height = 1.5,
       path = "Manuscript/Plots/Fig.3")



### proportion of 6 #####
sub_cells = sub_6@meta.data%>%
  group_by(sub.cluster, Status, individual)%>%
  subset(Status !='NRM')%>%
  summarize(cells = n())

total_cells = sub_6@meta.data%>%
  group_by(Status, individual)%>%
  subset(Status !='NRM')%>%
  summarize(total_cells = n())

full_data <- sub_cells%>%
  right_join(total_cells, by = 'individual')

full_data$diff = full_data$total_cells-full_data$cells

glm_output = data.frame()
for(i in 0 :3){
  subset_data = subset(full_data, sub.cluster  ==paste0('6_',i) )
  
  glm_matrix = matrix(c(subset_data$cells,
                          subset_data$diff),
                        nrow(subset_data), 2)
  
  glm_model = glm(glm_matrix~Status.x,
                      family = binomial('logit'),
                      data = subset_data)
  
  av= car::Anova(glm_model, type = 3)
  
  newd = data.frame(cluster = i,
                    av = av$`Pr(>Chisq)`[1])
  glm_output = rbind(glm_output, newd)
  
}
# 6_3, p = 0.07 but I am going to take that as significant
full_data$Status.x = as.character(full_data$Status.x )

full_data$Phase = ifelse(full_data$Status.x=='D', 'I',full_data$Status.x )
full_data$Phase = ifelse(full_data$Status.x=='E', 'LI',full_data$Phase )
full_data$Phase = factor(full_data$Phase, levels = c('M',
                                                     'I',
                                                     'LI',
                                                     'NF',
                                                     'F'))

  subset_data = subset(full_data, sub.cluster  ==paste0('6_',3) )
  
  glm_matrix = matrix(c(subset_data$cells,
                          subset_data$diff),
                        nrow(subset_data), 2)
  
  glm_model = glm(glm_matrix~Status.x,
                      family = binomial('logit'),
                      data = subset_data)
  
  av= car::Anova(glm_model, type = 3)
  
  newd = data.frame(cluster = i,
                    av = av$`Pr(>Chisq)`[1])
  glm_output = rbind(glm_output, newd)


pair3 = pairs(emmeans(glm_model, 'Status.x'), adjust = 'none')

prop_6_3 =ggplot(subset(full_data, sub.cluster=='6_3'), aes(x = Phase, y = cells/total_cells, fill = Phase))+
    geom_boxplot(data = subset(full_data, sub.cluster=='6_3'&Phase != "NF"),
                 aes(x = Phase, y =cells/total_cells ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
    scale_y_continuous(labels = scales::percent)+
    labs(x = 'Phase', y = '% of 6')+
    geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.17), annotation =c('p = 0.176'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.22), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
    geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.17), annotation =c('p = 0.226'), color = "black",tip_length = c(0,0), textsize = textsize)+
    theme(legend.position = 'none')

ggsave(plot = prop_6_3,
       file = "prop_6_3.svg",
       device = "svg",
       units = "in",
       width = 1.5,
       height = 1.5,
       path = "Manuscript/Plots/Fig.3")

######degs
degs_6 = subset(degs, cluster == 6)

p_value_annotation <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return(paste0("p = ", round(p, 3)))
  }
}
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')

deg = 'pgr'
data = degs

obj@meta.data$Phase = obj@meta.data$Status
obj@meta.data$Phase= ifelse(obj@meta.data$Phase == 'D', 'I', obj@meta.data$Phase)
obj@meta.data$Phase= ifelse(obj@meta.data$Phase == 'E', 'LI', obj@meta.data$Phase)

obj@meta.data$Phase = factor(obj@meta.data$Phase, levels = c('NRM',
                                                             'M',
                                                             'I',
                                                             'LI',
                                                             'NF',
                                                             'F'))

plotter_function = function(data, deg){
  
  deg_data = subset(data, gene == deg)
  d_f_p.value = deg_data$d_f_p.value
  d_m_p.value = deg_data$d_m_p.value
  f_m_p.value = deg_data$f_m_p.value
  
  d_f_annotation = p_value_annotation(d_f_p.value)
  d_m_annotation = p_value_annotation(d_m_p.value)
  f_m_annotation = p_value_annotation(f_m_p.value)
  
  meta.data = subset(obj@meta.data, final_clusters ==6)
  meta.data$expression = obj@assays$RNA$data[deg,obj@meta.data$final_clusters == 6]
  
  meta_grouped = meta.data%>%
    group_by(Phase, individual)%>%
    summarize(mean_expression = mean(expression))%>%
    subset(Phase != 'NRM')
  
  meta_grouped$Phase = factor(meta_grouped$Phase, levels = c(
                                                             'M',
                                                             'I',
                                                             'LI',
                                                             'NF',
                                                             'F'))

  
    r = ggplot(meta_grouped, aes(x = Phase, y = mean_expression))+
    geom_boxplot(data = subset(meta_grouped, Phase != "NF"),
                 aes(x = Phase, y =mean_expression ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
    labs(x = 'Phase', y = 'Normalized Expression')+
    geom_signif(xmin = c(1),
                xmax = c(1.9),
                y_position = c(max(meta_grouped$mean_expression*1.1)),
                annotation =c(d_m_annotation),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
          geom_signif(xmin = c(1),
                xmax = c(5),
                y_position = c(max(meta_grouped$mean_expression*1.3)),
                annotation =c(f_m_annotation),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
    geom_signif(xmin = c(2.1),
                xmax = c(5),
                y_position = c(max(meta_grouped$mean_expression*1.1)),
                annotation =c(d_f_annotation),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
      theme(legend.position = 'none', axis.title.y = element_blank())
    r
  return(r)
}

for(i in data$gene[data$cluster=='6']){
  plot = plotter_function(data, i)
  
  #ggsave(plot = plot,
  #     file = paste0('deg_', i, '.svg'),
   #    device = "svg",
    #   units = "in",
     #  width = 1.35,
      # height = 1.5,
       #path = "Manuscript/Plots/Fig.3/degs_6")

  
}

#### GnRH + and cckb+ neurons ####
sub_6$GnRH = ifelse(sub_6@assays$RNA$data['LOC111571064',]>0, T, F)

gnrh = sub_6@meta.data%>%
  group_by(individual, Status)%>%
  summarize(n_GnRH = sum(GnRH == T))

tot= sub_6@meta.data%>%
  group_by(individual, Status)%>%
  summarize(total_cells = n())

tog = gnrh%>%
  right_join(tot, by = c('individual'))

tog_6_gnrh= subset(tog, Status.x!='NRM')

gnrh_glm =glm(cbind(tog_6_gnrh$n_GnRH, tog_6_gnrh$total_cells-tog_6_gnrh$n_GnRH)~Status.x, 
              data = tog_6_gnrh,
              family = binomial('logit'))

anova(gnrh_glm, test = 'Chisq')

tog_6_gnrh$Phase = as.character(tog_6_gnrh$Status.x)
tog_6_gnrh$Phase = ifelse(tog_6_gnrh$Phase =='D','I', tog_6_gnrh$Phase )
tog_6_gnrh$Phase = ifelse(tog_6_gnrh$Phase =='E','LI', tog_6_gnrh$Phase )
tog_6_gnrh$Phase = factor(tog_6_gnrh$Phase, levels = c('M',
                                                       'I',
                                                       'LI',
                                                       'NF',
                                                       'F'))

gnrh_pairs = pairs(emmeans(gnrh_glm, 'Status.x'), adjust = 'none')

tog_6_gnrh$prop = tog_6_gnrh$n_GnRH/tog_6_gnrh$total_cells

gnrh_prop = ggplot(tog_6_gnrh, aes(x = Phase, y = n_GnRH/total_cells))+
    geom_boxplot(data = subset(tog_6_gnrh,Phase != "NF"),
                 aes(x = Phase, y =n_GnRH/total_cells ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
    scale_y_continuous(labels = scales::percent)+
    labs(x = 'Phase', y = '% GnRH of 6')+
    theme(legend.position = 'none')+
      geom_signif(xmin = c(1),
                xmax = c(1.9),
                y_position = c(max(tog_6_gnrh$prop*1.1)),
                annotation =c('p = 0.514'),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
          geom_signif(xmin = c(1),
                xmax = c(5),
                y_position = c(max(tog_6_gnrh$prop*1.3)),
                annotation =c('p = 0.127'),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
    geom_signif(xmin = c(2.1),
                xmax = c(5),
                y_position = c(max(tog_6_gnrh$prop*1.1)),
                annotation =c('p = 0.372'),
                color = "black",
                tip_length = c(0,0), textsize = 3)
  

#ggsave(plot = gnrh_prop,
#       file =  'gnrh_prop.svg',
#       device = "svg",
#       units = "in",
#       width = 1.5,
#       height = 1.5,
#       path = "Manuscript/Plots/Fig.3")


###
sub_6$cckb = ifelse(sub_6@assays$RNA$data['cckb',]>0, T, F)

cckb = sub_6@meta.data%>%
  group_by(individual, Status)%>%
  summarize(n_cckb = sum(cckb == T))

tot= sub_6@meta.data%>%
  group_by(individual, Status)%>%
  summarize(total_cells = n())

tog = cckb%>%
  right_join(tot, by = c('individual'))

tog_6_cckb= subset(tog, Status.x !=c('NRM'))

cckb_glm =glm(cbind(tog_6_cckb$n_cckb, tog_6_cckb$total_cells-tog_6_cckb$n_cckb)~Status.x, 
              data = tog_6_cckb,
              family = binomial('logit'))

anova(cckb_glm, test = 'Chisq')

tog_6_cckb$Phase = as.character(tog_6_cckb$Status.x)
tog_6_cckb$Phase = ifelse(tog_6_cckb$Phase =='D','I', tog_6_cckb$Phase )
tog_6_cckb$Phase = ifelse(tog_6_cckb$Phase =='E','LI', tog_6_cckb$Phase )
tog_6_cckb$Phase = factor(tog_6_cckb$Phase, levels = c('M',
                                                       'I',
                                                       'LI',
                                                       'NF',
                                                       'F'))

cckb_pairs = pairs(emmeans(cckb_glm, 'Status.x'), adjust = 'none')

tog_6_cckb$prop = tog_6_cckb$n_cckb/tog_6_cckb$total_cells

cckb_prop = ggplot(tog_6_cckb, aes(x = Phase, y = n_cckb/total_cells))+
    geom_boxplot(data = subset(tog_6_cckb,Phase != "NF"),
                 aes(x = Phase, y =n_cckb/total_cells ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
    scale_y_continuous(labels = scales::percent)+
    labs(x = 'Phase', y = '% cckb of 6')+
    theme(legend.position = 'none')+
      geom_signif(xmin = c(1),
                xmax = c(1.9),
                y_position = c(max(tog_6_cckb$prop*1.1)),
                annotation =c('p = 0.629'),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
          geom_signif(xmin = c(1),
                xmax = c(5),
                y_position = c(max(tog_6_cckb$prop*1.3)),
                annotation =c('*'),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
    geom_signif(xmin = c(2.1),
                xmax = c(5),
                y_position = c(max(tog_6_cckb$prop*1.1)),
                annotation =c('*'),
                color = "black",
                tip_length = c(0,0), textsize = 3)

ggsave(plot = cckb_prop,
       file =  'cckb_prop.svg',
       device = "svg",
       units = "in",
       width = 1.5,
       height = 1.5,
       path = "Manuscript/Plots/Fig.3")

#### how are the genes grouped ####
library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
  clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)
  library(factoextra)

    mecd = readRDS("Functions/mean_expression_cluster_data.rds")  


  obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_6 = FindSubCluster(obj, 6, graph.name='harmony.wsnn')
Idents(sub_6) <- 'sub.cluster'
sub_6 = subset(sub_6, final_clusters ==6)
sub_6$sub.cluster = factor(sub_6$sub.cluster, levels = c(paste0('6_',0:3)))
sub_6$Status = factor(sub_6$Status, levels = c('NRM','M',"D",'E','NF','F'))


  ### DEG enrichgment
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_6 = degs$gene[degs$cluster==6]
degs_6_df = subset(degs, cluster ==6)

deg_data = data.frame()
for(i in degs_6){
dat = mecd(sub_6, i, '6', 'final_clusters')
dat$gene = i
dat$mean = scale(dat$mean)
deg_data = rbind(dat, deg_data)

}
 labels = data.frame(gene = degs_6_df$gene, short_label = degs_6_df$short_label)
matching_vector= match(deg_data$gene, labels$gene)
 
deg_data$short_label = labels$short_label[matching_vector]

mat_data = deg_data%>%
  dplyr::select(individual, Status, gene, mean)%>%
  pivot_wider(names_from= gene,
              values_from = mean)

piv_mat = mat_data[,-c(1,2)]%>%t()
piv_mat <- as.data.frame(piv_mat)
hab <- labels$short_label[match(rownames(piv_mat), labels$gene)]  
pca = prcomp(piv_mat, scale = T)

fviz_pca_biplot(pca)

fviz_pca_ind(pca)
fviz_pca_ind(pca, habillage = hab)
# would be interesting to find unbiased clusters here









