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

sub_1 = FindSubCluster(obj,1, graph.name='harmony.wsnn')
Idents(sub_1) <- 'sub.cluster'
sub_1 = subset(sub_1, final_clusters ==1)
sub_1$Status = factor(sub_1$Status, levels = c('NRM','M',"D",'E','NF','F'))

dimplot = DimPlot(sub_1, label = F)+
  theme_void()+
  theme(legend.position = "none")
dimplot
DimPlot(sub_1, label = T)

#ggsave(plot = dimplot,
 #      file = "UMAP_sub1.tiff",
  #     device = "tiff",
   #    units = "in",
    #   width = 4,
     #  height = 4,
      # path = "Manuscript/Plots/Fig.4")

### markers would be good
marks_1 = FindAllMarkers(sub_1, only.pos = T)

markers_1_subclusters = DotPlot(sub_1, group.by = 'sub.cluster',
        features = c(
                  'sox2',
          's100b',
                  'gfap',
                  'elavl3',
        #'nkx2.2a',
        'gli1',
        'gli3',
        'foxg1a',
        'fgfr3',
        'LOC111577263', # aromatase
        'LOC111586965', #latent tgfb binding protein4
        'slc6a1b', # gaba transporter
        'slc6a11b', #gaba transporter
        'slc4a4a', # sodium bicarb transporter
        'LOC111571821', # erbb4 like
        'grik2',
        'asic2',
        'shha',
        'six6a', #shared marker for GnRH neurons by the way
        'efna5b',
        'slit3',
        'LOC111567448' ,#glypican 3 like
        'slc1a2b', # glutamate transporter
        'grid1a', #glutamate recepitr
        'LOC111576189', #cub and sushi domain containign 1
        'lmx1bb',
        'zic1',
        'smad9',
        'tgfbi',
        'bmp6',
        'wnt1',
        'LOC111582483', # lmx1.2
        'lmx1al',
        #'lmx1bb',
        'frzb',
        'ptch1',
        'sema3c'
        ),       dot.min = .05,
        scale=T,
        dot.scale = 5)+coord_flip()+
  labs(y = 'Subcluster')+
  theme(axis.title.y = element_blank())
markers_1_subclusters
# tbh this isnt very good so lets just do GO of markers and report that

marks_1_0 = FindMarkers(sub_1, '1_0', only.pos = T)
marks_1_1 = FindMarkers(sub_1, '1_1', only.pos = T)
marks_1_2 = FindMarkers(sub_1, '1_2', only.pos = T)
marks_1_3 = FindMarkers(sub_1, '1_3', only.pos = T)
marks_1_4 = FindMarkers(sub_1, '1_4', only.pos = T)
marks_1_5 = FindMarkers(sub_1, '1_5', only.pos = T)

strong_marks_1_0 = rownames(marks_1_0[marks_1_0$p_val_adj < 0.05, ])
strong_marks_1_1 = rownames(marks_1_1[marks_1_1$p_val_adj < 0.05, ])
strong_marks_1_2 = rownames(marks_1_2[marks_1_2$p_val_adj < 0.05, ])
strong_marks_1_3 = rownames(marks_1_3[marks_1_3$p_val_adj < 0.05, ])
strong_marks_1_4 = rownames(marks_1_4[marks_1_4$p_val_adj < 0.05, ])
strong_marks_1_5 = rownames(marks_1_5[marks_1_5$p_val_adj < 0.05, ])

make_unique = function(target, others) {
  all_others = unique(unlist(others))
  target[!target %in% all_others]
}

all_strong = list(strong_marks_1_0, strong_marks_1_1, strong_marks_1_2,
                  strong_marks_1_3, strong_marks_1_4, strong_marks_1_5)

unique_marks_1_0 = make_unique(strong_marks_1_0, all_strong[-1])
unique_marks_1_1 = make_unique(strong_marks_1_1, all_strong[-2])
unique_marks_1_2 = make_unique(strong_marks_1_2, all_strong[-3])
unique_marks_1_3 = make_unique(strong_marks_1_3, all_strong[-4])
unique_marks_1_4 = make_unique(strong_marks_1_4, all_strong[-5])
unique_marks_1_5 = make_unique(strong_marks_1_5, all_strong[-6])
library(enrichplot)

dotplot(clown_go(unique_marks_1_0))
dotplot(clown_go(unique_marks_1_1))
dotplot(clown_go(unique_marks_1_2))
dotplot(clown_go(unique_marks_1_3))
dotplot(clown_go(unique_marks_1_4))
dotplot(clown_go(unique_marks_1_5))

for(i in 0:5){
  dot = dotplot(clown_go(get(paste0('unique_marks_1_', i)))) +
    labs(title = paste0('Cluster 1_', i, ' Markers')) +
    theme(axis.text.y = element_text(size = 4),
          text = element_text(size = 4),
          axis.text.x = element_text(size = 4),
          axis.title.x = element_text(size = 4))
  
  #ggsave(plot = dot,
   #      file = paste0("go_1_", i, ".svg"),
    #     device = 'svg',
     #    system_fonts = list(sans = "Helvetica"),  # force font
      #   units = "in",
       #  width = 4,
        # height = 4,
         #path = "Manuscript/Plots/Fig.4/GO")
} 


markers_1_subclusters_better = DotPlot(sub_1, group.by = 'sub.cluster',
        features = c(unique_marks_1_0[1:5],
                     unique_marks_1_1[1:5],
                     unique_marks_1_2[1:5],
                     unique_marks_1_3[1:5],
                     unique_marks_1_4[1:5],
                     unique_marks_1_5[1:5]
        ),       dot.min = .05,
        scale=T,
        dot.scale = 5)+coord_flip()+
  labs(y = 'Subcluster')+
  theme(axis.title.y = element_blank())
markers_1_subclusters_better

# idk I am not convinced this is good either, lets just do the bread and butter proportions and cytotrace#

### CytoTRACE ####
library(CytoTRACE)
cyto = CytoTRACE(sub_1@assays$RNA$data%>%as.matrix())
sub_1$cyto = cyto$CytoTRACE
FeaturePlot(sub_1, 'cyto')
DimPlot(sub_1) #woah 1 and 3

cyto_plot = sub_1@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(mean_cyto = mean(cyto),
            se_cyto = sd(cyto)/sqrt(n()))
library(lme4)
library(multcomp)
av = lmer(mean_cyto~sub.cluster+(1|individual), data = subset(cyto_plot, Status != 'NRM'))
car::Anova(av,3)

cld_results = cld(emmeans(av, 'sub.cluster'), Letters = letters, adjust = "none", alpha = 0.05)
cld_df = as.data.frame(cld_results)

cytoplot = ggplot(subset(cyto_plot, Status != 'NRM'), aes(x = sub.cluster,y = mean_cyto))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(shape = 1, size =2)+
  theme_minimal()+
  labs(x = 'Subcluster', y = 'Mean CytoTRACE')+
    geom_text(data = cld_df, aes(x = sub.cluster, y = 1, label = .group), 
          size = 3, inherit.aes = FALSE)+
  ylim(0, 1)
  

#ggsave(plot = cytoplot,
##       file = "sub_1_cyto.svg",
#       device = 'svg',
#       units = "in",
#       width = 2,
#       height = 2,
#       path = "Manuscript/Plots/Fig.4")


### proportions ####
cells_ind = sub_1@meta.data%>%
  group_by(individual)%>%
  summarize(n_cells = n())

cells_sub_ind = sub_1@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(n_cells_in = n())

cells_total = cells_ind%>%
  right_join(cells_sub_ind, by = 'individual')
cells_total$prop = cells_total$n_cells_in/cells_total$n_cells

ggplot(subset(cells_total, Status %in% c('M','D',"F")), aes(x = sub.cluster, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(1))

ggplot(subset(cells_total, sub.cluster =='1_0'), aes(x = sub.cluster, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(1))
## i would bet there is a difference here between m and f

ggplot(subset(cells_total, sub.cluster =='1_1'), aes(x = sub.cluster, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(1))
#def something with dominants here

out_df = data.frame()
for(i in 0:5){
cells = subset(cells_total, sub.cluster == paste0('1_',i) & Status !='NRM')
matrix_1 = cbind(cells$n_cells_in,cells$n_cells-cells$n_cells_in)

matrix_1 = glm(matrix_1~Status, family = binomial(), data = cells)
av_ = car::Anova(matrix_1, 3)
pairs = pairs(emmeans(matrix_1, 'Status'), adjust ='none')%>%as.data.frame()

newd = data.frame(cluster = i,
                  av_p = av_$`Pr(>Chisq)`[1],
                 d_m_p =  pairs$p.value[pairs$contrast == 'M - D'],
                 f_m_p =  pairs$p.value[pairs$contrast == 'M - F'],
                 d_f = pairs$p.value[pairs$contrast == 'D - F'])
out_df = rbind(out_df, newd)

}
out_df$signif = ifelse(out_df$av_p<0.05, '*', NA) # not comparing for multiple comparisons because its only 6

# differences in 0:3
cells_total$Phase = as.character(cells_total$Status)
cells_total$Phase = ifelse(cells_total$Status=='D', 'I', cells_total$Phase)
cells_total$Phase = ifelse(cells_total$Phase=='E', 'LI', cells_total$Phase)
cells_total$Phase = factor(cells_total$Phase, levels = c('M', 'I', "LI", 'NF','F'))

# 1_0 
cells_total_1_0 = subset(cells_total, sub.cluster == paste0('1_',0) & Status !='NRM')
textsize=3

prop_1_0 = ggplot(cells_total_1_0,
       aes(x =Phase, y = n_cells_in/n_cells ))+
    geom_boxplot(data = subset(cells_total_1_0, Phase != "NF"),
                 aes(x = Phase, y =n_cells_in/n_cells ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
    scale_y_continuous(labels = scales::percent)+
    labs(x = 'Phase', y =  '% of 1_RGC')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.55), annotation =c('p=0.110'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.6), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.55), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.32)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')

#ggsave(plot = prop_1_0,
 #      file = "prop_1_0.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 1.5,
     #  height = 1.5,
      # path = "Manuscript/Plots/Fig.4/prop")

#1_1

cells_total_1_1 = subset(cells_total, sub.cluster == paste0('1_',1) & Status !='NRM')
prop_1_1 = ggplot(cells_total_1_1,
       aes(x =Phase, y = n_cells_in/n_cells ))+
    geom_boxplot(data = subset(cells_total_1_1, Phase != "NF"),
                 aes(x = Phase, y =n_cells_in/n_cells ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
    scale_y_continuous(labels = scales::percent)+
    labs(x = 'Phase', y =  '% of 1_RGC')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.4), annotation =c('p=0.110'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.45), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.4), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.32)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')
prop_1_1
# no pairwise differences probably LIs driving the anova

## 1_2 ###
cells_total_1_2 = subset(cells_total, sub.cluster == paste0('1_',2) & Status !='NRM')

prop_1_2 = ggplot(cells_total_1_2,
       aes(x =Phase, y = n_cells_in/n_cells ))+
    geom_boxplot(data = subset(cells_total_1_2, Phase != "NF"),
                 aes(x = Phase, y =n_cells_in/n_cells ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
    scale_y_continuous(labels = scales::percent)+
    labs(x = 'Phase', y =  '% of 1_RGC')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.32), annotation =c('p=0.662'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.38), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.32), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')
prop_1_2 # how the fuck is this true, I really dont want to report this...

#ggsave(plot = prop_1_2,
 #      file = "prop_1_2.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 1.5,
     #  height = 1.5,
      # path = "Manuscript/Plots/Fig.4/prop")

### 1_1 ### 
cells_total_1_1 = subset(cells_total, sub.cluster == paste0('1_',3) & Status !='NRM')

prop_1_1 = ggplot(cells_total_1_1,
       aes(x =Phase, y = n_cells_in/n_cells ))+
    geom_boxplot(data = subset(cells_total_1_1, Phase != "NF"),
                 aes(x = Phase, y =n_cells_in/n_cells ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
    scale_y_continuous(labels = scales::percent)+
    labs(x = 'Phase', y =  '% of 1_RGC')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.25), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.27), annotation =c('p=0.524'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.25), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')
prop_1_1 

ggsave(plot = prop_1_1,
       file = "prop_1_1.svg",
       device = "svg",
       units = "in",
       width = 1.5,
       height = 1.5,
       path = "Manuscript/Plots/Fig.4/prop")

#### Aromatase ####
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_arom = subset(degs, gene =='LOC111577263')

sub_1@meta.data$Phase = as.character(sub_1@meta.data$Status)
sub_1@meta.data$Phase= ifelse(sub_1@meta.data$Phase == 'D', 'I', sub_1@meta.data$Phase)
sub_1@meta.data$Phase= ifelse(sub_1@meta.data$Phase == 'E', 'LI', sub_1@meta.data$Phase)

sub_1@meta.data$Phase = factor(sub_1@meta.data$Phase, levels = c('NRM',
                                                             'M',
                                                             'I',
                                                             'LI',
                                                             'NF',
                                                             'F'))


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

plotter_function_final = function(obj,data, deg){
  
  deg_data = subset(data, gene == deg & cluster ==1)
  d_f_p.value = deg_data$d_f_p.value
  d_m_p.value = deg_data$d_m_p.value
  f_m_p.value = deg_data$f_m_p.value
  
  d_f_annotation = p_value_annotation(d_f_p.value)
  d_m_annotation = p_value_annotation(d_m_p.value)
  f_m_annotation = p_value_annotation(f_m_p.value)
  
  meta.data = subset(obj@meta.data, final_clusters ==1)
  meta.data$expression = obj@assays$RNA$data[deg,obj@meta.data$final_clusters == 1]
  
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
      theme(legend.position = 'none',)
    r
  return(r)
}

plotter_function_sub = function(obj,data, deg, cluster){
  
  deg_data = subset(data, gene == deg)
  d_f_p.value = deg_data$d_f_p.value
  d_m_p.value = deg_data$d_m_p.value
  f_m_p.value = deg_data$f_m_p.value
  
  d_f_annotation = p_value_annotation(d_f_p.value)
  d_m_annotation = p_value_annotation(d_m_p.value)
  f_m_annotation = p_value_annotation(f_m_p.value)
  
  meta.data = subset(obj@meta.data, sub.cluster == cluster)
  meta.data$expression = obj@assays$RNA$data[deg,obj@meta.data$sub.cluster == cluster]
  
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
      theme(legend.position = 'none',axis.title.y = element_blank())
    r
  return(r)
}


arom_whole = plotter_function_final(sub_1, degs_arom, 'LOC111577263')+
  labs( y = 'mean expression')

 # ggsave(plot = arom_whole,
  #     file = 'arom_all_1.svg',
   #    device = "svg",
    #   units = "in",
     #  width = 1.5,
      # height = 1.5,
       #path = "Manuscript/Plots/Fig.4/arom")


  ### ok now subclusters ###
  rgc_sub_degs =data.frame()

  for(i in c(paste0('1_', 0:5))){
    dat = read.csv(paste0('~/Desktop/multiome_poa/DEG Outputs/09_16_2025 1 Subclusters Neg Bin Anova First/cluster_', i,'.csv'))
    dat_sub = subset(dat, av_q.value < 0.05)
    dat_sub$cluster == i
    rgc_sub_degs= rbind(rgc_sub_degs, dat_sub)
    
  }
  
  arom_df = subset(rgc_sub_degs, gene =='LOC111577263')
  
  for(i in arom_df$cluster){
    
    arom_sub = plotter_function_sub(sub_1, arom_df[arom_df$cluster ==i,], 'LOC111577263', i)
     
   #  ggsave(plot = arom_sub,
  #     file = paste0('arom_',i,'.svg'),
  #     device = "svg",
  #     units = "in",
  #     width = 1.35,
  #     height = 1.5,
  #     path = "Manuscript/Plots/Fig.4/arom")

    
  }
  
  plotter_function_sub(sub_1, rgc_sub_degs[rgc_sub_degs$cluster =='1_5',], 'avp', '1_5')
#what  
  
  ### Transcription factors ####
  degs = read.csv('DEG Outputs/FINAL degs classified w singular.csv')

degs_1_concise = degs[degs$cluster==1, c('gene', 'gene_name','type', 'short_label', 'first_word', 'second_word')]
degs_1= subset(degs, cluster == 1)

eyes = table(degs_1_concise$type, degs_1_concise$first_word)%>%as.data.frame.matrix()

# transcription factors are predominantly late
#other gene types are less structured

# maybe a model where ERE signaling -> these txnal changes

# I think it would be good to plot these TFs too

tf_degs =degs_1_concise$gene[degs_1_concise$type == 'transcription factor']

for(i in c(tf_degs)){
  
  plok = plotter_function_final(sub_1, degs_1, i)+
    labs( y = i)
  
       ggsave(plot = plok,
       file = paste0('tf_',i,'.svg'),
       device = "svg",
       units = "in",
       width = 1.35,
       height = 1.5,
       path = "Manuscript/Plots/Fig.4/tfs")

  
}

#glceb another ECM gene, and we could discuss junctional adhesion molecule 3B-like
# and gap junction gamma-1 protein-like as well in that same sentence maybe 





  