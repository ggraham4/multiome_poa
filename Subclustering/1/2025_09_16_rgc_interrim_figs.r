### RGC Interrim Figs

{
  library(ggsignif)
  library(patchwork)
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
  library(glmGamPoi)
  library(scran)
  library(parallel)
    library(parallel)

  library(Seurat)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(emmeans)
  library(openxlsx)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
  library(biomaRt)
  library(hdWGCNA)
  
  clown_go = readRDS("Functions/clown_go2")  
  mecd = readRDS("Functions/mean_expression_cluster_data.rds")
    mecp = readRDS("Functions/mean_expression_cluster_plot.rds")

  `%notin%` <- Negate(`%in%`)
  
  #define functions
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

prop_cluster_plot=
function(object, gene, cluster, clustering = 'final_clusters'){
    library(stringr)
    library(forcats)
      options(dplyr.summarise.inform = FALSE)

    counts <- t(as.matrix(object@assays$RNA$data[,object@meta.data[[clustering]] == cluster]))
  Counts_of_interest <- as.data.frame(counts[,gene]>0) #binarize
    Counts_of_interest$expression <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
    Counts_of_interest$Status <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]

  results <-Counts_of_interest%>%
    group_by(individual, Status)%>%
    summarize(mean = mean(expression),
              se = sd(expression)/sqrt(n()))
  results$Sex <- results$Status
  
  
  results$Sex <- str_sub(results$individual, -1)
  results$Sex[results$individual == 'T17D'] = 'NF'
  results$Sex[results$individual == 'A12D'] = 'E'
  results$Sex[results$individual == 'T11D'] = 'E'
  results$Sex[results$individual == 'GH'] = 'NRM'
  
  results$factor <- ifelse(results$Sex == "NRM", 0, NA)
  results$factor <- ifelse(results$Sex == "M", 1, results$factor)
  results$factor <- ifelse(results$Sex == "D", 2, results$factor)
  results$factor <- ifelse(results$Sex == "E", 3, results$factor)
  results$factor <- ifelse(results$Sex == "NF", 4, results$factor)
  results$factor <- ifelse(results$Sex == "F", 5, results$factor)
  results$individual <- fct_reorder(results$individual, results$factor)
  
  results$Sex <- factor(results$Sex, levels = c('NRM','M','D','E','NF','F'))
  plot <- ggplot(results, aes(x = individual, y = mean, color = Sex))+
    geom_boxplot(aes(group = Sex, fill = Sex), alpha = 0.25, outlier.shape = NA)+
    geom_point()+
    geom_pointrange(aes(x = individual, y = mean, ymin = mean - se, ymax = mean+se))+
    theme_classic()+
    labs(x  ='FishID', y = '% Expressing', title = paste0(gene,'_cluster_',cluster))+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
}


}
#read in data
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

####  subcluster and subset ----
Idents(obj) = 'final_clusters'
rgc <- FindSubCluster(obj, 1, 'harmony.wsnn')
rgc = subset(rgc, final_clusters == 1)

Idents(rgc) <- 'sub.cluster'

DimPlot(rgc, label = T, reduction = 'harmony_wnn.umap')

rgc$Status = factor(rgc$Status, levels = c('NRM', 'M','D','E','NF','F'))

### dim ----
dim = DimPlot(rgc, label = T, reduction = 'harmony_wnn.umap')+
  theme_void()+
  theme(legend.position = 'none')
dim


ggsave(plot = dim,
       file = "dimplot.tif",
       device = "tif",
       units = "in",
       width = 4,
       height = 4,
       path = "Subclustering/1/interrim_figs")


#### CytoTRACE ----
cyto = CytoTRACE(rgc@assays$RNA$data%>%as.matrix())
rgc$cyto = cyto$CytoTRACE
FeaturePlot(rgc, 'cyto')
DimPlot(rgc) #woah 1 and 3


cyto = ggplot(rgc@meta.data, aes(x = sub.cluster,y = cyto))+
  geom_boxplot()+
  theme_minimal()+
  labs(x = 'Subcluster', y = 'CytoTRACE')+
  ylim(0,1)

ggsave(plot = cyto,
       file = "cyto.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")

rgc$sub.cluster <- factor(rgc$sub.cluster, levels = c('1_0','1_1','1_2','1_3','1_4','1_5'))

### dotplot ----
dot <- DotPlot(rgc, c(
                'pax6b',
               'dclk1a',
               'dclk2a',
               'fgfr3',
                     'fgfr4',
               'gfap',
                     'nfia',
                     'cd44b',
               'slc1a3a',
               'LOC111577263'
               ),
               group.by = 'sub.cluster',
               dot.min = .15,
                 dot.scale = 3)+
  coord_flip()+
  scale_x_discrete(labels = c(
                     'pax6b',
                     'dclk1a',
                     'dclk2a',
                     'fgfr3',
                     'fgfr4',
                     'gfap',
                     'nfia',
                     'cd44b',
                     'EAA1',
                     'cyp19a1b'))+
  labs(y = 'Subcluster')+
  theme(  axis.title.y = element_blank(), axis.text.x = element_text(angle = -15,
                                                                     hjust =.8,
                                                                     vjust = -0.75))
dot
ggsave(plot = dot,
       file = "dot.tif",
       device = "tif",
       units = "in",
       width = 5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")


### arom mean ----
arom_mean_se = rgc@meta.data%>%
  group_by(sub.cluster)%>%
  summarize(mean = mean(arom),
            sd_arom = sd(arom))

arom = ggplot(arom_mean_se, aes(x = sub.cluster,y = mean))+
  geom_pointrange(aes(ymin = mean - sd_arom, ymax =mean+sd_arom))+
  theme_minimal()+
  labs(x = 'Subcluster', y = 'Normalized cyp19a1b')
arom

ggsave(plot = arom,
       file = "arom_by_cluster.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")

rgc$arom_binom = ifelse(rgc@assays$RNA$data['LOC111577263',]>0, 1, 0)

#proportion cells ----
rgc$Phase  = unfactor(rgc$Status)
rgc$Phase = ifelse(rgc$Phase == 'D', 'I', rgc$Phase)
rgc$Phase = ifelse(rgc$Phase == 'E', 'LI', rgc$Phase )
rgc$Phase <- factor(rgc$Phase, levels = c('NRM',
                                          'M',
                                          'I',
                                          'LI',
                                          'NF',
                                          'F'))
cells_ind = rgc@meta.data%>%
  group_by(individual)%>%
  summarize(n_cells = n())

cells_sub_ind = rgc@meta.data%>%
  group_by(individual, Phase, sub.cluster)%>%
  summarize(n_cells_in = n())

cells_total = cells_ind%>%
  right_join(cells_sub_ind, by = 'individual')
cells_total$prop = cells_total$n_cells_in/cells_total$n_cells

cells_1_3 = subset(cells_total, sub.cluster == '1_3' & Phase %in% c('M','I','F'))
matrix_1_3 = cbind(cells_1_3$n_cells_in,cells_1_3$n_cells-cells_1_3$n_cells_in)

model_1_3 = glmer(matrix_1_3~Phase +(1|individual), family = 'binomial', data = cells_1_3)
car::Anova(model_1_3, 3)
#YES
pairs(emmeans(model_1_3, 'Phase'), adjust ='none') # dominants both


prop_1_3= ggplot(subset(cells_total, sub.cluster =='1_3' & Phase != 'NRM'), aes(x = Phase,
                                                                      y = prop, 
                                                                      fill = Phase))+
  geom_boxplot(alpha = 0.5)+
  geom_point(position = position_jitterdodge(1))+
  labs(y = 'Proportion of Cells')+
  ggtitle('1_3')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_signif(xmin = c(1),
              xmax = c(1.9),
              y_position = c(0.24),
              annotation =c("*"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
    geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(0.24),
              annotation =c("*"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(0.08, 1.15*0.24)+
  theme(legend.position = 'none')

ggsave(plot = prop_1_3,
       file = "prop_1_3.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")

### proportion_ arom+ ------

rgc$arom = rgc@assays$RNA$data['LOC111577263',]

arom_pos_neg = rgc@meta.data%>%
  group_by(Phase, sub.cluster, individual)%>%
  summarize(arom_pos = sum(arom_binom),
            n = n(),
            diff = n - arom_pos)

arom_prop = data.frame()
for(i in 0:5){
  print(i)
  model_dat = subset(arom_pos_neg, Phase %in% c('M','I','F')&
                       sub.cluster == paste0('1_',i))
  
  model_mat = cbind(model_dat$arom_pos, model_dat$diff)
  model= glmer(model_mat~Phase+(1|individual), data = model_dat, family = binomial('probit'))
  av = car::Anova(model, 3)
  av_p =av$`Pr(>Chisq)`[2]
  newd = data.frame(cluster =i,
                    av_p = av_p,
                    singular = isSingular(model),
                    error = length(model@optinfo$warnings))
  arom_prop = rbind(arom_prop, newd)
  pairs(emmeans(model, 'Phase'), adjust ='none')
  
}

prop_arom_1_3 = ggplot(subset(arom_pos_neg, sub.cluster =='1_3' & Phase != 'NRM'), 
                       aes(x = Phase, y = arom_pos/n,fill = Phase)
                       )+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_point(position = position_jitterdodge(1))+
  labs(y = 'Proportion of Cells')+
  ggtitle('cyp19a1b+ 1_3')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.05),
              xmax = c(1.9),
              y_position = c(1.05),
              annotation =c("*"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(1.05),
              annotation =c("**"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(0.45, 1.1)

prop_arom_1_3

ggsave(plot = prop_arom_1_3,
       file = "prop_arom_1_3.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")


prop_arom_1_0 = ggplot(subset(arom_pos_neg, sub.cluster =='1_0' & Phase != 'NRM'), 
                       aes(x = Phase, y = arom_pos/n,fill = Phase)
                       )+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_point(position = position_jitterdodge(1))+
  labs(y = 'Proportion of Cells')+
  ggtitle('cyp19a1b+ 1_0')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.05),
              xmax = c(1.9),
              y_position = c(.8),
              annotation =c("***"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(.8),
              annotation =c("**"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(0.2, 0.85)

prop_arom_1_0

ggsave(plot = prop_arom_1_0,
       file = "prop_arom_1_0.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")


prop_arom_1_1 = ggplot(subset(arom_pos_neg, sub.cluster =='1_1' & Phase != 'NRM'), 
                       aes(x = Phase, y = arom_pos/n,fill = Phase)
                       )+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_point(position = position_jitterdodge(1))+
  labs(y = 'Proportion of Cells')+
  ggtitle('cyp19a1b+ 1_1')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.05),
              xmax = c(1.9),
              y_position = c(.55),
              annotation =c("*"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(.55),
              annotation =c("***"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(0, 0.6)

prop_arom_1_1

ggsave(plot = prop_arom_1_1,
       file = "prop_arom_1_1.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")

prop_arom_1_4 = ggplot(subset(arom_pos_neg, sub.cluster =='1_4' & Phase != 'NRM'), 
                       aes(x = Phase, y = arom_pos/n,fill = Phase)
                       )+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_point(position = position_jitterdodge(1))+
  labs(y = 'Proportion of Cells')+
  ggtitle('cyp19a1b+ 1_4')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.05),
              xmax = c(1.9),
              y_position = c(.55),
              annotation =c("NS"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(.55),
              annotation =c("***"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(0, 0.6)

prop_arom_1_4

ggsave(plot = prop_arom_1_4,
       file = "prop_arom_1_4.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")


#### deg ----
rgc$arom

real_degs = read.csv( 'Subclustering/degs_1_defined_09_12_2025.csv')

rgc_arom_expression = rgc@meta.data%>%
  group_by(sub.cluster, individual, Phase)%>%
  summarize(mean_arom = mean(arom))

arom_1_0 = ggplot(subset(rgc_arom_expression, sub.cluster =='1_0' & Phase != 'NRM'), 
                       aes(x = Phase, y = mean_arom,fill = Phase)
                       )+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_point(position = position_jitterdodge(1))+
  labs(y = 'Normalized Expression')+
  ggtitle('cyp19a1b 1_0')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.05),
              xmax = c(1.9),
              y_position = c(2.55),
              annotation =c("***"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(2.55),
              annotation =c("***"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(0.5, 2.7)
arom_1_0

ggsave(plot = arom_1_0,
       file = "arom_1_0.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")

arom_1_1 = ggplot(subset(rgc_arom_expression, sub.cluster =='1_1' & Phase != 'NRM'), 
                       aes(x = Phase, y = mean_arom,fill = Phase)
                       )+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_point(position = position_jitterdodge(1))+
  labs(y = 'Normalized Expression')+
  ggtitle('cyp19a1b 1_1')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.05),
              xmax = c(1.9),
              y_position = c(.95),
              annotation =c("*"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(.95),
              annotation =c("***"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(0, 1.05)
arom_1_1

ggsave(plot = arom_1_1,
       file = "arom_1_1.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")

arom_1_3 = ggplot(subset(rgc_arom_expression, sub.cluster =='1_3' & Phase != 'NRM'), 
                       aes(x = Phase, y = mean_arom,fill = Phase)
                       )+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_point(position = position_jitterdodge(1))+
  labs(y = 'Normalized Expression')+
  ggtitle('cyp19a1b 1_3')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.05),
              xmax = c(1.9),
              y_position = c(3.85),
              annotation =c("***"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(3.85),
              annotation =c("***"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(0.95, 4.1)
arom_1_3

ggsave(plot = arom_1_3,
       file = "arom_1_3.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")

arom_1_4 = ggplot(subset(rgc_arom_expression, sub.cluster =='1_4' & Phase != 'NRM'), 
                       aes(x = Phase, y = mean_arom,fill = Phase)
                       )+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_point(position = position_jitterdodge(1))+
  labs(y = 'Normalized Expression')+
  ggtitle('cyp19a1b 1_4')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.05),
              xmax = c(5),
              y_position = c(2.85),
              annotation =c("***"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(2.4),
              annotation =c("***"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(0,3.1)
arom_1_4

ggsave(plot = arom_1_4,
       file = "arom_1_4.tif",
       device = "tif",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")






