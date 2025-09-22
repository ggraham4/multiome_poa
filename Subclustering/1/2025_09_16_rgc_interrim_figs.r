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
  
  geneNamer = function(gene){
  names = read.csv('Reference/genes updated.csv')
  
  name = names$NIH_description[names$NIH_accession==gene][1]
  
  if(is.na(name)){name = gene}
  return(name)
}

  
  #define functions
p_annotate <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else if (p_value < 0.1) {
    return(paste0("p = ", round(p_value, 3)))
  } else {
    return("NS")
  }
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
rgc$sub.cluster= factor(rgc$sub.cluster, levels = c(paste0('1_', 0:5)))
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
rgc$Phase  = (rgc$Status)
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

### wgcna
### WGCNA ----
clusters_list <- c('1_0',
                   '1_1',
                   '1_2',
                   '1_3',
                   '1_4',
                   '1_5')

rgc <- SetupForWGCNA(
  rgc,
  gene_select = "fraction",
  fraction = 0.05, 
  wgcna_name = "rgc" 
)

### Construct metacells
rgc <- MetacellsByGroups(
  seurat_obj = rgc,
  group.by = c("sub.cluster"),
  reduction = 'harmony_wnn.umap', 
  k = 25, 
  max_shared = 10,
  ident.group = 'sub.cluster',
  min_cells = 90 
)

# normalize metacell expression matrix:
rgc <- NormalizeMetacells(rgc)

### Coexpression Network analysis
rgc <- SetDatExpr(
  rgc,
  group_name = clusters_list,
  group.by='sub.cluster', 
  assay = 'RNA', 
  layer = 'data' 
)

### Select soft power threshold
rgc <- TestSoftPowers(
  rgc,
  networkType = 'signed' #not sure if this should be altered either
)

# plot the results:
plot_list <- PlotSoftPowers(rgc)
wrap_plots(plot_list, ncol=2)

rgc <- ConstructNetwork(
  rgc,
  tom_name = 'rgc', 
  tom_outdir = '/Users/ggraham/WGCNA_rgc/',
  overwrite_tom = T
)

PlotDendrogram(rgc, main='clust_1 hdWGCNA Dendrogram')

#model eigengenes
rgc <- ModuleEigengenes(
  rgc,
  group.by.vars="sub.cluster"
)

# harmonized module eigengenes:
hMEs <- GetMEs(rgc)
# module eigengenes:
MEs <- GetMEs(rgc, harmonized=FALSE)

rgc@misc[["rgc"]][["wgcna_net"]][["TOMFiles"]] <-'/Users/ggraham/WGCNA_rgc/rgc_TOM.rda'

# compute eigengene-based connectivity (kME):
rgc <- ModuleConnectivity(
  rgc,
  group.by = 'sub.cluster', group_name =  clusters_list
  
)

p <- PlotKMEs(rgc, ncol=5)
p

# get the module assignment table:
modules <- GetModules(rgc) %>% subset(module != 'grey')

### Get hub genes
hub_df <- GetHubGenes(rgc, n_hubs = 10)

#enrich hub genes
enrich_df = data.frame()
for(i in unique(hub_df$module)){
  #print(i)
  genes = hub_df$gene_name[hub_df$module==i]
  
  go =clown_go(genes)
  if(length(go$Description)>1){
  message(paste0(i),paste0('_', go$Description)) 
  newd = data.frame(module = i,
                    description = go$Description)
  enrich_df = rbind(newd, enrich_df)
  }
  else{
    message(paste0(i,' no enrichment'))
  }
}


MEs$Status = rgc$Status
MEs$individual = rgc$individual
MEs$sub.cluster = rgc$sub.cluster

plot_module('turquoise')
plot_module('yellow')

dmes = data.frame()
for(subcluster in unique(rgc$sub.cluster)){
for(i in unique(modules$color)){
  formula_string <- paste0(i, " ~ Status +(1|individual)")%>%as.formula()
  model = lmer(formula_string, data = subset(MEs, Status %in% c('M','D','F') & 
                                               sub.cluster == subcluster))
  p =car::Anova(model, 3)
  newd = data.frame(module = i,
                    subcluster = subcluster,
                    pval = p$`Pr(>Chisq)`[2])
  dmes = rbind(dmes, newd)
  
}
}

radar = ModuleRadarPlot(
  rgc,
  group.by = 'sub.cluster',
  axis.label.size=4,
  grid.label.size=4,
  combine = F
)
radar

clown_go(modules$gene_name[modules$color=='red'])%>%
  dotplot() # neurogenesis
clown_go(modules$gene_name[modules$color=='yellow'])%>%
  dotplot() # cell diff
clown_go(modules$gene_name[modules$color=='blue'])%>%
  dotplot() # translation


for(i in names(radar)){
  
  plot = radar[[i]]
  ggsave(plot = plot,
       file = paste0("radar_",i,'.svg'),
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")

}
rgc$blue = MEs$blue
grouped_and_sded_blue = subset(MEs, sub.cluster == '1_3'& Phase != 'NRM')%>%
  group_by(individual,Phase )%>%
  summarize(mean_blue = mean(blue),
            se_blue = sd(blue)/sqrt(n()))

model = lmer(blue~Status+(1|individual), data = subset(MEs, Status %in% c('M','D','F') & 
                                               sub.cluster == '1_3'))
car::Anova(model, 3)
pairs(emmeans(model, 'Status'), adjust = 'none')

set.seed(10)
blue_1_3 = ggplot(grouped_and_sded_blue, aes(x = Phase, y = mean_blue,fill = Phase))+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_pointrange(aes(x = Phase, y = mean_blue, ymin = mean_blue-se_blue,
                      ymax= mean_blue+se_blue),position = position_jitterdodge(1))+
  labs(y = 'Blue Mean +/- SE')+
  ggtitle('1_3')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.0),
              xmax = c(1.9),
              y_position = c(4.2),
              annotation =c("p=0.07"), 
              color = "black",
              tip_length = c(0,0),
              textsize=4)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(4.2),
              annotation =c("*"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(-2, 4.7)

  
ggsave(plot = blue_1_3,
       file = 'blue_1_3.tiff',
       device = "tiff",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")



### how are the red hub dfs correlated ------

ggplot(data = data.frame(), aes(x = rgc@assays$RNA$data['LOC111577263',],
                                y = rgc@assays$RNA$data['foxg1a',]))+
  geom_point()

ggplot(data = data.frame(), aes(x = rgc@assays$RNA$data['LOC111577263',],
                                y = rgc@assays$RNA$data['sox1a',]))+
  geom_point()

ggplot(data = data.frame(), aes(x = rgc@assays$RNA$data['LOC111577263',],
                                y = rgc@assays$RNA$data['lhx2b',]))+
  geom_point()

mecp(rgc, 'lhx2b', '1_3', 'sub.cluster')
mecp(rgc, 'sox1a', '1_3', 'sub.cluster')
mecp(rgc, 'foxg1a', '1_3', 'sub.cluster')
mecp(rgc, 'LOC111577260', '1_3', 'sub.cluster')

#### DNA damage module ====
#rip the GO term
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")


biomart_basic <-getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'entrezgene_description',
                   'go_id',
                   'name_1006',
                   'namespace_1003'))

dna_damage_genes = unique(
  biomart_basic$entrezgene_accession[
                                       grepl('break repair', biomart_basic$name_1006) |
                                       grepl('DNA repair', biomart_basic$name_1006) |
                                       grepl('DNA damage', biomart_basic$name_1006) |
                                       grepl('DNA damage', biomart_basic$entrezgene_description)
                                     ]
  )
dna_damage_genes

rgc = AddModuleScore(rgc, list(dna_damage_genes),name = 'Damage')

DotPlot(rgc, 'Damage1')

plotSeuratModule('Damage1', '1_0')
plotSeuratModule('Damage1', '1_1')
plotSeuratModule('Damage1', '1_2') #woooah
plotSeuratModule('Damage1', '1_3')
plotSeuratModule('Damage1', '1_4')
plotSeuratModule('Damage1', '1_5')

data_1_2 = subset(rgc@meta.data, sub.cluster == '1_2')
damage_model2 = lmer(Damage1 ~ Status+(1|individual), data =data_1_2 )
car::Anova(damage_model2, 3) # big money but what does it mean

pairs(emmeans(damage_model2,'Status'), adjust ='none')

grouped_and_sded_dmg = subset(rgc@meta.data, sub.cluster == '1_2'& Phase != 'NRM')%>%
  group_by(individual,Phase )%>%
  summarize(mean_dmg = mean(Damage1),
            se_dmg= sd(Damage1)/sqrt(n()))

set.seed(10)
dmg_1_2 = ggplot(grouped_and_sded_dmg, aes(x = Phase, y = mean_dmg,fill = Phase))+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_pointrange(aes(x = Phase, y = mean_dmg, ymin = mean_dmg-se_dmg,
                      ymax= mean_dmg+se_dmg),position = position_jitterdodge(1))+
  labs(y = 'DNA Damage Mean +/- SE')+
  ggtitle('1_2')+
    theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.0),
              xmax = c(1.9),
              y_position = c(.01),
              annotation =c("*"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(.01),
              annotation =c("**"), 
              color = "black",
              tip_length = c(0,0),
              textsize=6)+
  ylim(-0.02, 0.015)

ggsave(plot = dmg_1_2,
       file = 'dmg_1_2.tiff',
       device = "tiff",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")


###### growth factors and whatnot ----
lower_string_degs = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/cluster_1.csv')
degs = subset(lower_string_degs, av_q.value < 0.1)
degs$name = sapply(degs$gene, geneNamer)
sparse_degs = degs[,c('gene', 'name', 'singular')]   

all_degs_1_1 =  read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/09_16_2025 1 Subclusters Neg Bin Anova First/cluster_1_1.csv')

#define a function to plot DEGs the way I like
deg_plotter = function(object = rgc, 
                       gene, 
                       cluster , 
                       clustering='sub.cluster',
                       signif_dm  ,
                       signif_df ,
                       signif_mf ,
                       singular=F,
                       common_name){
  set.seed(10)
  singular = ifelse(singular == T, 'Singular', '')
  meta= object@meta.data
  meta$gene = object@assays$RNA$data[gene,]
  
  meta_grouped_and_sded = meta%>%
    filter(Phase != 'NRM' & !!sym(clustering) == cluster) %>%
    group_by(individual, Phase)%>%
    summarize(mean_gene = mean(gene),
              se_gene = sd(gene)/sqrt(n()))
  
  plot_lower_lim = min(meta_grouped_and_sded$mean_gene -meta_grouped_and_sded$se_gene )
  plot_upper_lim= max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.2
  plot_signif_lower = max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.05
  plot_signif_upper = max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.25
  

  
    signif_dm = p_annotate(signif_dm)
    signif_df = p_annotate(signif_df)
    signif_mf = p_annotate(signif_mf)
  
  textsize_dm = ifelse(grepl("\\*", signif_dm), 6, 3)  
  textsize_df = ifelse(grepl("\\*", signif_df), 6, 3)  
  textsize_mf = ifelse(grepl("\\*", signif_mf), 6, 3)    

  
      plot_upper_lim= ifelse(signif_mf!= 'NS', max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.4,plot_upper_lim )

  
plot = ggplot(meta_grouped_and_sded, aes(x = Phase, y = mean_gene,fill = Phase))+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_pointrange(aes(x = Phase,
                      y = mean_gene,
                      ymin = mean_gene-se_gene,
                      ymax= mean_gene+se_gene),
                  position = position_jitterdodge(1), 
                  size = 0.2
                  )+
  labs(y = 'Mean +/- SE Expression', subtitle = singular)+
  ggtitle(paste0(common_name, ': ', cluster))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size =12),
        plot.subtitle = element_text(hjust = 0.5, size =8))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.0),
              xmax = c(1.9),
              y_position = c(plot_signif_lower),
              annotation =c(signif_dm), 
              color = "black",
              tip_length = c(0,0),
              textsize=textsize_dm)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(plot_signif_lower),
              annotation =c(signif_df), 
              color = "black",
              tip_length = c(0,0),
              textsize=textsize_df)+
  ylim(plot_lower_lim, plot_upper_lim)
   
 if(signif_mf != 'NS'){
   plot  <- plot+
      geom_signif(xmin = c(1),
              xmax = c(5),
              y_position = c(plot_signif_upper),
              annotation =c(signif_mf), 
              color = "black",
              tip_length = c(0,0),
              textsize=textsize_mf)
    
  }

return(plot)
  
}


#rarb
rarb = deg_plotter(rgc, 
            'LOC111573403',
            '1_1',
            'sub.cluster',
            real_degs$d_m_p.value[real_degs$gene=='LOC111573403'],
            real_degs$d_f_p.value[real_degs$gene=='LOC111573403'],
            real_degs$f_m_p.value[real_degs$gene=='LOC111573403'],
            real_degs$singular[real_degs$gene=='LOC111573403'],
            'rarb'
            )

ggsave(plot = rarb,
       file = 'rarb_1_1.tiff',
       device = "tiff",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")



grb10 = deg_plotter(rgc, 
            'grb10b',
            '1_1',
            'sub.cluster',
            all_degs_1_1$d_m_p.value[all_degs_1_1$gene=='grb10b'],
            all_degs_1_1$d_f_p.value[all_degs_1_1$gene=='grb10b'],
            all_degs_1_1$f_m_p.value[all_degs_1_1$gene=='grb10b'],
            singular = all_degs_1_1$singular[all_degs_1_1$gene=='grb10b'],
            'grb10b'
            )

grb10

ggsave(plot = grb10,
       file = 'grb10_1_1.tiff',
       device = "tiff",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")


ccdc85cb = deg_plotter(rgc, 
            'ccdc85cb',
            '1_1',
            'sub.cluster',
            all_degs_1_1$d_m_p.value[all_degs_1_1$gene=='ccdc85cb'],
            all_degs_1_1$d_f_p.value[all_degs_1_1$gene=='ccdc85cb'],
            all_degs_1_1$f_m_p.value[all_degs_1_1$gene=='ccdc85cb'],
            singular = all_degs_1_1$singular[all_degs_1_1$gene=='ccdc85cb'],
            'ccdc85cb'
            )

ccdc85cb

ggsave(plot = ccdc85cb,
       file = 'ccdc85cb_1_1.tiff',
       device = "tiff",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")

 

fgfbp3 = deg_plotter(rgc, 
            'fgfbp3',
            '1',
            'final_clusters',
            lower_string_degs$d_m_p.value[lower_string_degs$gene=='fgfbp3'],
            lower_string_degs$d_f_p.value[lower_string_degs$gene=='fgfbp3'],
            lower_string_degs$f_m_p.value[lower_string_degs$gene=='fgfbp3'],
            singular = lower_string_degs$singular[lower_string_degs$gene=='fgfbp3'],
            'fgfbp3'
            )

fgfbp3

ggsave(plot = fgfbp3,
       file = 'fgfbp3_1.tiff',
       device = "tiff",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")

fgf12a = deg_plotter(rgc, 
            'fgf12a',
            '1_1',
            'sub.cluster',
            lower_string_degs$d_m_p.value[lower_string_degs$gene=='fgfbp3'],
            lower_string_degs$d_f_p.value[lower_string_degs$gene=='fgfbp3'],
            lower_string_degs$f_m_p.value[lower_string_degs$gene=='fgfbp3'],
            singular = lower_string_degs$singular[lower_string_degs$gene=='fgfbp3'],
            'fgf12a'
            )

fgf12a

ggsave(plot = fgf12a,
       file = 'fgf12a_1_1.tiff',
       device = "tiff",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")



dot2<- DotPlot(rgc, c('fgf2',
               'fgf3',
               'fgf7',
               'fgf10a',
                     'fgf10b',
               'fgf22'
               ),
               group.by = 'sub.cluster',
               dot.min = .01,
                 dot.scale = 3)+
  labs(y = 'Subcluster')+
  theme(  axis.title.y = element_blank(), axis.text.x = element_text(angle = -90,
                                                                     hjust =.8,
                                                                     vjust = 0.75))
dot2

ggsave(plot = dot2,
       file = 'dot2.tiff',
       device = "tif",
       units = "in",
       width = 4,
       height = 2.5,
       path = "Subclustering/1/interrim_figs")







