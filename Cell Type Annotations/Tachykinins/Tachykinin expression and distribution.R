### looking at tachykinin and tachykinin receptor expression in the obj

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
  library(glmGamPoi)
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)
  library(emmeans)
  library(CytoTRACE)
  library(ggrepel)
  
  library(Polychrome)
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
names(P40) <- NULL

mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
define_degs_prop<- readRDS('Functions/define_degs_prop.rds')
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
prop_deg_function.rds<- readRDS('Functions/DEG_functions/prop_deg_function.rds')
define_behavior_degs<- readRDS('Functions/define_behavior_degs')
clown_go<- readRDS('Functions/clown_go')

}

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

genes <- c('elavl3',
           'gad2', #GABA
           'LOC111588076', #GABA
           'LOC111584103', #GLU
           'slc17a6b', #glu
           'tac1',
           'tacr1a',
          'tacr1b',
          'tacr2',
           'tac3a',
          'tacr3a',
          'tacr3l',
          'npy',
          'kiss1',
          'kiss2',
          'kiss1ra',
          'kiss1rb',
          'esr1',
          'esr2a',
          'esr2b',
          'ar',
          'LOC111571064',
          'gnrh2',
          'gnrh3',
           'gnrhr1',
          'LOC111568468', #gnrhr2
          'gnrhr4'
          )

DotPlot(object = obj, 
                 group.by = "harmony.wnn_res0.4_clusters", 
                 features = genes
        ) + 
  coord_flip()+
  geom_hline(yintercept = 20, linetype = 2, alpha = 0.5)
          


mean_tacr3a <- mean_expression_cluster_data(obj,
                                            'tacr3a',
                                            19)
mean_kt <- data.frame(kt = obj@meta.data$Log_11KT,
                      individual = obj@meta.data$individual)%>%
  group_by(individual)%>%
  summarize(kt= mean(kt))

mean_tacr3a%>%
  right_join(mean_kt, by='individual')%>%
  ggplot(aes(x = kt, y = mean, color = Sex))+
  geom_smooth(method = 'lm')+
  geom_pointrange(aes(x = kt, y = mean, ymin =mean-se, ymax = mean+se))


mean_tac3a <- mean_expression_cluster_data(obj,
                                            'tac3a',
                                            19)
mean_tac3a%>%
  right_join(mean_kt, by='individual')%>%
  ggplot(aes(x = kt, y = mean, color = Sex))+
  geom_smooth(method = 'lm')+
  geom_pointrange(aes(x = kt, y = mean, ymin =mean-se, ymax = mean+se))


pgr_tacr3a <- mean_expression_cluster_data(obj,
                                            'pgr',
                                            19)%>%
  right_join(mean_tacr3a, by= 'individual')

    ggplot(pgr_tacr3a,aes(x = mean.y, y = mean.x, color = Sex.y))+
  geom_smooth(method = 'lm')+
      geom_point()+
      labs(x = 'tacr3a',y = 'pgr')
  





