library(ggplot2)
library(ggsignif)
library(tidyverse)
library(emmeans)
library(Seurat)

obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

obj_subset = obj[,obj@meta.data$harmony.wnn_res0.4_clusters!=15 & 
                   obj@meta.data$harmony.wnn_res0.4_clusters!=30&
                   obj@meta.data$Status != 'NRM']


gabaergic <- c('1',
               '5',
               '20',
               '25',
               '28',
               '31')

glutamatergic <- c('3',
                   '6',
                   '7',
                   '8',
                   '9',
                   '10',
                   '11',
                   '12',
                   '13',
                   '16',
                   '17',
                   '21',
                   '23',
                   '24',
                   '27')
mixed <- c('0','19')

obj_subset$type <- ifelse(obj_subset$harmony.wnn_res0.4_clusters %in% gabaergic, 'GABA',
                          ifelse(obj_subset$harmony.wnn_res0.4_clusters %in% glutamatergic,
                                 'GLUT',
                                 ifelse(obj_subset$harmony.wnn_res0.4_clusters %in% mixed, 'Mixed',
                                 'Non-Neuronal'))
)

type_plott <- DimPlot(obj_subset, group.by = 'type', label = F, reduction = 'harmony_wnn.umap')

arr <- list(x = -15, y = -15, x_len = 4, y_len = 4)

type_plott2 <- type_plott+
  theme_void()+
  theme(plot.title = element_blank())+
  annotate('text',
           x = -17, y = -12.5, label = 'UMAP_2', angle = 90, size =3)+
  annotate('text',
           x = -12.5, y = -17, label = 'UMAP_1', size =3)+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), 
           arrow = arrow(type = "closed", length = unit(10, 'pt'))) 
type_plott2

ggsave(plot = type_plott2,
       file = "type_plott2.tiff",
       device = "tiff",
       units = "in",
       width = 3.6,
       height = 2.5,
       path = "Bachelors Thesis/Plots/Figure 2",
       dpi = 300)

  
