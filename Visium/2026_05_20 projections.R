library(Seurat)
library(ggplot2)
library(tidyverse)

vis_neuron=  readRDS('/Users/ggraham/Desktop/Visium/vis_neuron.rds' )

vis_neuron_subset= subset(vis_neuron, Slice=='6P17')
new_names <- gsub("_\\d+$", "", colnames(vis_neuron_subset))
vis_neuron_subset <- RenameCells(vis_neuron_subset, new.names = new_names)

dissection = read.csv("~/Desktop/multiome_poa/Visium/Groupings/dissectioncsv.csv")

vis_neuron_subset$dissection = ifelse(colnames(vis_neuron_subset)%in%c(dissection$Barcode),
                                      T,
                                      F)

vis_neuron_subset= subset(vis_neuron_subset, dissection == T)

clusters = c(6,3,14,0
             ,7,24
             )

#9 10 14 15 17 18 26
colors <- c('#DC050C','#4EB265','#EE8026', '#1965B0', '#7BAFDE','#CAE0AB' ,'#882E72' ,
            'black')
named_colors <- setNames(colors[1:length(clusters)], clusters)



p =SpatialDimPlot(vis_neuron_subset%>%
                 subset(predicted.multiome_confidence_0.5%in%clusters), group.by = 'predicted.multiome_confidence_0.5',
               pt.size.factor = 0.5,
               images = 's_6P17.polygons',
               cols = named_colors)
p

ggsave(plot = p,
       file = "6p17_multione_0_3_14_7_24_dissection_aware.tiff",
       device = "tiff",
       units = "in",
       width = 5,
       height = 5,
       path = "Manuscript/Plots/Manuscript v1.2.1/visium/")


clusters = c(6,3,14,0
             ,7,24
             )

#9 10 14 15 17 18 26
colors <- c('#DC050C','#4EB265','#EE8026', '#1965B0', '#7BAFDE','#CAE0AB' ,'#882E72' ,
            'black', 'brown')
named_colors <- setNames(colors[1:length(clusters)], clusters)


clust2 = c(
  4,#3
  5,#6
  #8,#9
  17, #10
  18, #12
  22, #14
  10,#15
  16,#16
  23,#24
  0,#21
  24#26
)
colors2 <- c(
            '#F4A736',
            '#AE76A3',
            '#882E72',
            '#1965B0',
                                    '#DC050C',
            '#4EB265',
            '#E8601C',
                        '#90C987',
                         '#D1BBD7',
                        '#7BAFDE' 


            )
named_colors2 <- setNames(colors2[1:length(clust2)], clust2)


p2 =SpatialDimPlot(vis_neuron%>%
                 subset(predicted.multiome_confidence_0.5%in%clust2), group.by = 'predicted.multiome_confidence_0.5',
               pt.size.factor = 0.5,
              images = 's_4P5.polygons',
             cols = named_colors2
               )
p2

ggsave(plot = p2,
       file = "4p5_multiome_4_5_17_18_22_10_16_23_0_24.tiff",
       device = "tiff",
       units = "in",
       width = 5,
       height = 5,
      path = "Manuscript/Plots/Manuscript v1.2.1/visium/",
      dpi=320)


clust3 = c(
  0,
4,
3,
5,
22,
17,
24,
16,
10,
19,
25,
7
)

colors3 <- c(
            '#F4A736',
            '#AE76A3',
            '#882E72',
            '#1965B0',
             '#DC050C',
            '#4EB265',
            '#E8601C',
             '#D1BBD7',
             '#90C987',
           '#7BAFDE' ,
            '#72190E',
           'darkgray'
            )
named_colors3 <- setNames(colors3[1:length(clust3)], clust3)


p3 =SpatialDimPlot(vis_neuron%>%
                 subset(predicted.multiome_confidence_0.5%in%clust3), group.by = 'predicted.multiome_confidence_0.5',
               pt.size.factor = 0.5,
              images = 's_4P10.polygons',
             cols = named_colors3
               )
p3

ggsave(plot = p3,
       file = "4p10_multiome_0_10_4_3_5_22_17_6_24_25_7.tiff",
       device = "tiff",
       units = "in",
       width = 5,
       height = 5,
      path = "Manuscript/Plots/Manuscript v1.2.1/visium/",
      dpi=320)

