library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)

vis2 = readRDS("C:/Users/Gabe/Desktop/Visium/vis_better_barcodes_dissection_c2l05_projection_anatomical.rds")

heatmat = table(vis$Predicted_final_clusters,
                vis$Anatomical)%>%
  as.data.frame.matrix()%>%
  dplyr::select(!c('not sure 2',
           'not sure 1',
           'ventral diffuse',
           'No Region',
           'not sure 4',
           'not sure 3',
           'Diffuse',
           'What'))

# Define anchor colors
my_colors <- c("gray",'white', "orange", "red")

# Create a function, then evaluate for n = 50 colors
color_func <- colorRampPalette(my_colors)
continuous_palette <- color_func(50)

p =pheatmap((heatmat/rowSums(heatmat))%>%scale(),
         cluster_rows = T,
         cluster_cols = T,
         treeheight_row=0,
         treeheight_col = 0,
         color =continuous_palette,
         angle_col =90,
         number_color= 'black')

# Open an SVG device, draw, then close
svg(filename = 'C:/Users/Gabe/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2.3/region heatmap.svg',
    width = 6.5,
    height = 4.5)
p
dev.off()


