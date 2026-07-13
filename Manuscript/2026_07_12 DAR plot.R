library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)


colors = c("Early Downregulated" = "red",
           "Early Upregulated" = "#006400",
           "Late Downregulated" = "blue",
           "Late Upregulated"='#000000',
           "Progressively Downregulated" = 'purple',
           "Progressively Upregulated"='gray',
           "Transiently Downregulated"='brown',
           "Transiently Upregulated"='orange')

#Dars
dars =read.csv("Collaboration/all_clusters_DARs_peak_level_classified_with_support.csv")

plot_dars = dars %>%
  group_by(cluster_id, classification) %>%
  summarize(n = n(), .groups = "drop")

plot_dars$color = colors[plot_dars$classification] %>% unname()

plot_dars$cluster_id_num = as.numeric(as.character(plot_dars$cluster_id))

# add a spacer level before 0, using a level that has no matching data
plot_dars$cluster_id_num = factor(plot_dars$cluster_id_num, 
                                  levels = c(-1, 0:26))

dar_plot = ggplot(plot_dars, aes(x = cluster_id_num, y = n, fill = classification)) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  scale_fill_manual(values = colors) +
  labs(x = 'Cluster', y = 'DARs', fill = 'Classification') +
  scale_x_discrete(drop = FALSE, labels = c("", 0:26)) +  # hide the -1 label, show blank
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dar_plot

ggsave(plot = dar_plot,
      file = "DARs.svg",
     device = "svg",
    units = "in",
   width = 8.5,
  height = 2.5,
 path = "Manuscript/Plots/Manuscript v1.3/")
