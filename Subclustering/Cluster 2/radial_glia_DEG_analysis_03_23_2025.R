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
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(Polychrome)
  library(scCustomize)
  library(hdWGCNA)
  
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
  define_degs<- readRDS('Functions/define_degs')
  
}

radial_glia <- readRDS(file='A:/anemonefish_multiome_radial_glia_03_20_2025/radial_glia_object.rds')




### DEGS #####
neg_bin_mult_windows<- readRDS('Functions/DEG_functions/neg_bin_mult_windows.rds')

neg_bin <- data.frame()
for (i in 7:9) {
  cluster <- paste0('2_',i)
  print(cluster)
  output <- neg_bin_mult_windows(obj = radial_glia,
                                 clustering = 'sub',
                                 cluster = cluster)
  output$cluster = cluster
  neg_bin <- rbind(neg_bin, output)
}

neg_bin_defined <- define_degs(neg_bin)

neg_bin_defined_filtered <- neg_bin_defined%>%
  filter(!is.na(issignif)& is.na(warning))

write.csv(neg_bin_defined_filtered, 'DEG Outputs/subclusters_2_03_23_2025.csv')

neg_bin_defined_counts<- neg_bin_defined_filtered%>%
  group_by(cluster, class)%>%
  summarize(class_count = n())

neg_bin_defined_counts$colors <- NA
neg_bin_defined_counts$colors <- ifelse(neg_bin_defined_counts$class == 'Early Upregulated', 'blue', NA)
neg_bin_defined_counts$colors <- ifelse(neg_bin_defined_counts$class == 'Early Downregulated', 'red', neg_bin_defined_counts$colors )
neg_bin_defined_counts$colors <- ifelse(neg_bin_defined_counts$class == 'Late Downregulated', 'pink', neg_bin_defined_counts$colors )
neg_bin_defined_counts$colors <- ifelse(neg_bin_defined_counts$class == 'Late Upregulated', 'orange', neg_bin_defined_counts$colors )
neg_bin_defined_counts$colors <- ifelse(neg_bin_defined_counts$class == 'Progressively Upregulated', 'cyan', neg_bin_defined_counts$colors )
neg_bin_defined_counts$colors <- ifelse(neg_bin_defined_counts$class == 'Progressively Downregulated', 'hotpink', neg_bin_defined_counts$colors )
neg_bin_defined_counts$colors <- ifelse(neg_bin_defined_counts$class == 'Terminally Downregulated', 'maroon', neg_bin_defined_counts$colors )
neg_bin_defined_counts$colors <- ifelse(neg_bin_defined_counts$class == 'Terminally Upregulated', 'darkgreen', neg_bin_defined_counts$colors )
neg_bin_defined_counts$colors <- ifelse(neg_bin_defined_counts$class == 'Transiently Downregulated', 'yellow', neg_bin_defined_counts$colors )
neg_bin_defined_counts$colors <- ifelse(neg_bin_defined_counts$class == 'Transiently Upregulated', 'green', neg_bin_defined_counts$colors )


neg_bin_defined_counts$class <- factor(neg_bin_defined_counts$class, 
                                       levels = rev(c('Transiently Upregulated',
                                                      'Transiently Downregulated',
                                                      'Terminally Upregulated',
                                                      "Terminally Downregulated",
                                                      'Progressively Upregulated',
                                                      'Progressively Downregulated',
                                                      'Late Upregulated',
                                                      'Late Downregulated',
                                                      'Early Upregulated',
                                                      'Early Downregulated')))

class_colors <- c(
  'Transiently Upregulated' = 'green',
  'Transiently Downregulated' = 'yellow',
  'Terminally Upregulated' = 'darkgreen',
  'Terminally Downregulated' = 'maroon',
  'Progressively Upregulated' = 'cyan',
  'Progressively Downregulated' = 'hotpink',
  'Late Upregulated' = 'orange',
  'Late Downregulated' = 'pink',
  'Early Upregulated' = 'blue',
  'Early Downregulated' = 'red'
)

# Now update your plot
deg_counts <- ggplot(neg_bin_defined_counts, aes(x = as.factor(cluster), y = class_count, fill = class)) +
  labs(x = "Subcluster", y = "Number of DEGs", fill = 'Expression Pattern') +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0)) +
  scale_fill_manual(values = class_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 0.4, angle = -90),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.background = element_rect(color = 'white'),
        legend.direction = 'vertical',
        legend.title.position = 'top',
        legend.position = 'right'
  ) +
  theme(axis.title = element_text(size = 14),
        axis.text =element_text(size = 12),
        axis.title.y = element_text(hjust =1)
  )+
  ylim(0, 11)
deg_counts
