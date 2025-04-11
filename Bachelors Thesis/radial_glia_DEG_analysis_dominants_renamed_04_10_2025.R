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
for (i in 0:9) {
  cluster <- paste0('2_',i)
  print(cluster)
  output <- neg_bin_mult_windows(obj = radial_glia,
                                 clustering = 'sub',
                                 cluster = cluster)
  output$cluster = cluster
  neg_bin <- rbind(neg_bin, output)
}

neg_bin_defined <- define_degs(neg_bin)

neg_bin_defined$issignif <- ifelse(!is.na(neg_bin_defined$class), '*',NA)

#write.csv(neg_bin_defined, 'DEG Outputs/subclusters_2_03_23_2025.csv')

neg_bin_defined_filtered <- neg_bin_defined%>%
  filter(!is.na(issignif)& is.na(warning))

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

ggsave(plot = deg_counts,
       file = "deg_counts_2.svg",
       device = "svg",
       units = "in",
       width = 4,
       height = 3,
       path = "Bachelors Thesis/Plots/Figure 5")



#overall
clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='2_2'])%>%dotplot()## cytokine mediated??
clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='2_4'])%>%dotplot() ### important
clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='2_4'])$geneID

#late down
#clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='2_0' & neg_bin_defined_filtered$class =='Late Downregulated'])%>%dotplot() ### important
clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='2_1' & neg_bin_defined_filtered$class =='Late Downregulated'])%>%dotplot() ### important
#### oooooooooooohhh locked innnnnn

#early down
clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='2_1' & neg_bin_defined_filtered$class =='Early Downregulated'])%>%dotplot() ### important

clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='2_2' & neg_bin_defined_filtered$class =='Terminally Downregulated'])%>%dotplot() ### important

clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='2_2'])$geneID


## save plots

late_down_2_1 = clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='2_1' & neg_bin_defined_filtered$class =='Late Downregulated'])%>%
  dotplot()+
  theme(legend.position = 'none')+
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = -45),
        legend.position = 'none')+
  scale_size_continuous(range = c(1,3))
### important

ggsave(plot = late_down_2_1,
       file = "late_down_2_1.svg",
       device = "svg",
       units = "in",
       width = 6,
       height = 2.3,
       path = "Bachelors Thesis/Plots/Figure 5")


ggsave(plot = late_down_2_1,
       file = "late_down_2_1_short.svg",
       device = "svg",
       units = "in",
       width = 4,
       height = 2.3,
       path = "Bachelors Thesis/Plots/Figure 5")


all_2_4 = clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='2_4'])%>%dotplot() +### important
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = -45),
        legend.position = 'none')+
  scale_size_continuous(range = c(1,3))
all_2_4


ggsave(plot = all_2_4,
       file = "all_2_4_short.svg",
       device = "svg",
       units = "in",
       width = 4,
       height = 2.3,
       path = "Bachelors Thesis/Plots/Figure 5")



###plot brain aromatase DEG in 2_1
radial_glia_data <- radial_glia@assays$RNA$data[,radial_glia$sub =='2_1']%>%
  t()%>%
  as.data.frame()%>%
  dplyr::select('LOC111577263')

aromatase_data <- data.frame(aromatase_expression = radial_glia_data,
                             individual = radial_glia$individual[radial_glia$sub =='2_1'],
                             Sex = radial_glia$Status[radial_glia$sub =='2_1'])%>%
  group_by(individual, Sex)%>%
  summarize(mean_expression = mean(LOC111577263),
            se = sd(LOC111577263)/sqrt(n()))%>%
  subset(Sex != 'NRM')
library(ggsignif)


aromatase_data$Sex <- ifelse(aromatase_data$Sex =='D','I',aromatase_data$Sex)
aromatase_data$Sex <- ifelse(aromatase_data$Sex =='E','LI',aromatase_data$Sex)

aromatase_data$Sex <- factor(aromatase_data$Sex, levels = c('M','I','LI',"NF",'F'))


aromatase_plot <- ggplot(aromatase_data, aes(x = Sex, y = mean_expression))+
  geom_boxplot(alpha = 0.25, outlier.shape = NA, aes(color = Sex, fill = Sex))+
  geom_jitter(  shape = 1, color = 'black', size =2)+
  theme_classic()+
  theme(legend.position ='none',plot.title = element_text(hjust=0.5))+
  labs(title = '2_1', y = 'Mean cyp19a1b')+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(4), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(4), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)+
  
  ylim(1,4.5)
aromatase_plot

ggsave(plot = aromatase_plot,
       file = "aromatase_expression_2_1.svg",
       device = "svg",
       units = "in",
       width = 1.6,
       height = 2,
       path = "Bachelors Thesis/Plots_dominants_renamed/Figure 5")

#### in 2_3
radial_glia_data_2_3 <- radial_glia@assays$RNA$data[,radial_glia$sub =='2_3']%>%
  t()%>%
  as.data.frame()%>%
  dplyr::select('LOC111577263')

aromatase_data_2_3 <- data.frame(aromatase_expression = radial_glia_data_2_3,
                             individual = radial_glia$individual[radial_glia$sub =='2_3'],
                             Sex = radial_glia$Status[radial_glia$sub =='2_3'])%>%
  group_by(individual, Sex)%>%
  summarize(mean_expression = mean(LOC111577263),
            se = sd(LOC111577263)/sqrt(n()))%>%
  subset(Sex != 'NRM')
library(ggsignif)

aromatase_data_2_3$Sex <- ifelse(aromatase_data_2_3$Sex =='D','I',aromatase_data_2_3$Sex)
aromatase_data_2_3$Sex <- ifelse(aromatase_data_2_3$Sex =='E','LI',aromatase_data_2_3$Sex)

aromatase_data_2_3$Sex <- factor(aromatase_data_2_3$Sex, levels = c('M','I','LI',"NF",'F'))


aromatase_data_2_3_plot <- ggplot(aromatase_data_2_3, aes(x = Sex, y = mean_expression))+
  geom_boxplot(alpha = 0.25, outlier.shape = NA, aes(color = Sex, fill = Sex))+
  geom_jitter(  shape = 1, color = 'black', size =2)+
  theme_classic()+
  theme(legend.position ='none',plot.title = element_text(hjust=0.5))+
  labs(title = '2_3', y = 'Mean cyp19a1b')+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(4), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(4), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  theme(axis.text.y= element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())+
  ylim(1,4.5)
aromatase_data_2_3_plot

ggsave(plot = aromatase_data_2_3_plot,
       file = "aromatase_expression_2_3.svg",
       device = "svg",
       units = "in",
       width = 1.2,
       height = 2,
       path = "Bachelors Thesis/Plots_dominants_renamed/Figure 5")


#### in 2_7
radial_glia_data_2_7 <- radial_glia@assays$RNA$data[,radial_glia$sub =='2_7']%>%
  t()%>%
  as.data.frame()%>%
  dplyr::select('LOC111577263')

aromatase_data_2_7 <- data.frame(aromatase_expression = radial_glia_data_2_7,
                                 individual = radial_glia$individual[radial_glia$sub =='2_7'],
                                 Sex = radial_glia$Status[radial_glia$sub =='2_7'])%>%
  group_by(individual, Sex)%>%
  summarize(mean_expression = mean(LOC111577263),
            se = sd(LOC111577263)/sqrt(n()))%>%
  subset(Sex != 'NRM')
library(ggsignif)

aromatase_data_2_7$Sex <- ifelse(aromatase_data_2_7$Sex =='D','I',aromatase_data_2_7$Sex)
aromatase_data_2_7$Sex <- ifelse(aromatase_data_2_7$Sex =='E','LI',aromatase_data_2_7$Sex)

aromatase_data_2_7$Sex <- factor(aromatase_data_2_7$Sex, levels = c('M','I','LI',"NF",'F'))

aromatase_data_2_7_plot <- ggplot(aromatase_data_2_7, aes(x = Sex, y = mean_expression))+
  geom_boxplot(alpha = 0.25, outlier.shape = NA, aes(color = Sex, fill = Sex))+
  geom_jitter(  shape = 1, color = 'black', size =2)+
  theme_classic()+
  theme(legend.position ='none',plot.title = element_text(hjust=0.5))+
  labs(title = '2_7', y = 'Mean cyp19a1b')+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(1.7), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=6)+
  theme(
        axis.title.y = element_blank())+
  ylim(0,2)
aromatase_data_2_7_plot

ggsave(plot = aromatase_data_2_7_plot,
       file = "aromatase_expression_2_7.svg",
       device = "svg",
       units = "in",
       width = 1.35,
       height = 2,
       path = "Bachelors Thesis/Plots_dominants_renamed/Figure 5")


