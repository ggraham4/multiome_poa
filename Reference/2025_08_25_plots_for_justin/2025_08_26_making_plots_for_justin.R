###making plots for justin


library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(Polychrome)
library(Polychrome)
set.seed(935234)
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
names(P40) <- NULL

obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

setwd('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin')
#DimPlot

dim_labeled = DimPlot(obj, label = T, reduction = 'harmony_wnn.umap')+
  theme_void()+
  theme(legend.position="none")
dim_labeled

ggsave('dimplot_labeled.jpeg',
       dim_labeled, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )
setwd('/Users/ggraham/Desktop/multiome_poa')
gabe_labels = read.csv("Reference/2025_08_26_current_celltype_ids.csv")
gabe_labels$cluster = as.factor(gabe_labels$cluster)

labels = obj@meta.data%>%
  right_join(gabe_labels,by = join_by('final_clusters'=='cluster'))
obj$label = labels$label         

dim_named = DimPlot(obj, label = T, reduction = 'harmony_wnn.umap', group.by = 'label')+
  theme_void()+
  theme(legend.position="none")
dim_named
#not gonna happen

### DEG Plot ###
deg_file = read.csv("DEG Outputs/FINAL_degs_classified_08_11_2025.csv")
deg_file_summarized = deg_file%>%
  group_by(cluster, short_label)%>%
  summarize(n = n())
deg_file_summarized$signif = ifelse(deg_file_summarized$cluster %in% c(0,1,6,11), '*', NA)

missing_clusters = setdiff(0:26, deg_file_summarized$cluster)

deg_file_summarized_14 = data.frame(cluster = c(14, 18, 21), 
                                    short_label = rep('Early Downregulated',3),
                                    n = c(0,0,0),
                                    signif = c(NA, NA, NA))
deg_file_summarized = rbind(deg_file_summarized,deg_file_summarized_14 )

deg_file_total = deg_file%>%
  group_by(cluster)%>%
  summarize(n_total= n())

deg_file_summarized = deg_file_summarized%>%
  right_join(deg_file_total, by = 'cluster')

deg_plot = ggplot(deg_file_summarized, aes(x = as.factor(cluster), y = n, fill = short_label))+
  geom_bar(stat = 'identity', position = 'stack')+
  labs(x = 'Cluster', y = 'n_DEGs', fill = 'DEG Class')+
  scale_fill_manual(values =P40)+
  theme_minimal()+
  scale_y_continuous(breaks = c(0:40))+
  geom_text(aes(label =signif, y = n_total), size = 10)
deg_plot

"ggsave('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/degs_bar.jpeg',
       deg_plot, 
       device = 'jpeg',
       width = 20,
       height = 10,
       units = 'cm',
       dpi =600
       )
"

## 6 subclusters
obj_6 = FindSubCluster(obj, 6, graph.name ='harmony.wsnn')

labeled_6 = DimPlot(subset(obj_6, final_clusters ==6), group.by = 'sub.cluster', label =T, reduction = 'harmony_wnn.umap')+
  theme_void()+
  theme(legend.position="none")+
  labs(title = element_blank())
labeled_6

"ggsave('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/dimplot_6_subclusters.jpeg',
       labeled_6, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )
"
#interesting genes
for(i in c('ar',
           'esr2b',
           'pgr',
           'tacr3a',
           'drd3',
           'cckb',
           'cckbrb',
           'galr1a',
           'gal'
           
           )){
  
  plot = FeaturePlot(subset(obj_6, final_clusters ==6),feature = 'tacr3a', reduction = 'harmony_wnn.umap', pt.size = 0.5)+
      theme_void()+
    theme(title = element_blank())+
    labs(title = i)+
  theme(legend.position="none", plot.title =  element_text(hjust = .5))
  plot

 " ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/feautureplot_',i,'.jpeg'),
       plot, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )
"

}
#LOC111571064

  plot_gnrh = FeaturePlot(subset(obj_6, final_clusters ==6),feature = 'LOC111571064', reduction = 'harmony_wnn.umap', pt.size = 0.5)+
      theme_void()+
    theme(title = element_blank())+
    labs(title = 'gnrh1')+
  theme(legend.position="none", plot.title =  element_text(hjust = .5))
  plot_gnrh

  "ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/feautureplot_','gnrh1','.jpeg'),
       plot_gnrh, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )"

# tomorrow, I will do the same with MECP and also Proportion of cells in 6 expressing

mecp = function(object, gene, cluster, clustering = 'final_clusters', include_others = F){
    library(stringr)
    library(forcats)
      options(dplyr.summarise.inform = FALSE)

    counts <- t(as.matrix(object@assays$RNA$data[,object@meta.data[[clustering]] == cluster]))
  Counts_of_interest <- as.data.frame(counts[,gene])
    Counts_of_interest$expression <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
    Counts_of_interest$Sex <- object@meta.data$Sex[object@meta.data[[clustering]] == cluster]

  results <-Counts_of_interest%>%
    group_by(individual, Sex)%>%
    summarize(mean = mean(expression),
              se = sd(expression)/sqrt(n()))
    results$Sex= factor(results$Sex, levels = c('M', 'I','LI','NF', 'F'))


 if(include_others ==T){
    results= subset(results, Sex %in% c('M', 'I','LI','NF', 'F'))   
    jitter_width = 2.25


 }else{
   results= subset(results, Sex %in% c('M', 'I', 'F'))
   jitter_width = 1.75
}
 
  plot <- ggplot(results, aes(x = Sex, y = mean, color = Sex))+
    geom_boxplot(aes(group = Sex, fill = Sex), alpha = 0.25, outlier.shape = NA)+
    geom_pointrange(aes(x = Sex, y = mean, ymin = mean - se, ymax = mean+se), position = position_jitterdodge(jitter_width, seed = 1))+
    theme_classic()+
    labs(x  ='Status', y = 'Mean Normalized Expression +/- SE', title = paste0(gene,'_cluster_',cluster))
  return(plot)
}
mecd = readRDS('/Users/ggraham/Desktop/multiome_poa/Functions/mean_expression_cluster_data.RDS')

obj_6$Sex = ifelse(obj_6$Status == 'D', 'I', obj_6$Status)
obj_6$Sex = ifelse(obj_6$Sex == 'E', 'LI', obj_6$Sex)

mecp(obj_6, 'drd3' ,'6_0','sub.cluster', include_others = T)
mecp(obj_6, 'drd3' ,'6_0','sub.cluster', include_others = F)


#interesting DEGs
writer_function = function(gene, cluster,include_others, ...){
  
  plot = mecp(...)+
      labs(x = element_blank(), y = 'Mean +/- SE Normalized Counts', color = 'Status', fill = 'Status')

  if(include_others == T){
  path = paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/mean_expression_',gene,'_',cluster,'_others.jpeg')
  }else{
      path = paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/mean_expression_',gene,'_',cluster,'_mif.jpeg')

    
  }
  ggsave(path,
       plot, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )
  
}
writer_function('drd3', '6_0',T, obj_6, 'drd3', '6_0', 'sub.cluster',T)

writer_wrapper = function(object, gene, cluster, clustering = 'sub.cluster', include_others){
  writer_function(gene, cluster,include_others, object, gene, cluster, clustering,include_others)

}

c('ar',
           'esr2b',
           'pgr',
           'tacr3a',
           'drd3',
           'cckb',
           'cckbrb',
           'galr1a',
           'gal'
           
           )

"writer_wrapper(obj_6, 'pgr' ,'6','final_clusters',T)
writer_wrapper(obj_6, 'tacr3a' ,'6','final_clusters',T)
writer_wrapper(obj_6, 'esr2b' ,'6','final_clusters',T)
writer_wrapper(obj_6, 'ar' ,'6','final_clusters',T)
writer_wrapper(obj_6, 'drd3' ,'6','final_clusters',T)

writer_wrapper(obj_6, 'pgr' ,'6','final_clusters',F)
writer_wrapper(obj_6, 'tacr3a' ,'6','final_clusters',F)
writer_wrapper(obj_6, 'esr2b' ,'6','final_clusters',F)
writer_wrapper(obj_6, 'ar' ,'6','final_clusters',F)
writer_wrapper(obj_6, 'drd3' ,'6','final_clusters',F)"


#writer_wrapper(obj_6, 'tacr3a' ,'6_0','sub.cluster')
gnrh_6_3_others = mecp(obj_6, 'LOC111571064' ,'6_3','sub.cluster',T)+
    labs(title = 'gnrh1_cluster_6_3')



  "ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/mean_expression_gnrh1','_6_3_others','.jpeg'),
       gnrh_6_3_others, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )"
  
  gnrh_6_3_mif = mecp(obj_6, 'LOC111571064' ,'6_3','sub.cluster',F)+
        labs(title = 'gnrh1_cluster_6_3')


   " ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/mean_expression_gnrh1','_6_3_mif','.jpeg'),
       gnrh_6_3_mif, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )"
  

  gnrh_6_others = mecp(obj_6, 'LOC111571064' ,'6','final_clusters',T)+
  labs(title = 'gnrh1_cluster_6', x = element_blank())

 " ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/mean_expression_gnrh1','_6_others','.jpeg'),
       gnrh_6_others, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )"
  
    gnrh_6_mif = mecp(obj_6, 'LOC111571064' ,'6','final_clusters',F)+
  labs(title = 'gnrh1_cluster_6', x = element_blank())

"  ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/mean_expression_gnrh1','_6_mif','.jpeg'),
       gnrh_6_mif, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )
  
  
  writer_wrapper(obj_6, 'esr2b' ,'6_0','sub.cluster')

  writer_wrapper(obj_6, 'esr2b' ,'6_1','sub.cluster')
  writer_wrapper(obj_6, 'esr2b' ,'6_2','sub.cluster')
  writer_wrapper(obj_6, 'esr2b' ,'6_3','sub.cluster')
  
    writer_wrapper(obj_6, 'ar' ,'6_0','sub.cluster')

  writer_wrapper(obj_6, 'ar' ,'6_1','sub.cluster')
  writer_wrapper(obj_6, 'ar' ,'6_2','sub.cluster')
  writer_wrapper(obj_6, 'ar' ,'6_3','sub.cluster')
"


#proportions
obj_6$drd3 = ifelse(obj_6@assays$RNA$data['drd3',]>0, 'drd3', 'not')

dop_table = table(obj_6$Status, obj_6$drd3)%>%as.data.frame.matrix()
dop_table$prop = dop_table$drd3/(dop_table$drd3+dop_table$not)

  gene_prop_data = function(object, gene){
    if(!require(tidyverse)){library(tidyverse)}
    
    #classify cells as positive or negative for the gene
    object$pos_neg = ifelse(object@assays$RNA$data[gene,]>0, TRUE, FALSE)
    
    #make a table
    model_data = object@meta.data%>%
      group_by(individual, Sex)%>%
      subset(Sex %in% c('M','I', 'LI', 'NF', 'F'))%>%
      summarize(n_pos = sum(pos_neg==TRUE),
                n_neg = sum(pos_neg==FALSE))
    
    return(model_data)
  }
drd3_prop_data = gene_prop_data(subset(obj_6, final_clusters ==6), 'drd3')

drd3_prop_data$prop = drd3_prop_data$n_pos/(drd3_prop_data$n_pos+drd3_prop_data$n_neg)
drd3_prop_data$sum = (drd3_prop_data$n_pos+drd3_prop_data$n_neg)
drd3_prop_data$se = sqrt((drd3_prop_data$sum*drd3_prop_data$prop)*(1-drd3_prop_data$prop))

drd3_prop_data$Sex = factor(drd3_prop_data$Sex, levels = c('M','I','LI','NF','F'))

#plot
drd_prop_positive = ggplot(subset(drd3_prop_data, Sex %in% c('M','I','F')), aes(x = Sex, y = prop, color = Sex))+
  geom_jitter(position = position_jitterdodge(1), size =2)+
  geom_boxplot(alpha = 0.25, aes(fill = Sex), outlier.shape = NA)+
  labs(y = 'Proportion drd3+', title = 'Cluster 6', x = element_blank())+
  theme_classic()
   drd_prop_positive     
"
     ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/proportion_drd3_cluster6.jpeg'),
       drd_prop_positive, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )"
     
     drd_prop_positive2 = ggplot(subset(drd3_prop_data), aes(x = Sex, y = prop, color = Sex))+
  geom_jitter(position = position_jitterdodge(1), size =2)+
  geom_boxplot(alpha = 0.25, aes(fill = Sex), outlier.shape = NA)+
  labs(y = 'Proportion drd3+', title = 'Cluster 6', x = element_blank())+
  theme_classic()
   drd_prop_positive2   
   
        "ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/proportion_drd3_cluster6_others.jpeg'),
       drd_prop_positive2, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )
"

### proportion of cells
     
     sub_cells = obj_6@meta.data%>%
       subset(final_clusters==6)%>%
  group_by(sub.cluster, Sex, individual)%>%
  summarize(cells = n())

total_cells = obj_6@meta.data%>%
         subset(final_clusters==6)%>%
  group_by(Sex, individual)%>%
  summarize(total_cells = n())

full_data <- sub_cells%>%
  right_join(total_cells, by = 'individual')

full_data$diff = full_data$total_cells-full_data$cells

##plot
full_data$percent = full_data$cells/full_data$total_cells
full_data$Sex.x = factor(full_data$Sex.x, levels = c('NRM','M', 'I','LI','NF',"F"))

ggplot(full_data, aes(x = sub.cluster, y = percent, color = Sex.x))+
  geom_point(position = position_dodge(0.75))+
  geom_boxplot(alpha = 0)+
  theme_minimal()

proportion_6_3 = ggplot(subset(full_data, sub.cluster == '6_3' & Sex.x != 'NRM'), aes(x = Sex.x, y = percent, color = Sex.x, fill = Sex.x))+
  geom_point(position = position_jitterdodge(1), size =2)+
  geom_boxplot(alpha = 0.25)+
  theme_classic()+
  labs(x = element_blank(), y = 'Proportion of Nuclei',color = 'Sex', fill = 'Sex',
       title = 'Proportion of Cluster 6 Nuclei in 6_3')
proportion_6_3

     "ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/proportion_of_cells_in_6_3.jpeg'),
       proportion_6_3, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )"

     proportion_6_3_mif = ggplot(subset(full_data, sub.cluster == '6_3'& Sex.x %in% c('M','I','F')), aes(x = Sex.x, y = percent, color = Sex.x, fill = Sex.x))+
  geom_point(position = position_jitterdodge(1), size =2)+
  geom_boxplot(alpha = 0.25)+
  theme_classic()+
  labs(x = element_blank(), y = 'Proportion of Nuclei',color = 'Sex', fill = 'Sex',
       title = 'Proportion of Cluster 6 Nuclei in 6_3')
proportion_6_3_mif

     "ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/proportion_of_cells_in_6_3_mif.jpeg'),
       proportion_6_3_mif, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )"

     
     
obj$color = ifelse(obj$final_clusters %in% c(6,0,21,17,19), 'POA', NA)
obj$color = ifelse(obj$final_clusters %in% c(1,26), 'RG', obj$color)
obj$color = ifelse(obj$final_clusters %in% c(11), 'MG', obj$color)
obj$color = ifelse((is.na(obj$color) & obj$final_clusters %in% c(4,5,10,7,8,12,16,18,22)), 'GLUT', obj$color)
obj$color = ifelse((is.na(obj$color) & obj$final_clusters %in% c(3,23,25,14)), 'GABA', obj$color)
obj$color = ifelse((is.na(obj$color) & obj$final_clusters %in% c(9,24)), 'Mixed', obj$color)
obj$color = ifelse((is.na(obj$color) & obj$final_clusters %in% c(2)), 'OG', obj$color)
obj$color = ifelse((is.na(obj$color) & obj$final_clusters %in% c(13)), 'OPC', obj$color)
obj$color = ifelse((is.na(obj$color) & obj$final_clusters %in% c(15)), 'EG', obj$color)
obj$color = ifelse((is.na(obj$color) & obj$final_clusters %in% c(20)), 'Leuko', obj$color)

# Define cell type colors
cell_type_colors <- c(
  'POA' = '#FF6B6B',     # Red
  'RG' = '#4ECDC4',      # Teal
  'MG' = '#45B7D1',      # Blue
  'GLUT' = '#96CEB4',    # Green
  'GABA' = '#FFEAA7',    # Yellow
  'Mixed' = '#DDA0DD',   # Plum
  'OG' = '#FFB347',      # Orange
  'OPC' = '#98D8C8',     # Mint
  'EG' = '#F7DC6F',      # Light yellow
  'Leuko' = '#BB8FCE'    # Light purple
)

# Create mapping from cluster to color
cluster_to_color <- obj@meta.data %>%
  select(final_clusters, color) %>%
  distinct() %>%
  arrange(final_clusters) %>%
  mutate(actual_color = cell_type_colors[color])

# Create named vector of colors for each cluster (this is key!)
cluster_colors <- setNames(cluster_to_color$actual_color, cluster_to_color$final_clusters)

# Plot with cluster numbers as labels but colored by cell type
p <- DimPlot(obj, 
        reduction = 'harmony_wnn.umap', 
        group.by = 'final_clusters', 
        cols = cluster_colors,
        label = TRUE) +
    theme_void()+
  theme(legend.position = "none", plot.title = element_blank(), legend = element_blank())
p

# Create custom legend data
legend_data <- data.frame(
  cell_type = names(cell_type_colors),
  color = cell_type_colors,
  x = 1,
  y = seq_along(cell_type_colors)
)

# Create the legend plot
legend_plot <- ggplot(legend_data, aes(x = x, y = y, fill = cell_type)) +
  geom_point(size = 5, shape = 21) +
  scale_fill_manual(values = cell_type_colors) +
  geom_text(aes(x = 1.3, label = cell_type), hjust = 0, size = 3) +  # Fixed x position
  theme_void() +
  theme(legend.position = "none") +
  xlim(0.8, 2.5) +
  ggtitle("Cell Types")

# Combine plots
library(patchwork)
combined_plot <- p + legend_plot + plot_layout(widths = c(4, 1))

combined_plot

     ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/dimplot_major_cells_2.jpeg'),
       combined_plot, 
       device = 'jpeg',
       width = 15,
       height = 10,
       units = 'cm',
       dpi =600
       )


main_map = DimPlot(obj, reduction = 'harmony_wnn.umap', group.by = 'color', label = T)+
   theme_void()+
  theme(legend.position="none", title = element_blank())
main_map


     "ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/dimplot_major_cells_1.jpeg'),
       main_map, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )
     "


    # main_map2 = 
  DimPlot(obj, reduction = 'harmony_wnn.umap', group.by = 'final_clusters', cols = color_df$actual_color)
   theme_void()+
  theme(legend.position="none", title = element_blank())
main_map2

     "ggsave(paste0('/Users/ggraham/Desktop/multiome_poa/Reference/2025_08_25_plots_for_justin/dimplot_major_cells_2.jpeg'),
       main_map2, 
       device = 'jpeg',
       width = 10,
       height = 10,
       units = 'cm',
       dpi =600
       )"
     

