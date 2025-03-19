#Reclustering cluster 19 
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
  library('glmGamPoi')
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
  
  define_degs <- function(data, singular = TRUE) {
    if (!singular) {
      data <- data[data$singular == FALSE, ]
    }
    
    # Assign classes based on conditions
    data$class <- NA  # Initialize class column
    
    data$class[data$d_m_q.value < 0.05 & 
                 data$d_m_estimate > 0 & 
                 data$d_f_q.value > 0.05] <- 'Early Upregulated'
    
    
    data$class[data$d_m_q.value > 0.05 & 
                 data$d_f_q.value < 0.05 & 
                 data$d_f_estimate < 0] <- 'Late Upregulated'
    
    data$class[data$d_m_q.value < 0.05 & 
                 data$d_m_estimate < 0 & 
                 data$d_f_q.value > 0.05] <- 'Early Downregulated'
    
    
    data$class[data$d_m_q.value > 0.05 & 
                 data$d_f_q.value < 0.05 & 
                 data$d_f_estimate > 0] <- 'Late Downregulated'
    
    data$class[data$d_m_q.value < 0.05 & 
                 data$d_f_q.value < 0.05 & 
                 data$d_f_estimate > 0 & 
                 data$d_m_estimate > 0] <- 'Transiently Upregulated'
    
    data$class[data$d_m_q.value < 0.05 & 
                 data$d_f_q.value < 0.05 & 
                 data$d_f_estimate < 0 & 
                 data$d_m_estimate < 0] <- 'Transiently Downregulated'
    
    data$class[data$d_m_q.value < 0.05 & 
                 data$d_f_q.value < 0.05 & 
                 data$d_m_estimate < 0 & 
                 data$d_f_estimate > 0] <- 'Progressively Downregulated'
    
    data$class[data$d_m_q.value < 0.05 & 
                 data$d_f_q.value < 0.05 & 
                 data$d_m_estimate > 0 & 
                 data$d_f_estimate < 0] <- 'Progressively Upregulated'
    
    data$class[data$f_m_q.value < 0.05 & 
                 data$d_f_q.value > 0.05 & 
                 data$d_m_q.value > 0.05 & 
                 data$f_m_estimate > 0 ] <- 'Terminally Upregulated'
    
    data$class[data$f_m_q.value < 0.05 & 
                 data$d_f_q.value > 0.05 & 
                 data$d_m_q.value > 0.05 & 
                 data$f_m_estimate < 0  ] <- 'Terminally Downregulated'
    
    
    data$issignif <- NA
    data$issignif <- ifelse(data$f_m_q.value<0.05|
                              data$d_m_q.value<0.05|
                              data$d_f_q.value<0.05,
                            '*',NA)
    
    return(data)
  }
  
  mean_expression_cluster_data <- function(object, gene, cluster, clustering = 'harmony.wnn_res0.4_clusters'){
    counts <- t(object@assays$RNA$counts[,object@meta.data[[clustering]] == cluster])
    Counts_of_interest <- as.data.frame(counts[,gene])
    Counts_of_interest[[gene]] <- Counts_of_interest[,1]
    Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
    results <- data.frame()
    
    for (i in unique(object@meta.data$individual)) {
      Counts <- Counts_of_interest[[gene]][Counts_of_interest$individual==i]
      mean <- mean(Counts)
      mean_se <- sd(Counts) / sqrt(length(Counts))
      df <- data.frame(
        individual = i,
        mean = mean,
        se = mean_se
      )
      results <- rbind(results, df)
    }
    results$Sex <- str_sub(results$individual, -1)
    results$Sex[results$individual == 'T17D'] = 'NF'
    results$Sex[results$individual == 'A12D'] = 'E'
    results$Sex[results$individual == 'T11D'] = 'E'
    results$Sex[results$individual == 'GH'] = 'NRM'
    return(results)
  }
  
  
  
}

obj <- readRDS('C:/Users/Gabe/Desktop/RNA Object.rds')

cluster_19 <- subset(obj, harmony.wnn_res0.4_clusters ==19)
cluster_19 <- FindSubCluster(cluster_19, 19, 'harmony.wsnn',  subcluster.name= 'sub')

Idents(cluster_19) <- 'sub'

arr <- list(x = -3.5, y = -6, x_len = 2, y_len = 2)

dimplot_19 <- DimPlot(cluster_19, label = T, reduction = 'harmony_wnn.umap')+
  theme_void()+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), 
           arrow = arrow(type = "closed", length = unit(10, 'pt'))) +
  theme(legend.position = 'none')+
  annotate('text',
           x = -4.3, y = -5.2, label = 'UMAP_2', angle = 90, size =3)+
  annotate('text',
           x = -3, y = -6.5, label = 'UMAP_1', size =3)
dimplot_19

ggsave(plot = dimplot_19,
       file = "dimplot_19.svg",
       device = "svg",
       units = "in",
       width = 7,
       height = 7,
       path = "Bachelors Thesis/Plots/Figure 4")



### What are the differences between the subclusters ####
FeaturePlot(cluster_19, c('gad2','slc17a6b'))

gaba_glut_markers <- DotPlot(object = cluster_19, 
        group.by = "sub", 
        features = c('LOC111588076','gad2', 'LOC111584103','slc17a6b','slc17a7a'),
        dot.min = 0.1,
        col.min =-2.3)  +
  scale_x_discrete(labels = c('gad1',
                   'gad2',
                   'vglut2.1',
                   'slc17a6b',
                   'slc17a7a'))+
  coord_flip()+
  theme(axis.text.x = element_text(angle = -90), legend.position = 'none', axis.title.y = element_blank())+
  labs(y= 'Subcluster')

ggsave(plot = gaba_glut_markers,
       file = "gaba_glut_markers.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 4")


#split cells into glut, gaba , and mixed
cluster_19$primary_neurotransmitter <- ifelse(((cluster_19@assays$RNA$data['LOC111588076',]>0 
                                               | cluster_19@assays$RNA$data['gad2',]>0) &
                                                (cluster_19@assays$RNA$data['LOC111584103',]==0 & 
                                                   cluster_19@assays$RNA$data['slc17a6b',]==0 &  
                                                   cluster_19@assays$RNA$data['slc17a7a',]==0)), 'GABA',NA)

cluster_19$primary_neurotransmitter <- ifelse(((cluster_19@assays$RNA$data['LOC111588076',]==0 
                                                & cluster_19@assays$RNA$data['gad2',]==0) &
                                                 (cluster_19@assays$RNA$data['LOC111584103',]>0 | 
                                                    cluster_19@assays$RNA$data['slc17a6b',]>0 |  
                                                    cluster_19@assays$RNA$data['slc17a7a',]>0)), 'GLUT',cluster_19$primary_neurotransmitter)

ccluster_19$primary_neurotransmitter <- ifelse(((cluster_19@assays$RNA$data['LOC111588076',]>0 
                                                | cluster_19@assays$RNA$data['gad2',]>0) &
                                                 (cluster_19@assays$RNA$data['LOC111584103',]>0 | 
                                                    cluster_19@assays$RNA$data['slc17a6b',]>0 |  
                                                    cluster_19@assays$RNA$data['slc17a7a',]>0)), 'Mixed',cluster_19$primary_neurotransmitter)

cluster_19$primary_neurotransmitter <- ifelse(is.na(cluster_19$primary_neurotransmitter), 'Neither', cluster_19$primary_neurotransmitter)

table(cluster_19$primary_neurotransmitter)

 prim_neuro <-DimPlot(cluster_19, label = F, reduction = 'harmony_wnn.umap', group.by ='primary_neurotransmitter')+
   theme_void()+
   annotate("segment", 
            x = arr$x, xend = arr$x + c(arr$x_len, 0), 
            y = arr$y, yend = arr$y + c(0, arr$y_len), 
            arrow = arrow(type = "closed", length = unit(10, 'pt'))) +
   theme(legend.position=c(.5,0.3))+
   annotate('text',
            x = -4.3, y = -5.2, label = 'UMAP_2', angle = 90, size =3)+
   annotate('text',
            x = -3, y = -6.5, label = 'UMAP_1', size =3)
 
 ggsave(plot = prim_neuro,
        file = "prim_neuro.svg",
        device = "svg",
        units = "in",
        width = 7,
        height = 7,
        path = "Bachelors Thesis/Plots/Figure 4")
 
 ##markers by primary neurotransmitter type
 Idents(cluster_19) <- 'primary_neurotransmitter'
 markers_prim_type <- FindAllMarkers(cluster_19 )

 markers_mixed <-   markers_prim_type$gene[markers_prim_type$cluster=='Mixed'&markers_prim_type$p_val_adj<0.05 & markers_prim_type$pct.1>markers_prim_type$pct.2]
 markers_deficient <-   markers_prim_type$gene[markers_prim_type$cluster=='Mixed'&markers_prim_type$p_val_adj<0.05 & markers_prim_type$pct.1<markers_prim_type$pct.2]
 
 
 clown_go<- readRDS('Functions/clown_go')
clown_go(markers_mixed)%>%dotplot()

markers_GABA <-   markers_prim_type$gene[markers_prim_type$cluster=='GABA' & markers_prim_type$p_val_adj<0.05]

clown_go(markers_GABA)%>%dotplot()

markers_GLUT <-   markers_prim_type$gene[markers_prim_type$cluster=='GLUT'& markers_prim_type$p_val_adj<0.05] 

clown_go(markers_GLUT)%>%dotplot()

markers_neither <-   markers_prim_type[markers_prim_type$cluster=='Neither',] #top marker is amyloid beta??? what on earth
markers_neither$gene[markers_neither$p_val_adj<0.05 & markers_neither$pct.1>markers_neither$pct.2] # really no good markers weird

clown_go(markers_neither$gene[markers_neither$p_val_adj<0.05])%>%dotplot() #this is really something huh, the top markers all seem related to metabolism but they are all goig to be not expressed in this celltype
# my primary theory here is these are normal cells and just sequencing didn't happent to pickup a glut or gaba transcript

### cytotrace
cluster_19_matrix <- as.matrix(cluster_19@assays$RNA$counts) #not necessary to normalize
cluster_19_cyto <- CytoTRACE(mat = cluster_19_matrix
)

cluster_19$cyto <-cluster_19_cyto$CytoTRACE

cluster_19_data <- data.frame(individual = cluster_19@meta.data$individual,
                              status = cluster_19@meta.data$Status,
                              cluster = cluster_19@meta.data$sub,
                              cyto = cluster_19@meta.data$cyto)%>%
  subset(status != 'NRM')
cluster_19_data$status <- factor(cluster_19_data$status, levels = c('M','D',"E",'NF','F'))

 ggplot(cluster_19_data, aes(x = fct_reorder(prim_neuro, cyto, .desc = T), y = cyto, group = interaction(prim_neuro, status), color = status))+
  geom_boxplot(aes(fill = status), alpha = 0.25)

 cell_counts <- table(cluster_19$primary_neurotransmitter, cluster_19$Status)%>%as.data.frame()
 cell_counts$prop[cell_counts$Var2 == 'D'] = cell_counts$Freq[cell_counts$Var2 == 'D']/ nrow(cluster_19@meta.data[cluster_19@meta.data$Status == 'D',])
 cell_counts$prop[cell_counts$Var2 == 'M'] = cell_counts$Freq[cell_counts$Var2 == 'M']/ nrow(cluster_19@meta.data[cluster_19@meta.data$Status == 'M',])
 cell_counts$prop[cell_counts$Var2 == 'F'] = cell_counts$Freq[cell_counts$Var2 == 'F']/ nrow(cluster_19@meta.data[cluster_19@meta.data$Status == 'F',])
 cell_counts$prop[cell_counts$Var2 == 'E'] = cell_counts$Freq[cell_counts$Var2 == 'E']/ nrow(cluster_19@meta.data[cluster_19@meta.data$Status == 'E',])
 cell_counts$prop[cell_counts$Var2 == 'NF'] = cell_counts$Freq[cell_counts$Var2 == 'NF']/ nrow(cluster_19@meta.data[cluster_19@meta.data$Status == 'NF',])
 
 cell_counts$Var2 <- factor(cell_counts$Var2, levels = c('M','D',"E",'NF','F'))
 
 ggplot(cell_counts, aes(x = Var2, y = prop, group = interaction(Var2, Var1), fill = Var1))+
   geom_bar(stat='identity')
 
 cluster_prim_type <- table(cluster_19$primary_neurotransmitter, cluster_19$sub)%>%as.data.frame()
 cluster_prim_type$prop[cluster_prim_type$Var2 == '19_0'] = cluster_prim_type$Freq[cluster_prim_type$Var2 == '19_0']/ nrow(cluster_19@meta.data[cluster_19@meta.data$sub == '19_0',])
 cluster_prim_type$prop[cluster_prim_type$Var2 == '19_1'] = cluster_prim_type$Freq[cluster_prim_type$Var2 == '19_1']/ nrow(cluster_19@meta.data[cluster_19@meta.data$sub == '19_1',])
 cluster_prim_type$prop[cluster_prim_type$Var2 == '19_2'] = cluster_prim_type$Freq[cluster_prim_type$Var2 == '19_2']/ nrow(cluster_19@meta.data[cluster_19@meta.data$sub == '19_2',])
 
 subcluster_comp <- ggplot(cluster_prim_type, aes(x = Var2, y = prop, group = interaction(Var2, Var1), fill = Var1))+
   geom_bar(stat='identity')+
   labs(x = 'Subcluster', y = 'Proportion')+
   theme_classic()+
   theme(legend.position = 'none', axis.title = element_text(size = 14), axis.text =element_text(size = 12), axis.text.x = element_text(angle =-90))
 subcluster_comp
 
 ggsave(plot = subcluster_comp,
        file = "subcluster_comp.svg",
        device = "svg",
        units = "in",
        width = 1.65,
        height = 2,
        path = "Bachelors Thesis/Plots/Figure 4")
 
 ggplot(cluster_19_data, aes(x = cluster, y = cyto, fill = prim_neuro, group = prim_neuro))+
   stat_summary(geom = 'bar', fun = 'mean', position  ='dodge')+
   stat_summary(geom = 'errorbar', position  ='dodge', color = 'black')
 
 
 linear_model <- lmer(cyto~prim_neuro + cluster+ (1|individual), data = subset(cluster_19_data, status %in% c('M','D','F')))
 car::Anova(linear_model)
 pairs(emmeans(linear_model, 'prim_neuro', by = 'cluster'))
 #Alright I think more of the cytotrace score is explained by subcluster than NT type, so I am going to move on from this idea.
 
 linear_model_subcluster <- lmer(cyto~ status*cluster+ (1|individual), data = subset(cluster_19_data, status %in% c('M','D','F')))
 car::Anova(linear_model_subcluster, type ='III')

 
  pairs(emmeans(linear_model_subcluster, 'status', by = 'cluster'), adjust ='none')
  #using unadjusted because tukey is v conservative
 
  pairs(emmeans(linear_model_subcluster, 'cluster'), adjust ='none')
  
  data_cyto_plot <- subset(cluster_19_data, status %in% c('M','D','F'))%>%
    na.omit()%>%
    group_by(status, cluster)%>%
    summarize(mean_cyto = mean(cyto),
              se = stats::sd(cyto)/sqrt(n()))
  
 cyto_plot_19 <-  ggplot(data_cyto_plot, aes(x = cluster, y = mean_cyto, color = status, fill = status))+
   geom_bar(stat = 'identity', position = 'dodge')+
    geom_errorbar(aes(x = cluster, y = mean_cyto, ymin = mean_cyto-se, ymax =mean_cyto+se), position = position_dodge(0.9), width = 0.7, color ='black', size =1)+
   labs(x ='Subcluster', y = 'Mean CytoTRACE +/- SE')+
   geom_signif(xmin = 1, xmax = 1.33,
               y_position = 0.8,
               annotation = '*' , 
               color = "black",
               tip_length = c(0.02,.15), textsize = 8,
               size = 1)+
    geom_signif(xmin = 1.66, xmax = 2,
                y_position = 0.55,
                annotation = '**' , 
                color = "black",
                tip_length = c(0.2,0.02), textsize = 8,
                size = 1)+
    geom_signif(xmin = 1.66, xmax = 2.33,
                y_position = 0.73,
                annotation = '*' , 
                color = "black",
                tip_length = c(0,.25), textsize = 8,
                size = 1)+
    geom_signif(xmin = 3, xmax = 3.33,
                y_position = 0.55,
                annotation = '*' , 
                color = "black",
                tip_length = c(0.02,0.20), textsize = 8,
                size = 1)+
    geom_signif(xmin = 0.66, xmax = 2.33,
                y_position = 1,
                annotation = '***' , 
                color = "black",
                tip_length = c(0,0), textsize = 8,
                size = 1)+
    geom_signif(xmin = 0.66, xmax = 3.33,
                y_position = 1.2,
                annotation = '***' , 
                color = "black",
                tip_length = c(0,0), textsize = 8,
                size = 1)+
   theme_minimal()+
   geom_text(label = "M", aes(x = 0.66, y = -0.03), color= 'black', size =3)+
   geom_text(label = "D", aes(x = 1, y = -0.03), color= 'black', size =3)+
   geom_text(label = "F", aes(x = 1.33, y = -0.03), color= 'black', size =3)+
   geom_text(label = "M", aes(x = 1.66, y = -0.03), color= 'black', size =3)+
   geom_text(label = "D", aes(x = 2, y = -0.03), color= 'black', size =3)+
   geom_text(label = "F", aes(x = 2.33, y = -0.03), color= 'black', size =3)+
   geom_text(label = "M", aes(x = 2.66, y = -0.03), color= 'black', size =3)+
   geom_text(label = "D", aes(x = 3, y = -0.03), color= 'black', size =3)+
   geom_text(label = "F", aes(x = 3.33, y = -0.03), color= 'black', size =3)+
   ylim(-0.04, 1.35)+
   theme(legend.position = 'none', axis.title = element_text(size = 14), axis.text =element_text(size = 12), axis.text.x =element_text(angle =-90))
 cyto_plot_19
 
 ggsave(plot = cyto_plot_19,
        file = "cyto_plot_19.svg",
        device = "svg",
        units = "in",
        width = 2,
        height = 3,
        path = "Bachelors Thesis/Plots/Figure 4")
 
 ### proportion differences ####
 
sub_cells = cluster_19@meta.data%>%
   group_by(sub, Status, individual)%>%
   subset(Status %in% c('M','D','F'))%>%
   summarize(cells = n())
 
 total_cells = cluster_19@meta.data%>%
   group_by(Status, individual)%>%
   subset(Status %in% c('M','D','F'))%>%
   summarize(total_cells = n())
 
 full_data <- sub_cells%>%
   right_join(total_cells, by = 'individual')
 
 glmer_matrix <- matrix(NA, nrow(full_data),2)
 glmer_matrix[,1] <- full_data$cells
 glmer_matrix[,2] <- full_data$total_cells -full_data$cells
 
 glmer_model <- glmer(glmer_matrix~sub*Status.x +(1|individual),data = full_data, family = binomial('logit'))
 car::Anova(glmer_model, type = 'III')
 
 
 ### DEGS #####
 neg_bin_mult_windows<- readRDS('Functions/DEG_functions/neg_bin_mult_windows.rds')
 
 neg_bin <- data.frame()
 for (i in 0:2) {
   cluster <- paste0('19_',i)
   print(cluster)
   output <- neg_bin_mult_windows(obj = cluster_19,
                          clustering = 'sub',
                          cluster = cluster)
   output$cluster = cluster
   neg_bin <- rbind(neg_bin, output)
 }
 
 neg_bin$f_m_q.value <- p.adjust(neg_bin$f_m_p.value, 'fdr', nrow(neg_bin))
 neg_bin$d_m_q.value <- p.adjust(neg_bin$d_m_p.value, 'fdr', nrow(neg_bin))
 neg_bin$d_f_q.value <- p.adjust(neg_bin$d_f_p.value, 'fdr', nrow(neg_bin))
 
 neg_bin_defined <- define_degs(neg_bin)
 
 neg_bin_defined_filtered <- neg_bin_defined%>%
   filter(!is.na(issignif)& is.na(warning))
 
 #write.csv(neg_bin_defined_filtered, 'DEG Outputs/subclusters_19_03_18_2025.csv')
 
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
        file = "19_degs.svg",
        device = "svg",
        units = "in",
        width = 3,
        height = 2,
        path = "Bachelors Thesis/Plots/Figure 4")
 
 

 
 neg_bin_defined_filtered[neg_bin_defined_filtered$cluster=='19_0' ,]
 neg_bin_defined_filtered[neg_bin_defined_filtered$cluster=='19_1' ,]
 neg_bin_defined_filtered[neg_bin_defined_filtered$cluster=='19_2' ,]
 
 clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_0'])%>%dotplot()
 clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_1'])%>%dotplot()
 clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_2'])%>%dotplot()
 
 clown_go(neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_0' & neg_bin_defined_filtered$class == 'Late Upregulated'])%>%dotplot()
 
 
 ### plotting hsd17b14 ####
 clust_19_1_data <- cluster_19@assays$RNA$data[,cluster_19$sub =='19_1']%>%
   t()%>%
   as.data.frame()%>%
   dplyr::select('hsd17b14')
 
 hsd17b14_data <- data.frame(hsd17b14_expression = clust_19_1_data,
                              individual = cluster_19$individual[obj$sub =='19_1'],
                              Sex = cluster_19$Status[cluster_19$sub =='19_1'])%>%
   group_by(individual, Sex)%>%
   summarize(mean_expression = mean(hsd17b14),
             se = sd(hsd17b14)/sqrt(n()))%>%
   subset(Sex != 'NRM')%>%
   na.omit()
 library(ggsignif)
 
 hsd17b14_data$Sex <- factor(hsd17b14_data$Sex, levels = c('M','D','E',"NF",'F'))
 hsd17b14_data_plot <- ggplot(hsd17b14_data, aes(x = Sex, y = mean_expression))+
   geom_boxplot(alpha = 0.25, outlier.shape = NA, aes(color = Sex, fill = Sex))+
   geom_jitter(  shape = 1, color = 'black', size =2)+
   theme_classic()+
   theme(legend.position ='none',plot.title = element_text(hjust=0.5))+
   labs(title = '19_1', y = 'Mean hsd17b14')+
   geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(.17), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)+
   ylim(-0.02,.2)
 hsd17b14_data_plot
 
 ggsave(plot = hsd17b14_data_plot,
        file = "hsd17b14_data_plot.svg",
        device = "svg",
        units = "in",
        width = 1.6,
        height = 2,
        path = "Bachelors Thesis/Plots/Figure 4")
 
 ### plot pgr
 clust_19_1_data_pgr <- cluster_19@assays$RNA$data[,cluster_19$sub =='19_2']%>%
   t()%>%
   as.data.frame()%>%
   dplyr::select('pgr')
 
 pgr_data <- data.frame(pgr_expression = clust_19_1_data_pgr,
                             individual = cluster_19$individual[obj$sub =='19_2'],
                             Sex = cluster_19$Status[cluster_19$sub =='19_2'])%>%
   group_by(individual, Sex)%>%
   summarize(mean_expression = mean(pgr),
             se = sd(pgr)/sqrt(n()))%>%
   subset(Sex != 'NRM')%>%
   na.omit()
 library(ggsignif)
 
 pgr_data$Sex <- factor(pgr_data$Sex, levels = c('M','D','E',"NF",'F'))
 pgr_data_plot <- ggplot(pgr_data, aes(x = Sex, y = mean_expression))+
   geom_boxplot(alpha = 0.25, outlier.shape = NA, aes(color = Sex, fill = Sex))+
   geom_jitter(  shape = 1, color = 'black', size =2)+
   theme_classic()+
   theme(legend.position ='none',plot.title = element_text(hjust=0.5))+
   labs(title = '19_2', y = 'Mean pgr')+
   geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(1.8), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=6)+
   ylim(0,2)
 pgr_data_plot
 
 ggsave(plot = pgr_data_plot,
        file = "pgr_data_plott.svg",
        device = "svg",
        units = "in",
        width = 1.6,
        height = 2,
        path = "Bachelors Thesis/Plots/Figure 4")
 
 
 ##subcluster markers
 Idents(cluster_19) <- 'sub'
 sub_markers <- FindAllMarkers(cluster_19)
 
markers_19_0 <- sub_markers$gene[sub_markers$cluster == '19_0' & sub_markers$pct.1>sub_markers$pct.2 & sub_markers$p_val_adj <0.05]
clown_go(markers_19_0)%>%dotplot()

go_19_0 <- clown_go(markers_19_0)
go_19_0$geneID[go_19_0$Description =='nervous system development']

markers_19_1 <- sub_markers$gene[sub_markers$cluster == '19_1' & sub_markers$pct.1>sub_markers$pct.2 & sub_markers$p_val_adj <0.05]
clown_go(markers_19_1)%>%dotplot()

markers_19_2 <- sub_markers$gene[sub_markers$cluster == '19_2' & sub_markers$pct.1>sub_markers$pct.2 & sub_markers$p_val_adj <0.05]
clown_go(markers_19_2)%>%dotplot()


### differences in proportion ###
total_cells_individual <- cluster_19@meta.data%>%
  group_by(individual)%>%
  summarize(total_cells = n())

total_cells_cluster <- cluster_19@meta.data%>%
  group_by(sub,individual, Status)%>%
  summarize(cluster_cells = n())

full_data <- total_cells_cluster%>%
  right_join(total_cells_individual, by  = 'individual')%>%
  subset(Status != 'NRM')

full_data$success = full_data$cluster_cells
full_data$failure = full_data$total_cells - full_data$cluster_cells

### all 3 have prop differences, so I'm going to add in pairwise comparisons
out <- data.frame()
for(i in unique(cluster_19$sub)){
  data <- full_data%>%
    subset(sub == i & Status %in% c('M','D','F'))
  
logistic <- glmer(cbind(data$success, data$failure)~Status +(1|individual), 
                  data = data, 
                  family = binomial('logit'))
type_iii <- car::Anova(logistic, type = 'III')

pairs <- pairs(emmeans(logistic, "Status", type = 'log'), adjust = 'none')%>%
  as.data.frame() 

newd <- data.frame(cluster = i,
                   anova_p = type_iii$`Pr(>Chisq)`[1],
                   singular = isSingular(logistic),
                   d_f_p.value = pairs$p.value[pairs$contrast == 'D - F'],
                   d_m_p.value = pairs$p.value[pairs$contrast == 'D - M'],
                   f_m_p.value = pairs$p.value[pairs$contrast == 'F - M']
)
out <- rbind(out, newd)
}
full_data$Status <- factor(full_data$Status, levels = c('M','D',"E","NF",'F'))
prop_19_0 <- ggplot(subset(full_data, sub=='19_0'), aes(y = cluster_cells/total_cells, x = Status, fill = Status))+
  geom_boxplot(alpha = 0.25, outlier.shape = NA,aes(color = Status))+
  geom_jitter(  shape = 1, color = 'black', size =2)+
  theme_classic()+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size =12),
        axis.text = element_text(size = 10))+
  labs(x = 'Sex', y = 'Proportion of Cells', title = '19_0')+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(.65), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=6)+
  ylim(0.03,0.75)
prop_19_0

ggsave(plot = prop_19_0,
       file = "prop_19_0.svg",
       device = "svg",
       units = "in",
       width = 1.6,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 4")

prop_19_1 <- ggplot(subset(full_data, sub=='19_1'), aes(y = cluster_cells/total_cells, x = Status, fill = Status))+
  geom_boxplot(alpha = 0.25, outlier.shape = NA,aes(color = Status))+
  geom_jitter(  shape = 1, color = 'black', size =2)+
  theme_classic()+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size =12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10))+
  labs(x = 'Sex', y = 'Proportion of Cells', title = '19_1')+
  ylim(0.03,0.75)
prop_19_1

ggsave(plot = prop_19_1,
       file = "prop_19_1.svg",
       device = "svg",
       units = "in",
       width = 1.2,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 4")

windowsFonts(Arial=windowsFont("Arial"))
windowsFonts(Calibri=windowsFont("Calibri"))

prop_19_2 <- ggplot(subset(full_data, sub=='19_2'), aes(y = cluster_cells/total_cells, x = Status, fill = Status))+
  geom_boxplot(alpha = 0.25, outlier.shape = NA,aes(color = Status))+
  geom_jitter(  shape = 1, color = 'black', size =2)+
  theme_classic()+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size =12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10))+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(.65), annotation =c(""), color = "black", tip_length = c(0,0), textsize=4)+
  labs(x = 'Sex', y = 'Proportion of Cells', title = '19_2')+
  ylim(0.03,0.75)
prop_19_2

ggsave(plot = prop_19_2,
       file = "prop_19_2.svg",
       device = "svg",
       units = "in",
       width = 1.2,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 4")


DimPlot(cluster_19, reduction = 'atacUMAP')


#### slingshot trajectory analysis ####

sce= SingleCellExperiment(assays = list(
  counts = GetAssayData(cluster_19, slot='counts'),
  logcounts = GetAssayData(cluster_19, slot='data')
))

reducedDims(sce)= list(
  PCA = Embeddings(cluster_19),
  UMAP = Embeddings(cluster_19, 'harmony_wnn.umap')
)

colData(sce)$cell_type = cluster_19$sub
colData(sce)$cluster = cluster_19$sub


set.seed(123)
library(slingshot)

sce <- slingshot(sce, 
                 clusterLabels = 'cluster',
                 reducedDim = 'UMAP',
                 start.clus = '19_0')


pseudotime_slingshot <- slingPseudotime(sce)
avg_pseudotime <- rowMeans(pseudotime_slingshot, na.rm = T)
avg_pseudotime[is.na(avg_pseudotime)] <- min(avg_pseudotime, na.rm =T)
avg_pseudotime_norm <- 100* (avg_pseudotime - min(avg_pseudotime))/
  (max(avg_pseudotime)- min(avg_pseudotime))

cluster_19$slingshot_pseudotime <- avg_pseudotime_norm

curves = slingCurves(sce)
umap_coords <- reducedDims(sce)$UMAP

plot_data <- data.frame(UMAP1 =umap_coords[,1],
                        UMAP2 = umap_coords[,2],
                        Pseudotime = cluster_19$slingshot_pseudotime,
                        CellType = cluster_19$harmony.wnn_res0.4_clusters)

slingshot_plot <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Pseudotime))+
  geom_point(size =0.8, alpha = 0.7)+
  scale_color_viridis_c()+
  labs(title = 'Slingshot Trajectories')+
  theme_minimal()

#add curves
for(i in seq_along(curves)){
  curve_data <- data.frame(UMAP1 = curves[[i]]$s[,1],
                           UMAP2 = curves[[i]]$s[,2])
  slingshot_plot <- slingshot_plot+ 
    geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
              size = 1, color = 'black', alpha = 0.8)
  
}

slingshot_plot

std_umap <-DimPlot(cluster_19, reduction = 'harmony_wnn.umap', label=T)

print(slingshot_plot + std_umap)
#alright well that's i guess a result

#### Diffusion pseudotime #####
library(destiny)
sce <- as.SingleCellExperiment(cluster_19)
var_genes = VariableFeatures(cluster_19)
dm = DiffusionMap(sce, k =30, n_pcs =20, verbose =T)
dc <- eigenvectors(dm)[,1:2]
cluster_19[['diffmap']]<- CreateDimReducObject(
  embeddings = dc,
  key = 'DC_',
  assay = 'RNA'
)
root_cell_type <- '19_0'
#find radial glia indicies
root_cell_indicies <- which(cluster_19$sub == root_cell_type)

#choose median radial glia as root
root_cell_dc <- dc[root_cell_indicies,]
root_cell_center <- colMeans(root_cell_dc)

#calculate distances for all cells
dists <- apply(root_cell_dc, 1, function(x) sum((x- root_cell_center)^2))

#set smallest distance as root cell
root_cell <- root_cell_indicies[which.min(dists)]
dpt <- DPT(dm, tips = root_cell)
cluster_19$diffusion_pseudotime <- rank(dpt$dpt) #why do rank also what is rank 
cluster_19$diffusion_pseudotime <- 100* (cluster_19$diffusion_pseudotime-
                                    min(cluster_19$diffusion_pseudotime))/
  (max(cluster_19$diffusion_pseudotime)-
     min(cluster_19$diffusion_pseudotime)) ##have no idea what this does, I think it normalizes?
# ok so its 100 * (score above min)/range(score) so its a percentage I think
p1 <- ggplot(data.frame(DC1 = dc[,1],
                        DC2 = dc[,2],
                        CellType = cluster_19$sub))+
  geom_point(aes(x = DC1, y = DC2, color = CellType), size = 1)+
  theme_minimal()+
  guides(color = guide_legend(override.aes =list(size =3)))+ #waht does this do??
  labs(title= "Diffusion Map", x = 'DC1',y = 'Dc2')

p2 <- FeaturePlot(cluster_19, features = 'diffusion_pseudotime',
                  reduction = 'harmony_wnn.umap',
                  label = T,
                  pt.size = 1)+
  scale_color_viridis_c()+
  labs(title = 'Diffusion Pseudotime')

p1+p2

## theyre all going to predict that 190-191-192, but that'ts just kind of cause that's how the graph is setup
# i want to try rna velocity

diff_score_data <- data.frame(individual = cluster_19$individual,
                              Status = cluster_19$Status, 
                              cluster = cluster_19$sub,
                              pseudotime = cluster_19$diffusion_pseudotime)%>%
  group_by(individual, Status, cluster)%>%
  summarize(pseudotime = mean(pseudotime),
            se_pseudotime = sd(pseudotime)/sqrt(n()))
library(forcats)
diff_score_data$Status <- factor(diff_score_data$Status, levels = c('NRM','M','D',"E",'NF','F'))

plot <- ggplot(diff_score_data, aes(x = fct_reorder(cluster, pseudotime), y = pseudotime))+
  geom_jitter(size = 1, alpha =0.5, aes(color = Status))+
  geom_boxplot(fill = NA, outlier.shape = NA)+
  theme_minimal()+
  labs(x = 'Cluster',y = 'Mean +/- SE Diff Pseudotime')#+
 # stat_summary(geom = 'line', fun = 'mean', aes(group = Status, color= Status)) # conceptually this doesnt make sense
plot

phase_data <- data.frame(Phase =cluster_19@meta.data$Phase,
                         sub = cluster_19$sub,
                         Status = cluster_19$Status,
                         individual = cluster_19$individual)%>%
  group_by(individual, Status, sub, Phase)%>%
  summarize(n_phase = n())%>%
  right_join(total_cells_cluster, by = c('sub', 'Status', 'individual'))%>%
  mutate(percent_n_phase = n_phase/cluster_cells)%>%
  subset(Status != 'NRM')

phase_data$Status <- factor(phase_data$Status, levels = c('M','D','E','NF','F'))

 ggplot(phase_data, aes( x = Status, y = percent_n_phase, color = Status, fill = Status))+
   geom_boxplot(alpha =0.25,
                outlier.shape = NA)+
   geom_jitter(size = 2, color ='black', shape = 1)+
   facet_grid(Phase~sub)
 

 ggplot(subset(phase_data, sub == '19_2'), aes( x = Status, y = percent_n_phase, color = Status, fill = Status))+
   geom_boxplot(alpha =0.25,
                outlier.shape = NA)+
   geom_jitter(size = 2, color ='black', shape = 1)+
   facet_grid(~Phase)
 
 #ok when I have more time I will do a mixed logistic I think there is somethhing here
 #except it makes no sense, neurons shouldnt be in S phase
 
 #####################################################







