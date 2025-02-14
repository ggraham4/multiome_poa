library(ggplot2)
library(ggsignif)
library(tidyverse)

data <- read_csv("Measures/all_data.csv")

data$Status_dummy = data$Status
data$Status_dummy <- factor(data$Status_dummy, levels = c('M',
                                                          'F',
                                                          'D',
                                                          'S',
                                                          'E',
                                                          'EP',
                                                          'NM',
                                                          'NF'
))

data$Status <- ifelse(data$Status %in% c('NM','NF'), NA, data$Status)

##Behaviors day 2 ####
data$Behaviors_Day_2 <- as.numeric(data$Behaviors_Day_2)


beh_model <- lm(Behaviors_Day_2~Status_dummy, data = data)
beh_pairs <- as.data.frame(pairs(emmeans(beh_model, 'Status_dummy'), adjust = 'none'))
beh_pairs$issignif <- ifelse(beh_pairs$p.value<0.05, '*', NA)


behavior_plot <- ggplot(data, aes(x = Status_dummy, y = Behaviors_Day_2))+
  geom_boxplot(data = subset(data, !is.na(Status)),aes(group = Status, fill = Status, , color = Status), alpha = 0.25, outlier.shape = NA)+
  geom_point(size = 2, color = 'black', shape = 1, position = position_dodge2(0.5))+
  theme_classic()+
  geom_signif(xmin = c(1), xmax = c(4), y_position = c(max(data$Behaviors_Day_2)*1.2), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(3), xmax = c(4), y_position = c(max(data$Behaviors_Day_2)*1.05), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  labs(x  ='Sex', y = 'Parental Care Behaviors')+
  ylim(c(0, max(data$Behaviors_Day_2)*1.25))+
  theme(axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0, size = 10), legend.position = 'none')
behavior_plot

ggsave(plot = behavior_plot,
      file = "behavior_plot.tiff",
      device = "tiff",
      units = "in",
      width = 2.5,
      height = 2.5,
      path = "Bachelors Thesis/Plots/Figure 1")

##Time day 2 ###
data$Time_Day_2 <- as.numeric(data$Time_Day_2)

time_model <- lm(Time_Day_2~Status_dummy, data = data)
time_pairs <- as.data.frame(pairs(emmeans(time_model, 'Status_dummy'), adjust = 'none'))
time_pairs$issignif <- ifelse(time_pairs$p.value<0.05, '*', NA)

time_plot <- ggplot(data, aes(x = Status_dummy, y = Time_Day_2))+
  geom_boxplot(data = subset(data, !is.na(Status)),aes(group = Status, fill = Status, , color = Status), alpha = 0.25, outlier.shape = NA)+
  geom_point(size = 2, color = 'black', shape = 1, position = position_dodge2(0.5))+
  theme_classic()+
  labs(x  ='Sex', y = 'Time in Nest (s)')+
  geom_signif(xmin = c(1), xmax = c(2), y_position = c(max(data$Time_Day_2)*1.05), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(1), xmax = c(4), y_position = c(max(data$Time_Day_2)*1.2), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(3), xmax = c(4), y_position = c(max(data$Time_Day_2)*1.05), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=5)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0), legend.position = 'none')+
  ylim(c(0, max(data$Time_Day_2)*1.25))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
  
time_plot

ggsave(plot = time_plot,
      file = "time_plot.tiff",
      device = "tiff",
      units = "in",
      width = 2.5,
      height = 2.5,
      path = "Bachelors Thesis/Plots/Figure 1")

###Testosterone ####
data$Log_11KT <- as.numeric(data$Log_11KT)

kt_model <- lm(Log_11KT~Status_dummy, data = data)
kt_pairs <- as.data.frame(pairs(emmeans(kt_model, 'Status_dummy'), adjust = 'none'))
kt_pairs$issignif <- ifelse(kt_pairs$p.value<0.05, '*', NA)

kt_data <- subset(data, !is.na(Log_11KT))
kt_plot <- ggplot(kt_data, aes(x = Status_dummy, y = Log_11KT))+
  geom_boxplot(data = subset(data, !is.na(Status)),aes(group = Status, fill = Status, , color = Status), alpha = 0.25, outlier.shape = NA)+
  geom_point(size = 2, color = 'black', shape = 1, position = position_dodge2(0.5))+
  theme_classic()+
  labs(x  ='Sex', y = 'Log10 11KT (pg/ml)')+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(max(kt_data$Log_11KT)*1.05), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(1), xmax = c(4), y_position = c(max(kt_data$Log_11KT)*1.5), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(2), xmax = c(4), y_position = c(max(kt_data$Log_11KT)*1.3), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
   geom_signif(xmin = c(2.1), xmax = c(3), y_position = c(max(kt_data$Log_11KT)*1.05), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
   theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0), legend.position = 'none')+
  ylim(c(0, max(kt_data$Log_11KT)*1.6))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
  
kt_plot

ggsave(plot = kt_plot,
      file = "kt_plot.tiff",
      device = "tiff",
      units = "in",
      width = 2.5,
      height = 2.5,
      path = "Bachelors Thesis/Plots/Figure 1")

##percent testicuar ##
data$Percent_Testicular <- as.numeric(data$Percent_Testicular)
test_model <- lm(Percent_Testicular~Status_dummy, data = data)
test_pairs <- as.data.frame(pairs(emmeans(test_model, 'Status_dummy'), adjust = 'none'))
test_pairs$issignif <- ifelse(test_pairs$p.value<0.05, '*', NA)


`%notin%` <- Negate(`%in%`)

test_data <- data[data$Status %notin% c('F','NF'),]
test_data <- test_data[test_data$Status_dummy %notin% c('F','NF'),]
test_data$Status_dummy <- factor(test_data$Status_dummy, 
                                 levels = c('M',
                                            'D',
                                            'S',
                                            'E',
                                            'EP',
                                            'NM'))
test_data$Status <- factor(test_data$Status, 
                                 levels = c('M',
                                            'D',
                                            'S',
                                            'E',
                                            'EP',
                                            'NM'))

test_data<- test_data[!is.na(test_data$Percent_Testicular),]

test_plot <- ggplot(test_data, aes(x = Status_dummy, y = Percent_Testicular))+
  geom_boxplot(data = subset(test_data, !is.na(Status)),aes(group = Status, fill = Status, , color = Status), alpha = 0.25, outlier.shape = NA)+
  geom_point(size = 2, color = 'black', shape = 1, position = position_dodge2(0.5))+
  theme_classic()+
  labs(x  ='Sex', y = 'Percent Testicular')+
  geom_signif(xmin = c(1), xmax = c(2), y_position = c(max(test_data$Percent_Testicular)*1.05), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(1), xmax = c(3), y_position = c(max(test_data$Percent_Testicular)*1.2), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0), legend.position = 'none')+
  ylim(c(0, max(test_data$Percent_Testicular)*1.25))+
  scale_color_manual(values =c('#619CFF', '#F8766D','#DB72FB','#D39200','#00BA38'))+
  scale_fill_manual(values =c('#619CFF', '#F8766D','#DB72FB','#D39200','#00BA38'))

  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
  
test_plot

ggsave(plot = test_plot,
      file = "test_plot.tiff",
      device = "tiff",
      units = "in",
      width = 2.5,
      height = 2.5,
      path = "Bachelors Thesis/Plots/Figure 1")

##Percent Ovarian ###

data$Percent_Ovarian <- as.numeric(data$Percent_Ovarian)
ov_model <- lm(Percent_Ovarian~Status_dummy, data = data)
ov_pairs <- as.data.frame(pairs(emmeans(ov_model, 'Status_dummy'), adjust = 'none'))
ov_pairs$issignif <- ifelse(ov_pairs$p.value<0.05, '*', NA)

ov_data <- data[data$Status %notin% c('F','NF'),]
ov_data <- ov_data[ov_data$Status_dummy %notin% c('F','NF'),]
ov_data$Status_dummy <- factor(ov_data$Status_dummy, 
                                 levels = c('M',
                                            'D',
                                            'S',
                                            'E',
                                            'EP',
                                            'NM'))
ov_data$Status <- factor(ov_data$Status, 
                           levels = c('M',
                                      'D',
                                      'S',
                                      'E',
                                      'EP',
                                      'NM'))

ov_data<- ov_data[!is.na(ov_data$Percent_Ovarian),]

ov_plot <- ggplot(ov_data, aes(x = Status_dummy, y = Percent_Ovarian))+
  geom_boxplot(data = subset(ov_data, !is.na(Status)),aes(group = Status, fill = Status, , color = Status), alpha = 0.25, outlier.shape = NA)+
  geom_point(size = 2, color = 'black', shape = 1, position = position_dodge2(0.5))+
  theme_classic()+
  labs(x  ='Sex', y = 'Percent Ovarian')+
  geom_signif(xmin = c(1), xmax = c(2), y_position = c(max(ov_data$Percent_Ovarian)*1.05), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(1), xmax = c(3), y_position = c(max(ov_data$Percent_Ovarian)*1.2), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0), legend.position = 'none')+
  ylim(c(0, max(ov_data$Percent_Ovarian)*1.25))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+

  scale_color_manual(values =c('#619CFF', '#F8766D','#DB72FB','#D39200','#00BA38'))+
  scale_fill_manual(values =c('#619CFF', '#F8766D','#DB72FB','#D39200','#00BA38'))
ov_plot

ggsave(plot = ov_plot,
      file = "ov_plot.tiff",
      device = "tiff",
      units = "in",
      width = 2.5,
      height = 2.5,
      path = "Bachelors Thesis/Plots/Figure 1")

## volume estimate ####

data$Log10_Volume <- as.numeric(data$Log10_Volume)
ev_2.5x_model <- lm(Log10_Volume~Status_dummy, data = data)
ev_2.5x_pairs <- as.data.frame(pairs(emmeans(ev_2.5x_model, 'Status_dummy'), adjust = 'none'))
ev_2.5x_pairs$issignif <- ifelse(ev_2.5x_pairs$p.value<0.05, '*', NA)

ev_2.5x_data <- data[data$Status %notin% c('F','NF'),]
ev_2.5x_data <- ev_2.5x_data[ev_2.5x_data$Status_dummy %notin% c('F','NF'),]
ev_2.5x_data$Status_dummy <- factor(ev_2.5x_data$Status_dummy, 
                               levels = c('M',
                                          'D',
                                          'S',
                                          'E',
                                          'EP',
                                          'NM'))
ev_2.5x_data$Status <- factor(ev_2.5x_data$Status, 
                         levels = c('M',
                                    'D',
                                    'S',
                                    'E',
                                    'EP',
                                    'NM'))

ev_2.5x_plot <- ggplot(ev_2.5x_data, aes(x = Status_dummy, y = Log10_Volume))+
  geom_boxplot(data = subset(ev_2.5x_data, !is.na(Status)),aes(group = Status, fill = Status, , color = Status), alpha = 0.25, outlier.shape = NA)+
  geom_point(size = 2, color = 'black', shape = 1, position = position_dodge2(0.5))+
  theme_classic()+
  labs(x  ='Sex', y = 'Log10 Gonadal Volume (px)')+
  geom_signif(xmin = c(2), xmax = c(3), y_position = c(max(ev_2.5x_data$Log10_Volume)*1.03), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(1), xmax = c(3), y_position = c(max(ev_2.5x_data$Log10_Volume)*1.1), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0), legend.position = 'none')+
  ylim(c(5.5, max(ev_2.5x_data$Log10_Volume)*1.13))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+
  scale_color_manual(values =c('#619CFF', '#F8766D','#DB72FB','#D39200','#00BA38'))+
  scale_fill_manual(values =c('#619CFF', '#F8766D','#DB72FB','#D39200','#00BA38'))
ev_2.5x_plot

ggsave(plot = ev_2.5x_plot,
      file = "log10 volume.tiff",
      device = "tiff",
      units = "in",
      width = 2.5,
      height = 2.5,
      path = "Bachelors Thesis/Plots/Figure 1")

###UMAP #### --- BEGIN FIG 2 #####
obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

obj_subset = obj[,obj@meta.data$harmony.wnn_res0.4_clusters!=15 & obj@meta.data$harmony.wnn_res0.4_clusters!=30]

umap <- DimPlot(obj_subset, 
                reduction = 'harmony_wnn.umap',
                group.by = 'harmony.wnn_res0.4_clusters', label = T, label.size = 3, pt.size = 0.0001)+
  theme(legend.position = 'none')+
  labs(x = 'UMAP_1', y = 'UMAP_2', title = NULL)+
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
  
umap

ggsave(plot = umap,
      file = "umap.tiff",
      device = "tiff",
      units = "in",
      width = 2.5,
      height = 2.5,
      path = "Bachelors Thesis/Plots/Figure 2")

##RNA ####
rna_umap <- DimPlot(obj_subset, 
                reduction = 'harmony_rnaUMAP',
                group.by = 'harmony.wnn_res0.4_clusters', label = T, label.size = 3, pt.size = 0.0001)+
  theme(legend.position = 'none')+
  labs(x = 'UMAP_1', y = 'UMAP_2', title = NULL)+
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
  
rna_umap

ggsave(plot = rna_umap,
      file = "rna_umap.tiff",
      device = "tiff",
      units = "in",
      width = 2.5,
      height = 2.5,
      path = "Bachelors Thesis/Plots/Figure 2")

###ATAC ####
atac_umap <- DimPlot(obj_subset, 
                reduction = 'atacUMAP',
                group.by = 'harmony.wnn_res0.4_clusters', label = T, label.size = 3, pt.size = 0.0001)+
  theme(legend.position = 'none')+
  labs(x = 'UMAP_1', y = 'UMAP_2', title = NULL)+
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))

atac_umap

ggsave(plot = atac_umap,
      file = "atac_umap.tiff",
      device = "tiff",
      units = "in",
      width = 2.5,
      height = 2.5,
      path = "Bachelors Thesis/Plots/Figure 2")


Idents(obj_subset) <- "harmony.wnn_res0.4_clusters"
markers <- unique(c(
  'elavl3',#neuron
  'gad2',#gABA 
  'LOC111588076', #gad1
  'LOC111584103', #vgat2.1
  'slc17a6b', #vglut
  'slc17a7a', #vglut1
  'sst1.1', #interneuron marker
  'LOC111577263', #brain aromatase - radial glia
  'gfap', #astrocyte marker
  'crocc2', #ependymal cell marker
  'mbpa', #oligo marker
  'cspg4', #OPC marker
  'p2ry12', #microglia marker
  'ptprc', #leukocyte marker
  'oxt',
  'avp',
  'slc18a3b', #ach marker
  'kiss1',
  'kiss1ra',
  'kiss1rb',
  'tac1',
  'tacr1a',
  'tac3a',
  'tacr3a',
  'tacr3l',
  'npy',
  'esr1',
  'esr2a',
  'esr2b',
  'ar',
  'LOC111562384', #ccka
  'cckar',
  'cckb',
  'cckbra',
  'cckbrb',
  'gal',
  'galr1a',
  'rbm47',
  'gdpd5a',
  'col15a1b',
  'flvcr2b'
) )


marker_plot  <- DotPlot(object = obj_subset, 
                 group.by = "harmony.wnn_res0.4_clusters", 
                 features = markers,
                 dot.min=0.1
) + 
  coord_flip()+
  theme(axis.text.x = element_text(angle = -45),
        dot.scale =3)+
  scale_x_discrete(labels= c(
    'elavl3',
    'gad2',
    'gad1',
    'vglut2.1',
    'slc17a6b (GLUT)',
    'slc17a7a (GLUT)',
    'sst1.1',
    'cyp19a1b',
    'gfap',
    'crocc2',
    'mbpa',
    'cspg4',
    'p2ry12',
    'ptprc',
    'oxt',
    'avp',
    'slc18a3b (ach)',
    'kiss1',
    'kiss1ra',
    'kiss1rb',
    'tac1',
    'tacr1a',
    'tac3a',
    'tacr3a',
    'tacr3L',
    'npy',
    'esr1',
    'esr2a',
    'esr2b',
    'ar',
    'ccka',
    'cckar',
    'cckb',
    'cckbra',
    'cckbrb',
    'gal',
    'galr1a',
    'rbm47',
    'gdpd5a',
    'col15a1b',
    'flvcr2b'                 
    ))+
  labs(y = 'Cluster')+
  theme(axis.title.y =  element_blank())+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
  scale_radius(range = c(0,3))
marker_plot

ggsave(plot = marker_plot,
       file = "marker_plot.tiff",
       device = "tiff",
       units = "in",
       width = 7.5,
       height = 5,
       path = "Bachelors Thesis/Plots/Figure 2")

