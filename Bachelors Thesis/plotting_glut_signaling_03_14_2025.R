### plotting differential communication 
library(lme4)
library(tidyverse)
library(Seurat)
library(car)
library(emmeans)

human_named = readRDS("C:/Users/Gabe/Desktop/RNA_object_human_names.rds")


individuals = unique(human_named$individual)
individuals = individuals[individuals!= 'GH']

individuals_meta_data <- data.frame(individual = human_named$individual,
                                    Status = human_named$Status)%>%
  distinct()


plot_communication_cluster_pair = function(interacting_pair, cluster_pair, individuals_meta_data){
  `%notin%` <- Negate(`%in%`)
  individuals = individuals_meta_data$individual
  individuals = individuals[individuals!= 'GH']
  
  means_together = data.frame()
  for(i in individuals){
    message(i)
    data = read_tsv(paste0("A:/CellPhoneDB 030225/",i,"/degs_analysis_means_",i,".txt"))%>%suppressMessages()
    if(cluster_pair %notin% colnames(data)){next}
    mean = as.list(data[data$interacting_pair == interacting_pair, cluster_pair])
    mean= unname(unlist(mean))
    if(length(mean)==0){next}
    newd = data.frame(individual= i, 
                      mean = mean)
    means_together <- rbind(means_together, newd)
  }
  
  #add sex 
  full_data <- means_together%>%
    right_join(individuals_meta_data, by = 'individual')%>%
    na.omit()
  
  ## plot
  full_data$Status <- factor(full_data$Status, levels = c('M','D','E','NF',"F"))
  plot <- ggplot(full_data, aes(x = Status, y = mean, color = Status))+
    geom_boxplot(outlier.shape = NA, alpha = 0.25, aes(fill= Status))+
    labs(y = paste0(interacting_pair), title = paste0('Clusters ',cluster_pair))
  plot
  return(plot)
  
  
}

sig_data <- read_csv("A:/CellPhoneDB 030225/signif_data_bound.csv")

glut_12_19 <- plot_communication_cluster_pair('Glutamate_byGLS_and_SLC1A6_GRM1', '12|19',individuals_meta_data )
library(ggsignif)

g_12_19 <- glut_12_19+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'Glutamate Signaling',, title = '12>19')+
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
  ylim(0, 1.5)+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(1.2), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(1.2), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)

g_12_19

ggsave(plot = g_12_19,
       file = "g_12_19.svg",
       device = "svg",
       units = "in",
       width = 1.6,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 3")

glut_27_14 <- plot_communication_cluster_pair('Glutamate_byGLS_and_SLC17A7_GRM7', '27|14',individuals_meta_data )
g_27_14 <- glut_27_14+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'DHEA Signaling', title = '27>14')+
  theme(legend.position = 'none', axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  ylim(0, 1.5)+
  geom_signif(xmin = c(2), xmax = c(5), y_position = c(0.9), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(5), y_position = c(1.2), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)

g_27_14

ggsave(plot = g_27_14,
       file = "g_27_14.svg",
       device = "svg",
       units = "in",
       width = 1.35,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 3")






