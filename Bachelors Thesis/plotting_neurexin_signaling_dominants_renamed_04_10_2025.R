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
  full_data$Status <- ifelse(full_data$Status=='D','I',full_data$Status)
  full_data$Status <- ifelse(full_data$Status=='E','LI',full_data$Status)
  
  full_data$Status <- factor(full_data$Status, levels = c('M','I','LI','NF',"F"))
  plot <- ggplot(full_data, aes(x = Status, y = mean, color = Status))+
    geom_boxplot(outlier.shape = NA, alpha = 0.25, aes(fill= Status))+
    labs(y = paste0(interacting_pair), title = paste0('Clusters ',cluster_pair))
  plot
  return(plot)
  
  
}

sig_data <- read_csv("A:/CellPhoneDB 030225/signif_data_bound.csv")

nrxn_19_4 <- plot_communication_cluster_pair('LRRTM4_NRXN2', '19|4',individuals_meta_data )
library(ggsignif)

nr_19_4 <- nrxn_19_4+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'Neurexin Signaling',, title = '19>4')+
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
  ylim(1, 4)+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(3), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(3), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(5), y_position = c(3.6), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)

nr_19_4

ggsave(plot = nr_19_4,
       file = "nr_19_4.svg",
       device = "svg",
       units = "in",
       width = 1.6,
       height = 2,
       path = "Bachelors Thesis/Plots_dominants_renamed/Figure 3")

nrxn_10_24 <- plot_communication_cluster_pair('LRRTM4_NRXN2', '10|24',individuals_meta_data )
nr_10_24 <- nrxn_10_24+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'DHEA Signaling', title = '10>24')+
  theme(legend.position = 'none', axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5))+
  ylim(2, 5)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(3.9), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(5), y_position = c(4.5), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)

nr_10_24

ggsave(plot = nr_10_24,
       file = "nr_10_24.svg",
       device = "svg",
       units = "in",
       width = 1.35,
       height = 2,
       path = "Bachelors Thesis/Plots_dominants_renamed/Figure 3")







