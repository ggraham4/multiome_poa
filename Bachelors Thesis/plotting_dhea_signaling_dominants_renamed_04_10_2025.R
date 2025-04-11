### plotting differential communication 
library(lme4)
library(tidyverse)
library(Seurat)
library(car)
library(emmeans)

human_named = readRDS("C:/Users/Gabe/Desktop/RNA_object_human_names.rds")

unique_cell_cell_pairs = colnames(example_data[,14:ncol(example_data)])

individuals = unique(human_named$individual)
individuals = individuals[individuals!= 'GH']


##lets do an example first 

t13d_data = read_tsv("A:/CellPhoneDB 030225/T13D/degs_analysis_means_T13D.txt")

t13_mean_cyp19_esr2_2_29 = as.list(t13d_data[t13d_data$interacting_pair == 'Estradiol_byCYP19A1_ESR2', '2|19'])%>%unname()%>%unlist()


means_together = data.frame()
for(i in individuals){
  message(i)
  data = read_tsv(paste0("A:/CellPhoneDB 030225/",i,"/degs_analysis_means_",i,".txt"))%>%suppressMessages()
  mean = as.list(data[data$interacting_pair == 'Estradiol_byCYP19A1_ESR2', '2|19'])%>%unname()%>%unlist()
  
  if(length(data)==0){next}
  newd = data.frame(individual= i, 
                    mean = mean)
  means_together <- rbind(means_together, newd)
  
}
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

dhea_24_20 <- plot_communication_cluster_pair('Dehydroepiandrosterone_bySTS_PPARA', '24|20',individuals_meta_data )
library(ggsignif)

dh24_20 <- dhea_24_20+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'DHEA Signaling',, title = '24>20')+
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
  ylim(0, 0.8)+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(0.7), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(0.7), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)
dh24_20

ggsave(plot = dh24_20,
       file = "dhea_24_20.svg",
       device = "svg",
       units = "in",
       width = 1.6,
       height = 2,
       path = "Bachelors Thesis/Plots_dominants_renamed/Figure 3")

dhea_24_21 <- plot_communication_cluster_pair('Dehydroepiandrosterone_bySTS_PPARA', '24|21',individuals_meta_data )
dh24_21 <- dhea_24_21+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'DHEA Signaling', title = '24>21')+
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5))+
  ylim(0, 0.8)+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(0.7), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(0.7), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)
dh24_21

ggsave(plot = dh24_21,
       file = "dhea_24_21.svg",
       device = "svg",
       units = "in",
       width = 1.2,
       height = 2,
       path = "Bachelors Thesis/Plots_dominants_renamed/Figure 3")







