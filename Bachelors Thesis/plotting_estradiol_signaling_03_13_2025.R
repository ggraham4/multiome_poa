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
  full_data$Status <- factor(full_data$Status, levels = c('M','D','E','NF',"F"))
  plot <- ggplot(full_data, aes(x = Status, y = mean, color = Status))+
    geom_boxplot(outlier.shape = NA, alpha = 0.25, aes(fill= Status))+
    labs(y = paste0(interacting_pair), title = paste0('Clusters ',cluster_pair))
  plot
  return(plot)
  
  
}

sig_data <- read_csv("A:/CellPhoneDB 030225/signif_data_bound.csv")

estradiol_2_19 <- plot_communication_cluster_pair('Estradiol_byCYP19A1_ESR2', '2|19',individuals_meta_data )
library(ggsignif)
es2 <- estradiol_2_19+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'Estradiol Signaling', title = '2>19')+
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5))+
  ylim(0.45, 1.8)+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(1.6), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)
es2

ggsave(plot = es2,
       file = "estradiol_2_19.svg",
       device = "svg",
       units = "in",
       width = 1.2,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 3")


estradiol_2_2 <- plot_communication_cluster_pair('Estradiol_byCYP19A1_ESR2', '2|2',individuals_meta_data )
es2_2 <- estradiol_2_2+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'Estradiol Signaling', title = '2>2')+
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
  ylim(0.45, 1.8)+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(1.6), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(1.6), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=6)
es2_2

ggsave(plot = es2_2,
       file = "estradiol_2_2.svg",
       device = "svg",
       units = "in",
       width = 1.6,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 3")

estradiol_2_11 <- plot_communication_cluster_pair('Estradiol_byCYP19A1_ESR2', '2|11',individuals_meta_data )
es2_11 <- estradiol_2_11+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'Estradiol Signaling', title = '2>11')+
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5))+
  ylim(0.45, 1.8)+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(1.6), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(1.6), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=6)
es2_11

ggsave(plot = es2_11,
       file = "estradiol_2_11.svg",
       device = "svg",
       units = "in",
       width = 1.2,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 3")


estradiol_2_8 <- plot_communication_cluster_pair('Estradiol_byCYP19A1_ESR2', '2|8',individuals_meta_data )
es2_8 <- estradiol_2_8+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'Estradiol Signaling', title = '2>8')+
  theme(legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5))+
  ylim(0.45, 1.8)+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(1.6), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(1.6), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)
es2_8

ggsave(plot = es2_8,
       file = "estradiol_2_8.svg",
       device = "svg",
       units = "in",
       width = 1.2,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 3")





  
