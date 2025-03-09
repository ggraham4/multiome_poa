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
  plot <- ggplot(full_data, aes(x = Status, y = mean, color = Status, shape = Status))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(size = 3)+
    scale_shape_manual(values = c(1,2,3,4,5))+
    theme_classic()+
    labs(y = paste0(interacting_pair), title = paste0('Clusters ',cluster_pair))
  plot
  return(plot)

  
}

plot_communication_cluster_pair('Estradiol_byCYP19A1_ESR2', '2|19',individuals_meta_data )
  
plot_communication_cluster_pair('Glutamate_byGLS_and_SLC1A6_GRM1', '12|19',individuals_meta_data )

plot_communication_cluster_pair('Glutamate_byGLS_and_SLC17A7_GRM7', '27|14',individuals_meta_data )

plot_communication_cluster_pair('Dehydroepiandrosterone_bySTS_PPARA', '24|21',individuals_meta_data )

  
  
  

