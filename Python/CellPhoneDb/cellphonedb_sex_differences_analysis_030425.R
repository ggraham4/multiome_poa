### cellphone db differential communication analysis

library(lme4)
library(tidyverse)
library(Seurat)

human_named = readRDS("C:/Users/Gabe/Desktop/RNA_object_human_names.rds")

## first, figure out all pathways is analyzed
example_data = read_tsv("A:/CellPhoneDB 030225/simple_analysis_means_result_03_03_2025_194650.txt")

unique_interactions = unique(example_data$interacting_pair)
unique_cell_cell_pairs = colnames(example_data[,14:ncol(example_data)])

base_path = 'A:/CellPhoneDB 030225/'
#read in all sample data
individual_data = list()
for(i in unique(human_named$individual)){
  if(i =='GH'){next}
  print(i)
  
  mean_data = read_tsv(paste0(base_path, i,'/degs_analysis_means_',i,'.txt'),show_col_types = FALSE)
  individual_data[[i]]$means = mean_data
  
}

# test the distribution of an example first
head(example_data)

dist_data = data.frame()
for(i in unique(human_named$individual)){
  if(i =='GH'){next}
  data =  individual_data[[i]]$means
  
  mean = data[data$interacting_pair =='CDH2_CDH2', "0|0"]
  
  newd = data.frame(individual =i, 
                         mean = mean)
  dist_data = rbind(dist_data, newd)
}

hist(dist_data$X0.0)
#normal

### test a few more examples
individuals = unique(human_named$individual)
individuals = individuals[individuals!= 'GH']

`%notin%` <- Negate(`%in%`)
for(interacting_pair in sample(unique_interactions,10)){
  interacting_pair_data = data.frame()
  for(i in individuals){
    data =  individual_data[[i]]$means
    
    columns = colnames(data)
    if("19|27" %notin% columns){next}
    
    mean = data[data$interacting_pair == interacting_pair, "19|27"]
    
    if(nrow(mean)==0){next}
    
    newd = data.frame(individual =i, 
                      mean = mean)
    interacting_pair_data = rbind(interacting_pair_data, newd)
  }
  hist(interacting_pair_data[,2])
}
### if more than 4 == 0, use fisher, else, use lm

#first, make a loop to coalesce all the data

cluster_pair_list = list()  
for(interacting_pair in unique_interactions){
  print(interacting_pair)
  
  full_cluster_pair_data = data.frame()
 
  for(cluster_pair in unique_cell_cell_pairs){
    
    cluster_pair_data = data.frame()
    for(i in individuals){
      
      data =  individual_data[[i]]$means
      
      columns = colnames(data)
      
      if(cluster_pair %notin% columns){next}
      
      mean = data[data$interacting_pair == interacting_pair, cluster_pair]
      colnames(mean) = 'mean'
      
      if(nrow(mean)==0){next}
      
      newd = data.frame(individual =i,
                        interacting_pair = interacting_pair,
                        cluster_pair = cluster_pair,
                       mean = mean)
      cluster_pair_data = rbind(cluster_pair_data, newd)
      
      
    }
    full_cluster_pair_data =rbind(cluster_pair_data, full_cluster_pair_data)
  }
  cluster_pair_list[[interacting_pair]] =full_cluster_pair_data
  

}

saveRDS(cluster_pair_list, 'A:/CellPhoneDB 030225/coalesced_list_030525.RDS')




