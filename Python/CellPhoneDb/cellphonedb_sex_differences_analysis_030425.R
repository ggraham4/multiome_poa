### cellphone db differential communication analysis

library(lme4)
library(tidyverse)
library(Seurat)

human_named = readRDS("C:/Users/Gabe/Desktop/RNA_object_human_names.rds")

## first, figure out all pathways is analyzed
example_data = read_tsv("A:/CellPhoneDB 030225/degs_analysis_means_03_03_2025_224407.txt")

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
test_pull = readRDS('A:/CellPhoneDB 030225/coalesced_list_030525.RDS')


hist(cluster_pair_list[["Dihydrotestosterone_bySRD5A1_AR"]]$mean)
#ok so these might actually be negative binomially distributed
hist(cluster_pair_list[["CCK_CCKBR"]]$mean) #this is a neuropeptide so it should have high expression
hist(cluster_pair_list[["CCK_CCKAR"]]$mean) #ok these seem to be gamma distributed
hist(cluster_pair_list[["TAC1_TACR1"]]$mean) 
hist(cluster_pair_list[["TAC1_TACR3"]]$mean) 
#i wonder if model like glmer.nb(interacting_pair~Status*cluster_pair) might be better as a broad model and then 
#each subsequent model only evaluates what is significnat in the first model

## add sex into the df
individual_by_sex_data = data.frame(human_named$individual, human_named$Status)%>%distinct()
ccka_cckar_sex = cluster_pair_list[["CCK_CCKAR"]]%>%
  right_join(individual_by_sex_data, by =join_by('individual'== 'human_named.individual'))%>%
  na.omit()


###Ok so the data is actually not normally distributed, it is gamma distributed
library(glmmTMB)

test_model = glmmTMB(mean+1~human_named.Status*cluster_pair,
                     data = ccka_cckar_sex,
                     family = Gamma())

summary(test_model)

# a way to cut down this data might be to pull in "relevant_interactions.txt" and eliminate any interactions that are 0

