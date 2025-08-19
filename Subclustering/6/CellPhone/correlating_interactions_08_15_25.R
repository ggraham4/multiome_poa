### Analyzing CellPhoneDB Data with 6 subclusters, whole and sub ###
### GJG 08_15_2025 ###
## libaries ####
{
  library(lme4)
  library(tidyverse)
  library(Seurat)
  library(car)
  library(emmeans)
  `%notin%` <- Negate(`%in%`)
  
}

##read in data #####
#seurat
human_named = readRDS("A:/optimal_clustering_05_06_2025/RNA_object_human_names.rds")

#example lists
whole_example_data = read_tsv("A:/cellphone_6_08_14_2025/whole_output/degs_analysis_means_08_14_2025_124528.txt")
sub_example_data = read_tsv("A:/cellphone_6_08_14_2025/sub_output/degs_analysis_means_08_14_2025_133918.txt")

#significance data 
whole_significant = read.csv('A:/cellphone_6_08_14_2025/whole_output/signif_data_bound.csv')
sub_significant = read.csv('A:/cellphone_6_08_14_2025/sub_output/signif_data_bound.csv')

#individuals
individuals = unique(human_named$individual[!human_named$individual %in% 'GH'])

##whole
#individual data
whole_base_path = "A:/cellphone_6_08_14_2025/whole_output/"

whole_individual_data = list()
for(i in individuals){
  message(paste0("Individual  ", which(individuals ==i), ' of ', length(individuals)))
  file_name = paste0(whole_base_path,i,"/degs_analysis_means_",i,".txt")
  whole_means = read_tsv(file_name)
  whole_individual_data[[i]]= whole_means
}

##sub
#individual data
sub_base_path = "A:/cellphone_6_08_14_2025/sub_output/"

sub_individual_data = list()
for(i in individuals){
  message(paste0("Individual  ", which(individuals ==i), ' of ', length(individuals)))
  file_name = paste0(sub_base_path,i,"/degs_analysis_means_",i,".txt")
  sub_means = read_tsv(file_name)
  sub_individual_data[[i]]= sub_means
}

#info about individuals
measures =read.csv("Measures/all_data.csv")

### define functions ####
coalesce_interaction = function(LRpair, CCpair, whole_sub){
  "
  the goal here is to coalesce the means of each individual for the given
  cluster cluster pair and ligand receptor pair
  LRpair : Ligand - receptor pair, string e.g., 'TAC1_TACR3'
  CCpair : Cluster - cluster pair, sring e.g., '0|10'
  whole_sub: which list to reference, should be of level 'whole' or 'sub'
  "
  
  if(whole_sub == 'whole'){
    ref_signif = whole_significant
    ref_means = whole_individual_data
  }
  if(whole_sub == 'sub'){
    ref_signif = sub_significant
    ref_means = sub_individual_data
    
  }
  
  coalesced_data = data.frame()
  for(ind in individuals){
    individual_df =  ref_means[[ind]]
    CCpair_column = which(colnames(individual_df)==CCpair)
    LRpair_row = which(individual_df$interacting_pair== LRpair)
    
    individual_mean = (individual_df[LRpair_row, CCpair_column]) %>% as.numeric()
    if(length(individual_mean) == 0 || is.na(individual_mean)){
      next  # Skip this individual if mean is NA or empty
    }    
    newd = data.frame(individual = ind,
                      mean = individual_mean)
    coalesced_data =rbind(coalesced_data, newd)
  }
  full_data = coalesced_data%>%
    right_join(measures, by = join_by('individual'=='Fish'))%>%
    subset(Status %in% c('M','D','E','NF','F'))
  full_data$Status = factor(full_data$Status, levels = c('M','D','E','NF','F'))
  return(full_data)
}

plot_interaction = function(LRpair, CCpair, whole_sub){
  if(whole_sub=='whole'){
    status_q_value = whole_significant$main_effect_q.value[whole_significant$cluster_pair == CCpair & 
                                                             whole_significant$interacting_pair == LRpair]
    status_q_value = round(status_q_value, 5)
    
  } else if(whole_sub=='sub'){
    status_q_value = sub_significant$main_effect_q.value[sub_significant$cluster_pair == CCpair & 
                                                           sub_significant$interacting_pair == LRpair]
    status_q_value = round(status_q_value, 5)
  } else {
    message('Plot ERROR: whole_sub not "whole" or "sub"')
    status_q_value = 'Invalid whole_sub parameter'
  }
  
  if(length(status_q_value) == 0 || is.na(status_q_value)){
    status_q_value = ' Interaction not relevant'
  }  
  
  data_for_plot = coalesce_interaction(LRpair, CCpair, whole_sub)
  
  data_for_plot = data_for_plot[!is.na(data_for_plot$mean),]
  plot = ggplot(data_for_plot, aes(x = Status, y = mean, 
                                   color = Status, 
                                   fill = Status, 
                                   shape = Status))+
    geom_point(position = position_jitterdodge(2))+
    geom_boxplot(alpha = 0.25, outlier.shape = NA)+
    theme_minimal()+
    labs(x = 'Status', y = 'Interaction Mean',
         title = paste0(LRpair, ' ', CCpair),
         subtitle = paste0('Status q.value =',status_q_value ))
  
  return(plot)
  
}
n_interactions <- length(unique(sub_significant$interacting_pair[sub_significant$main_effect_q.value<0.05]))
n_clusters <- length(unique(sub_significant$cluster_pair))

# Adjusted calculation - roughly half since we're avoiding duplicates
interactions <- ((n_interactions^2 - n_interactions) / 2) * n_clusters

pairwise_corr = function(){
  i = 0
  data_new = data.frame()
  
  interactions_list <- unique(sub_significant$interacting_pair[sub_significant$main_effect_q.value<0.05])
  
  for(idx1 in 1:(length(interactions_list)-1)){
    interaction_1 <- interactions_list[idx1]
    
    for(idx2 in (idx1+1):length(interactions_list)){
      interaction_2 <- interactions_list[idx2]
      
      # Skip if first 4 letters are the same OR last 4 letters are the same
      if(substr(interaction_1, 1, 4) == substr(interaction_2, 1, 4) ||
         substr(interaction_1, nchar(interaction_1)-3, nchar(interaction_1)) == 
         substr(interaction_2, nchar(interaction_2)-3, nchar(interaction_2))){
        next
      }
      
      for(cluster in unique(sub_significant$cluster_pair)){
        i = i + 1
        message(paste0(i, ' out of ', interactions))
        
        data_to_corr_1 = coalesce_interaction(interaction_1, cluster, 'sub')
        data_to_corr_2 = coalesce_interaction(interaction_2, cluster, 'sub')
        
        joined <- data_to_corr_1 %>%
          right_join(data_to_corr_2, by = 'individual')
        
        test = cor.test(joined$mean.x, joined$mean.y)
        newd = data.frame(interaction_1 = interaction_1,
                          interaction_2 = interaction_2,
                          cluster_pair = cluster,
                          p = test$p.value,
                          correlation = test$estimate)
        data_new = rbind(data_new, newd)
      }
    }
  }
  return(data_new)
}

corrs = pairwise_corr()

# Now do FDR adjustment
corrs$p_adjusted <- p.adjust(corrs$p, method = "fdr")
corrs_filtered= subset(corrs, p_adjusted<0.05)
