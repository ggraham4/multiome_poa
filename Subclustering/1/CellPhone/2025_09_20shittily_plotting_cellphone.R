human_named = readRDS("A:/optimal_clustering_05_06_2025/RNA_object_human_names.rds")


whole_base_path = "Subclustering/1/CellPhone/Whole/"
whole_significant = read.csv('Subclustering/1/CellPhone/Whole/signif_data_bound.csv')

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


### lets look at some interactions  ####
whole_significant[whole_significant$main_effect_q.value<0.05,]

### estradiol
whole_significant[whole_significant$main_effect_q.value<0.05 &
                    whole_significant$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA',]
#individuals
individuals = unique(human_named$individual[!human_named$individual %in% 'GH'])

##whole
#individual data

whole_individual_data = list()
for(i in individuals){
  message(paste0("Individual  ", which(individuals ==i), ' of ', length(individuals)))
  file_name = paste0(whole_base_path,i,"/degs_analysis_means_",i,".txt")
  whole_means = read_tsv(file_name)
  whole_individual_data[[i]]= whole_means
}

measures =read.csv("Measures/all_data.csv")

plot_interaction('Dehydroepiandrosterone_bySTS_PPARA', '1_0|0', 'whole') 

plot_interaction('RTN4R_ADGRB1', '3|1_3', 'whole') 

plot_interaction('KISS1_KISS1R', '10|1_3', 'whole') 

plot_interaction('NRG3_ERBB4', '1_5|3', 'whole') 

plot_interaction('Glutamate_byGLS_and_SLC17A7_GRM3', '22|1_2', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC17A7_GRM3', '24|1_2', 'whole') 
plot_interaction('SEMA3A_NRP1', '1_1|14', 'whole') 
plot_interaction('SEMA3A_PlexinA2_complex1', '1_1|15', 'whole') 
plot_interaction('WNT5B_FZD3_LRP5', '14|1_1', 'whole') 
plot_interaction('NRG3_ERBB4', '1_5|1_1', 'whole') 

plot_interaction('RTN4R_ADGRB1', '5|1_3', 'whole') 
plot_interaction('RTN4R_ADGRB1', '16|1_3', 'whole') 
plot_interaction('GABA_byGAD1_and_SLC6A11_GABA-A_a4b3d_complex', '1_3|10', 'whole') 


signif_data_bound[signif_data_bound$cluster_pair=='10|6' &signif_data_bound$main_effect_p.value<0.05, ]

plot_interaction('Glutamate_byGLS_and_SLC1A6_GRM8', '10|6', 'whole') 
# increase in glutamate to potentially gnrh+ cells does not support my hypothesis

#it seems like the cells are at least touching?
plot_interaction('GABA_byGAD1_and_SLC6A1_GABA-A_a5b3g2S_complex', '10|6', 'whole') 
# but also increase in gaba?
plot_interaction('GABA_byGAD1_and_SLC6A6_GABA-A_a5b3g2S_complex', '10|6', 'whole') 


signif_data_bound[signif_data_bound$cluster_pair=='6|10' &signif_data_bound$main_effect_p.value<0.05, ]
#they def seem to be touchiing
plot_interaction('Glutamate_byGLS_and_SLC17A6_GRIA3', '6|10', 'whole') 

# confused why are they sending both to eachother. 10 is like fully glutamatergic I think
