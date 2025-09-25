### cellphone db differential communication analysis

library(lme4)
library(tidyverse)
library(Seurat)
`%notin%` <- Negate(`%in%`)

human_named = readRDS("A:/optimal_clustering_05_06_2025/RNA_object_human_names.rds")

## first, figure out all pathways is analyzed
example_data = read_tsv("Subclustering/1/CellPhone/Whole/degs_analysis_significant_means_09_19_2025_155229.txt")

unique_interactions = unique(example_data$interacting_pair)
unique_cell_cell_pairs = colnames(example_data[,14:ncol(example_data)])

base_path = "Subclustering/1/CellPhone/Whole/"
#read in all sample data
individual_data = list()
for(i in unique(human_named$individual)){
  if(i =='GH'){next}
  print(i)
  
  mean_data = read_tsv(paste0(base_path, i,'/degs_analysis_means_',i,'.txt'),show_col_types = FALSE)
  individual_data[[i]]$means = mean_data
  
}

#first, make a loop to coalesce all the data
individuals = unique(human_named$individual)

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
### breakpoint ###

saveRDS(cluster_pair_list, paste0(base_path,'coalesced_list_09_19_2025.RDS'))

test_pull = readRDS(paste0(base_path,'coalesced_list_09_19_2025.RDS'))
### ####
## read in reference data
example_data = read_tsv("A:/cellphone_1_09_19_2025/whole_output/degs_analysis_means_09_19_2025_155229.txt")
coalesced_data = readRDS('A:/cellphone_1_09_19_2025/whole_output/coalesced_list_09_19_2025.RDS')
base_path = 'A:/cellphone_1_09_19_2025/whole_output/'

unique_interactions = unique(example_data$interacting_pair)
unique_cell_cell_pairs = colnames(example_data[,14:ncol(example_data)])

individuals = unique(human_named$individual)
individuals = individuals[individuals!= 'GH']

## read in significance data
significant_interactions_reference = read_tsv("A:/cellphone_1_09_19_2025/whole_output/degs_analysis_relevant_interactions_09_19_2025_155229.txt")

significant_interactions_reference_filtered = significant_interactions_reference[,-c(1,3:13)]

significant_interactions_reference_filtered_pivoted = significant_interactions_reference_filtered%>%
  pivot_longer(cols = colnames(significant_interactions_reference_filtered[,-1]), values_to = 'signif', names_to = 'cluster_pair')

significant_interactions_reference_filtered_pivoted_signif_only = subset(significant_interactions_reference_filtered_pivoted, signif ==1)

nrow(significant_interactions_reference_filtered_pivoted_signif_only)/nrow(significant_interactions_reference_filtered_pivoted)
## only 19% are significant this is great

#Find the unique significant ligand receptor interactions
significant_interactions = unique(significant_interactions_reference_filtered_pivoted_signif_only$interacting_pair)

#subset the BIG list to that
coalesced_data_subset = coalesced_data[which(names(coalesced_data) %in% significant_interactions)]

## also add sex to the list
status_data = data.frame(individual = human_named$individual,
                         Status = human_named$Status)%>%
  distinct()


#for loop for subsetting out the cluster interactions that are significant
new_list = list()
for(interacting_pair in unique(names(coalesced_data_subset))){
  print(interacting_pair)
  query_data= coalesced_data_subset[[interacting_pair]]
  
  reference_data = significant_interactions_reference_filtered_pivoted_signif_only[significant_interactions_reference_filtered_pivoted_signif_only$interacting_pair == interacting_pair,]
  
  query_data_signif = query_data%>%
    right_join(reference_data, by = 'cluster_pair')%>%
    subset(signif == 1)%>%
    right_join(status_data, by = 'individual')
  
  new_list[[interacting_pair]] = query_data_signif
  
}

saveRDS(new_list, 'A:/cellphone_1_09_19_2025/whole_output/signif_only_list_09_19_2025.RDS')

### ####

unique_cell_cell_pairs = colnames(example_data[,14:ncol(example_data)])

individuals = unique(human_named$individual)
individuals = individuals[individuals!= 'GH']

signif_list = readRDS( 'A:/cellphone_1_09_19_2025/whole_output/signif_only_list_09_19_2025.RDS')


unique_interactions = unique(names(signif_list))

library(emmeans)

signif_data = list()
for(interacting_pair in unique_interactions){
  message(paste0('Interaction ',which(interacting_pair ==unique_interactions), ' of ', length(unique_interactions)))
  
  full_data = signif_list[[interacting_pair]]
  full_data =full_data[full_data$Status %in% c('M','D','F'),]%>%
    na.omit()
  
  interacting_pair_data = data.frame()
  for(cluster_pair in unique(full_data$cluster_pair)){
    
    temp_data = full_data[full_data$cluster_pair == cluster_pair,]
    
if(sum(temp_data$mean==0)>=4){
  # run a fisher's exact test, first, binarize the data
  temp_data$binary = ifelse(temp_data$mean == 0, 0, 1)
  
  test_type = 'Fisher'
  
  # Create contingency table for each comparison separately
  
  ##main effect - all three groups
  fisher_data_main = temp_data %>%
    group_by(Status) %>%
    summarize(success = sum(binary),
              failure = n() - success, .groups = 'drop')
  
  fisher_matrix_main <- fisher_data_main[,-1] %>%
    as.matrix()
  
  whole_fish = fisher.test(fisher_matrix_main)
  main_effect_p.value = whole_fish$p.value
  
  ### d_m comparison
  fisher_data_d_m = temp_data %>%
    filter(Status %in% c('D', 'M')) %>%
    group_by(Status) %>%
    summarize(success = sum(binary),
              failure = n() - success, .groups = 'drop')
  
  d_m_matrix = fisher_data_d_m[,-1] %>%
    as.matrix()
  
  d_m_fish = fisher.test(d_m_matrix)
  d_m_p.value = d_m_fish$p.value
  
  ###d_f comparison
  fisher_data_d_f = temp_data %>%
    filter(Status %in% c('D', 'F')) %>%
    group_by(Status) %>%
    summarize(success = sum(binary),
              failure = n() - success, .groups = 'drop')
  
  d_f_matrix = fisher_data_d_f[,-1] %>%
    as.matrix()
  
  d_f_fish = fisher.test(d_f_matrix)
  d_f_p.value = d_f_fish$p.value
  
  ###f_m comparison
  fisher_data_f_m = temp_data %>%
    filter(Status %in% c('F', 'M')) %>%
    group_by(Status) %>%
    summarize(success = sum(binary),
              failure = n() - success, .groups = 'drop')
  
  f_m_matrix = fisher_data_f_m[,-1] %>%
    as.matrix()
  
  f_m_fish = fisher.test(f_m_matrix)
  f_m_p.value = f_m_fish$p.value
  
      new_data = data.frame(
        interacting_pair = interacting_pair, 
        cluster_pair = cluster_pair,
        test_type = test_type,
        main_effect_p.value = main_effect_p.value,
        d_m_p.value = d_m_p.value,
        d_f_p.value = d_f_p.value,
        f_m_p.value = f_m_p.value
      )
      interacting_pair_data = rbind(interacting_pair_data, new_data)
    }
    else{
      
      test_type = 'Linear Regression'
      
      model = lm(mean ~ Status, data = temp_data)
      av = car::Anova(model, type = 'III')%>%as.data.frame()
      
      main_effect_p.value = av$`Pr(>F)`[2]
      
      pairs = pairs(emmeans(model, 'Status'), adjust = 'none')%>%as.data.frame()
      
      d_m_p.value = pairs$p.value[pairs$contrast == 'D - M']
      f_m_p.value = pairs$p.value[pairs$contrast == 'F - M']
      d_f_p.value = pairs$p.value[pairs$contrast == 'D - F']
      
      new_data = data.frame(
        interacting_pair = interacting_pair, 
        cluster_pair = cluster_pair,
        test_type = test_type,
        main_effect_p.value = main_effect_p.value,
        d_m_p.value = d_m_p.value,
        d_f_p.value = d_f_p.value,
        f_m_p.value = f_m_p.value
      )
      interacting_pair_data = rbind(interacting_pair_data, new_data)
      
    }
    
    
  }
  signif_data[[interacting_pair]] =interacting_pair_data
}

signif_data2 =signif_data

for(i in names(signif_data2)){
  signif_data2[[i]]$main_effect_q.value = p.adjust(signif_data2[[i]]$main_effect_p.value, 'fdr',nrow(signif_data2[[i]]))
  
  
}

# do p adjust within interacting pair not across, 72K will kill any p values

signif_data_bound = do.call(rbind, signif_data2)
rownames(signif_data_bound) =(1:nrow(signif_data_bound))

write_csv(signif_data_bound, 'A:/cellphone_1_09_19_2025/whole_output/signif_data_bound.csv')


for(i in names(signif_data2)){
  if(all(signif_data2[[i]]$test_type=='Fisher')){
    signif_data2[[i]]$main_effect_q.value_no_fisher = 1
  }
  else{
    signif_data2[[i]]$main_effect_q.value_no_fisher = p.adjust(signif_data2[[i]]$main_effect_p.value, 'fdr',nrow(signif_data2[[i]]))
  }
  
}

signif_data_no_fisher = signif_data2
for(i in names(signif_data_no_fisher)){
  
  signif_data_no_fisher[[i]] <- subset(signif_data_no_fisher[[i]], test_type != 'Fisher')
  
  signif_data_no_fisher[[i]]$main_effect_q.value = p.adjust(signif_data_no_fisher[[i]]$main_effect_p.value, 'fdr',nrow(signif_data_no_fisher[[i]]))
  
  
}

signif_data_no_fisher_bound = do.call(rbind, signif_data_no_fisher)
rownames(signif_data_no_fisher_bound) =(1:nrow(signif_data_no_fisher_bound))

write_csv(signif_data_no_fisher_bound, 'A:/cellphone_1_09_19_2025/whole_output/signif_data_bound_no_fisher.csv')



