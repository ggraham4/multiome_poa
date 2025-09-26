### cellphone db differential communication analysis - OPTIMIZED

library(lme4)
library(tidyverse)
library(Seurat)
`%notin%` <- Negate(`%in%`)

human_named = readRDS("A:/optimal_clustering_05_06_2025/RNA_object_human_names.rds")

## first, figure out all pathways is analyzed
example_data = read_tsv("A:/cellphone_1_09_19_2025/whole_output/degs_analysis_means_09_19_2025_155229.txt")

unique_interactions = unique(example_data$interacting_pair)
unique_cell_cell_pairs = colnames(example_data[,14:ncol(example_data)])

base_path = "A:/cellphone_1_09_19_2025/whole_output/"
#read in all sample data
individual_data = list()
for(i in unique(human_named$individual)){
  if(i =='GH'){next}
  print(i)
  
  mean_data = read_tsv(paste0(base_path, i,'/degs_analysis_means_',i,'.txt'),show_col_types = FALSE)
  individual_data[[i]]$means = mean_data
  
}

#first, make a loop to coalesce all the data - OPTIMIZED
individuals = unique(human_named$individual)

cluster_pair_list = list()  
for(interacting_pair in unique_interactions){
  print(interacting_pair)
  
  # Collect all results for this interacting pair in one go
  all_results = vector("list", length = length(unique_cell_cell_pairs) * length(individuals))
  result_idx = 1
  
  for(cluster_pair in unique_cell_cell_pairs){
    for(i in individuals){
      
      data = individual_data[[i]]$means
      columns = colnames(data)
      
      if(cluster_pair %notin% columns){next}
      
      mean = data[data$interacting_pair == interacting_pair, cluster_pair]
      
      if(length(mean)==0 || nrow(mean)==0){next}
      
      all_results[[result_idx]] = data.frame(
        individual = i,
        interacting_pair = interacting_pair,
        cluster_pair = cluster_pair,
        mean = mean[[1]]  # Extract the actual value
      )
      result_idx = result_idx + 1
    }
  }
  
  # Combine all results for this interaction
  all_results = all_results[!sapply(all_results, is.null)]
  if(length(all_results) > 0){
    cluster_pair_list[[interacting_pair]] = do.call(rbind, all_results)
    rownames(cluster_pair_list[[interacting_pair]]) = NULL
  }
}
### breakpoint ###

saveRDS(cluster_pair_list, 'A:/cellphone_1_09_19_2025/whole_output/coalesced_list_09_19_2025.RDS')
test_pull = readRDS('A:/cellphone_1_09_19_2025/whole_output/coalesced_list_09_19_2025.RDS')

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

#for loop for subsetting out the cluster interactions that are significant - OPTIMIZED
new_list = list()
for(interacting_pair in unique(names(coalesced_data_subset))){
  print(interacting_pair)
  query_data = coalesced_data_subset[[interacting_pair]]
  
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

# OPTIMIZED STATISTICAL ANALYSIS
signif_data = list()

for(interacting_pair in unique_interactions){
  message(paste0('Interaction ',which(interacting_pair ==unique_interactions), ' of ', length(unique_interactions)))
  
  full_data = signif_list[[interacting_pair]]
  full_data = full_data[full_data$Status %in% c('M','D','F'),]%>%
    na.omit()
  
  # Get unique cluster pairs once
  unique_cluster_pairs = unique(full_data$cluster_pair)
  
  # Pre-allocate list for results
  results_list = vector("list", length = length(unique_cluster_pairs))
  names(results_list) = unique_cluster_pairs
  
  for(idx in seq_along(unique_cluster_pairs)){
    cluster_pair = unique_cluster_pairs[idx]
    temp_data = full_data[full_data$cluster_pair == cluster_pair,]
    
    if(sum(temp_data$mean==0)>=4){
      #run a fisher's exact test, first , binarize the data
      temp_data$binary = ifelse(temp_data$mean == 0, 0, 1)
      
      test_type = 'Fisher'
      
      # ok so I need to collect main_effect_p.value and pairwise comparisons, I expect I will have no power
      
      ##main effect
      fisher_data = temp_data%>%
        group_by(Status)%>%
        summarize(success = sum(binary),
                  failure = n() - success, .groups = 'drop')
      
      fisher_matrix <- fisher_data[,-1]%>%
        as.matrix()
      
      whole_fish = fisher.test(fisher_matrix)
      main_effect_p.value = whole_fish$p.value
      
      ### d_m
      d_m_matrix = fisher_data%>%
        subset(Status != 'F')%>%
        dplyr::select(-1)%>%
        as.matrix()
      
      d_m_fish = fisher.test(d_m_matrix)
      d_m_p.value = d_m_fish$p.value
      
      ###d_f
      d_f_matrix = fisher_data%>%
        subset(Status != 'M')%>%
        dplyr::select(-1)%>%
        as.matrix()
      
      d_f_fish = fisher.test(d_f_matrix)
      d_f_p.value = d_f_fish$p.value
      
      ###f_m
      f_m_matrix = fisher_data%>%
        subset(Status != 'D')%>%
        dplyr::select(-1)%>%
        as.matrix()
      
      f_m_fish = fisher.test(f_m_matrix)
      f_m_p.value = f_m_fish$p.value
      
      results_list[[idx]] = data.frame(
        interacting_pair = interacting_pair, 
        cluster_pair = cluster_pair,
        test_type = test_type,
        main_effect_p.value = main_effect_p.value,
        d_m_p.value = d_m_p.value,
        d_f_p.value = d_f_p.value,
        f_m_p.value = f_m_p.value
      )
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
      
      results_list[[idx]] = data.frame(
        interacting_pair = interacting_pair, 
        cluster_pair = cluster_pair,
        test_type = test_type,
        main_effect_p.value = main_effect_p.value,
        d_m_p.value = d_m_p.value,
        d_f_p.value = d_f_p.value,
        f_m_p.value = f_m_p.value
      )
    }
  }
  
  # Combine results for this interacting pair
  results_list = results_list[!sapply(results_list, is.null)]
  if(length(results_list) > 0){
    signif_data[[interacting_pair]] = do.call(rbind, results_list)
  }
}

signif_data2 = signif_data

# Vectorized p-adjustment
for(i in names(signif_data2)){
  signif_data2[[i]]$main_effect_q.value = p.adjust(signif_data2[[i]]$main_effect_p.value, 'fdr', nrow(signif_data2[[i]]))
}

# do p adjust within interacting pair not across, 72K will kill any p values

signif_data_bound = do.call(rbind, signif_data2)
rownames(signif_data_bound) = (1:nrow(signif_data_bound))

write_csv(signif_data_bound, 'A:/cellphone_1_09_19_2025/whole_output/signif_data_bound.csv')