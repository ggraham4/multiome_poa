library(lme4)
library(tidyverse)
library(Seurat)

human_named = readRDS("A:/optimal_clustering_05_06_2025/RNA_object_human_names.rds")

## read in reference data
example_data = read_tsv("A:/cellphonedb_05_12_2025/degs_analysis_means_05_12_2025_222054.txt")
coalesced_data = readRDS('A:/cellphonedb_05_12_2025/coalesced_list_05_20_2025.RDS')
base_path = 'A:/cellphonedb_05_12_2025/'

unique_interactions = unique(example_data$interacting_pair)
unique_cell_cell_pairs = colnames(example_data[,14:ncol(example_data)])

individuals = unique(human_named$individual)
individuals = individuals[individuals!= 'GH']

## read in significance data
significant_interactions_reference = read_tsv("A:/cellphonedb_05_12_2025/degs_analysis_relevant_interactions_05_12_2025_222054.txt")

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

saveRDS(new_list, 'A:/cellphonedb_05_12_2025/signif_only_list_05_21_2025.RDS')






