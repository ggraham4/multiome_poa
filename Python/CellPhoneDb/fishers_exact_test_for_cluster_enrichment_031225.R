### Fisher's exact test

data = read_csv('A:/CellPhoneDB 030225/signif_data_bound.csv')

data[c("sender", "receiver")] <- do.call(rbind, strsplit(data$cluster_pair, "\\|"))


signif = data$interacting_pair[data$main_effect_q.value<0.05]
all_interacting = data$interacting_pair

#### senders #####

fisher_test_pair = function(signif_data, cluster) {
  contingency_table = matrix(c(
    sum(signif_data$main_effect_q.value < 0.05 & signif_data$sender == cluster),  # Success in significant
    sum(signif_data$main_effect_q.value >= 0.05 & signif_data$sender == cluster), # Success in non-significant
    sum(signif_data$main_effect_q.value < 0.05 & signif_data$sender != cluster),  # Failure in significant
    sum(signif_data$main_effect_q.value >= 0.05 & signif_data$sender != cluster)  # Failure in non-significant
  ), nrow=2, byrow=TRUE)
  
  test = fisher.test(contingency_table, alternative = "greater")$p.value
  return(test)
}

for(i in unique(data$sender[data$main_effect_q.value<0.05])){
  print(i)
  print(fisher_test_pair(data, i))
  
}

# 2 and 12 are overrepresented as a sender
#### receiver ####
fisher_test_pair_receiver = function(signif_data, cluster) {
  contingency_table = matrix(c(
    sum(signif_data$main_effect_q.value < 0.05 & signif_data$receiver == cluster),  # Success in significant
    sum(signif_data$main_effect_q.value >= 0.05 & signif_data$receiver == cluster), # Success in non-significant
    sum(signif_data$main_effect_q.value < 0.05 & signif_data$receiver != cluster),  # Failure in significant
    sum(signif_data$main_effect_q.value >= 0.05 & signif_data$receiver != cluster)  # Failure in non-significant
  ), nrow=2, byrow=TRUE)
  
  test = fisher.test(contingency_table, alternative = "greater")$p.value
  return(test)
}


for(i in unique(data$receiver[data$main_effect_q.value<0.05])){
  print(i)
  print(fisher_test_pair_receiver(data, i))
  
}


# 8 overrepresented as receiver:

### combination #####

fisher_test_pair_combination = function(signif_data, cluster_pair) {
  contingency_table = matrix(c(
    sum(signif_data$main_effect_q.value < 0.05 & signif_data$cluster_pair == cluster_pair),  # Success in significant
    sum(signif_data$main_effect_q.value >= 0.05 & signif_data$cluster_pair == cluster_pair), # Success in non-significant
    sum(signif_data$main_effect_q.value < 0.05 & signif_data$cluster_pair != cluster_pair),  # Failure in significant
    sum(signif_data$main_effect_q.value >= 0.05 & signif_data$cluster_pair != cluster_pair)  # Failure in non-significant
  ), nrow=2, byrow=TRUE)
  
  test = fisher.test(contingency_table, alternative = "greater")$p.value
  return(test)
}

for(i in unique(data$cluster_pair[data$main_effect_q.value<0.05])){
  print(i)
  print(fisher_test_pair_combination(data, i))
  
}

# they are all overreporesented because there are 900 cluster pair combinations

