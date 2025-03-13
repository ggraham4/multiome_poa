### Fisher's exact test

data = read_csv('A:/CellPhoneDB 030225/signif_data_bound.csv')

signif = data$interacting_pair[data$main_effect_q.value<0.05]
all_interacting = data$interacting_pair

successes =  length(data$interacting_pair[data$main_effect_q.value<0.05 & data$interacting_pair == 'Estradiol_byCYP19A1_ESR2'])
failures = length(data$interacting_pair[data$main_effect_q.value<0.05 & data$interacting_pair != 'Estradiol_byCYP19A1_ESR2'])
all_successes_possible = length(data$interacting_pair[data$interacting_pair == 'Estradiol_byCYP19A1_ESR2'])
all_failures_possible = length(data$interacting_pair[data$interacting_pair != 'Estradiol_byCYP19A1_ESR2'])

fisher_test_pair = function(signif_data, term) {
  contingency_table = matrix(c(
    sum(signif_data$main_effect_q.value < 0.05 & signif_data$interacting_pair == term),  # Success in significant
    sum(signif_data$main_effect_q.value >= 0.05 & signif_data$interacting_pair == term), # Success in non-significant
    sum(signif_data$main_effect_q.value < 0.05 & signif_data$interacting_pair != term),  # Failure in significant
    sum(signif_data$main_effect_q.value >= 0.05 & signif_data$interacting_pair != term)  # Failure in non-significant
  ), nrow=2, byrow=TRUE)
  
  test = fisher.test(contingency_table, alternative = "greater")$p.value
  return(test)
}

fisher_test_pair(data, 'Estradiol_byCYP19A1_ESR2')

fisher_test_pair(data, 'SLITRK4_PTPRS') # not significant

fisher_test_pair(data, 'WNT2_SFRP1')

fisher_test_pair(data, 'LRRTM4_NRXN2')

fisher_test_pair(data, 'Glutamate_byGLS_and_SLC1A6_GRM1') # not significant

fisher_test_pair(data, 'Glutamate_byGLS_and_SLC17A7_GRM7')
