### hypergeometric test of cellphonedb terms
data = read_csv('A:/CellPhoneDB 030225/signif_data_bound.csv')

signif = data$interacting_pair[data$main_effect_q.value<0.05]
all_interacting = data$interacting_pair

successes =  length(data$interacting_pair[data$main_effect_q.value<0.05 & data$interacting_pair == 'Estradiol_byCYP19A1_ESR2'])
failures = length(data$interacting_pair[data$main_effect_q.value<0.05 & data$interacting_pair != 'Estradiol_byCYP19A1_ESR2'])
all_successes_possible = length(data$interacting_pair[data$interacting_pair == 'Estradiol_byCYP19A1_ESR2'])
all_failures_possible = length(data$interacting_pair[data$interacting_pair != 'Estradiol_byCYP19A1_ESR2'])

  
phyper(q = successes, m = all_successes_possible, n =  all_failures_possible, k = failures+successes, lower.tail = F, log.p = F)



hypergeom_pair = function(signif_data, term){
  
  all_interacting = signif_data$interacting_pair
  
  successes =  length(signif_data$interacting_pair[signif_data$main_effect_q.value<0.05 & signif_data$interacting_pair == term])
  failures = length(signif_data$interacting_pair[signif_data$main_effect_q.value<0.05 & signif_data$interacting_pair != term])
  all_successes_possible = length(signif_data$interacting_pair[signif_data$interacting_pair == term])
  all_failures_possible = length(signif_data$interacting_pair[signif_data$interacting_pair != term])
  
  
  test = phyper(q = successes, m = all_successes_possible, n =  all_failures_possible, k = failures+successes, lower.tail = F, log.p =F)
  
  return(test)
  
}

hypergeom_pair(data, 'Estradiol_byCYP19A1_ESR2')

hypergeom_pair(data, 'SLITRK4_PTPRS')

hypergeom_pair(data, 'WNT2_SFRP1')

hypergeom_pair(data, 'LRRTM4_NRXN2')

hypergeom_pair(data, 'Glutamate_byGLS_and_SLC1A6_GRM1')

hypergeom_pair(data, 'Glutamate_byGLS_and_SLC17A7_GRM7')

