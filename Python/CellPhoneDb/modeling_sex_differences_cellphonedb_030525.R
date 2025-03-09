library(lme4)
library(tidyverse)
library(Seurat)
library(car)
library(emmeans)

human_named = readRDS("C:/Users/Gabe/Desktop/RNA_object_human_names.rds")

## read in reference data
example_data = read_tsv("A:/CellPhoneDB 030225/degs_analysis_means_03_03_2025_224407.txt")
base_path = 'A:/CellPhoneDB 030225/'

unique_cell_cell_pairs = colnames(example_data[,14:ncol(example_data)])

individuals = unique(human_named$individual)
individuals = individuals[individuals!= 'GH']

signif_list = readRDS('A:/CellPhoneDB 030225/signif_only_list_030525.RDS')


unique_interactions = unique(names(signif_list))

signif_data = list()
for(interacting_pair in unique_interactions){
  message(paste0('Interaction ',which(interacting_pair ==unique_interactions), ' of ', length(unique_interactions)))
  
  full_data = signif_list[[interacting_pair]]
  full_data =full_data[full_data$Status %in% c('M','D','F'),]%>%
    na.omit()
  
  interacting_pair_data = data.frame()
  for(cluster_pair in unique(full_data$cluster_pair)){
    if(grepl(cluster_pair, '30') ==T | grepl(cluster_pair, '15')==T){next} ### skip clusters 15 and 31, they are no goof
    
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
                  failure = n() - success)
      
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

#theres a way to do this in parallel but I dont want to think to hard, i'll just let it run, make sure to save it
signif_data2 =signif_data

for(i in names(signif_data2)){
  signif_data2[[i]]$main_effect_q.value = p.adjust(signif_data2[[i]]$main_effect_p.value, 'fdr',nrow(signif_data2[[i]]))
  
  
}

# do p adjust within interacting pair not across, 72K will kill any p values

signif_data_bound = do.call(rbind, signif_data2)
rownames(signif_data_bound) =(1:nrow(signif_data_bound))

write_csv(signif_data_bound, 'A:/CellPhoneDB 030225/signif_data_bound.csv')

