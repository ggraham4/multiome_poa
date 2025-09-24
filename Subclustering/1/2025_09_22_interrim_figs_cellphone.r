
#Everything
{
  library(ggsignif)
  library(patchwork)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(CytoTRACE)
  library(BiocManager)
  library(lme4)
  library(enrichR)
  `%notin%` = Negate(`%in%`)
  library(car)
  library(emmeans)
  library(glmGamPoi)
  library(scran)
  library(parallel)
    library(parallel)

  library(Seurat)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(emmeans)
  library(openxlsx)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
  library(biomaRt)
  library(hdWGCNA)
  
  clown_go = readRDS("Functions/clown_go2")  
  mecd = readRDS("Functions/mean_expression_cluster_data.rds")
    mecp = readRDS("Functions/mean_expression_cluster_plot.rds")

  `%notin%` <- Negate(`%in%`)
  
  geneNamer = function(gene){
  names = read.csv('Reference/genes updated.csv')
  
  name = names$NIH_description[names$NIH_accession==gene][1]
  
  if(is.na(name)){name = gene}
  return(name)
}

  
  #define functions
p_annotate <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else if (p_value < 0.1) {
    return(paste0("p = ", round(p_value, 3)))
  } else {
    return("NS")
  }
}

  plot_module <-  function(module){
  #ur supposed to use ME not hME
  me_subset=MEs%>%
    group_by(individual, Status)%>%
    summarize(mean_module = mean(.data[[module]]),
              se_module = sd(.data[[module]])/sqrt(n()))
  
  me_subset$Status = factor(me_subset$Status, levels = c('NRM','M',"D",'E','NF','F'))
  
  plot <- ggplot(me_subset, aes(x = Status, y = mean_module, color = Status))+
    geom_boxplot(aes(group = Status, fill = Status), alpha = 0.25, outlier.shape = NA)+
    geom_pointrange(aes(x = Status, y = mean_module,
                        ymin = mean_module - se_module,
                        ymax = mean_module +se_module),
                    position = position_jitterdodge(3))+
    theme_classic()+
    labs(x  ='Status', y = module)+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
  
}

plot_module_cluster <-  function(module){
  #ur supposed to use ME not hME
  me_subset=MEs%>%
    group_by(individual, sub.cluster)%>%
    summarize(mean_module = mean(.data[[module]]),
              se_module = sd(.data[[module]])/sqrt(n()))
  
  plot <- ggplot(me_subset, aes(x = sub.cluster, y = mean_module, color = sub.cluster))+
    geom_boxplot(alpha = 0.25, outlier.shape = NA)+
    geom_pointrange(aes(x = sub.cluster, y = mean_module,
                        ymin = mean_module - se_module,
                        ymax = mean_module +se_module),
                    position = position_jitterdodge(1))+
    theme_classic()+
    labs(x  ='sub.cluster', y = module)+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
  
}

prop_cluster_plot=
function(object, gene, cluster, clustering = 'final_clusters'){
    library(stringr)
    library(forcats)
      options(dplyr.summarise.inform = FALSE)

    counts <- t(as.matrix(object@assays$RNA$data[,object@meta.data[[clustering]] == cluster]))
  Counts_of_interest <- as.data.frame(counts[,gene]>0) #binarize
    Counts_of_interest$expression <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
    Counts_of_interest$Status <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]

  results <-Counts_of_interest%>%
    group_by(individual, Status)%>%
    summarize(mean = mean(expression),
              se = sd(expression)/sqrt(n()))
  results$Sex <- results$Status
  
  
  results$Sex <- str_sub(results$individual, -1)
  results$Sex[results$individual == 'T17D'] = 'NF'
  results$Sex[results$individual == 'A12D'] = 'E'
  results$Sex[results$individual == 'T11D'] = 'E'
  results$Sex[results$individual == 'GH'] = 'NRM'
  
  results$factor <- ifelse(results$Sex == "NRM", 0, NA)
  results$factor <- ifelse(results$Sex == "M", 1, results$factor)
  results$factor <- ifelse(results$Sex == "D", 2, results$factor)
  results$factor <- ifelse(results$Sex == "E", 3, results$factor)
  results$factor <- ifelse(results$Sex == "NF", 4, results$factor)
  results$factor <- ifelse(results$Sex == "F", 5, results$factor)
  results$individual <- fct_reorder(results$individual, results$factor)
  
  results$Sex <- factor(results$Sex, levels = c('NRM','M','D','E','NF','F'))
  plot <- ggplot(results, aes(x = individual, y = mean, color = Sex))+
    geom_boxplot(aes(group = Sex, fill = Sex), alpha = 0.25, outlier.shape = NA)+
    geom_point()+
    geom_pointrange(aes(x = individual, y = mean, ymin = mean - se, ymax = mean+se))+
    theme_classic()+
    labs(x  ='FishID', y = '% Expressing', title = paste0(gene,'_cluster_',cluster))+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
}




#obj = readRDS("A:/optimal_clustering_05_06_2025/RNA_object_human_names.rds")
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")


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
  
  data_for_plot$Phase  = (data_for_plot$Status)
  data_for_plot$Phase = ifelse(data_for_plot$Phase == 'D', 'I', data_for_plot$Phase)
  data_for_plot$Phase = ifelse(data_for_plot$Phase == 'E', 'LI', data_for_plot$Phase )
  data_for_plot$Phase <- factor(data_for_plot$Phase, levels = c('NRM',
                                          'M',
                                          'I',
                                          'LI',
                                          'NF',
                                          'F'))
  

  plot = ggplot(data_for_plot, aes(x = Phase, y = mean, 
                                   fill = Phase))+
    geom_point(position = position_jitterdodge(1))+
    geom_boxplot(alpha = 0.25, outlier.shape = NA)+
    theme_minimal()+
    labs(x = 'Status', y = 'Interaction Mean',
         title = paste0(LRpair, ' ', CCpair),
         subtitle = paste0('Phase q.value =',status_q_value ))
  
  return(plot)
  
}



#individuals
individuals = unique(obj$individual[!obj$individual %in% 'GH'])

##whole
#individual data

whole_individual_data = list()
for(i in individuals){
  message(paste0("Individual  ", which(individuals ==i), ' of ', length(individuals)))
  file_name = paste0(whole_base_path,i,"/degs_analysis_means_",i,".txt")
  whole_means = read_tsv(file_name)
  whole_individual_data[[i]]= whole_means
}



signif_data_bound  =read.csv('Subclustering/1/CellPhone/Whole/signif_data_bound.csv')
coalesced_list = readRDS("Subclustering/1/CellPhone/Whole/coalesced_list_09_19_2025.RDS")
signif_only_list = readRDS("Subclustering/1/CellPhone/Whole/signif_only_list_09_19_2025.RDS") 
relevant_interaction = read.csv("Subclustering/1/CellPhone/Whole/degs_analysis_relevant_interactions_09_19_2025_155229.txt", sep ='\t')

# 1_3_ signif 
signif_1_3 = signif_data_bound[signif_data_bound$main_effect_q.value<0.1 & 
                                 grepl('1_3',signif_data_bound$cluster_pair),]
plot_interaction('RTN4R_ADGRB1', '16|1_3','whole')
plot_interaction('RTN4R_ADGRB1', '3|1_3','whole')
plot_interaction('RTN4R_ADGRB1', '5|1_3','whole')

signif_1_2 = signif_data_bound[signif_data_bound$main_effect_q.value<0.1 & 
                                 grepl('1_2',signif_data_bound$cluster_pair),]

mecp(rgc, 'fgfbp3', '1_2', 'sub.cluster')

plot_interaction('KISS1_KISS1R', '10|1_3','whole')
plot_interaction('GABA_byGAD1_and_SLC6A13_GABA-A_a4b3g2S_complex', '1_3|10','whole')




measures =read.csv("Measures/all_data.csv")
}

plot_interaction2 = function(LRpair, CCpair, whole_sub, 
                          signif_dm = NULL,
                          signif_df = NULL,
                          signif_mf = NULL){
  
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
  
  # Transform status to phase
  data_for_plot$Phase = case_when(
    data_for_plot$Status == 'D' ~ 'I',
    data_for_plot$Status == 'E' ~ 'LI',
    TRUE ~ as.character(data_for_plot$Status)
  )
  
  data_for_plot$Phase <- factor(data_for_plot$Phase, 
                               levels = c('NRM', 'M', 'I', 'LI', 'NF', 'F'))
  
  # Calculate plot limits
  plot_lower_lim = min(data_for_plot$mean, na.rm = TRUE) * 0.9
  plot_upper_lim = max(data_for_plot$mean, na.rm = TRUE) * 1.2
  plot_signif_lower = max(data_for_plot$mean, na.rm = TRUE) * 1.05
  plot_signif_upper = max(data_for_plot$mean, na.rm = TRUE) * 1.25
  
  # Process significance values
  if(!is.null(signif_dm)) signif_dm = p_annotate(signif_dm)
  if(!is.null(signif_df)) signif_df = p_annotate(signif_df)
  if(!is.null(signif_mf)) signif_mf = p_annotate(signif_mf)
  
  # Calculate text sizes
  textsize_dm = if(!is.null(signif_dm)) ifelse(grepl("\\*", signif_dm), 6, 3) else 3
  textsize_df = if(!is.null(signif_df)) ifelse(grepl("\\*", signif_df), 6, 3) else 3
  textsize_mf = if(!is.null(signif_mf)) ifelse(grepl("\\*", signif_mf), 6, 3) else 3
  
  # Adjust upper limit if mf is significant
  if(!is.null(signif_mf) && signif_mf != 'NS'){
    plot_upper_lim = max(data_for_plot$mean, na.rm = TRUE) * 1.4
  }
  
  # Create plot
  plot = ggplot(data_for_plot, aes(x = Phase, y = mean, fill = Phase)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(1), size = 1) +
    labs(y = 'Mean Signaling') +
    ggtitle(paste0(LRpair)) +
    labs(subtitle = CCpair)+
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 8)) +
    theme(legend.position = 'none') +
    ylim(plot_lower_lim, plot_upper_lim)
  
  # Add significance annotations
  if(!is.null(signif_dm)){
    plot = plot +
      geom_signif(xmin = c(1.0),
                  xmax = c(1.9),
                  y_position = c(plot_signif_lower),
                  annotation = c(signif_dm), 
                  color = "black",
                  tip_length = c(0,0),
                  textsize = textsize_dm)
  }
  
  if(!is.null(signif_df)){
    plot = plot +
      geom_signif(xmin = c(2.1),
                  xmax = c(5),
                  y_position = c(plot_signif_lower),
                  annotation = c(signif_df), 
                  color = "black",
                  tip_length = c(0,0),
                  textsize = textsize_df)
  }
  
  if(!is.null(signif_mf) && signif_mf != 'NS'){
    plot = plot +
      geom_signif(xmin = c(1),
                  xmax = c(5),
                  y_position = c(plot_signif_upper),
                  annotation = c(signif_mf), 
                  color = "black",
                  tip_length = c(0,0),
                  textsize = textsize_mf)
  }
  
  return(plot)
}
test_ =plot_interaction2('KISS1_KISS1R', '10|1_3','whole', 0.1, 0.05, 0.001)

#ggsave(plot = test_,
#       file = 'test.tiff',
#       device = "tif",
#       units = "in",
#       width = 2.5,
#       height = 2.5,
#       path = "Subclustering/1/interrim_figs")

kiss_10_1_3 =plot_interaction2('KISS1_KISS1R',
                               '10|1_3',
                               'whole',
                               signif_1_3$d_m_p.value[signif_1_3$interacting_pair=='KISS1_KISS1R'],
                                signif_1_3$d_f_p.value[signif_1_3$interacting_pair=='KISS1_KISS1R'],
                                signif_1_3$f_m_p.value[signif_1_3$interacting_pair=='KISS1_KISS1R'])

ggsave(plot = kiss_10_1_3,
       file = 'kiss_10_1_3.tiff',
       device = "tif",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")

WNT5B_FZD3_LRP5_14_1_1 =plot_interaction2('WNT5B_FZD3_LRP5',
                               '14|1_1',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='WNT5B_FZD3_LRP5' & 
                                                               signif_data_bound$cluster_pair== '14|1_1'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='WNT5B_FZD3_LRP5'& 
                                                                signif_data_bound$cluster_pair== '14|1_1'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='WNT5B_FZD3_LRP5'& 
                                                                signif_data_bound$cluster_pair== '14|1_1'])

ggsave(plot = WNT5B_FZD3_LRP5_14_1_1,
       file = 'WNT5B_FZD3_LRP5_14_1_1.tiff',
       device = "tif",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")


NRG3_ERBB4_1_5_1_1 =plot_interaction2('NRG3_ERBB4',
                               '1_5|1_1',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='NRG3_ERBB4' & 
                                                               signif_data_bound$cluster_pair== '1_5|1_1'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='NRG3_ERBB4'& 
                                                                signif_data_bound$cluster_pair== '1_5|1_1'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='NRG3_ERBB4'& 
                                                                signif_data_bound$cluster_pair== '1_5|1_1'])

ggsave(plot = NRG3_ERBB4_1_5_1_1,
       file = 'NRG3_ERBB4_1_5_1_1.tiff',
       device = "tif",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")

DHEA_1_8_1_3 =plot_interaction2('Dehydroepiandrosterone_bySTS_PPARA',
                               '18|1_3',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA' & 
                                                               signif_data_bound$cluster_pair== '18|1_3'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA'& 
                                                                signif_data_bound$cluster_pair== '18|1_3'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA'& 
                                                                signif_data_bound$cluster_pair== '18|1_3'])+
  labs(title ='DHEA_bySTS_PPARA' )

ggsave(plot = DHEA_1_8_1_3,
       file = 'DHEA_1_8_1_3.tiff',
       device = "tif",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")

Nogo_16_1_3 =plot_interaction2('RTN4R_ADGRB1',
                               '16|1_3',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='RTN4R_ADGRB1' & 
                                                               signif_data_bound$cluster_pair== '16|1_3'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='RTN4R_ADGRB1'& 
                                                                signif_data_bound$cluster_pair== '16|1_3'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='RTN4R_ADGRB1'& 
                                                                signif_data_bound$cluster_pair== '16|1_3'])

ggsave(plot = Nogo_16_1_3,
       file = 'Nogo_16_1_3',
       device = "tif",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")

Nogo_3_1_3 =plot_interaction2('RTN4R_ADGRB1',
                               '3|1_3',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='RTN4R_ADGRB1' & 
                                                               signif_data_bound$cluster_pair== '3|1_3'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='RTN4R_ADGRB1'& 
                                                                signif_data_bound$cluster_pair== '3|1_3'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='RTN4R_ADGRB1'& 
                                                                signif_data_bound$cluster_pair== '3|1_3'])

ggsave(plot = Nogo_3_1_3,
       file = 'Nogo_3_1_3',
       device = "tif",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")

Nogo_5_1_3 =plot_interaction2('RTN4R_ADGRB1',
                               '5|1_3',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='RTN4R_ADGRB1' & 
                                                               signif_data_bound$cluster_pair== '5|1_3'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='RTN4R_ADGRB1'& 
                                                                signif_data_bound$cluster_pair== '5|1_3'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='RTN4R_ADGRB1'& 
                                                                signif_data_bound$cluster_pair== '5|1_3'])

ggsave(plot = Nogo_5_1_3,
       file = 'Nogo_5_1_3',
       device = "tif",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")

KISS1_KISS1R_10_1_3 =plot_interaction2('KISS1_KISS1R',
                               '10|1_3',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='KISS1_KISS1R' & 
                                                               signif_data_bound$cluster_pair== '10|1_3'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='KISS1_KISS1R'& 
                                                                signif_data_bound$cluster_pair== '10|1_3'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='KISS1_KISS1R'& 
                                                                signif_data_bound$cluster_pair== '10|1_3'])

ggsave(plot = KISS1_KISS1R_10_1_3,
       file = 'KISS1_KISS1R_10_1_3',
       device = "tif",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")

GABA_1_3_10 =plot_interaction2('GABA_byGAD1_and_SLC6A11_GABA-A_a4b3d_complex',
                               '1_3|10',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='GABA_byGAD1_and_SLC6A11_GABA-A_a4b3d_complex' & 
                                                               signif_data_bound$cluster_pair== '1_3|10'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='GABA_byGAD1_and_SLC6A11_GABA-A_a4b3d_complex'& 
                                                                signif_data_bound$cluster_pair== '1_3|10'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='GABA_byGAD1_and_SLC6A11_GABA-A_a4b3d_complex'& 
                                                                signif_data_bound$cluster_pair== '1_3|10'])+
  labs(title = 'GABA_GABA-A_a4b3d')

ggsave(plot = GABA_1_3_10,
       file = 'GABA_1_3_10',
       device = "tif",
       units = "in",
       width = 2.5,
       height = 3,
       path = "Subclustering/1/interrim_figs")

plot_interaction2('Glutamate_byGLS_and_SLC17A7_GRM3',
                               '22|1_2',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='Glutamate_byGLS_and_SLC17A7_GRM3' & 
                                                               signif_data_bound$cluster_pair== '22|1_2'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='Glutamate_byGLS_and_SLC17A7_GRM3'& 
                                                                signif_data_bound$cluster_pair== '22|1_2'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='Glutamate_byGLS_and_SLC17A7_GRM3'& 
                                                                signif_data_bound$cluster_pair== '22|1_2'])



dhea_18_14 = plot_interaction2('Dehydroepiandrosterone_bySTS_PPARA',
                               '18|14',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA' & 
                                                               signif_data_bound$cluster_pair== '18|14'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA'& 
                                                                signif_data_bound$cluster_pair== '18|14'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA'& 
                                                                signif_data_bound$cluster_pair== '18|14'])

dhea_18_16 = plot_interaction2('Dehydroepiandrosterone_bySTS_PPARA',
                               '18|16',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA' & 
                                                               signif_data_bound$cluster_pair== '18|16'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA'& 
                                                                signif_data_bound$cluster_pair== '18|16'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA'& 
                                                                signif_data_bound$cluster_pair== '18|16'])

dhea_18_6 = plot_interaction2('Dehydroepiandrosterone_bySTS_PPARA',
                               '18|6',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA' & 
                                                               signif_data_bound$cluster_pair== '18|6'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA'& 
                                                                signif_data_bound$cluster_pair== '18|6'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='Dehydroepiandrosterone_bySTS_PPARA'& 
                                                                signif_data_bound$cluster_pair== '18|6'])

LRFN5_PTPRF_5_6 = plot_interaction2('LRFN5_PTPRF',
                               '5|6',
                               'whole',
                               signif_data_bound$d_m_p.value[signif_data_bound$interacting_pair=='LRFN5_PTPRF' & 
                                                               signif_data_bound$cluster_pair== '5|6'],
                                signif_data_bound$d_f_p.value[signif_data_bound$interacting_pair=='LRFN5_PTPRF'& 
                                                                signif_data_bound$cluster_pair== '5|6'],
                                signif_data_bound$f_m_p.value[signif_data_bound$interacting_pair=='LRFN5_PTPRF'& 
                                                                signif_data_bound$cluster_pair== '5|6'])
LRFN5_PTPRF_5_6




