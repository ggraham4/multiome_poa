
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

deg_plotter = function(object = rgc, 
                       gene, 
                       cluster , 
                       clustering='sub.cluster',
                       signif_dm  ,
                       signif_df ,
                       signif_mf ,
                       singular=F,
                       common_name){
  set.seed(10)
  singular = ifelse(singular == T, 'Singular', '')
  meta= object@meta.data
  meta$gene = object@assays$RNA$data[gene,]
  
  meta_grouped_and_sded = meta%>%
    filter(Phase != 'NRM' & !!sym(clustering) == cluster) %>%
    group_by(individual, Phase)%>%
    summarize(mean_gene = mean(gene),
              se_gene = sd(gene)/sqrt(n()))
  
  plot_lower_lim = min(meta_grouped_and_sded$mean_gene -meta_grouped_and_sded$se_gene )
  plot_upper_lim= max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.2
  plot_signif_lower = max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.05
  plot_signif_upper = max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.25
  

  
    signif_dm = p_annotate(signif_dm)
    signif_df = p_annotate(signif_df)
    signif_mf = p_annotate(signif_mf)
  
  textsize_dm = ifelse(grepl("\\*", signif_dm), 6, 3)  
  textsize_df = ifelse(grepl("\\*", signif_df), 6, 3)  
  textsize_mf = ifelse(grepl("\\*", signif_mf), 6, 3)    

  
      plot_upper_lim= ifelse(signif_mf!= 'NS', max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.4,plot_upper_lim )

  
plot = ggplot(meta_grouped_and_sded, aes(x = Phase, y = mean_gene,fill = Phase))+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_pointrange(aes(x = Phase,
                      y = mean_gene,
                      ymin = mean_gene-se_gene,
                      ymax= mean_gene+se_gene),
                  position = position_jitterdodge(1), 
                  size = 0.2
                  )+
  labs(y = 'Mean +/- SE Expression', subtitle = singular)+
  ggtitle(paste0(common_name, ': ', cluster))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size =12),
        plot.subtitle = element_text(hjust = 0.5, size =8))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.0),
              xmax = c(1.9),
              y_position = c(plot_signif_lower),
              annotation =c(signif_dm), 
              color = "black",
              tip_length = c(0,0),
              textsize=textsize_dm)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(plot_signif_lower),
              annotation =c(signif_df), 
              color = "black",
              tip_length = c(0,0),
              textsize=textsize_df)+
  ylim(plot_lower_lim, plot_upper_lim)
   
 if(signif_mf != 'NS'){
   plot  <- plot+
      geom_signif(xmin = c(1),
              xmax = c(5),
              y_position = c(plot_signif_upper),
              annotation =c(signif_mf), 
              color = "black",
              tip_length = c(0,0),
              textsize=textsize_mf)
    
  }

return(plot)
  
}


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


}
#read in data
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
whole_base_path = "/Users/ggraham/Desktop/multiome_poa/Subclustering/6/CellPhone/whole_output/"
whole_significant = read.csv("Subclustering/6/CellPhone/whole_output/signif_data_bound.csv")
measures = read.xlsx("Reference/Complete Data Frame (Hormones, Behavior, Size, Gonads) GG.xlsx")
whole_6_signif = read.csv("Subclustering/6/CellPhone/whole_output/signif_data_bound.csv")
whole_6_signif_linear = subset(whole_6_signif, test_type =='Linear Regression' & 
                                 main_effect_q.value <0.1)

individuals = unique(obj$individual)
whole_individual_data = list()
for(i in individuals){
  if(i =='GH'){next}
  message(paste0("Individual  ", which(individuals ==i), ' of ', length(individuals)))
  file_name = paste0(whole_base_path,i,"/degs_analysis_means_",i,".txt")
  whole_means = read_tsv(file_name)
  whole_individual_data[[i]]= whole_means
}

# analysis ----
# for now, I ignore AR signaling because not thought to be bioactive in teleosts

plot_interaction('BDNF_NTRK2', '17|6_0', 'whole')
#bdnf decreases in dominants and females

plot_interaction('NECTIN2_NECTIN3', '6_0|8', 'whole')
plot_interaction('NECTIN2_NECTIN3', '6_0|22', 'whole')
# in physical contact with 8 and 22, strongly increased in dominants

plot_interaction('LRRTM4_NRXN2', '6_0|2', 'whole')
# decreased interaction with oligodendrocytes

plot_interaction('NRG2_ERBB3', '6_0|2', 'whole')
plot_interaction('NRG2_ERBB3', '6_2|2', 'whole') #interesting how it remains increased in 6_0 but not 6_2
#increased this pro differentiation signal in oligos huh

### all of the signaling is happening to 6_0

#nectin2 ensures cell survival
plot_interaction('CD226_NECTIN2', '1|6_0', 'whole') 
plot_interaction('CD226_NECTIN2', '2|6_0', 'whole') 

# ptprf promotes cell adhesion, axon regeneration, neurite outgrowth

#6_0
plot_interaction('LRFN5_PTPRF', '9|6_0', 'whole') 
plot_interaction('LRFN5_PTPRF', '10|6_0', 'whole') 
plot_interaction('LRFN5_PTPRF', '12|6_0', 'whole') 
plot_interaction('LRFN5_PTPRF', '13|6_0', 'whole') 

plot_interaction('LRRC4C_PTPRF', '16|6_0', 'whole') 
plot_interaction('LRRC4C_PTPRF', '17|6_0', 'whole') 
plot_interaction('LRRC4C_PTPRF', '18|6_0', 'whole') ### DHEA cluster
plot_interaction('LRRC4C_PTPRF', '5|6_0', 'whole') 
plot_interaction('LRRC4C_PTPRF', '9|6_0', 'whole') 

plot_interaction('IL1RAP_PTPRF', '7|6_0', 'whole') 


plot_interaction('IL1RAP_PTPRF', '6_2|6_0', 'whole') 

#6_2
plot_interaction('LRFN5_PTPRF', '13|6_2', 'whole') 
plot_interaction('LRFN5_PTPRF', '6_1|6_2', 'whole') 
plot_interaction('LRRC4C_PTPRF', '0|6_2', 'whole') 
plot_interaction('LRRC4C_PTPRF', '13|6_2', 'whole') 
plot_interaction('LRRC4C_PTPRF', '3|6_2', 'whole') 
plot_interaction('LRRC4C_PTPRF', '8|6_2', 'whole') 

plot_interaction('IL1RAP_PTPRF', '7|6_2', 'whole') 

# GLUTAMATE
#6_1
plot_interaction('Glutamate_byGLS_and_SLC1A6_GRM1', '5|6_1', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC1A1_GRM1', '0|6_1', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC17A7_GRM1', '24|6_1', 'whole') 

#6_0
plot_interaction('Glutamate_byGLS_and_SLC1A6_Glutamate_Kainate_2_5_complex', '11|6_0', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC1A6_Glutamate_Kainate_2_5_complex', '5|6_0', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC1A1_Glutamate_Kainate_2_5_complex', '19|6_0', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC1A1_Glutamate_Kainate_2_5_complex', '6_0|6_0', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC1A1_Glutamate_Kainate_2_5_complex', '8|6_0', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC1A6_GRIK2', '5|6_0', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC17A7_GRIA4', '24|6_0', 'whole') # unique


# GABA
plot_interaction('GABA_byGAD2_and_SLC6A6_GABA-A_a1b3g2S_complex', '6_0|25', 'whole')
# increased


### now lets look at the fishers
whole_6_signif_fisher = subset(whole_6_signif, test_type =='Fisher' & 
                                 main_effect_q.value <0.1)
#christ thats a lot of GABA
plot_interaction('PENK_OPRM1', '6_3|17', 'whole') 
plot_interaction('SST_SSTR5', '6_2|24', 'whole') # probably not real, but this may also reflect hte loss in females


## who does 10 talk to

#not a lot going on with 10 - 6 BUT, 10 does communicate
whole_6_signif$sender = sapply(strsplit(whole_6_signif$cluster_pair, "\\|"), `[`, 1)
whole_6_signif$receiver = sapply(strsplit(whole_6_signif$cluster_pair, "\\|"), `[`, 2)

table(whole_6_signif$receiver[whole_6_signif$sender==10])%>%sort()
# bro is yapping to cluster 8, stronger evidence for 6_0 than 6_3, so my theory is taking a hit there

table(whole_6_signif$sender[whole_6_signif$receiver==10])%>%sort()
#likewuse 0 is talking to it, but 61 is making itself known as well as 5

# FINDING 1 ----
"Sex change leads to an increase in glut signaling to 6_1, and a decrease to 6_0
"
plot_interaction('SLITRK5_PTPRS', '25|24', 'whole') 

plot_interaction('LRRTM4_NRXN2', '7|18', 'whole') 

plot_interaction('EFNA3_EPHA3', '13|4', 'whole') 

plot_interaction('BTC_ERBB3', '19|2', 'whole') 
plot_interaction('NRG2_ERBB3', '9|2', 'whole') 

plot_interaction('TENM4_ADGRL3', '9|16', 'whole') 

plot_interaction('GABA_byGAD2_and_SLC6A11_GABA-A_a1b2_complex', '0|3', 'whole') 
plot_interaction('GABA_byGAD2_and_SLC6A6_GABA-A_a1b3g2S_complex', '19|14', 'whole') 

plot_interaction('GABA_byGAD2_and_SLC6A6_GABA-A_a3b3g2S_complex', '8|17', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC17A7_Glutamate_Kainate_1_4_complex', '24|4', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC17A7_GRM1', '24|3', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC17A7_GRM1', '24|24', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC17A7_GRM1', '24|13', 'whole') 

plot_interaction('Glutamate_byGLS_and_SLC17A7_GRM5', '24|12', 'whole') 
plot_interaction('Glutamate_byGLS_and_SLC17A7_GRM7', '22|11', 'whole') 

plot_interaction('IL1RAP_PTPRF', '7|9', 'whole') 
plot_interaction('LRRC4C_PTPRF', '9|6_0', 'whole') 

plot_interaction('LRRC4C_PTPRF', '8|9', 'whole') 
plot_interaction('LRRC4C_PTPRF', '8|6_2', 'whole') 
plot_interaction('LRRC4C_PTPRF', '5|6_0', 'whole') 
plot_interaction('LRRC4C_PTPRF', '3|9', 'whole') 
plot_interaction('LRRC4C_PTPRF', '3|6_2', 'whole') 
plot_interaction('LRRC4C_PTPRF', '18|6_0', 'whole') 
plot_interaction('LRRC4C_PTPRF', '17|6_0', 'whole') 
plot_interaction('LRRC4C_PTPRF', '16|6_0', 'whole') 

plot_interaction('LRRC4C_PTPRF', '13|9', 'whole') 
plot_interaction('LRRC4C_PTPRF', '13|6_2', 'whole') 

plot_interaction('LRRC4C_PTPRF', '0|9', 'whole') 

plot_interaction('LRRC4C_PTPRF', '0|6_2', 'whole') 

plot_interaction('LRFN5_PTPRF', '6_1|6_2', 'whole') 

plot_interaction('LRRC4C_PTPRF', '5|9', 'whole') 
plot_interaction('LRRC4C_PTPRF', '5|6_2', 'whole') 
plot_interaction('LRFN5_PTPRF', '13|6_0', 'whole') 
plot_interaction('LRFN5_PTPRF', '12|6_0', 'whole') 

plot_interaction('atRetinoicAcid_byALDH1A2_RAreceptor_RARB', '11|19', 'whole') 
plot_interaction('WNT2_SFRP1', '10|17', 'whole') 

plot_interaction('GABA_byGAD1_and_SLC6A6_GABA-A_a3b3g2S_complex', '6_0|6_3', 'whole') 

plot_interaction('GABA_byGAD1_and_SLC6A6_GABA-A_a3b3g2S_complex', '6_1|6_3', 'whole') 
plot_interaction('GABA_byGAD1_and_SLC6A6_GABA-A_a3b3g2S_complex', '14|6_3', 'whole') 



table(whole_6_signif_fisher$sender, whole_6_signif_fisher$receiver)
table(whole_6_signif_fisher$receiver)
table(whole_6_signif_fisher$sender)

