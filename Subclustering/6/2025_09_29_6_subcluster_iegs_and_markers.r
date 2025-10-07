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
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")


calculateCoexpressedGenes = function(seurat_obj2, geneList, clustering = 'sub.cluster'){
  ### the goal is to split cells in a cluster into + and negative for the gene list,
  ### then find other genes coexpressed with those
  
  `%notin%` = Negate(`%in%`)
  seurat_obj <- seurat_obj2
  seurat_obj$clustering = seurat_obj[[clustering]]
  
  seurat_obj$positive = ifelse(seurat_obj@assays$RNA$data[geneList[1],]>0,T, NA)
  for(i in 2:length(geneList)){
    if(geneList[i] %in% rownames(seurat_obj@assays$RNA$data)){
      seurat_obj$positive = ifelse(seurat_obj@assays$RNA$data[geneList[i],]>0,T, seurat_obj$positive)
    }}
    seurat_obj$positive = ifelse(is.na(seurat_obj$positive), F, seurat_obj$positive)

  
  marker_all<- data.frame()
  for(cluster in unique(seurat_obj$clustering)){
    temp_obj = subset(seurat_obj, clustering == cluster)
    Idents(temp_obj) = 'positive'
    
    if(sum(temp_obj$positive)>3){
    
    markers <- FindAllMarkers(temp_obj,
                              assay = "RNA",
                              group.by = 'positive',
                              logfc.threshold = 0,
                              min.pct = 1/nrow(temp_obj@meta.data))
    markers$sub.cluster= cluster
    marker_all = rbind(markers, marker_all)
    }
  }
  #extract significant   genes
  marker_all_signif <- marker_all[marker_all$p_val_adj<0.05 &marker_all$cluster==T ,]
  
  #genes in over half of clusters
  marker_genes_half = marker_all_signif%>%
    group_by(gene)%>%
    summarize(n_clusters = n())%>%
  subset(n_clusters >= length(unique(seurat_obj$clustering))/2
  ) 
  
  message('Found ', nrow(marker_genes_half) - length(geneList), ' new markers')
  print(marker_genes_half$gene[marker_genes_half$gene %notin% geneList])
  
  all_markers = marker_genes_half$gene
  #module_score time
  #seurat_obj = AddModuleScore(seurat_obj, 
  #                            features = list(all_markers),
  #                            name = 'coexModule')
  
  return(all_markers)
}


sub_6 = FindSubCluster(obj, 6, graph.name='harmony.wsnn')
Idents(sub_6) <- 'sub.cluster'
sub_6 = subset(sub_6, final_clusters ==6)
sub_6$sub.cluster = factor(sub_6$sub.cluster, levels = c(paste0('6_',0:3)))
sub_6$Status = factor(sub_6$Status, levels = c('NRM','M',"D",'E','NF','F'))

DimPlot(sub_6)

marks_6 = FindAllMarkers(sub_6)

#iegs look like degs for 6_0

iegs <- c('npas4a', 'fosb','egr1')

sub_6 = AddModuleScore(sub_6, list(iegs), name = 'ieg')
DotPlot(sub_6, 'ieg1', group.by = 'sub.cluster')+
  coord_flip()

## lets do the coexpression module
iegs2 =calculateCoexpressedGenes(sub_6, iegs)
sub_6 = AddModuleScore(sub_6, list(iegs2), name = 'iegs2')
DotPlot(sub_6, 'iegs21', group.by = 'sub.cluster')+
  coord_flip()

# ok, it seems like 6_1 and 6_2 are much more IEG + than 6_0 and 6_3 both in terms of prop and expression
module_plot <- sub_6@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(mean_ieg = mean(iegs21),
            se_ieg = sd(iegs21)/sqrt(n()))

ggplot(module_plot, aes(x = sub.cluster, y = mean_ieg, fill = Status, shape = Status))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
    geom_point(position = position_jitterdodge(1))
## looks like females might be less active than males in 6_1

ieg_module_lmer = lmer(iegs21~Status+(1|individual), data = subset(sub_6@meta.data, sub.cluster =='6_1' & 
                                                                     Status %in% c('M','D','F')))
car::Anova(ieg_module_lmer, 3) #nothing, we can still say 6_1 and 6_2 are undergoing more plasticity, but not due to sex change?
# I wonder if the plasticity is due to parental care assay

module_corr_plot <- sub_6@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(mean_ieg = mean(iegs21),
            mean_behaviors_day_2 = mean(Behaviors_Day_2),
            mean_time = mean(Time_Day_2))

ggplot(
  subset(module_corr_plot, sub.cluster =='6_1' &  Status %in% c('M','D','F')),
  aes(x = mean_behaviors_day_2, y = mean_ieg, color = Status, shape = Status))+
    geom_point()+
  geom_smooth(method = 'lm')

ggplot(
  subset(module_corr_plot, sub.cluster =='6_2' &  Status %in% c('M','D','F')),
  aes(x = mean_behaviors_day_2, y = mean_ieg, color = Status, shape = Status))+
    geom_point()+
  geom_smooth(method = 'lm')
# seems like nothing to me

ggplot(
  subset(module_corr_plot, sub.cluster =='6_1' &  Status %in% c('M','D','F')),
  aes(x = mean_time, y = mean_ieg, color = Status, shape = Status))+
    geom_point()+
  geom_smooth(method = 'lm')

ggplot(
  subset(module_corr_plot, sub.cluster =='6_2' &  Status %in% c('M','D','F')),
  aes(x = mean_time, y = mean_ieg, color = Status, shape = Status))+
    geom_point()+
  geom_smooth(method = 'lm')
# could be something

time_ieg = lm(mean_ieg~mean_time+Status, data =  subset(module_corr_plot, sub.cluster =='6_2' &  Status %in% c('M','D','F')))
anova(time_ieg, test ='Chisq')
# nothing

#then why do I see this IEG signal? It seems unrelated to sex or behavior
# maybe its not that 6_1 and 6_2 are special, but that 6_0 and 6_3 are?


DotPlot(sub_6, 
        c('ar',
          'esr2a',
          'esr2b',
          'esrrga',
          'gal',
          'npy',
          'kiss1ra',
          'cckbra',
          'cckbrb',
          'LOC111571064'),
        group.by = 'sub.cluster')+
  coord_flip()

##cckbrba and gnrh1 are in the same cells


sub_6$cckbrb_gnrh3 = ifelse(
  (sub_6@assays$RNA$data['cckbrb',] > 0 & 
    ! sub_6@assays$RNA$data['LOC111571064',] > 0  ),
    'cckbrb_only',
  NA
)

sub_6$cckbrb_gnrh3 = ifelse(
  (!sub_6@assays$RNA$data['cckbrb',] > 0 & 
     sub_6@assays$RNA$data['LOC111571064',] > 0  ),
    'gnrh1_only',
  sub_6$cckbrb_gnrh3 
)

sub_6$cckbrb_gnrh3 = ifelse(
  (sub_6@assays$RNA$data['cckbrb',] > 0 & 
     sub_6@assays$RNA$data['LOC111571064',] > 0  ),
    'both',
  sub_6$cckbrb_gnrh3 
)

table(sub_6$sub.cluster, sub_6$cckbrb_gnrh3)

#what happens if I subcluster 6_3 farther

sub_6 = FindSubCluster(sub_6, '6_3', subcluster.name = 'sub_sub.cluster', graph.name = 'harmony.wsnn')
DimPlot(sub_6, group.by ='sub_sub.cluster')
# I guess its not a fair thing to do

marks_6_3 = FindMarkers(sub_6, '6_3')

#maybe there is a difference in 6_3 excitability

ieg_module_lmer_3 = lmer(iegs21~Status+(1|individual), data = subset(sub_6@meta.data, sub.cluster =='6_3' & 
                                                                     Status %in% c('M','D','F')))
car::Anova(ieg_module_lmer_3, 3)#nein

# do the cck or gnrh neurons differ in excitability
gonadotroph_ieg = sub_6@meta.data%>%
  group_by(sub.cluster, cckbrb_gnrh3, individual, Status)%>%
  summarize(ieg = mean(iegs21),
            se_ieg = sd(iegs21)/sqrt(n()))%>%
  subset(!is.na(cckbrb_gnrh3))

ggplot(subset(gonadotroph_ieg, sub.cluster == '6_3' & Status %in% c('M','D','F') & cckbrb_gnrh3 != 'both'), aes(x = Status, y = ieg, color = cckbrb_gnrh3))+
  geom_boxplot(alpha = 0.2)+
  geom_pointrange(aes(x =Status, y = ieg, ymin = ieg - se_ieg, ymax  =ieg+se_ieg ), position = position_jitterdodge(.75))

cck_gnrh_model_6_3 = lmer(iegs21~cckbrb_gnrh3*Status +(1|individual), data = subset(sub_6@meta.data, sub.cluster =='6_3' & 
                                                                     Status %in% c('M','D','F')&
                                                                       cckbrb_gnrh3 %in% c('cckbrb_only', 'gnrh1_only')))
car::Anova(cck_gnrh_model_6_3, 3)
pairs(emmeans(cck_gnrh_model_6_3, c('cckbrb_gnrh3', 'Status')), adjust ='none')

library(sjPlot)
plot_model(cck_gnrh_model_6_3,type = 'eff', terms = c('cckbrb_gnrh3', 'Status') )

cck_gnrh_dat =as.data.frame.array(table(sub_6@meta.data$cckbrb_gnrh3, sub_6@meta.data$sub.cluster, sub_6@meta.data$Status))
# might be glmer time, do they differ in prop my conclusion is yes

#fishers exact test
fish_cck = data.frame()
for(i in paste0('6_', 0:3)){
  
  subset_dat = subset(sub_6@meta.data, sub.cluster == i & Status %in% c('M',"D","F")) 
  subset_dat$Status = as.character(subset_dat$Status)
  
  cck_table = as.data.frame.matrix(table(subset_dat$Status, subset_dat$cckbrb_gnrh3 ))
  cck_matrix = cbind(cck_table[,'cckbrb_only'], rowSums(cck_table)-cck_table[,'cckbrb_only'])  
  d_m_cck_test = fisher.test(cck_matrix[c(1,3),])
  d_m_cck_p = d_m_cck_test$p.value
 
  d_f_cck_test = fisher.test(cck_matrix[c(1,2),])
  d_f_cck_p = d_f_cck_test$p.value
  
    f_m_cck_test = fisher.test(cck_matrix[c(2,3),])
  f_m_cck_p = f_m_cck_test$p.value
  
  newd = data.frame(cluster =i,
                    d_m_cck_p= d_m_cck_p,
                    d_f_cck_p= d_f_cck_p,
                   f_m_cck_p = f_m_cck_p )
  
  fish_cck = rbind(fish_cck, newd)
}
# nuffink

prop_cck = sub_6@meta.data%>%
  group_by(individual,Status, sub.cluster)%>%
  summarise(n_cck = sum(cckbrb_gnrh3 == 'cckbrb_only', na.rm = T),
            n_cells = n())%>%
  mutate(prop_cck = n_cck/n_cells)

ggplot(subset(prop_cck,sub.cluster =='6_3'), aes(x = Status, y = prop_cck))+
  geom_boxplot()+
  geom_point()

fish_gnrh = data.frame()
for(i in paste0('6_', 0:3)){
  
  subset_dat = subset(sub_6@meta.data, sub.cluster == i & Status %in% c('M',"D","F")) 
  subset_dat$Status = as.character(subset_dat$Status)
  
  gnrh_table = as.data.frame.matrix(table(subset_dat$Status, subset_dat$cckbrb_gnrh3 ))
  gnrh_matrix = cbind(gnrh_table[,'gnrh1_only'], rowSums(gnrh_table)-gnrh_table[,'gnrh1_only'])
  
  d_m_gnrh_test = fisher.test(gnrh_matrix[c(1,3),])
  d_m_cck_p = d_m_gnrh_test$p.value
 
  d_f_gnrh_test = fisher.test(gnrh_matrix[c(1,2),])
  d_f_cck_p = d_f_gnrh_test$p.value
  
    f_m_gnrh_test = fisher.test(gnrh_matrix[c(2,3),])
  f_m_cck_p = f_m_gnrh_test$p.value
  
  newd = data.frame(cluster =i,
                    d_m_cck_p= d_m_cck_p,
                    d_f_cck_p= d_f_cck_p,
                   f_m_cck_p = f_m_cck_p )
  
  fish_gnrh = rbind(fish_gnrh, newd)
}

#kinda surprised its othing, could argue df difference in 6_#

prop_gnrh = sub_6@meta.data%>%
  group_by(individual,Status, sub.cluster)%>%
  summarise(n_gnrh = sum(cckbrb_gnrh3 == 'gnrh1_only', na.rm = T),
            n_cells = n())%>%
  mutate(prop_gnrh = n_gnrh/n_cells)

ggplot(subset(prop_gnrh,sub.cluster =='6_3'), aes(x = Status, y = prop_gnrh))+
  geom_boxplot()+
  geom_point()


marks_cck <- FindMarkers(subset(sub_6, sub.cluster =='6_3'), ident.1 = 'gnrh1_only', ident.2 = 'cckbrb_only', group.by = 'cckbrb_gnrh3')

### is there a way I can find sets of neurons that are selectively active or inactive in males and females
FeaturePlot(subset(sub_6, sub.cluster == '6_3'), 'iegs21')

hist(sub_6$iegs21[sub_6$sub.cluster=='6_3'])
bounds = quantile(sub_6$iegs21[sub_6$sub.cluster=='6_3'], probs = c(0.05, 0.95))

sub_6$activity = ifelse(sub_6$iegs21<= bounds[1], 'inhibited', '95th')
sub_6$activity = ifelse(sub_6$iegs21>= bounds[2], 'stimulated', sub_6$activity)

activity_marks = FindMarkers(sub_6, group.by = 'activity', ident.1 = 'stimulated', test.use ='LR')
#?dclk1a

clown_go(rownames(activity_marks[activity_marks$p_val_adj<0.05 &
                                   activity_marks$pct.1>activity_marks$pct.2,]))%>%
  dotplot()

#avp, oxt enriched in stimulated cells, ar enriched in inhibited cells
#hsd17b12a - e1 to e2 conversion enriched in stimulated cells

terms = clown_go(rownames(activity_marks[activity_marks$p_val_adj<0.05 &
                                   activity_marks$pct.1>activity_marks$pct.2,]))@result
terms$geneID[terms$Description=='integrated stress response signaling']
terms$geneID[terms$Description=='response to steroid hormone']

clown_go(rownames(activity_marks[activity_marks$p_val_adj<0.05 &
                                   activity_marks$pct.1<activity_marks$pct.2,]))%>%
  dotplot()


cyto = CytoTRACE(sub_6@assays$RNA$data%>%as.matrix())
sub_6$cyto = cyto$CytoTRACE

ggplot(subset(sub_6@meta.data, activity != '95th' &
                Status %in% c('M','D','F')&
                sub.cluster=='6_3'), aes(x = activity, y = cyto))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~Status)
# in males, inhibited neurons seem to be more immature in 6_0

ggplot(subset(sub_6@meta.data, 
                Status %in% c('M','D','F')&
                sub.cluster=='6_3'), aes(x = activity, y = cyto))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~Status)

#lets do an unbiased search of IEG+ genes
# First, get the cells from cluster 6_3 with Status M/D/F
cells_to_use = sub_6@meta.data$sub.cluster == '6_3' & 
               sub_6@meta.data$Status %in% c('M','D','F')

data_mat = sub_6@assays$RNA$data[, cells_to_use] %>% as.matrix() %>% t()
data_mat = data_mat > 0
genes = colnames(data_mat[, colSums(data_mat) > 0])

IEG_test = function(gene){
  meta_copy = sub_6@meta.data
  
  if(gene %in% genes){
    
    message(which(genes == gene))
    
    # Initialize positive column with NA for all cells
    meta_copy$positive = NA
    
    # Get gene data for all cells in sub_6 (not just the subset)
    gene_expr = sub_6@assays$RNA$data[gene, ] > 0
    meta_copy$positive = ifelse(gene_expr, 'positive', 'negative')
    
    # Filter to only the cells we want for the model
    meta_copy = meta_copy[cells_to_use, ]
    
    # Wrap model fitting and analysis in tryCatch
    result = tryCatch({
      model = lmer(iegs21~positive*Status +(1|individual), data = meta_copy)
      singular = isSingular(model)
      
      av = car::Anova(model, 3)
      positive_p = av$`Pr(>Chisq)`[2]
      av_p = av$`Pr(>Chisq)`[3]
      interaction_p = av$`Pr(>Chisq)`[4]
      
      # Main effect pairwise comparisons
      pair = pairs(emmeans(model, 'Status'), adjust = 'none') %>% as.data.frame()
      
      d_m_p.value = pair$p.value[pair$contrast== 'M - D']
      d_f_p.value = pair$p.value[pair$contrast== 'D - F']
      f_m_p.value = pair$p.value[pair$contrast== 'M - F']
      d_m_estimate = pair$estimate[pair$contrast== 'M - D']
      d_f_estimate = pair$estimate[pair$contrast== 'D - F']
      f_m_estimate = pair$estimate[pair$contrast== 'M - F']
      
      # Interaction pairwise comparisons
      interaction_pairs = pairs(emmeans(model, ~ positive|Status), adjust = 'none') %>% as.data.frame()
      
      # Extract interaction comparisons for each Status level
      int_d_p.value = interaction_pairs$p.value[interaction_pairs$Status == 'D']
      int_f_p.value = interaction_pairs$p.value[interaction_pairs$Status == 'F']
      int_m_p.value = interaction_pairs$p.value[interaction_pairs$Status == 'M']
      int_d_estimate = interaction_pairs$estimate[interaction_pairs$Status == 'D']
      int_f_estimate = interaction_pairs$estimate[interaction_pairs$Status == 'F']
      int_m_estimate = interaction_pairs$estimate[interaction_pairs$Status == 'M']
      
      newd = data.frame(
        gene = gene,
        av_p = av_p,
        positive_p = positive_p,
        interaction_p = interaction_p, 
        d_m_p.value = d_m_p.value,
        d_f_p.value = d_f_p.value, 
        f_m_p.value = f_m_p.value,
        d_m_estimate = d_m_estimate,
        d_f_estimate = d_f_estimate,
        f_m_estimate = f_m_estimate,
        int_d_p.value = int_d_p.value,
        int_f_p.value = int_f_p.value,
        int_m_p.value = int_m_p.value,
        int_d_estimate = int_d_estimate,
        int_f_estimate = int_f_estimate,
        int_m_estimate = int_m_estimate,
        warning = ifelse(length(model@optinfo$conv$lme4$code) != 0, 
                        substr(model@optinfo$conv$lme4$messages, 1, 50), 
                        NA),
        error = NA
      )
      
      return(newd)
      
    }, error = function(e){
      # Return a data frame with NA values if there's an error
      newd = data.frame(
        gene = gene,
        av_p = NA,
        positive_p = NA,
        interaction_p = NA, 
        d_m_p.value = NA,
        d_f_p.value = NA, 
        f_m_p.value = NA,
        d_m_estimate = NA,
        d_f_estimate = NA,
        f_m_estimate = NA,
        int_d_p.value = NA,
        int_f_p.value = NA,
        int_m_p.value = NA,
        int_d_estimate = NA,
        int_f_estimate = NA,
        int_m_estimate = NA,
        warning = NA,
        error = substr(e$message, 1, 100)
      )
      
      message(paste("Error for gene", gene, ":", e$message))
      return(newd)
    })
    
    return(result)
    
  } else {
    message('Gene not in genes')
    return(NULL)
  }
}

ieg_degs = lapply(genes, IEG_test)

ieg_degs_df = do.call(rbind, ieg_degs)
ieg_degs_df = ieg_degs_df[complete.cases(ieg_degs_df[1:(ncol(ieg_degs_df)-2)]),]

ieg_degs_df$av_q = p.adjust(ieg_degs_df$av_p, 'fdr', nrow(ieg_degs_df))
ieg_degs_df$interaction_q = p.adjust(ieg_degs_df$interaction_p, 'fdr', nrow(ieg_degs_df))
ieg_degs_df$positive_q = p.adjust(ieg_degs_df$positive_p, 'fdr', nrow(ieg_degs_df))

interesting_genes_1 = ieg_degs_df[ieg_degs_df$av_q<0.05 
             # & ieg_degs_df$interaction_q<0.05 
              &ieg_degs_df$positive_q<0.05
            ,]

clown_go(interesting_genes_1$gene)%>%dotplot() 


interesting_genes_2 = ieg_degs_df[
             ieg_degs_df$interaction_q<0.05 
            ,]


interesting_genes_3 = ieg_degs_df[
             ieg_degs_df$positive_q<0.05
            ,]

clown_go(interesting_genes_3$gene)%>%dotplot() 

mecp(sub_6,'gabrg2', '6_3', 'sub.cluster' )
prop_cluster_plot(sub_6,'gabrg2', '6_3', 'sub.cluster' )
sub_6$gabrg2 = ifelse(sub_6@assays$RNA$data['gabrg2',]>0, T,F)

sub_6_gabrg2_plot  =sub_6@meta.data%>%
  group_by(individual, Status, gabrg2, sub.cluster)%>%
  summarize(mean_ieg = mean(iegs21))%>%
  subset(Status %in% c('M','D','F'))

ggplot(sub_6_gabrg2_plot, aes(x = Status, y = mean_ieg, color = gabrg2))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(.5))+
  facet_wrap(~sub.cluster, scales ='free')

plot_gene_ieg = function(gene_name, seurat_obj = sub_6, statuses = c('M','D','F')) {
  
  # Check if gene exists in the data
  if (!gene_name %in% rownames(seurat_obj@assays$RNA$data)) {
    stop(paste("Gene", gene_name, "not found in the data"))
  }
  
  # Create binary gene expression variable
  seurat_obj@meta.data[[gene_name]] = ifelse(
    seurat_obj@assays$RNA$data[gene_name,] > 0, 
    TRUE, 
    FALSE
  )
  
  # Prepare plot data
  plot_data = seurat_obj@meta.data %>%
    group_by(individual, Status, !!sym(gene_name), sub.cluster) %>%
    summarize(mean_ieg = mean(iegs21), .groups = 'drop') %>%
    subset(Status %in% statuses)
  
  # Create plot
  p = ggplot(plot_data, aes(x = Status, y = mean_ieg, color = !!sym(gene_name))) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(.5)) +
    facet_wrap(~sub.cluster, scales = 'free') +
    labs(title = paste("IEG expression by", gene_name, "status"),
         y = "Mean IEG21",
         color = gene_name) +
    theme_minimal()
  
  return(p)
}

plot_gene_ieg('cckbrb')
plot_gene_ieg('esr2b')
plot_gene_ieg('npy')
plot_gene_ieg('tacr3a')
plot_gene_ieg('kiss1ra')
plot_gene_ieg('npy8ar')
plot_gene_ieg('tac1')
plot_gene_ieg('gal')
plot_gene_ieg('tac1')

plot_gene_ieg('eif3ea')
plot_gene_ieg('nav3')

#write.csv(ieg_degs_df, 'DEG Outputs/6_3_ieg_genes.csv')

clown_go(ieg_degs_df$gene[ieg_degs_df$interaction_q<0.1])%>%dotplot()
clown_go(ieg_degs_df$gene[ieg_degs_df$positive_q<0.1])%>%dotplot()

plot_gene_ieg('fbxo25')

#### unbiased ----
#cells_to_use = sub_6@meta.data$sub.cluster == '6_0' & 
            #   sub_6@meta.data$Status %in% c('M','D','F')

#data_mat = sub_6@assays$RNA$data[, cells_to_use] %>% as.matrix() %>% t()
#data_mat = data_mat > 0
#genes = colnames(data_mat[, colSums(data_mat) > 0])

IEG_test = function(gene){
  meta_copy = sub_6@meta.data
  
  if(gene %in% genes){
    
    message(which(genes == gene))
    
    # Initialize positive column with NA for all cells
    meta_copy$positive = NA
    
    # Get gene data for all cells in sub_6 (not just the subset)
    gene_expr = sub_6@assays$RNA$data[gene, ] > 0
    meta_copy$positive = ifelse(gene_expr, 'positive', 'negative')
    
    # Filter to only the cells we want for the model
    meta_copy = meta_copy[cells_to_use, ]
    
    # Wrap model fitting and analysis in tryCatch
    result = tryCatch({
      model = lmer(iegs21~positive*Status +(1|individual), data = meta_copy)
      singular = isSingular(model)
      
      av = car::Anova(model, 3)
      positive_p = av$`Pr(>Chisq)`[2]
      av_p = av$`Pr(>Chisq)`[3]
      interaction_p = av$`Pr(>Chisq)`[4]
      
      # Main effect pairwise comparisons
      pair = pairs(emmeans(model, 'Status'), adjust = 'none') %>% as.data.frame()
      
      d_m_p.value = pair$p.value[pair$contrast== 'M - D']
      d_f_p.value = pair$p.value[pair$contrast== 'D - F']
      f_m_p.value = pair$p.value[pair$contrast== 'M - F']
      d_m_estimate = pair$estimate[pair$contrast== 'M - D']
      d_f_estimate = pair$estimate[pair$contrast== 'D - F']
      f_m_estimate = pair$estimate[pair$contrast== 'M - F']
      
      # Interaction pairwise comparisons
      interaction_pairs = pairs(emmeans(model, ~ positive|Status), adjust = 'none') %>% as.data.frame()
      
      # Extract interaction comparisons for each Status level
      int_d_p.value = interaction_pairs$p.value[interaction_pairs$Status == 'D']
      int_f_p.value = interaction_pairs$p.value[interaction_pairs$Status == 'F']
      int_m_p.value = interaction_pairs$p.value[interaction_pairs$Status == 'M']
      int_d_estimate = interaction_pairs$estimate[interaction_pairs$Status == 'D']
      int_f_estimate = interaction_pairs$estimate[interaction_pairs$Status == 'F']
      int_m_estimate = interaction_pairs$estimate[interaction_pairs$Status == 'M']
      
      newd = data.frame(
        gene = gene,
        av_p = av_p,
        positive_p = positive_p,
        interaction_p = interaction_p, 
        d_m_p.value = d_m_p.value,
        d_f_p.value = d_f_p.value, 
        f_m_p.value = f_m_p.value,
        d_m_estimate = d_m_estimate,
        d_f_estimate = d_f_estimate,
        f_m_estimate = f_m_estimate,
        int_d_p.value = int_d_p.value,
        int_f_p.value = int_f_p.value,
        int_m_p.value = int_m_p.value,
        int_d_estimate = int_d_estimate,
        int_f_estimate = int_f_estimate,
        int_m_estimate = int_m_estimate,
        warning = ifelse(length(model@optinfo$conv$lme4$code) != 0, 
                        substr(model@optinfo$conv$lme4$messages, 1, 50), 
                        NA),
        error = NA
      )
      
      return(newd)
      
    }, error = function(e){
      # Return a data frame with NA values if there's an error
      newd = data.frame(
        gene = gene,
        av_p = NA,
        positive_p = NA,
        interaction_p = NA, 
        d_m_p.value = NA,
        d_f_p.value = NA, 
        f_m_p.value = NA,
        d_m_estimate = NA,
        d_f_estimate = NA,
        f_m_estimate = NA,
        int_d_p.value = NA,
        int_f_p.value = NA,
        int_m_p.value = NA,
        int_d_estimate = NA,
        int_f_estimate = NA,
        int_m_estimate = NA,
        warning = NA,
        error = substr(e$message, 1, 100)
      )
      
      message(paste("Error for gene", gene, ":", e$message))
      return(newd)
    })
    
    return(result)
    
  } else {
    message('Gene not in genes')
    return(NULL)
  }
}

#ieg_degs_6_0 = lapply(genes, IEG_test)

#ieg_degs_6_0_df = do.call(rbind, ieg_degs_6_0)
#ieg_degs_6_0_df = ieg_degs_6_0_df[complete.cases(ieg_degs_6_0_df[1:(ncol(ieg_degs_6_0_df)-2)]),]

#ieg_degs_6_0_df$av_q = p.adjust(ieg_degs_6_0_df$av_p, 'fdr', nrow(ieg_degs_6_0_df))
#ieg_degs_6_0_df$interaction_q = p.adjust(ieg_degs_6_0_df$interaction_p, 'fdr', nrow(ieg_degs_6_0_df))
#ieg_degs_6_0_df$positive_q = p.adjust(ieg_degs_6_0_df$positive_p, 'fdr', nrow(ieg_degs_6_0_df))

#write.csv(ieg_degs_6_0_df, 'DEG Outputs/6_0_ieg_genes.csv')
plot_gene_ieg('sts')
plot_gene_ieg('oxtrb')
plot_gene_ieg('trh')

plot_gene_ieg('pparaa')
plot_gene_ieg('pparab')

plot_gene_ieg('cckb')
plot_gene_ieg('galr1a')

#loop it

IEG_test = function(gene, min_individuals = 4){
  meta_copy = sub_6@meta.data
  
  if(gene %in% genes){
    
    message(which(genes == gene))
    
    # Initialize positive column with NA for all cells
    meta_copy$positive = NA
    
    # Get gene data for all cells in sub_6 (not just the subset)
    gene_expr = sub_6@assays$RNA$data[gene, ] > 0
    meta_copy$positive = ifelse(gene_expr, 'positive', 'negative')
    
    # Filter to only the cells we want for the model
    meta_copy = meta_copy[cells_to_use, ]
    
    # Check representation: at least min_individuals from each Status in both positive and negative
    representation_check = meta_copy %>%
      group_by(Status, positive) %>%
      summarise(n_individuals = n_distinct(individual), .groups = 'drop')
    
    # Check if each Status has at least min_individuals in both positive and negative
    status_levels = unique(meta_copy$Status)
    sufficient_representation = TRUE
    
    for(status in status_levels){
      n_pos = representation_check$n_individuals[representation_check$Status == status & 
                                                   representation_check$positive == 'positive']
      n_neg = representation_check$n_individuals[representation_check$Status == status & 
                                                   representation_check$positive == 'negative']
      
      # Handle cases where a combination doesn't exist
      if(length(n_pos) == 0) n_pos = 0
      if(length(n_neg) == 0) n_neg = 0
      
      if(n_pos < min_individuals | n_neg < min_individuals){
        sufficient_representation = FALSE
        break
      }
    }
    
    # If insufficient representation, return NA results with informative message
    if(!sufficient_representation){
      message(paste("Gene", gene, "has insufficient representation across Status groups"))
      
      newd = data.frame(
        gene = gene,
        av_p = NA,
        positive_p = NA,
        interaction_p = NA, 
        d_m_p.value = NA,
        d_f_p.value = NA, 
        f_m_p.value = NA,
        d_m_estimate = NA,
        d_f_estimate = NA,
        f_m_estimate = NA,
        int_d_p.value = NA,
        int_f_p.value = NA,
        int_m_p.value = NA,
        int_d_estimate = NA,
        int_f_estimate = NA,
        int_m_estimate = NA,
        singular = NA,  # Added this
        warning = "Insufficient representation",
        error = NA
      )
      
      return(newd)
    }
    
    # Wrap model fitting and analysis in tryCatch
    result = tryCatch({
      model = lmer(iegs21~positive*Status +(1|individual), data = meta_copy)
      singular = isSingular(model)
      
      av = car::Anova(model, 3)
      positive_p = av$`Pr(>Chisq)`[2]
      av_p = av$`Pr(>Chisq)`[3]
      interaction_p = av$`Pr(>Chisq)`[4]
      
      # Main effect pairwise comparisons
      pair = pairs(emmeans(model, 'Status'), adjust = 'none') %>% as.data.frame()
      
      d_m_p.value = pair$p.value[pair$contrast== 'M - D']
      d_f_p.value = pair$p.value[pair$contrast== 'D - F']
      f_m_p.value = pair$p.value[pair$contrast== 'M - F']
      d_m_estimate = pair$estimate[pair$contrast== 'M - D']
      d_f_estimate = pair$estimate[pair$contrast== 'D - F']
      f_m_estimate = pair$estimate[pair$contrast== 'M - F']
      
      # Interaction pairwise comparisons
      interaction_pairs = pairs(emmeans(model, ~ positive|Status), adjust = 'none') %>% as.data.frame()
      
      # Extract interaction comparisons for each Status level
      int_d_p.value = interaction_pairs$p.value[interaction_pairs$Status == 'D']
      int_f_p.value = interaction_pairs$p.value[interaction_pairs$Status == 'F']
      int_m_p.value = interaction_pairs$p.value[interaction_pairs$Status == 'M']
      int_d_estimate = interaction_pairs$estimate[interaction_pairs$Status == 'D']
      int_f_estimate = interaction_pairs$estimate[interaction_pairs$Status == 'F']
      int_m_estimate = interaction_pairs$estimate[interaction_pairs$Status == 'M']
      
      newd = data.frame(
        gene = gene,
        av_p = av_p,
        positive_p = positive_p,
        interaction_p = interaction_p, 
        d_m_p.value = d_m_p.value,
        d_f_p.value = d_f_p.value, 
        f_m_p.value = f_m_p.value,
        d_m_estimate = d_m_estimate,
        d_f_estimate = d_f_estimate,
        f_m_estimate = f_m_estimate,
        int_d_p.value = int_d_p.value,
        int_f_p.value = int_f_p.value,
        int_m_p.value = int_m_p.value,
        int_d_estimate = int_d_estimate,
        int_f_estimate = int_f_estimate,
        int_m_estimate = int_m_estimate,
        singular = singular,  # Using the variable calculated above
        warning = ifelse(length(model@optinfo$conv$lme4$code) != 0, 
                        substr(model@optinfo$conv$lme4$messages, 1, 50), 
                        NA),
        error = NA
      )
      
      return(newd)
      
    }, error = function(e){
      # Return a data frame with NA values if there's an error
      newd = data.frame(
        gene = gene,
        av_p = NA,
        positive_p = NA,
        interaction_p = NA, 
        d_m_p.value = NA,
        d_f_p.value = NA, 
        f_m_p.value = NA,
        d_m_estimate = NA,
        d_f_estimate = NA,
        f_m_estimate = NA,
        int_d_p.value = NA,
        int_f_p.value = NA,
        int_m_p.value = NA,
        int_d_estimate = NA,
        int_f_estimate = NA,
        int_m_estimate = NA,
        singular = NA,  # Added this
        warning = NA,
        error = substr(e$message, 1, 100)
      )
      
      message(paste("Error for gene", gene, ":", e$message))
      return(newd)
    })
    
    return(result)
    
  } else {
    message('Gene not in genes')
    return(NULL)
  }
}

for(i in paste0('6_', c(
  2,
  1,
  3,
  0
  ))){
  cells_to_use = sub_6@meta.data$sub.cluster == i & 
                 sub_6@meta.data$Status %in% c('M','D','F')
  data_mat = sub_6@assays$RNA$data[, cells_to_use] %>% as.matrix() %>% t()
  data_mat = data_mat > 0
  genes = colnames(data_mat[, colSums(data_mat) > 0])
  
  ieg_degs_list = lapply(genes, IEG_test)
  
  # Remove NULL entries before rbind
  ieg_degs_list = ieg_degs_list[!sapply(ieg_degs_list, is.null)]
  
  # Check if there are any results left
  if(length(ieg_degs_list) == 0){
    message(paste("No valid results for cluster", i))
    next
  }
  
  ieg_degs_df = do.call(rbind, ieg_degs_list)
  ieg_degs_df = ieg_degs_df[complete.cases(ieg_degs_df[1:(ncol(ieg_degs_df)-2)]),]
  ieg_degs_df = subset(ieg_degs_df, is.na(warning))
  ieg_degs_df$av_q = p.adjust(ieg_degs_df$av_p, 'fdr', nrow(ieg_degs_df))
  ieg_degs_df$interaction_q = p.adjust(ieg_degs_df$interaction_p, 'fdr', nrow(ieg_degs_df))
  ieg_degs_df$positive_q = p.adjust(ieg_degs_df$positive_p, 'fdr', nrow(ieg_degs_df))
  
  write.csv(ieg_degs_df, 
            paste0('DEG Outputs/6 subcluster iegs/', i, '_ieg_genes.csv'),
            row.names = FALSE)
  
  message(paste("Completed cluster", i, "with", nrow(ieg_degs_df), "genes"))
}

plot_gene_ieg('nav3')

ieg_interactions_6_0 = read.csv("DEG Outputs/6 subcluster iegs/6_0_ieg_genes.csv")
clown_go(ieg_interactions_6_0$gene[ieg_interactions_6_0$interaction_q<0.05])%>%dotplot()

ieg_interactions_6_1 = read.csv("DEG Outputs/6 subcluster iegs/6_1_ieg_genes.csv")
#clown_go(ieg_interactions_6_1$gene[ieg_interactions_6_1$interaction_q<0.05])%>%dotplot()

ieg_interactions_6_2 = read.csv("DEG Outputs/6 subcluster iegs/6_2_ieg_genes.csv")
#clown_go(ieg_interactions_6_2$gene[ieg_interactions_6_2$interaction_q<0.5])%>%dotplot()

ieg_interactions_6_3 = read.csv("DEG Outputs/6 subcluster iegs/6_3_ieg_genes.csv")
clown_go(ieg_interactions_6_3$gene[ieg_interactions_6_3$interaction_q<0.05])%>%dotplot()
clown_go(ieg_interactions_6_3$gene[ieg_interactions_6_3$positive_q<0.05])%>%dotplot()

ieg_interactions_6_3$gene[ieg_interactions_6_3$interaction_q<0.1]
ieg_interactions_6_3$gene[ieg_interactions_6_3$positive_q<0.05]
plot_gene_ieg('arxa')
plot_gene_ieg('vwc2l')


### i wonder if IEG module is the wrong way to do it, and if IEG + vs - might be more appropriate


