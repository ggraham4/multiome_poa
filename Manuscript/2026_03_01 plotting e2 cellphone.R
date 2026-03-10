# CellPhone Plots
library(lme4)
library(tidyverse)
library(Seurat)
library(car)
library(emmeans)
p_value_annotation <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return(paste0("p = ", round(p, 3)))
  }
}


reference = read.csv('/Users/ggraham/Desktop/multiome_poa/Measures/2025_12_26 all_data.csv')
reference$individual=reference$Fish

individuals =reference$Fish[reference$Status %in% c('M','D','E','NF',"F")]
example_data = read_tsv('/Users/ggraham/Desktop/multiome_poa/Python/CellPhoneDb/degs_analysis_means_05_12_2025_222054.txt')
significant_interactions = read_csv("/Users/ggraham/Desktop/multiome_poa/Python/CellPhoneDb/signif_data_bound.csv")

significant_interactions_0.1 = subset(significant_interactions, main_effect_q.value < 0.1)

unique_cell_cell_pairs = colnames(example_data[,14:ncol(example_data)])

all_individual_means <- list()
for(i in individuals){
  message("Loading: ", i)
  filepath <- paste0('/Users/ggraham/Desktop/multiome_poa/Python/CellPhoneDb/Analysis Means/degs_analysis_means_', i, ".txt")
  if(file.exists(filepath)){
    all_individual_means[[i]] <- read_tsv(filepath) %>% suppressMessages()
  }
}

plot_communication_cluster_pair <- function(interacting_pair, cluster_pair, individuals_meta_data, 
                                             significant_interactions_df = significant_interactions_0.1,
                                             preloaded_means = all_individual_means){
  `%notin%` <- Negate(`%in%`)
  
  # --- Build means dataframe from preloaded list ---
  means_together <- data.frame()
  for(i in names(preloaded_means)){
    data <- preloaded_means[[i]]
    if(cluster_pair %notin% colnames(data)) next
    mean_val <- data[data$interacting_pair == interacting_pair, cluster_pair]
    mean_val <- unname(unlist(mean_val))
    if(length(mean_val) == 0) next
    means_together <- rbind(means_together, data.frame(individual = i, mean = mean_val))
  }
  
  # --- Join metadata ---
  full_data <- means_together %>%
    right_join(individuals_meta_data, by = 'individual') %>%
    dplyr::select(individual, mean, Status) %>%
    na.omit()
  
  full_data$Status <- ifelse(full_data$Status == 'D', 'I', full_data$Status)
  full_data$Status <- ifelse(full_data$Status == 'E', 'LI', full_data$Status)
  full_data$Status <- factor(full_data$Status, levels = c('M', 'I', 'LI', 'NF', 'F'))
  
  # --- Pull significance annotations for this interaction + cluster pair ---
  sig_row <- significant_interactions_df[significant_interactions_df$interacting_pair == interacting_pair & 
                    significant_interactions_df$cluster_pair == cluster_pair,]
                   
  
  # Build annotations — adjust column names to match your significant_interactions df
  get_annot <- function(p){
    if(is.null(p) || length(p) == 0 || is.na(p)) return("ns")
    p_value_annotation(p)
  }
  
  d_m_annotation <- get_annot(sig_row$d_m_p.value[1])
  d_f_annotation <- get_annot(sig_row$d_f_p.value[1])
  f_m_annotation <- get_annot(sig_row$f_m_p.value[1])
  
  y_max <- max(full_data$mean, na.rm = TRUE)
  
  # --- Plot ---
  plot <- ggplot(full_data, aes(x = Status, y = mean, fill = Status)) +
    geom_boxplot(data = subset(full_data, Status != "NF"),
                 aes(x = Status, y = mean, fill = Status),
                 outlier.shape = NA) +
    geom_point(size = 0.5) +
    theme_classic() +
    scale_x_discrete(drop = FALSE) +
    labs(x = 'Phase', 
         y = paste0(interacting_pair), 
         title = paste0('Clusters ', cluster_pair)) +
    geom_signif(xmin = c(1),
                xmax = c(1.9),
                y_position = c(y_max * 1.1),
                annotation = c(d_m_annotation),
                color = "black",
                tip_length = c(0, 0), textsize = 3) +
    geom_signif(xmin = c(1),
                xmax = c(5),
                y_position = c(y_max * 1.3),
                annotation = c(f_m_annotation),
                color = "black",
                tip_length = c(0, 0), textsize = 3) +
    geom_signif(xmin = c(2.1),
                xmax = c(5),
                y_position = c(y_max * 1.1),
                annotation = c(d_f_annotation),
                color = "black",
                tip_length = c(0, 0), textsize = 3) +
    theme(legend.position = 'none', 
          axis.title.y = element_blank(),
          axis.title.x = element_blank())+
    scale_fill_manual(values = c('#F8766D','#7cae00','#00bfc4','#c77cff','black'))

  
  return(plot)
}

plot_communication_cluster_pair(interacting_pair = 'Estradiol_byCYP19A1_ESR2', 
                                cluster_pair = '1|8', 
                                individuals_meta_data = reference)

significant_interactions_estradiol = significant_interactions_0.1$cluster_pair[significant_interactions_0.1$interacting_pair=='Estradiol_byCYP19A1_ESR2']

for(i in significant_interactions_estradiol){
  
  p = plot_communication_cluster_pair(interacting_pair = 'Estradiol_byCYP19A1_ESR2', 
                                cluster_pair = i, 
                                individuals_meta_data = reference)
  
      ggsave(plot = p,
       file = paste0('clusters_', i, '_e2.svg'),
       device = "svg",
       units = "in",
       width = 1.35,
       height = 1.45,
       path = "Manuscript/Plots/Fig.4/cellphone_E2")


}

