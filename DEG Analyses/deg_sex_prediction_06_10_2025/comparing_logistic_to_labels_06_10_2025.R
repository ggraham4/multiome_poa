library(Seurat)
library(dplyr)
mecp = readRDS("Functions/mean_expression_cluster_plot.rds")
define_degs2 = readRDS("Functions/define_degs2")
### read in data

obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
logistic_data = read.csv("DEG Analyses/deg_sex_prediction_06_10_2025/logistic_probabilities_06_10_2026.csv")

deg_csvs <- dir('DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering')
deg_path <- 'DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/'

# read in DEGs list, I think first I will do all singular and nonsingular
all_degs <- data.frame()
for(i in deg_csvs){
  data <- read.csv(paste0(deg_path,i))
  data = define_degs2(data)
  all_degs <- rbind(all_degs, data)
}

all_degs_subset <- subset(all_degs, av_q.value <0.05 & singular ==F)

joined = all_degs_subset%>%
  right_join(logistic_data, by = c('cluster','gene'))

hist(joined$mean_probability_female[joined$first_word=='Transiently']) # it seems like it really struggles with the transient genes but that's interesting th
hist(joined$mean_probability_female[joined$first_word=='Early'])
hist(joined$mean_probability_female[joined$first_word=='Late']) ## ok so this and early tracks with my classifications

mecp(obj,'LOC111577263',1)

mecp(obj,'tacr3a',6)

table(joined$short_label[joined$first_word=='Transiently'])
