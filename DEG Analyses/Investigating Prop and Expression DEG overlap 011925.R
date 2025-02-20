{
  library(parallel)
  library(clusterProfiler)
  library(blme)
  library(Seurat)
  library(tidyverse)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(SeuratObject)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)
  library(emmeans)
  library(CytoTRACE)
  library(ggrepel)
  
  mean_expression_cluster_plot <- function(object, gene, cluster, clustering = 'harmony.wnn_res0.4_clusters'){
    counts <- t(object@assays$RNA$data[,object@meta.data[[clustering]] == cluster])
  Counts_of_interest <- as.data.frame(counts[,gene])
    Counts_of_interest[[gene]] <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
  results <- data.frame()
  for (i in unique(object@meta.data$individual)) {
    Counts <- Counts_of_interest[[gene]][Counts_of_interest$individual==i]
        mean <- mean(Counts)
        mean_se <- sd(Counts) / sqrt(length(Counts))
        df <- data.frame(
          individual = i,
          mean = mean,
          se = mean_se
        )
        results <- rbind(results, df)
  }
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
    labs(x  ='FishID', y = 'Mean Counts +/- SE of Counts', title = paste0(gene,'_cluster_',cluster))+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
  }
  
  mean_expression_cluster_plot(obj, 
                               'utrn',
                               19)

  
  prop_cluster_plot <- function(object,
                         gene,
                         cluster,
                         clustering = 'harmony.wnn_res0.4_clusters'
){
  options(dplyr.summarise.inform = FALSE)


  counts <- suppressWarnings(t(as.matrix(object@assays$RNA$counts[, object@meta.data[[clustering]] == cluster ]))[,gene])
  
  binary_counts <- (counts>0)+0
  
  n_cells <- nrow(binary_counts)
  n_genes <- ncol(binary_counts)
  
  

  #make meta data column
  full_data <- data.frame(
    cells = rownames(object@meta.data[object@meta.data[[clustering]] == cluster,]),
    Status = object@meta.data$Status[object@meta.data[[clustering]] == cluster],
    individual =object@meta.data$individual[object@meta.data[[clustering]] == cluster]
  )
    full_data$gene_expression <- binary_counts
    
    full_data <- full_data%>%
      group_by(individual, Status)%>%
    summarize(prop = sum(gene_expression)/n(),
              se = sqrt(prop*(1-prop)/n())
              )
  
    full_data$factor <- ifelse(full_data$Status == "NRM", 0, NA)
  full_data$factor <- ifelse(full_data$Status == "M", 1, full_data$factor)
  full_data$factor <- ifelse(full_data$Status == "D", 2, full_data$factor)
  full_data$factor <- ifelse(full_data$Status == "E", 3, full_data$factor)
  full_data$factor <- ifelse(full_data$Status == "NF", 4, full_data$factor)
  full_data$factor <- ifelse(full_data$Status == "F", 5, full_data$factor)

  full_data$Status <- factor(full_data$Status, levels <- c('NRM','M','D','E','NF',"F"))
  
  
  plot <- ggplot(full_data, aes(x = fct_reorder(individual, factor), y = prop, color = Status))+
    geom_boxplot(aes(group = Status, fill = Status), alpha = 0.25, outlier.shape = NA)+
    geom_point()+
    geom_pointrange(aes(x = individual, y = prop, ymin = prop - se, ymax = prop+se))+
    theme_classic()+
    labs(x  ='FishID', y = ' Proportion +/- SE ', title = paste0(gene,'_cluster_',cluster))+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
  
}
  clown_go <- readRDS('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/Gabe/clown_go.rds')

  
  library(Polychrome)
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
names(P40) <- NULL

  
        define_degs <- function(data, singular = TRUE) {
  if (!singular) {
    data <- data[data$singular == FALSE, ]
  }
  
  # Assign classes based on conditions
  data$class <- NA  # Initialize class column
  
  data$class[data$d_m_q.value < 0.05 & 
             data$prop_expressing_D >data$prop_expressing_M  & 
             data$d_f_q.value > 0.05] <- 'Early Upregulated'
  
  
  data$class[data$d_m_q.value > 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$prop_expressing_D <data$prop_expressing_F] <- 'Late Upregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
               data$prop_expressing_D <data$prop_expressing_M & 
             data$d_f_q.value > 0.05] <- 'Early Downregulated'
  
  
  data$class[data$d_m_q.value > 0.05 & 
             data$d_f_q.value < 0.05 & 
               data$prop_expressing_D >data$prop_expressing_F] <- 'Late Downregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
            data$prop_expressing_D >data$prop_expressing_F& 
            data$prop_expressing_D >data$prop_expressing_M] <- 'Transiently Upregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$prop_expressing_D <data$prop_expressing_F & 
             data$prop_expressing_D <data$prop_expressing_M] <- 'Transiently Downregulated'
  
    data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$prop_expressing_D <data$prop_expressing_M & 
             data$prop_expressing_D >data$prop_expressing_F] <- 'Progressively Downregulated'
    
      data$class[data$d_m_q.value < 0.05 & 
                data$d_f_q.value < 0.05 & 
            data$prop_expressing_D >data$prop_expressing_M & 
              data$prop_expressing_D <data$prop_expressing_F] <- 'Progressively Upregulated'
      
      data$class[data$f_m_q.value < 0.05 & 
                data$d_f_q.value > 0.05 & 
                data$d_m_q.value > 0.05 & 
              data$prop_expressing_F >data$prop_expressing_M] <- 'Terminally Upregulated'
      
  data$class[data$f_m_q.value < 0.05 & 
                data$d_f_q.value > 0.05 & 
                data$d_m_q.value > 0.05 & 
             data$prop_expressing_F <data$prop_expressing_M ] <- 'Terminally Downregulated'
  

  data$issignif <- NA
  data$issignif <- ifelse(data$f_m_q.value<0.05|
                            data$d_m_q.value<0.05|
                            data$d_f_q.value<0.05,
                          '*',NA)

  return(data)
        }
  
  mean_expression_cluster_data <- function(object, gene, cluster, clustering = 'harmony.wnn_res0.4_clusters'){
    counts <- t(object@assays$RNA$counts[,object@meta.data[[clustering]] == cluster])
    Counts_of_interest <- as.data.frame(counts[,gene])
    Counts_of_interest[[gene]] <- Counts_of_interest[,1]
    Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
    results <- data.frame()
    
    for (i in unique(object@meta.data$individual)) {
      Counts <- Counts_of_interest[[gene]][Counts_of_interest$individual==i]
      mean <- mean(Counts)
      mean_se <- sd(Counts) / sqrt(length(Counts))
      df <- data.frame(
        individual = i,
        mean = mean,
        se = mean_se
      )
      results <- rbind(results, df)
    }
    results$Sex <- str_sub(results$individual, -1)
    results$Sex[results$individual == 'T17D'] = 'NF'
    results$Sex[results$individual == 'A12D'] = 'E'
    results$Sex[results$individual == 'T11D'] = 'E'
    results$Sex[results$individual == 'GH'] = 'NRM'
    return(results)
  }
  
  
  
}

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

##Read in Prop DEGs ###
for(i in c(0:31)){
  print(i)
  data <- read.csv(paste0('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/Prop DEGs 011525/prop_degs_cluster_',i,'.csv'))
  data <-define_degs(data) 
   assign(paste0('prop_degs_cluster_', i), data, envir = .GlobalEnv)
  }

together_data <- data.frame()
for(i in 1:31){
  print(i)
  data <- get(paste0('prop_degs_cluster_',i))
  if(ncol(data)==25){
  data <- data%>%
    dplyr::select(-X.1)
  }
  together_data <- rbind(together_data, data)
}

#read in expression degs
for(i in c(0:31)){
  print(i)
  data <- read.csv(paste0('//Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/122324 Neg Bin with Doms/cluster_',i,'.csv'))
  data <-define_degs(data) 
   assign(paste0('exp_degs_cluster_', i), data, envir = .GlobalEnv)
  }

together_data_exp <- data.frame()
for(i in 0:31){
  data <- get(paste0('exp_degs_cluster_',i))
  together_data_exp <- rbind(together_data_exp, data)
}

prop_data_genes <- together_data$gene[together_data$issignif=='*' & is.na(together_data$warning) & !is.na(together_data$issignif)]
unique_prop_data_genes <-unique(prop_data_genes)

exp_data_genes <- together_data_exp$gene[together_data_exp$issignif=='*' & is.na(together_data_exp$warning) & !is.na(together_data_exp$issignif)]
unique_exp_data_genes <-unique(exp_data_genes)

all_genes <- unique(append(unique_exp_data_genes,unique_prop_data_genes))

together_genes <- data.frame(genes = all_genes)

together_genes$type <- ifelse(together_genes$genes%in%prop_data_genes & !(together_genes$genes%in%exp_data_genes),
                                           'prop',
                                           NA)

together_genes$type <- ifelse(together_genes$genes%in%exp_data_genes & !(together_genes$genes%in%prop_data_genes),
                                           'exp',
                                           together_genes$type)

together_genes$type <- ifelse(together_genes$genes%in%exp_data_genes & together_genes$genes%in%prop_data_genes,
                                           'both',
                                           together_genes$type)
together_genes_plot <- together_genes%>%
  group_by(type)%>%
  summarize(count = n())

ggplot(together_genes_plot, aes(x = type, y = count))+
  geom_bar(stat = 'identity')

overlapping_genes_go <- clown_go(together_genes$genes[together_genes$type=='both'])
dotplot(overlapping_genes_go)

prop_genes_go <- clown_go(together_genes$genes[together_genes$type=='prop'])
dotplot(prop_genes_go)
prop_genes_go$geneID

exp_genes_go <- clown_go(together_genes$genes[together_genes$type=='exp'])
dotplot(exp_genes_go)

prop_both_genes_go <- clown_go(together_genes$genes[together_genes$type=='prop'|together_genes$type=='both'])
dotplot(prop_both_genes_go)

prop_both_genes_go$geneID[prop_both_genes_go$Description == 'chromatin remodeling']
prop_both_genes_go_enriched <- unlist(strsplit(prop_both_genes_go$geneID[prop_both_genes_go$Count == 25], "/"))

together_data$cluster[together_data$gene %in% prop_both_genes_go_enriched &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()%>% plot()

together_data$cluster[together_data$gene %in% c('kdm3b') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()

prop_cluster_plot(obj,
                  'kdm3b',
                  5)

together_data$cluster[together_data$gene %in% c('kdm5c') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()

prop_cluster_plot(obj,
                  'kdm5c',
                  3)

together_data$cluster[together_data$gene %in% c('kdm4b') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()

prop_cluster_plot(obj,
                  'kdm4b',
                  1)

together_data$cluster[together_data$gene %in% c('kdm6ba') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()
prop_cluster_plot(obj,
                  'kdm6ba',
                  8)

together_data$cluster[together_data$gene %in% c('kmt2ba') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()
prop_cluster_plot(obj,
                  'kmt2ba',
                  8)

together_data$cluster[together_data$gene %in% c('ezh1') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()
prop_cluster_plot(obj,
                  'ezh1',
                  9)

together_data$cluster[together_data$gene %in% c('kmt2bb') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()

prop_cluster_plot(obj,
                  'kmt2bb',
                  5)

together_data$cluster[together_data$gene %in% c('kmt5c') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()
prop_cluster_plot(obj,
                  'kmt5c',
                  13)

together_data$cluster[together_data$gene %in% c('hdac5') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()
prop_cluster_plot(obj,
                  'hdac5',
                  24)

prop_both_genes_go_enriched_signaling <- unlist(strsplit(prop_both_genes_go$geneID[prop_both_genes_go$Description =='chemical synaptic transmission'], "/"))
#lots of discs large homologs here


together_data$cluster[together_data$gene %in% c('dlg3') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()
prop_cluster_plot(obj,
                  'dlg3',
                  9)

together_data$cluster[together_data$gene %in% c('dlg2') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()

prop_cluster_plot(obj,
                  'dlg2',
                  15)

 
together_data$cluster[together_data$gene %in% c('LOC111563200') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()

prop_cluster_plot(obj,
                  'LOC111563200',
                  19)


together_data$cluster[together_data$gene %in% c('LOC111580207') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()
prop_cluster_plot(obj,
                  'LOC111563200',
                  4)

together_data$cluster[together_data$gene %in% c('grm1b') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()
prop_cluster_plot(obj,
                  'grm1b',
                  3)


#### ok what about cell comm, both only ####

overlapping_genes_go_enriched <- unlist(strsplit(overlapping_genes_go$geneID[overlapping_genes_go$Description =='cell communication'], "/"))


together_data$cluster[together_data$gene %in% c('LOC111569772') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()
prop_cluster_plot(obj,
                  'LOC111569772',
                  4)+
  mean_expression_cluster_plot(obj,
                               'LOC111569772',
                               4)

together_data$cluster[together_data$gene %in% c('LOC111564152') &
                        is.na(together_data$warning)&
                        !is.na(together_data$issignif)]%>% table()

prop_cluster_plot(obj,
                  'LOC111564152',
                  8)+
  mean_expression_cluster_plot(obj,
                               'LOC111564152',
                               8)
prop_cluster_plot(obj,
                  'LOC111564152',
                  14)+
  mean_expression_cluster_plot(obj,
                               'LOC111564152',
                               14)






















