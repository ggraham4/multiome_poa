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

  
  prop_cluster_plot <- function(object,
                         gene,
                         cluster,
                         clustering = 'harmony.wnn_res0.4_clusters'
){
  options(dplyr.summarise.inform = FALSE)


  counts <- suppressWarnings(t(as.matrix(object@assays$RNA$data[, object@meta.data[[clustering]] == cluster ]))[,gene])
  
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
#  clown_go <- readRDS('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/Gabe/clown_go.rds')

  
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
    return(results)
  }
  
  
  
}

obj <- readRDS('C:/Users/Gabe/Desktop/RNA Object.rds')

prop_deg_function <- function(object = obj, cluster = 19, clustering = 'harmony.wnn_res0.4_clusters'){
  start <- Sys.time()
  library(lme4)
  library(dplyr)
  library(parallel)
  library(car)
  
  options(dplyr.summarise.inform = FALSE)
  
  #extract counts matrix 
  counts <- t(as.matrix(object@assays$RNA$data[, object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")]))
  
  binary_counts <- (counts>0)+0
  
  n_cells <- nrow(binary_counts)
  n_genes <- ncol(binary_counts)
  
  #convert to matrix
  counts_matrix <- matrix(data = binary_counts, 
                          nrow =n_cells,
                          ncol= n_genes)
  
  rownames(counts_matrix) <- rownames(binary_counts)
  colnames(counts_matrix) <- colnames(binary_counts)
  
  ##Remove 0 genes 
  filtered_cols_matrix <- counts_matrix[,which(colSums(counts_matrix) != 0)]
  
  #make meta data column
  meta_data <- data.frame(
    cells = rownames(object@meta.data[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D"),]),
    Status = object@meta.data$Status[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")],
    individual =object@meta.data$individual[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")]
  )
  
  
  n_cells_individual <- meta_data%>%
    group_by(individual)%>%
    summarize(n_cells = n())
  
  genes <- colnames(filtered_cols_matrix)
  genes <- genes[!is.na(genes)]
  
  
  deg_function <- function(gene){
    message(paste0('gene ', which(genes == gene), ' of ', n_genes))
    
    gene_expression <- filtered_cols_matrix[, gene, drop = FALSE]
    
    meta_data$gene <- gene_expression
    
    joined_data_restrictions <- meta_data%>%
      group_by(Status)%>%
      summarize(cells_expressing = sum(gene))
    
    if (any(joined_data_restrictions$cells_expressing < 1)) return(NULL)
    
    data_for_matrix <- meta_data %>%
      group_by(individual, Status) %>%
      summarise(
        cells_expressing = sum(gene),
        cells_total = n()
      )
    
    model_matrix <- with(data_for_matrix, cbind(cells_expressing, cells_total - cells_expressing))
    
    
    if (length(unique(model_matrix[, 1])) == 1) return(NULL)
    
    
    glmer_model <- tryCatch({suppressMessages(glmer(model_matrix~Status + (1|individual), family = binomial('logit'), data = data_for_matrix))
    }, error = function(e) {
    return(NULL)
  })
    if (is.null(glmer_model)) return(NULL)
    
    pairwise_comps <- as.data.frame(pairs(emmeans(glmer_model, 'Status'), adjust ='none'))
    
    data_for_output <- data_for_matrix %>%
      group_by(Status) %>%
      summarise(
        prop_expressing = sum(cells_expressing) / sum(cells_total),
        cells_expressing = sum(cells_expressing),
        cells_total = sum(cells_total)
      ) %>%
      pivot_wider(names_from = 'Status', values_from = c(cells_expressing,cells_total, prop_expressing))
    
    data_for_output$anova_p.value <- car::Anova(glmer_model, type = 'III')[2,3]
    data_for_output$d_f_p.value <- pairwise_comps$p.value[pairwise_comps$contrast=='D - F']
    data_for_output$d_m_p.value <- pairwise_comps$p.value[pairwise_comps$contrast=='D - M']
    data_for_output$f_m_p.value <- pairwise_comps$p.value[pairwise_comps$contrast=='F - M']
    data_for_output$singular = isSingular(glmer_model)
    data_for_output$warning = ifelse(length(glmer_model@optinfo$conv$lme4$code) != 0, substr(glmer_model@optinfo$conv$lme4$messages, 1, 50), NA)
    data_for_output$gene <- gene
    data_for_output$cluster <- cluster
    
    return(data_for_output)
    
  }
  
  library(parallelsugar)
  deg_output <- parallelsugar::mclapply(X=genes, FUN=deg_function
                       ,mc.cores= detectCores()-1
  )
  
  deg_output2 <- do.call(rbind, deg_output)
  
  deg_output2$anova_q.value <- p.adjust(deg_output2$anova_p.value, 'fdr', nrow(deg_output2))
  
  deg_output2$d_f_q.value <- p.adjust(deg_output2$d_f_p.value, 'fdr', nrow(deg_output2))
  deg_output2$d_m_q.value <- p.adjust(deg_output2$d_m_p.value, 'fdr', nrow(deg_output2))
  deg_output2$f_m_q.value <- p.adjust(deg_output2$f_m_p.value, 'fdr', nrow(deg_output2))
  
  end <- Sys.time()
  print(end-start)
  return(deg_output2)
  
}


for (i in c(0:31)) {
  print(i)
  start <-Sys.time()
   prop_deg_data<-  prop_deg_function(object = obj, cluster = i)

  assign(paste0('prop_degs_cluster_', i), prop_deg_data, envir = .GlobalEnv)
  
    write.csv(prop_deg_data, 
              paste0('X:/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/Prop DEGs 011525/prop_degs_cluster_', i, '.csv'))
end <- Sys.time()
print(end-start)
    }

###Analysis ####
for(i in c(6:31)){
  print(i)
  data <- read.csv(paste0('X:/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/Prop DEGs 011525/prop_degs_cluster_',i,'.csv'))
  data <-define_degs(data) 
   assign(paste0('prop_degs_cluster_', i), data, envir = .GlobalEnv)
  }

together_data <- data.frame()
for(i in 0:31){
  data <- get(paste0('prop_degs_cluster_',i))
  together_data <- rbind(together_data, data)
}

together_data_summed <- together_data%>%
    subset(!is.na(class)& is.na(warning))%>%
  group_by(cluster, class)%>%
  summarize(class_count = n())

ggplot(together_data_summed, aes(x = cluster, y = class_count, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DEGs") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
    scale_x_continuous(breaks = c(0:31))+
  scale_y_continuous()+
  scale_fill_manual(values = P40)

go_8_late_down <- clown_go(prop_degs_cluster_8$gene[!is.na(prop_degs_cluster_8$class) & is.na(prop_degs_cluster_8$warning)& prop_degs_cluster_8$class == 'Late Downregulated'])
dotplot(go_8_late_down)

go_8_early_up <- clown_go(prop_degs_cluster_8$gene[!is.na(prop_degs_cluster_8$class) & is.na(prop_degs_cluster_8$warning)& prop_degs_cluster_8$class == 'Early Upregulated'])
dotplot(go_8_early_up)

go_8_trans_up <- clown_go(prop_degs_cluster_8$gene[!is.na(prop_degs_cluster_8$class) & is.na(prop_degs_cluster_8$warning)& prop_degs_cluster_8$class == 'Transiently Upregulated'])
dotplot(go_8_trans_up)

go_4_late_down <- clown_go(prop_degs_cluster_4$gene[!is.na(prop_degs_cluster_4$class) & is.na(prop_degs_cluster_4$warning)& prop_degs_cluster_4$class == 'Late Downregulated'])
dotplot(go_4_late_down)

go_1_early_up<- clown_go(prop_degs_cluster_1$gene[!is.na(prop_degs_cluster_1$class) & is.na(prop_degs_cluster_1$warning)& prop_degs_cluster_1$class == 'Early Upregulated'])
dotplot(go_1_early_up)

go_1_trans_up <- clown_go(prop_degs_cluster_1$gene[!is.na(prop_degs_cluster_1$class) & is.na(prop_degs_cluster_1$warning)& prop_degs_cluster_1$class == 'Transiently Upregulated'])
dotplot(go_1_trans_up)

go_19_all <-  clown_go(prop_degs_cluster_19$gene[!is.na(prop_degs_cluster_19$class) & is.na(prop_degs_cluster_19$warning)])
dotplot(go_19_all)

go_14_all <-  clown_go(prop_degs_cluster_14$gene[!is.na(prop_degs_cluster_14$class) & is.na(prop_degs_cluster_14$warning)])
dotplot(go_14_all)

go_29_all <-  clown_go(prop_degs_cluster_29$gene[!is.na(prop_degs_cluster_29$class) & is.na(prop_degs_cluster_29$warning)])
dotplot(go_29_all)

prop_degs_cluster_8$gene[!is.na(prop_degs_cluster_8$class) & is.na(prop_degs_cluster_8$warning)]

together_data[together_data$gene == 'gfap' & !is.na(together_data$issignif),]

###Plot DEGs
prop_cluster_plot(obj, 
                  'tacr3a',
                  19)

mean_expression_cluster_plot(obj,
                             'tacr3a',
                             19)

mean_expression_cluster_plot(obj,
                             'pgr',
                             19)

tac_plots <- prop_cluster_plot(obj, 
                  'tacr3a',
                  19)+
  mean_expression_cluster_plot(obj,
                             'tacr3a',
                             19)

pgr_plots <- prop_cluster_plot(obj, 
                  'pgr',
                  19)+
  mean_expression_cluster_plot(obj,
                             'pgr',
                             19)


  

