#Reclustering cluster 19 
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

        define_degs <- function(data, singular = TRUE) {
  if (!singular) {
    data <- data[data$singular == FALSE, ]
  }
  
  # Assign classes based on conditions
  data$class <- NA  # Initialize class column
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_m_estimate > 0 & 
             data$d_f_q.value > 0.05] <- 'Early Upregulated'
  
  
  data$class[data$d_m_q.value > 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$d_f_estimate < 0] <- 'Late Upregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_m_estimate < 0 & 
             data$d_f_q.value > 0.05] <- 'Early Downregulated'
  
  
  data$class[data$d_m_q.value > 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$d_f_estimate > 0] <- 'Late Downregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$d_f_estimate > 0 & 
             data$d_m_estimate > 0] <- 'Transiently Upregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$d_f_estimate < 0 & 
             data$d_m_estimate < 0] <- 'Transiently Downregulated'
  
    data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$d_m_estimate < 0 & 
             data$d_f_estimate > 0] <- 'Progressively Downregulated'
    
      data$class[data$d_m_q.value < 0.05 & 
                data$d_f_q.value < 0.05 & 
             data$d_m_estimate > 0 & 
             data$d_f_estimate < 0] <- 'Progressively Upregulated'
      
      data$class[data$f_m_q.value < 0.05 & 
                data$d_f_q.value > 0.05 & 
                data$d_m_q.value > 0.05 & 
             data$f_m_estimate > 0 ] <- 'Terminally Upregulated'
      
  data$class[data$f_m_q.value < 0.05 & 
                data$d_f_q.value > 0.05 & 
                data$d_m_q.value > 0.05 & 
             data$f_m_estimate < 0  ] <- 'Terminally Downregulated'
  

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

cluster_19 <- subset(obj, harmony.wnn_res0.4_clusters ==19)
cluster_19 <- FindSubCluster(cluster_19, 19, 'harmony.wsnn',  subcluster.name= 'sub')

Idents(cluster_19) <- 'sub'
DimPlot(cluster_19, label = T)

### What are the differences between the subclusters ####
FeaturePlot(cluster_19, c('gad2','slc17a6b'))

DotPlot(object = cluster_19, 
        group.by = "sub", 
        features = c('gad2','LOC111584103','slc17a6b','slc17a7a'))  +
  coord_flip()

markers_19 <- FindAllMarkers(cluster_19)

top_markers <- markers_19 %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 1, with_ties = FALSE) %>%
  ungroup()

DotPlot(object = cluster_19, 
        group.by = "sub", 
        features = top_markers$gene) +
  coord_flip()

clown_go <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/clown_go.rds')

markers_19_0 <- markers_19$gene[markers_19$p_val_adj<0.05& markers_19$cluster=='19_0']
markers_19_1 <- markers_19$gene[markers_19$p_val_adj<0.05& markers_19$cluster=='19_1']
markers_19_2 <- markers_19$gene[markers_19$p_val_adj<0.05& markers_19$cluster=='19_2']

unique_markers_19_0 <- markers_19_0[!(markers_19_0%in% markers_19_1 |markers_19_0 %in%markers_19_2)]
unique_markers_19_1 <- markers_19_1[!(markers_19_1%in% markers_19_0 |markers_19_1 %in%markers_19_2)]
unique_markers_19_2 <- markers_19_2[!(markers_19_2%in% markers_19_0 |markers_19_2 %in%markers_19_1)]
#test if there are any overlaps
unique(unique_markers_19_0 %in%unique_markers_19_1)
unique(unique_markers_19_1 %in%unique_markers_19_2)
unique(unique_markers_19_0 %in%unique_markers_19_2)
#Swag

go_19_0 <- clown_go(unique_markers_19_0)
dotplot(go_19_0)+ labs(title = '19_0')

go_19_1 <- clown_go(unique_markers_19_1)
dotplot(go_19_1)+ labs(title = '19_1')

go_19_2 <- clown_go(unique_markers_19_2)
dotplot(go_19_2)+ labs(title = '19_2')


#### CytoTRACE ######
cluster_19_matrix <- as.matrix(cluster_19@assays$RNA$counts)
cluster_19_cyto <- CytoTRACE(mat = cluster_19_matrix
)

cluster_19$cyto <-cluster_19_cyto$CytoTRACE

FeaturePlot(cluster_19, 'cyto')

cluster_19_data <- data.frame(individual = cluster_19@meta.data$individual,
                              status = cluster_19@meta.data$Status,
                              cluster = cluster_19@meta.data$sub,
                              cyto = cluster_19@meta.data$cyto)

cluster_19_data <- subset(cluster_19_data, status == 'F' |status == 'M' | status == 'D')

cluster_19_cyto_gross <- lmer(cyto~cluster *status+ (1|individual), data = cluster_19_data)

car::Anova(cluster_19_cyto_gross, type = 'III')

cyto_gross_plot <- ggplot(cluster_19_data, aes(x = cluster, y = cyto, group = interaction(cluster, status), color = status))+
  geom_violin(alpha = 0)+
  geom_point(position = position_dodge(0.9))
cyto_gross_plot

library(sjPlot)
plot_model(cluster_19_cyto_gross,terms = c('cluster','status'), type = 'eff')

cyto_gross_pairs <- pairs(emmeans(cluster_19_cyto_gross, c('cluster','status')),adjust = 'none', by = 'cluster')
cyto_gross_pairs

#### Negative binomial ####
neg.bin.mult <- function(obj, clustering = 'sub', cluster, n_cores = detectCores() - 1) {
  start_time <- Sys.time()  # Start timing
  
  message('Extracting Counts')
  counts <- obj@assays$RNA$counts[, obj@meta.data[[clustering]] == cluster & (obj@meta.data$Status == "M" | obj@meta.data$Status == "F" | obj@meta.data$Status == "D")]
  combined_counts <- counts
  
  df_counts <- data.frame(t(combined_counts))
  colnames(df_counts) <- rownames(obj@assays$RNA)
  
  n_genes = ncol(df_counts)
  n_cells = nrow(df_counts)
  
  message("Making Counts Data Frame...")
  df_counts_meta <- data.frame(rownames(df_counts))
  df_counts_meta$id <- df_counts_meta$rownames.df_counts.
  df_counts_meta$rownames.df_counts. = NULL
  df_counts_meta$individual = obj$individual[obj@meta.data[[clustering]] == cluster & (obj@meta.data$Status == "M" | obj@meta.data$Status == "F" | obj@meta.data$Status == "D")]
  df_counts_meta$Status = obj$Status[obj@meta.data[[clustering]] == cluster & (obj@meta.data$Status == "M" | obj@meta.data$Status == "F" | obj@meta.data$Status == "D")]
  
  message("Removing Genes with 0 Counts...")
  df_counts_no_0 <- df_counts[, which(colSums(df_counts) != 0)]
  
  message("Making New Counts Data Frame Without 0s...")
  n_genes_no_0 = ncol(df_counts_no_0)
  
  df_counts_no_0 <- cbind(df_counts_no_0, df_counts_meta)
  df_counts_no_0_split_by_subject <- split(df_counts_no_0, f = df_counts_no_0$individual)
  
  message("Finding Good Genes for Subject...")
  # REMOVE GENES WITH ZERO COUNTS IN EACH SUBJECT 
  for (l in 1:length(df_counts_no_0_split_by_subject)) {
    correct_gene_names <- colnames(df_counts_no_0)
    
    temp_subject_l <- data.frame(df_counts_no_0_split_by_subject[[l]]) ### AND HERE THEY GET FUCKED UP 
    colnames(temp_subject_l) <- correct_gene_names
    
    temp_subject_l_counts <- temp_subject_l[, 1:n_genes_no_0]
    temp_subject_l_counts_no_0 <- temp_subject_l_counts[, which(colSums(temp_subject_l_counts) != 0)]
    out <- data.frame(colnames(temp_subject_l_counts_no_0))
    assign(x = paste0("gene_list_subject_", l), value = get("out"))
  }
  
  # GENERATE A LIST OF GENES FOR EACH SUBJECT
  good_gene_list <- gene_list_subject_1$colnames.temp_subject_l_counts_no_0.
  
  for (m in 2:length(df_counts_no_0_split_by_subject)) {
    temp_good_gene_list_m <- data.frame(value = get(paste0("gene_list_subject_", m)))
    temp_good_gene_list_m <- temp_good_gene_list_m$colnames.temp_subject_l_counts_no_0.
    good_gene_list <- intersect(good_gene_list, temp_good_gene_list_m)
  }
  
  p <- length(good_gene_list)
  
  message('Making Gene Data Frame for Each Subject...')
  valid_genes <- good_gene_list[good_gene_list %in% colnames(df_counts_no_0)]
  v <- length(valid_genes)
  
  df_counts_no_0_all_subjects <- df_counts_no_0[, valid_genes]
  count_matrix_final <- as.matrix(df_counts_no_0_all_subjects)
  count_matrix_final <- as.data.frame(t(count_matrix_final))
  df_counts_no_0_all_subjects <- cbind(df_counts_no_0_all_subjects, df_counts_meta)
  
  Status <- as.factor(df_counts_no_0_all_subjects$Status)
  subject <- as.factor(df_counts_no_0_all_subjects$individual)
  
  message("Estimating Dispersion Using Gamma-Poisson...")
  cluster_size <- ncol(count_matrix_final)
  
  size_factors <- calculateSumFactors(count_matrix_final,
                                      clusters = NULL,
                                      ref.clust = NULL,
                                      max.cluster.size = cluster_size,
                                      positive = TRUE,
                                      scaling = NULL,
                                      min.mean = NULL,
                                      subset.row = NULL)
  
  coldata <- data.frame(Status)
  fit <- glm_gp(as.matrix(count_matrix_final), col_data = coldata, size_factors = size_factors, design = ~ Status, on_disk = FALSE)
  dispersions.RAW <- fit$overdispersion_shrinkage_list$ql_disp_estimate
  log.sizeFactors.RAW <- log(size_factors)
  
  offset <- log.sizeFactors.RAW
  index <- v
  
  results <- mclapply(1:index, function(i) {
    message('Calculating Gene', paste0(i, ' of ', index, '...'))
    
    dispersion <- dispersions.RAW[i]
    outcome <- df_counts_no_0_all_subjects[, i]
    
    glmer_model <- glmer(outcome ~ Status + (1 | subject),
                         offset = offset,
                         family = MASS::negative.binomial(theta = 1 / dispersion))
    
    pairs <- pairs(emmeans(glmer_model, 'Status'), adjust = 'none')
    
    output_df <- data.frame(
      gene = valid_genes[i],
      f_m_estimate = as.data.frame(pairs[pairs@grid$contrast == 'F - M'])[, 2], # Ensure the correct reference to pairs columns
      f_m_p.value = as.data.frame(pairs[pairs@grid$contrast == 'F - M'])[, 6],
      d_m_estimate = as.data.frame(pairs[pairs@grid$contrast == 'D - M'])[, 2],
      d_m_p.value = as.data.frame(pairs[pairs@grid$contrast == 'D - M'])[, 6],
      d_f_estimate = as.data.frame(pairs[pairs@grid$contrast == 'D - F'])[, 2],
      d_f_p.value = as.data.frame(pairs[pairs@grid$contrast == 'D - F'])[, 6],
      warning = ifelse(length(glmer_model@optinfo$conv$lme4$code) != 0, substr(glmer_model@optinfo$conv$lme4$messages, 1, 50), NA),
      singular = ifelse(isSingular(glmer_model), TRUE, FALSE)
    )
    return(output_df)
  }, mc.cores = n_cores)
  
  results <- do.call(rbind, results)
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  results$f_m_p.value <- ifelse(results$f_m_p.value == 0, 1, results$f_m_p.value)
  results$d_m_p.value <- ifelse(results$d_m_p.value == 0, 1, results$d_m_p.value)
  results$d_f_p.value <- ifelse(results$d_f_p.value == 0, 1, results$d_f_p.value)
  
  results$f_m_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$f_m_p.value), method = "fdr"), "NA")
  results$d_m_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$d_m_p.value), method = "fdr"), "NA")
  results$d_f_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$d_f_p.value), method = "fdr"), "NA")
  
  assign(paste0("results_", "cluster", cluster), results, envir = .GlobalEnv)
  
  message('Complete')
  end_time <- Sys.time()  # End timing
  message(end_time - start_time)  # Print the time difference
  return(results)
}

neg_bin <- data.frame()
for (i in 0:2) {
  cluster <- paste0('19_',i)
  print(cluster)
  output <- neg.bin.mult(cluster_19,
                         'sub',
                         cluster)
  output$cluster = cluster
  neg_bin <- rbind(neg_bin, output)
}

 neg_bin$f_m_q.value <- p.adjust(neg_bin$f_m_p.value, 'fdr', nrow(neg_bin))
 neg_bin$d_m_q.value <- p.adjust(neg_bin$d_m_p.value, 'fdr', nrow(neg_bin))
 neg_bin$d_f_q.value <- p.adjust(neg_bin$d_f_p.value, 'fdr', nrow(neg_bin))

 neg_bin_defined <- define_degs(neg_bin)

 neg_bin_defined_filtered <- neg_bin_defined%>%
   filter(!is.na(issignif))
 
 neg_bin_defined_counts<- neg_bin_defined_filtered%>%
   group_by(cluster, class)%>%
   summarize(class_count = n())
 
deg_counts <-ggplot(neg_bin_defined_counts, aes(x = cluster, y = class_count, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
    geom_bar(position="stack", stat="identity")
 deg_counts
 
  neg_bin_defined_counts_no_singular<- neg_bin_defined_filtered%>%
    filter(singular !=T)%>%
   group_by(cluster, class)%>%
   summarize(class_count = n())
  
deg_counts_no_singular<-ggplot(neg_bin_defined_counts_no_singular, aes(x = cluster, y = class_count, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
    geom_bar(position="stack", stat="identity")
 deg_counts_no_singular


neg_bin_defined_filtered[neg_bin_defined_filtered$cluster=='19_0' ,]
neg_bin_defined_filtered[neg_bin_defined_filtered$cluster=='19_1' ,]
neg_bin_defined_filtered[neg_bin_defined_filtered$cluster=='19_2' ,]

### GO of DEGs####
genes_19_0 <- neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_0']
go_19_0 <- clown_go(genes_19_0)
dotplot(go_19_0)

genes_19_0_late_up <- neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_0'& neg_bin_defined_filtered$class=='Late Upregulated']
go_19_0_late_up <- clown_go(genes_19_0_late_up)
dotplot(go_19_0_late_up)

genes_19_0_term_down <- neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_0'& neg_bin_defined_filtered$class=='Terminally Downregulated']
go_19_0_term_down <- clown_go(genes_19_0_term_down)
dotplot(go_19_0_term_down)


genes_19_1 <- neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_1']
go_19_1 <- clown_go(genes_19_1)
dotplot(go_19_1)

genes_19_1_term_down <- neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_1'& neg_bin_defined_filtered$class=='Terminally Downregulated']
go_19_1_term_down <- clown_go(genes_19_1_term_down)
dotplot(go_19_1_term_down)


genes_19_2 <- neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_2']
go_19_2 <- clown_go(genes_19_2)
dotplot(go_19_2)

genes_19_2_late_up<- neg_bin_defined_filtered$gene[neg_bin_defined_filtered$cluster=='19_2'& neg_bin_defined_filtered$class=='Late Upregulated']
go_19_2_late_up <- clown_go(genes_19_2_late_up)
dotplot(go_19_2_late_up)


#### Proportion Differences ####
total_cells_m <- nrow(cluster_19@meta.data[cluster_19@meta.data$Status=='M',])
total_cells_f <- nrow(cluster_19@meta.data[cluster_19@meta.data$Status=='F',])
total_cells_d <- nrow(cluster_19@meta.data[cluster_19@meta.data$Status=='D',])

cluster_19_cells <- cluster_19@meta.data%>%
  filter(Status == 'F' | Status == 'M'| Status == 'D')%>%
  group_by(sub, Status)%>%
  summarize(n())

cluster_19_cells$total <- NA
cluster_19_cells$total[cluster_19_cells$Status=='D'] <- total_cells_d
cluster_19_cells$total[cluster_19_cells$Status=='F'] <- total_cells_f
cluster_19_cells$total[cluster_19_cells$Status=='M'] <- total_cells_m

cluster_19_cells$in_cluster <- cluster_19_cells[,3]
cluster_19_cells$non <- cluster_19_cells$total-cluster_19_cells$in_cluster
cluster_19_cells$non <- unlist(cluster_19_cells$non)
cluster_19_cells$in_cluster <- unlist(cluster_19_cells$in_cluster)

counts_19_0 <- subset(cluster_19_cells, sub == '19_0')

cluster_19_0_matrix <- matrix(NA,nrow(counts_19_0) ,2)
cluster_19_0_matrix[,1] <- counts_19_0$in_cluster
cluster_19_0_matrix[,2] <- counts_19_0$non

cluster_19_0_glm <- glm(cluster_19_0_matrix~Status, family = binomial('logit'), data=counts_19_0)
summary(cluster_19_0_glm)
pairs(emmeans(cluster_19_0_glm, 'Status'),adjust = 'none')

counts_19_0$prop <- counts_19_0$n/counts_19_0$total
counts_19_0$Status <- factor(counts_19_0$Status, levels = c('M','D','F'))
ggplot(counts_19_0, aes(x = Status, y = prop, fill = Status))+
  geom_bar(stat = 'identity')


counts_19_1 <- subset(cluster_19_cells, sub == '19_1')

cluster_19_1_matrix <- matrix(NA,nrow(counts_19_1) ,2)
cluster_19_1_matrix[,1] <- counts_19_1$in_cluster
cluster_19_1_matrix[,2] <- counts_19_1$non

cluster_19_1_glm <- glm(cluster_19_1_matrix~Status, family = binomial('logit'), data=counts_19_1)
summary(cluster_19_1_glm)
pairs(emmeans(cluster_19_1_glm, 'Status'),adjust = 'none')

counts_19_1$prop <- counts_19_1$n/counts_19_1$total
counts_19_1$Status <- factor(counts_19_1$Status, levels = c('M','D','F'))
ggplot(counts_19_1, aes(x = Status, y = prop, fill = Status))+
  geom_bar(stat = 'identity')


counts_19_2 <- subset(cluster_19_cells, sub == '19_2')

cluster_19_2_matrix <- matrix(NA,nrow(counts_19_2) ,2)
cluster_19_2_matrix[,1] <- counts_19_2$in_cluster
cluster_19_2_matrix[,2] <- counts_19_2$non

cluster_19_2_glm <- glm(cluster_19_2_matrix~Status, family = binomial('logit'), data=counts_19_2)
summary(cluster_19_2_glm)
pairs(emmeans(cluster_19_2_glm, 'Status'),adjust = 'none')

counts_19_2$prop <- counts_19_2$n/counts_19_2$total
counts_19_2$Status <- factor(counts_19_2$Status, levels = c('M','D','F'))
ggplot(counts_19_2, aes(x = Status, y = prop, fill = Status))+
  geom_bar(stat = 'identity')

### Lasso Regression ####
        lasso_scorer <- function(deg_data, cluster, multiome_object = cluster_19){
          library(glmnet)
#subset data
  data<- subset(deg_data, cluster == cluster & !is.na(gene))
  #remove clusters with no DEGs
  if(nrow(data)<1){return(NULL)}
  #coerce to numeric
  data$d_f_q.value <- as.numeric(data$d_f_q.value)
  data$d_m_q.value <- as.numeric(data$d_m_q.value)
  data$f_m_q.value <- as.numeric(data$f_m_q.value)
  
  #list degs

  degs <- data$gene[data$f_m_q.value<0.05|
                      data$d_f_q.value<0.05|
                      data$d_m_q.value<0.05]
  degs <- degs[!is.na(degs)]
  
  if(length(degs) ==0){return(NULL)}

degs_expression <- data.frame()
for(i in degs){
  gene_expression <-mean_expression_cluster_data(
    object = cluster_19,
    i,
    cluster,
    clustering ='sub'
    )

  new_data <-data.frame(
    cluster= cluster,
    gene = i,
    mean = gene_expression$mean,
    sex = gene_expression$Sex,
    individual = gene_expression$individual
    )
  degs_expression <- rbind(new_data, degs_expression)
}

#Define sexes as binomial  
degs_expression$fish.class <- ifelse(degs_expression$sex == 'F',1,NA)
degs_expression$fish.class <- ifelse(degs_expression$sex == 'M',0,degs_expression$fish.class)
#remove the other weird sexes
degs_expression <- subset(degs_expression, sex == 'F'|
                            sex=='M'|
                            sex == 'D')

#pivot to make matrix
pivoted_data<- degs_expression%>%
  pivot_wider(names_from = gene, 
              values_from = mean)

#training should only be males and females
training_data <- pivoted_data[!is.na(pivoted_data$fish.class),]
training_data <- training_data[complete.cases(training_data[, 5:ncol(training_data)]), ]




x.training <- as.matrix(na.omit(training_data[,5:ncol(training_data)]))
  if(ncol(x.training)<2){return(NULL)}

y.training <- training_data$fish.class

#calculate lambda
lasso <-cv.glmnet(y = y.training, x = x.training, family = "binomial", alpha = 1, lambda = NULL)


#use as lambda in final trainer
min <- lasso$lambda.min

#train
lasso.final <- glmnet(x=x.training, y=y.training, alpha = 1, family = "binomial",
                      lambda = min)

#now predict dominants
test_data <- pivoted_data[is.na(pivoted_data$fish.class),]
test_data <- test_data[complete.cases(test_data[, 5:ncol(test_data)]), ]


x.test <- as.matrix(na.omit(test_data[,5:ncol(test_data)]))
  if(ncol(x.test)<2){return(NULL)}


#predict probabilities
probabilities <- lasso.final %>% predict(newx = x.test, type = 'response')


#make a dataframe with results
predicted.classes <- ifelse(probabilities > 0.5, 'f', "m")

predicted_data <- as.data.frame(probabilities)
made.data <- as.data.frame(probabilities)
made.data$predicted <- predicted.classes
made.data$fish <- test_data$individual
made.data$probabilities <- probabilities

return(made.data)

}


score_data <- data.frame()
for(i in 0:2){
  
  cluster = paste0('19_',i)
probs <- lasso_scorer(deg_data =results_cluster19_0, cluster = cluster, multiome_object = cluster_19)
probs$cluster = cluster
score_data <- rbind(score_data, probs)

}
library(forcats)

ggplot(subset(score_data, cluster == '19_0'), aes(x = fct_reorder(.f=fish, .x = probabilities), y = probabilities, color = predicted))+
  geom_point()+
  geom_text(data = subset(score_data, cluster == '19_0'),aes(label = fish, vjust = 2))



ggplot(subset(score_data, cluster == '19_1'), aes(x = fct_reorder(.f=fish, .x = probabilities), y = probabilities, color = predicted))+
  geom_point()+
  geom_text(data = subset(score_data, cluster == '19_1'),aes(label = fish, vjust = 2))

ggplot(subset(score_data, cluster == '19_2'), aes(x = fct_reorder(.f=fish, .x = probabilities), y = probabilities, color = predicted))+
  geom_point()+
  geom_text(data = subset(score_data, cluster == '19_2'),aes(label = fish, vjust = 2))

cell_clusters <- Idents(cluster_19)

cluster_labels <- score_data%>%
  group_by(cluster)%>%
  summarize(mean = mean(probabilities))

cluster_to_predicted_sex <- setNames(cluster_labels$mean, cluster_labels$cluster)

metadata <- cluster_to_predicted_sex[as.character(cell_clusters)]
names(metadata) <- names(cell_clusters)
cluster_19 <- AddMetaData(cluster_19, metadata, col.name = "probabilities")

FeaturePlot(cluster_19, feature = 'probabilities')


mean(score_data$probabilities[score_data$cluster=='19_0'])
mean(score_data$probabilities[score_data$cluster=='19_1'])
mean(score_data$probabilities[score_data$cluster=='19_2'])


##### Continuity #####
Idents(cluster_19) <- 'harmony.wnn_res0.4_clusters'
DimPlot(cluster_19)

#First, I want to make a PCA based on the DEGs
cluster_19_degs <- read.csv('/Users/ggraham/Desktop/snRNA-seq R Files 122524/Seurat Outputs/122324 Neg Bin with Doms/cluster_19.csv')
cluster_19_degs <- cluster_19_degs%>%
  subset(d_m_q.value <0.05|
           d_f_q.value<0.05|
           f_m_q.value<0.05)
deg_list <- cluster_19_degs$gene

pca_data <- data.frame()
for(i in deg_list){
  data <- mean_expression_cluster_data(cluster_19,
                                       i,
                                       19)
  data$gene <- i
  pca_data <- rbind(pca_data, data)
}

pca_data_pivoted <- pca_data %>%
  dplyr::select(individual, mean, Sex, gene)%>%
  dplyr::filter(Sex == 'D'| Sex =='F' | Sex == 'M')%>%
  pivot_wider(names_from = gene, 
              values_from = mean)

cluster_19_prcomp <- prcomp(pca_data_pivoted[,3:ncol(pca_data_pivoted)])
fviz_pca_biplot(cluster_19_prcomp)+
  geom_text(label = pca_data_pivoted$individual)+
  geom_point(aes(color = pca_data_pivoted$Sex))

cluster_19_pca_loadings <- cluster_19_prcomp$x
       
cluster_19_pca_loadings <- as.data.frame(cluster_19_pca_loadings[,1:2])
cluster_19_pca_loadings$Sex <- pca_data_pivoted$Sex

grouped_means <- cluster_19_pca_loadings%>%
  group_by(Sex)%>%
    summarize(across(starts_with("PC"), mean))

mean_m <- grouped_means[grouped_means$Sex=='M',2:3]

mean_f <- grouped_means[grouped_means$Sex=='F',2:3]

mean_d <- grouped_means[grouped_means$Sex=='D',2:3]

f_m_distance <- stats::dist(rbind(as.numeric(mean_m), as.numeric(mean_f)))

d_m_distance <- stats::dist(rbind(as.numeric(mean_m), as.numeric(mean_d)))

d_f_distance <- stats::dist(rbind(as.numeric(mean_d), as.numeric(mean_f)))

continuum_score <- f_m_distance/(d_m_distance+d_f_distance)

ggplot(grouped_means, aes(x = PC1, y = PC2, color = Sex))

##### How old is cluster 19 really ###
neuron_only <- subset(obj, harmony.wnn_res0.4_clusters !=2&
                        harmony.wnn_res0.4_clusters !=14&
                        harmony.wnn_res0.4_clusters !=29&
                        harmony.wnn_res0.4_clusters !=22&
                        harmony.wnn_res0.4_clusters !=26&
                        harmony.wnn_res0.4_clusters !=28&
                        harmony.wnn_res0.4_clusters !=18&
                        harmony.wnn_res0.4_clusters !=4)

DotPlot(object = neuron_only, 
        group.by = "harmony.wnn_res0.4_clusters", 
        features = c('gad2','LOC111584103','slc17a6b','elavl3'),
        cols = c("#D2B4DE", "#8E44AD", "#6C3483")
) + 
  coord_flip()

#alright I'm convinced theyre all neurons


neuron_only_matrix <- as.matrix(neuron_only@assays$RNA$counts)
neuron_only_cyto <- CytoTRACE(mat = neuron_only_matrix
)

neuron_only$cyto <-neuron_only_cyto$CytoTRACE

FeaturePlot(neuron_only, 'cyto')

cyto_data <- data.frame(cluster = neuron_only$harmony.wnn_res0.4_clusters,
                        cyto = neuron_only$cyto,
                        individual = neuron_only$individual,
                        status = neuron_only$Status)

cyto_data$color <- ifelse(cyto_data$cluster==19, '19', NA)

ggplot(cyto_data, aes(x = cluster, y = cyto, fill = color))+
  geom_boxplot()+
  geom_hline(yintercept = mean(cyto_data$cyto))

cyto_model <- lmer(cyto~cluster+(1|individual), data = cyto_data)
summary(cyto_model)
car::Anova(cyto_model, type = 'III')

t.test(x=cyto_data$cyto[cyto_data$cluster ==19], y= rep(1/2, times = length(cyto_data$cyto[cyto_data$cluster ==19])))

random <- rnorm(n = length(cyto_data$cyto[cyto_data$cluster ==19]),mean= 1/2,sd = sd(cyto_data$cyto))

nll <-t.test(x=cyto_data$cyto[cyto_data$cluster ==19], y= random)
nll

#I want to look at this cyto data more later
#write.csv(cyto_data, 'X:/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/011125 all neurons cyto.csv')

##### I want to look at proportion degs in this cluster, it will be fast because its a smaller matrix so I can be inefficient

### Proportion DEGs ####

prop_deg_function <- function(object = obj, cluster = 19, clustering = 'harmony.wnn_res0.4_clusters'){
  start <- Sys.time()
  library(lme4)
  library(dplyr)
  library(parallel)
  library(car)
  
  options(dplyr.summarise.inform = FALSE)
  
   #extract counts matrix 
  counts <- t(as.matrix(object@assays$RNA$counts[, object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")]))
  
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
  

  deg_function <- function(gene){
    message(paste0('gene ', which(genes == gene), ' of ', n_genes))
    
    gene_expression <- filtered_cols_matrix[, gene, drop = FALSE]
    
    meta_data$gene <- gene_expression
    
    data_to_analyze$gene <- data_to_analyze[,1]
    
      joined_data_restrictions <- joined_data_to_analyze%>%
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

    
    glmer_model <- suppressMessages(glmer(model_matrix~Status + (1|individual), family = binomial('logit'), data = data_for_matrix))
    
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
  
  
    deg_output <- mclapply(X=genes, FUN=deg_function
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

deg_output2 <-prop_deg_function(object = cluster_19, cluster = 19, clustering = 'harmony.wnn_res0.4_clusters') 

deg_output2<- deg_output2[deg_output2$gene!='Error : Response is constant',]
deg_output2$gene[deg_output2$anova_q.value<0.05 & is.na(deg_output2$warning)]

write.csv(deg_output2, '/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/cluster 19 prop degs 011525.csv')

deg_output3 <- deg_output2[is.na(deg_output2$warning),]
deg_output3$d_f_q.value <- p.adjust(deg_output3$d_f_p.value, 'fdr', nrow(deg_output3))
deg_output3$d_m_q.value <- p.adjust(deg_output3$d_m_p.value, 'fdr', nrow(deg_output3))
deg_output3$f_m_q.value <- p.adjust(deg_output3$f_m_p.value, 'fdr', nrow(deg_output3))
deg_output3$anova_q.value <- p.adjust(deg_output3$anova_p.value, 'fdr', nrow(deg_output3))

deg_output3$color <- ifelse(deg_output3$anova_q.value<0.05 | ####### Using P values not q values because wtf is this data
                              deg_output3$d_f_q.value<0.05|
                              deg_output3$f_m_q.value<0.05| 
                              deg_output3$d_m_q.value<0.05,
                            'signif',
                            NA)

deg_output3$label <- ifelse( is.na(deg_output3$warning)&(
                              deg_output3$anova_q.value<0.1 | 
                              deg_output3$d_f_q.value<0.1|
                              deg_output3$f_m_q.value<0.1| 
                              deg_output3$d_m_q.value<0.1),
                            deg_output3$gene,
                            NA)


ggplot(deg_output3, aes(x = as.numeric(prop_expressing_F),y =  as.numeric(prop_expressing_M), color = color))+
  geom_point(alpha = 0.25, aes(size = as.numeric(prop_expressing_D)))+
  geom_text(data = deg_output3, aes(label = label), color ='blue', size =3)

FeaturePlot(data, 'pgr', split.by = 'Status')


                   

