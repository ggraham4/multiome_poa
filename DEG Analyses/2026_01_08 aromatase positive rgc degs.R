#Negative binomial lower stringency 
{
  library(parallel)
  library(factoextra)
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
  library(Polychrome)
  P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
names(P40) <- NULL

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")

a <- listAttributes(ensembl)
biomart_basic <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'entrezgene_description'))


namer = function(gene){
  gene_name = biomart_basic$entrezgene_description[biomart_basic$entrezgene_accession==gene]
  return(gene_name[1])
  
}


  mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
clown_go<- readRDS('Functions/clown_go')
define_degs<- readRDS('Functions/define_degs')

latent_plotter = function(gene,  cluster, clustering ='res0.8_50nn_40PC_45LSI'){
  
  temp = obj@meta.data
  temp$exp = obj@assays$RNA$data[gene,]
  
  temp_subset = temp[temp[clustering] == cluster,]
  temp_out = temp_subset%>%
    group_by(Status.x, individual)%>%
    summarize(mean_exp = mean(exp),
              SexState = mean(SexState))
  
  p = ggplot(temp_out, aes(x = SexState, y = mean_exp))+
    geom_point(aes( color = Status.x))+
    labs(title = paste0(gene,'_', cluster))+
    geom_smooth()
  
  return(p)
}

}

obj <- readRDS('~/Desktop/optimal_clustering_rna_only.rds')


obj$arom = ifelse(obj@meta.data$res0.8_50nn_40PC_45LSI == 1 & 
                    obj@assays$RNA$data['LOC111577263',]>0, 
                  T,
                  F)

DimPlot(obj, group.by = 'arom')


 cluster=T
 clustering = 'arom'
 n_cores = detectCores() - 1
 
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
  
  
  ##### ok here is where my changes are going to have to be #####
  
  df_counts_no_0 <- cbind(df_counts_no_0, df_counts_meta)
  df_counts_no_0_split_by_subject <- split(df_counts_no_0, f = df_counts_no_0$individual)
  
  message("Finding Good Genes for Subject...")
  # REMOVE GENES WITH ZERO COUNTS IN EACH SUBJECT 
  for (l in 1:length(df_counts_no_0_split_by_subject)) {
    correct_gene_names <- colnames(df_counts_no_0)
    
    temp_subject_l <- data.frame(df_counts_no_0_split_by_subject[[l]]) ### AND HERE THEY GET FUCKED UP 
    colnames(temp_subject_l) <- correct_gene_names
    
    temp_subject_l_counts <- temp_subject_l[, 1:n_genes_no_0]
    ###
    temp_subject_l_counts_no_0 <- temp_subject_l_counts #<- temp_subject_l_counts[, which(colSums(temp_subject_l_counts) != 0)]
    #ok here I am making the real code a comment to stop the filtering without fucking up the rest of the code
    #I'm realizing there does need to be some way I test for a gene being missing in several sexes huh
    out <- data.frame(colnames(temp_subject_l_counts_no_0))
    assign(x = paste0("gene_list_subject_", l), value = get("out"))
  }
  
  # GENERATE A LIST OF GENES FOR EACH SUBJECT -- I believe everything else should still work
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
    
    
     tryCatch({suppressMessages(
    glmer_model <- glmer(outcome ~ Status + (1 | subject),
                         offset = offset,
                         family = MASS::negative.binomial(theta = 1 / dispersion)))
     }, error = function(e) {
      return(NULL)
    })
    if(!exists('glmer_model')){
      glmer_model = NULL
    }
    
    if(!is.null(glmer_model)){
    
    pairs <- pairs(emmeans(glmer_model, 'Status'), adjust = 'none')
    
  av = car::Anova(glmer_model, type = 3)
    output_df <- data.frame(
      gene = valid_genes[i],
      f_m_estimate = as.data.frame(pairs[pairs@grid$contrast == 'F - M'])[, 2], # Ensure the correct reference to pairs columns
      f_m_p.value = as.data.frame(pairs[pairs@grid$contrast == 'F - M'])[, 6],
      d_m_estimate = as.data.frame(pairs[pairs@grid$contrast == 'D - M'])[, 2],
      d_m_p.value = as.data.frame(pairs[pairs@grid$contrast == 'D - M'])[, 6],
      d_f_estimate = as.data.frame(pairs[pairs@grid$contrast == 'D - F'])[, 2],
      d_f_p.value = as.data.frame(pairs[pairs@grid$contrast == 'D - F'])[, 6],
      av_p.value = av$`Pr(>Chisq)`[2],    
      warning = ifelse(length(glmer_model@optinfo$conv$lme4$code) != 0, substr(glmer_model@optinfo$conv$lme4$messages, 1, 50), NA),
      singular = ifelse(isSingular(glmer_model), TRUE, FALSE)
    )
    return(output_df)
    }
    
  }, mc.cores = n_cores
  )
  
  results <- do.call(rbind, results)
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  results$f_m_p.value <- ifelse(results$f_m_p.value == 0, 1, results$f_m_p.value)
  results$d_m_p.value <- ifelse(results$d_m_p.value == 0, 1, results$d_m_p.value)
  results$d_f_p.value <- ifelse(results$d_f_p.value == 0, 1, results$d_f_p.value)
  
  results$f_m_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$f_m_p.value), method = "fdr", n = nrow(results)), "NA")
  results$d_m_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$d_m_p.value), method = "fdr", n = nrow(results)), "NA")
  results$d_f_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$d_f_p.value), method = "fdr", n = nrow(results)), "NA")
  
  results$av_q.value <- ifelse(test = is.na(results$warning),p.adjust(as.numeric(results$av_p.value), method = 'fdr',nrow(results)), "NA")
  
  assign(paste0("results_", "cluster", cluster), results, envir = .GlobalEnv)
  
  message('Complete')
  end_time <- Sys.time()  # End timing
  message(end_time - start_time)  # Print the time difference


out = results
out$cluster ='Aromatase+ RGC'
#write.csv(out, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2026_01_08 aromatase positive RGC'))


define_degs3 = define_degs_2 = function(data, alpha = 0.05){
  ref_data = read.csv("DEG Analyses/Expression DEGs/pairwise_patterns.csv")
  
  significant = subset(data, av_q.value < alpha)
  if(nrow(significant)<1){return(NULL)}
  significant$full_label = NA
  significant$first_word = NA
  significant$second_word = NA
  significant$short_label = NA

  newd <- data.frame()
  
  for(index in 1:nrow(significant)){
    data_of_interest = significant[index, ]
    
    if(data_of_interest$d_m_p.value < 0.05){
      ref_locs_d_m = which(ref_data$d_m != 'NS')
    } else{
      ref_locs_d_m = which(ref_data$d_m == 'NS')
    }
    
    if(data_of_interest$f_m_p.value < 0.05){
      ref_locs_f_m = which(ref_data$f_m != 'NS')
    } else{
      ref_locs_f_m = which(ref_data$f_m == 'NS')
    }
    
    if(data_of_interest$d_f_p.value < 0.05){
      ref_locs_d_f = which(ref_data$d_f != 'NS')
    } else{
      ref_locs_d_f = which(ref_data$d_f == 'NS')
    }
    
    locs_with_correct_signif_pattern = Reduce(intersect, list(ref_locs_d_f, ref_locs_f_m, ref_locs_d_m))
    
    if(data_of_interest$d_m_p.value < 0.05){
      if(data_of_interest$d_m_estimate < 0){
        dir_locs_d_m = which(ref_data$d_m == '<')
      } else{
        dir_locs_d_m = which(ref_data$d_m == '>')
      }
    } else{
      dir_locs_d_m = which(ref_data$d_m == 'NS')
    }
    
    if(data_of_interest$f_m_p.value < 0.05){
      if(data_of_interest$f_m_estimate < 0){
        dir_locs_f_m = which(ref_data$f_m == '<')
      } else{
        dir_locs_f_m = which(ref_data$f_m == '>')
      }
    } else{
      dir_locs_f_m = which(ref_data$f_m == 'NS')
    }
    
    if(data_of_interest$d_f_p.value < 0.05){
      if(data_of_interest$d_f_estimate < 0){
        dir_locs_d_f = which(ref_data$d_f == '<')
      } else{
        dir_locs_d_f = which(ref_data$d_f == '>')
      }
    } else{
      dir_locs_d_f = which(ref_data$d_f == 'NS')
    }
    
    dir_locs = Reduce(intersect, list(dir_locs_d_f, dir_locs_d_m, dir_locs_f_m))
    
    final_loc = intersect(dir_locs, locs_with_correct_signif_pattern)
    
    if(length(final_loc) > 0){
      full_label = ref_data[final_loc, ]$Full_label
      first_word = ref_data[final_loc, ]$First_word
      second_word = ref_data[final_loc, ]$Second_word
      short_label = ref_data[final_loc, ]$Short_label

      data_of_interest$full_label = full_label
      data_of_interest$first_word = first_word
      data_of_interest$second_word = second_word
      data_of_interest$short_label = short_label

      newd = rbind(newd, data_of_interest)
    }
  }
  return(newd)     
}

out_defined = define_degs3(out)
gene_namer = readRDS('Functions/gene_namer.rds')
out_defined$name = sapply(out_defined$gene, namer)

#write.csv(out_defined, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2026_01_08 defined aromatase positive RGC'))
out_defined = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2026_01_08 defined aromatase positive RGC')
out = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2026_01_08 aromatase positive RGC')

clown_go2 = readRDS('Functions/clown_go2')
clown_go2(out$gene[out$av_q.value<0.05],p=.05)%>%dotplot()

clown_go2(out_defined$gene[out_defined$av_q.value<0.05 & out_defined$second_word == 'Downregulated'],p=.05)%>%dotplot()+
  labs(title = 'Downregulated')

clown_go2(out_defined$gene[out_defined$av_q.value<0.05 & out_defined$second_word == 'Upregulated'],p=.05)%>%dotplot()+
  labs(title = 'Upregulated')

clown_go2(out_defined$gene[out_defined$av_q.value<0.05 & out_defined$first_word == 'Transiently'],p=.05)%>%dotplot()+
  labs(title = 'Transiently')

clown_go2(out_defined$gene[out_defined$av_q.value<0.05 & out_defined$first_word == 'Early'],p=.05)%>%dotplot()+
  labs(title = 'Early')

clown_go2(out_defined$gene[out_defined$av_q.value<0.05 & out_defined$first_word == 'Late'],p=.05)%>%dotplot()+
  labs(title = 'Late')

table(out_defined$full_label)

clown_gsea= readRDS('Functions/clown_gsea.rds')
sorter = readRDS('Functions/gsea_sorter.rds')

rownames(out) = out$gene
clown_gsea(sorter(out, 'd_m_estimate'))%>%dotplot()
clown_gsea(sorter(out, 'f_m_estimate'))%>%dotplot()
clown_gsea(sorter(out, 'd_f_estimate'))%>%dotplot() # nothing

### 
mean_expression_cluster_plot(obj, 'nkx2.2a', T, 'arom') # promotes differentiation, acts with pdgrfa

#together
mean_expression_cluster_plot(obj, 'gli1', T, 'arom') # promotes opc
mean_expression_cluster_plot(obj, 'sox1a', T, 'arom')#promotes neuronal
mean_expression_cluster_plot(obj, 'LOC111568198', T, 'arom') # promotes neuronal
#LOC111568198  pbx1 


all_degs = read.csv('DEG Outputs/FINAL degs classified w singular.csv')
all_degs_1 = subset(all_degs, cluster ==1)

mean_expression_cluster_plot(obj, 'nkx2.2a', 1, 'res0.8_50nn_40PC_45LSI') 
mean_expression_cluster_plot(obj, 'sox1a', 1, 'res0.8_50nn_40PC_45LSI') 
mean_expression_cluster_plot(obj, 'gli1', 1, 'res0.8_50nn_40PC_45LSI') 
mean_expression_cluster_plot(obj, 'LOC111568198', 1, 'res0.8_50nn_40PC_45LSI') 
# consistent with the aro+ cells

# is there a difference in OPC proportion
obj@meta.data$Status = factor(obj@meta.data$Status, levels =c('NRM',
                                                              'M',
                                                              'D',
                                                              'E',
                                                              'NF',
                                                              "F"))
n_cells_clust = obj@meta.data%>%
  group_by(res0.8_50nn_40PC_45LSI, individual, Status)%>%
  summarize(n =n())

n_cells = obj@meta.data%>%
  group_by(individual, Status)%>%
  summarize(nTotal =n())

joined = n_cells_clust%>%
  left_join(n_cells, by = 'individual')%>%
  mutate(prop_cells = n/nTotal)

# do they have more glia? what is the evidence of this increased glial differentiation?
ggplot(subset(joined, res0.8_50nn_40PC_45LSI == 1), aes(x = Status.x, y = prop_cells))+
    geom_boxplot()+
  geom_point()
# not really a convincing difference

ggplot(subset(joined, res0.8_50nn_40PC_45LSI == 26), aes(x = Status.x, y = prop_cells))+
    geom_boxplot()+
  geom_point()
# no

ggplot(subset(joined, res0.8_50nn_40PC_45LSI == 13), aes(x = Status.x, y = prop_cells))+
    geom_boxplot()+
  geom_point()
# doesnt seem like there are more OPCs

ggplot(subset(joined, res0.8_50nn_40PC_45LSI == 2), aes(x = Status.x, y = prop_cells))+
    geom_boxplot()+
  geom_point()
#MAYBE more Oligos but like thats not convincing

ggplot(subset(joined, res0.8_50nn_40PC_45LSI == 15), aes(x = Status.x, y = prop_cells))+
    geom_boxplot()+
  geom_point()
# no difference in ependymal

# what the hell is happening?

sub_1 = FindSubCluster(obj, 1, 'harmony.wsnn')
sub_1 = subset(sub_1, res0.8_50nn_40PC_45LSI==1)
Idents(sub_1) ='sub.cluster'
DimPlot(sub_1)

n_cells_subclust = sub_1@meta.data%>%
  group_by(sub.cluster, individual, Status)%>%
  summarize(n =n())

n_cells_sub = sub_1@meta.data%>%
  group_by(individual, Status)%>%
  summarize(nTotal =n())

joined_sub = n_cells_subclust%>%
  left_join(n_cells_sub, by = 'individual')%>%
  mutate(prop_cells = n/nTotal)

ggplot(joined_sub, aes(x = Status.x, y = prop_cells))+
    geom_boxplot()+
  geom_point()+
  facet_wrap(~sub.cluster, scales = 'free')
 
# very much struggling to find a glial cluster where females are higher, 1_3 is a progenitor
# could it be that 


# lets look at all of the 1 TFs

#up
mean_expression_cluster_plot(obj, 'nkx2.2a', 1, 'res0.8_50nn_40PC_45LSI') 
mean_expression_cluster_plot(obj, 'foxp2', 1, 'res0.8_50nn_40PC_45LSI') 
mean_expression_cluster_plot(obj, 'sox6', 1, 'res0.8_50nn_40PC_45LSI') 

#transient up 
mean_expression_cluster_plot(obj, 'LOC111587547', 1, 'res0.8_50nn_40PC_45LSI') #homeobox protein DBX1-B-like
mean_expression_cluster_plot(obj, 'znf710b', 1, 'res0.8_50nn_40PC_45LSI')

#down
mean_expression_cluster_plot(obj, 'sox1a', 1, 'res0.8_50nn_40PC_45LSI') 
mean_expression_cluster_plot(obj, 'gli1', 1, 'res0.8_50nn_40PC_45LSI') 
mean_expression_cluster_plot(obj, 'LOC111568198', 1, 'res0.8_50nn_40PC_45LSI') #pbx1
mean_expression_cluster_plot(obj, 'skida1', 1, 'res0.8_50nn_40PC_45LSI') 
mean_expression_cluster_plot(obj, 'LOC111563338', 1, 'res0.8_50nn_40PC_45LSI') # COUP tf2 like
# implicated in sex hormone regulation

#epigenetic
mean_expression_cluster_plot(obj, 'LOC111584294', 1, 'res0.8_50nn_40PC_45LSI') 
mean_expression_cluster_plot(obj, 'dnmt3bb.1', 1, 'res0.8_50nn_40PC_45LSI') 


# latent var
latent = read.csv('Measures/2025_12_26 all_data.csv')
latent$individual = latent$Fish

obj@meta.data= obj@meta.data%>%
  left_join(latent, by = 'individual')

#epigenetic
latent_plotter('dnmt3bb.1', 1)
latent_plotter('LOC111584294', 1)

# tfs
latent_plotter('LOC111563338', 1) # coup2 tf like , this onw is super interesting
latent_plotter('gli1', 1) 
latent_plotter('sox1a', 1) 

# i would like to put all of these together
transcription_factors =c(
  'nkx2.2a',
  'foxp2',
  'sox6',
  'LOC111587547', #dbx1-like
  'znf710b',
  'sox1a',
  'gli1',
  'LOC111568198', #pbx1
  'dnmt3bb.1', #methylation
  'LOC111584294', #bromodomain containing
  'LOC111563338' # couptf2
    )

Idents(sub_1)='res0.8_50nn_40PC_45LSI'
avg_expr <- AverageExpression(
  sub_1, 
  features = transcription_factors,
  group.by = "individual",
  slot = "data" 
)

gex = avg_expr[['RNA']]%>%t()%>%as.data.frame()
gex$individual = rownames(gex)

renamed = list(
    'nkx2.2a' ='nkx2.2a' ,
  'foxp2'='foxp2',
  'sox6' = 'sox6',
  'LOC111587547'= 'dbx-1 like',
  'znf710b' = 'znf710b',
  'sox1a' = 'sox1a',
  'gli1' = 'gli1',
  'LOC111568198'= 'pbx1',
  'dnmt3bb.1' = 'dnmt3bb.1',
  'LOC111584294' = 'brd3',
  'LOC111563338' = 'COUP TFII like'
)


gex_lat = latent %>%
  left_join(gex, by = 'individual') %>%
  mutate(across(all_of(transcription_factors), scale)) %>%  # Scale each TF
  pivot_longer(cols = transcription_factors)

gex_lat$name2 = renamed[gex_lat$name]%>%unlist()

ggplot(gex_lat, aes(x = SexState, y = value))+
  geom_smooth(aes(color = name2), se = F)+
  geom_rug(aes(x = SexState), sides = 'b')+
  labs(x = 'Latent Sex Score', y = 'Z-Score Expression', color ='Gene')
# foxp1 and nkx2.2 covary # and maybe sox6
# sox1 and gli1 and dnmt
# dbx1 and znf10b


##correlation matirx
b = avg_expr[['RNA']]%>%t()%>%as.data.frame()
colnames(b) = renamed[colnames(b)]%>%unlist()
mat = cor(b )
heatmap(mat)

prcomp(mat)%>%fviz_pca_ind()

mat2 = cov( avg_expr[['RNA']]%>%t()%>%as.data.frame())
heatmap(mat2)

