# IEGs based on PC1

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
    summarize(mean_ieg = mean(ieg_PC1), .groups = 'drop') %>%
    subset(Status %in% statuses)
  
  # Create plot
  p = ggplot(plot_data, aes(x = Status, y = mean_ieg, color = !!sym(gene_name))) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(.5)) +
    facet_wrap(~sub.cluster, scales = 'free') +
    labs(title = paste("IEG expression by", gene_name, "status"),
         y = "IEG PC1",
         color = gene_name) +
    theme_minimal()
  
  return(p)
}

  #define functions

prop_cluster_plot=function(object, gene, cluster, clustering = 'final_clusters'){
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




}
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_6 = FindSubCluster(obj, 6, graph.name='harmony.wsnn')
Idents(sub_6) <- 'sub.cluster'
sub_6 = subset(sub_6, final_clusters ==6)
sub_6$sub.cluster = factor(sub_6$sub.cluster, levels = c(paste0('6_',0:3)))
sub_6$Status = factor(sub_6$Status, levels = c('NRM','M',"D",'E','NF','F'))

DimPlot(sub_6)

# Find IEGs ---
iegs <- c('npas4a', 'fosb','egr1') #using fosb instead of cfos , I think it should be more or less the same

iegs_= calculateCoexpressedGenes(sub_6, iegs) #c-fos like is included anyway so its good

### Define a function to return PC1 for the module
whole_matrix = sub_6@assays$RNA$data%>%as.matrix()%>%t()
module_pc1 = function(module){
  genes = colnames(whole_matrix)
  module_matrix = whole_matrix[,which(genes %in% module)]
  
  pc = prcomp(module_matrix, scale = T)
  scores <- pc$x[,1]

  if(mean(pc$rotation[c("fosb","egr1","npas4a"),1], na.rm=TRUE) < 0){
  scores <- -scores
}
return(scores)

}

sub_6$ieg_PC1 = module_pc1(iegs_)

# check to make sure it makes sense
DotPlot(sub_6, 'ieg_PC1')

ggplot(subset(sub_6@meta.data, Status != 'NRM'), aes(x = Status, y = ieg_PC1))+
  geom_boxplot()+
  facet_wrap(~sub.cluster, scales ='free')

ggplot(subset(sub_6@meta.data, Status != 'NRM'), aes(x = sub.cluster, y = ieg_PC1))+
  geom_boxplot()

DotPlot(sub_6, iegs)
# ok the PC score is correct I believe

## 6_3 looks interesting lets look at it rq
sub_6_3_pc1 = sub_6@meta.data%>%
  subset(sub.cluster=='6_3' & Status%in%c('M','D','F'))%>%
  group_by(individual, Status)%>%
  summarize(mean_pc1 = mean(ieg_PC1))
  
ggplot(sub_6_3_pc1, aes(x = Status, y = mean_pc1))+
  geom_boxplot()+
  geom_jitter()
# ok this is something for sure

ieg_6_3_model = lmer(ieg_PC1~Status+(1|individual), data = subset(sub_6@meta.data, sub.cluster=='6_3'&
                                                                    Status %in% c('M',"D","F")))

car::Anova(ieg_6_3_model, 3) #not quite :(

### ok, now lets do our unbiased search of 6_3 subpopulations that are IEG different ----

mat_6_3 = sub_6@assays$RNA$data[,sub_6@meta.data$sub.cluster =='6_3' &
                                  sub_6@meta.data$Status %in%
                                  c('M','D','F')]%>%
  as.matrix()%>%
  t()

meta_6_3 = subset(sub_6@meta.data, sub.cluster =='6_3' & Status %in% c('M','D','F'))

# 1) helper: genes with at least N positive cells per Status
min_cells_expr <- 3
statuses <- c('M','D','F')

nonZero_genes <- colnames(mat_6_3)[
  sapply(colnames(mat_6_3), function(gene){
    expr_counts <- tapply(mat_6_3[,gene] > 0, meta_6_3$Status, sum)
    all(expr_counts[statuses] >= min_cells_expr)
  })
]
geneIEG_lmer <- function(gene){
  gene_data = data.frame(
    individual = meta_6_3$individual,
    Status = factor(meta_6_3$Status, levels = c('M','D','F')), 
    ieg_PC1 = meta_6_3$ieg_PC1,
    expression = as.numeric(mat_6_3[,gene] > 0) 
  )

  if(sum(gene_data$expression) < 3){   
    stop("Too few expressing cells")
  }

  model = lmer(ieg_PC1 ~ Status * expression + (1|individual), data = gene_data)

  return(model)
}

geneIEG_av = function(model){
  return(car::Anova(model, type = 3))
}

geneIEG_Status_pairs = function(model){
  return(as.data.frame(pairs(emmeans(model, 'Status'), adjust = 'none')))
}

# 4) Unbiased search wrapper, with tryCatch to avoid crashes
unbiasedSearch_geneIEG = function(gene){
  model = geneIEG_lmer(gene)        
  av_ = geneIEG_av(model)
  
  out_data = data.frame(
    gene = gene,
    singular = isSingular(model),
  warning = ifelse(length(model@optinfo$conv$lme4$code) != 0, 
       substr(model@optinfo$conv$lme4$messages, 1, 50), NA),
    av_Status_p = av_$`Pr(>Chisq)`[2],
    av_expression_p = av_$`Pr(>Chisq)`[3],
    av_interaction_p = av_$`Pr(>Chisq)`[4]
  )
  return(out_data)
}

safe_unbiasedSearch_geneIEG <- function(gene){
  message(paste0(which(nonZero_genes == gene ), ' / ', length(nonZero_genes), " : ", gene))
  tryCatch({
    unbiasedSearch_geneIEG(gene)
  }, error = function(e) {
    message("ERROR: ", gene, " -> ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("WARNING (treated as NULL): ", gene, " -> ", conditionMessage(w))
    return(NULL)
  })
}

results <- lapply(nonZero_genes, safe_unbiasedSearch_geneIEG)
results <- do.call(rbind, results)

results = results[complete.cases(results[,4:6]),]
results$av_Status_q = p.adjust(results$av_Status_p, 'fdr', nrow(results))
results$av_interaction_q = p.adjust(results$av_interaction_p, 'fdr', nrow(results))
results$av_expression_q = p.adjust(results$av_expression_p, 'fdr', nrow(results))

interaction_genes = results[results$av_interaction_q < 0.1, ]
Status_genes = results[results$av_Status_q < 0.1, ]
expression_genes = results[results$av_expression_q < 0.1, ]

# i guess this is a dead end

plot_gene_ieg('kdm6bb')
plot_gene_ieg('LOC111571064')
plot_gene_ieg('cckb')
plot_gene_ieg('LOC111562384')


### trying whole cluster ----
mat_6 = sub_6@assays$RNA$data[,
                                  sub_6@meta.data$Status %in%
                                  c('M','D','F')]%>%
  as.matrix()%>%
  t()

meta_6 = subset(sub_6@meta.data, Status %in% c('M','D','F'))

# 1) helper: genes with at least N positive cells per Status
min_cells_expr <- 3
statuses <- c('M','D','F')

nonZero_genes <- colnames(mat_6)[
  sapply(colnames(mat_6), function(gene){
    expr_counts <- tapply(mat_6[,gene] > 0, meta_6$Status, sum)
    all(expr_counts[statuses] >= min_cells_expr)
  })
]
geneIEG_lmer <- function(gene){
  gene_data = data.frame(
    individual = meta_6$individual,
    Status = factor(meta_6$Status, levels = c('M','D','F')), 
    ieg_PC1 = meta_6$ieg_PC1,
    expression = as.numeric(mat_6[,gene] > 0) 
  )

  if(sum(gene_data$expression) < 3){   
    stop("Too few expressing cells")
  }

  model = lmer(ieg_PC1 ~ Status * expression + (1|individual), data = gene_data)

  return(model)
}

geneIEG_av = function(model){
  return(car::Anova(model, type = 3))
}

geneIEG_Status_pairs = function(model){
  return(as.data.frame(pairs(emmeans(model, 'Status'), adjust = 'none')))
}

# 4) Unbiased search wrapper, with tryCatch to avoid crashes
unbiasedSearch_geneIEG = function(gene){
  model = geneIEG_lmer(gene)        
  av_ = geneIEG_av(model)
  
  out_data = data.frame(
    gene = gene,
    singular = isSingular(model),
  warning = ifelse(length(model@optinfo$conv$lme4$code) != 0, 
       substr(model@optinfo$conv$lme4$messages, 1, 50), NA),
    av_Status_p = av_$`Pr(>Chisq)`[2],
    av_expression_p = av_$`Pr(>Chisq)`[3],
    av_interaction_p = av_$`Pr(>Chisq)`[4]
  )
  return(out_data)
}

safe_unbiasedSearch_geneIEG <- function(gene){
  message(paste0(which(nonZero_genes == gene ), ' / ', length(nonZero_genes), " : ", gene))
  tryCatch({
    unbiasedSearch_geneIEG(gene)
  }, error = function(e) {
    message("ERROR: ", gene, " -> ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("WARNING (treated as NULL): ", gene, " -> ", conditionMessage(w))
    return(NULL)
  })
}

results <- lapply(nonZero_genes, safe_unbiasedSearch_geneIEG)
results <- do.call(rbind, results)

results = results[complete.cases(results[,4:6]),]
results$av_Status_q = p.adjust(results$av_Status_p, 'fdr', nrow(results))
results$av_interaction_q = p.adjust(results$av_interaction_p, 'fdr', nrow(results))
results$av_expression_q = p.adjust(results$av_expression_p, 'fdr', nrow(results))

interaction_genes = results[results$av_interaction_q < 0.1, ]
Status_genes = results[results$av_Status_q < 0.1, ]
expression_genes = results[results$av_expression_q < 0.1, ]

clown_go(expression_genes$gene)%>%dotplot()
clown_go(interaction_genes$gene)%>%dotplot()


### kind of curious for the whole object
whole_matrix = obj@assays$RNA$data%>%as.matrix()%>%t()

obj$ieg_PC1 = module_pc1(iegs_)

DotPlot(obj, 'ieg_PC1')

ieg_grouped = obj@meta.data%>%
  subset(Status %in% c('M','D','F'))%>%
  group_by(final_clusters, individual, Status)%>%
  summarize(mean_ieg = mean(ieg_PC1))

ggplot(ieg_grouped, aes(x = final_clusters, y = mean_ieg, color =Status))+
  geom_boxplot()

## clusters of interest 
#> 5
#> 17?
#> 18?

ggplot(subset(ieg_grouped, final_clusters == 5), aes(x = Status, y = mean_ieg, color =Status))+
  geom_boxplot()+
  geom_point()

model_5 = lmer(ieg_PC1~Status+(1|individual), data = subset(obj@meta.data, final_clusters == 5 & Status %in% c('M','D','F')))
car::Anova(model_5, 3) # damn

model_6 = lmer(ieg_PC1~Status+(1|individual), data = subset(obj@meta.data, final_clusters == 6 & Status %in% c('M','D','F')))
car::Anova(model_6, 3) # damn

for(i in 0:26){
  model = lmer(ieg_PC1~Status+(1|individual), data = subset(obj@meta.data, final_clusters == i & Status %in% c('M','D','F')))
m = car::Anova(model, 3) 
print(i)
print(m$`Pr(>Chisq)`[2])
  
}
# i wonder if this may not be the right approach lol, maybe the binary is better after all

### cursory glance at binary approach #### 
mat <- sub_6@assays$RNA$data

# genes of interest
iegs_genes <- iegs_  # from calculateCoexpressedGenes

# subset to IEG genes that exist in the matrix
iegs_genes <- iegs_genes[iegs_genes %in% rownames(mat)]

# create a logical matrix: TRUE if gene expressed in cell
ieg_detected <- mat[iegs_genes, ] > 0

# count per cell how many IEGs are detected
sub_6$IEG_count <- colSums(ieg_detected)

# check
head(sub_6$IEG_count)
summary(sub_6$IEG_count)

### binary ====
sub_6$IEG_pos = as.numeric(sub_6$IEG_count>0)

ieg_prop_pos = sub_6@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarise(n_pos = sum(IEG_pos),
            n_cells = n())%>%
  mutate(prop_ieg = n_pos/n_cells)

ggplot(ieg_prop_pos, aes(x = Status, y = prop_ieg))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~sub.cluster)
# 6_3 again seems less active in females
ieg_prop_pos_6_3 = subset(ieg_prop_pos, sub.cluster=='6_3')

prop_glmer = glmer(cbind(ieg_prop_pos_6_3$n_pos, ieg_prop_pos_6_3$n_cells)~Status+(1|individual), family = binomial(), data = ieg_prop_pos_6_3)
car::Anova(prop_glmer, 3)

# will probably have to look at my beta binomial code to analyze some of these higher ones

###
sub_6_pc1 = sub_6@meta.data%>%
  subset( Status%in%c('M','D','F'))%>%
  group_by(individual, Status)%>%
  summarize(mean_pc1 = mean(ieg_PC1))
  
ggplot(sub_6_pc1, aes(x = Status, y = mean_pc1))+
  geom_boxplot()+
  geom_jitter()

sub_6_2_pc1 = sub_6@meta.data%>%
  subset(sub.cluster =='6_2'& Status%in%c('M','D','F'))%>%
  group_by(individual, Status)%>%
  summarize(mean_pc1 = mean(ieg_PC1))
  
ggplot(sub_6_2_pc1, aes(x = Status, y = mean_pc1))+
  geom_boxplot()+
  geom_jitter()



