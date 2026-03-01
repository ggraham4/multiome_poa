{
  library(scater)
  library(Seurat)
  library(condiments)
  library(emmeans)
  library(forcats)
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
  library(CytoTRACE)
  library(glmGamPoi)
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)
  library(BiocManager)
  library(slingshot)
  library(SingleCellExperiment)
  library(SeuratWrappers)
}

# read in data
obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds") 
# need this obj because it has the pca
FeaturePlot(obj, 'cckb',reduction = 'harmony_wnn.umap' )
FeaturePlot(obj, 'LOC111562384',reduction = 'harmony_wnn.umap' )
DimPlot(obj, reduction = 'harmony_wnn.umap')

#find rgc subclusters - this tells me the source node (1_1 is the most immature I know)
Idents(obj) = 'res0.8_50nn_40PC_45LSI'
obj <- FindSubCluster(obj, 1, 'harmony.wsnn')
DimPlot(obj, group.by = 'sub.cluster', reduction = 'harmony_wnn.umap')

# factorize status to have plots properly ordered
obj@meta.data$Status = factor(obj@meta.data$Status , levels = c('NRM',
                                                                              'M',
                                                                              'D',
                                                                              'E',
                                                                              'NF',
                                                                              "F"))


###subset to only radial glia and neurons
subset_obj <- subset(obj, 
                     #oligos
                     res0.8_50nn_40PC_45LSI!=2&
                     #microglia
                    res0.8_50nn_40PC_45LSI!=11&
                    #opcs
                    res0.8_50nn_40PC_45LSI!=13&
                    #dividing glia
                    res0.8_50nn_40PC_45LSI!=26&
                    #leuko
                    res0.8_50nn_40PC_45LSI!=20&
                    #ependymal
                    res0.8_50nn_40PC_45LSI!=15
                    
                    )
sce =readRDS( '/Users/ggraham/Desktop/slingshot.rds')

pseudotime_slingshot <- slingPseudotime(sce)
sce$Status = as.factor(subset_obj$Status)

#### Generalize to all lineages ####
num_lineages = ncol(pseudotime_slingshot)

# put all lineages in one place
num_lineages_pt_df = list()
for(i in 1:num_lineages){
  lin <- i
pt_df <- data.frame(
  cell = colnames(sce),
  pseudotime = pseudotime_slingshot[, lin],
  Status = sce$Status,
  individual = subset_obj$individual
) %>%
  filter(!is.na(pseudotime),
         Status %in% c("M","D","F"))

  pt_df$lineage = i
  num_lineages_pt_df[[i]] =  pt_df
}

### rolling window setup
window_width <- 10

all_eval_points = list() 
for(i in 1:length(num_lineages_pt_df)){

  minmax = num_lineages_pt_df[[i]]  

step_size <- (max(minmax$pseudotime) -  min(minmax$pseudotime)) / 50

eval_points <- seq(
  min(minmax$pseudotime),
  max(minmax$pseudotime),
  by = step_size
)

all_eval_points[[i]] = eval_points
}

### total cells per individual 
total_cells <- subset_obj@meta.data %>%
  filter(Status %in% c("M", "D", "F")) %>%
  group_by(individual, Status) %>%
  summarize(n_cells = n(), .groups = "drop")


# rolling function 
rolling_calc <- function(window, lineage){
  
  window_start <- window - (window_width/2)
  window_end   <- window + (window_width/2)
  
  lineage_data = num_lineages_pt_df[[lineage]]
  
  meta_window <- subset(
    lineage_data,
    pseudotime >= window_start & pseudotime <= window_end
  )
  
  cells_in_window <- meta_window %>%
    group_by(individual) %>%
    summarize(cells_in = n(), .groups = "drop")
  
  joint <- cells_in_window %>%
    right_join(total_cells, by = "individual") %>%
    mutate(
      cells_in = ifelse(is.na(cells_in), 0, cells_in),
      proportion_cells = cells_in / n_cells,
      window = window
    )
  
  return(joint)
}


big_boy <- list()
for(i in 1:length(num_lineages_pt_df)){
  points_to_evaluate_here <- all_eval_points[[i]]
  
  big_boy[[i]] <- list()  # initialize inner list for this lineage
  
  for(j in seq_along(points_to_evaluate_here)){
    point <- points_to_evaluate_here[j]
    big_boy[[i]][[j]] <- rolling_calc(point, i)
  }
}

big_boy_collapsed <- lapply(big_boy, bind_rows)

# now analysis time # 
rolling_stats_out = data.frame()
for(lineage in 1:length(big_boy)){
  print(lineage)
  
  relevant = big_boy_collapsed[[lineage]]
  
  for(point in all_eval_points[[lineage]]){
    
    relevant_eval_points = subset(relevant,window == point & Status %in% c("M","D","F"))
    
      model_mat <- cbind(
    relevant_eval_points$cells_in,
    relevant_eval_points$n_cells - relevant_eval_points$cells_in
  )
  
  model <- glm(model_mat ~ Status,
               data = relevant_eval_points,
               family = "binomial")
  
  av_ <- anova(model, test = "Chisq")
  
  pair <- pairs(
    emmeans(model, "Status"),
    adjust = "none"
  ) %>% as.data.frame()
  
  newd <- data.frame(
    lineage =lineage, 
    window = point,
    av_p = av_$`Pr(>Chi)`[2],
    m_d_p = pair$p.value[pair$contrast == "M - D"],
    m_f_p = pair$p.value[pair$contrast == "M - F"],
    d_f_p = pair$p.value[pair$contrast == "D - F"],
    m_d_estimate = pair$estimate[pair$contrast == "M - D"],
    m_f_estimate = pair$estimate[pair$contrast == "M - F"],
    d_f_estimate = pair$estimate[pair$contrast == "D - F"]
  )
  
  rolling_stats_out <- rbind(rolling_stats_out, newd)

  }
}

rolling_stats_out$av_q <- p.adjust(
  rolling_stats_out$av_p,
  method = "fdr",
  nrow(rolling_stats_out)
)

rolling_stats_out$signif <- ifelse(rolling_stats_out$av_q < 0.05, "*", NA)

table(rolling_stats_out$signif, rolling_stats_out$lineage)
# 2, 3, and 7 are interesting

lin_2=ggplot(big_boy_collapsed[[2]],
       aes(x = window, y = proportion_cells, color = Status)) +
  stat_summary(geom = 'line', fun = 'mean', aes(color = Status), size =2)+
  geom_point(alpha = 0.4) +
  geom_vline(
    data = subset(subset(rolling_stats_out, lineage == 2), signif == "*"),
    aes(xintercept = window),
    linetype = "dashed"
  )
# more dominants at earlier windows and fewer at later windows

lin_3 = ggplot(big_boy_collapsed[[3]],
       aes(x = window, y = proportion_cells, color = Status)) +
  stat_summary(geom = 'line', fun = 'mean', aes(color = Status), size =2)+
  geom_point(alpha = 0.4) +
  geom_vline(
    data = subset(subset(rolling_stats_out, lineage == 3), signif == "*"),
    aes(xintercept = window),
    linetype = "dashed"
  )
# not really sure what to make of this one

lin_7 = ggplot(big_boy_collapsed[[7]],
       aes(x = window, y = proportion_cells, color = Status)) +
  stat_summary(geom = 'line', fun = 'mean', aes(color = Status), size =2)+
  geom_point(alpha = 0.4) +
  geom_vline(
    data = subset(subset(rolling_stats_out, lineage == 7), signif == "*"),
    aes(xintercept = window),
    linetype = "dashed"
  )
# dominants enriched in this lineage at all points


#### some plotting ####
plot_sig_pseudotime <- function(lineage_num, rolling_stats_out, sce, window_width = 10) {
  
  half_window <- window_width / 2
  pt_col <- paste0("slingPseudotime_", lineage_num)
  pt_vals <- sce[[pt_col]]
  
  sig_windows <- rolling_stats_out$window[
    rolling_stats_out$signif == "*" & 
    rolling_stats_out$lineage == lineage_num
  ]
  
  in_sig_window <- rep(FALSE, ncol(sce))
  for(sig_win in sig_windows) {
    in_window <- !is.na(pt_vals) & pt_vals >= (sig_win - half_window) & pt_vals <= (sig_win + half_window)
    in_sig_window <- in_sig_window | in_window
  }
  
  sce$pseudotime_sig <- ifelse(in_sig_window, pt_vals, NA)
  plotUMAP(sce, colour_by = "pseudotime_sig") + ggtitle(paste("Lineage", lineage_num))
}

# usage
lin_2
plot_sig_pseudotime(2, rolling_stats_out, sce)
# woah this is interesante, more doms in the radial glia, fewer in these dorsal cells

lin_7
plot_sig_pseudotime(7, rolling_stats_out, sce)
# cluster 0
# again same enrichment in RGC and I think enrichment in 22? interesting 

lin_3
plot_sig_pseudotime(3, rolling_stats_out, sce)
DimPlot(subset_obj, reduction = 'harmony_wnn.umap')
#basically all of 0 and 24, like robustly 0 and 3
# what is the actual pairwise comparison in the later stages they all look the same to me

lin_3_stats = rolling_stats_out[rolling_stats_out$lineage==3 & !is.na(rolling_stats_out$signif),]
# i think a very small effect size 

lin_3+xlim(60, 75)
# dominants have more cells at the end which is 24 kinda

plot_sig_pseudotime(18, rolling_stats_out, sce)
plotUMAP(sce, color = 'slingPseudotime_18')



### which one includes 6 
plot_sig_pseudotime(9, rolling_stats_out, sce)


lin_9 = ggplot(big_boy_collapsed[[9]],
       aes(x = window, y = proportion_cells, color = Status)) +
  stat_summary(geom = 'line', fun = 'mean', aes(color = Status), size =2)+
  geom_point(alpha = 0.4) +
  geom_vline(
    data = subset(subset(rolling_stats_out, lineage == 9), signif == "*"),
    aes(xintercept = window),
    linetype = "dashed"
  )
lin_9

lin_9+xlim(65, 75)
# more in females, this tracks with the 6 population results

# i think it would be a very cool plot to be able to make a plot with this lineage
embedded_all <- embedCurves(sce, "UMAP")
embedded_all <-  slingCurves(embedded_all)

emb_all = list()
for(i in 1:length(embedded_all)){
  emb_all[[paste0(i)]] = (embedded_all)[[i]] 
}


  dat = emb_all[[paste0(9)]][["s"]]%>%as.data.frame()
  dat$curve_id <- paste0(9)   # evaluate NOW
  
p=plotUMAP(sce, color = 'slingPseudotime_9') +
    geom_path(
      data = dat,
      aes(x = harmonywnnUMAP_1,
          y = harmonywnnUMAP_2),
      size = 1.2
    )+
   theme_void()+
      theme(legend.position = 'top')

     #ggsave(plot = p,
    #   file = "2026_02_18_pseudotime_by_sex_lineage_9.tiff",
    #   device = "tiff",
    #   units = "in",
    #   width = 4,
    #   height = 4,
     #  path = "Manuscript/Plots/Fig.4")
     
     
lin_9 = ggplot(big_boy_collapsed[[9]],
       aes(x = window, y = proportion_cells, color = Status)) +
  stat_summary(geom = 'line', fun = 'mean', aes(color = Status), size =2)+
  geom_point(alpha = 0.4, size =0.5)+
  scale_x_continuous(breaks = seq(0, max((big_boy_collapsed[[9]]$window)), by =10))+
  theme_minimal()+
geom_label(
    data = subset(subset(rolling_stats_out, lineage == 9), signif == "*"),
    aes(x = window, label = '*', y = .1), inherit.aes = F
  )+
  theme(legend.position = 'top')

  #   ggsave(plot = lin_9,
   #    file = "2026_02_18_pseudotime_by_sex_lineage_9_proportions.svg",
    #   device = "svg",
     #  units = "in",
      # width = 3,
       #height = 3,
       #path = "Manuscript/Plots/Fig.4")

     
#### DEGs along the pseudotime ####
sce_subset = sce[,sce$Status %in% c('M','D','F')]

     library(TSCAN)
     
### gam ###
pseudo <- testPseudotime(sce_subset, pseudotime=sce_subset$slingPseudotime_9)
pseudo$SYMBOL <- rowData(sce_subset)$SYMBOL
pseudo_ordered =pseudo[order(pseudo$p.value),]

plotExpression(sce_subset, features='sox6',
    x="slingPseudotime_9", colour_by="Status")

plotExpression(sce_subset, features=c(head(rownames(pseudo_ordered))),
    x="slingPseudotime_9", colour_by="Status")


sorted <- pseudo[order(pseudo$p.value),]
up_genes <- sorted[sorted$logFC > 0,]
down_genes <-  sorted[sorted$logFC < 0,]
     
plotExpression(sce_subset, features=c(head(rownames(up_genes))),
    x="slingPseudotime_9", colour_by="Status")

plotExpression(sce_subset, features=c(head(rownames(down_genes))),
    x="slingPseudotime_9", colour_by="Status")

# like this is just entirely unconvincing

## lets look at the RGC lineage it should be all TFs ###
pseudo_18 <- testPseudotime(sce_subset, pseudotime=sce_subset$slingPseudotime_18)
pseudo_18$SYMBOL <- rowData(sce_subset)$SYMBOL
pseudo_18_ordered =pseudo_18[order(pseudo_18$p.value),]

# yeah I dont know man Im not really convinced, maybe I should do a more complex model like
# degs by Status along the pseudotime, though that could be very computationally expensive,
#though in theory testLinearModel could do it, but I think its not worth my time


### heatmap of lineages ####

cell_clusters <- sce$cluster
pt_matrix <- slingPseudotime(sce)
num_lineages <- ncol(pt_matrix)

### Order clusters by mean pseudotime


mean_pt_per_cluster <- data.frame(
  cell = rownames(pt_matrix),
  cluster = cell_clusters,
  mean_pt = rowMeans(pt_matrix, na.rm = TRUE)
) %>%
  group_by(cluster) %>%
  summarize(mean_pt = mean(mean_pt, na.rm = TRUE)) %>%
  arrange(mean_pt)

cluster_order <- mean_pt_per_cluster$cluster


##  Order lineages by mean pseudotime

mean_lineage_pt <- as.data.frame(pt_matrix) %>%
  mutate(cell = rownames(pt_matrix)) %>%
  pivot_longer(cols = starts_with("Lineage")) %>%
  group_by(name) %>%
  summarize(mean_pseudotime = mean(value, na.rm = TRUE)) %>%
  arrange(mean_pseudotime)

lineage_order <- mean_lineage_pt$name


lineage_cluster_mat <- matrix(
  0,
  nrow = num_lineages,
  ncol = length(unique(cell_clusters)),
  dimnames = list(
    colnames(pt_matrix),
    sort(unique(cell_clusters))
  )
)

for (lin in 1:num_lineages) {
  cells_in_lin <- !is.na(pt_matrix[, lin])
  clusters_hit <- cell_clusters[cells_in_lin]
  tab <- table(clusters_hit)
  lineage_cluster_mat[lin, names(tab)] <- as.numeric(tab)
}


### 4. Normalize to proportions


lineage_cluster_mat_norm <- lineage_cluster_mat / 
                            rowSums(lineage_cluster_mat)


row_means <- rowMeans(lineage_cluster_mat_norm)
col_means <- colMeans(lineage_cluster_mat_norm)

lineage_cluster_mat_norm <- cbind(
  lineage_cluster_mat_norm,
  Mean = row_means
)

lineage_cluster_mat_norm <- rbind(
  lineage_cluster_mat_norm,
  Mean = c(col_means, mean(row_means))
)


lineage_cluster_mat_norm <- 
  lineage_cluster_mat_norm[lineage_order, cluster_order, drop = FALSE]

# Re-add Mean row/column at bottom/right
lineage_cluster_mat_norm <- rbind(
  lineage_cluster_mat_norm,
  Mean = c(col_means[cluster_order], mean(row_means))
)

lineage_cluster_mat_norm <- cbind(
  lineage_cluster_mat_norm,
  Mean = c(row_means[lineage_order], mean(row_means))
)



heatmap = pheatmap(
  lineage_cluster_mat_norm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 45,
  display_numbers = round(lineage_cluster_mat_norm, 2),
  color = colorRampPalette(c("white", "blue"))(100),
  fontsize = 7
)
heatmap
     
#ggsave(plot = heatmap,
 #      file = "2026_02_18_pseudotime_by_sex lineage heatmap.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 6,
     #  height = 4,
      # path = "Manuscript/Plots/RGC supplementary")

mean_pseudotime = pt_matrix %>%
  as.data.frame() %>%
  mutate(cell = rownames(pt_matrix),
         cluster = sce$cluster) %>%
  pivot_longer(cols = starts_with('Lineage')) %>%
  group_by(name, cluster) %>%
  summarize(mean_pseudotime = mean(value, na.rm = T)) %>%
  pivot_wider(names_from = cluster, values_from = mean_pseudotime)

mean_pseudo_matrix = as.matrix(mean_pseudotime[,-1])
rownames(mean_pseudo_matrix) = mean_pseudotime$name

# Sort rows and columns by their means
mean_pseudo_matrix = mean_pseudo_matrix[order(rowMeans(mean_pseudo_matrix, na.rm = T)), 
                                        order(colMeans(mean_pseudo_matrix, na.rm = T))]

# Add mean row and column
mean_pseudo_matrix = cbind(mean_pseudo_matrix, Row_Mean = rowMeans(mean_pseudo_matrix, na.rm = T))
mean_pseudo_matrix = rbind(mean_pseudo_matrix, Col_Mean = colMeans(mean_pseudo_matrix, na.rm = T))

heatmap2 = pheatmap(
  mean_pseudo_matrix,
  na_col = "gray",
  cluster_cols = F,
  cluster_rows = F,
    fontsize = 7,
    display_numbers = round(mean_pseudo_matrix, 2),
  angle_col = 45
)
heatmap2

#ggsave(plot = heatmap2,
 #      file = "2026_02_18_pseudotime_by_sex lineage heatmap mean pseudo.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 6,
     #  height = 4,
      # path = "Manuscript/Plots/RGC supplementary")


