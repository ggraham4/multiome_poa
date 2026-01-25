#> sex differences in pseudotime
#> here, I want to see if 
#> 1) The sexes differ in cell density along the pseudotime (do they get hung up)
#> 2) do they have differential expression along the pseudotime
#> 3) Are there sex specific branches in the pseudotime
#> 
#> here, I am using slingshot instead of monocle
#remove.packages(c('SeuratObject',
#                  'Seurat'))
#install.packages('SeuratObject')
#packageVersion('SeuratObject')
#install.packages('Seurat')


{
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


### Learn the trajectory ####
sce= SingleCellExperiment(assays = list(
  counts = GetAssayData(subset_obj, layer='counts'),
  logcounts = GetAssayData(subset_obj, layer='data')
))

subset_obj = FindVariableFeatures(subset_obj, layer  = "data")
subset_obj = RunPCA(subset_obj, dim = 50, verbose = TRUE, assay = "RNA",               
 features = VariableFeatures(object = subset_obj), reduction.name = "pca_slingshot",
 reduction.key = "pca_slingshot_")

reducedDims(sce)= list(
  PCA = Embeddings(subset_obj,'pca_slingshot' ),
  UMAP = Embeddings(subset_obj, 'harmony_wnn.umap')
)
#add these reductions into the cds object
reducedDims(sce)$PCA <- Embeddings(subset_obj, "pca_slingshot")

colData(sce)$cell_type = subset_obj$sub.cluster
colData(sce)$cluster = subset_obj$sub.cluster

set.seed(0)
sce <- slingshot(sce, 
                 clusterLabels = 'cluster',
                 reducedDim = 'PCA',
                 start.clus = '1_1')

pseudotime_slingshot <- slingPseudotime(sce)
#what to do with this matrix

sce$Status = as.factor(subset_obj$Status)
prog_results <- progressionTest(sce, conditions = sce$Status)

lin <- 1
pt_df <- data.frame(
  cell = colnames(sce),
  pseudotime = pseudotime_slingshot[, lin],
  Status = sce$Status
) %>%
  filter(!is.na(pseudotime),
         Status %in% c("M","D","F"))

# Plot
ggplot(pt_df, aes(x = pseudotime, fill = Status)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  facet_wrap(~Status) # looks like there is a missing male peak at 40 that is present in doms and females
# and a male and dominant peak at 20 that is missing in females

#### statistics ####
### rolling window setup
window_width <- 10

step_size <- (max(pt_df$pseudotime) -  min(pt_df$pseudotime)) / 50

eval_points <- seq(
  min(pt_df$pseudotime),
  max(pt_df$pseudotime),
  by = step_size
)

meta <- pt_df
meta$individual <- subset_obj$individual[match(meta$cell, colnames(subset_obj))]

### total cells per individual
total_cells <- meta %>%
  group_by(individual, Status) %>%
  summarize(n_cells = n(), .groups = "drop")

### rolling calculation
rolling_calc <- function(window){
  
  window_start <- window - (window_width/2)
  window_end   <- window + (window_width/2)
  
  meta_window <- subset(
    meta,
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

rolling_proportion <- lapply(eval_points, rolling_calc)
rolling_proportion_bound <- do.call(rbind, rolling_proportion)

ggplot(rolling_proportion_bound,
       aes(x = window, y = proportion_cells, color = Status)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE)

rolling_stats_out <- data.frame()

for(i in eval_points){
  subset_eval <- subset(
    rolling_proportion_bound,
    window == i & Status %in% c("M","D","F")
  )
  
  model_mat <- cbind(
    subset_eval$cells_in,
    subset_eval$n_cells - subset_eval$cells_in
  )
  
  model <- glm(model_mat ~ Status,
               data = subset_eval,
               family = "binomial")
  
  av_ <- anova(model, test = "Chisq")
  
  pair <- pairs(
    emmeans(model, "Status"),
    adjust = "none"
  ) %>% as.data.frame()
  
  newd <- data.frame(
    window = i,
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

rolling_stats_out$av_q <- p.adjust(
  rolling_stats_out$av_p,
  method = "fdr",
  nrow(rolling_stats_out)
)

rolling_stats_out$signif <- ifelse(rolling_stats_out$av_q < 0.05, "*", NA)

ggplot(rolling_proportion_bound,
       aes(x = window, y = proportion_cells, color = Status)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE) +
  geom_vline(
    data = subset(rolling_stats_out, signif == "*"),
    aes(xintercept = window),
    linetype = "dashed"
  )

# need to do this for several lineages ...

#saveRDS(sce, '/Users/ggraham/Desktop/slingshot.rds')
sce =readRDS( '/Users/ggraham/Desktop/slingshot.rds')
library(patchwork)
plot(as.SlingshotDataSet(sce))

embedded <- embedCurves(sce, "UMAP")
em <- slingCurves(embedded)[[1]] # only 1 path.
em <- as.data.frame(embedded$s[embedded$ord,])

DimPlot(subset_obj)+
    geom_path(data=em, aes(x=UMAP1, y=UMAP2))
