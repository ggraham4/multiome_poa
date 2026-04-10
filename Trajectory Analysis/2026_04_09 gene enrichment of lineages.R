# analyzing which genes are enricheed along my lineages of interest
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
  library(slingshot)
}

clown_go = readRDS('Functions/clown_go2')

#### read in data #### 
obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds") 

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
# read in precomputed SCE object from 2026_01_17 slingshot by sex.R
sce =readRDS( '/Users/ggraham/Desktop/slingshot.rds')
pseudotime_slingshot <- slingPseudotime(sce)
sce$Status = as.factor(subset_obj$Status)

### some analysis ####
### extract cells that are positive for my lineages
# so cluster 0 has several lineages going through it, so lets look at cluster 6 
# for 0, I might do lineage 7

# put all lineages in one place
num_lineages = ncol(pseudotime_slingshot)
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

# add lineage 8 cells to object
cells_lin_8 = num_lineages_pt_df[[8]]$cell
subset_obj$lineage_8_pos = ifelse(rownames(subset_obj@meta.data)%in%cells_lin_8, T, F)
subset_obj$sox2_pos = ifelse(subset_obj@assays$RNA$data['sox2',]>0, T, F)

FeaturePlot(subset_obj, 'lineage_8_pos', reduction = 'harmony_wnn.umap')
FeaturePlot(subset_obj, 'sox2_pos', reduction = 'harmony_wnn.umap')

table(subset_obj@meta.data$sox2_pos,subset_obj@meta.data$lineage_8_pos )
table(subset_obj@meta.data$sox2_pos,subset_obj@meta.data$lineage_8_pos )%>%fisher.test()
# insanely significant, I still think bootstrapping is better

# so we could also do this with a bootstrap, its also worth considering it by individual 
# though I suspect it will still be signif, also we could do find markers to find 
# enriched genes I realize

enriched_lin_8 = FindAllMarkers(subset_obj, group.by = 'lineage_8_pos', only.pos = T)

# ok lets do by individual 
dat_sox_ind = subset_obj@meta.data%>%
  group_by(individual, Status, lineage_8_pos)%>%
  summarize(mean_sox = mean(sox2_pos))

ggplot(dat_sox_ind, aes(x = lineage_8_pos, y=mean_sox ))+
  geom_boxplot()+
  geom_point()
# yeah its not close

dat_sox_ind$Status = factor(dat_sox_ind$Status, levels = c('NRM', 'M','D','E','NF','F'))
ggplot(dat_sox_ind, aes(x = lineage_8_pos, y=mean_sox, color = Status))+
  geom_boxplot()+
  geom_point(position = position_dodge(0.75))+
  labs(y = 'Prop Sox2+', x = 'Lineage 8 Status')

##### ok lets do a bootstrap anyway #####
observed_table = table(subset_obj@meta.data$sox2_pos,
                       subset_obj@meta.data$lineage_8_pos)%>%as.data.frame.matrix()

observed_data = observed_table[2,2]/(observed_table[1,2]+observed_table[2,2])
sample_size = sum(observed_table[,2])

n_perm = 10000
out_list = c()
for(i in 1:n_perm){
  print(i)
  set.seed(i)
  cells = sample(rownames(subset_obj@meta.data), sample_size)
  
  temp_data = subset_obj@meta.data[rownames(subset_obj@meta.data)%in%cells,]
  
  sox2_count = sum(temp_data$sox2_pos)
  
  sox2_prop = sox2_count/sample_size
  
  out_list = append(out_list, sox2_prop)
}

sum(out_list>=observed_data)/n_perm

hist(out_list)
# permutation analysis also agrees, p value = 0

sum(subset_obj$res0.8_50nn_40PC_45LSI==1 & subset_obj$sox2_pos==T)/sum(subset_obj$res0.8_50nn_40PC_45LSI==1 )
# only 18% of RGC are sox2+ so I think this effect is not just driven by RGCs

# need to do this for other lineages I think that is the true null distribution ####
# function to get per-individual sox2 proportions for a given lineage
get_sox2_by_lineage <- function(lineage_num) {
  cells_lin = num_lineages_pt_df[[lineage_num]]$cell
  subset_obj$lineage_pos = ifelse(rownames(subset_obj@meta.data) %in% cells_lin, T, F)
  
  dat = subset_obj@meta.data %>%
    group_by(individual, lineage_pos) %>%
    summarize(mean_sox = mean(sox2_pos), .groups = 'drop') %>%
    filter(lineage_pos == TRUE) %>%
    mutate(lineage = lineage_num)
  
  return(dat)
}

# function to get fisher p value for a given lineage
get_fisher_p <- function(lineage_num) {
  cells_lin = num_lineages_pt_df[[lineage_num]]$cell
  subset_obj$lineage_pos = ifelse(rownames(subset_obj@meta.data) %in% cells_lin, T, F)
  
  p = table(subset_obj@meta.data$sox2_pos, subset_obj@meta.data$lineage_pos) %>%
    fisher.test() %>%
    .$p.value
  
  return(data.frame(lineage = lineage_num, p_value = p))
}

# run across all lineages
all_sox2 = lapply(1:num_lineages, get_sox2_by_lineage) %>% bind_rows()
all_fisher = lapply(1:num_lineages, get_fisher_p) %>% bind_rows()

# plot
ggplot(all_sox2, aes(x = factor(lineage), y = mean_sox)) +
  geom_boxplot() +
  geom_point() +
  labs(x = 'Lineage', y = 'Prop Sox2+')

### correct for sampe size using logistic regression ##
get_corrected_p <- function(lineage_num) {
  cells_lin = num_lineages_pt_df[[lineage_num]]$cell
  
  dat = subset_obj@meta.data %>%
    mutate(
      lineage_pos = rownames(subset_obj@meta.data) %in% cells_lin,
      res0.8_50nn_40PC_45LSI = factor(res0.8_50nn_40PC_45LSI)
    )
  
  fit = glm(sox2_pos ~ lineage_pos , # before, I had a covariate of cluster ID here, should that be the case?
            data = dat,
            family = binomial)
  
  z = summary(fit)$coefficients['lineage_posTRUE', 'z value']
  p = summary(fit)$coefficients['lineage_posTRUE', 'Pr(>|z|)']
  
  return(data.frame(lineage = lineage_num, z_value = z,
                    p_value = p))
}

all_corrected = lapply(1:num_lineages, get_corrected_p) %>% bind_rows()

mecp = readRDS('Functions/mean_expression_cluster_plot.rds')
mecp(subset_obj, 'LOC111574271', 6, 'res0.8_50nn_40PC_45LSI')
mecp(subset_obj, 'LOC111574271', 9, 'res0.8_50nn_40PC_45LSI')
mecp(subset_obj, 'pcna', 6, 'res0.8_50nn_40PC_45LSI')

### is 8 actually the most enriched for 6 ####

big_data = data.frame()
for(i in 1:num_lineages){
      cells_lin = num_lineages_pt_df[[i]]
      big_data = rbind(big_data, cells_lin)
}

subset_obj@meta.data$cell = rownames(subset_obj@meta.data)

big_joint = big_data%>%
  right_join(subset_obj@meta.data, by = 'cell')

cells_per_cluster =big_joint%>%
  group_by(res0.8_50nn_40PC_45LSI, lineage)%>%
  summarize(n = n())%>%
  na.omit()

total_cells =  big_joint%>%
  group_by(res0.8_50nn_40PC_45LSI)%>%
  summarize(n = n())%>%
  na.omit()

cells_per_joint = cells_per_cluster%>%right_join(total_cells, by = 'res0.8_50nn_40PC_45LSI')%>%
  mutate(prop_cells = n.x/n.y)
cells_per_joint$prop_cells = round(cells_per_joint$prop_cells, 4)

# i guess I was looking at the wrong lineage the whole time???

### LINEAGE 9 ####
# add lineage 8 cells to object
cells_lin_9 = num_lineages_pt_df[[9]]$cell
subset_obj$lineage_9_pos = ifelse(rownames(subset_obj@meta.data)%in%cells_lin_9, T, F)
subset_obj$sox2_pos = ifelse(subset_obj@assays$RNA$data['sox2',]>0, T, F)

FeaturePlot(subset_obj, 'lineage_9_pos', reduction = 'harmony_wnn.umap')
FeaturePlot(subset_obj, 'lineage_9_pos', reduction = 'rnaPCA')

FeaturePlot(subset_obj, 'sox2_pos', reduction = 'harmony_wnn.umap')
Marks_9 = FindAllMarkers(subset_obj, group.by = 'lineage_9_pos', only.pos = T)

Marks_rgc =  FindMarkers(subset_obj, 1, group.by = 'res0.8_50nn_40PC_45LSI', only.pos = T)
Marks_rgc$gene = rownames(Marks_rgc)

Marks_6=  FindMarkers(subset_obj, 6, group.by = 'res0.8_50nn_40PC_45LSI', only.pos = T)
Marks_6$gene = rownames(Marks_6)

Marks_clust_9=  FindMarkers(subset_obj, 9, group.by = 'res0.8_50nn_40PC_45LSI', only.pos = T)
Marks_clust_9$gene = rownames(Marks_clust_9)

# elements in 1st object that are not in second
genes_notrgc = setdiff(Marks_9$gene, Marks_rgc$gene )
genes_not_6 = setdiff(Marks_9$gene, Marks_6$gene )
genes_not_1_6 = intersect(genes_notrgc, genes_not_6)
genes_not_1_6_9 = setdiff(genes_not_1_6,Marks_clust_9$gene )

clown_go(genes_not_1_6_9)%>%dotplot()


# its not goin well for ol gabo ###