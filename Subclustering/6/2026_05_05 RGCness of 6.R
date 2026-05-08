library(Seurat)
library(tidyverse)
library(AUCell)

RGC_markers = FindMarkers(obj, 1, group.by = 'final_clusters', only.pos = T)
RGC_markers_strong = head(rownames(RGC_markers), 100)
  
counts_mat = obj@assays$RNA$counts%>%as.matrix()
cell_rankings <- AUCell_buildRankings(counts_mat, plotStats = FALSE)

auc_scores <- AUCell_calcAUC(RGC_markers_strong, cell_rankings)

auc_mat <- t(auc_scores@assays@data[['AUC']]%>%as.matrix())  # now cells × cell types
colnames(auc_mat) <- paste0("AUC_", colnames(auc_mat))

obj <- AddMetaData(obj, metadata = as.data.frame(auc_mat))


DotPlot(obj, 'AUC_geneSet')

dat_6_rgclike = obj@meta.data%>%
  subset(final_clusters==6& Status != "NRM")%>%
  group_by(individual, Status)%>%
  summarize(mean_rgc =mean(AUC_geneSet))

ggplot(dat_6_rgclike, aes(x = Status, y = mean_rgc))+
  geom_boxplot()+
  geom_point()

#lets fucking gooooo


mod = lm(mean_rgc~Status, data = dat_6_rgclike)
anova(mod, test = 'Chisq')

pairs(emmeans(mod, 'Status'), adjust = 'none')

moder = lmer(AUC_geneSet~Status+(1|individual), data = obj@meta.data%>%
               subset(final_clusters == 6 & Status != 'NRM'))
car::Anova(moder, 3) # sad!
pairs(emmeans(moder, 'Status'), adjust = 'none')


# ok well this isnt significant but like look at it