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

module_pc1 = function(module, object){
  whole_matrix = object@assays$RNA$data%>%as.matrix()%>%t()
  genes = colnames(whole_matrix)
  module_matrix = whole_matrix[,which(genes %in% module)]
  
  pc = prcomp(module_matrix, scale = T)
  scores <- pc$x[,1]

  if(mean(pc$rotation[c("fosb","egr1","npas4a"),1], na.rm=TRUE) < 0){
  scores <- -scores
}
return(scores)

}

sub_10 = FindSubCluster(obj, 10, graph.name='harmony.wsnn')
Idents(sub_10) <- 'sub.cluster'
sub_10 = subset(sub_10, final_clusters ==10)
sub_10$sub.cluster = factor(sub_10$sub.cluster, levels = c(paste0('10_',0:3)))
sub_10$Status = factor(sub_10$Status, levels = c('NRM','M',"D",'E','NF','F'))

DimPlot(sub_10)

marks_10 = FindAllMarkers(sub_10)

cyto = CytoTRACE(sub_10@assays$RNA$data%>%as.matrix())
sub_10$cyto = cyto$CytoTRACE

ggplot(sub_10@meta.data, aes(x = sub.cluster, color = Status, y = cyto))+
  geom_boxplot()

ggplot(sub_10@meta.data, aes(x = sub.cluster, y = cyto))+
  geom_boxplot()

iegs <- c('npas4a', 'fosb','egr1')

sub_10 = AddModuleScore(sub_10, list(iegs), name = 'ieg')
DotPlot(sub_10, 'ieg1', group.by = 'sub.cluster')+
  coord_flip()

ggplot(sub_10@meta.data, aes(x = sub.cluster, y = ieg1, color = Status))+
  geom_boxplot()

iegs2 =calculateCoexpressedGenes(sub_10, iegs)
sub_10 = AddModuleScore(sub_10, list(iegs2), name = 'iegs2')



module_plot <- sub_10@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(mean_ieg = mean(iegs21),
            se_ieg = sd(iegs21)/sqrt(n()))

ggplot(module_plot, aes(x = sub.cluster, y = mean_ieg, fill = Status, shape = Status))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
    geom_point(position = position_jitterdodge(1))

ggplot(module_plot, aes(x = Status, y = mean_ieg, fill = Status, shape = Status))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
    geom_point(position = position_jitterdodge(1))


iegs3 = module_pc1(iegs2, sub_10)
sub_10$iegs3 = iegs3


module_plot2 <- sub_10@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(mean_ieg = mean(iegs3),
            se_ieg = sd(iegs3)/sqrt(n()))

ggplot(subset(module_plot2, Status != 'NRM'), aes(x = sub.cluster, y = mean_ieg, fill = Status, shape = Status))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
    geom_point(position = position_jitterdodge(1))+
  labs(y = 'IEG_PC1')

ggplot(subset(module_plot2, Status != 'NRM'), aes(x = Status, y = mean_ieg, fill = Status, shape = Status))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
    geom_point(position = position_jitterdodge(1))

dat_10_0 = subset(sub_10@meta.data, sub.cluster == '10_0' & Status != 'NRM')

l = lmer(iegs3~Status + (1|individual), data = dat_10_0)
car::Anova(l, 3) #0.2

pairs(emmeans(l, 'Status'), adjust ='none')

DotPlot(sub_10, c('kiss1',
                  'slc18a3b'), 
        group.by = 'sub.cluster')+
  coord_flip()

DimPlot(sub_10, reduction = 'harmony_wnn.umap', label = T)

sub_10$kiss = ifelse(sub_10@assays$RNA$data['kiss1',]>0, T, F)
sub_10$Ach = ifelse(sub_10@assays$RNA$data['slc18a3b',]>0, T, F)


module_plot3 <- sub_10@meta.data%>%
  subset(sub.cluster== '10_0')%>%
  group_by(individual, Status, kiss, Ach)%>%
  summarize(mean_ieg = mean(iegs3),
            se_ieg = sd(iegs3)/sqrt(n()))

module_plot3$NT = ifelse(module_plot3$kiss == T & module_plot3$Ach == T, 'both', NA)
module_plot3$NT = ifelse(module_plot3$kiss == T & module_plot3$Ach == F, 'kiss_only', module_plot3$NT )
module_plot3$NT = ifelse(module_plot3$kiss == F & module_plot3$Ach == T, 'ach_only', module_plot3$NT )

module_plot3$NT  = factor(module_plot3$NT, levels  = c('ach_only', 'kiss_only', 'both'))

ggplot(subset(module_plot3, Status != 'NRM' & !is.na(NT)), aes(x = Status, y = mean_ieg, fill = Status, shape = Status))+
  geom_boxplot(alpha = 0.5, outlier.shape = NA)+
    geom_point(position = position_jitterdodge(1))+
  facet_wrap(~NT)+
  labs(y = 'IEG PC1')
