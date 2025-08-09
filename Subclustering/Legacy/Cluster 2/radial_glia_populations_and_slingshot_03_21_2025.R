library(here)

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
  library(emmeans)
  library(CytoTRACE)
  library(ggrepel)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(Polychrome)
  library(scCustomize)
  
  P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
  swatch(P40)
  names(P40) <- NULL
  
  mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
  prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
  define_degs_prop<- readRDS('Functions/define_degs_prop.rds')
  mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
  prop_deg_function.rds<- readRDS('Functions/DEG_functions/prop_deg_function.rds')
  define_behavior_degs<- readRDS('Functions/define_behavior_degs')
  clown_go<- readRDS('Functions/clown_go')
  define_degs<- readRDS('Functions/define_degs')
  
}

obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")
radial_glia <- subset(obj, harmony.wnn_res0.4_clusters==2)
radial_glia <- FindSubCluster(radial_glia, cluster = 2, subcluster.name = 'sub', graph.name = 'harmony.wsnn')

glia_markers <- readRDS('Subclustering/Cluster 2/glia subcluster markers.rds')

Idents(radial_glia)<- 'sub'

##### Marker Genes #####
glia_markers

overexpressed_genes <- list()
for(i in unique(glia_markers$cluster)){
  overexpressed_genes[[i]]$genes=   glia_markers$gene[glia_markers$p_val_adj <0.05 & 
                                                        (glia_markers$pct.1>glia_markers$pct.2)&
                                                        glia_markers$cluster==i]
  
}

#clown_go(overexpressed_genes$`2_0`$genes)%>%dotplot()
clown_go(overexpressed_genes$`2_1`$genes)%>%dotplot() 
#axon guidance, ephrins, rna pol II

#clown_go(overexpressed_genes$`2_2`$genes)%>%dotplot()
clown_go(overexpressed_genes$`2_3`$genes)%>%dotplot()
#cell adhesion, NS development, cell development

clown_go(overexpressed_genes$`2_4`$genes)%>%dotplot()
#wnt, notch, neural crest

clown_go(overexpressed_genes$`2_5`$genes)%>%dotplot()
#system dev, cell differentiation

clown_go(overexpressed_genes$`2_6`$genes)%>%dotplot()
#angiogenesis, hindbrain dev, signal trans

clown_go(overexpressed_genes$`2_7`$genes)%>%dotplot()
# wnt and hindbrain again

clown_go(overexpressed_genes$`2_8`$genes)%>%dotplot()

#clown_go(overexpressed_genes$`2_9`$genes)%>%dotplot()
#really?

##### Lineages #####
# I want to do slingshot on this especially
DimPlot(radial_glia, reduction = 'harmony_wnn.umap', label = T)

library(slingshot)

sce= SingleCellExperiment(assays = list(
  counts = GetAssayData(radial_glia, slot='counts'),
  logcounts = GetAssayData(radial_glia, slot='data')
))

reducedDims(sce)= list(
  PCA = Embeddings(radial_glia),
  UMAP = Embeddings(radial_glia, 'harmony_wnn.umap')
)

colData(sce)$cell_type = radial_glia$sub
colData(sce)$cluster = radial_glia$sub

sce <- slingshot(sce, 
                 clusterLabels = 'cluster',
                 reducedDim = 'UMAP',
                 start.clus = '2_9')

pseudotime_slingshot <- slingPseudotime(sce)
avg_pseudotime <- rowMeans(pseudotime_slingshot, na.rm = T)
avg_pseudotime[is.na(avg_pseudotime)] <- min(avg_pseudotime, na.rm =T)
#avg_pseudotime_norm <- 100* (avg_pseudotime - min(avg_pseudotime))/
#  (max(avg_pseudotime)- min(avg_pseudotime))

radial_glia$slingshot_pseudotime <- avg_pseudotime

curves = slingCurves(sce)
umap_coords <- reducedDims(sce)$UMAP

arr <- list(x = -3.5, y = 4, x_len = 2, y_len = 2)

std_umap <-DimPlot(radial_glia, reduction = 'harmony_wnn.umap', label=T, repel = T,   label.size = 4)+
  xlim(-5,5)+
  ylim(3,10)+
  theme_void()+
  theme(legend.position = 'none')+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), 
           arrow = arrow(type = "closed", length = unit(10, 'pt'))) +
  theme(legend.position = 'none')+
  annotate('text',
           x = -4, y = 5, label = 'UMAP_2', angle = 90, size =3)+
  annotate('text',
           x = -2.5, y = 3.5, label = 'UMAP_1', size =3)


std_umap

#ggsave(plot = std_umap,
#       file = "dimplot_2.svg",
#       device = "svg",
#       units = "in",
#       width = 3,
#       height = 3,
#       path = "Bachelors Thesis/Plots/Figure 4")

plot_data <- data.frame(UMAP1 =umap_coords[,1],
                        UMAP2 = umap_coords[,2],
                        Pseudotime = radial_glia$slingshot_pseudotime,
                        CellType = radial_glia$sub)

slingshot_plot <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Pseudotime))+
  geom_point(size =0.5)+
  scale_color_viridis_c()+
  theme_void()+
  xlim(-5,3)+
  ylim(3,10)
#add curves
for(i in seq_along(curves)){
  curve_data <- data.frame(UMAP1 = curves[[i]]$s[,1],
                           UMAP2 = curves[[i]]$s[,2])
  slingshot_plot <- slingshot_plot+ 
    geom_path(data = curve_data, aes(x = UMAP1, y = UMAP2),
              size = 1, color = 'black', alpha = 0.8)
  
}


arr2 <- list(x = -3.5, y = 4, x_len = 1, y_len = 2)
slingshot_plot2 <- slingshot_plot+
  annotate('text',
           x = -4, y = 5, label = 'UMAP_2', angle = 90, size =3)+
  annotate('text',
           x = -4, y = 3.5, label = 'UMAP_1', size =3)+
  theme_void()+
  annotate("segment", 
           x = arr2$x, xend = arr2$x + c(arr2$x_len, 0), 
           y = arr2$y, yend = arr2$y + c(0, arr2$y_len), 
           arrow = arrow(type = "closed", length = unit(10, 'pt'))) 
slingshot_plot2

#ggsave(plot = slingshot_plot2,
#       file = "slingshot_2.svg",
#       device = "svg",
#       units = "in",
#       width = 3.3,
#       height = 3,
#       path = "Bachelors Thesis/Plots/Figure 5")




sling_score_data <- data.frame(individual = radial_glia$individual,
                              Status = radial_glia$Status, 
                              cluster = radial_glia$sub,
                              pseudotime = radial_glia$slingshot_pseudotime)
library(forcats)

plot <- ggplot(sling_score_data, aes(x = fct_reorder(cluster, pseudotime), y = pseudotime))+
  geom_point(size = 1, shape = 1, alpha = 0.5, position = position_jitter(width = 0))+
  geom_boxplot(fill = NA, outlier.shape = NA)+
  theme_minimal()+
  labs(x = 'Subcluters',y = 'Pseudotime')+
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size =12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = -90))+
  ylim(20,27.5)
  
plot

#ggsave(plot = plot,
#       file = "slingshot_box.svg",
#       device = "svg",
#       units = "in",
#       width = 3,
#       height = 1.4,
#       path = "Bachelors Thesis/Plots/Figure 5")


cluster_2_matrix <- as.matrix(radial_glia@assays$RNA$counts) #not necessary to normalize
cluster_2_cyto <- CytoTRACE(mat = cluster_2_matrix
)
radial_glia$cyto <-cluster_2_cyto$CytoTRACE

cluster_2_data <- data.frame(individual = radial_glia@meta.data$individual,
                              status = radial_glia@meta.data$Status,
                              cluster = radial_glia@meta.data$sub,
                              cyto = radial_glia@meta.data$cyto)%>%
  subset(status != 'NRM')
cluster_2_data$status <- factor(cluster_2_data$status, levels = c('M','D',"E",'NF','F'))

cluster_cyto = ggplot(cluster_2_data, aes(x = fct_reorder(cluster, cyto, .desc = T), y = cyto))+
  geom_boxplot(fill = NA, outlier.shape = NA)+
  theme_minimal()+
  labs(x = 'Subcluster', y = 'CytoTRACE')+
  theme(axis.text.x = element_text(angle= -90),
        axis.text = element_text(size = 10),
        axis.title = element_text(size =12))+
  geom_point(size = 1, shape = 1, alpha = 0.2, position = position_jitter(width = 0))
  
cluster_cyto

#ggsave(plot = cluster_cyto,
#       file = "cluster_cyto.svg",
#       device = "svg",
#       units = "in",
#       width = 3,
#       height = 1.4,
#       path = "Bachelors Thesis/Plots/Figure 5")


###looking at manually curated markers
#https://www.nature.com/articles/s41593-020-00794-1#Sec11

ugomma_markers <- read_excel('A:/anemonefish_multiome_radial_glia_03_20_2025/markers_ugomma_2021.xlsx',9)

sapiens_to_ocellaris <- read.csv('C:/Users/Gabe/Desktop/multiome_poa/Reference/hsapiens_to_aocellaris.csv')

ugomma_markers_translated <- ugomma_markers%>%
  right_join(sapiens_to_ocellaris, by  = join_by('gene'=='hsapiens_name'))%>%
  filter(!is.na(p_val))

module <- list()
for(i in unique(ugomma_markers_translated$cluster)){
  b = 1

  marker_genes <- ugomma_markers_translated$aocellaris_name[ugomma_markers_translated$cluster %in% i]
  marker_gene_final <- unlist(marker_genes)
  marker_gene_final <- unique(marker_gene_final[marker_gene_final != ''])
  
  module[[paste0(i)]] <- c(marker_gene_final) #this keeps overwriting
  
  b = b+1
}


radial_glia <- AddModuleScore(
  object = radial_glia,
  features = module,
  name = 'new_modules'
  
)

DotPlot(radial_glia, 
        features = c(paste0(rep('new_modules',10),1:10)),
        dot.min = .5,
        col.min =  .5*2.5)+
  coord_flip()+
  scale_x_discrete(labels = names(module))+
  scale_color_gradientn(colors = c('white','red','blue'))

FeaturePlot(radial_glia, 'ntrk3b', reduction = 'harmony_wnn.umap' )
FeaturePlot(radial_glia, 'ntrk3a', reduction = 'harmony_wnn.umap' )

DotPlot(radial_glia, 
        features = c('lum',
                     'alx1',
                     #'sox2',
                     'top2a',
                     'dlx6a',
                     'rbm47',
                     'nog1',
                     'elavl3'))+
coord_flip()
  
#sox2 marks progenitors so 2_8 and 2_9 are likely progenitors; this corresponds with the module score
radial_glia$sub = factor(radial_glia$sub, levels = c('2_0',
                                                     '2_1',
                                                     '2_2',
                                                     '2_3',
                                                     '2_4',
                                                     '2_5',
                                                     '2_6',
                                                     '2_7',
                                                     '2_8',
                                                     '2_9'))
radial_glia_markers <- DotPlot(radial_glia, 
        group.by = 'sub',
        features = c(
                     'elavl3',
                     'ntrk3a',
                     'sox2',
                     'gfap',
                     's100b',
                     'LOC111577263', # 2_3 seems like it expresses a lot of neuron like markers suggesting to me it is a neuronal lineag
                     'kcnj10a' #tripartite synapse marker
                     ),
        dot.min = 0.1)+
  coord_flip()+
  theme(axis.text.x = element_text(angle=-90, vjust =-0))+
  scale_x_discrete(labels =c('elavl3',
                     'ntrk3a',
                     'sox2',
                     'gfap',
                     's100b',
                     'cyp19a1b',
                     'kcnj10a'))+
  theme(axis.title.y = element_blank())+
  labs(y = 'Subcluster')+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size=10),
        legend.title  = element_text(size=8),
        legend.position = 'right')
radial_glia_markers

ggsave(plot = radial_glia_markers,
       file = "radial_glia_markers.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Bachelors Thesis/Plots/Figure 5")



###########
DotPlot(radial_glia, features = c('LOC111563038' # pankx1
                                  ))+
  coord_flip()

FeaturePlot(radial_glia, 'elavl3', reduction = 'harmony_wnn.umap' , label = T, repel = T)
