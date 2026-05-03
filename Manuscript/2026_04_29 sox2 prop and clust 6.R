obj$'sox2_pos'=ifelse(obj@assays$RNA$data['sox2',]>0, T, F)

d = FeaturePlot(obj, 
            'sox2_pos',
            reduction = 'harmony_wnn.umap',
            raster = T,
            pt.size = 2
            )+
  theme_void()+
  theme(legend.position = 'none')
d 

ggsave(plot = d,
       file = paste0('sox2_umap.svg'),
       device = "svg",
       units = "in",
       width =2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/"))

DotPlot(obj, 'sox2')

obj$clust_6=ifelse(obj$final_clusters==6, T, F)


six = FeaturePlot(obj, 
            'clust_6',
            reduction = 'harmony_wnn.umap',
            raster = T,
            pt.size = 2
            )+
  theme_void()+
  theme(legend.position = 'none',
        title = element_blank())
six

ggsave(plot = six,
       file = paste0('clust_6_umap.svg'),
       device = "svg",
       units = "in",
       width =2.33,
       height = 2.33,
       path = paste0("Manuscript/Plots/Manuscript v1.2/"))

obj$six_sox2 = ifelse(obj$clust_6 & obj$sox2_pos==T, T, F)
six_sox2_pos = FindAllMarkers(obj, 
            group.by = 'six_sox2', only.pos = T)

