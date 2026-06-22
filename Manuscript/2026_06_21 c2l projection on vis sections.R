library(Seurat)

vis = readRDS("C:/Users/Gabe/Desktop/Visium/vis_better_barcodes_dissection_c2l05_projection_anatomical.rds")

setwd('C:/Users/Gabe/Desktop/multiome_poa')

spatial_figure_function = function(slice = '6P17', clusters = c()) {
  
  # Subset once and reuse
  vis_sub = subset(vis, !is.na(Predicted_final_clusters) & 
                     Predicted_final_clusters %in% clusters)
  
  image_name  = paste0('s_', slice, '.polygons')
  save_base   = paste0('Manuscript/Plots/Manuscript v1.2.3/', slice, 
                       '_projected_clusters/cluster_', 
                       paste(clusters, collapse = "_"))
  
  # Shared plot args
  plot_args = list(
    object   = vis_sub,
    group.by = "Predicted_final_clusters",
    images   = image_name,
    image.scale = 'hires'
  )
  
  # Tissue plot (with image)
  p = do.call(SpatialDimPlot, c(plot_args, image.alpha = 0.5))
  ggsave(paste0(save_base, '_tissue.jpg'), p, 
         width = 10, height = 8, dpi = 300, device = 'jpeg')
  
  # Points-only plot (no image)
  d = do.call(SpatialDimPlot, c(plot_args, image.alpha = 0))
  ggsave(paste0(save_base, '_points.jpg'), d, 
         width = 10, height = 8, dpi = 300, device = 'jpeg')
}

spatial_figure_function(clusters = c(1, 11, 14, 15, 3, 6))
spatial_figure_function('3P5',clusters = c(0,1,2,4,5,8,10,22,23,26))
spatial_figure_function('4P5',clusters = c(0,1,4,7,23,6))
spatial_figure_function('4P5',clusters = c(0,1,4,23,6))
for(i in 0:26){
  print(i)
  for(f in c('6P17','4P5','4P10','3P5')){
    print(f)
    p = SpatialDimPlot(
      object = subset(vis, !is.na(Predicted_final_clusters)&
                        Predicted_final_clusters==i), 
      group.by = "Predicted_final_clusters", 
      images= paste0('s_',f,'.polygons'),
      image.alpha = 0.5,
      image.scale = 'hires'
    )
    ggsave(paste0(
      'C:/Users/Gabe/Desktop/Visium/plots/',f,'_projected_clusters/cluster_',i,'.png'
    ), p, width = 10, height = 8, dpi = 300)
    
  }
}
