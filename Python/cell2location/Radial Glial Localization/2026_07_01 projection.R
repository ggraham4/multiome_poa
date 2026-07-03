library(Seurat)
library(tidyverse)
library(Matrix)

# ── Paths ──────────────────────────────────────────────────────────────────
results_dir   <-"C:/Users/Gabe/Desktop/Visium/RGC/cell2location_map"
groupings_dir <- "C:/Users/Gabe/Desktop/multiome_poa/Visium/Groupings"
plot_dir      <- "C:/Users/Gabe/Desktop/Visium/plots"
dir.create(plot_dir, showWarnings = FALSE)

vis = readRDS("C:/Users/Gabe/Desktop/Visium/vis_better_barcodes_dissection_c2l05_projection_anatomical.rds")
#have to perform this subsetting to make it run
vis = subset(vis,seurat_clusters.projected0.4%in%c(1,26, 9) )

cat("Loading c2l high-confidence (q05) abundance...\n")
abundance <- read.csv(
  file.path(results_dir, "abundance_q05.csv"),
  row.names = 1,
  check.names = FALSE
)
cat("Abundance dims:", dim(abundance), "\n")
cat("Abundance barcode example:", head(rownames(abundance), 3), "\n")

# ── Determine Dominant Cell Types ──────────────────────────────────────────
cat("\nIdentifying dominant cell type per sketch spot...\n")

# Find the column index with the highest q05 value for each cell/row
max_idx <- max.col(abundance, ties.method = "first")
dominant_labels <- colnames(abundance)[max_idx]
names(dominant_labels) <- rownames(abundance)

# Add this categorical label to the Seurat object's metadata
vis$dominant_c2l_cluster <- NA
vis$dominant_c2l_cluster[names(dominant_labels)] <- dominant_labels

cat("Unique cell types assigned to sketch:", length(unique(na.omit(vis$dominant_c2l_cluster))), "\n")

# ── Project Labels from Sketch to Full Object ──────────────────────────────
cat("\nProjecting dominant classifications to full object...\n")

vis <- ProjectData(
  object             = vis,
  assay              = "Spatial.Polygons",
  full.reduction     = "full.pca.sketch",
  sketched.assay     = "sketch",          # Targets the original sketch assay
  sketched.reduction = "pca.sketch",
  dims               = 1:15,
  # Pass the metadata column string so Seurat maps the categorical predictions
  refdata            = list(c2l_predicted_id = "dominant_c2l_cluster")
)

cat("Projected metadata column 'c2l_predicted_id' added successfully.\n")

# Verify predictions populated
cat("Assigned spots in full object:", sum(!is.na(vis$c2l_predicted_id)), "\n")

# ── Replot with Projected Categories ───────────────────────────────────────
# Make sure available_fovs is defined or grab them from the object
if (!exists("available_fovs")) {
  available_fovs <- Images(vis)
}

for (fov in available_fovs) {
  fov_clean <- gsub("s_|\\.polygons", "", fov)
  fov_dir   <- file.path(plot_dir, paste0('RGC/',fov_clean, "_projected_clusters"))
  dir.create(fov_dir, showWarnings = FALSE)
  DefaultFOV(vis) <- fov
  
  cat("Plotting projected clusters for FOV:", fov, "\n")
  
  # Generate a discrete spatial plot of our projected classifications
  tryCatch({
    p <- ImageDimPlot(
      object = vis, 
      group.by = "c2l_predicted_id", 
      size = .75,
      cols = "alphabet" # Or provide a custom color vector if desired
    ) + ggtitle(paste0(fov_clean, " - Projected C2L Clusters"))
    
    outfile <- file.path(fov_dir, "projected_c2l_map.png")
    ggsave(outfile, p, width = 10, height = 8, dpi = 300)
    cat("  Saved cluster map to:", outfile, "\n")
  }, error = function(e) cat("  Failed plotting FOV map:", fov, "-", conditionMessage(e), "\n"))
}

cat("\nPipeline complete!\n")

#saveRDS(vis, file.path("C:/Users/Gabe/Desktop/Visium/RGC_vis_better_barcodes_dissection_c2l95_projection.rds"))
vis = readRDS( file.path("C:/Users/Gabe/Desktop/Visium/RGC_vis_better_barcodes_dissection_c2l95_projection.rds"))

vis$Predicted_RGC <- sub("^.*sf-", "", vis$c2l_predicted_id)


library(pheatmap)
heatmat = table(vis$Predicted_RGC,
                vis$Anatomical)%>%
  as.data.frame.matrix()%>%
  dplyr::select(!c('not sure 2',
                   'not sure 1',
                   'ventral diffuse',
                   'No Region',
                   'not sure 4',
                   'not sure 3',
                   'Diffuse',
                   'What'))

# Define anchor colors
my_colors <- c("gray",'white', "orange", "red")

# Create a function, then evaluate for n = 50 colors
color_func <- colorRampPalette(my_colors)
continuous_palette <- color_func(50)

p =pheatmap((heatmat/rowSums(heatmat))%>%scale(),
            cluster_rows = T,
            cluster_cols = T,
            treeheight_row=0,
            treeheight_col = 0,
            color =continuous_palette,
            angle_col =90,
            number_color= 'black')
p

#svg(filename = 'C:/Users/Gabe/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2.3/RGC heatmap.svg',
 #   width = 6.5,
#    height = 4.5)
#p
#dev.off()

spatial_figure_function = function(slice = '6P17', clusters = c()) {
  
  # Subset once and reuse
  vis_sub = subset(vis, !is.na(Predicted_RGC) & 
                     Predicted_RGC %in% clusters)
  
  image_name  = paste0('s_', slice, '.polygons')
  save_base   = paste0('Manuscript/Plots/Manuscript v1.2.3/', slice, 
                       '_projected_clusters/cluster_', 
                       paste(clusters, collapse = "_"))
  
  # Shared plot args
  plot_args = list(
    object   = vis_sub,
    group.by = "Predicted_RGC",
    images   = image_name,
    image.scale = 'hires',
    crop =F
  )
  
  # Tissue plot (with image)
  p = do.call(SpatialDimPlot, c(plot_args, image.alpha = 0.5))
  
  ggsave(paste0(save_base, '_tissue.jpg'), p, 
         width = 10, height = 8, dpi = 300, device = 'jpeg')
  cat("File size:", file.info(outfile)$size, "bytes\n")
  cat("Modified:", format(file.info(outfile)$mtime), "\n")
  Sys.sleep(1)  # let any AV/OneDrive lock release
  cat("Still says:", format(file.info(outfile)$mtime), "\n")
  print(paste0('saved to ',save_base ))
  # Points-only plot (no image)
  d = do.call(SpatialDimPlot, c(plot_args, image.alpha = 0))
  ggsave(paste0(save_base, '_points.jpg'), d, 
         width = 10, height = 8, dpi = 300, device = 'jpeg')
}
setwd('C:/Users/Gabe/Desktop/multiome_poa')
spatial_figure_function(clusters = unique(vis$Predicted_RGC))
spatial_figure_function('4P5',clusters = unique(vis$Predicted_RGC))
spatial_figure_function('3P5',clusters = unique(vis$Predicted_RGC))
spatial_figure_function('4P10',clusters = unique(vis$Predicted_RGC))


vis_sub = subset(vis, !is.na(Predicted_RGC) & 
                   Predicted_RGC %in% unique(vis$Predicted_RGC))

SpatialDimPlot( object   = vis_sub,
                group.by = "Predicted_RGC",
                images   = 's_6P17.polygons',
                image.scale = 'hires',
                crop =F)


