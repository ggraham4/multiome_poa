library(Seurat)
library(tidyverse)
library(Matrix)
library(pheatmap)

# ── Paths ──────────────────────────────────────────────────────────────────
results_dir   <- "C:/Users/Gabe/Desktop/Visium/POA/cell2location_map"
groupings_dir <- "C:/Users/Gabe/Desktop/multiome_poa/Visium/Groupings"
plot_dir      <- "C:/Users/Gabe/Desktop/Visium/plots"
dir.create(plot_dir, showWarnings = FALSE)

parker_2024 <- readRDS("A:/Anemonefish POA Legacy R Objects/parker_2024_realigned_ocellaris.rds")
parker_2024 <- subset(parker_2024, clusters49 != '35')

vis <- readRDS("C:/Users/Gabe/Desktop/Visium/vis_better_barcodes_c2l05_projection_parker_2024_clusters49.rds")
vis$Predicted <- sub("^.*sf-", "", vis$c2l_predicted_id)

images <- unique(vis$Slice)
barcodes_included <- data.frame()
for (i in images) {
  print(i)
  data <- read.csv(paste0("C:/Users/Gabe/Desktop/multiome_poa/Visium/Groupings/Maruska Regions Data Informed/2026_05_17 ",
                          tolower(i), " Anatomical.csv")) %>%
    subset(!is.na(Anatomical))
  data$Barcode_updated <- paste0(data$Barcode, i)
  barcodes_included <- rbind(data, barcodes_included)
}

vis$Anatomical <- "No Region"
idx <- match(colnames(vis), barcodes_included$Barcode_updated)
vis$Anatomical[!is.na(idx)] <- barcodes_included$Anatomical[idx[!is.na(idx)]]

renamed <- c(
  'could be nPPa'                    = 'What',
  'could be nMMp'                    = 'nMMp',
  'could be nPMP'                    = 'nPMp',
  'dorsolateral nucleus'              = 'Dln-r',
  'dorsolateral nucleus caudal part'  = 'Dln-c',
  'Dl-v2 dorsal part'                 = 'Dl-v2d',
  'TGn'                               = 'TGN'
)
renaming_col <- ifelse(vis$Anatomical %in% names(renamed),
                       renamed[vis$Anatomical],
                       vis$Anatomical)
vis$Anatomical <- renaming_col

vis$Anatomical <- ifelse(vis$Anatomical %in% c('not sure 21', 'OT'),
                         NA,
                         vis$Anatomical)

vis$Predicted_final_clusters <- sub(".*_", "", vis$dominant_c2l_cluster)

# ── Cluster label mapping ───────────────────────────────────────────────────
# This vector as originally written is new -> old:
#   names  = NEW cluster labels (e.g. '1a','1b',...,'18')
#   values = OLD Seurat/c2l cluster IDs (e.g. '0','7',...,'40')
new_to_old <- c('1a'='0',
                '1b'='7',
                '1c'='11',
                '1d'='13',
                '1e'='27',
                '1f'='29',
                '1g'='47',
                '2a'='2',
                '2b'='6',
                '2c'='17',
                '2d'='20',
                '2e'='23',
                '2f'='30',
                '2g'='31',
                '2h'='41',
                '2i'='43',
                '3'='9',
                '4'='24',
                '5'='26',
                '6a'='34',
                '6b'='38',
                '7'='37',
                '8'='39',
                '9a'= '1',
                '9b'='4',
                '9c'='5',
                '9d'='12',
                '9e'='16',
                '9f'='18',
                '10'='21',
                '11'='36',
                '12'='42',
                '13a'='8',
                '13b'='10',
                '13c'='19',
                '13d'='25',
                '13e'='45',
                '13f'='48',
                '14a'='3',
                '14b'='15',
                '14c'='32',
                '14d'='44',
                '15a'='14',
                '15b'='33',
                '16'='28',
                '17'='22',
                '18'='40'
)

# Invert to get old -> new: names = OLD cluster IDs, values = NEW labels
old_to_new <- setNames(names(new_to_old), new_to_old)

# Sanity check - should be no duplicated old IDs, i.e. inversion is 1:1
stopifnot(!anyDuplicated(new_to_old))

# vis$Predicted_renamed is now a CHARACTER vector (labels like "1a", "9c", etc.)
# NOT numeric, since new labels aren't all purely numeric
vis$Predicted_renamed <- unname(old_to_new[as.character(vis$Predicted)])

vis <- subset(vis, Predicted != '35')

heatmat <- table(vis$Predicted_renamed,
                 vis$Anatomical) %>%
  as.data.frame.matrix() %>%
  dplyr::select(-c('not sure 2',
                   'not sure 3',
                   'No Region',
                   'ventral diffuse',
                   'not sure 1',
                   'What',
                   'Diffuse',
                   'not sure 4'))

table(vis$Anatomical) %>% sort()
table(vis$Predicted_renamed) %>% sort()

# Define anchor colors
my_colors <- c("gray", 'white', "orange", "red")
color_func <- colorRampPalette(my_colors)
continuous_palette <- color_func(50)

# Scale by row proportions
piv <- (heatmat / rowSums(heatmat)) %>% scale()

# Sort rows by the order the labels appear in new_to_old (i.e. "1a","1b",...,"18")
# instead of as.numeric(), since labels like "1a" aren't numeric
piv <- piv[order(match(rownames(piv), names(new_to_old))), ]

p <- pheatmap(piv,
              cluster_rows = FALSE,
              cluster_cols = TRUE,
              treeheight_row = 0,
              treeheight_col = 0,
              color = continuous_palette,
              angle_col = 90,
              number_color = 'black')
p



# Apply the same old -> new relabeling used for vis, so labels match between plots
parker_2024$clusters49_renamed <- unname(old_to_new[as.character(parker_2024$clusters49)])

DimPlot(parker_2024, label = TRUE, group.by = 'clusters49_renamed')
DimPlot(parker_2024, label = TRUE, group.by = 'clusters49')

SpatialDimPlot(subset(vis, Predicted_renamed=='13d'),
               group.by ='Predicted_renamed',
               images = 's_6P17.polygons',
               image.scale = "hires")

SpatialFeaturePlot(vis,
                     'insm1b',
                     images = 's_6P17.polygons',
                   #  image.scale = "hires", 
                   pt.size.factor = 0.05)

svg(filename = 'C:/Users/Gabe/Desktop/multiome_poa/Manuscript/Plots/clusters49 c2l heatmap.svg',
    width = 8,
    height = 6.5)
p
dev.off()



