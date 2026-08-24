library(Seurat)
library(tidyverse)
library(Matrix)
library(pheatmap)

# ── Paths ──────────────────────────────────────────────────────────────────
results_dir   <- "C:/Users/Gabe/Desktop/Visium/POA/cell2location_map"
groupings_dir <- "C:/Users/Gabe/Desktop/multiome_poa/Visium/Groupings"
plot_dir      <- "C:/Users/Gabe/Desktop/Visium/plots"
dir.create(plot_dir, showWarnings = FALSE)

vis <- readRDS("C:/Users/Gabe/Desktop/Visium/vis_better_barcodes_c2l05_projection_parker_2024_clusters20.rds")
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

# ── Old cluster ID -> new cluster label mapping ─────────────────────────────
# names = OLD Seurat/c2l cluster IDs, values = NEW labels
# old cluster 13 is intentionally excluded (no mapping) - it's not used at all
# old cluster 3 maps to new label 13 (previously missing/left as NULL)
old_to_new <- c('1'='1',
                '2'='2',
                '6'='3',
                '9'='4',
                '10'='5',
                '12'='6',
                '15'='7',
                '16'='8',
                '0'='9',
                '7'='10',
                '14'='11',
                '18'='12',
                '3'='13',
                '4'='14',
                '5'='15',
                '11'='16',
                '8'='17',
                '17'='18',
                '19'='19')

vis$Predicted_renamed <- as.numeric(old_to_new[as.character(vis$Predicted)])

heatmat <- table(vis$Predicted_renamed[!is.na(vis$Predicted_renamed)],
                 vis$Anatomical[!is.na(vis$Predicted_renamed)]) %>%
  as.data.frame.matrix() %>%
  dplyr::select(-c('not sure 2',
                   'not sure 3',
                   'No Region',
                   'ventral diffuse',
                   'not sure 4',
                   'not sure 1',
                   'What',
                   'Diffuse'))

table(vis$Anatomical) %>% sort()

# Define anchor colors
my_colors <- c("gray", 'white', "orange", "red")
color_func <- colorRampPalette(my_colors)
continuous_palette <- color_func(50)

# Scale by row proportions
piv <- (heatmat / rowSums(heatmat)) %>% scale()

# Sort rows by the new numeric cluster label
piv <- piv[order(as.numeric(rownames(piv)), decreasing = FALSE), ]

p <- pheatmap(piv,
              cluster_rows = FALSE,
              cluster_cols = TRUE,
              treeheight_row = 0,
              treeheight_col = 0,
              color = continuous_palette,
              angle_col = 90,
              number_color = 'black')
p

parker_2024 <- readRDS("A:/Anemonefish POA Legacy R Objects/parker_2024_realigned_ocellaris.rds")
DimPlot(parker_2024, label = TRUE, group.by = 'clusters20')


svg(filename = 'C:/Users/Gabe/Desktop/multiome_poa/Manuscript/Plots/clusters20 c2l heatmap.svg',
    width =8,
    height = 4.5)
p
dev.off()
