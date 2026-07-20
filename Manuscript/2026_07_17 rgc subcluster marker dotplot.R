library(Seurat)
library(ggplot2)
library(tidyverse)
library(emmeans)
library(ggsignif)
library(lme4)


obj  = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

obj = FindSubCluster(obj,
                     1, 'harmony.wsnn', resolution = 0.2, subcluster.name = 'sub_res0.2')

sub_1 = subset(obj, res0.8_50nn_40PC_45LSI==1)

# Set the identity to the subclusters
Idents(sub_1) <- 'sub_res0.2'

# Make sure you're pulling markers from the RNA assay (not ATAC)
DefaultAssay(sub_1) <- 'RNA'

# Find positive markers for every subcluster
markers <- FindAllMarkers(
  sub_1,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# --- Keep only non-redundant markers (genes unique to a single subcluster) ---
non_redundant_markers <- markers %>%
  group_by(gene) %>%
  filter(n() == 1) %>%   # gene appears as a marker for only one cluster
  ungroup()

# --- Top 2 per subcluster, ranked by fold-change ---
top2_markers <- non_redundant_markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_max(order_by = avg_log2FC, n = 2, with_ties = FALSE) %>%
  ungroup()

top2_markers

# label 
labels  = c(
  '1_1' = '1_NSC',
  '1_0' = '1_NP',
  '1_2' = '1_shha',
  '1_3' = '1_gfap'
)
sub_1@meta.data$rgc_label = unname(labels[sub_1@meta.data$sub_res0.2])

Idents(sub_1)='rgc_label'

genes = c(
  'gja1b',
  'LOC111577263',
  'shha',
  'gfap',
  'angpt1'
)

p=DotPlot(sub_1, c(genes,
          top2_markers$gene))+
  coord_flip()+
  scale_x_discrete(labels = c(
    genes,
    top2_markers$gene
    ))+
  theme(
    axis.text.x = element_text(size = 6, angle =45),       # X-axis (feature) text size
    axis.text.y = element_text(size = 6),       # Y-axis (cell identity) text size
    axis.title = element_blank(), # Axis titles
    legend.text = element_text(size = 6),       # Legend labels
    legend.title = element_text(size=6),       # Legend title
    legend.position = 'top',
    legend.justification = "center"  
    
  )

p     

ggsave(plot = p,
      file = 'marker dotplot rgc.svg',
     device = "svg",
   units = "in",
  width = 2.19,
 height = 7.5,
path = 'Manuscript/Plots/Manuscript v1.3')
