### Weighted Nearest Neighbor (WNN) Analysis
library(ggpubr) 

nemo <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/Seurat_Data/seurat_objects/nemo_gex.sct.atac_integration_10.27.24.rds')

DefaultAssay(nemo) <- 'RNA'
Idents(nemo) <- 'individual'

n_pc <- c(30, 35, 40, 45, 50)
n_lsi <- c(30, 35, 40, 45, 50)
knn <- c(35, 45)
resolutions <- c(0.4, 0.6, 0.8, 1, 1.4, 1.6)

for (pc in n_pc) {
  for (lsi in n_lsi) {
    for (k in knn) {
      
      nemo <- FindMultiModalNeighbors(
        object = nemo,
        reduction.list = list("harmony.rna", 'atacLSI'),
        dims.list = list(1:pc, 1:lsi),        
        verbose = TRUE,
        prune.SNN = 0,
        knn.graph.name = "wknn",
        snn.graph.name = "wsnn",
        weighted.nn.name = "wnn",
        modality.weight.name = c("RNA.weight", "ATAC.weight"),
        k.nn = k,
      )
      
      nemo <- RunUMAP(nemo, 
                      nn.name = "harmony.wnn", 
                      verbose = TRUE, 
                      n.neighbors = k, 
                      n.epochs = 500,
                      min.dist = 0.4, 
                      metric = "euclidean", 
                      reduction.name = "wnn.umap", 
                      reduction.key = "wnnUMAP_") 
      
      
      DefaultAssay(nemo) <- "RNA"
      a1 <- FeaturePlot_scCustom(nemo, "s1pr1", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("s1pr1 (RG)"))))) & NoAxes() & NoLegend() 
      
      a2 <- FeaturePlot_scCustom(nemo, "myrf", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("Myrf (Oligo)"))))) & NoAxes() & NoLegend() 
      
      a3 <- FeaturePlot_scCustom(nemo, "gpr17", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) +
        ggtitle(expression(atop(paste(italic("gpr17 (NFO)"))))) & NoAxes() & NoLegend()
      
      a4 <- FeaturePlot_scCustom(nemo, "csf1rb", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("csf1r (MG)"))))) & NoAxes() & NoLegend() 
      
      a5 <- FeaturePlot_scCustom(nemo, "aplnra", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) +
        ggtitle(expression(atop(paste(italic("aplnra (OPC)"))))) & NoAxes() & NoLegend()
      
      a6 <- FeaturePlot_scCustom(nemo, "tbx21", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) +
        ggtitle(expression(atop(paste(italic("tbx21 (Macrophages)"))))) & NoAxes() & NoLegend()
      
      a7 <- FeaturePlot_scCustom(nemo, "osr1", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("osr1 (Pericytes)"))))) & NoAxes() & NoLegend() 
      
      a8 <- FeaturePlot_scCustom(nemo, "pou4f1", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("pou4f1 (Excitatory Projection Neurons"))))) & NoAxes() & NoLegend() 
      
      a9 <- FeaturePlot_scCustom(nemo, "fsip1", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("fsip1 (Ependymal)"))))) & NoAxes() & NoLegend() 
      
      a10 <- FeaturePlot_scCustom(nemo, "slc17a6b", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("excitatory neurons (slc17a6b+)"))))) & NoAxes() & NoLegend() 
      
      a11 <- FeaturePlot_scCustom(nemo, "nrxn3a", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("inhibitory neurons (nrxn3a)"))))) & NoAxes() & NoLegend() 
      
      a12 <- FeaturePlot_scCustom(nemo, "pitx2", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("pitx2 (Ventrolateral POA Exciatory Neuroendocrine Neurons)"))))) & NoAxes() & NoLegend() 
      
      a13 <- FeaturePlot_scCustom(nemo, "sst1.1", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("sst1.1 (Somatostatinergic neurons)"))))) & NoAxes() & NoLegend() 
      
      a14 <- FeaturePlot_scCustom(nemo, "six3a", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("six3a (Neural Proginator Cells)"))))) & NoAxes() & NoLegend() 
      
      a15 <- FeaturePlot_scCustom(nemo, "shox2", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("shox2 (Ventromedial POA Excitatory Neuroendocrine Neurons)"))))) & NoAxes() & NoLegend() 
      
      a16 <- FeaturePlot_scCustom(nemo, "dcn", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("dcn (Fibroblasts)"))))) & NoAxes() & NoLegend() 
      
      a17 <- FeaturePlot_scCustom(nemo, "lhx6a", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("lhx6a (GABAergic interneuron progenitors)"))))) & NoAxes() & NoLegend() 
      
      a18 <- FeaturePlot_scCustom(nemo, "calb2a", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("calb2a (Mature excitatory high-activity neurons)"))))) & NoAxes() & NoLegend() 
      
      a19 <- FeaturePlot_scCustom(nemo, "npy", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("npy (NPY Neurons)"))))) & NoAxes() & NoLegend() 
      
      a20 <- FeaturePlot_scCustom(nemo, "hmx3a", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("hmx3a (Neurosecretory neuronal progenitors)"))))) & NoAxes() & NoLegend() 
      
      a21 <- FeaturePlot_scCustom(nemo, "tac1", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("tac1 (Mature excitatory neuroendocrine-modulatory neurons
)"))))) & NoAxes() & NoLegend() 
      
      a22 <- FeaturePlot_scCustom(nemo, "crhr1", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("crhr1 (CRH neurons)"))))) & NoAxes() & NoLegend() 
      
      a23 <- FeaturePlot_scCustom(nemo, "gata2b", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("gata2b (GABAergic Committed Neuronal Proginators)"))))) & NoAxes() & NoLegend() 
      
      a24 <- FeaturePlot_scCustom(nemo, "fat2", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("fat2 (Inhibitory Projection Neurons)"))))) & NoAxes() & NoLegend() 
      
      a25 <- FeaturePlot_scCustom(nemo, "gal", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 2) + 
        ggtitle(expression(atop(paste(italic("gal (Galaninergic neurons)"))))) & NoAxes() & NoLegend() 
      
      a26 <- FeaturePlot_scCustom(nemo, "trh", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("trh (TRH neurons)"))))) & NoAxes() & NoLegend() 
      
      a27 <- FeaturePlot_scCustom(nemo, "oxtrb", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("oxtrb (Oxytocinergic neurons)"))))) & NoAxes() & NoLegend() 
      
      a28 <- FeaturePlot_scCustom(nemo, "avpr2aa", order = T, colors_use = colors, reduction = "harmony_wnn.umap", min.cutoff = 'q10', max.cutoff = 'q98', na_cutoff = 1.5) + 
        ggtitle(expression(atop(paste(italic("avpr2aa (Vasopressinergic neurons)"))))) & NoAxes() & NoLegend() 
      
      
      
      a <- ggarrange(a1, a2, a3, a4, ncol = 2, nrow = 2, legend = "none")
      feature_plot <- annotate_figure(a, bottom = text_grob(paste0("n.neighbors = ", k, ", n.pc = ", pc, ", n.lsi = ", lsi),
                                                            color = "#227CDE", hjust = 1.2, x = 1, face = "italic", size = 12))
      
      print(feature_plot)
      
      a <- ggarrange(a5, a6, a7, a8, ncol = 2, nrow = 2, legend = "none")
      feature_plot <- annotate_figure(a, bottom = text_grob(paste0("n.neighbors = ", k, ", n.pc = ", pc, ", n.lsi = ", lsi),
                                                            color = "#227CDE", hjust = 1.2, x = 1, face = "italic", size = 12))
      print(feature_plot)
      
      a <- ggarrange(a9, a10, a11, a12, ncol = 2, nrow = 2, legend = "none")
      feature_plot <- annotate_figure(a, bottom = text_grob(paste0("n.neighbors = ", k, ", n.pc = ", pc, ", n.lsi = ", lsi),
                                                            color = "#227CDE", hjust = 1.2, x = 1, face = "italic", size = 12))
      print(feature_plot)
      
      a <- ggarrange(a13, a14, a15, a16, ncol = 2, nrow = 2, legend = "none")
      feature_plot <- annotate_figure(a, bottom = text_grob(paste0("n.neighbors = ", k, ", n.pc = ", pc, ", n.lsi = ", lsi),
                                                            color = "#227CDE", hjust = 1.2, x = 1, face = "italic", size = 12))
      print(feature_plot)
      
      a <- ggarrange(a17, a18, a19, a20, ncol = 2, nrow = 2, legend = "none")
      feature_plot <- annotate_figure(a, bottom = text_grob(paste0("n.neighbors = ", k, ", n.pc = ", pc, ", n.lsi = ", lsi),
                                                            color = "#227CDE", hjust = 1.2, x = 1, face = "italic", size = 12))
      print(feature_plot)
      
      a <- ggarrange(a21, a22, a23, a24, ncol = 2, nrow = 2, legend = "none")
      feature_plot <- annotate_figure(a, bottom = text_grob(paste0("n.neighbors = ", k, ", n.pc = ", pc, ", n.lsi = ", lsi),
                                                            color = "#227CDE", hjust = 1.2, x = 1, face = "italic", size = 12))
      print(feature_plot)
      
      a <- ggarrange(a25, a26, a27, a28, ncol = 2, nrow = 2, legend = "none")
      feature_plot <- annotate_figure(a, bottom = text_grob(paste0("n.neighbors = ", k, ", n.pc = ", pc, ", n.lsi = ", lsi),
                                                            color = "#227CDE", hjust = 1.2, x = 1, face = "italic", size = 12))
      print(feature_plot)
      
      for (res in resolutions) {
        
        cluster_name <- paste0("res", res, "_", k, "nn_", pc, "PC_", lsi, "LSI")
        
        nemo <- FindClusters(nemo, resolution = res, graph.name = "wsnn", algorithm = 3, cluster.name = cluster_name)
        
        b <- DimPlot(nemo, label = TRUE, reduction = "harmony_wnn.umap", group.by = cluster_name)
        dim_plot <- annotate_figure(b, bottom = text_grob(paste0("n.neighbors = ", k, ", n.pc = ", pc, ", n.lsi = ", lsi, ", res = ", res),
                                                          color = "#227CDE", hjust = 1.2, x = 1, face = "italic", size = 12))
        print(dim_plot)
        
        cat(paste("Completed clustering with ", pc, "PCs ", lsi, "LSIs, and ", k, "nearest neighbors at resolution = ", res, "\n"))
      }
    }
  }
}
dev.off()

gc()
gc()

combos <- expand.grid(Res = resolutions, KNN = knn, PC = n_pc, LSI = n_lsi)
combos <- combos[order(combos$Res, combos$KNN, combos$PC, combos$LSI), ]
combos$ResStr <- gsub("\\.$", "", as.character(combos$Res))
res_names <- paste0("res", combos$ResStr, "_", combos$KNN, "nn_", combos$PC, "PC_", combos$LSI, "LSI")

# computing Adjusted Rand Index (ARI) with Dune
clusMat <- nemo@meta.data %>% select(res_names)

ARI <- apply(clusMat, 2, function(x) {
  apply(clusMat, 2, function(y) {
    x_unc <- x
    y_unc <- y
    aricode::ARI(x_unc, y_unc)
  })
})
rownames(ARI) <- colnames(ARI)
row_means <- rowMeans(ARI)
max_mean <- names(row_means)[which.max(row_means)]

pdf("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/harmony_WNN_clustering_RandIndex_nemo.orig_all_cells.pdf", width = 20, height = 20)
df <- ARI %>% as.data.frame() %>%
  dplyr::mutate(label1 = rownames(ARI)) %>%
  tidyr::gather(key = label2, value = metric, -(ncol(ARI) + 1))

row_index <- which(levels(factor(df$label2)) == max_mean)

ht <- ggplot(df, aes(x = label1, y = label2, fill = metric)) +
  geom_tile() +
  scale_fill_gradientn(colours = brewer.pal(9, "Spectral"), 
                       values = scales::rescale(c(0, 0.25, 0.35, 0.45, 0.5, 0.75, 1)), 
                       limits = c(0, 1)) +
  theme_classic() +
  theme(axis.line = element_blank()) +
  geom_text(aes(label = round(metric, 2)), size = 4) +
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "Clusters from Tested Parameters", y = "Clusters from Tested Parameters", 
       title = "Adjusted Rand Index for joint snRNA-seq and snATAC-seq clustering nemoerparameters (snMultiome)") + 
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = row_index - 0.5,
                ymax = row_index + 0.5), 
            inherit.aes = FALSE, color = "#FF4777", fill = NA, linewidth = 1) 
annotate_figure(ht, bottom = text_grob(paste0("Optimal clustering parameters = ", max_mean, " with mean Rand Index = ", max(row_means)),
                                       color = "#227CDE", hjust = 1.2, x = 1, face = "italic", size = 12))
dev.off()
# the most stable clustering parameters are 30 nearest neighbors, 25 principal components, 30 latent semantic indexing, and resolution 1.3

ht

max_mean

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/Seurat_Data/seurat_objects/nemo.orig_harmony.integration_optimal_clusters.rds')
