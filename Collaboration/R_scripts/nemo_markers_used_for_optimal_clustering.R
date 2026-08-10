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