library(Seurat)
library(Signac)
library(paletteer)
library(ggtext)
library(tidyverse)
library(plotly)
library(patchwork)
library(ggplot2)
library(knitr)
library(tibble)
library(ggridges)
library(glmGamPoi)
library(rtracklayer)
library(biomartr)
library(igraph)
library(circlize)
library(webshot)
library(gplots)
library(ggrepel)
library(plotly)
library(kableExtra)
library(scCustomize)
library(DropletUtils)
library(DropletQC)
library(HDF5Array)
library(scater)
library(qlcMatrix)
library(hdf5r)
library(parallel)
library(parallelly)
library(future)
library(paletteer)
library(ggtext)
library(RColorBrewer)
library(viridis)
library(wesanderson)

options(future.globals.maxSize = 50000 * 1024^2)
mem.maxVSize(500000000)

options(parallelly.fork.enable = TRUE)
options(parallelly.supportsMulticore.unstable = 'quiet')
options(mc.cores = parallel::detectCores())

# options(future.globals.maxSize = 50000 * 1024^2)
# mem.maxVSize(500000000)

annotations <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/clownfish_genome_files/clownfish_genome_annotations.rds')
nemo_gene.table <- as.data.frame(annotations)

nemo_fasta <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/clownfish_genome_files/clownfish_fasta.rds')

nemo_id_info <- read_csv('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/clownfish_genome_files/clownfish_experiment_id_info.csv')

nemo_gene_info <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/clownfish_genome_files/clownfish_gene_info_homologs.rds')

ambient_markers <- read_csv('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/Katie/MZ_genome_files/M_zebra_UMD2a/ambient_markers_E.Caglayan2022.csv')

# nemo <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/clownfish_atac.gex_preQC_8.26.24.rds')


s1 <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/peaks_assay/s1_clownfish_gex.atac_common.peaks.rds')
s2 <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/peaks_assay/s2_clownfish_gex.atac_common.peaks.rds')
s3 <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/peaks_assay/s3_clownfish_gex.atac_common.peaks.rds')
s4 <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/peaks_assay/s4_clownfish_gex.atac_common.peaks.rds')
s5 <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/peaks_assay/s5_clownfish_gex.atac_common.peaks.rds')
s6 <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/peaks_assay/s6_clownfish_gex.atac_common.peaks.rds')

nemo <- merge(x=s1, y=c(s2,s3,s4,s5,s6), merge.data = TRUE)
DefaultAssay(nemo) <- 'RNA'
Idents(nemo) <- 'individual'

mito <- rownames(nemo)[grep('^KEG47-', rownames(nemo))]
nemo[['pct.mt']] <- PercentageFeatureSet(nemo, pattern = '^KEG47-')

ribogenes <- c(rownames(nemo)[grep('^rpl', rownames(nemo))], rownames(nemo)[grep('^rps', rownames(nemo))])
ribogenes <- rownames(nemo)[rownames(nemo) %in% ribogenes]
nemo[['pct.rb']] <- PercentageFeatureSet(nemo, features = ribogenes)

ribosomal_mito_genes <- c(rownames(nemo)[grep('^mrpl', rownames(nemo))], rownames(nemo)[grep('^mrps', rownames(nemo))])
ribosomal_mito_genes <- rownames(nemo)[rownames(nemo) %in% ribosomal_mito_genes]
nemo[['pct.mt.rb']] <- PercentageFeatureSet(nemo, features = ribosomal_mito_genes)

coding <- unique(nemo_gene.table$gene_id[nemo_gene.table$gene_biotype == 'protein_coding'])
coding <- coding[which(! coding %in% c(mito, ribosomal_mito_genes, ribogenes))]
cgenes <- rownames(nemo)[rownames(nemo) %in% coding]
nemo[['pct.cd']] <- PercentageFeatureSet(nemo, features = cgenes)

lncRNA <- nemo_gene.table$gene_id[nemo_gene.table$gene_biotype == 'lncRNA']
snRNA <- nemo_gene.table$gene_id[nemo_gene.table$gene_biotype == 'snRNA']
snoRNA <- nemo_gene.table$gene_id[nemo_gene.table$gene_biotype == 'snoRNA']
ncRNA <- nemo_gene.table$gene_id[nemo_gene.table$gene_biotype == 'ncRNA']
nuclear <- unique(c(lncRNA, snRNA, snoRNA, ncRNA))
nuclear <- rownames(nemo)[rownames(nemo) %in% nuclear]
nemo[['pct.nuclear']] <- PercentageFeatureSet(nemo, features = nuclear)

ambient_markers$Gene <- tolower(ambient_markers$Gene)
ambient_markers$nemo_homolog <- annotations$gene_id[match(ambient_markers$Gene, annotations$human_greatest_homolog, incomparables = T)]
ambient <- ambient_markers %>% dplyr::filter(PValue < 1e-3)
known_ambient <- rownames(nemo)[which(ambient$nemo_homolog %in% rownames(nemo))]
ambient <- unique(c(known_ambient))
nemo[['pct.ambient']] <- PercentageFeatureSet(nemo, features = ambient)

nemo <- Add_Cell_Complexity(object = nemo, assay = 'RNA', overwrite = T)
nemo$log10PeaksPerUMI <- log10(nemo$nFeature_ATAC)/log10(nemo$nCount_ATAC)
nemo$genes_per_umi <- nemo$nFeature_RNA/nemo$nCount_RNA
nemo$peaks_per_umi <- nemo$nFeature_ATAC/nemo$nCount_ATAC

nemo_gene_info <- nemo_gene.table %>% dplyr::filter(gbkey == 'Gene')
nemo_genes <- as.data.frame(nemo_gene_info$gene_id)
colnames(nemo_genes)[1] <- 'seurat_gene_name'
nemo_gene_info <- cbind(nemo_genes, nemo_gene_info)
nemo_gene_info$phase<-NULL
nemo_gene_info$transcript_id<-NULL
nemo_gene_info$NCBI<-NULL
nemo_gene_info[14:19] <- NULL
nemo_gene_info$gene <- NULL
nemo_gene_info$gene_name <- NULL

s.genes <- tolower(cc.genes$s.genes)
g2m.genes <- tolower(cc.genes$g2m.genes)

s_genes <- nemo_gene_info$seurat_gene_name[which(s.genes %in% nemo_gene_info$human_greatest_homolog)]
nomatch_s.genes <- s.genes[which(! s.genes %in% nemo_gene_info$human_greatest_homolog)]
g2m_genes <- nemo_gene_info$seurat_gene_name[which(g2m.genes %in% nemo_gene_info$human_greatest_homolog)]
nomatch_g2m.genes <- nemo_gene_info$seurat_gene_name[which(! g2m.genes %in% nemo_gene_info$human_greatest_homolog)]

nemo[['pct_s.genes']] <- PercentageFeatureSet(nemo, features = s_genes)
nemo[['pct_g2m.genes']] <- PercentageFeatureSet(nemo, features = g2m_genes)

nemo <- Add_Top_Gene_Pct_Seurat(seurat_object = nemo, num_top_genes = 50, assay = 'RNA', meta_col_name = 'percent_top50_genes', overwrite = T)
nemo <- Add_Top_Gene_Pct_Seurat(seurat_object = nemo, num_top_genes = 100, assay = 'RNA', meta_col_name = 'percent_top100_genes', overwrite = T)

######################

clown.s1_demulti <- read_csv('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/deconvolution/souporcell_demuxlet_selected/s1_barcodes_sampleIDs.csv')
clown.s1_demulti <- clown.s1_demulti %>%
  mutate(id_barcode = paste0(Sample, sep = '_', barcode))
same_cells <- colnames(nemo)[which(colnames(nemo) %in% clown.s1_demulti$id_barcode)]
clown.s1_demulti <- clown.s1_demulti %>% dplyr::filter(id_barcode %in% same_cells)

clown.s2_demulti <- read_csv('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/deconvolution/souporcell_demuxlet_selected/s2_barcodes_sampleIDs.csv')
clown.s2_demulti <- clown.s2_demulti %>%
  mutate(id_barcode = paste0(Sample, sep = '_', barcode))
same_cells <- colnames(nemo)[which(colnames(nemo) %in% clown.s2_demulti$id_barcode)]
clown.s2_demulti <- clown.s2_demulti %>% dplyr::filter(id_barcode %in% same_cells)

clown.s3_demulti <- read_csv('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/deconvolution/souporcell_demuxlet_selected/s3_barcodes_sampleIDs.csv')
clown.s3_demulti <- clown.s3_demulti %>%
  mutate(id_barcode = paste0(Sample, sep = '_', barcode))
same_cells <- colnames(nemo)[which(colnames(nemo) %in% clown.s3_demulti$id_barcode)]
clown.s3_demulti <- clown.s3_demulti %>% dplyr::filter(id_barcode %in% same_cells)

clown.s4_demulti <- read_csv('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/deconvolution/souporcell_demuxlet_selected/s4_barcodes_sampleIDs.csv')
clown.s4_demulti <- clown.s4_demulti %>%
  mutate(id_barcode = paste0(Sample, sep = '_', barcode))
same_cells <- colnames(nemo)[which(colnames(nemo) %in% clown.s4_demulti$id_barcode)]
clown.s4_demulti <- clown.s4_demulti %>% dplyr::filter(id_barcode %in% same_cells)

clown.s5_demulti <- read_csv('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/deconvolution/souporcell_demuxlet_selected/s5_barcodes_sampleIDs.csv')
clown.s5_demulti <- clown.s5_demulti %>%
  mutate(id_barcode = paste0(Sample, sep = '_', barcode))
same_cells <- colnames(nemo)[which(colnames(nemo) %in% clown.s5_demulti$id_barcode)]
clown.s5_demulti <- clown.s5_demulti %>% dplyr::filter(id_barcode %in% same_cells)

clown.s6_demulti <- read_csv('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/deconvolution/souporcell_demuxlet_selected/s6_barcodes_sampleIDs.csv')
clown.s6_demulti <- clown.s6_demulti %>%
  mutate(id_barcode = paste0(Sample, sep = '_', barcode))
same_cells <- colnames(nemo)[which(colnames(nemo) %in% clown.s6_demulti$id_barcode)]
clown.s6_demulti <- clown.s6_demulti %>% dplyr::filter(id_barcode %in% same_cells)

deconv_ids <- rbind(clown.s1_demulti, clown.s2_demulti, clown.s3_demulti, clown.s4_demulti, clown.s5_demulti, clown.s6_demulti)

##################


nemo$individual <- deconv_ids$Sample[match(colnames(nemo), deconv_ids$id_barcode)]
Idents(nemo) <- 'orig.ident'

nemo$Tank <- nemo_id_info$Tank[match(nemo$individual, nemo_id_info$Fish)]
nemo$Trigger <- nemo_id_info$Trigger[match(nemo$individual, nemo_id_info$Fish)]
nemo$Prev_tank <- nemo_id_info$Prev_tank[match(nemo$individual, nemo_id_info$Fish)]
nemo$Condition <- nemo_id_info$Condition[match(nemo$individual, nemo_id_info$Fish)]
nemo$Status <- nemo_id_info$Status[match(nemo$individual, nemo_id_info$Fish)]
nemo$Status_Long <- nemo_id_info$Status_Long[match(nemo$individual, nemo_id_info$Fish)]
nemo$Time_Day_2 <- nemo_id_info$Time_Day_2[match(nemo$individual, nemo_id_info$Fish)]
nemo$Behaviors_Day_2 <- nemo_id_info$Behaviors_Day_2[match(nemo$individual, nemo_id_info$Fish)]
nemo$Change_Length <- nemo_id_info$Change_Length[match(nemo$individual, nemo_id_info$Fish)]
nemo$Change_Mass <- nemo_id_info$Change_Mass[match(nemo$individual, nemo_id_info$Fish)]
nemo$Log_11KT <- nemo_id_info$Log_11KT[match(nemo$individual, nemo_id_info$Fish)]
nemo$Log_E2 <- nemo_id_info$Log_E2[match(nemo$individual, nemo_id_info$Fish)]
nemo$Average_Area_2.5x <- nemo_id_info$Average_Area_2.5x[match(nemo$individual, nemo_id_info$Fish)]
nemo$Total_Slides_With_Gonad <- nemo_id_info$Total_Slides_With_Gonad[match(nemo$individual, nemo_id_info$Fish)]
nemo$Estimated_Volume_2.5x <- nemo_id_info$Estimated_Volume_2.5x[match(nemo$individual, nemo_id_info$Fish)]
nemo$Log10_Volume <- nemo_id_info$Log10_Volume[match(nemo$individual, nemo_id_info$Fish)]
nemo$Percent_Testicular <- nemo_id_info$Percent_Testicular[match(nemo$individual, nemo_id_info$Fish)]
nemo$Testicular_Estimate <- nemo_id_info$Testicular_Estimate[match(nemo$individual, nemo_id_info$Fish)]
nemo$Log10_Testicular_Estimate <- nemo_id_info$Log10_Testicular_Estimate[match(nemo$individual, nemo_id_info$Fish)]
nemo$Log_E2 <- NULL

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_merged.pools_gex.atac_preQC_NEW.rds')

nemo <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/clownfish_merged.pools_gex.atac_preQC.rds')

gc()
gc()

# hex codes <- c("#334666", "#7484A6", "#9DABC8", "#E8B7C3", "#C99FA9", "#81BEAE", "#B3DDD0", "#D1DCE2", "#F5B993", "#ED9C6B", "#FF7B89", "#8A5082", "#6F5F90", "#A5CAD2", "#F3CCDB", "#E5E5F3")

perPool_colorPallet <- paletteer_d("ggthemes::excel_Parallax")
perIndividual_colorPallet <- c('#1170AA', '#1BA3C6FF', '#2CB5C0FF', '#30BCADFF', '#21B087FF', '#5FBB68', '#57A337FF', '#A2B627FF', '#F4DE3A', '#F8B620FF', '#F89217FF', '#E76618', 'red2', '#C8133B', '#F64971FF', '#FC719EFF', '#EB73B3FF', '#CE69BEFF', "#8A5082", '#A26DC2FF', '#7873C0FF', "#6F5F90")

######################## 


test <- nemo
meta <- nemo@meta.data
DefaultAssay(nemo) <- 'RNA'
Idents(nemo) <- 'orig.ident'

dim(nemo)
table(nemo$orig.ident)
table(nemo$individual)
table(nemo$Status_Long)

nemo_qc.metrics <- as_tibble(
  nemo[[]],
  rownames="Cell.Barcode"
)

nemo_qc.metrics %>%
  ggplot(aes(pct.mt)) + 
  geom_histogram(binwidth = 0.5, fill="darkorange", colour="black") +
  ggtitle("sample (afternoon) Distribution of Percentage Mitochondrion") +
  geom_vline(xintercept = 20, color = 'red3')
nemo_qc.metrics %>%
  arrange(pct.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=pct.mt)) + 
  geom_point(size=0.7) + 
  scale_color_gradientn(colors=wes_palette("Zissou1", 100, type = "continuous")) +
  ggtitle("Lognormalized QC metrics: testfish POA % MT genes highlighted") +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 30000) +
  geom_vline(xintercept = 300) +
  geom_vline(xintercept = 30000) +
  stat_smooth(method="lm", color = 'blue3') +
  scale_x_log10() + scale_y_log10()
nemo_qc.metrics %>%
  arrange(pct.cd) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=pct.cd)) + 
  geom_point(size=0.7) + 
  scale_color_gradientn(colors=wes_palette("Zissou1", 100, type = "continuous")) +
  ggtitle("Lognormalized QC metrics: testfish POA % coding genes highlighted") +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 30000) +
  geom_vline(xintercept = 300) +
  geom_vline(xintercept = 30000) +
  stat_smooth(method="lm", color = 'blue3') +
  scale_x_log10() + scale_y_log10()
nemo_qc.metrics %>%
  arrange(pct.rb) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=pct.rb)) + 
  geom_point(size=0.7) + 
  scale_color_gradientn(colors=wes_palette("Zissou1", 100, type = "continuous")) +
  ggtitle("Lognormalized QC metrics: testfish POA % coding genes highlighted") +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 30000) +
  geom_vline(xintercept = 300) +
  geom_vline(xintercept = 30000) +
  stat_smooth(method="lm", color = 'blue3') +
  scale_x_log10() + scale_y_log10()
nemo_qc.metrics %>%
  arrange(pct.nuclear) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=pct.nuclear)) + 
  geom_point(size=0.7) + 
  scale_color_gradientn(colors=wes_palette("Zissou1", 100, type = "continuous")) +
  ggtitle("Lognormalized QC metrics: testfish POA % coding genes highlighted") +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 30000) +
  geom_vline(xintercept = 300) +
  geom_vline(xintercept = 30000) +
  stat_smooth(method="lm", color = 'blue3') +
  scale_x_log10() + scale_y_log10()



a <- VlnPlot(nemo, 
             features = c('nCount_RNA'), 
             group.by = "orig.ident", 
             ncol = 1,
             log = T, 
             pt.size = 0.1,
             cols = paletteer_d("ggthemes::excel_Parallax"),
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 25000, color = 'red') +
  geom_hline(yintercept = 200, color = 'red')
b <- VlnPlot(nemo, 
             features = c('nFeature_RNA'), 
             group.by = "orig.ident", 
             ncol = 1,
             log = T, 
             pt.size = 0.1, 
             cols = paletteer_d("ggthemes::excel_Parallax"),
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 9000, color = 'red') +
  geom_hline(yintercept = 200, color = 'red')
c <- VlnPlot(nemo, 
             features = c('nCount_ATAC'), 
             group.by = "orig.ident", 
             ncol = 1,
             log = T, 
             cols = paletteer_d("ggthemes::excel_Parallax"),
             pt.size = 0.1, 
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 25000, color = 'red')
d <- VlnPlot(nemo, 
             features = c('nFeature_ATAC'), 
             group.by = "orig.ident", 
             ncol = 1,
             log = T, 
             cols = paletteer_d("ggthemes::excel_Parallax"),
             pt.size = 0.1, 
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 15000, color = 'red') 
p1 <- wrap_plots(a,b,c,d, ncol =2)
p1


a <- VlnPlot(nemo, 
             features = c('pct.mt'), 
             group.by = "orig.ident", 
             ncol = 1,
             log = F, 
             cols = paletteer_d("ggthemes::excel_Parallax"),
             pt.size = 0.1, 
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 20, color = 'red')
b <- VlnPlot(nemo, 
             features = c('pct.cd'), 
             group.by = "orig.ident", 
             ncol = 1,
             log = F,  
             cols = paletteer_d("ggthemes::excel_Parallax"),
             pt.size = 0.1, 
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 65, color = 'red')
c <- VlnPlot(nemo, 
             features = c('pct.rb'), 
             group.by = "orig.ident", 
             ncol = 1,
             log = F,  
             cols = paletteer_d("ggthemes::excel_Parallax"),
             pt.size = 0.1, 
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 4, color = 'red')
d <- VlnPlot(nemo, 
             features = c('pct.mt.rb'), 
             group.by = "orig.ident", 
             ncol = 1,
             log = F,  
             cols = paletteer_d("ggthemes::excel_Parallax"),
             pt.size = 0.1, 
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 0.5, color = 'red')
p2 <- wrap_plots(a,b,c,d, ncol =2)
p2

var_feats_df <- FindVariableFeatures(nemo, nfeatures = 28004, verbose = T, assay = 'RNA')
var_feats_df <- var_feats_df@assays[["RNA"]]@meta.data
rownames(var_feats_df) <- var_feats_df$var.features


var_feats_df$log_variance <- log10(var_feats_df$vf_vst_counts_variance)
var_feats_df$log_variance.expected <- log10(var_feats_df$vf_vst_counts_variance.expected)

var_feats_df <- Add_Pct_Diff(var_feats_df, pct.1_name = 'vf_vst_counts_variance.expected', pct.2_name = 'vf_vst_counts_variance', overwrite = T)
colnames(var_feats_df)[colnames(var_feats_df) == 'pct_diff'] = 'variance_pct.diff'
var_feats_df <- Add_Pct_Diff(var_feats_df, pct.1_name = 'log_variance.expected', pct.2_name = 'log_variance', overwrite = T)
colnames(var_feats_df)[colnames(var_feats_df) == 'pct_diff'] = 'log.variance_pct.diff'

good_nuc <- var_feats_df %>% dplyr::filter(vf_vst_counts_variance > 0 & vf_vst_counts_variance < 5)
good_nuc <- good_nuc$var.features
good_nuc <- good_nuc[which(! good_nuc %in% c(mito, ribogenes, ribosomal_mito_genes))]

bad_nuc <- var_feats_df %>% dplyr::filter(vf_vst_counts_variance > 5)
bad_nuc <- bad_nuc$var.features
bad_nuc <- bad_nuc[which(! bad_nuc %in% good_nuc)]

nemo[['pct.bad_nuc']] <- PercentageFeatureSet(nemo, features = bad_nuc)
good_nuc <- good_nuc[which(! good_nuc %in% bad_nuc)]
nemo[['pct.good_nuc']] <- PercentageFeatureSet(nemo, features = good_nuc)

a <- QC_Plot_UMIvsGene(seurat_object = nemo, meta_gradient_name = "pct.good_nuc", meta_gradient_low_cutoff = 55, low_cutoff_gene = 200, low_cutoff_UMI = 200) +
  scale_x_log10() + 
  scale_y_log10() + labs(title = 'Percentage of Variable "good" Genes per Nuc')
a

b <- QC_Plot_UMIvsGene(seurat_object = nemo, meta_gradient_name = "pct.bad_nuc", meta_gradient_low_cutoff = 45, low_cutoff_gene = 200, low_cutoff_UMI = 200) +
  scale_x_log10() + 
  scale_y_log10() + labs(title = 'Percentage of Variable Ambient Genes per Nuc')
b

p3 <- wrap_plots(a,b,ncol=1)
p3

c <- QC_Plot_UMIvsGene(seurat_object = nemo, meta_gradient_name = "percent_top50_genes", meta_gradient_low_cutoff = 50, low_cutoff_gene = 200, low_cutoff_UMI = 200) +
  scale_x_log10() + 
  scale_y_log10() + labs(title = 'Percentage of Top 50 Genes per Nuc')
d <- QC_Plot_UMIvsGene(seurat_object = nemo, meta_gradient_name = "percent_top100_genes", meta_gradient_low_cutoff = 60,low_cutoff_gene = 200, low_cutoff_UMI = 200) +
  scale_x_log10() + 
  scale_y_log10() + labs(title = 'Percentage of Top 100 Genes per Nuc')
p4 <- wrap_plots(c, d, ncol=1)
p4

p5 <- wrap_plots(p3,p4, ncol = 2)
p5 

cells <- WhichCells(nemo, expression = pct.bad_nuc < 45 & 
                      pct.good_nuc > 55 & 
                      percent_top50_genes < 50 &
                      percent_top100_genes < 60 )
nemo <- subset(nemo, cells = cells)


nemo@meta.data %>%
  arrange(pct.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=pct.mt)) + 
  geom_point(size=0.7) + 
  scale_color_gradientn(colors=wes_palette("Zissou1", 100, type = "continuous")) +
  ggtitle("Lognormalized QC metrics: testfish POA % MT genes highlighted") +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 30000) +
  geom_vline(xintercept = 300) +
  geom_vline(xintercept = 30000) +
  stat_smooth(method="lm", color = 'blue3') +
  scale_x_log10() + scale_y_log10()
nemo@meta.data %>%
  arrange(pct.cd) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=pct.cd)) + 
  geom_point(size=0.7) + 
  scale_color_gradientn(colors=wes_palette("Zissou1", 100, type = "continuous")) +
  ggtitle("Lognormalized QC metrics: testfish POA % coding genes highlighted") +
  geom_hline(yintercept = 300) +
  geom_hline(yintercept = 30000) +
  geom_vline(xintercept = 300) +
  geom_vline(xintercept = 30000) +
  stat_smooth(method="lm", color = 'blue3') +
  scale_x_log10() + scale_y_log10()


a <- VlnPlot(nemo, 
             features = c('pct.mt'), 
             group.by = "individual", 
             ncol = 1,
             log = F, 
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.1) + NoLegend() + 
  geom_hline(yintercept = 20, color = 'red')
b <- VlnPlot(nemo, 
             features = c('pct.cd'), 
             group.by = "individual", 
             ncol = 1,
             log = F,  
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.1) + NoLegend() + 
  geom_hline(yintercept = 65, color = 'red')
c <- VlnPlot(nemo, 
             features = c('pct.rb'), 
             group.by = "individual", 
             ncol = 1,
             log = F,  
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.1) + NoLegend() + 
  geom_hline(yintercept = 4, color = 'red')
d <- VlnPlot(nemo, 
             features = c('pct.mt.rb'), 
             group.by = "individual", 
             ncol = 1,
             log = F,  
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.1) + NoLegend() + 
  geom_hline(yintercept = 0.5, color = 'red')
p4 <- wrap_plots(a,b,c,d, ncol =2)
p4


nemo <- subset(nemo, subset = pct.mt < 15 &
                 pct.cd > 85 &
                 pct.rb < 2 &
                 pct.mt.rb < 0.4)
table(nemo$individual)


a <- VlnPlot(nemo, 
             features = c('nCount_RNA'), 
             group.by = "individual", 
             ncol = 1,
             log = T, 
             pt.size = 0.1,
             cols = perIndividual_colorPallet,
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 25000, color = 'red') +
  geom_hline(yintercept = 200, color = 'red')
b <- VlnPlot(nemo, 
             features = c('nFeature_RNA'), 
             group.by = "individual", 
             ncol = 1,
             log = T, 
             pt.size = 0.1, 
             cols = perIndividual_colorPallet,
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 8500, color = 'red') +
  geom_hline(yintercept = 200, color = 'red')
c <- VlnPlot(nemo, 
             features = c('nCount_ATAC'), 
             group.by = "individual", 
             ncol = 1,
             log = T, 
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 25000, color = 'red')
d <- VlnPlot(nemo, 
             features = c('nFeature_ATAC'), 
             group.by = "individual", 
             ncol = 1,
             log = T, 
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 16000, color = 'red') 
p5 <- wrap_plots(a,b,c,d, ncol =2)
p5

nemo <- subset(nemo, subset = nCount_RNA < 25000 &
                 nCount_RNA > 200 &
                 nFeature_RNA < 8500 &
                 nFeature_RNA > 200 &
                 nCount_ATAC < 25000 &
                 nFeature_ATAC < 16000)
table(nemo$individual)

a <- VlnPlot(nemo, 
             features = c('pct.mt'), 
             group.by = "individual", 
             ncol = 1,
             log = F, 
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.1) + NoLegend() + 
  geom_hline(yintercept = 20, color = 'red')
b <- VlnPlot(nemo, 
             features = c('pct.cd'), 
             group.by = "individual", 
             ncol = 1,
             log = F,  
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.1) + NoLegend() + 
  geom_hline(yintercept = 65, color = 'red')
c <- VlnPlot(nemo, 
             features = c('pct.rb'), 
             group.by = "individual", 
             ncol = 1,
             log = F,  
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.1) + NoLegend() + 
  geom_hline(yintercept = 4, color = 'red')
d <- VlnPlot(nemo, 
             features = c('pct.mt.rb'), 
             group.by = "individual", 
             ncol = 1,
             log = F,  
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.1) + NoLegend() + 
  geom_hline(yintercept = 0.5, color = 'red')
p6 <- wrap_plots(a,b,c,d, ncol =2)
p6


a <- VlnPlot(nemo, 
             features = c('nCount_RNA'), 
             group.by = "individual", 
             ncol = 1,
             log = T, 
             pt.size = 0.1,
             cols = perIndividual_colorPallet,
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 25000, color = 'red') +
  geom_hline(yintercept = 200, color = 'red')
b <- VlnPlot(nemo, 
             features = c('nFeature_RNA'), 
             group.by = "individual", 
             ncol = 1,
             log = T, 
             pt.size = 0.1, 
             cols = perIndividual_colorPallet,
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 8500, color = 'red') +
  geom_hline(yintercept = 200, color = 'red')
c <- VlnPlot(nemo, 
             features = c('nCount_ATAC'), 
             group.by = "individual", 
             ncol = 1,
             log = T, 
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 25000, color = 'red')
d <- VlnPlot(nemo, 
             features = c('nFeature_ATAC'), 
             group.by = "individual", 
             ncol = 1,
             log = T, 
             cols = perIndividual_colorPallet,
             pt.size = 0.1, 
             alpha = 0.05) + NoLegend() + 
  geom_hline(yintercept = 16000, color = 'red') 
p7 <- wrap_plots(a,b,c,d, ncol =2)
p7



g <- QC_Plot_UMIvsGene(seurat_object = nemo, meta_gradient_name = "pct.mt", meta_gradient_low_cutoff = 15, low_cutoff_gene = 250, low_cutoff_UMI = 250) +
  scale_x_log10() + 
  scale_y_log10()
h <- QC_Plot_UMIvsGene(seurat_object = nemo, meta_gradient_name = "pct.cd", meta_gradient_low_cutoff = 85,low_cutoff_gene = 250, low_cutoff_UMI = 250) +
  scale_x_log10() + 
  scale_y_log10()



VlnPlot(nemo, 
        features = c('pct.mt'), 
        group.by = "orig.ident", 
        ncol = 1,log = F, 
        pt.size = 0.1, alpha = 0.1) + NoLegend() + geom_hline(yintercept = 12, color = 'red')
nemo <- subset(nemo, subset = pct.mt < 12)
VlnPlot(nemo, 
        features = c('pct.mt', 'pct.cd'), 
        group.by = "orig.ident", 
        ncol = 2,log = F, 
        pt.size = 0.1, alpha = 0.1) + NoLegend()
table(nemo$individual)



density_nCountATAC_tss <- DensityScatter(nemo, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = c(1,5,95,99)) + 
  scale_color_viridis(option = "B", 
                      breaks = c(0.05, 0.15, 0.3, 0.5, 0.7, 1.0), 
                      labels = c(0.05, 0.15, 0.3, 0.5, 0.7, 1.0))
density_nCountATAC_tss <- density_nCountATAC_tss + labs(title = 'Quantiles: Density of Cells by TSS Enrichment')
density_nCountATAC_tss

a <- VlnPlot(nemo, 
             features = c('TSS.enrichment'), 
             group.by = "orig.ident", 
             ncol = 1,log = F, 
             pt.size = 0.1, alpha = 0.1) + NoLegend() + geom_hline(yintercept = 0.8, color = 'red')
b <- VlnPlot(nemo, 
             features = c('nucleosome_signal'), 
             group.by = "orig.ident", 
             ncol = 1,log = F, 
             pt.size = 0.1, alpha = 0.1) + NoLegend() + geom_hline(yintercept = 2, color = 'red')
wrap_plots(a,b)

nemo <- subset(nemo, subset = TSS.enrichment > 0.8 &
                 TSS.enrichment < 2 &
                 nucleosome_signal < 2)
a <- VlnPlot(nemo, 
             features = c('TSS.enrichment'), 
             group.by = "orig.ident", 
             ncol = 1,log = F, 
             pt.size = 0.1, alpha = 0.1) + NoLegend() + geom_hline(yintercept = 1, color = 'red')
b <- VlnPlot(nemo, 
             features = c('nucleosome_signal'), 
             group.by = "orig.ident", 
             ncol = 1,log = F, 
             pt.size = 0.1, alpha = 0.1) + NoLegend() + geom_hline(yintercept = 2, color = 'red')
wrap_plots(a,b)

table(nemo$individual)


tmp_nemo <- nemo

gc()
gc()

tmp_nemo[['ATAC']]<-NULL
tmp_nemo[['SCT']]<-NULL
tmp_nemo[['RNA']]@layers$scale.data <- NULL

gc()
gc()

tmp_nemo <- NormalizeData(tmp_nemo, assay = 'RNA')
tmp_nemo <- FindVariableFeatures(tmp_nemo, nfeatures = 4000)
tmp_nemo <- ScaleData(tmp_nemo)

gc()
gc()

tmp_nemo <- RunPCA(tmp_nemo, 
                   npcs = 50, 
                   features = rownames(tmp_nemo[['RNA']]))
tmp_nemo <- CellCycleScoring(tmp_nemo, s.features = s_genes, g2m.features = g2m_genes)

nemo[['S.Score']] <- tmp_nemo[['S.Score']]
nemo[['G2M.Score']] <- tmp_nemo[['G2M.Score']]

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/nemo_postQC_testing_10.17.24_v2.rds')


#############################


clown <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/nemo_postQC_testing_10.18.24_v2.rds')

# clown[["RNA"]] <- JoinLayers(clown[["RNA"]])

gc()
gc()

nemo = clown

table(nemo$individual)
clown_median_stats <- Median_Stats(seurat_object = nemo, group_by_var = "individual", median_var = c('pct.mt', 'pct.cd', 'nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC'))
clown_median_stats
nemo$median_UMI_count <- clown_median_stats$Median_nCount_RNA[match(nemo$individual, clown_median_stats$individual)]
nemo$median_gene_count <- clown_median_stats$Median_nFeature_RNA[match(nemo$individual, clown_median_stats$individual)]
nemo$median_peak_count <- clown_median_stats$Median_nCount_ATAC[match(nemo$individual, clown_median_stats$individual)]
nemo$median_fragment_count <- clown_median_stats$Median_nFeature_ATAC[match(nemo$individual, clown_median_stats$individual)]
nemo$median_pct.mt <- clown_median_stats$Median_pct.mt[match(nemo$individual, clown_median_stats$individual)]


##################



nemo_tmp <- nemo
nemo_tmp[['ATAC']]<-NULL
nemo_tmp[['SCT']]<-NULL
nemo_tmp <- NormalizeData(nemo_tmp)
nemo_tmp <- FindVariableFeatures(nemo_tmp, nfeatures = 4000)
nemo_tmp <- ScaleData(nemo_tmp, 
                  vars.to.regress = 'nuc_prep_batch',
                  scale.max = 50,
                  assay = 'RNA')
nemo <- CellCycleScoring(nemo_tmp, s.features = s_genes, g2m.features = g2m_genes)



##################


gc()
gc()

DefaultAssay(nemo) <- 'RNA'
Idents(nemo) <- 'orig.ident'
nemo[["RNA"]] <- split(nemo[["RNA"]], f = nemo$orig.ident)

nemo <- NormalizeData(nemo)
nemo <- FindVariableFeatures(nemo, nfeatures = 4000)

gc()
gc()

nemo <- ScaleData(nemo, 
                       vars.to.regress = 'nuc_prep_batch',
                       scale.max = 50,
                       assay = 'RNA')

nemo <- RunPCA(nemo, dim = 50, verbose = TRUE, assay = 'RNA', features = rownames(nemo[['RNA']]), reduction.name = 'rnaPCA', reduction.key = 'rnaPCA_')


ElbowPlot(nemo, ndims = ncol(Embeddings(nemo, "rnaPCA")), reduction = 'rnaPCA')
DepthCor(nemo, reduction = 'rnaPCA', assay = 'RNA', n=50)

p1 <- DimPlot(object = nemo, reduction = "rnaPCA", pt.size = .1, group.by = "orig.ident", shuffle = T) + 
  labs(title = 'PCA (Pre-Integration)', subtitle = 'Dimensionality of Merged Clownfish Pools \n') + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 15, vjust = 0), plot.subtitle = 
          element_text(hjust = 0.5, 
                       size = 12, 
                       face="italic", 
                       vjust = -3))
p1
p2 <- VlnPlot(object = nemo, features = "rnaPCA_1", group.by = "orig.ident", pt.size = .1)
p1 + p2

gc()
gc()

# integration with Harmony
nemo <- IntegrateLayers(
  object = nemo, method = HarmonyIntegration,
  orig.reduction = "rnaPCA", new.reduction = "harmony.rna",
  verbose = TRUE, assay = 'RNA'
)


harmony.rna <- DimPlot(object = nemo, reduction = "harmony.rna", pt.size = .1, group.by = "orig.ident", shuffle = T) + 
  labs(title = 'Harmony Integration', subtitle = 'Dimensionality of Merged Clownfish Pools \n') + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 15, vjust = 0), plot.subtitle = 
          element_text(hjust = 0.5, 
                       size = 12, 
                       face="italic", 
                       vjust = -3))
harmony.rna
wrap_plots(p1, harmony.rna)

gc()
gc()

# integration with reciprocal PCA
nemo <- IntegrateLayers(
  object = nemo, method = RPCAIntegration,
  orig.reduction = "rnaPCA", new.reduction = "rpca.rna",
  verbose = TRUE
)


rpca.rna <- DimPlot(object = nemo, reduction = "rpca.rna", pt.size = .1, group.by = "orig.ident", shuffle = T) + 
  labs(title = 'RPCA Integration', subtitle = 'Dimensionality of Merged Clownfish Pools \n') + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 15, vjust = 0), plot.subtitle = 
          element_text(hjust = 0.5, 
                       size = 12, 
                       face="italic", 
                       vjust = -3))
rpca.rna 
wrap_plots(p1, rpca.rna)

gc()
gc()

# integration with CCA
nemo <- IntegrateLayers(
  object = nemo, method = CCAIntegration,
  normalization.method = "LogNormalize",
  orig.reduction = "rnaPCA", new.reduction = "cca.rna",
  verbose = TRUE
)


cca.rna <- DimPlot(object = nemo, reduction = "cca.rna", pt.size = .1, group.by = "orig.ident", shuffle = T) + 
  labs(title = 'CCA Integration', subtitle = 'Dimensionality of Merged Clownfish Pools \n') + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 15, vjust = 0), plot.subtitle = 
          element_text(hjust = 0.5, 
                       size = 12, 
                       face="italic", 
                       vjust = -3))
cca.rna
wrap_plots(p1, cca.rna)

gc()
gc()

rna_integration <- wrap_plots(p1, harmony.rna, rpca.rna, cca.rna, ncol = 2)
rna_integration


nemo <- JoinLayers(nemo)

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_integration_testing_10.18.24.rds')


####################


nemo <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_integration_testing_10.21.24.rds')

nemo[['RNA_orig']] <- nemo[['RNA']]
nemo[['RNA_orig']]@layers$scale.data <- NULL
all_genes <- rownames(nemo)
DefaultAssay(nemo) <- 'RNA_orig'

nemo <- SCTransform(nemo,
                    vst.flavor = 'v2',
                    verbose = T,
                    variable.features.n = 4000,
                    return.only.var.genes = F) %>%
  RunPCA(dim = 50, verbose = TRUE, features = all_genes, reduction.name = 'sctPCA', reduction.key = 'sctPCA_')

gc()
gc()

nemo[['RNA_orig']] <- NULL

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_sct.integration_testing_10.27.24.rds')

gc()
gc()

DefaultAssay(nemo) <- 'SCT'

# integration with Harmony
nemo <- IntegrateLayers(
  object = nemo, 
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "sctPCA", 
  new.reduction = "harmony.sct",
  verbose = TRUE
)

p3 <- DimPlot(object = nemo, reduction = "harmony.sct", pt.size = .1, group.by = "orig.ident", shuffle = T)
p4 <- VlnPlot(object = nemo, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
p3 + p4

gc()
gc()

# integration with reciprocal PCA
nemo <- IntegrateLayers(
  object = nemo, 
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "sctPCA", 
  new.reduction = "rpca.sct",
  verbose = TRUE
)


p3 <- DimPlot(object = nemo, reduction = "rpca.sct", pt.size = .1, group.by = "orig.ident", shuffle = T)
p4 <- VlnPlot(object = nemo, features = "rpcasct_1", group.by = "orig.ident", pt.size = .1)
p3 + p4

gc()
gc()

# integration with CCA
nemo <- IntegrateLayers(
  object = nemo, method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "sctPCA", new.reduction = "cca.sct",
  verbose = TRUE
)

p3 <- DimPlot(object = nemo, reduction = "cca.sct", pt.size = .1, group.by = "orig.ident", shuffle = T)
p4 <- VlnPlot(object = nemo, features = "ccasct_1", group.by = "orig.ident", pt.size = .1)
p3 + p4

gc()
gc()

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_sct.integration_testing_10.27.24.rds')



####################



nemo.atac <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_atac_integration_testing_10.21.24.rds')
nemo <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_integration_testing_10.21.24.rds')

atac <- nemo
DefaultAssay(atac) <- 'ATAC'
atac[['RNA']] <- NULL
nemo.atac_list <- SplitObject(atac, split.by = 'orig.ident')

s1.atac <- nemo.atac_list[["pool_1"]]
s2.atac <- nemo.atac_list[["pool_2"]]
s3.atac <- nemo.atac_list[["pool_3"]]
s4.atac <- nemo.atac_list[["pool_4"]]
s5.atac <- nemo.atac_list[["pool_5"]]
s6.atac <- nemo.atac_list[["pool_6"]]

# compute LSI
DefaultAssay(s1.atac) <- "ATAC"
s1.atac <- FindTopFeatures(s1.atac, min.cutoff = 10)
s1.atac <- RunTFIDF(s1.atac)
s1.atac <- RunSVD(s1.atac, 
                  n= 50,
                  features = VariableFeatures(s1.atac))

gc()
gc()

# compute LSI
DefaultAssay(s2.atac) <- "ATAC"
s2.atac <- FindTopFeatures(s2.atac, min.cutoff = 10)
s2.atac <- RunTFIDF(s2.atac)
s2.atac <- RunSVD(s2.atac, 
                  n= 50,
                  features = VariableFeatures(s2.atac))

gc()
gc()

# compute LSI
DefaultAssay(s3.atac) <- "ATAC"
s3.atac <- FindTopFeatures(s3.atac, min.cutoff = 10)
s3.atac <- RunTFIDF(s3.atac)
s3.atac <- RunSVD(s3.atac, 
                  n= 50,
                  features = VariableFeatures(s3.atac))

gc()
gc()

# compute LSI
DefaultAssay(s4.atac) <- "ATAC"
s4.atac <- FindTopFeatures(s4.atac, min.cutoff = 10)
s4.atac <- RunTFIDF(s4.atac)
s4.atac <- RunSVD(s4.atac, 
                  n= 50,
                  features = VariableFeatures(s4.atac))

gc()
gc()

# compute LSI
DefaultAssay(s5.atac) <- "ATAC"
s5.atac <- FindTopFeatures(s5.atac, min.cutoff = 10)
s5.atac <- RunTFIDF(s5.atac)
s5.atac <- RunSVD(s5.atac, 
                  n= 50,
                  features = VariableFeatures(s5.atac))

gc()
gc()

# compute LSI
DefaultAssay(s6.atac) <- "ATAC"
s6.atac <- FindTopFeatures(s6.atac, min.cutoff = 10)
s6.atac <- RunTFIDF(s6.atac)
s6.atac <- RunSVD(s6.atac, 
                  n= 50,
                  features = VariableFeatures(s6.atac))
ElbowPlot(s6.atac, ndims = 50, reduction = "lsi")

gc()
gc()

# merge
nemo.atac <- merge(x=s1.atac, y=c(s2.atac, s3.atac, s4.atac, s5.atac, s6.atac), merge.dr = T))

nemo.atac_list <- SplitObject(nemo.atac, split.by = 'orig.ident')

saveRDS(nemo.atac_list, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_atac_integration_testing_10.21.24.rds')

gc()
gc()


saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_integration_testing_10.21.24.rds')


# process the combined dataset
nemo.atac <- FindTopFeatures(nemo.atac, min.cutoff = 'q5')
nemo.atac <- RunTFIDF(nemo.atac)
nemo.atac <- RunSVD(nemo.atac, 
                    n= 50,
                    features = VariableFeatures(nemo.atac))
ElbowPlot(nemo.atac, ndims = 50, reduction = "lsi")

nemo.atac <- RunUMAP(nemo.atac, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(nemo.atac, group.by = "orig.ident")
p1

# atac.feats <- rownames(nemo.atac)
atac_feats <- VariableFeatures(nemo.atac)

gc()
gc()


# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = nemo.atac_list,
  anchor.features = atac_feats,
  reduction = "rlsi",
  dims = 2:30,
  n.trees = 100,
  l2.norm = F,
  scale = F,
  max.features = 500
)


gc()
gc()

# integrate LSI embeddings
nemo <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = nemo.atac[["lsi"]],
  new.reduction.name = "atacLSI",
)


###########################


DefaultAssay(nemo) <- 'ATAC'

head(rownames(nemo@assays$ATAC@data))
head(rownames(nemo@assays$ATAC@counts))
head(rownames(nemo@assays$ATAC@meta.features))
head(nemo@assays[["ATAC"]]@var.features)

rownames(nemo@assays$ATAC@data)   = stringr::str_replace(rownames(nemo@assays$ATAC@data),   "_", "-")
rownames(nemo@assays$ATAC@counts) = stringr::str_replace(rownames(nemo@assays$ATAC@counts), "_", "-")
rownames(nemo@assays$ATAC@meta.features) = stringr::str_replace(rownames(nemo@assays$ATAC@meta.features), "_", "-")
nemo@assays[["ATAC"]]@var.features = stringr::str_replace(nemo@assays[["ATAC"]]@var.features, "_", "-")

head(rownames(nemo@assays$ATAC@data))
head(rownames(nemo@assays$ATAC@counts))
head(rownames(nemo@assays$ATAC@meta.features))
head(nemo@assays[["ATAC"]]@var.features)

nemo <- RegionStats(nemo, genome = nemo_fasta)

rownames(nemo@assays$ATAC@data)   = stringr::str_replace(rownames(nemo@assays$ATAC@data),   "-", "_")
rownames(nemo@assays$ATAC@counts) = stringr::str_replace(rownames(nemo@assays$ATAC@counts), "-", "_")
rownames(nemo@assays$ATAC@meta.features) = stringr::str_replace(rownames(nemo@assays$ATAC@meta.features), "-", "_")
nemo@assays[["ATAC"]]@var.features = stringr::str_replace(nemo@assays[["ATAC"]]@var.features, "-", "_")

head(rownames(nemo@assays$ATAC@data))
head(rownames(nemo@assays$ATAC@counts))
head(rownames(nemo@assays$ATAC@meta.features))
head(nemo@assays[["ATAC"]]@var.features)

gc()
gc()

# create a new UMAP using the integrated embeddings
nemo <- RunUMAP(nemo, 
                  reduction = 'atacLSI', 
                  dims = 2:30, 
                  reduction.key = "atacUMAP_", 
                  reduction.name = 'atacUMAP', 
                  assay = 'ATAC')

nemo <- FindNeighbors(nemo, 
                      reduction = 'atacLSI', 
                      dims = 2:30,
                      assay = 'ATAC',                      
                      k.param = 50, 
                      prune.SNN = 0,
                      graph.name = c('atac_nn','atac_snn'))

nemo <- FindClusters(nemo, 
                     algorithm = 3,
                     resolution = 0.1, 
                     cluster.name = "lsi_res0.1_atacUMAP",
                     graph.name = 'atac_snn')

atacUMAP <- DimPlot(nemo, 
                        reduction = "atacUMAP", 
                        pt.size = 0.1, 
                        label = T, 
                        order = T, 
                        group.by = c('lsi_res0.1_atacUMAP'))
atacUMAP

gc()
gc()

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_gex.sct.atac_integration_10.27.24.rds')

# saveRDS(nemo.atac, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_atac_integration_rlsi_10.22.24.rds')

# nemo.atac<- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_atac_integration_rlsi_10.22.24.rds')

nemo[['ATAC']] <- NULL
nemo[['ATAC']] <- nemo.atac[['ATAC']]
nemo@reductions[['atacLSI']] <- nemo.atac@reductions[['lsi']]
nemo@reductions[['atacUMAP']] <- nemo.atac@reductions[['umap']]

nemo[['RNA']] <- JoinLayers(nemo[['RNA']])

#############################


gc()
gc()

DefaultAssay(nemo) <- 'RNA'

nemo <- RunUMAP(nemo,
                     dims = 1:50, min.dist = 0.4,
                     n.neighbors = 50,
                     metric = "euclidean",
                     reduction = 'harmony.rna',
                     reduction.name = "harmony_rnaUMAP",
                     assay = 'RNA')

nemo <- FindNeighbors(nemo, 
                      reduction = "harmony_rnaUMAP", 
                      k.param = 50 , 
                      dims = 1:2, 
                      prune.SNN = 0, 
                      assay = 'RNA',
                      graph.name = c('harmony_rna_nn','harmony_rna_snn'))

nemo <- FindClusters(nemo, 
                     resolution = 0.1, 
                     algorithm = 2, 
                     assay = 'RNA', 
                     cluster.name = "harmony_res0.1_rnaUMAP",
                     graph.name = 'harmony_rna_snn')

harmony_rnaUMAP <- DimPlot(nemo, 
                           reduction = "harmony_rnaUMAP", 
                           pt.size = 0.1, 
                           label = T, 
                           shuffle = T, 
                           order = T, 
                           group.by = c('harmony_res0.1_rnaUMAP'))
harmony_rnaUMAP 

gc()
gc()

nemo <- RunUMAP(nemo,
                     dims = 1:50, min.dist = 0.4,
                     n.neighbors = 50,
                     metric = "euclidean",
                     reduction = 'rpca.rna',
                     reduction.name = "rpca_rnaUMAP",
                     assay = 'RNA')

nemo <- FindNeighbors(nemo, 
                      reduction = "rpca_rnaUMAP", 
                      k.param = 50 , 
                      dims = 1:2, 
                      prune.SNN = 0, 
                      assay = 'RNA',
                      graph.name = c('rpca_rna_nn','rpca_rna_snn'))
nemo <- FindClusters(nemo, 
                     resolution = 0.1, 
                     algorithm = 2, 
                     assay = 'RNA', 
                     cluster.name = "rpca_res0.1_rnaUMAP",
                     graph.name = 'rpca_rna_snn')

rpca_rnaUMAP <- DimPlot(nemo, 
                        reduction = "rpca_rnaUMAP", 
                        pt.size = 0.1, 
                        label = T, 
                        shuffle = T, 
                        order = T, 
                        group.by = c('rpca_res0.1_rnaUMAP'))
rpca_rnaUMAP

gc()
gc()

nemo <- RunUMAP(nemo,
                     dims = 1:50, min.dist = 0.4,
                     n.neighbors = 50,
                     metric = "euclidean",
                     reduction = 'cca.rna',
                     reduction.name = "cca_rnaUMAP",
                     assay = 'RNA')

nemo <- FindNeighbors(nemo, 
                      reduction = "cca_rnaUMAP", 
                      k.param = 50 , 
                      dims = 1:2, 
                      prune.SNN = 0, 
                      assay = 'RNA',
                      graph.name = c('cca_rna_nn','cca_rna_snn'))

nemo <- FindClusters(nemo, 
                     resolution = 0.1, 
                     algorithm = 2, 
                     assay = 'RNA', 
                     cluster.name = "cca_res0.1_rnaUMAP",
                     graph.name = 'cca_rna_snn')

cca_rnaUMAP <- DimPlot(nemo, 
                       reduction = "cca_rnaUMAP", 
                       pt.size = 0.1, 
                       label = T, 
                       shuffle = T, 
                       order = T, 
                       group.by = c('cca_res0.1_rnaUMAP'))
cca_rnaUMAP

rna_UMAP <- wrap_plots(harmony_rnaUMAP, rpca_rnaUMAP, cca_rnaUMAP, ncol = 3)
rna_UMAP

gc()
gc()

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_gex.sct.atac_integration_10.27.24.rds')



#############################



DefaultAssay(nemo) <- 'SCT'

nemo <- RunUMAP(nemo,
                dims = 1:50, min.dist = 0.4,
                n.neighbors = 50,
                metric = "euclidean",
                reduction = 'harmony.sct',
                reduction.name = "harmony_sctUMAP",
                assay = 'SCT')

nemo <- FindNeighbors(nemo, 
                      reduction = "harmony_sctUMAP", 
                      k.param = 50 , 
                      dims = 1:2, 
                      prune.SNN = 0, 
                      assay = 'SCT',
                      graph.name = c('harmony_sct_nn','harmony_sct_snn'))

nemo <- FindClusters(nemo, 
                     resolution = 0.1, 
                     algorithm = 2, 
                     assay = 'SCT', 
                     cluster.name = "harmony_res0.1_sctUMAP",
                     graph.name = 'harmony_sct_snn')

harmony_sctUMAP <- DimPlot(nemo, 
                           reduction = "harmony_sctUMAP", 
                           pt.size = 0.1, 
                           label = T, 
                           shuffle = T, 
                           order = T, 
                           group.by = c('harmony_res0.1_sctUMAP'))
harmony_sctUMAP 

gc()
gc()

nemo <- RunUMAP(nemo,
                dims = 1:50, min.dist = 0.4,
                n.neighbors = 50,
                metric = "euclidean",
                reduction = 'rpca.sct',
                reduction.name = "rpca_sctUMAP",
                assay = 'SCT')

nemo <- FindNeighbors(nemo, 
                      reduction = "rpca_sctUMAP", 
                      k.param = 50 , 
                      dims = 1:2, 
                      prune.SNN = 0, 
                      assay = 'SCT',
                      graph.name = c('rpca_sct_nn','rpca_sct_snn'))
nemo <- FindClusters(nemo, 
                     resolution = 0.1, 
                     algorithm = 2, 
                     assay = 'SCT', 
                     cluster.name = "rpca_res0.1_sctUMAP",
                     graph.name = 'rpca_sct_snn')

rpca_sctUMAP <- DimPlot(nemo, 
                        reduction = "rpca_sctUMAP", 
                        pt.size = 0.1, 
                        label = T, 
                        shuffle = T, 
                        order = T, 
                        group.by = c('rpca_res0.1_sctUMAP'))
rpca_sctUMAP

gc()
gc()

nemo <- RunUMAP(nemo,
                dims = 1:50, min.dist = 0.4,
                n.neighbors = 50,
                metric = "euclidean",
                reduction = 'cca.sct',
                reduction.name = "cca_sctUMAP",
                assay = 'SCT')

nemo <- FindNeighbors(nemo, 
                      reduction = "cca_sctUMAP", 
                      k.param = 50 , 
                      dims = 1:2, 
                      prune.SNN = 0, 
                      assay = 'SCT',
                      graph.name = c('cca_sct_nn','cca_sct_snn'))

nemo <- FindClusters(nemo, 
                     resolution = 0.1, 
                     algorithm = 2, 
                     assay = 'SCT', 
                     cluster.name = "cca_res0.1_sctUMAP",
                     graph.name = 'cca_sct_snn')

cca_sctUMAP <- DimPlot(nemo, 
                       reduction = "cca_sctUMAP", 
                       pt.size = 0.1, 
                       label = T, 
                       shuffle = T, 
                       order = T, 
                       group.by = c('cca_res0.1_sctUMAP'))
cca_sctUMAP

sct_UMAP <- wrap_plots(harmony_sctUMAP, rpca_sctUMAP, cca_sctUMAP, ncol = 3)
sct_UMAP

gc()
gc()

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_gex.sct.atac_integration_10.27.24.rds')



############################


gc()
gc()

DefaultAssay(nemo) <- 'RNA'
Idents(nemo) <- 'orig.ident'

#Joint UMAP visualization
# build a joint neighbor graph using both assays - WNN
nemo <- FindMultiModalNeighbors(
  object = nemo,
  reduction.list = list("harmony.rna", 'atacLSI'),
  dims.list = list(1:50, 2:50),
  verbose = TRUE,
  prune.SNN = 0,
  knn.graph.name = "harmony.wknn",
  snn.graph.name = "harmony.wsnn",
  weighted.nn.name = "harmony.wnn"
)

# build a joint UMAP visualization
nemo <- RunUMAP(
  object = nemo,
  nn.name = "harmony.wnn",
  reduction.name = "harmony_wnn.umap",
  reduction.key = "harmony_wnnUMAP_",
  verbose = TRUE,
  min.dist = 0.4,
  metric = "euclidean"
)

nemo <- FindClusters(nemo,
                     graph.name = "harmony.wsnn",
                     algorithm = 3,
                     verbose = T,
                     resolution = 0.1,
                     cluster.name = 'harmony.wnn_res0.1_clusters')



p1 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'harmony_wnn.umap', group.by = 'harmony.wnn_res0.1_clusters') +
  ggplot2::ggtitle("Clownfish Joint WNN UMAP: RNA + ATAC")
p2 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'harmony_rnaUMAP', group.by = c('harmony.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (RNA assay ONLY)")
p3 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'atacUMAP', group.by = c('harmony.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (ATAC assay ONLY)")
harmony_rna.atac_wnnUMAP <- wrap_plots(p1,p2,p3)
harmony_rna.atac_wnnUMAP

DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'harmony_wnn.umap', split.by = 'Status_Long') +
  ggplot2::ggtitle("Clownfish Joint RNA + ATAC WNN UMAP: Harmony")

gc()
gc()

gc()
gc()

DefaultAssay(nemo) <- 'RNA'
Idents(nemo) <- 'orig.ident'

#Joint UMAP visualization
# build a joint neighbor graph using both assays - WNN
nemo <- FindMultiModalNeighbors(
  object = nemo,
  reduction.list = list("rpca.rna", 'atacLSI'),
  dims.list = list(1:50, 2:50),
  verbose = TRUE,
  prune.SNN = 0,
  knn.graph.name = "rpca.wknn",
  snn.graph.name = "rpca.wsnn",
  weighted.nn.name = "rpca.wnn"
)

# build a joint UMAP visualization
nemo <- RunUMAP(
  object = nemo,
  nn.name = "rpca.wnn",
  reduction.name = "rpca_wnn.umap",
  reduction.key = "rpca_wnnUMAP_",
  verbose = TRUE,
  min.dist = 0.4,
  metric = "euclidean"
)

nemo <- FindClusters(nemo,
                     graph.name = "rpca.wsnn",
                     algorithm = 3,
                     verbose = T,
                     resolution = 0.1,
                     cluster.name = 'rpca.wnn_res0.1_clusters')



p1 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'rpca_wnn.umap', group.by = 'rpca.wnn_res0.1_clusters') +
  ggplot2::ggtitle("Clownfish Joint WNN UMAP: RNA + ATAC")
p2 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'rpca_rnaUMAP', group.by = c('rpca.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (RNA assay ONLY)")
p3 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'atacUMAP', group.by = c('rpca.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (ATAC assay ONLY)")
rpca_rna.atac_wnnUMAP <- wrap_plots(p1,p2,p3)
rpca_rna.atac_wnnUMAP

DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'rpca_wnn.umap', split.by = 'Status_Long') +
  ggplot2::ggtitle("Clownfish Joint RNA + ATAC WNN UMAP: rpca")

gc()
gc()

gc()
gc()

DefaultAssay(nemo) <- 'RNA'
Idents(nemo) <- 'orig.ident'

#Joint UMAP visualization
# build a joint neighbor graph using both assays - WNN
nemo <- FindMultiModalNeighbors(
  object = nemo,
  reduction.list = list("cca.rna", 'atacLSI'),
  dims.list = list(1:50, 2:50),
  verbose = TRUE,
  prune.SNN = 0,
  knn.graph.name = "cca.wknn",
  snn.graph.name = "cca.wsnn",
  weighted.nn.name = "cca.wnn"
)

# build a joint UMAP visualization
nemo <- RunUMAP(
  object = nemo,
  nn.name = "cca.wnn",
  reduction.name = "cca_wnn.umap",
  reduction.key = "cca_wnnUMAP_",
  verbose = TRUE,
  min.dist = 0.4,
  metric = "euclidean"
)

nemo <- FindClusters(nemo,
                     graph.name = "cca.wsnn",
                     algorithm = 3,
                     verbose = T,
                     resolution = 0.1,
                     cluster.name = 'cca.wnn_res0.1_clusters')



p1 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'cca_wnn.umap', group.by = 'cca.wnn_res0.1_clusters') +
  ggplot2::ggtitle("Clownfish Joint WNN UMAP: RNA + ATAC")
p2 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'cca_rnaUMAP', group.by = c('cca.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (RNA assay ONLY)")
p3 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'atacUMAP', group.by = c('cca.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (ATAC assay ONLY)")
cca_rna.atac_wnnUMAP <- wrap_plots(p1,p2,p3)
cca_rna.atac_wnnUMAP

DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'cca_wnn.umap', split.by = 'Status_Long') +
  ggplot2::ggtitle("Clownfish Joint RNA + ATAC WNN UMAP: cca")

gc()
gc()

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_gex.sct.atac_integration_10.27.24.rds')


############################


gc()
gc()

DefaultAssay(nemo) <- 'SCT'
Idents(nemo) <- 'orig.ident'

#Joint UMAP visualization
# build a joint neighbor graph using both assays - WNN
nemo <- FindMultiModalNeighbors(
  object = nemo,
  reduction.list = list("harmony.sct", 'atacLSI'),
  dims.list = list(1:50, 2:50),
  verbose = TRUE,
  prune.SNN = 0,
  knn.graph.name = "harmony.sct.wknn",
  snn.graph.name = "harmony.sct.wsnn",
  weighted.nn.name = "harmony.sct.wnn"
)

# build a joint UMAP visualization
nemo <- RunUMAP(
  object = nemo,
  nn.name = "harmony.sct.wnn",
  reduction.name = "harmony.sct_wnn.umap",
  reduction.key = "harmony.sct_wnnUMAP_",
  verbose = TRUE,
  min.dist = 0.4,
  metric = "euclidean"
)

nemo <- FindClusters(nemo,
                     graph.name = "harmony.sct.wsnn",
                     algorithm = 3,
                     verbose = T,
                     resolution = 0.1,
                     cluster.name = 'harmony.sct.wnn_res0.1_clusters')



p1 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'harmony.sct_wnn.umap', group.by = 'harmony.sct.wnn_res0.1_clusters') +
  ggplot2::ggtitle("Clownfish Joint WNN UMAP: SCT + ATAC")
p2 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'harmony_sctUMAP', group.by = c('harmony.sct.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (SCT assay ONLY)")
p3 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'atacUMAP', group.by = c('harmony.sct.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (ATAC assay ONLY)")
harmony_sct.atac_wnnUMAP <- wrap_plots(p1,p2,p3)
harmony_sct.atac_wnnUMAP

DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'harmony.sct_wnn.umap', split.by = 'Status_Long') +
  ggplot2::ggtitle("Clownfish Joint SCT + ATAC WNN UMAP: Harmony")

gc()
gc()

#Joint UMAP visualization
# build a joint neighbor graph using both assays - WNN
nemo <- FindMultiModalNeighbors(
  object = nemo,
  reduction.list = list("rpca.sct", 'atacLSI'),
  dims.list = list(1:50, 2:50),
  verbose = TRUE,
  prune.SNN = 0,
  knn.graph.name = "rpca.sct.wknn",
  snn.graph.name = "rpca.sct.wsnn",
  weighted.nn.name = "rpca.sct.wnn"
)

# build a joint UMAP visualization
nemo <- RunUMAP(
  object = nemo,
  nn.name = "rpca.sct.wnn",
  reduction.name = "rpca.sct_wnn.umap",
  reduction.key = "rpca.sct_wnnUMAP_",
  verbose = TRUE,
  min.dist = 0.4,
  metric = "euclidean"
)

nemo <- FindClusters(nemo,
                     graph.name = "rpca.sct.wsnn",
                     algorithm = 3,
                     verbose = T,
                     resolution = 0.1,
                     cluster.name = 'rpca.sct.wnn_res0.1_clusters')



p1 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'rpca.sct_wnn.umap', group.by = 'rpca.sct.wnn_res0.1_clusters') +
  ggplot2::ggtitle("Clownfish Joint WNN UMAP: SCT + ATAC")
p2 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'sct_rpcaUMAP', group.by = c('rpca.sct.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (SCT assay ONLY)")
p3 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'atacUMAP', group.by = c('rpca.sct.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (ATAC assay ONLY)")
rpca_sct.atac_wnnUMAP <- wrap_plots(p1,p2,p3)
rpca_sct.atac_wnnUMAP

DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'rpca.sct_wnn.umap', split.by = 'Status_Long') +
  ggplot2::ggtitle("Clownfish Joint SCT + ATAC WNN UMAP: rpca")

gc()
gc()


#Joint UMAP visualization
# build a joint neighbor graph using both assays - WNN
nemo <- FindMultiModalNeighbors(
  object = nemo,
  reduction.list = list("cca.sct", 'atacLSI'),
  dims.list = list(1:50, 2:50),
  verbose = TRUE,
  prune.SNN = 0,
  knn.graph.name = "cca.sct.wknn",
  snn.graph.name = "cca.sct.wsnn",
  weighted.nn.name = "cca.sct.wnn"
)

# build a joint UMAP visualization
nemo <- RunUMAP(
  object = nemo,
  nn.name = "cca.sct.wnn",
  reduction.name = "cca.sct_wnn.umap",
  reduction.key = "cca.sct_wnnUMAP_",
  verbose = TRUE,
  min.dist = 0.4,
  metric = "euclidean"
)

nemo <- FindClusters(nemo,
                     graph.name = "cca.sct.wsnn",
                     algorithm = 3,
                     verbose = T,
                     resolution = 0.1,
                     cluster.name = 'cca.sct.wnn_res0.1_clusters')



p1 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'cca.sct_wnn.umap', group.by = 'cca.sct.wnn_res0.1_clusters') +
  ggplot2::ggtitle("Clownfish Joint WNN UMAP: SCT + ATAC")
p2 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'cca_sctUMAP', group.by = c('cca.sct.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (SCT assay ONLY)")
p3 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'atacUMAP', group.by = c('cca.sct.wnn_res0.1_clusters')) + ggplot2::ggtitle("Clownfish WNN Clusters (ATAC assay ONLY)")
cca_sct.atac_wnnUMAP <- wrap_plots(p1,p2,p3)
cca_sct.atac_wnnUMAP

DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'cca.sct_wnn.umap', split.by = 'Status_Long') +
  ggplot2::ggtitle("Clownfish Joint SCT + ATAC WNN UMAP: cca")

gc()
gc()

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_gex.sct.atac_integration_10.27.24.rds')



wnnUMAP1 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'harmony_wnn.umap', group.by = 'harmony.wnn_res0.1_clusters') +
  ggplot2::ggtitle("Harmony RNA + ATAC WNN UMAP")

wnnUMAP2 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'rpca_wnn.umap', group.by = 'rpca.wnn_res0.1_clusters') +
  ggplot2::ggtitle("RPCA RNA + ATAC WNN UMAP")
  
wnnUMAP3 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'cca_wnn.umap', group.by = c('cca.wnn_res0.1_clusters')) + ggplot2::ggtitle("CCA RNA + ATAC WNN UMAP")

wnnUMAP4 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'harmony.sct_wnn.umap', group.by = 'harmony.sct.wnn_res0.1_clusters') +
  ggplot2::ggtitle("Harmony SCT + ATAC WNN UMAP")

wnnUMAP5 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'rpca.sct_wnn.umap', group.by = 'rpca.sct.wnn_res0.1_clusters') +
  ggplot2::ggtitle("RPCA SCT + ATAC WNN UMAP")

wnnUMAP6 <- DimPlot(nemo, pt.size = 0.1, label = T, shuffle = T, reduction = 'cca.sct_wnn.umap', group.by = c('cca.sct.wnn_res0.1_clusters')) + ggplot2::ggtitle("CCA SCT + ATAC WNN UMAP")
wnnUMAPs <- wrap_plots(wnnUMAP1, wnnUMAP2, wnnUMAP3, wnnUMAP4, wnnUMAP5, wnnUMAP6)
wnnUMAPs


