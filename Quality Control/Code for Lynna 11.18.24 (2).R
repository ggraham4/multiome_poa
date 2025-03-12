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


cells <- WhichCells(nemo, expression = pct.bad_nuc < 45 & 
                      pct.good_nuc > 55 & 
                      percent_top50_genes < 50 &
                      percent_top100_genes < 60 )
nemo <- subset(nemo, cells = cells)

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

gc()
gc()

# integration with Harmony
nemo <- IntegrateLayers(
  object = nemo, method = HarmonyIntegration,
  orig.reduction = "rnaPCA", new.reduction = "harmony.rna",
  verbose = TRUE, assay = 'RNA'
)


gc()
gc()

nemo <- JoinLayers(nemo)

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_integration_testing_10.18.24.rds')


####################


nemo <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_integration_testing_10.21.24.rds')

nemo[['RNA_orig']] <- nemo[['RNA']]
nemo[['RNA_orig']]@layers$scale.data <- NULL
all_genes <- rownames(nemo)
DefaultAssay(nemo) <- 'RNA_orig'

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
nemo.atac <- merge(x=s1.atac, y=c(s2.atac, s3.atac, s4.atac, s5.atac, s6.atac), merge.dr = T)

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

nemo.atac <- RunUMAP(nemo.atac, reduction = "lsi", dims = 2:30)

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

nemo <- RegionStats(nemo, genome = nemo_fasta)

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

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_gex.sct.atac_integration_10.27.24.rds')



#############################

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


gc()
gc()

DefaultAssay(nemo) <- 'RNA'
Idents(nemo) <- 'orig.ident'

gc()
gc()

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/clownfish/QC/nemo_gex.sct.atac_integration_10.27.24.rds')
