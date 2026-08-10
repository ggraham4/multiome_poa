library(rhdf5)
library(data.table)

nemo <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/Seurat_Data/seurat_objects/nemo.orig_harmony.integration_optimal_clusters_v2.rds')

annotations <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/clownfish_genome_files/main_files/clownfish_genome_annotations.rds')


nemo$pool <- nemo$orig.ident
nemo$pool[nemo$pool == 'pool_1'] = 's1'
nemo$pool[nemo$pool == 'pool_2'] = 's2'
nemo$pool[nemo$pool == 'pool_3'] = 's3'
nemo$pool[nemo$pool == 'pool_4'] = 's4'
nemo$pool[nemo$pool == 'pool_5'] = 's5'
nemo$pool[nemo$pool == 'pool_6'] = 's6'
table(nemo$pool)

nemo$barcode_id <- rownames(nemo@meta.data)
nemo$barcode_orig <- sub(".*?_", "", nemo$barcode_id)
head(nemo$barcode_orig)

samples <- paste0("s", 1:6)


#################################


# for loop code to run for each sample provided in the 'samples' character list

#### NOTE: the code within the for loop can be replaced with any code as long as the seurat object is named ' seu ' and the plots you want to generate figures for in the pdf output are written as ' print(plot) '

for (sample in samples) {
  message("Processing ", sample, " …")
  
  # 1. Read data + metrics + EDM + linkages
  dat     <- Read10X_h5(paste0('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/cellranger_outs/', sample, '_clownfish/raw_feature_bc_matrix.h5'))

  seu_linkages <- fread(paste0("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/cellranger_outs/", sample, "_clownfish/analysis/feature_linkage/feature_linkage.bedpe"), header = FALSE)
  
  metrics <- read_csv(paste0("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/cellranger_outs/", sample, "_clownfish/per_barcode_metrics.csv"))
  
  # 2. Build Seurat + ATAC assay
  seu <- CreateSeuratObject(dat$`Gene Expression`)
  
  frags <- CreateFragmentObject(
    path = paste0('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/cellranger_outs/', sample, '_clownfish/atac_fragments.tsv.gz'),
    validate.fragments = TRUE)
  
  seu[['ATAC']] <- CreateChromatinAssay(
    counts = dat$Peaks, 
    fragments = frags, 
    annotation = annotations,
    sep = c(":", "-"))
  
  seu$pool <- sample
  seu.cells <- as.data.frame(nemo$barcode_orig[which(nemo$pool %in% seu$pool)])
  colnames(seu.cells) <- 'barcode'
  cells <- seu.cells$barcode
  seu <- subset(seu, cells = cells)
  dim(seu)
  seu.cells$barcode_id <- rownames(seu.cells)
  head(seu.cells)
  seu$barcode <- rownames(seu@meta.data)
  seu$barcode_id <- seu.cells$barcode_id[match(seu.cells$barcode, seu$barcode)]
  head(seu$barcode_id)
  head(seu$barcode)
  
  barcodes_keep <- seu.cells$barcode
  head(barcodes_keep)
  
  metrics <- metrics %>% dplyr::filter(barcode %in% barcodes_keep)

  # Assign column names based on 10x format
  # Assign column names based on 10x format
  colnames(seu_linkages)[1:13] <- c(
    "chrom1", "start1", "end1",
    "chrom2", "start2", "end2",
    "name", "score", "strand1", "strand2",
    "significance", "distance", "seu_linkages_type"
  )
  seu_linkages$fragment1 <- paste0(seu_linkages$chrom1, ":", seu_linkages$start1, "-", seu_linkages$end1)
  seu_linkages$fragment2 <- paste0(seu_linkages$chrom2, ":", seu_linkages$start2, "-", seu_linkages$end2)
  
  # Extract value1 (gene)
  seu_linkages$feature1 <- sub("^<(.*?)><.*$", "\\1", seu_linkages$name)
  
  seu_linkages$feature2 <- sub("^<.*?><(.*?)>$", "\\1", seu_linkages$name)
  
  split_column <- function(df, col_name) {
    # Find rows that contain a semicolon
    has_semicolon <- grepl(";", df[[col_name]])
    rows_to_change <- which(has_semicolon)
    
    # Duplicate those rows and keep second part
    duplicates <- df[rows_to_change, ]
    duplicates[[col_name]] <- sub(".*?;", "", duplicates[[col_name]])
    
    # In the original rows, keep the first part
    df[[col_name]][rows_to_change] <- sub(";.*", "", df[[col_name]][rows_to_change])
    
    # Combine original and new rows
    rbind(df, duplicates)
  }
  
  # Apply the function to both columns
  seu_linkages <- split_column(seu_linkages, "feature1")
  seu_linkages <- split_column(seu_linkages, "feature2")
  
  
  # Create match indices for the suffix
  matcheseu <- regexpr("_(distal|promoter|intergenic)", seu_linkages$feature1)
  
  matches2 <- regexpr("_(distal|promoter|intergenic)", seu_linkages$feature2)
  
  
  # Extract th# Extract th# Extract the matched strings (will be "" for no match, not dropped)
  seu_linkages[, feature1_seu_linkages_location := ifelse(matcheseu > 0, regmatches(feature1, matcheseu), NA)]
  
  seu_linkages[, feature2_seu_linkages_location := ifelse(matches2 > 0, regmatches(feature2, matches2), NA)]
  
  # Remove the matched suffix from value1
  seu_linkages$feature1 <- gsub("_(distal|promoter|intergenic)", "", seu_linkages$feature1)
  seu_linkages$feature1_seu_linkages_location <- gsub("_", "", seu_linkages$feature1_seu_linkages_location)
  
  seu_linkages$feature2 <- gsub("_(distal|promoter|intergenic)", "", seu_linkages$feature2)
  seu_linkages$feature2_seu_linkages_location <- gsub("_", "", seu_linkages$feature2_seu_linkages_location)
  
  seu_linkages$feature1[seu_linkages$feature1 == ""] <- NA
  seu_linkages$feature2[seu_linkages$feature2 == ""] <- NA
  
  table(seu_linkages$seu_linkages_type)
  
  peak_gene_seu_linkages <- seu_linkages %>% dplyr::filter(seu_linkages_type == 'peak-gene')
  peak_gene_seu_linkages$linked_feature1 <- peak_gene_seu_linkages$fragment1
  peak_gene_seu_linkages$fragment1 <- NULL
  peak_gene_seu_linkages$fragment2 <- NULL
  peak_gene_seu_linkages$feature1_associated_gene_symbol <- peak_gene_seu_linkages$feature1
  peak_gene_seu_linkages$feature1 <- NULL
  peak_gene_seu_linkages$feature1_seu_linkages_type <- peak_gene_seu_linkages$feature1_seu_linkages_location
  peak_gene_seu_linkages$feature1_seu_linkages_location <- NULL
  peak_gene_seu_linkages$linked_feature2 <- peak_gene_seu_linkages$feature2
  peak_gene_seu_linkages$feature2 <- NULL
  peak_gene_seu_linkages$feature2_associated_gene_symbol <- peak_gene_seu_linkages$linked_feature2
  peak_gene_seu_linkages$feature2_seu_linkages_type <- peak_gene_seu_linkages$feature2_seu_linkages_location
  peak_gene_seu_linkages$feature2_seu_linkages_location <- NULL
  
  peak_gene_seu_linkages$linked_feature1 <- gsub(":", "-", peak_gene_seu_linkages$linked_feature1)
  
  p_g_linked_featseu <- unique(peak_gene_seu_linkages$linked_feature1)
  p_g_linked_feats2 <- unique(peak_gene_seu_linkages$linked_feature2)
  
  gene_peak_seu_linkages <- seu_linkages %>% dplyr::filter(seu_linkages_type == 'gene-peak')
  gene_peak_seu_linkages$linked_feature1 <- gene_peak_seu_linkages$feature1
  gene_peak_seu_linkages$feature1 <- NULL
  gene_peak_seu_linkages$fragment1 <- NULL
  gene_peak_seu_linkages$feature1_associated_gene_symbol <- gene_peak_seu_linkages$linked_feature1
  gene_peak_seu_linkages$feature1_seu_linkages_type <- gene_peak_seu_linkages$feature1_seu_linkages_location
  gene_peak_seu_linkages$feature1_seu_linkages_location <- NULL
  gene_peak_seu_linkages$linked_feature2 <- gene_peak_seu_linkages$fragment2
  gene_peak_seu_linkages$fragment2 <- NULL
  gene_peak_seu_linkages$feature2_associated_gene_symbol <- gene_peak_seu_linkages$feature2
  gene_peak_seu_linkages$feature2 <- NULL
  gene_peak_seu_linkages$feature2_seu_linkages_type <- gene_peak_seu_linkages$feature2_seu_linkages_location
  gene_peak_seu_linkages$feature2_seu_linkages_location <- NULL
  
  gene_peak_seu_linkages$linked_feature2 <- gsub(":", "-", gene_peak_seu_linkages$linked_feature2)
  
  g_p_linked_featseu <- unique(gene_peak_seu_linkages$linked_feature1)
  g_p_linked_feats2 <- unique(gene_peak_seu_linkages$linked_feature2)
  
  peak_peak_seu_linkages <- seu_linkages %>% dplyr::filter(seu_linkages_type == 'peak-peak')
  peak_peak_seu_linkages$linked_feature1 <- peak_peak_seu_linkages$fragment1
  peak_peak_seu_linkages$fragment1 <- NULL
  peak_peak_seu_linkages$feature1_associated_gene_symbol <- peak_peak_seu_linkages$feature1
  peak_peak_seu_linkages$feature1 <- NULL
  peak_peak_seu_linkages$feature1_seu_linkages_type <- peak_peak_seu_linkages$feature1_seu_linkages_location
  peak_peak_seu_linkages$feature1_seu_linkages_location <- NULL
  peak_peak_seu_linkages$linked_feature2 <- peak_peak_seu_linkages$fragment2
  peak_peak_seu_linkages$fragment2 <- NULL
  peak_peak_seu_linkages$feature2_associated_gene_symbol <- peak_peak_seu_linkages$feature2
  peak_peak_seu_linkages$feature2 <- NULL
  peak_peak_seu_linkages$feature2_seu_linkages_type <- peak_peak_seu_linkages$feature2_seu_linkages_location
  peak_peak_seu_linkages$feature2_seu_linkages_location <- NULL
  
  peak_peak_seu_linkages$linked_feature1 <- gsub(":", "-", peak_peak_seu_linkages$linked_feature1)
  
  peak_peak_seu_linkages$linked_feature2 <- gsub(":", "-", peak_peak_seu_linkages$linked_feature2)
  
  p_p_linked_featseu <- unique(peak_peak_seu_linkages$linked_feature1)
  p_p_linked_feats2 <- unique(peak_peak_seu_linkages$linked_feature2)
  
  linked_peaks <- unique(c(p_g_linked_featseu, p_p_linked_featseu, p_p_linked_feats2, g_p_linked_feats2))
  linked_genes <- unique(c(p_g_linked_feats2,g_p_linked_featseu))
  
  seu_linkages <- rbind(peak_gene_seu_linkages, gene_peak_seu_linkages, peak_peak_seu_linkages)
  
  seu_linkages$seu_linkages_score <- seu_linkages$score
  seu_linkages$seu_linkages_significance <- seu_linkages$significance
  seu_linkages$seu_linkages_distance <- seu_linkages$distance
  colnames(seu_linkages)
  seu_linkages <- seu_linkages[, -c(1:12)]
  colnames(seu_linkages)
  
  
  # identify peak-gene, gene-peak, and peak-peak seu_linkages
  peak_linked_genes <- seu_linkages %>% dplyr::filter(
    seu_linkages_type == 'gene-peak' &
      seu_linkages_score > 0)
  genes_linked_to_peaks <- unique(peak_linked_genes$linked_feature1)
  genes_linked_to_peaks <- genes_linked_to_peaks[which(genes_linked_to_peaks %in% rownames(seu[['RNA']]))]
  
  gene_linked_peaks <- seu_linkages %>% dplyr::filter(
    seu_linkages_type == 'peak-gene' &
      seu_linkages_score > 0)
  peaks_linked_to_genes <- unique(c(gene_linked_peaks$linked_feature1))
  peaks_linked_to_genes <- gsub("_", "-", peaks_linked_to_genes)
  
  peak_linked_peaks <- seu_linkages %>% dplyr::filter(
    seu_linkages_type == 'peak-peak' &
      seu_linkages_score > 0)
  peaks_linked_to_peaks <- unique(c(peak_linked_peaks$linked_feature1))
  peaks_linked_to_peaks <- gsub("_", "-", peaks_linked_to_peaks)
  
  
  # add peak linkages as a QC metric into the metadata of seurat object
  seu[['pct.peaks_linked_to_genes']] <- PercentageFeatureSet(seu, features = peaks_linked_to_genes, assay = 'ATAC')
  seu[['pct.genes_linked_to_peaks']] <- PercentageFeatureSet(seu, features = genes_linked_to_peaks, assay = 'RNA')
  seu[['pct.peaks_linked_to_peaks']] <- PercentageFeatureSet(seu, features = peaks_linked_to_peaks, assay = 'ATAC')
  
  seu_metadata <- seu@meta.data
  
  metrics$pct.peaks_linked_to_genes <- seu_metadata$pct.peaks_linked_to_genes[match(metrics$barcode, seu_metadata$barcode, nomatch = NA)]
  
  metrics$pct.genes_linked_to_peaks <- seu_metadata$pct.genes_linked_to_peaks[match(metrics$barcode, seu_metadata$barcode, nomatch = NA)]
  
  metrics$pct.peaks_linked_to_peaks <- seu_metadata$pct.peaks_linked_to_peaks[match(metrics$barcode, seu_metadata$barcode, nomatch = NA)]
  
  metrics$pct.peaks_linked_to_genes[is.nan(metrics$pct.peaks_linked_to_genes)] <- 0
  metrics$pct.genes_linked_to_peaks[is.nan(metrics$pct.genes_linked_to_peaks)] <- 0
  metrics$pct.peaks_linked_to_peaks[is.nan(metrics$pct.peaks_linked_to_peaks)] <- 0
  
  
  metrics <- metrics %>% dplyr::filter(barcode %in% rownames(seu@meta.data))
  
  
  # 6. Merge raw counts + compute additional percentages
  md <- seu@meta.data
  seu$gex_raw_reads  <- metrics$gex_raw_reads[ match(metrics$barcode, rownames(md)) ]
  seu$atac_raw_reads <- metrics$atac_raw_reads[ match(metrics$barcode, rownames(md)) ]
  
  metrics <- metrics %>% mutate(
    pct.fragments_in_peaks = atac_peak_region_fragments / atac_fragments * 100,
    pct.gex_mapped         = gex_mapped_reads    / gex_raw_reads    * 100,
    gex_unmapped_reads     = gex_raw_reads - gex_mapped_reads,
    pct.gex_unmapped       = gex_unmapped_reads   / gex_raw_reads    * 100,
    atac_mapped_reads      = atac_raw_reads - atac_unmapped_reads,
    pct.atac_mapped        = atac_mapped_reads    / atac_raw_reads   * 100,
    pct.atac_unmapped      = atac_unmapped_reads  / atac_raw_reads   * 100,
    pct.mito_atac          = atac_mitochondrial_reads / atac_mapped_reads * 100
  )
  
  metrics$barcode_id <- seu$barcode_id[match(seu$barcode, metrics$barcode)]
  metrics$pct.peaks_linked_to_genes <- seu$pct.peaks_linked_to_genes[match(seu$barcode, metrics$barcode)]
  
  metrics$pct.genes_linked_to_peaks <- seu$pct.genes_linked_to_peaks[match(seu$barcode, metrics$barcode)]
  
  metrics$pct.peaks_linked_to_peaks <- seu$pct.peaks_linked_to_peaks[match(seu$barcode, metrics$barcode)]
  
  seu$pct.fragments_in_peaks <- metrics$pct.fragments_in_peaks[ match(metrics$barcode, rownames(md)) ]
  seu$pct.gex_mapped        <- metrics$pct.gex_mapped[        match(metrics$barcode, rownames(md)) ]
  seu$pct.gex_unmapped      <- metrics$pct.gex_unmapped[      match(metrics$barcode, rownames(md)) ]
  seu$pct.atac_mapped       <- metrics$pct.atac_mapped[       match(metrics$barcode, rownames(md)) ]
  seu$pct.atac_unmapped     <- metrics$pct.atac_unmapped[     match(metrics$barcode, rownames(md)) ]
  seu$pct.mito_atac         <- metrics$pct.mito_atac[         match(metrics$barcode, rownames(md)) ]
  
  gc(); gc()
  
  seu_metadata <- seu@meta.data
  
  assign(paste0(sample, "_metrics"), metrics)
  assign(paste0(sample, "_linkages"), seu_linkages)
  
  message("✅ All samples done, with QC plots saved.") 
  
}


nemo_linkages <- rbind(s1_linkages, s2_linkages, s3_linkages, s4_linkages, s5_linkages, s6_linkages)

View(nemo_linkages)

nemo@misc$gene_peak_linkages <- nemo_linkages

saveRDS(nemo_linkages, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/nemo_gene_peak_linkages_df.rds')

nemo_metrics <- rbind(s1_metrics, s2_metrics, s3_metrics, s4_metrics, s5_metrics, s6_metrics)
View(nemo_metrics)

meta <- nemo@meta.data

meta$pct.peaks_linked_to_genes <- nemo$pct.peaks_linked_to_genes[match(meta$barcode_id, nemo$barcode_id)]

meta$pct.genes_linked_to_peaks <- nemo$pct.genes_linked_to_peaks[match(meta$barcode_id, nemo$barcode_id)]

meta$pct.peaks_linked_to_peaks <- nemo$pct.peaks_linked_to_peaks[match(meta$barcode_id, nemo$barcode_id)]

meta$gex_raw_reads <- nemo_metrics$gex_raw_reads[match(meta$barcode_id, nemo_metrics$barcode_id)]

meta$atac_raw_reads <- nemo_metrics$atac_raw_reads[match(meta$barcode_id, nemo_metrics$barcode_id)]

meta$pct.fragments_in_peaks <- nemo_metrics$pct.fragments_in_peaks[match(meta$barcode_id, nemo_metrics$barcode_id)]

meta$pct.gex_mapped <- nemo_metrics$pct.gex_mapped[match(meta$barcode_id, nemo_metrics$barcode_id)]

meta$pct.gex_unmapped <- nemo_metrics$pct.gex_unmapped[match(meta$barcode_id, nemo_metrics$barcode_id)]

meta$pct.atac_mapped <- nemo_metrics$pct.atac_mapped[match(meta$barcode_id, nemo_metrics$barcode_id)]

meta$pct.atac_unmapped <- nemo_metrics$pct.atac_unmapped[match(meta$barcode_id, nemo_metrics$barcode_id)]

meta$pct.mito_atac <- nemo$pct.mito_atac[match(meta$barcode_id, nemo_metrics$barcode_id)]


nemo@meta.data <- meta

nemo@misc$nemo_barcode_metrics <- nemo_metrics

saveRDS(metrics, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/nemo_cell_barcode_summary_metrics.rds')

meta2 <- meta[-(76:225)]
meta3 <- meta[(76:225)]
saveRDS(meta3, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/optimal_cluster_parameters_dataframe.rds')

meta3 <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/optimal_cluster_parameters_dataframe.rds')

nemo$optimal_clust_res0.8_40nn_45PC_40LSI <- meta3$res0.8_40nn_45PC_40LSI

nemo@meta.data <- meta2

saveRDS(nemo, '/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/Seurat_Data/seurat_objects/nemo.orig_harmony.integration_optimal_clusters_v2.rds')




