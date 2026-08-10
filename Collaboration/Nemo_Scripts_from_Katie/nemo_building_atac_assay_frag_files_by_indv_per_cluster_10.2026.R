library(future)
library(furrr)
library(parallel)
library(parallelly)

# Set workers to 75–80% of physical cores
workers <- max(2, floor(24 * 0.75))  # → 8 workers
# plan("multisession", workers = workers)
options(future.globals.maxSize = 50000 * 1024^2)

RunParallel <- function(
    input_list,
    FUN, # FUN = function(clust) is the function to apply to each element in input_list.
    ...,
    method = c("furrr", "bplapply"),
    workers = parallel::detectCores(logical = FALSE),
    progress = TRUE,
    bpparam = NULL
) {
  method <- match.arg(method)
  
  if (method == "furrr") {
    # Use multisession plan for RStudio / Mac safety
    future::plan("multisession", workers = workers)
    if (!requireNamespace("furrr", quietly = TRUE)) install.packages("furrr")
    library(furrr)
    library(progressr)
    if (progress) handlers(global = TRUE)
    
    result <- future_map(
      input_list,
      .f = FUN,
      ...,
      .options = furrr_options(seed = TRUE),
      .progress = progress
    )
    future::plan("sequential") # Reset plan
    return(result)
    
  } else if (method == "bplapply") {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) BiocManager::install("BiocParallel")
    library(BiocParallel)
    if (is.null(bpparam)) {
      bpparam <- BiocParallel::SnowParam(workers = workers, type = "SOCK", progressbar = progress)
    }
    return(BiocParallel::bplapply(input_list, FUN, ..., BPPARAM = bpparam))
  }
}

library(Signac)
library(igraph)
library(Seurat)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(scales)
library(scCustomize)
library(patchwork)
library(readr)
library(glue)
library(Biostrings)
library(biomartr)
library(tools)
library(rtracklayer)
library(Rsamtools)
library(BiocParallel)
library(fastmatch)
library(stringr)
library(utils)
library(GenomeInfoDb)
library(data.table)
library(DropletUtils)
library(rhdf5)

mem.maxVSize(30 * 1024^3)
options(future.globals.maxSize = 30 * 1024^3)

nemo <- readRDS("~/scratch/nemo_multi/seurat/seurat_objects/nemo_atac_fixed_10.22.25.rds")

annotations_works <- Annotation(nemo[['peaks']])
ann_df <- as.data.frame(annotations_works)
View(ann_df)

# Load genome fasta
nemo_genome_file <- '/storage/home/hcoda1/8/kleatherbury3/scratch/nemo_multi/nemo_ref_genome/Nemo/fasta/genome.fa'
nemo_genome <- read_genome(format = 'fasta', file = nemo_genome_file)

# nemo_genome is your DNAStringSet
seqnames    <- names(nemo_genome)
seqlengths  <- width(nemo_genome)
# if you know any scaffolds are circular, set TRUE; otherwise FALSE (or leave as NA)
is_circular <- rep(FALSE, length(seqnames))

# create the Seqinfo
nemo_seqinfo <- Seqinfo(
  seqnames    = seqnames,
  seqlengths  = seqlengths,
  isCircular  = is_circular,
  genome      = "Amphiprion_ocellaris"    # or whatever genome build you want to label
)

# annotations <- readRDS("~/scratch/nemo_multi/seurat/genome_files/nemo_genome_annotaton_coverageplot.rds")
# ann_df <- as.data.frame(annotations)
# View(ann_df)

nemo_seqinfo_works <- readRDS("~/scratch/nemo_multi/seurat/genome_files/nemo_seqinfo_coverageplot.rds")


##############################################################


gc()

colnames(nemo@meta.data)[colnames(nemo@meta.data) == 'barcode_orig'] = 'CellRanger_barcode'
colnames(nemo@meta.data)[colnames(nemo@meta.data) == 'barcode_id'] = 'individual_barcode'
nemo@meta.data$pool_barcode <- paste0(nemo@meta.data$pool, '_', nemo@meta.data$CellRanger_barcode)
head(nemo$CellRanger_barcode)
head(nemo$individual_barcode)
head(nemo$pool_barcode)


nemo <- RenameCells(nemo, new.names = paste0(nemo$pool, '_', nemo$CellRanger_barcode))

# Define pools
pool_list <- unique(nemo$pool)

main_frag_dir <- "/storage/home/hcoda1/8/kleatherbury3/scratch/nemo_multi/cellranger_seurat_files"
main_output_dir <- "/storage/home/hcoda1/8/kleatherbury3/scratch/nemo_multi/seurat/peakCalling_macs2_v2/renamedCellRangerBarcodeFeatureFiles"

RunParallel(
  input_list = pool_list,
  FUN = function(pool) {
    
    fragments_dir <- file.path(main_frag_dir, pool)
    output_dir <- file.path(main_output_dir, pool)
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    frag_file <- file.path(fragments_dir, "atac_fragments.tsv.gz")
    renamed_frag <- file.path(output_dir, paste0(pool, "_fragments_renamed.tsv"))
    renamed_frag_gz <- paste0(renamed_frag, ".gz")
    
    # Ensure working in output dir and clean
    setwd(output_dir)
    if (file.exists(renamed_frag_gz)) file.remove(renamed_frag_gz)
    
    # Rewrite using awk and bgzip
    awk_cmd <- glue::glue(
      "gunzip -c \"{frag_file}\" | ",
      "awk '{{ $4 = \"{pool}_\" $4; print }}' OFS='\t' > \"{basename(renamed_frag)}\""
    )
    stopifnot(system(awk_cmd) == 0)
    stopifnot(system(paste("bgzip -f", basename(renamed_frag))) == 0)
    stopifnot(system(paste("tabix -p bed", basename(renamed_frag_gz))) == 0)
    
    # Confirm valid bgzip
    message(glue("✅ Fragment file created and indexed: {renamed_frag_gz}"))
    
    # system(glue::glue("gunzip -c '{renamed_frag_gz}' | head -n 5"))
    
    ## --- Rename barcodes in .h5 ---
    h5_file <- normalizePath(file.path(fragments_dir, "raw_feature_bc_matrix.h5"))
    if (!file.exists(h5_file)) stop(glue("❌ Missing .h5 file for {pool}: {h5_file}"))
    file.copy(from = h5_file, to = output_dir, overwrite = TRUE)
    
    data <- Read10X_h5(file.path(output_dir, paste0('raw_feature_bc_matrix.h5')))
    
    # Gene expression
    genes.data <- data$`Gene Expression`
    colnames(genes.data) <- paste0(pool, "_", colnames(genes.data))
    output_h5_dir_gex <- normalizePath(file.path(output_dir, paste0(pool, "_GEX_raw_feature_bc_matrix_renamed.h5")))
    if (dir.exists(output_h5_dir_gex)) unlink(output_h5_dir_gex, recursive = TRUE)
    write10xCounts(
      path = output_h5_dir_gex,
      x = genes.data,
      type = "HDF5",
      version = "3",
      overwrite = TRUE
    )
    
    # ATAC
    atac.data <- data$Peaks
    colnames(atac.data) <- paste0(pool, "_", colnames(atac.data))
    output_h5_file_atac <- normalizePath(file.path(output_dir, paste0(pool, "_ATAC_raw_feature_bc_matrix_renamed.h5")), mustWork = FALSE)
    if (file.exists(output_h5_file_atac)) file.remove(output_h5_file_atac)
    
    write10xCounts(
      path = output_h5_file_atac,
      x = atac.data,
      type = "HDF5",
      version = "3",
      overwrite = TRUE
    )
    message(glue("✅ Saved renamed .h5 for {pool} at {output_h5_dir_gex} and {output_h5_file_atac}"))
    message(glue("✅ Completed renaming for pool: {pool}"))
    
  },
  method = "furrr",
  workers = 4
)

gc()
gc()

annotations <- nemo_poa@assays[["peaks"]]@annotation
nemo_seqinfo <- nemo_poa@assays[["peaks"]]@seqinfo
nemo$seq_pool <- nemo$pool
seq_pool_list <- unique(nemo$seq_pool)

allseq_pool_objectList <- list()

for (seq_pool in seq_pool_list) {
  
  main_dir = "/storage/home/hcoda1/8/kleatherbury3/scratch/nemo_multi/seurat/peakCalling_macs2_v2"
  
  meta <- nemo@meta.data
  
  pool_list <- meta$pool[meta$seq_pool == seq_pool]
  pool_list <- unique(pool_list)
  pool_list
  
  pool_ranges <- list()
  pool_objects <- list()
  median_frag_lengths <- list()
  
  for (pool in pool_list) {
    
    pool_dir <- file.path(main_dir, 'renamedCellRangerBarcodeFeatureFiles', pool)
    fragFile_suffix <- "_fragments_renamed.tsv.gz"
    fragPath <- file.path(pool_dir, paste0(pool, fragFile_suffix))
    
    meta <- nemo@meta.data
    individual_list <- unique(meta$individual[meta$seq_pool == seq_pool])
    meta <- meta %>% dplyr::filter(individual %in% individual_list)
    
    cells_to_keep <- rownames(meta)
    meta$nCount_RNA <- NULL
    meta$nFeature_RNA <- NULL
    meta$nCount_ATAC <- NULL
    meta$nFeature_ATAC <- NULL
    
    frag_obj <- CreateFragmentObject(path = fragPath, validate.fragments = TRUE)
    
    seu <- Read10X_h5(file.path(pool_dir, paste0(pool, '_GEX_raw_feature_bc_matrix_renamed.h5')))
    seu <- CreateSeuratObject(seu)
    
    atac.counts <- Read10X_h5(file.path(pool_dir, paste0(pool, '_ATAC_raw_feature_bc_matrix_renamed.h5')))
    seu[['ATAC']] <- CreateChromatinAssay(
      counts = atac.counts,
      fragments = frag_obj,
      annotation = annotations,
      genome = nemo_seqinfo,
      validate.fragments = TRUE,
      sep = c(":", "-")
    )
    
    seu <- subset(seu, cells = cells_to_keep)
    seu <- AddMetaData(seu, metadata = meta)
    
    seqlevels(seu[['ATAC']]@ranges) <- seqlevels(nemo_seqinfo)
    seqinfo(seu[['ATAC']]@ranges) <- nemo_seqinfo
    
    seu <- RenameCells(seu, new.names = paste0(seu$individual, '_', seu$CellRanger_barcode))
    seu_split <- SplitObject(seu, split.by = 'individual')
    individual_list <- names(seu_split)
    
    safe_fragment_qc <- function(seu, pool, main_dir, seq_pool, nemo_seqinfo) {
      individual <- unique(seu$individual)
      
      if (length(individual) != 1 || is.na(individual) || individual == "") {
        warning(glue::glue("⚠️ Skipping: invalid individual for pool {pool}"))
        return(NULL)
      }
      
      indv_pool <- unique(seu$pool)
      indv_seq_pool <- unique(seu$seq_pool)
      
      mcols(seu[['ATAC']]@ranges)$individual <- individual
      mcols(seu[['ATAC']]@ranges)$pool <- indv_pool
      mcols(seu[['ATAC']]@ranges)$seq_pool <- indv_seq_pool
      mcols(seu[['ATAC']]@ranges)$peak_id <- paste0(individual, '_', indv_seq_pool, '_', indv_pool, "_peak_", seq_along(seu[['ATAC']]@ranges))
      names(seu[['ATAC']]@ranges) <- mcols(seu[['ATAC']]@ranges)$peak_id
      
      fragPath <- normalizePath(seu[["ATAC"]]@fragments[[1]]@path)
      frags <- data.table::fread(
        cmd = glue::glue("gunzip -c '{fragPath}' | grep -v '^#' | head -n 5000000"),
        col.names = c("chr", "start", "end", "barcode", "readcount"),
        sep = "\t",
        fill = TRUE
      )
      frags[, length := end - start]
      median_frag_length <- median(frags$length)
      
      binned <- frags %>%
        dplyr::filter(length <= 1000) %>%
        dplyr::mutate(bin = floor(length / 20) * 20) %>%
        dplyr::count(bin)
      
      plot <- ggplot(binned, aes(x = bin, y = n)) +
        geom_col(width = 20, fill = "steelblue", color = "black", size = 0.8) +
        geom_rect(aes(xmin = 0, xmax = 130, ymin = 0, ymax = Inf), fill = "palegreen3", alpha = 0.002, inherit.aes = FALSE) +
        annotate("text", x = 90, y = Inf, label = "NFR\n(<125 bp)", vjust = 1.5, size = 4) +
        geom_rect(aes(xmin = 130, xmax = 275, ymin = 0, ymax = Inf), fill = "orange", alpha = 0.002, inherit.aes = FALSE) +
        annotate("text", x = 202.5, y = Inf, label = "Mononucleosome\n(130–250 bp)", vjust = 1.5, size = 4) +
        geom_rect(aes(xmin = 275, xmax = 1000, ymin = 0, ymax = Inf), fill = "orangered", alpha = 0.002, inherit.aes = FALSE) +
        annotate("text", x = 600, y = Inf, label = "Dinucleosome\n(300–450 bp)", vjust = 1.5, size = 4) +
        scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 1000, by = 100)) +
        theme_classic() +
        labs(
          title = paste0('Pool ', pool, ": Fragment Length Distribution \n"),
          x = "Fragment Length (bp) \n",
          y = "\n Count"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 15, vjust = -1.5, color = 'navy'),
          axis.title.x = element_text(face = "bold", size = 12, vjust = -1, color = 'black'),
          axis.title.y = element_text(face = "bold", size = 12, vjust = 1, color = 'black')
        )
      
      seq_pool_dir <- file.path(main_dir, "peakCallingPerseq_pool", seq_pool)
      indv_dir <- file.path(seq_pool_dir, individual)
      if (!dir.exists(indv_dir)) dir.create(indv_dir, recursive = TRUE)
      
      ggsave(
        filename = file.path(indv_dir, paste0(individual, "_CellRanger_fragLengths.png")),
        plot = plot,
        height = 8, width = 15
      )
      
      indv_rna_counts <- GetAssayData(seu, assay = 'RNA', layer = 'counts')
      indv_atac_counts <- GetAssayData(seu, assay = 'ATAC', layer = 'counts')
      seu_barcodes <- seu@assays[["ATAC"]]@fragments[[1]]@cells
      
      # Save barcode list (as .txt for readability)
      writeLines(seu_barcodes, con = file.path(indv_dir, paste0(individual, "_barcode_list.txt")))
      
      # Save objects to disk, without overwriting the in-memory objects
      saveRDS(seu, file = file.path(indv_dir, paste0(individual, "_CellRanger_SeuratObject.rds")))
      saveRDS(seu[['ATAC']]@ranges, file = file.path(indv_dir, paste0(individual, "_CellRanger_GRanges.rds")))
      
      return(list(
        seu = seu,
        ranges = seu[['ATAC']]@ranges,
        median_frag_length = median_frag_length
      ))
    }
    
    # Run in parallel
    results <- RunParallel(
      input_list = seu_split,
      FUN = function(seu) safe_fragment_qc(
        seu = seu,
        pool = pool,
        main_dir = main_dir,
        seq_pool = seq_pool,
        nemo_seqinfo = nemo_seqinfo
      ),
      method = "furrr",
      workers = 4
    )
    
    # Remove skipped entries
    results <- purrr::compact(results)
    
    # Split into separate outputs
    seu_result_list <- purrr::map(results, "seu")
    ranges_result_list <- purrr::map(results, "ranges")
    frag_length_result <- purrr::map_dbl(results, "median_frag_length")
    
    individual_ids <- purrr::map_chr(seu_result_list, ~ unique(.x$individual))
    names(seu_result_list) <- individual_ids
    names(ranges_result_list) <- individual_ids
    names(frag_length_result) <- individual_ids
    
    pool_objects[[pool]] <- seu_result_list
    pool_ranges[[pool]] <- ranges_result_list
    median_frag_lengths[[pool]] <- frag_length_result
    
    rm(results)
    gc()
    
  }
  
  seq_pool_result <- list(
    seq_pool_objects = purrr::flatten(pool_objects),
    seq_pool_ranges = purrr::flatten(pool_ranges),
    seq_pool_medianFragLengths = purrr::flatten(median_frag_lengths)
  )
  
  rm(pool_objects)
  rm(pool_ranges)
  rm(median_frag_lengths)
  
  gc()
  
  allseq_pool_objectList[[seq_pool]] <- seq_pool_result
  
}


saveRDS(allseq_pool_objectList, file.path(main_peak_dir, 'allseq_pool_objectList.rds'))


future::plan("sequential")  # reset parallel backend before calling again


reduced_rangesPerseq_pool <- list()
all_rangesPerseq_pool <- list()
indvObjectsPerseq_pool <- list()

gc()
gc()

for (seq_pool in names(allseq_pool_objectList)) {
  
  seq_pool_dir <- file.path(main_dir, "peakCallingPerseq_pool", seq_pool)
  if (!dir.exists(seq_pool_dir)) dir.create(seq_pool_dir, recursive = TRUE)
  
  pool_ranges <- allseq_pool_objectList[[seq_pool]][['seq_pool_ranges']]
  pool_objects <- allseq_pool_objectList[[seq_pool]][['seq_pool_objects']]
  
  all_ranges <- unlist(as(c(pool_ranges[1:length(pool_ranges)]), "GRangesList"))
  all_ranges <- trim(all_ranges)
  
  reduce_peaks <- GenomicRanges::reduce(
    unlist(as(c(pool_ranges[1:length(pool_ranges)]), "GRangesList")),
    with.revmap = TRUE, min.gapwidth = 0L
  )
  reduce_peaks <- reduce_peaks[!is.na(seqnames(reduce_peaks))]
  reduce_peaks <- trim(reduce_peaks)
  
  reduced_rangesPerseq_pool[[seq_pool]] <- reduce_peaks
  all_rangesPerseq_pool[[seq_pool]] <- all_ranges
  
  saveRDS(reduce_peaks, file.path(seq_pool_dir, paste0(seq_pool, "_CellRanger_reducedPeaks.rds")))
  saveRDS(all_ranges, file.path(seq_pool_dir, paste0(seq_pool, "_CellRanger_allPeaksMerged.rds")))
  
  individual_list <- names(pool_objects)
  indv_objects <- list()
  
  for (individual in individual_list) {
    seu <- pool_objects[[individual]]
    
    peak_counts <- FeatureMatrix(
      fragments = Fragments(seu[['ATAC']]), 
      features = reduce_peaks, 
      cells = colnames(seu)
    )
    
    seu[['ATAC']] <- CreateChromatinAssay(
      counts = peak_counts,
      fragments = Fragments(seu[['ATAC']]),
      annotation = annotations,
      genome = nemo_seqinfo,
      validate.fragments = TRUE
    )
    
    gc()
    seqlevels(seu[['ATAC']]@ranges) <- seqlevels(nemo_seqinfo)
    seqinfo(seu[['ATAC']]@ranges) <- nemo_seqinfo
    seu[['ATAC']]@ranges <- trim(seu[['ATAC']]@ranges)
    
    indv_dir <- file.path(seq_pool_dir, individual)
    if (!dir.exists(indv_dir)) dir.create(indv_dir, recursive = TRUE)
    saveRDS(seu, file.path(indv_dir, paste0(individual, "_CellRanger_reducedPeaks_seuratObject.rds")))
    
    indv_objects[[individual]] <- seu
  }
  
  # Merge individual Seurat objects into a seq_pool-level object
  if (length(indv_objects) > 1) {
    seq_pool_object <- merge(
      x = indv_objects[[1]], 
      y = indv_objects[2:length(indv_objects)], 
      merge.data = TRUE
    )
  } else {
    seq_pool_object <- indv_objects[[1]]
  }
  
  saveRDS(seq_pool_object, file.path(seq_pool_dir, paste0(seq_pool, '_mergedObject.rds')))
  
  indvObjectsPerseq_pool[[seq_pool]] <- seq_pool_object
  
}


gc()
gc()

saveRDS(indvObjectsPerseq_pool, file.path(main_dir, "indvObjectsPerseq_pool.rds"))


ObjectsPerGroup <- list()
group_list <- names(indvObjectsPerseq_pool)
peak_list <- list()

main_peak_dir <- "/storage/home/hcoda1/8/kleatherbury3/scratch/nemo_multi/seurat/peakCalling_macs2_v2"

macs3_path <- "/storage/home/hcoda1/8/kleatherbury3/r-js585-0/opt/miniconda3/envs/macs_env/bin/macs3"
py_path    <- "/storage/home/hcoda1/8/kleatherbury3/r-js585-0/opt/miniconda3/envs/macs_env/bin/python3.10"

stopifnot(file.exists(macs3_path))
stopifnot(file.exists(py_path))
system(paste(shQuote(py_path), shQuote(macs3_path), "--version"))

macs3_wrapper <- file.path(main_peak_dir, "macs3_wrapper.sh")
writeLines(
  c(
    "#!/usr/bin/env bash",
    paste(shQuote(py_path), shQuote(macs3_path), '"$@"')
  ),
  con = macs3_wrapper
)
Sys.chmod(macs3_wrapper, mode = "0755")
stopifnot(file.exists(macs3_wrapper))
stopifnot(file.access(macs3_wrapper, mode = 1) == 0)

for (group in group_list) {
  
  group_dir <- file.path(main_peak_dir, "peakCallingPerseq_pool", group)
  group_object <- indvObjectsPerseq_pool[[group]]
  
  peak_dir <- file.path(group_dir, "clusters")
  if (!dir.exists(peak_dir)) dir.create(peak_dir, recursive = TRUE)
  
  additional.args <- c("--slocal 500 --llocal 3000 --max-gap 25 --keep-dup all --tsize 150 --min-length 100 --bw 200")
  
  peaks <- CallPeaks(
    object = group_object,
    macs2.path = macs3_wrapper,
    effective.genome.size = 8.5e8,
    outdir = peak_dir,
    cleanup = FALSE,
    assay = "ATAC",
    fragment.tempdir = peak_dir,
    group.by = "final_clusters",
    name = gsub(" ", "_", Project(group_object)),
    combine.peaks = TRUE,
    nomodel = TRUE,
    qvalue = 0.05,
    call_summits = TRUE,
    extsize = 200,
    shift = -100,
    additional.args = additional.args
  )
  
  saveRDS(peaks, file.path(group_dir, paste0(group, "_MACS2_allPeaksPerGroup_preReduce.rds")))
  
  peak_list[[group]] <- peaks
  ObjectsPerGroup[[group]] <- group_object
}

saveRDS(peak_list, file.path(main_peak_dir, 'macs2_allPeaksMerged_preReduce.rds'))
saveRDS(ObjectsPerGroup, file.path(main_peak_dir, 'macs2_ObjectsPerGroup_preReduce.rds'))

for (pool in names(peak_list)){
  granges <- peak_list[[pool]]
  granges@seqinfo <- nemo_seqinfo
  peak_list[[pool]] <- granges
}
all_ranges <- unlist(as(c(peak_list[1:length(peak_list)]), "GRangesList"))

## test ranges
# Idents(nemo) <- 'final_clusters'
# DefaultAssay(nemo) <- 'ATAC'
# 
# a <- CoveragePlot(
#   object = nemo,
#   region = "gfap",        
#   features = "gfap",
#   extend.upstream = 5000,
#   extend.downstream = 5000,
#   # expression.assay = 'RNA',
#   assay = 'ATAC',
#   sep = c('-', '-'),
#   peaks = T,
#   annotation = T,
#   ymax = 'q95',
#   ranges = all_ranges,
#   links = F
# )
# 
# b <- CoveragePlot(
#   object = nemo,
#   region = "cga",        
#   features = "cga",
#   extend.upstream = 5000,
#   extend.downstream = 5000,
#   # expression.assay = 'RNA',
#   assay = 'ATAC',
#   sep = c('-', '-'),
#   peaks = T,
#   annotation = T,
#   ymax = 'q95',
#   ranges = all_ranges,
#   links = F
# )
# a + b

all_ranges <- trim(all_ranges)

macs2_reduce_peaks <- GenomicRanges::reduce(all_ranges, with.revmap = TRUE)
macs2_reduce_peaks <- macs2_reduce_peaks[!is.na(seqnames(macs2_reduce_peaks))]
mcols(macs2_reduce_peaks)$peak_id <- paste0("peak_", seq_along(macs2_reduce_peaks))

cr_gr_list <- list()
for (group in group_list) {
  cr_gr <- granges(indvObjectsPerseq_pool[[group]][["ATAC"]])
  cr_gr@elementMetadata@listData[["group_peak_ids"]] <- NULL
  cr_gr_list[[group]] <- cr_gr
}
cr_gr <- unlist(as(c(cr_gr_list[1:length(cr_gr_list)]), "GRangesList"))
cr_df <- as.data.frame(mcols(cr_gr))

macs_gr <- macs2_reduce_peaks

gc()

############



# --- pick likely column names directly (no helpers) ---
nm <- tolower(names(cr_df))
col_indv <- if (any(nm %in% c("individual","individuals","indv","sample_id","donor"))) names(cr_df)[which(nm %in% c("individual","individuals","indv","sample_id","donor"))[1]] else NA_character_
col_pool <- if (any(nm %in% c("pool","pool_id","pools")))                                           names(cr_df)[which(nm %in% c("pool","pool_id","pools"))[1]]                 else NA_character_
col_Status  <- if (any(nm %in% c("Status","Status","seq_pool")))                                 names(cr_df)[which(nm %in% c("Status","Status","seq_pool"))[1]]       else NA_character_
col_indv_pid  <- if (any(nm %in% c('indv_peak_ids','indv_peak_id')))                                     names(cr_df)[which(nm %in% c('indv_peak_ids','indv_peak_id'))[1]]           else NA_character_
col_gpid <- if (any(nm %in% c("Status_peak_id","Status_peak_ids")))                        names(cr_df)[which(nm %in% c("Status_peak_id","Status_peak_ids"))[1]] else NA_character_

INDV <- if (!is.na(col_indv)) as.character(cr_df[[col_indv]]) else rep(NA_character_, nrow(cr_df))
POOL <- if (!is.na(col_pool)) as.character(cr_df[[col_pool]]) else rep(NA_character_, nrow(cr_df))
Status  <- if (!is.na(col_Status )) as.character(cr_df[[col_Status ]]) else rep(NA_character_, nrow(cr_df))
PID  <- if (!is.na(col_indv_pid )) as.character(cr_df[[col_indv_pid ]]) else rep(NA_character_, nrow(cr_df))
GPID <- if (!is.na(col_gpid)) as.character(cr_df[[col_gpid]]) else rep(NA_character_, nrow(cr_df))

# --- overlaps (any bp; simple) ---
hits <- findOverlaps(macs_gr, cr_gr, ignore.strand = TRUE)

# initialize columns on macs_gr
n <- length(macs_gr)
mcols(macs_gr)$individuals        <- rep("", n)
mcols(macs_gr)$pools              <- rep("", n)
mcols(macs_gr)$Status             <- rep("", n)
mcols(macs_gr)$cr_peak_ids        <- rep("", n)
mcols(macs_gr)$cr_Status_peak_ids  <- rep("", n)
mcols(macs_gr)$num_unique_indvs   <- integer(n)
mcols(macs_gr)$num_unique_pools   <- integer(n)
mcols(macs_gr)$num_unique_Status  <- integer(n)
mcols(macs_gr)$num_cr_peaks       <- integer(n)


if (length(hits) > 0L) {
  DT <- data.table(
    macs_idx = as.integer(queryHits(hits)),
    cr_idx   = as.integer(subjectHits(hits)),
    indv     = INDV[subjectHits(hits)],
    pool     = POOL[subjectHits(hits)],
    status      = Status [subjectHits(hits)],
    pid      = PID [subjectHits(hits)],
    gpid     = GPID[subjectHits(hits)]
  )
  
  # collapse unique values per MACS peak (like your barcode example)
  DT[, indv := as.character(indv)]
  DT[, pool := as.character(pool)]
  DT[, status  := as.character(status )]
  DT[, pid  := as.character(pid )]
  DT[, gpid := as.character(gpid)]
  
  agg <- DT[, .(
    individuals        = paste(sort(unique(indv [!is.na(indv ) & nzchar(indv )])), collapse = ","),
    pools              = paste(sort(unique(pool [!is.na(pool ) & nzchar(pool )])), collapse = ","),
    Status             = paste(sort(unique(status  [!is.na(status  ) & nzchar(status  )])), collapse = ","),
    cr_peak_ids        = paste(sort(unique(pid  [!is.na(pid  ) & nzchar(pid  )])), collapse = ","),
    cr_Status_peak_ids  = paste(sort(unique(gpid [!is.na(gpid ) & nzchar(gpid )])), collapse = ","),
    num_unique_indvs   = uniqueN(indv [nzchar(indv )]),
    num_unique_pools   = uniqueN(pool [nzchar(pool )]),
    num_unique_Status  = uniqueN(status  [nzchar(status  )]),
    num_cr_peaks       = .N
  ), by = macs_idx]
  
  mcols(macs_gr)$individuals       [agg$macs_idx] <- agg$individuals
  mcols(macs_gr)$pools             [agg$macs_idx] <- agg$pools
  mcols(macs_gr)$Status            [agg$macs_idx] <- agg$Status
  mcols(macs_gr)$cr_peak_ids       [agg$macs_idx] <- agg$cr_peak_ids
  mcols(macs_gr)$cr_Status_peak_ids [agg$macs_idx] <- agg$cr_Status_peak_ids
  mcols(macs_gr)$num_unique_indvs  [agg$macs_idx] <- agg$num_unique_indvs
  mcols(macs_gr)$num_unique_pools  [agg$macs_idx] <- agg$num_unique_pools
  mcols(macs_gr)$num_unique_Status [agg$macs_idx] <- agg$num_unique_Status
  mcols(macs_gr)$num_cr_peaks      [agg$macs_idx] <- agg$num_cr_peaks
}

allPeaks_info <- as.data.frame(macs_gr)
write_csv(allPeaks_info, file.path(main_peak_dir, paste0("MACS2_peakInfoDF_reducedPeaks.csv")))

gc()

macs2_reduce_peaks <- macs_gr

saveRDS(macs2_reduce_peaks, file.path(main_peak_dir, 'macs2_allGroupsMerged_reducedPeaks.rds'))


#############################

group_list <- names(ObjectsPerGroup)
final_peaks <- list()

#### MAKE SURE TO DOUBLE CHECK CORRECT ANNOTATION AND SEQINFO OBJECT ARE SET PRIOR TO RUNNING
for (group in group_list) {
  group_dir <- file.path(main_peak_dir, "peakCallingPerseq_pool", group)
  seu <- ObjectsPerGroup[[group]]
  
  peak_counts <- FeatureMatrix(
    fragments = Fragments(seu[['ATAC']]),
    features = macs2_reduce_peaks,
    cells = colnames(seu)
  )
  
  seu[['peaks']] <- CreateChromatinAssay(
    counts = peak_counts,
    fragments = Fragments(seu[['ATAC']]),
    annotation = annotations,
    genome = nemo_seqinfo,
    validate.fragments = TRUE
  )
  
  gc()
  
  saveRDS(
    seu,
    paste0(group_dir, '/', group, "_gex.atac.peaks_seuratObject.rds")
  )
  
  seu[['ATAC']] <- NULL
  final_peaks[[group]] <- seu
}

# Merge individual grouprat objects into a group-level object
if (length(final_peaks) > 1) {
  final_object <- merge(
    x = final_peaks[[1]], 
    y = final_peaks[2:length(final_peaks)], 
    merge.data = TRUE
  )
} else {
  final_object <- final_peaks[[1]]
}

final_object[['RNA']] <- JoinLayers(final_object[['RNA']])

saveRDS(final_object, file.path(main_peak_dir, paste0('nemo_gex.atac.peaks.mergedObject.rds')))
saveRDS(final_peaks, file.path(main_peak_dir, paste0('nemo_final_reduced_peaks_list_per_group.rds')))

gc()


#################################

# nemo_orig <- nemo
final_object <- readRDS(file.path(main_peak_dir, paste0('nemo_gex.atac.peaks.mergedObject.rds')))

main_dir = "/storage/home/hcoda1/8/kleatherbury3/scratch/nemo_multi/seurat/peakCalling_macs2_v2"

gc()

nemo <- final_object
frags_list <- Fragments(nemo[['peaks']])
names(frags_list) <- as.character(seq_along(frags_list))
names(frags_list)

for (name in names(frags_list)) {
  barcode <- names(frags_list[[name]]@cells)[[1]]
  prefix <- sub("_.*", "", barcode)
  names(frags_list)[names(frags_list) == name] <- prefix
}

out_dir <- file.path(main_dir, 'fragment_files')
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# paths + names derived once (used by worker)
frag_paths <- vapply(frags_list, function(f) GetFragmentData(f, slot = "path"), character(1))
frag_names <- names(frags_list)

# Require external tools (run once, not in workers)
need <- c("gunzip", "bgzip", "tabix", "sort")
miss <- need[Sys.which(need) == ""]
if (length(miss)) stop("Missing tools on PATH: ", paste(miss, collapse = ", "))

## --- NEW: decide which indices actually need processing (missing folder OR incomplete/corrupt outputs) ---
.is_complete_and_ok <- function(frag_dir, frag_name) {
  out_gz  <- file.path(frag_dir, paste0(frag_name, ".tsv.gz"))
  out_tbi <- paste0(out_gz, ".tbi")
  
  # Must exist and be non-empty
  if (!file.exists(out_gz) || !file.exists(out_tbi)) return(FALSE)
  if (is.na(file.size(out_gz)) || file.size(out_gz) <= 0) return(FALSE)
  if (is.na(file.size(out_tbi)) || file.size(out_tbi) <= 0) return(FALSE)
  
  # Must be readable gzip
  ok_gz <- system(sprintf("gunzip -t %s", shQuote(out_gz))) == 0
  if (!ok_gz) return(FALSE)
  
  # Must be readable by tabix using its index
  ok_tabix <- system(sprintf("tabix -H %s > /dev/null 2>&1", shQuote(out_gz))) == 0
  if (!ok_tabix) return(FALSE)
  
  TRUE
}

to_run <- which(vapply(seq_along(frag_names), function(i) {
  frag_dir <- file.path(out_dir, frag_names[i])
  # If folder missing, run
  if (!dir.exists(frag_dir)) return(TRUE)
  # If outputs incomplete/corrupt, run
  !.is_complete_and_ok(frag_dir, frag_names[i])
}, logical(1)))

message("Total fragments: ", length(frag_names))
message("Will process (missing/incomplete/corrupt): ", length(to_run))
message("Will skip (complete): ", length(frag_names) - length(to_run))

## --- worker: EXACTLY your per-fragment loop body, scoped to one index i ---
.process_one_fragment <- function(i) {
  in_path <- frag_paths[i]
  
  # Make subfolder for this fragment
  frag_dir <- file.path(out_dir, frag_names[i])
  if (!dir.exists(frag_dir)) dir.create(frag_dir, recursive = TRUE)
  
  # Define output paths inside this folder
  out_gz  <- file.path(frag_dir, paste0(frag_names[i], ".tsv.gz"))
  out_tbi <- paste0(out_gz, ".tbi")
  
  message("Processing: ", in_path, "  ->  ", out_gz)
  
  # Temp files inside fragment folder
  tmp_header <- file.path(frag_dir, paste0(frag_names[i], "_header.txt"))
  tmp_data   <- file.path(frag_dir, paste0(frag_names[i], "_data.tsv"))
  tmp_sorted <- file.path(frag_dir, paste0(frag_names[i], "_sorted.tsv"))
  tmp_comb   <- file.path(frag_dir, paste0(frag_names[i], "_combined.tsv"))
  
  # 1) Extract header and data separanemoy
  stopifnot(system(sprintf("gunzip -c %s | grep '^#' > %s",
                           shQuote(in_path), shQuote(tmp_header))) == 0)
  stopifnot(system(sprintf("gunzip -c %s | grep -v '^#' > %s",
                           shQuote(in_path), shQuote(tmp_data))) == 0)
  
  # 2) Sort data by chr and start
  stopifnot(system(sprintf("LC_ALL=C sort -k1,1 -k2,2n %s > %s",
                           shQuote(tmp_data), shQuote(tmp_sorted))) == 0)
  
  # 3) Recombine header + sorted data, compress, and index
  stopifnot(system(sprintf("cat %s %s > %s",
                           shQuote(tmp_header), shQuote(tmp_sorted), shQuote(tmp_comb))) == 0)
  stopifnot(system(sprintf("bgzip -c %s > %s",
                           shQuote(tmp_comb), shQuote(out_gz))) == 0)
  stopifnot(system(sprintf("tabix -c '#' -s 1 -b 2 -e 3 %s",
                           shQuote(out_gz))) == 0)
  
  # 4) Clean temp files
  file.remove(tmp_header, tmp_data, tmp_sorted, tmp_comb)
  
  message("Done: ", out_gz, "  (index: ", out_tbi, ")")
  # return something useful if you want:
  data.frame(name = frag_names[i], input = in_path, output = out_gz, index = out_tbi, stringsAsFactors = FALSE)
}

## --- parallelize your original for-loop with RunParallel()
RunParallelPACE <- function(
    input_list,
    FUN,
    ...,
    workers    = NULL,   # NULL = use current future::plan(); otherwise override
    progress   = TRUE,
    globals    = TRUE,   # TRUE = auto-detect; or pass list(obj1=..., obj2=...)
    batch_size = NULL    # e.g., 6–8 for heavy I/O; NULL=no batching
) {
  if (!requireNamespace("furrr", quietly = TRUE)) stop("Package 'furrr' is required.")
  have_progressr <- requireNamespace("progressr", quietly = TRUE)
  
  # Discover effective workers from your current plan / options
  eff_from_plan  <- tryCatch(future::nbrOfWorkers(), error = function(e) 1L)
  eff_from_arg   <- if (is.null(workers)) eff_from_plan else as.integer(workers)
  eff_from_cores <- tryCatch(parallelly::availableCores(logical = FALSE), error = function(e) eff_from_arg)
  eff_workers    <- max(1L, min(eff_from_arg, eff_from_plan, eff_from_cores))
  
  # Runner builds the global set, including progressor when used
  runner <- function(ix) {
    if (isTRUE(progress) && have_progressr) {
      progressr::with_progress({
        p <- progressr::progressor(along = ix)
        globs <- if (isTRUE(globals)) TRUE else c(list(FUN = FUN, p = p), globals)
        
        furrr::future_map(
          input_list[ix],
          .f = function(x, ...) { p(sprintf("processing: %s", x)); FUN(x, ...) },
          ...,
          .options = furrr::furrr_options(seed = TRUE, globals = globs)
        )
      })
    } else {
      globs <- if (isTRUE(globals)) TRUE else c(list(FUN = FUN), globals)
      furrr::future_map(
        input_list[ix],
        .f = function(x, ...) FUN(x, ...),
        ...,
        .options = furrr::furrr_options(seed = TRUE, globals = globs)
      )
    }
  }
  
  # Respect your plan; only throttle via batching if requested
  if (is.null(batch_size) || batch_size <= 0L) {
    out <- runner(seq_along(input_list))
  } else {
    bs   <- as.integer(batch_size)
    bs   <- max(1L, min(bs, eff_workers))          # don’t exceed available workers
    idx  <- seq_along(input_list)
    bats <- split(idx, ceiling(seq_along(idx) / bs))
    out  <- lapply(bats, runner)
    out  <- unlist(out, recursive = FALSE, use.names = FALSE)
  }
  out
}

res_tbl_list <- RunParallelPACE(
  input_list = to_run,  # NEW: only run missing/incomplete/corrupt
  FUN        = .process_one_fragment,
  progress   = TRUE,
  globals    = list(  # make sure workers see these objects used inside FUN
    frag_paths = frag_paths,
    frag_names = frag_names,
    out_dir    = out_dir
  ),
  batch_size = 6      # good starting point for bgzip/tabix/sort-heavy jobs
)
# res_tbl <- do.call(rbind, res_tbl_list)


main_dir = "/storage/home/hcoda1/8/kleatherbury3/scratch/nemo_multi/seurat/peakCalling_macs2_v2"

# nemo <- readRDS("~/scratch/nemo_multi/seurat/peakCalling_macs2_v2/nemo_gex.atac.peaks.mergedObject.rds")

# check that frags have individual ids
names(nemo[['peaks']]@fragments)

frags_list <- Fragments(nemo[['peaks']])
names(frags_list) <- as.character(seq_along(frags_list))
names(frags_list)

names(nemo[['peaks']]@fragments) <- names(frags_list)

for (name in names(frags_list)) {
  barcode <- names(frags_list[[name]]@cells)[[1]]
  prefix <- sub("_.*", "", barcode)
  names(frags_list)[names(frags_list) == name] = prefix
}

names(frags_list)

for (name in names(frags_list)) {
  cells <- frags_list[[name]]@cells
  frag_dir <- file.path(main_dir, 'fragment_files', name, paste0(name, '.tsv.gz'))
  frag <- CreateFragmentObject(frag_dir, cells = cells, validate.fragments = TRUE)
  frags_list[[name]] <- frag
}

frag_paths <- list()
for (name in names(frags_list)) {
  path <- frags_list[[name]]@path
  frag_paths[[name]] <- path
}

nemo[['peaks']]@fragments <- frags_list

for (name in names(nemo[['peaks']]@fragments)) {
  print(nemo[['peaks']]@fragments[[name]]@path)
  print(head(nemo[['peaks']]@fragments[[name]]@cells))
}

total_fragments <- CountFragments(frag_paths)

# main_dir = "/storage/home/hcoda1/8/kleatherbury3/scratch/nemo_multi/seurat/peakCalling_macs2_v2"

saveRDS(frags_list, file.path(main_dir, 'nemo_fragmentObjectList.rds'))
saveRDS(total_fragments, file.path(main_dir, 'nemo_total_fragment_counts_per_individual.rds'))
saveRDS(nemo, '~/scratch/nemo_multi/seurat/seurat_objects/nemo_atac_new_03.04.26.rds')


######################

a <- CoveragePlot(
  object = nemo_orig,
  region = "gfap",        
  features = "gfap",
  extend.upstream = 5000,
  extend.downstream = 5000,
  expression.assay = 'RNA',
  assay = 'ATAC',
  sep = c('-', '-'),
  peaks = T,
  annotation = T,
  idents = 15,
  ymax = 'q95',
  ranges = nemo[['peaks']]@ranges,
  links = F
)
a

c <- CoveragePlot(
  object = nemo_orig,
  region = "cga",        
  features = "cga",
  extend.upstream = 5000,
  extend.downstream = 5000,
  expression.assay = 'RNA',
  assay = 'ATAC',
  sep = c('-', '-'),
  peaks = T,
  annotation = T,
  idents = 22,
  ymax = 'q95',
  ranges = nemo[['peaks']]@ranges,
  links = F
)
c

######################



nemo <- NucleosomeSignal(object = nemo, assay = 'peaks')

nemo_tss <- import('~/scratch/nemo_multi/nemo_ref_genome/Nemo/regions/tss.bed', format = "BED")

nemo <- TSSEnrichment(nemo, 
                      cells = colnames(nemo), 
                      assay = 'peaks', 
                      verbose = T, 
                      tss.positions = nemo_tss,
                      region_extension = 2000)

gc()

head(colnames(nemo))

nemo <- RenameCells(nemo, new.names = paste0(nemo$pool, '_', nemo$CellRanger_barcode))
head(colnames(nemo))

rownames(total_fragments) <- total_fragments$CB
nemo$peak_region_fragments <- total_fragments[colnames(nemo), "frequency_count"]
nemo$peak_reads_count <- total_fragments[colnames(nemo), "reads_count"]
nemo$nucleosome_free_peak_count <- total_fragments[colnames(nemo), "nucleosome_free"]
nemo$mononucleosomal_peak_count <- total_fragments[colnames(nemo), "mononucleosomal"]
nemo$log10nfragments <- log10(nemo$peak_region_fragments) 
nemo$log10PeaksPerUMI <- log10(nemo$nFeature_peaks)/log10(nemo$nCount_peaks)
nemo$peaks_per_umi <- nemo$nFeature_peaks/nemo$nCount_peaks

nemo <- FRiP(
  object = nemo,
  assay = 'peaks',
  total.fragments = 'peak_region_fragments',
  col.name = 'pct_reads_in_peaks'
)

nemo <- RenameCells(nemo, new.names = paste0(nemo$individual, '_', nemo$CellRanger_barcode))
head(colnames(nemo))

saveRDS(nemo, file.path(main_dir, paste0('nemo_gex.atac.peaks.mergedObject.rds')))


#############

nemo_poa <- readRDS("~/scratch/nemo_multi/seurat/transfer_mapping/nemo_poa_only_mapped_to_hypo_celltypes_02.20.26.rds")

annotations <- readGFFAsGRanges('~/scratch/nemo_multi/nemo_ref_genome/GCF_022539595.1_ASM2253959v1_genomic.gtf')
ann_df <- as.data.frame(annotations)

# Load genome fasta
nemo_genome_file <- '~/scratch/nemo_multi/nemo_ref_genome/GCF_022539595.1_ASM2253959v1_genomic.fna'
nemo_genome <- read_genome(format = 'fasta', file = nemo_genome_file)
names(nemo_genome) = reshape2::colsplit(names(nemo_genome), " ", c("1", '2'))[,1]


# nemo_genome is your DNAStringSet
seqnames    <- names(nemo_genome)
seqlengths  <- width(nemo_genome)
# if you know any scaffolds are circular, set TRUE; otherwise FALSE (or leave as NA)
is_circular <- rep(FALSE, length(seqnames))

# create the Seqinfo
nemo_seqinfo <- Seqinfo(
  seqnames    = seqnames,
  seqlengths  = seqlengths,
  isCircular  = is_circular,
  genome      = "Nemo"    # or whatever genome build you want to label
)

# colnames(annotations@elementMetadata)[colnames(annotations@elementMetadata) == 'gene'] = 'gene_name'

# ensure the "transcript" feature types all have a value for "gene_biotype"
ann_df <- ann_df %>%
  group_by(gene_id) %>%
  mutate(
    gene_biotype = if_else(
      type == "transcript",
      # explicitly call dplyr::first()
      dplyr::first(gene_biotype[type == "gene"], 
                   default = gene_biotype[1]),
      gene_biotype
    )
  ) %>%
  ungroup()

# remove all columns that are not required for cellranger
colnames(ann_df)
columns_keep <- c("seqnames", "start", "end", "width", "strand", "NCBI", "type", "phase", "gene_biotype", "gbkey", "gene_id", "transcript_id", "gene")
ann_df <- ann_df %>% dplyr::select(any_of(columns_keep))

ann_df$row_order <- rownames(ann_df)
nemo_df_gene <- ann_df %>% filter(type == c('gene'))
nemo_df_trans <- ann_df %>% filter(type == c('transcript'))
nemo_df_exon <- ann_df %>% filter(type == c('exon'))
nemo_df_start <- ann_df %>% filter(type == c('start_codon'))
nemo_df_stop <- ann_df %>% filter(type == c('stop_codon'))
ann_df <- rbind(nemo_df_trans, 
                nemo_df_exon, 
                nemo_df_start,
                nemo_df_stop,
                nemo_df_gene)
ann_df <- ann_df %>%
  mutate(row_order = as.numeric(row_order)) %>%
  arrange(row_order)
ann_df$gene_name <- ann_df$gene_id
ann_df$row_order <- NULL
# Convert to GRanges and confirm with df
gtf_gr <- makeGRangesFromDataFrame(
  ann_df,
  keep.extra.columns = TRUE,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)
gtf_gr_df <- as.data.frame(gtf_gr)

counts <- GetAssayData(nemo, assay = 'peaks', layer = 'counts')
cells_target <- colnames(nemo_poa)
cells_common <- intersect(cells_target, colnames(counts))

peak_counts2 <- counts[, cells_common, drop = FALSE]

stopifnot(identical(colnames(nemo_poa), colnames(peak_counts2)))

merged_assay <- CreateChromatinAssay(
  counts = peak_counts2,
  validate.fragments = TRUE,
  genome = nemo_seqinfo,
  annotation = gtf_gr,
  fragments = nemo[["peaks"]]@fragments,
  ranges = nemo[['peaks']]@ranges
)
nemo2 <- nemo_poa
# Attach back to your Seurat object
nemo_poa[["peaks"]] <- merged_assay

# Optional: verify annotation alignment
Annotation(nemo_poa[["peaks"]])

head(rownames(nemo_poa@assays$peaks@data))
head(rownames(nemo_poa@assays$peaks@counts))
head(rownames(nemo_poa@assays$peaks@meta.features))

rownames(nemo_poa@assays$peaks@data)   = stringr::str_replace(rownames(nemo_poa@assays$peaks@data),   "-", "_")
rownames(nemo_poa@assays$peaks@counts) = stringr::str_replace(rownames(nemo_poa@assays$peaks@counts), "-", "_")
rownames(nemo_poa@assays$peaks@meta.features) = stringr::str_replace(rownames(nemo_poa@assays$peaks@meta.features), "-", "_")

head(rownames(nemo_poa@assays$peaks@data))
head(rownames(nemo_poa@assays$peaks@counts))
head(rownames(nemo_poa@assays$peaks@meta.features))

DefaultAssay(nemo_poa) <- 'RNA'
CoveragePlot(
  object = nemo_poa,
  region = "cckb",
  features = "cckb",
  assay = "peaks",
  extend.upstream = 5000,
  extend.downstream = 5000,
  # expression.assay = "RNA",
  annotation = TRUE,
  ymax = "q95",
  group.by = "final_clusters",
  window = 100
)

CoveragePlot(
  object = nemo_poa,
  region = "cckb",
  features = "cckb",
  assay = "peaks",
  extend.upstream = 5000,
  extend.downstream = 5000,
  # expression.assay = "RNA",
  annotation = TRUE,
  ymax = "q95",
  group.by = "final_clusters",
  idents = 6,
  split.by = 'Status',
  window = 200,
  downsample.rate = 0.05,
  heights = c(3,1)
)

CoveragePlot(
  object = nemo_poa,
  region = "lhb",
  features = "lhb",
  extend.upstream = 5000,
  extend.downstream = 5000,
  expression.assay = "RNA",
  group.by = "gabe_celltype_annotations",
  assay = "peaks",
  sep = c("-", "-"),
  peaks = TRUE,
  annotation = TRUE,
  ymax = "q95",
  links = FALSE,
  ranges = nemo[['peaks']]@ranges
  window = 100
)

DefaultAssay(nemo_poa) <- 'RNA'
nemo_poa[['GeneActivity']] <- NULL
gene.activities <- GeneActivity(nemo_poa, 
                                assay = "peaks",
                                extend.upstream = 2000,
                                extend.downstream = 50,
                                max.width = NULL,
                                biotypes = NULL)

# add gene activities as a new assay
nemo_poa[["GeneActivity"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(nemo_poa) <- "GeneActivity"
nemo_poa <- NormalizeData(object = nemo_poa, normalization.method = "LogNormalize", assay = 'GeneActivity', scale.factor = median(nemo_poa$nCount_GeneActivity))

gc()

# Set your default assay to GeneActivity
DefaultAssay(nemo_poa) <- "GeneActivity"

# Compute average accessibility per cluster (replace with your cluster column)
avg_access <- AverageExpression(
  nemo_poa,
  assays = "GeneActivity",
  group.by = "final_clusters",   # or "seurat_clusters"
  slot = "data"                  # "data" = normalized, "counts" = raw
)

# Inspect results
head(avg_access$GeneActivity[, 1:5])

avg_df <- as.data.frame(avg_access$GeneActivity)
colnames(avg_df) <- sub("^g", "", colnames(avg_df))

avg_df$gene <- rownames(avg_df)
saveRDS(avg_df, '/storage/scratch1/8/kleatherbury3/nemo_multi/seurat/seurat_objects/predictedRNA_expression_from_chromatin_assesibility_per_final_cluster_df_03.04.26.rds')

# Compute average accessibility per cluster (replace with your cluster column)
avg_access <- AverageExpression(
  nemo_poa,
  assays = "GeneActivity",
  group.by = "subcluster_L3_numeric",   # or "seurat_clusters"
  slot = "data"                  # "data" = normalized, "counts" = raw
)

# Inspect results
head(avg_access$GeneActivity[, 1:5])

avg_df <- as.data.frame(avg_access$GeneActivity)
colnames(avg_df) <- sub("^g", "", colnames(avg_df))
avg_df$Gene <- rownames(avg_df)

saveRDS(avg_df, '/storage/scratch1/8/kleatherbury3/nemo_multi/seurat/seurat_objects/predictedRNA_expression_from_chromatin_assesibility_per_subcluster_L3_numeric_df_03.04.26.rds')

colnames(nemo_poa@meta.data)[c(71,72,73,76,77,85,86, 78, 90:95, 106:111)]
nemo_poa@meta.data <- nemo_poa@meta.data[, -c(71,72,73,76,77,85,86, 78, 90:95, 106:111)]
colnames(nemo_poa@meta.data)

colnames(nemo_poa@meta.data)[92:107] <- paste0("MbunaHypo_", colnames(nemo_poa@meta.data)[92:107])
colnames(nemo_poa@meta.data)[108:119] <- paste0("MbunaTelOB_", colnames(nemo_poa@meta.data)[108:119])
colnames(nemo_poa@meta.data)[120:121] <- paste0("MbunaHypo_", colnames(nemo_poa@meta.data)[120:121])
colnames(nemo_poa@meta.data)

saveRDS(nemo_poa, '/storage/scratch1/8/kleatherbury3/nemo_multi/seurat/seurat_objects/nemo_ATAC_corrected_SCT.RNA.GA_final_object_03.04.26.rds')






############# PLOTTING AND VISUALIZATION #############

## -------------------------------------------
## density-based two-tailed thresholds
## -------------------------------------------
bounds_by_density <- function(x, frac_of_peak = 0.05, log10_scale = TRUE,
                              bw = "nrd0", adjust = 1, na.rm = TRUE) {
  x <- x[is.finite(x)]
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) < 2L) return(c(lb = NA_real_, ub = NA_real_))
  
  if (log10_scale) {
    x <- x[x > 0]
    if (length(x) < 2L) return(c(lb = NA_real_, ub = NA_real_))
    z <- log10(x)
  } else {
    z <- x
  }
  
  d <- stats::density(z, bw = bw, adjust = adjust, na.rm = TRUE)
  if (length(d$y) < 3L || !is.finite(max(d$y))) return(c(lb = NA_real_, ub = NA_real_))
  
  level <- frac_of_peak * max(d$y)
  if (!is.finite(level) || level <= min(d$y) || level >= max(d$y))
    return(c(lb = NA_real_, ub = NA_real_))
  
  ip <- which.max(d$y)
  
  ## left crossing: find k with y[k] <= level < y[k+1]
  iL <- NA_integer_
  for (k in seq(ip - 1L, 1L, by = -1L)) {
    if (isTRUE(d$y[k] <= level && d$y[k + 1L] > level)) { iL <- k; break }
  }
  xL <- if (is.na(iL)) NA_real_ else
    approx(x = d$y[c(iL, iL + 1L)], y = d$x[c(iL, iL + 1L)],
           xout = level, ties = "ordered")$y
  
  ## right crossing: find k with y[k] >= level > y[k+1]
  iR <- NA_integer_
  for (k in seq(ip, length(d$y) - 1L, by = 1L)) {
    if (isTRUE(d$y[k] >= level && d$y[k + 1L] < level)) { iR <- k; break }
  }
  xR <- if (is.na(iR)) NA_real_ else
    approx(x = d$y[c(iR, iR + 1L)], y = d$x[c(iR, iR + 1L)],
           xout = level, ties = "ordered")$y
  
  res <- if (log10_scale) c(lb = 10^xL, ub = 10^xR) else c(lb = xL, ub = xR)
  attr(res, "level") <- level
  res
}



## -------------------------
## plotting
## -------------------------

sample <- readRDS('/storage/home/hcoda1/8/kleatherbury3/scratch/nemo_multi/seurat/peakCalling_macs2_v2/nemo_gex.atac.peaks.mergedObject.rds')
qc_dir <- file.path('/storage/home/hcoda1/8/kleatherbury3/scratch/nemo_multi/seurat', 'QC')
if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

seq_poolQC_dir <- file.path(qc_dir, 'seq_pool')
if (!dir.exists(seq_poolQC_dir)) dir.create(seq_poolQC_dir, recursive = TRUE)

ind_col <- if ("individual" %in% names(sample@meta.data)) "individual" else "Individual"
table(sample@meta.data[[ind_col]])

# Reorder individuals
new_levels <- mixedsort(unique(sample@meta.data[[ind_col]]))
sample@meta.data[[ind_col]] <- factor(
  sample@meta.data[[ind_col]], 
  levels = new_levels)

indv_colors <- c(
  'lightskyblue',
  'palegreen',
  'yellow1',
  'hotpink',
  'sandybrown',
  'plum3',
  'turquoise',
  "tomato"
)

parameters <- c('nCount_RNA', 
                'nFeature_RNA',
                'nCount_peaks', 
                'nFeature_peaks', 
                'log10GenesPerUMI', 
                'log10PeaksPerUMI',
                'pct.gex_mapped',
                'pct.atac_mapped')

filtered_violin_plots1 <- list()

for (parameter in parameters) {
  e <- VlnPlot_scCustom(sample, 
                        features = parameter, 
                        group.by = "individual", 
                        num_columns = 1,
                        log = F, 
                        pt.size = 0.1, 
                        alpha = 0.1, 
                        plot_median = T, 
                        colors_use = indv_colors,
                        layer = 'counts') + 
    xlab(label = '') +
    labs(title = paste(parameter)) +
    theme(
      plot.title = element_text(size = 12, 
                                color = 'black', 
                                hjust = 0.5, 
                                face = "bold",
                                vjust = -1)) +
    NoLegend()
  
  filtered_violin_plots1[[parameter]] <- e
}

c <- wrap_plots(filtered_violin_plots1, ncol = 2) +
  plot_annotation(
    title = paste0(toupper(seq_pool), " Tel + OB Merged RNA:  Post-QC Metrics"),
    theme = theme(
      plot.title = element_text(size = 15, 
                                color = 'blue3', 
                                hjust = 0.5,
                                vjust = -0.8,
                                face = "bold")))
print(c)


parameters <- c('pct.mt', 
                'pct.cd', 
                'pct.mito_atac', 
                'pct_reads_in_peaks',
                'peak_region_fragments',
                'nucleosome_signal',
                'TSS.enrichment',
                'peak_reads_count'
)

filtered_violin_plots2 <- list()

for (parameter in parameters) {
  e <- VlnPlot_scCustom(sample, 
                        features = parameter, 
                        group.by = "individual", 
                        num_columns = 1,
                        log = F, 
                        pt.size = 0.1, 
                        alpha = 0.1, 
                        plot_median = T, 
                        colors_use = indv_colors,
                        layer = 'counts') + 
    xlab(label = '') +
    labs(title = paste(parameter)) +
    theme(
      plot.title = element_text(size = 12, 
                                color = 'black', 
                                hjust = 0.5, 
                                face = "bold",
                                vjust = -1)) +
    NoLegend()
  filtered_violin_plots2[[parameter]] <- e
}

d <- wrap_plots(filtered_violin_plots2,ncol = 4) +
  plot_annotation(
    title = paste0(toupper(sample), " Tel + OB Merged RNA:  Post-QC Metrics"),
    theme = theme(
      plot.title = element_text(size = 15, 
                                color = 'blue3', 
                                hjust = 0.5,
                                vjust = -0.8,
                                face = "bold")))
print(d)

parameters <- c('nCount_RNA', 
                'nFeature_RNA',
                'nCount_peaks', 
                'nFeature_peaks', 
                'log10GenesPerUMI', 
                'log10PeaksPerUMI',
                'pct.gex_mapped',
                'pct.atac_mapped')

frac_of_peak <- 0.05 # if bounds look too tight/strict then increase this value

# which params use log10 on x (others stay linear)
log_params <- c("nCount_RNA","nFeature_RNA","nCount_peaks","nFeature_peaks")

filtered_density_plots1 <- list()

for (parameter in parameters) {
  x <- sample@meta.data[[parameter]]
  median_val <- median(x, na.rm = TRUE)
  
  use_log10 <- parameter %in% log_params
  
  # density-based bounds (same as before)
  b  <- bounds_by_density(x, frac_of_peak = frac_of_peak, log10_scale = use_log10)
  lb <- unname(b["lb"]); ub <- unname(b["ub"])
  
  # label function for vline annotations (keeps <1 as decimals; no scientific)
  lbl <- function(v) trimws(formatC(v, format = "f", digits = 3, drop0trailing = TRUE, big.mark = ","))
  
  p <- sample@meta.data %>%
    ggplot(aes(x = .data[[parameter]],
               color = .data[[ind_col]],
               fill  = .data[[ind_col]])) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    scale_color_manual(values = indv_colors, name = "Individual Key") +
    scale_fill_manual(values  = indv_colors, name = "Individual Key") +
    {
      if (use_log10)
        scale_x_log10(
          breaks = scales::log_breaks(),
          labels = function(x) lbl(x)   # <- no scientific, preserves decimals
        )
      else
        scale_x_continuous(
          labels = function(x) lbl(x)   # <- no scientific, preserves decimals
        )
    } +
    theme_classic() +
    theme(
      plot.margin = margin(15, 15, 15, 15),
      plot.title  = element_text(size = 12, 
                                 color = 'black',
                                 hjust = 0.5, 
                                 face = "bold", 
                                 vjust = -1)
    ) +
    ylab("Cell density") +
    
    # median line + label
    geom_vline(xintercept = median_val, 
               color = 'black', 
               linewidth = 0.6, 
               linetype = 'solid') +
    annotate('label', 
             x = median_val, 
             y = 0,
             label = paste("md=", lbl(median_val)),
             size = 3, 
             color = 'black', 
             vjust = -7, 
             hjust = 0.5) +
    
    # lower / upper density-threshold lines + labels
    { if (is.finite(lb)) geom_vline(xintercept = lb, 
                                    color = "firebrick", 
                                    linewidth = 0.7, 
                                    linetype = 'dashed') else NULL } +
    { if (is.finite(lb)) annotate('label', 
                                  x = lb, 
                                  y = 0, 
                                  label = paste('LQ\n', lbl(lb)),
                                  size = 3, 
                                  color = "firebrick",
                                  vjust = -4.2, 
                                  hjust = 0.5) else NULL } +
    { if (is.finite(ub)) geom_vline(xintercept = ub, 
                                    color = "firebrick", 
                                    linewidth = 0.9) else NULL } +
    { if (is.finite(ub)) annotate('label', 
                                  x = ub, 
                                  y = 0, 
                                  label = paste('LQ\n', lbl(ub)),
                                  size = 3, 
                                  color = "firebrick",
                                  vjust = -4.2, 
                                  hjust = 0.5) else NULL } +
    labs(title = paste("\n Density of", parameter, '\n'))
  
  filtered_density_plots1[[parameter]] <- p
}

a <- wrap_plots(filtered_density_plots1, ncol = 2, guides = "collect") +
  plot_annotation(
    title = paste0(
      toupper(seq_pool),
      " Tel + OB Merged RNA:  Density Plots of Filtered QC Metrics"
    ),
    theme = theme(plot.title = element_text(
      size = 15, color = 'blue3', hjust = 0.5, vjust = -1, face = "bold"
    ))
  ) &
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.key.width  = unit(18, "pt"),
    legend.spacing.x  = unit(10, "pt"),
    legend.title      = element_text(face = "bold")
  ) &
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE, title.position = "top", title.hjust = 0.5),
    fill  = guide_legend(nrow = 1, byrow = TRUE, title.position = "top", title.hjust = 0.5)
  )
print(a)


parameters <- c('pct.mt', 
                'pct.cd', 
                'pct.mito_atac', 
                'pct_reads_in_peaks', 
                'peak_region_fragments', 
                'nucleosome_signal', 
                'TSS.enrichment', 
                'peak_reads_count')

filtered_density_plots2 <- list()

for (parameter in parameters) {
  x <- sample@meta.data[[parameter]]
  median_val <- median(x, na.rm = TRUE)
  
  use_log10 <- parameter %in% log_params
  
  # density-based bounds (same as before)
  b  <- bounds_by_density(x, frac_of_peak = frac_of_peak, log10_scale = use_log10)
  lb <- unname(b["lb"]); ub <- unname(b["ub"])
  
  # label function for vline annotations (keeps <1 as decimals; no scientific)
  lbl <- function(v) trimws(formatC(v, format = "f", digits = 3, drop0trailing = TRUE, big.mark = ","))
  
  p <- sample@meta.data %>%
    ggplot(aes(x = .data[[parameter]],
               color = .data[[ind_col]],
               fill  = .data[[ind_col]])) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    scale_color_manual(values = indv_colors, name = "Individual Key") +
    scale_fill_manual(values  = indv_colors, name = "Individual Key") +
    {
      if (use_log10)
        scale_x_log10(
          breaks = scales::log_breaks(),
          labels = function(x) lbl(x)   # <- no scientific, preserves decimals
        )
      else
        scale_x_continuous(
          labels = function(x) lbl(x)   # <- no scientific, preserves decimals
        )
    } +
    theme_classic() +
    theme(
      plot.margin = margin(15, 15, 15, 15),
      plot.title  = element_text(size = 12, 
                                 color = 'black',
                                 hjust = 0.5, 
                                 face = "bold", 
                                 vjust = -1)
    ) +
    ylab("Cell density") +
    
    # median line + label
    geom_vline(xintercept = median_val, 
               color = 'black', 
               linewidth = 0.6, 
               linetype = 'solid') +
    annotate('label', 
             x = median_val, 
             y = 0,
             label = paste("md=", lbl(median_val)),
             size = 3, 
             color = 'black', 
             vjust = -7, 
             hjust = 0.5) +
    
    # lower / upper density-threshold lines + labels
    { if (is.finite(lb)) geom_vline(xintercept = lb, 
                                    color = "firebrick", 
                                    linewidth = 0.7, 
                                    linetype = 'dashed') else NULL } +
    { if (is.finite(lb)) annotate('label', 
                                  x = lb, 
                                  y = 0, 
                                  label = paste('LQ\n', lbl(lb)),
                                  size = 3, 
                                  color = "firebrick",
                                  vjust = -4.2, 
                                  hjust = 0.5) else NULL } +
    { if (is.finite(ub)) geom_vline(xintercept = ub, 
                                    color = "firebrick", 
                                    linewidth = 0.9) else NULL } +
    { if (is.finite(ub)) annotate('label', 
                                  x = ub, 
                                  y = 0, 
                                  label = paste('LQ\n', lbl(ub)),
                                  size = 3, 
                                  color = "firebrick",
                                  vjust = -4.2, 
                                  hjust = 0.5) else NULL } +
    labs(title = paste("\n Density of", parameter, '\n'))
  
  filtered_density_plots2[[parameter]] <- p
}
b <- wrap_plots(filtered_density_plots2, ncol = 2, guides = "collect") +
  plot_annotation(
    title = paste0(
      toupper(seq_pool),
      " Tel + OB Merged RNA:  Density Plots of Filtered QC Metrics"
    ),
    theme = theme(plot.title = element_text(
      size = 15, color = 'blue3', hjust = 0.5, vjust = -1, face = "bold"
    ))
  ) &
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.key.width  = unit(18, "pt"),
    legend.spacing.x  = unit(10, "pt"),
    legend.title      = element_text(face = "bold")
  ) &
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE, title.position = "top", title.hjust = 0.5),
    fill  = guide_legend(nrow = 1, byrow = TRUE, title.position = "top", title.hjust = 0.5)
  )
print(b)

a <- CoveragePlot(
  object = nemo,
  region = "gfap",        
  features = "gfap",
  extend.upstream = 5000,
  extend.downstream = 5000,
  expression.assay = 'RNA',
  assay = 'ATAC',
  sep = c('-', '-'),
  peaks = T,
  annotation = T,
  idents = 15,
  ymax = 'q95',
  ranges = final_object[['peaks']]@ranges,
  links = F
)

b <- CoveragePlot(
  object = final_object,
  region = "gfap",        
  features = "gfap",
  extend.upstream = 5000,
  extend.downstream = 5000,
  expression.assay = 'RNA',
  assay = 'peaks',
  sep = c('-', '-'),
  peaks = T,
  annotation = T,
  ymax = 'q95',
  # ranges = nemo[['ATAC']]@ranges,
  links = F
)

p1 <- wrap_plots(a,c, ncol =1)
p1

c <- CoveragePlot(
  object = nemo,
  region = "cga",        
  features = "cga",
  extend.upstream = 5000,
  extend.downstream = 5000,
  expression.assay = 'RNA',
  assay = 'ATAC',
  sep = c('-', '-'),
  peaks = T,
  annotation = T,
  idents = 22,
  ymax = 'q95',
  ranges = final_object[['peaks']]@ranges,
  links = F
)

d <- CoveragePlot(
  object = final_object,
  region = "cga",        
  features = "cga",
  extend.upstream = 5000,
  extend.downstream = 5000,
  expression.assay = 'RNA',
  assay = 'ATAC',
  sep = c('-', '-'),
  peaks = T,
  annotation = T,
  ymax = 'q95',
  ranges = all_ranges,
  links = F
)

p2 <- wrap_plots(c,d, ncol =1)
p2


