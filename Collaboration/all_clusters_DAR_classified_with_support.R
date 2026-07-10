# ============================================================
# All-cluster pairwise DAR analysis with classification
# Assumes Seurat object `obj` is already loaded in memory
#
# Adds peak support metrics to master tables:
#   total_counts_all_individuals
#   nonzero_individuals_all
#   mean_counts_all_individuals
#   max_counts_any_individual
#   total_counts_test_individuals
#   nonzero_individuals_test
#   mean_counts_test_individuals
#   max_counts_test_individual
#
# Designed for RStudio/OOD / Campus Cluster:
#   - does NOT readRDS()
#   - does NOT create obj$final_clusters unless needed
#   - writes outputs to file.path(getwd(), "all_clusters_DAR_classified_with_support")
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(Signac)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(GenomeInfoDb)
  library(ggplot2)
})

############################################################
# USER SETTINGS
############################################################

# Object should already be loaded as `obj`
if (!exists("obj")) {
  stop("Object `obj` does not exist. Load your Seurat object first, then run this script.")
}

if (!inherits(obj, "Seurat")) {
  stop("`obj` exists but is not a Seurat object.")
}

# Cluster column to use
cluster_column <- "res0.8_50nn_40PC_45LSI"

# Preferred assay
atac_assay <- "ATAC"

# Groups to compare for DAR testing
groups_to_test <- c("M", "F", "D")

# Filtering thresholds
# Change these here if desired.
min_individuals_with_peak <- 10
min_total_counts <- 80

# Significance threshold
fdr_cutoff <- 0.20

# Output directory
outdir <- file.path(getwd(), "all_clusters_DAR_classified_with_support")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile_peak_csv     <- file.path(outdir, "all_clusters_DARs_peak_level_classified_with_support.csv")
outfile_pairwise_csv <- file.path(outdir, "all_clusters_DARs_pairwise_with_classification_and_support.csv")
outfile_summary_csv  <- file.path(outdir, "all_clusters_DAR_summary_by_cluster.csv")
outfile_plot_png     <- file.path(outdir, "all_clusters_DAR_classification_stacked_barplot.png")

cat("\nSaving outputs to:\n", outdir, "\n\n", sep = "")

# Check required metadata
required_meta <- c(cluster_column, "individual", "Status")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))

if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

if (!(atac_assay %in% names(obj@assays))) {
  stop("ATAC assay not found. Available assays: ", paste(names(obj@assays), collapse = ", "))
}

# Clusters to analyze
cluster_ids <- sort(unique(as.numeric(as.character(obj@meta.data[[cluster_column]]))))
cluster_ids <- as.character(cluster_ids[!is.na(cluster_ids)])

cat("Clusters to analyze:\n")
print(cluster_ids)
cat("\n")

############################################################
# HELPER FUNCTIONS
############################################################

get_counts_matrix <- function(seurat_obj, assay_name) {
  out <- tryCatch(
    GetAssayData(seurat_obj, assay = assay_name, layer = "counts"),
    error = function(e) {
      GetAssayData(seurat_obj, assay = assay_name, slot = "counts")
    }
  )

  if (!inherits(out, "Matrix")) {
    out <- as(out, "dgCMatrix")
  }

  out
}

subset_by_cluster <- function(seurat_obj, cluster_column, cluster_id) {
  cluster_keep <- rownames(seurat_obj@meta.data)[
    as.character(seurat_obj@meta.data[[cluster_column]]) == as.character(cluster_id)
  ]

  if (length(cluster_keep) == 0) {
    return(NULL)
  }

  subset(seurat_obj, cells = cluster_keep)
}

parse_peak <- function(x) {
  x <- as.character(x)

  chromosome <- NA_character_
  start <- NA_real_
  end <- NA_real_

  if (grepl("^.+:[0-9]+-[0-9]+$", x)) {
    chromosome <- sub("^(.+):[0-9]+-[0-9]+$", "\\1", x)
    start <- as.numeric(sub("^.+:([0-9]+)-([0-9]+)$", "\\1", x))
    end   <- as.numeric(sub("^.+:([0-9]+)-([0-9]+)$", "\\2", x))
  } else if (grepl("^.+-[0-9]+-[0-9]+$", x)) {
    chromosome <- sub("^(.*)-[0-9]+-[0-9]+$", "\\1", x)
    start <- as.numeric(sub("^.*-([0-9]+)-([0-9]+)$", "\\1", x))
    end   <- as.numeric(sub("^.*-([0-9]+)-([0-9]+)$", "\\2", x))
  } else if (grepl("^.+_[0-9]+_[0-9]+$", x)) {
    chromosome <- sub("^(.*)_[0-9]+_[0-9]+$", "\\1", x)
    start <- as.numeric(sub("^.*_([0-9]+)_([0-9]+)$", "\\1", x))
    end   <- as.numeric(sub("^.*_([0-9]+)_([0-9]+)$", "\\2", x))
  }

  data.frame(
    peak = x,
    chromosome = chromosome,
    start = start,
    end = end,
    stringsAsFactors = FALSE
  )
}

normalize_seqname <- function(x) {
  x <- as.character(x)
  x <- sub("^([A-Za-z]{2,})-([0-9].*)$", "\\1_\\2", x)
  x
}

safe_num_fmt <- function(x, digits = 15) {
  ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))
}

empty_pairwise_result_table <- function() {
  data.frame(
    cluster_id = character(),
    peak = character(),
    chromosome = character(),
    chromosome_norm = character(),
    start = numeric(),
    end = numeric(),
    nearest_gene = character(),
    distance_to_gene = numeric(),
    nearest_feature_type = character(),
    comparison = character(),
    logFC = numeric(),
    logCPM = numeric(),
    PValue = numeric(),
    FDR = numeric(),
    classification = character(),
    mean_logCPM_M = numeric(),
    mean_logCPM_F = numeric(),
    mean_logCPM_D = numeric(),
    total_counts_all_individuals = numeric(),
    nonzero_individuals_all = integer(),
    mean_counts_all_individuals = numeric(),
    max_counts_any_individual = numeric(),
    total_counts_test_individuals = numeric(),
    nonzero_individuals_test = integer(),
    mean_counts_test_individuals = numeric(),
    max_counts_test_individual = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_peak_result_table <- function() {
  data.frame(
    cluster_id = character(),
    peak = character(),
    chromosome = character(),
    chromosome_norm = character(),
    start = numeric(),
    end = numeric(),
    nearest_gene = character(),
    distance_to_gene = numeric(),
    nearest_feature_type = character(),
    mean_logCPM_M = numeric(),
    mean_logCPM_F = numeric(),
    mean_logCPM_D = numeric(),
    M_vs_F_logFC = numeric(),
    M_vs_F_FDR = numeric(),
    M_vs_D_logFC = numeric(),
    M_vs_D_FDR = numeric(),
    F_vs_D_logFC = numeric(),
    F_vs_D_FDR = numeric(),
    sig_M_vs_F = logical(),
    sig_M_vs_D = logical(),
    sig_F_vs_D = logical(),
    D_position_between_M_and_F = numeric(),
    classification = character(),
    total_counts_all_individuals = numeric(),
    nonzero_individuals_all = integer(),
    mean_counts_all_individuals = numeric(),
    max_counts_any_individual = numeric(),
    total_counts_test_individuals = numeric(),
    nonzero_individuals_test = integer(),
    mean_counts_test_individuals = numeric(),
    max_counts_test_individual = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_cluster_summary_table <- function() {
  data.frame(
    cluster_id = character(),
    cells_in_cluster = integer(),
    individuals_total = integer(),
    individuals_M = integer(),
    individuals_F = integer(),
    individuals_D = integer(),
    peaks_before_filter = integer(),
    peaks_after_filter = integer(),
    n_M_vs_F_DARs = integer(),
    n_M_vs_D_DARs = integer(),
    n_F_vs_D_DARs = integer(),
    n_unique_DARs = integer(),
    status = character(),
    stringsAsFactors = FALSE
  )
}

classify_one_peak <- function(mean_M, mean_F, mean_D,
                              sig_MF, sig_MD, sig_FD,
                              progressive_lower = 0.25,
                              progressive_upper = 0.75) {

  if (is.na(mean_M) || is.na(mean_F) || is.na(mean_D)) {
    return("Unclassified")
  }

  if (mean_M == mean_F) {
    if (mean_D > mean_M && sig_MD && sig_FD) {
      return("Transiently Upregulated")
    }
    if (mean_D < mean_M && sig_MD && sig_FD) {
      return("Transiently Downregulated")
    }
    return("Unclassified")
  }

  direction <- ifelse(mean_F > mean_M, "Upregulated", "Downregulated")

  # Position of D along the M -> F interval:
  #   0 = exactly at M
  #   1 = exactly at F
  p <- (mean_D - mean_M) / (mean_F - mean_M)

  # D lies outside the M-F interval
  if (p < 0 || p > 1) {

    # Only call transient if D differs significantly from BOTH endpoints
    if (sig_MD && sig_FD) {
      if (p > 1) {
        return(paste("Transiently", direction))
      }
      if (p < 0) {
        opposite_direction <- ifelse(direction == "Upregulated",
                                     "Downregulated",
                                     "Upregulated")
        return(paste("Transiently", opposite_direction))
      }
    }

    # Otherwise assign to nearest endpoint
    if (abs(mean_D - mean_F) <= abs(mean_D - mean_M)) {
      return(paste("Early", direction))
    } else {
      return(paste("Late", direction))
    }
  }

  # D lies between M and F
  if (p >= progressive_lower && p <= progressive_upper) {
    return(paste("Progressively", direction))
  }

  if (p > progressive_upper) {
    return(paste("Early", direction))
  }

  if (p < progressive_lower) {
    return(paste("Late", direction))
  }

  return("Unclassified")
}

annotate_peak_table <- function(cluster_obj, peak_level_df, atac_assay = "ATAC") {

  if (nrow(peak_level_df) == 0) return(peak_level_df)

  peak_info <- dplyr::bind_rows(lapply(unique(peak_level_df$peak), parse_peak)) %>%
    dplyr::distinct(peak, .keep_all = TRUE)

  peak_level_df <- peak_level_df %>%
    dplyr::left_join(peak_info, by = "peak")

  peak_level_df$chromosome_norm <- normalize_seqname(peak_level_df$chromosome)

  if (!(atac_assay %in% names(cluster_obj@assays)) ||
      !inherits(cluster_obj[[atac_assay]], "ChromatinAssay")) {
    warning(atac_assay, " assay is not a ChromatinAssay. Returning table without gene annotation.")
    peak_level_df$nearest_gene <- NA_character_
    peak_level_df$distance_to_gene <- NA_real_
    peak_level_df$nearest_feature_type <- NA_character_
    return(peak_level_df)
  }

  annot_gr <- Annotation(cluster_obj[[atac_assay]])

  if (is.null(annot_gr) || length(annot_gr) == 0) {
    warning("No annotation attached to cluster_obj[[", atac_assay, "]]. Returning table without gene annotation.")
    peak_level_df$nearest_gene <- NA_character_
    peak_level_df$distance_to_gene <- NA_real_
    peak_level_df$nearest_feature_type <- NA_character_
    return(peak_level_df)
  }

  old_seqlevels <- seqlevels(annot_gr)
  new_seqlevels <- normalize_seqname(old_seqlevels)
  names(new_seqlevels) <- old_seqlevels
  annot_gr <- renameSeqlevels(annot_gr, new_seqlevels)

  valid_idx <- !is.na(peak_level_df$chromosome_norm) &
    !is.na(peak_level_df$start) &
    !is.na(peak_level_df$end)

  peak_valid <- peak_level_df[valid_idx, , drop = FALSE]

  if (nrow(peak_valid) == 0) {
    peak_level_df$nearest_gene <- NA_character_
    peak_level_df$distance_to_gene <- NA_real_
    peak_level_df$nearest_feature_type <- NA_character_
    return(peak_level_df)
  }

  dar_gr <- GRanges(
    seqnames = peak_valid$chromosome_norm,
    ranges = IRanges(start = peak_valid$start, end = peak_valid$end)
  )

  common_seqlevels <- intersect(
    unique(as.character(seqnames(dar_gr))),
    unique(as.character(seqnames(annot_gr)))
  )

  if (length(common_seqlevels) == 0) {
    peak_level_df$nearest_gene <- NA_character_
    peak_level_df$distance_to_gene <- NA_real_
    peak_level_df$nearest_feature_type <- NA_character_
    return(peak_level_df)
  }

  keep_dar <- as.character(seqnames(dar_gr)) %in% common_seqlevels
  dar_gr <- dar_gr[keep_dar]
  peak_valid <- peak_valid[keep_dar, , drop = FALSE]

  if (nrow(peak_valid) == 0) {
    peak_level_df$nearest_gene <- NA_character_
    peak_level_df$distance_to_gene <- NA_real_
    peak_level_df$nearest_feature_type <- NA_character_
    return(peak_level_df)
  }

  annot_gr <- keepSeqlevels(annot_gr, common_seqlevels, pruning.mode = "coarse")

  annot_cols <- colnames(mcols(annot_gr))

  gene_col_candidates <- c(
    "gene_name", "gene", "symbol", "gene_symbol",
    "Name", "GENE", "SYMBOL", "transcript_name", "tx_name"
  )
  gene_col <- gene_col_candidates[gene_col_candidates %in% annot_cols][1]

  feature_col_candidates <- c(
    "type", "gene_biotype", "transcript_biotype",
    "feature_type", "annotation", "tx_type"
  )
  feature_col <- feature_col_candidates[feature_col_candidates %in% annot_cols][1]

  if (!is.na(gene_col) && length(gene_col) > 0) {
    annot_df <- as.data.frame(annot_gr)
    annot_df$gene_use <- as.character(annot_df[[gene_col]])
    annot_df <- annot_df[!is.na(annot_df$gene_use) & annot_df$gene_use != "", , drop = FALSE]

    if (nrow(annot_df) > 0) {
      gene_ranges <- annot_df %>%
        dplyr::group_by(seqnames, gene_use) %>%
        dplyr::summarise(
          start = min(start),
          end = max(end),
          .groups = "drop"
        )

      gene_gr <- GRanges(
        seqnames = gene_ranges$seqnames,
        ranges = IRanges(start = gene_ranges$start, end = gene_ranges$end)
      )
      mcols(gene_gr)$gene_name <- gene_ranges$gene_use
      mcols(gene_gr)$feature_type <- "gene_body"
      target_gr <- gene_gr
    } else {
      target_gr <- annot_gr
    }
  } else {
    target_gr <- annot_gr
  }

  hits <- distanceToNearest(dar_gr, target_gr, ignore.strand = TRUE)

  q_idx <- queryHits(hits)
  s_idx <- subjectHits(hits)
  dists <- mcols(hits)$distance
  nearest_target <- target_gr[s_idx]

  target_cols <- colnames(mcols(nearest_target))

  if ("gene_name" %in% target_cols) {
    gene_vals <- as.character(mcols(nearest_target)$gene_name)
  } else if (!is.na(gene_col) && gene_col %in% target_cols) {
    gene_vals <- as.character(mcols(nearest_target)[[gene_col]])
  } else {
    gene_vals <- rep(NA_character_, length(nearest_target))
  }

  if ("feature_type" %in% target_cols) {
    feature_vals <- as.character(mcols(nearest_target)$feature_type)
  } else if (!is.na(feature_col) && feature_col %in% target_cols) {
    feature_vals <- as.character(mcols(nearest_target)[[feature_col]])
  } else {
    feature_vals <- rep(NA_character_, length(nearest_target))
  }

  peak_level_df$row_merge_tmp <- seq_len(nrow(peak_level_df))
  peak_valid$row_merge_tmp <- q_idx
  peak_valid$nearest_gene <- gene_vals
  peak_valid$distance_to_gene <- dists
  peak_valid$nearest_feature_type <- feature_vals

  peak_level_df <- peak_level_df %>%
    dplyr::left_join(
      peak_valid %>%
        dplyr::select(row_merge_tmp, nearest_gene, distance_to_gene, nearest_feature_type),
      by = "row_merge_tmp"
    ) %>%
    dplyr::select(-row_merge_tmp)

  peak_level_df
}

############################################################
# MAIN FUNCTION
############################################################

analyze_one_cluster <- function(cluster_id) {

  cat("============================================================\n")
  cat("Running cluster ", cluster_id, "\n", sep = "")
  cat("============================================================\n")

  cluster_obj <- subset_by_cluster(obj, cluster_column, cluster_id)

  if (is.null(cluster_obj)) {
    cat("No cells found for cluster ", cluster_id, ". Skipping.\n\n", sep = "")
    return(list(
      pairwise = empty_pairwise_result_table(),
      peak_level = empty_peak_result_table(),
      summary = data.frame(
        cluster_id = as.character(cluster_id),
        cells_in_cluster = 0,
        individuals_total = 0,
        individuals_M = 0,
        individuals_F = 0,
        individuals_D = 0,
        peaks_before_filter = 0,
        peaks_after_filter = 0,
        n_M_vs_F_DARs = 0,
        n_M_vs_D_DARs = 0,
        n_F_vs_D_DARs = 0,
        n_unique_DARs = 0,
        status = "no cells",
        stringsAsFactors = FALSE
      )
    ))
  }

  meta <- cluster_obj@meta.data[, c("individual", "Status"), drop = FALSE]

  peak_counts <- get_counts_matrix(cluster_obj, atac_assay)

  stopifnot(all(colnames(peak_counts) == rownames(meta)))

  cat("Cells in cluster:", ncol(peak_counts), "\n")
  cat("Peaks in cluster object:", nrow(peak_counts), "\n")

  group_matrix <- sparse.model.matrix(~ 0 + meta$individual)
  colnames(group_matrix) <- gsub("meta\\$individual", "", colnames(group_matrix))

  pb_counts <- peak_counts %*% group_matrix
  pb_mat <- t(pb_counts)  # individuals x peaks

  cat("Pseudo-bulk matrix dimensions (individuals x peaks): ",
      nrow(pb_mat), " x ", ncol(pb_mat), "\n", sep = "")

  ind_meta <- meta %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(individual, Status) %>%
    dplyr::distinct()

  ind_meta <- as.data.frame(ind_meta)
  rownames(ind_meta) <- ind_meta$individual
  ind_meta <- ind_meta[rownames(pb_mat), , drop = FALSE]

  stopifnot(all(rownames(pb_mat) == rownames(ind_meta)))

  cat("Individuals per status in this cluster:\n")
  print(table(ind_meta$Status))
  cat("\n")

  # ----------------------------------------------------------
  # Peak support metrics across ALL individuals in this cluster
  # ----------------------------------------------------------
  peak_support_all <- data.frame(
    peak = colnames(pb_mat),
    total_counts_all_individuals = colSums(pb_mat),
    nonzero_individuals_all = colSums(pb_mat > 0),
    mean_counts_all_individuals = colMeans(pb_mat),
    max_counts_any_individual = apply(pb_mat, 2, max),
    stringsAsFactors = FALSE
  )

  keep <- peak_support_all$nonzero_individuals_all >= min_individuals_with_peak &
    peak_support_all$total_counts_all_individuals >= min_total_counts

  pb_mat_filt <- pb_mat[, keep, drop = FALSE]
  n_tests_per_comparison <- ncol(pb_mat_filt)

  cat("Peaks before filtering:", ncol(pb_mat), "\n")
  cat("Peaks after filtering:", ncol(pb_mat_filt), "\n")
  cat("Number of peaks tested per pairwise comparison: ",
      n_tests_per_comparison, "\n\n", sep = "")

  if (ncol(pb_mat_filt) == 0) {
    cat("No peaks passed filtering for cluster ", cluster_id, ". Skipping.\n\n", sep = "")
    return(list(
      pairwise = empty_pairwise_result_table(),
      peak_level = empty_peak_result_table(),
      summary = data.frame(
        cluster_id = as.character(cluster_id),
        cells_in_cluster = ncol(peak_counts),
        individuals_total = nrow(ind_meta),
        individuals_M = sum(ind_meta$Status == "M", na.rm = TRUE),
        individuals_F = sum(ind_meta$Status == "F", na.rm = TRUE),
        individuals_D = sum(ind_meta$Status == "D", na.rm = TRUE),
        peaks_before_filter = ncol(pb_mat),
        peaks_after_filter = 0,
        n_M_vs_F_DARs = 0,
        n_M_vs_D_DARs = 0,
        n_F_vs_D_DARs = 0,
        n_unique_DARs = 0,
        status = "no peaks after filter",
        stringsAsFactors = FALSE
      )
    ))
  }

  keep_groups <- ind_meta$Status %in% groups_to_test

  pb_test <- pb_mat_filt[keep_groups, , drop = FALSE]
  meta_test <- ind_meta[keep_groups, , drop = FALSE]

  meta_test$Status <- factor(meta_test$Status, levels = groups_to_test)

  cat("Individuals used for DAR testing:\n")
  print(table(meta_test$Status))
  cat("\n")

  if (any(table(meta_test$Status) < 2)) {
    cat("One or more groups has fewer than 2 individuals in cluster ",
        cluster_id, ". Skipping cluster.\n\n", sep = "")
    return(list(
      pairwise = empty_pairwise_result_table(),
      peak_level = empty_peak_result_table(),
      summary = data.frame(
        cluster_id = as.character(cluster_id),
        cells_in_cluster = ncol(peak_counts),
        individuals_total = nrow(ind_meta),
        individuals_M = sum(ind_meta$Status == "M", na.rm = TRUE),
        individuals_F = sum(ind_meta$Status == "F", na.rm = TRUE),
        individuals_D = sum(ind_meta$Status == "D", na.rm = TRUE),
        peaks_before_filter = ncol(pb_mat),
        peaks_after_filter = ncol(pb_mat_filt),
        n_M_vs_F_DARs = 0,
        n_M_vs_D_DARs = 0,
        n_F_vs_D_DARs = 0,
        n_unique_DARs = 0,
        status = "too few individuals in one or more groups",
        stringsAsFactors = FALSE
      )
    ))
  }

  # ----------------------------------------------------------
  # Peak support metrics restricted to M/F/D test individuals
  # ----------------------------------------------------------
  peak_support_test <- data.frame(
    peak = colnames(pb_test),
    total_counts_test_individuals = colSums(pb_test),
    nonzero_individuals_test = colSums(pb_test > 0),
    mean_counts_test_individuals = colMeans(pb_test),
    max_counts_test_individual = apply(pb_test, 2, max),
    stringsAsFactors = FALSE
  )

  peak_support <- peak_support_all %>%
    dplyr::left_join(peak_support_test, by = "peak")

  # ----------------------------------------------------------
  # edgeR pairwise DAR analysis
  # ----------------------------------------------------------
  y <- DGEList(counts = t(pb_test))
  y <- calcNormFactors(y)

  design <- model.matrix(~ 0 + Status, data = meta_test)
  colnames(design) <- levels(meta_test$Status)

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  contr <- makeContrasts(
    MvsF = M - F,
    MvsD = M - D,
    FvsD = F - D,
    levels = design
  )

  qlf_MvsF <- glmQLFTest(fit, contrast = contr[, "MvsF"])
  qlf_MvsD <- glmQLFTest(fit, contrast = contr[, "MvsD"])
  qlf_FvsD <- glmQLFTest(fit, contrast = contr[, "FvsD"])

  res_MvsF <- topTags(qlf_MvsF, n = Inf)$table %>%
    tibble::rownames_to_column("peak") %>%
    dplyr::mutate(comparison = "M_vs_F")

  res_MvsD <- topTags(qlf_MvsD, n = Inf)$table %>%
    tibble::rownames_to_column("peak") %>%
    dplyr::mutate(comparison = "M_vs_D")

  res_FvsD <- topTags(qlf_FvsD, n = Inf)$table %>%
    tibble::rownames_to_column("peak") %>%
    dplyr::mutate(comparison = "F_vs_D")

  all_pairwise <- dplyr::bind_rows(res_MvsF, res_MvsD, res_FvsD)

  sig_MvsF <- res_MvsF %>% dplyr::filter(FDR < fdr_cutoff)
  sig_MvsD <- res_MvsD %>% dplyr::filter(FDR < fdr_cutoff)
  sig_FvsD <- res_FvsD %>% dplyr::filter(FDR < fdr_cutoff)

  sig_peaks <- union(union(sig_MvsF$peak, sig_MvsD$peak), sig_FvsD$peak)

  cluster_summary <- data.frame(
    cluster_id = as.character(cluster_id),
    cells_in_cluster = ncol(peak_counts),
    individuals_total = nrow(ind_meta),
    individuals_M = sum(ind_meta$Status == "M", na.rm = TRUE),
    individuals_F = sum(ind_meta$Status == "F", na.rm = TRUE),
    individuals_D = sum(ind_meta$Status == "D", na.rm = TRUE),
    peaks_before_filter = ncol(pb_mat),
    peaks_after_filter = ncol(pb_mat_filt),
    n_M_vs_F_DARs = nrow(sig_MvsF),
    n_M_vs_D_DARs = nrow(sig_MvsD),
    n_F_vs_D_DARs = nrow(sig_FvsD),
    n_unique_DARs = length(sig_peaks),
    status = "ok",
    stringsAsFactors = FALSE
  )

  # ----------------------------------------------------------
  # Mean logCPM per group for classification
  # ----------------------------------------------------------
  logcpm_mat <- cpm(y, log = TRUE, prior.count = 2)

  mean_logcpm_by_group <- data.frame(
    peak = rownames(logcpm_mat),
    mean_logCPM_M = rowMeans(logcpm_mat[, meta_test$Status == "M", drop = FALSE]),
    mean_logCPM_F = rowMeans(logcpm_mat[, meta_test$Status == "F", drop = FALSE]),
    mean_logCPM_D = rowMeans(logcpm_mat[, meta_test$Status == "D", drop = FALSE]),
    stringsAsFactors = FALSE
  )

  peak_level <- all_pairwise %>%
    dplyr::select(peak, comparison, logFC, FDR) %>%
    tidyr::pivot_wider(
      names_from = comparison,
      values_from = c(logFC, FDR),
      names_glue = "{comparison}_{.value}"
    ) %>%
    dplyr::left_join(mean_logcpm_by_group, by = "peak") %>%
    dplyr::left_join(peak_support, by = "peak") %>%
    dplyr::mutate(
      sig_M_vs_F = M_vs_F_FDR < fdr_cutoff,
      sig_M_vs_D = M_vs_D_FDR < fdr_cutoff,
      sig_F_vs_D = F_vs_D_FDR < fdr_cutoff
    )

  peak_level_sig <- peak_level %>%
    dplyr::filter(sig_M_vs_F | sig_M_vs_D | sig_F_vs_D)

  if (nrow(peak_level_sig) == 0) {
    cat("No significant DARs found for cluster ", cluster_id, ".\n\n", sep = "")
    return(list(
      pairwise = empty_pairwise_result_table(),
      peak_level = empty_peak_result_table(),
      summary = cluster_summary
    ))
  }

  peak_level_sig$D_position_between_M_and_F <- (
    peak_level_sig$mean_logCPM_D - peak_level_sig$mean_logCPM_M
  ) / (
    peak_level_sig$mean_logCPM_F - peak_level_sig$mean_logCPM_M
  )

  peak_level_sig$classification <- mapply(
    classify_one_peak,
    mean_M = peak_level_sig$mean_logCPM_M,
    mean_F = peak_level_sig$mean_logCPM_F,
    mean_D = peak_level_sig$mean_logCPM_D,
    sig_MF = peak_level_sig$sig_M_vs_F,
    sig_MD = peak_level_sig$sig_M_vs_D,
    sig_FD = peak_level_sig$sig_F_vs_D,
    MoreArgs = list(
      progressive_lower = 0.25,
      progressive_upper = 0.75
    )
  )

  peak_level_sig$cluster_id <- as.character(cluster_id)

  peak_level_sig <- annotate_peak_table(cluster_obj, peak_level_sig, atac_assay = atac_assay)

  peak_level_sig <- peak_level_sig %>%
    dplyr::select(
      cluster_id,
      peak,
      chromosome, chromosome_norm, start, end,
      nearest_gene, distance_to_gene, nearest_feature_type,
      mean_logCPM_M, mean_logCPM_F, mean_logCPM_D,
      M_vs_F_logFC, M_vs_F_FDR,
      M_vs_D_logFC, M_vs_D_FDR,
      F_vs_D_logFC, F_vs_D_FDR,
      sig_M_vs_F, sig_M_vs_D, sig_F_vs_D,
      D_position_between_M_and_F,
      classification,
      total_counts_all_individuals,
      nonzero_individuals_all,
      mean_counts_all_individuals,
      max_counts_any_individual,
      total_counts_test_individuals,
      nonzero_individuals_test,
      mean_counts_test_individuals,
      max_counts_test_individual
    ) %>%
    dplyr::arrange(classification, pmin(M_vs_F_FDR, M_vs_D_FDR, F_vs_D_FDR, na.rm = TRUE))

  cat("Peak-level classified DARs kept:", nrow(peak_level_sig), "\n")
  print(table(peak_level_sig$classification))
  cat("\n")

  # Keep ALL significant pairwise rows for those classified peaks
  all_pairwise_sig <- all_pairwise %>%
    dplyr::filter(FDR < fdr_cutoff) %>%
    dplyr::filter(peak %in% peak_level_sig$peak) %>%
    dplyr::left_join(
      peak_level_sig %>%
        dplyr::select(
          cluster_id, peak,
          chromosome, chromosome_norm, start, end,
          nearest_gene, distance_to_gene, nearest_feature_type,
          classification,
          mean_logCPM_M, mean_logCPM_F, mean_logCPM_D,
          total_counts_all_individuals,
          nonzero_individuals_all,
          mean_counts_all_individuals,
          max_counts_any_individual,
          total_counts_test_individuals,
          nonzero_individuals_test,
          mean_counts_test_individuals,
          max_counts_test_individual
        ),
      by = "peak"
    ) %>%
    dplyr::select(
      cluster_id,
      peak,
      chromosome, chromosome_norm, start, end,
      nearest_gene, distance_to_gene, nearest_feature_type,
      comparison, logFC, logCPM, PValue, FDR,
      classification,
      mean_logCPM_M, mean_logCPM_F, mean_logCPM_D,
      total_counts_all_individuals,
      nonzero_individuals_all,
      mean_counts_all_individuals,
      max_counts_any_individual,
      total_counts_test_individuals,
      nonzero_individuals_test,
      mean_counts_test_individuals,
      max_counts_test_individual
    ) %>%
    dplyr::arrange(FDR)

  list(
    pairwise = all_pairwise_sig,
    peak_level = peak_level_sig,
    summary = cluster_summary
  )
}

############################################################
# RUN ALL CLUSTERS
############################################################

all_pairwise_results <- list()
all_peak_results <- list()
all_cluster_summaries <- list()

for (cl in cluster_ids) {
  res <- tryCatch(
    analyze_one_cluster(cl),
    error = function(e) {
      cat("ERROR in cluster ", cl, ": ", conditionMessage(e), "\n\n", sep = "")
      list(
        pairwise = empty_pairwise_result_table(),
        peak_level = empty_peak_result_table(),
        summary = data.frame(
          cluster_id = as.character(cl),
          cells_in_cluster = NA_integer_,
          individuals_total = NA_integer_,
          individuals_M = NA_integer_,
          individuals_F = NA_integer_,
          individuals_D = NA_integer_,
          peaks_before_filter = NA_integer_,
          peaks_after_filter = NA_integer_,
          n_M_vs_F_DARs = 0,
          n_M_vs_D_DARs = 0,
          n_F_vs_D_DARs = 0,
          n_unique_DARs = 0,
          status = paste0("error: ", conditionMessage(e)),
          stringsAsFactors = FALSE
        )
      )
    }
  )

  all_pairwise_results[[as.character(cl)]] <- res$pairwise
  all_peak_results[[as.character(cl)]] <- res$peak_level
  all_cluster_summaries[[as.character(cl)]] <- res$summary
}

pairwise_table <- dplyr::bind_rows(all_pairwise_results) %>%
  dplyr::distinct(
    cluster_id, peak, comparison, logFC, logCPM, PValue, FDR,
    .keep_all = TRUE
  )

peak_table <- dplyr::bind_rows(all_peak_results) %>%
  dplyr::distinct(cluster_id, peak, .keep_all = TRUE)

cluster_summary_table <- dplyr::bind_rows(all_cluster_summaries)

cat("============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n")
cat("Total pairwise significant DAR rows across all clusters:", nrow(pairwise_table), "\n")
cat("Total unique classified DAR peaks across all clusters:", nrow(peak_table), "\n")
cat("Clusters contributing at least one classified DAR:",
    length(unique(peak_table$cluster_id[!is.na(peak_table$cluster_id) & peak_table$cluster_id != ""])),
    "\n\n")

cat("Overall classification counts:\n")
print(table(peak_table$classification))
cat("\n")

############################################################
# EXPORT TABLES
############################################################

pairwise_export <- pairwise_table
peak_export <- peak_table
summary_export <- cluster_summary_table

# Numeric formatting helper for export
format_cols <- function(df, cols, digits) {
  cols <- intersect(cols, colnames(df))
  for (cc in cols) {
    df[[cc]] <- safe_num_fmt(df[[cc]], digits)
  }
  df
}

pairwise_export <- format_cols(
  pairwise_export,
  c("PValue", "FDR", "logFC", "logCPM",
    "mean_logCPM_M", "mean_logCPM_F", "mean_logCPM_D",
    "mean_counts_all_individuals", "mean_counts_test_individuals"),
  10
)

peak_export <- format_cols(
  peak_export,
  c("mean_logCPM_M", "mean_logCPM_F", "mean_logCPM_D",
    "M_vs_F_logFC", "M_vs_D_logFC", "F_vs_D_logFC",
    "M_vs_F_FDR", "M_vs_D_FDR", "F_vs_D_FDR",
    "D_position_between_M_and_F",
    "mean_counts_all_individuals", "mean_counts_test_individuals"),
  10
)

# Use extra precision for FDR/PValue columns
pairwise_export <- format_cols(pairwise_export, c("PValue", "FDR"), 15)
peak_export <- format_cols(peak_export, c("M_vs_F_FDR", "M_vs_D_FDR", "F_vs_D_FDR"), 15)

write.table(
  peak_export,
  file = outfile_peak_csv,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "",
  eol = "\n"
)

write.table(
  pairwise_export,
  file = outfile_pairwise_csv,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "",
  eol = "\n"
)

write.csv(
  summary_export,
  file = outfile_summary_csv,
  row.names = FALSE
)

cat("Peak-level classified DAR table saved to:\n")
cat(normalizePath(outfile_peak_csv, mustWork = FALSE), "\n\n")

cat("Pairwise DAR table with classification saved to:\n")
cat(normalizePath(outfile_pairwise_csv, mustWork = FALSE), "\n\n")

cat("Cluster summary table saved to:\n")
cat(normalizePath(outfile_summary_csv, mustWork = FALSE), "\n\n")

############################################################
# STACKED BAR PLOT
############################################################

classification_levels <- c(
  "Early Downregulated",
  "Early Upregulated",
  "Late Downregulated",
  "Late Upregulated",
  "Progressively Downregulated",
  "Progressively Upregulated",
  "Transiently Downregulated",
  "Transiently Upregulated",
  "Unclassified"
)

classification_colors <- c(
  "Early Downregulated" = "red",
  "Early Upregulated" = "green4",
  "Late Downregulated" = "blue",
  "Late Upregulated" = "black",
  "Progressively Downregulated" = "darkorchid",
  "Progressively Upregulated" = "grey70",
  "Transiently Downregulated" = "brown3",
  "Transiently Upregulated" = "orange",
  "Unclassified" = "pink"
)

if (nrow(peak_table) > 0) {

  plot_df <- peak_table %>%
    dplyr::filter(!is.na(classification)) %>%
    dplyr::mutate(
      cluster_id_num = as.numeric(as.character(cluster_id)),
      classification = factor(classification, levels = classification_levels)
    ) %>%
    dplyr::count(cluster_id, cluster_id_num, classification, name = "n_DARs") %>%
    dplyr::filter(n_DARs > 0) %>%
    dplyr::arrange(cluster_id_num)

  clusters_with_dars <- sort(unique(plot_df$cluster_id_num))

  plot_df$cluster_id <- factor(
    plot_df$cluster_id,
    levels = as.character(clusters_with_dars)
  )

  p <- ggplot(plot_df, aes(x = cluster_id, y = n_DARs, fill = classification)) +
    geom_bar(stat = "identity", width = 0.85) +
    scale_fill_manual(values = classification_colors, drop = FALSE) +
    labs(
      title = "DAR classification across clusters",
      subtitle = paste0("Filter: nonzero individuals >= ", min_individuals_with_peak,
                        ", total counts >= ", min_total_counts,
                        ", FDR < ", fdr_cutoff),
      x = "Cluster ID",
      y = "DARs",
      fill = "Classification"
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  ggsave(
    filename = outfile_plot_png,
    plot = p,
    width = 14,
    height = 6,
    dpi = 300
  )

  print(p)

  cat("Stacked DAR classification plot saved to:\n")
  cat(normalizePath(outfile_plot_png, mustWork = FALSE), "\n\n")

} else {
  cat("No DARs found across clusters, so no stacked bar plot was created.\n\n")
}

cat("All done. Files in output directory:\n")
print(list.files(outdir, full.names = TRUE))
