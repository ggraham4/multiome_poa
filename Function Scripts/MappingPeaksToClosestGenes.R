MappingPeaksToClosestGenes <- function(
    peak_fragments,  # Vector of peak fragment names (e.g., "chr1_1000-2000")
    seurat_object,
    gene_homolog_df,
    genome_annotation_df,
    max_distance = 100000  # Maximum distance to search for genes (default 100kb)
) {
  
  # Extract gene coordinates
  gene_coordinates <- genome_annotation_df %>% dplyr::filter(type == "gene")
  gene_coordinates <- as.data.frame(gene_coordinates)
  
  # Get fragment names directly from ATAC assay rownames
  all_fragments <- rownames(seurat_object[['ATAC']])
  
  # Extract ATAC peak coordinates from Seurat object
  peak_coordinates <- as.data.frame(unique(seurat_object[['ATAC']]@ranges))
  peak_coordinates$fragment <- all_fragments[1:nrow(peak_coordinates)]
  
  cat("First 5 input fragments:", paste(head(peak_fragments, 5), collapse = ", "), "\n")
  cat("First 5 Seurat fragments:", paste(head(peak_coordinates$fragment, 5), collapse = ", "), "\n")
  
  # Filter peak coordinates to only include the requested fragments
  peak_coordinates <- peak_coordinates[peak_coordinates$fragment %in% peak_fragments, ]
  
  if (nrow(peak_coordinates) == 0) {
    cat("Available fragments (first 10):\n")
    print(head(unique(as.data.frame(seurat_object[['ATAC']]@ranges))$fragment, 10))
    stop("No matching peak fragments found in the Seurat object. Check fragment naming format.")
  }
  
  # Function to find the closest gene for a single peak
  find_closest_gene <- function(peak_row) {
    peak_chr <- peak_row[["seqnames"]]
    peak_start <- as.numeric(peak_row[["start"]])
    peak_end <- as.numeric(peak_row[["end"]])
    peak_center <- round((peak_start + peak_end) / 2)
    
    # Filter genes on the same chromosome
    chr_genes <- gene_coordinates[gene_coordinates$seqnames == peak_chr, ]
    
    if (nrow(chr_genes) == 0) {
      return(data.frame(
        closest_gene_id = NA,
        closest_gene_name = NA,
        distance_to_gene = NA,
        gene_start = NA,
        gene_end = NA,
        gene_strand = NA,
        gene_biotype = NA,
        position_relative_to_gene = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    # Calculate distances from peak center to gene boundaries
    distances <- pmin(
      abs(peak_center - chr_genes$start),  # Distance to gene start
      abs(peak_center - chr_genes$end)     # Distance to gene end
    )
    
    # Find genes that overlap with the peak
    overlapping_genes <- which(
      chr_genes$start <= peak_end & chr_genes$end >= peak_start
    )
    
    if (length(overlapping_genes) > 0) {
      # If peak overlaps with genes, choose the one with maximum overlap
      overlaps <- pmin(chr_genes$end[overlapping_genes], peak_end) - 
                 pmax(chr_genes$start[overlapping_genes], peak_start)
      closest_idx <- overlapping_genes[which.max(overlaps)]
      distance <- 0  # Overlapping distance is 0
      position <- "overlapping"
    } else {
      # Find the closest non-overlapping gene
      closest_idx <- which.min(distances)
      distance <- distances[closest_idx]
      
      # Filter by maximum distance if specified
      if (!is.null(max_distance) && distance > max_distance) {
        return(data.frame(
          closest_gene_id = NA,
          closest_gene_name = NA,
          distance_to_gene = distance,
          gene_start = NA,
          gene_end = NA,
          gene_strand = NA,
          gene_biotype = NA,
          position_relative_to_gene = "too_far",
          stringsAsFactors = FALSE
        ))
      }
      
      # Determine position relative to gene
      gene_start <- chr_genes$start[closest_idx]
      gene_end <- chr_genes$end[closest_idx]
      
      if (peak_center < gene_start) {
        position <- "upstream"
      } else if (peak_center > gene_end) {
        position <- "downstream"
      } else {
        position <- "within_gene"  # Should not happen for non-overlapping
      }
    }
    
    closest_gene <- chr_genes[closest_idx, ]
    
    return(data.frame(
      closest_gene_id = closest_gene$gene_id,
      closest_gene_name = closest_gene$gene_name,
      distance_to_gene = distance,
      gene_start = closest_gene$start,
      gene_end = closest_gene$end,
      gene_strand = closest_gene$strand,
      gene_biotype = closest_gene$gene_biotype,
      position_relative_to_gene = position,
      stringsAsFactors = FALSE
    ))
  }
  
  # Apply the function to each peak
  cat("Finding closest genes for", nrow(peak_coordinates), "peaks...\n")
  
  closest_genes_list <- apply(peak_coordinates, 1, find_closest_gene)
  closest_genes_df <- do.call(rbind, closest_genes_list)
  
  # Combine with original peak information
  result_df <- cbind(peak_coordinates, closest_genes_df)
  
  # Add homolog information
  result_df$m.zebra_description <- gene_homolog_df$mzebra_description[
    match(result_df$closest_gene_id, gene_homolog_df$seurat_name)
  ]
  result_df$greatest_human_homolog <- gene_homolog_df$human[
    match(result_df$closest_gene_id, gene_homolog_df$seurat_name)
  ]
  result_df$one_to_one_homologs <- gene_homolog_df$human_agree[
    match(result_df$closest_gene_id, gene_homolog_df$seurat_name)
  ]
  
  # Add gene coordinates for convenience
  result_df$gene_coordinates <- paste(result_df$seqnames, 
                                     result_df$gene_start, 
                                     result_df$gene_end, 
                                     sep = c("-", "-"))
  
  # Reorder columns for better readability - only include columns that exist
  available_cols <- colnames(result_df)
  desired_cols <- c(
    "fragment", "seqnames", "start", "end", "width", "strand",
    "closest_gene_id", "closest_gene_name", "distance_to_gene", 
    "position_relative_to_gene", "gene_start", "gene_end", "gene_strand", 
    "gene_biotype", "gene_coordinates", "m.zebra_description", 
    "greatest_human_homolog", "one_to_one_homologs"
  )
  
  # Only select columns that actually exist
  cols_to_select <- desired_cols[desired_cols %in% available_cols]
  result_df <- result_df[, cols_to_select, drop = FALSE]
  
  # Create summary statistics
  summary_stats <- list(
    total_peaks = nrow(result_df),
    peaks_with_genes = sum(!is.na(result_df$closest_gene_id)),
    peaks_overlapping_genes = sum(result_df$position_relative_to_gene == "overlapping", na.rm = TRUE),
    peaks_upstream = sum(result_df$position_relative_to_gene == "upstream", na.rm = TRUE),
    peaks_downstream = sum(result_df$position_relative_to_gene == "downstream", na.rm = TRUE),
    peaks_too_far = sum(result_df$position_relative_to_gene == "too_far", na.rm = TRUE),
    median_distance = median(result_df$distance_to_gene[result_df$distance_to_gene > 0], na.rm = TRUE),
    mean_distance = mean(result_df$distance_to_gene[result_df$distance_to_gene > 0], na.rm = TRUE)
  )
  
  # Cleanup
  gc()
  
  cat("Mapping completed!\n")
  cat("Summary:\n")
  cat("- Total peaks:", summary_stats$total_peaks, "\n")
  cat("- Peaks with assigned genes:", summary_stats$peaks_with_genes, "\n")
  cat("- Overlapping with genes:", summary_stats$peaks_overlapping_genes, "\n")
  cat("- Upstream of genes:", summary_stats$peaks_upstream, "\n")
  cat("- Downstream of genes:", summary_stats$peaks_downstream, "\n")
  if (summary_stats$peaks_too_far > 0) {
    cat("- Too far from genes:", summary_stats$peaks_too_far, "\n")
  }
  if (!is.na(summary_stats$median_distance)) {
    cat("- Median distance to genes:", round(summary_stats$median_distance), "bp\n")
  }
  
  return(list(
    peak_to_gene_mapping = result_df,
    summary_statistics = summary_stats
  ))
}

# Example usage:
# Load required data
# gene_homolog_df <- readRDS("Collaboration/clownfish_gene_info_homologs.rds")
# genome_annotation_df <- readRDS("Collaboration/clownfish_genome_annotations.rds")

# Example peak fragments
# peak_fragments <- c("chr1_1000-2000", "chr2_5000-6000", "chr3_10000-11000")

# Run the function
# result <- MappingPeaksToClosestGenes(
#   peak_fragments = peak_fragments,
#   seurat_object = your_seurat_object,
#   gene_homolog_df = gene_homolog_df,
#   genome_annotation_df = genome_annotation_df,
#   max_distance = 50000  # Only consider genes within 50kb
# )

# Access results
# peak_mappings <- result$peak_to_gene_mapping
# summary_info <- result$summary_statistics

# Save the function
saveRDS(MappingPeaksToClosestGenes, 'Functions/MappingPeaksToClosestGenes.rds')
