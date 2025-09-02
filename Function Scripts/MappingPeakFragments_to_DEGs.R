MappingPeakFragments_to_DEGs <- function(
    testing_variable_name, # cluster/cell-type name or name of differential testing variable (example: testing_variable_name = "neuroendocrine_clust6")
    deg_output, 
    seurat_object,
    gene_homolog_df,
    genome_annotation_df, 
    range_extensions = c(1000, 5000, 10000, 25000, 50000)  # Customizable range extensions
) {
  
  
  
  colnames(deg_output)[colnames(deg_output) == 'gene'] = 'gene_id'
  # Load the GLMMSeq DEG outputs
  gene_coordinates <- genome_annotation_df %>% dplyr::filter(type =="gene")
  
  deg_output$seqnames <- gene_coordinates$seqnames[match(deg_output$gene_id, gene_coordinates$gene_id)]
  deg_output$gene_start <- gene_coordinates$start[match(deg_output$gene_id, gene_coordinates$gene_id)]
  deg_output$gene_end <- gene_coordinates$end[match(deg_output$gene_id, gene_coordinates$gene_id)]
  deg_output$gene_width <- gene_coordinates$width[match(deg_output$gene_id, gene_coordinates$gene_id)]
  deg_output$gene_strand <- gene_coordinates$strand[match(deg_output$gene_id, gene_coordinates$gene_id)]
  deg_output$gene_biotype <- gene_coordinates$gene_biotype[match(deg_output$gene_id, gene_coordinates$gene_id)]
  deg_output$LinkageGroup <- gene_coordinates$standard_name[match(deg_output$gene_id, gene_coordinates$gene_id)]
  
  deg_output$gene_coordinates <- paste(deg_output$seqnames, deg_output$gene_start, deg_output$gene_end, sep = c("-", "-"))
  
  deg_output$gene <- deg_output$gene_id
  deg_output$m.zebra_description <- gene_homolog_df$mzebra_description[match(deg_output$gene_id, gene_homolog_df$seurat_name)]
  deg_output$greatest_human_homolog <- gene_homolog_df$human[match(deg_output$gene_id, gene_homolog_df$seurat_name)]
  deg_output$one_to_one_homologs <- gene_homolog_df$human_agree[match(deg_output$gene_id, gene_homolog_df$seurat_name)]
  
  
  gc()
  gc()
  
  
  # Extract ATAC peak coordinates
  max_splits <- 3
  peak_coordinates <- as.data.frame(unique(seurat_object[['ATAC']]@ranges))
  peak_coordinates2 <- as.data.frame(unique(seurat_object[['ATAC']]@ranges))
  peak_coordinates2$fragment <- paste0(peak_coordinates2$seqnames, '_',peak_coordinates2$start, '-', peak_coordinates2$end)
  peak_coordinates2$seqnames <- NULL
  peak_coordinates2$start <- NULL
  peak_coordinates2$end <- NULL
  peak_coordinates2$width <- NULL
  peak_coordinates2$strand <- NULL
  peak_coordinates <- cbind(peak_coordinates2, peak_coordinates)
  colnames(peak_coordinates)[1] <- 'fragment'
  
  
  # Apply function for peaks within DEG coordinates (ie. when range extension = 0)
  add_peak_overlap_columns <- function(deg_output) {
    # Helper function to find matching peaks within a given distance
    find_matching_peaks <- function(gene_row, start, end, range_extension) {
      start <- as.numeric(gene_row[["gene_start"]]) - range_extension
      end <- as.numeric(gene_row[["gene_end"]]) + range_extension
      
      matching_peaks <- peak_coordinates$fragment[
        (peak_coordinates$seqnames == gene_row["seqnames"] & peak_coordinates$start >= start & peak_coordinates$start <= end) |
          (peak_coordinates$seqnames == gene_row["seqnames"] & peak_coordinates$end >= start & peak_coordinates$end <= end) |
          (peak_coordinates$seqnames == gene_row["seqnames"] & peak_coordinates$start <= start & peak_coordinates$end >= end)
      ]
      
      gc()  # Perform garbage collection
      paste(matching_peaks, collapse = ",")
    }
    
    # Apply the function to each row and add two new columns
    deg_output$peaks_within_DEG_coordinates <- apply(deg_output, 1, find_matching_peaks, start = 'gene_start', end = 'gene_end', range_extension = 0)
    
    deg_output$peakCount_within_DEG_coordinates <- sapply(strsplit(as.character(deg_output$peaks_within_DEG_coordinates), ","), function(x) ifelse(all(x == ""), 0, length(x)))
    
    gc()
    gc()
    
    return(deg_output)
  }
  
  deg_output <- add_peak_overlap_columns(deg_output)
  
  peaks_within_DEG_coordinates <- as.data.frame(unique(unlist(strsplit(deg_output$peaks_within_DEG_coordinates, ","))))
  colnames(peaks_within_DEG_coordinates)[[1]] <- 'fragment'
  
  find_matching_peaks <- function(gene_row, range_extension) {
    start <- as.numeric(gene_row[["gene_start"]]) - range_extension
    end <- as.numeric(gene_row[["gene_end"]]) + range_extension
    
    matching_peaks <- peak_coordinates$fragment[
      (peak_coordinates$seqnames == gene_row["seqnames"] & peak_coordinates$start >= start & peak_coordinates$start <= end) |
        (peak_coordinates$seqnames == gene_row["seqnames"] & peak_coordinates$end >= start & peak_coordinates$end <= end) |
        (peak_coordinates$seqnames == gene_row["seqnames"] & peak_coordinates$start <= start & peak_coordinates$end >= end)
    ]
    
    gc()  # Perform garbage collection
    paste(matching_peaks, collapse = ",")
  }
  
  # Apply function for each range extension
  for (range_extension in range_extensions) {
    column_name <- paste0("peaks_within_", range_extension, "bp")
    count_column_name <- paste0("peakCount_within_", range_extension, "bp")
    
    # Apply function to find matching peaks
    deg_output[[column_name]] <- apply(deg_output, 1, find_matching_peaks, range_extension = range_extension)
    deg_output[[count_column_name]] <- sapply(strsplit(as.character(deg_output[[column_name]]), ","), function(x) ifelse(all(x == ""), 0, length(x)))
  }
  
  # Organize output
  peak_results <- lapply(range_extensions, function(ext) {
    peaks <- unique(unlist(strsplit(deg_output[[paste0("peaks_within_", ext, "bp")]], ",")))
    peaks <- peaks[peaks != "NA"]
    data.frame(fragment = peaks)
  })
  
  names(peak_results) <- paste0('ALL_', testing_variable_name,  "_peaks_within_", range_extensions, "bp")
  
  result_list <- c(list(Fragments_Mapped_to_DEGs = deg_output, peaks_within_DEG_coordinates = peaks_within_DEG_coordinates), peak_results)
  
  names(result_list)[1] <- paste0(testing_variable_name, '_', names(result_list)[1])
  
  names(result_list)[2] <- paste0('ALL_', testing_variable_name, '_peaks_in_DEG_coordinates')

  return(result_list)
  
}

gene_homolog_df = readRDS("Collaboration/clownfish_gene_info_homologs.rds")
genome_annotation_df = readRDS("Collaboration/clownfish_genome_annotations.rds")
gene_coordinates <- genome_annotation_df[genome_annotation_df$type == "gene"]
gene_coordinates <- as.data.frame(gene_coordinates)

saveRDS(MappingPeakFragments_to_DEGs, 'Functions/MappingPeakFragments_to_DEGs.rds')

