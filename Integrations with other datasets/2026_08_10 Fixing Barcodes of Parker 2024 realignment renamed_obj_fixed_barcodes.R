library(Seurat)

# ============================================================
# Step 0: Load inputs
# ============================================================
parker_2024 <- readRDS("A:/Anemonefish POA Legacy R Objects/Old object.rds")
renamed_obj <- readRDS("A:/Anemonefish POA Legacy R Objects/allcells_filtered_2025-11-13.Rds")

# ============================================================
# Step 1: Determine the sample-to-sample mapping empirically
#   (renamed_obj's orig.ident labels -> parker's $sample labels)
#   via raw 10x barcode overlap -- this doesn't rely on any naming
#   convention lining up, only on the physical barcode sequences.
# ============================================================
parker_barcodes <- rownames(parker_2024@meta.data)
parker_core <- sub("^[a-zA-Z]_", "", parker_barcodes)   # drop leading sex_ prefix
parker_core <- sub("_[0-9]+$", "", parker_core)          # drop trailing merge suffix
parker_core_by_sample <- split(parker_core, parker_2024$sample)

renamed_core <- colnames(renamed_obj)
# renamed_obj's own "-1/-2/-3/-4" suffix reflects aggregation order, not sample
# identity once we're grouping by orig.ident -- normalize it away before matching
renamed_core_fixed <- sub("-[0-9]+$", "-1", renamed_core)
renamed_core_fixed_by_sample <- split(renamed_core_fixed, renamed_obj$orig.ident)
renamed_core_by_sample <- split(renamed_core, renamed_obj$orig.ident)

samples_new <- names(renamed_core_by_sample)
samples_old <- names(parker_core_by_sample)

overlap_mat <- matrix(0, nrow = length(samples_new), ncol = length(samples_old),
                      dimnames = list(samples_new, samples_old))
for (i in samples_new) {
  for (j in samples_old) {
    overlap_mat[i, j] <- length(intersect(renamed_core_fixed_by_sample[[i]],
                                          parker_core_by_sample[[j]]))
  }
}

cat("Sample overlap matrix (rows = renamed_obj orig.ident, cols = parker $sample):\n")
print(overlap_mat)

best_match <- apply(overlap_mat, 1, function(row) names(which.max(row)))
cat("\nInferred sample mapping:\n")
print(best_match)

# dominance check -- confirms each row's top match isn't a coin-flip
top2 <- apply(overlap_mat, 1, function(row) sort(row, decreasing = TRUE)[1:2])
cat("\nTop-2 overlap counts per sample (want a big gap between rows):\n")
print(top2)

# ============================================================
# Step 2: Build a cell-level lookup: renamed_obj colname -> parker colname
# ============================================================
matched_list <- lapply(samples_new, function(ns) {
  os <- best_match[[ns]]
  new_bcs_fixed  <- renamed_core_fixed_by_sample[[ns]]
  new_bcs_actual <- renamed_core_by_sample[[ns]]
  old_bcs_full   <- parker_barcodes[parker_2024$sample == os]
  old_core       <- parker_core_by_sample[[os]]
  
  idx <- match(new_bcs_fixed, old_core)
  data.frame(new_barcode     = new_bcs_actual[!is.na(idx)],
             parker_barcode  = old_bcs_full[idx[!is.na(idx)]],
             stringsAsFactors = FALSE)
})
matched_df <- do.call(rbind, matched_list)

cat("\nTotal matched cells:", nrow(matched_df), "of", ncol(renamed_obj), "renamed_obj cells\n")
cat("Duplicate new_barcodes (want 0):", sum(duplicated(matched_df$new_barcode)), "\n")
cat("Duplicate parker_barcodes (want 0):", sum(duplicated(matched_df$parker_barcode)), "\n")

# ============================================================
# Step 3: Subset renamed_obj to matched cells and rename to Parker's barcodes
#   (this makes the object trivially joinable to parker_2024 downstream --
#   see script 02)
# ============================================================
renamed_obj_fixed <- subset(renamed_obj, cells = matched_df$new_barcode)

# reorder the lookup to match the subset object's actual cell order
matched_df <- matched_df[match(colnames(renamed_obj_fixed), matched_df$new_barcode), ]
stopifnot(identical(colnames(renamed_obj_fixed), matched_df$new_barcode))

# keep the pre-rename barcode as metadata for provenance/traceability
renamed_obj_fixed$original_barcode <- matched_df$new_barcode

renamed_obj_fixed <- RenameCells(renamed_obj_fixed, new.names = matched_df$parker_barcode)

cat("\nrenamed_obj_fixed:", ncol(renamed_obj_fixed), "cells, barcodes now match Parker's naming:\n")
print(head(colnames(renamed_obj_fixed)))

# ============================================================
# Step 4: Save
# ============================================================
saveRDS(renamed_obj_fixed, "A:/Anemonefish POA Legacy R Objects/renamed_obj_fixed_barcodes.rds")
