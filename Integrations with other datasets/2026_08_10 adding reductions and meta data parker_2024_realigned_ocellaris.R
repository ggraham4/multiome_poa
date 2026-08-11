library(Seurat)

# ============================================================
# Step 0: Load inputs
#   parker_2024:        original 2021 object (a. percula-aligned counts,
#                        published UMAP, clusters, and metadata)
#   renamed_obj_fixed:   output of script 01 -- the 2025 a. ocellaris
#                        realignment, subset + renamed to Parker's barcodes
# ============================================================
parker_2024       <- readRDS("A:/Anemonefish POA Legacy R Objects/Old object.rds")
renamed_obj_fixed <- readRDS( "A:/Anemonefish POA Legacy R Objects/renamed_obj_fixed_barcodes.rds")

# ============================================================
# Step 1: Restrict Parker's object to cells that also survived the
#   2025 realignment/QC. This intentionally keeps Parker's original
#   cell set (and therefore its UMAP/cluster identities) as the frame
#   of reference -- see conversation for why this is preferred over
#   trying to fold in the realignment's "new" cells.
# ============================================================
common_cells <- intersect(colnames(parker_2024), colnames(renamed_obj_fixed))

cat("Cells in Parker 2024:            ", ncol(parker_2024), "\n")
cat("Cells in realigned (fixed) object:", ncol(renamed_obj_fixed), "\n")
cat("Cells common to both:             ", length(common_cells), "\n")

parker_filtered <- subset(parker_2024, cells = common_cells)

# ============================================================
# Step 2: Pull the matching ocellaris-aligned counts/data, in
#   parker_filtered's exact cell order
# ============================================================
new_counts <- GetAssayData(renamed_obj_fixed, assay = "RNA", layer = "counts")[, colnames(parker_filtered)]
new_data   <- GetAssayData(renamed_obj_fixed, assay = "RNA", layer = "data")[, colnames(parker_filtered)]

stopifnot(identical(colnames(new_counts), colnames(parker_filtered)))
stopifnot(identical(colnames(new_data),   colnames(parker_filtered)))

cat("\nGenes -- old (a. percula ref):", nrow(parker_filtered),
    "| new (a. ocellaris ref):", nrow(new_counts), "\n")

# ============================================================
# Step 3: Swap in the realigned counts/data. UMAP, PCA, and cluster
#   identities in parker_filtered are left untouched -- they still
#   describe these exact cells, just now backed by ocellaris-aligned
#   expression values instead of percula-aligned ones.
# ============================================================
parker_filtered[["RNA"]] <- CreateAssayObject(counts = new_counts)
parker_filtered <- SetAssayData(parker_filtered, assay = "RNA", layer = "data", new.data = new_data)

# Optional: original nCount_RNA/nFeature_RNA reflect the OLD percula
# alignment and are now stale for these cells. Uncomment to add
# ocellaris-based versions alongside (without overwriting the originals):
# parker_filtered$nCount_RNA_ocellaris   <- colSums(new_counts)
# parker_filtered$nFeature_RNA_ocellaris <- colSums(new_counts > 0)

# sanity check -- UMAP should look identical to the original paper's plot
DimPlot(parker_filtered, reduction = "umap")

# ============================================================
# Step 4: Save
# ============================================================
saveRDS(parker_filtered,  "A:/Anemonefish POA Legacy R Objects/parker_2024_realigned_ocellaris.rds")
DimPlot(parker_filtered, group.by = 'clusters49')
FeaturePlot(parker_filtered, 'esr2b')
