### old object peaks
library(GenomicRanges)

obj = readRDS("~/Desktop/legacy_multiome.rds")

obj@assays$peaks@fragments

### here, make sure you run this every time you connec tto the UIUC vpn and open it
for(i in 1:6){
obj@assays[["peaks"]]@fragments[[i]]@path = paste0("/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/cellranger_outs_atac_fragment_files/s",i,"_clownfish/atac_fragments.tsv.gz")
}

b = (obj@assays[['peaks']]@annotation)


CoveragePlot(
  object = obj,
  region = "syt1a",  assay = 'peaks', )

obj = JoinLayers(obj)

obj$neuron = ifelse(obj$seurat_clusters %in% c('3','27'),'Glia', 'other')

obj$neuron = ifelse(obj$seurat_clusters %in% c('0'),'neuron', obj$neuron)

sum(obj$neuron =='neuron')
sum(obj$neuron =='Glia')

Idents(obj) <- 'neuron'
CoveragePlot(
  object = obj,
  region = "NC_072786.1-4122587-4124587",   #syt1a   TSS
  #features = "syt1a",
  expression.assay = "RNA",
  #extend.upstream = 0,
  #extend.downstream = 100000,
  assay = 'peaks', 
  sep = c('-', '-'), 
  peaks = T,
  annotation = F)

specific_range <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = 1000000, end = 1005000)
)


FeaturePlot(obj, reduction = 'harmony_wnn.umap', feature = 'syt1a')


obj@assays[["ATAC"]]@fragments[[1]]@path <- "~/Desktop/ATAC Fragment Files/cellranger_outs_atac_fragment_files/s1_clownfish/atac_fragments.tsv.gz"
obj@assays[["ATAC"]]@fragments[[2]]@path <- "~/Desktop/ATAC Fragment Files/cellranger_outs_atac_fragment_files/s2_clownfish/atac_fragments.tsv.gz"
obj@assays[["ATAC"]]@fragments[[3]]@path <- "~/Desktop/ATAC Fragment Files/cellranger_outs_atac_fragment_files/s3_clownfish/atac_fragments.tsv.gz"
obj@assays[["ATAC"]]@fragments[[4]]@path <- "~/Desktop/ATAC Fragment Files/cellranger_outs_atac_fragment_files/s4_clownfish/atac_fragments.tsv.gz"
obj@assays[["ATAC"]]@fragments[[5]]@path <- "~/Desktop/ATAC Fragment Files/cellranger_outs_atac_fragment_files/s5_clownfish/atac_fragments.tsv.gz"
obj@assays[["ATAC"]]@fragments[[6]]@path <- "~/Desktop/ATAC Fragment Files/cellranger_outs_atac_fragment_files/s6_clownfish/atac_fragments.tsv.gz"


head(rownames(obj@assays$ATAC@data))
head(rownames(obj@assays$ATAC@counts))
head(rownames(obj@assays$ATAC@meta.features))
head(obj[['ATAC']]@var.features)

rownames(obj@assays$ATAC@data)   = stringr::str_replace(rownames(obj@assays$ATAC@data),   "-", "_")
rownames(obj@assays$ATAC@counts) = stringr::str_replace(rownames(obj@assays$ATAC@counts), "-", "_")
rownames(obj@assays$ATAC@meta.features) = stringr::str_replace(rownames(obj@assays$ATAC@meta.features), "-", "_")
obj[['ATAC']]@var.features = stringr::str_replace(obj[['ATAC']]@var.features, "-", "_")

head(rownames(obj@assays$ATAC@data))
head(rownames(obj@assays$ATAC@counts))
head(rownames(obj@assays$ATAC@meta.features))
head(obj[['ATAC']]@var.features)

DefaultAssay(obj) <- 'ATAC'

obj@assays$ATAC@fragments

CoveragePlot(object = obj, region = "syt1a", assay = 'ATAC')
Idents(obj) <- 'neuron'

CoveragePlot(
  object = obj,
  region = "syt1a",      
  features = "syt1a",
  expression.assay = "RNA",
  extend.upstream = 100,
  extend.downstream = 100,
  assay = 'ATAC', 
  sep = c('-', '-'), 
  peaks = T,
  annotation = F)
