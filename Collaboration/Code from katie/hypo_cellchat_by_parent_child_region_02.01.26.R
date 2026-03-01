

group_list <- unique(tel$group)
group_list

gc()

#################################################
## MASTER LEVELS (match L9_Child_CellType)
#################################################

levels_L9 <- levels(factor(as.character(tel$L9_Child_CellType)))

for (grp in group_list) {
  
  #################################################
  ## L9_Child_CellType
  #################################################
  
  obj <- tel
  # obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay") 
  
  level <- 'L9_Child_CellType'
  Idents(obj) <- 'L9_Child_CellType'
  levels(obj)
  
  data.input <- GetAssayData(obj, assay = 'SCT', layer = 'data') # normalized data matrix
  
  data.input.tel <- as.matrix(data.input[homologs$mzebra_seurat_gene_id[which(!is.na(homologs$one_to_one_human))],])
  rownames(data.input.tel) <- homologs$one_to_one_human[which(!is.na(homologs$one_to_one_human))]
  head(rownames(data.input.tel))
  
  obj@assays$SCT@data <- data.input.tel
  head(rownames(obj@assays$SCT@data))
  
  gc()
  
  meta <- obj@meta.data
  cell.use = rownames(meta)[meta$group == grp] # extract the cell names from disease data
  
  # Prepare input data for CelChat analysis
  data.input.tel = data.input.tel[, cell.use]
  meta = meta[cell.use, ]
  unique(meta$L9_Child_CellType)
  
  meta$samples <- as.factor(meta$individual)
  
  #################################################
  ## IMPORTANT: DROP unused levels BEFORE computeCommunProb()
  ## (computeCommunProb fails if idents contain zero-cell levels)
  #################################################
  
  meta$L9_Child_CellType <- droplevels(factor(as.character(meta$L9_Child_CellType)))
  
  cellchat.obj <- createCellChat(object = data.input.tel, group.by = "L9_Child_CellType", assay = "SCT", meta = meta)
  
  cellchat.obj <- addMeta(cellchat.obj, meta = meta)
  cellchat.obj <- setIdent(cellchat.obj, ident.use = "L9_Child_CellType") # set default cell identity
  
  gc()
  
  # Database (human)
  CellChatDB <- CellChatDB.human
  cellchat.obj@DB <- CellChatDB
  
  # Subset signaling genes to reduce compute
  cellchat.obj <- subsetData(cellchat.obj)
  
  ############################################
  ## Preprocessing + overexpression
  ############################################
  
  # Conservative overexpression settings to reduce run/noise effects
  cellchat.obj_varFeats <- identifyOverExpressedGenes(
    cellchat.obj,
    thresh.p   = 0.05,
    min.cells  = 5,
    return.object = FALSE
  )
  length(unique(cellchat.obj_varFeats$features))
  
  outdir <- file.path("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/tel_ob/cellchat", grp, "child_clusters")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  write_csv(cellchat.obj_varFeats, file.path(outdir, paste0(level, "_cellchat.", grp, "_OverExpressedGenes.rds")))
  
  cellchat.obj <- identifyOverExpressedGenes(
    cellchat.obj,
    thresh.p   = 0.05,
    min.cells  = 5
  )
  
  cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
  
  ############################################
  ## Optional smoothing (CellChat-documented)
  ############################################
  
  data(PPI.human)
  PPI.human <- as(PPI.human, "sparseMatrix")
  cellchat.obj <- smoothData(cellchat.obj, adj = PPI.human)
  
  ############################################
  ## Compute communication probability
  ############################################
  
  cellchat.obj <- computeCommunProb(
    cellchat.obj,
    population.size = FALSE,
    raw.use = FALSE
  )
  
  # Aggregate to pathway level (more robust biologically)
  cellchat.obj <- computeCommunProbPathway(cellchat.obj)
  cellchat.obj <- aggregateNet(cellchat.obj)
  
  ############################################
  ## FILTERING INSIDE THE CELLCHAT OBJECT
  ############################################
  
  cellchat.obj_subset <- filterCommunication(
    cellchat.obj,
    min.cells     = 2,
    min.samples   = 8,
    rare.keep     = TRUE,
    nonFilter.keep = FALSE
  )
  
  #################################################
  ## PLOTTING FIX:
  ## Re-introduce missing master celltypes ONLY FOR PLOTTING
  ## by padding the net matrices + groupSize + color vector.
  #################################################
  
  pad_net_matrix <- function(mat, all_nodes) {
    out <- matrix(
      0,
      nrow = length(all_nodes),
      ncol = length(all_nodes),
      dimnames = list(all_nodes, all_nodes)
    )
    common <- intersect(rownames(mat), all_nodes)
    out[common, common] <- mat[common, common, drop = FALSE]
    out
  }
  
  idents_ordered <- c(
    "1.1.1-PV1/CAMK2G",
    "1.1.2-PV1/CAMK2G",
    "1.1.3-PV1/CACNB2",
    "1.2-PV1/SHOX2",
    "1.3-PV1/PCSK2",
    "2.1.1-LYPD1/FOXJ3",
    "2.1.2-LYPD1/RGS11",
    "2.2.1-LYPD1/PCSK2",
    "2.2.2-CCK",
    "2.3-PDYN/AgRP",
    "3.1.1-AgRP/NPY",
    "3.1.2-AgRP/NPY",
    "3.1.3-NPY",
    "3.2.1-SF1/BDNF",
    "3.2.2-PNOC/NPY",
    "3.3.1-GABA/TRH",
    "3.3.2-AgRP/NPY",
    "4-ZEB2/ST18",
    "5.1-CHRM4/FOXP4",
    "5.2-CHRM4/DGKZ",
    "6.1-FOXP4/PROX1",
    "6.2-FOXP4/PROX1",
    "6.3-FOXP4/FOSB",
    "7.1-ZFHX4/KCNH4",
    "7.2.1-ZFHX4/FOSB",
    "7.2.2-ZFHX4/ERBB4",
    "8.1-LMX1A",
    "8.2-LMX1A/vglut2",
    "8.3-LMX1A/vgat/vglut",
    "8.4-LMX1A/nos1",
    "9.1-Oligo/APOE",
    "9.1.1-Oligo/NFASC",
    "9.1.2-MF-Oligo-MPZ",
    "9.2-Oligo/CAMK2G",
    "9.3-Oligo-pnOls-PTPRD",
    "10.1-PNN",
    "10.2-PNN",
    "11.1.1-IRX1/CALN1",
    "11.1.2-IRX1/CALN1",
    "11.2-IRX1/PLA2G4A",
    "12.1-Magno-OXT",
    "12.2.1-Magno-AVP",
    "12.2.2-Magno-AVP",
    "13.1.1-Astro-Proto",
    "13.1.2-Astro-Trans",
    "13.2.1-a-Tany",
    "13.2.2-b-Tany",
    "13.3-Pituicyte-Pit",
    "14.1-MG",
    "14.2-MG",
    "15.1-CALB2/TSHZ2",
    "15.2-NXPH1/BCL11B",
    "16.1-SF1/PACAP",
    "16.2-NPY/CBLN1",
    "17.1.1-PACAP/ZNF536",
    "17.1.2-PACAP/ZNF536",
    "17.2-PRDM8/EBF1",
    "18-Astro-Fibro",
    "19.1-SATB1/PROX1",
    "19.2-ESR1",
    "19.3-LHX6/ARX",
    "20.1.1-GHRH",
    "20.1.2-FOSB",
    "20.2-KNDy",
    "20.3-FOXP1/ISL1",
    "20.4-LepR",
    "21.1.1-RORA/LEF1",
    "21.1.2-QKI/FAT2",
    "21.1.3-RORA/LEF1",
    "21.2-QKI/KCNC3",
    "21.3-GATA2",
    "22-TIDA",
    "23.1-Gonado-Pit",
    "23.2-Folliculo-Pit",
    "24-Thermo",
    "25-OPC",
    "26-FOXP4/RP3V",
    "27.1-Somato-Pit",
    "27.2.1-ptThyro-Pit",
    "27.2.2-pdThyro-Pit",
    "28-Lacto-Pit",
    "29.1-SHOX2/TCF7L2",
    "29.2-SHOX2/TCF7L2",
    "30-q-hmRG",
    "31.1-CP",
    "31.2-BMEC",
    "31.3-PVF",
    "32-ParvoPre-TRH",
    "33-Cortico-Pit",
    "34-CARTPT/CRHR1",
    "35-MEIS2",
    "36.1-Pit-SOX2/PROP1",
    "36.2-Melano-Pit",
    "37-TH",
    "38.1-NFOL",
    "38.2-NFOL",
    "39-Chol"
  )
  
  #################################################
  ## FIX: prevent factor coercion from turning unknown labels into NA
  #################################################
  
  idents_present <- sort(unique(as.character(cellchat.obj_subset@idents)))
  missing_in_order <- setdiff(idents_present, idents_ordered)
  if (length(missing_in_order) > 0) {
    message("Appending missing cell types not in idents_ordered: ", paste(missing_in_order, collapse = ", "))
    idents_ordered_use <- c(idents_ordered, missing_in_order)
  } else {
    idents_ordered_use <- idents_ordered
  }
  
  # Apply to CellChat object (safe ordering)
  cellchat.obj@idents <- factor(cellchat.obj@idents, levels = idents_ordered_use)
  
  all_nodes <- idents_ordered_use
  
  color_key_child <- c(
    "1.1.1-PV1/CAMK2G"            = "#F2BEE1",
    "1.1.2-PV1/CAMK2G"            = "#E993CB",
    "1.1.3-PV1/CACNB2"            = "#ff9acd",
    "1.2-PV1/SHOX2"               = "#E275BD",
    "1.3-PV1/PCSK2"               = "#ffbfe2",
    "2.1.1-LYPD1/FOXJ3"           = "#DB1F83",
    "2.1.2-LYPD1/RGS11"           = "#E6589F",
    "2.2.1-LYPD1/PCSK2"           = "#EE88BA",
    "2.2.2-CCK"                   = "#F7BFDA",
    "2.3-PDYN/AgRP"               = "#F6A5C5",
    "3.1.1-AgRP/NPY"              = "#BBC7FF",
    "3.1.2-AgRP/NPY"              = "#91A7FF",
    "3.1.3-NPY"                   = "#5D7CFA",
    "3.2.1-SF1/BDNF"              = "#4C6EF5",
    "3.2.2-PNOC/NPY"              = "#364FC7",
    "3.3.1-GABA/TRH"              = "#4087FF",
    "3.3.2-AgRP/NPY"              = "#6CA3FF",
    "4-ZEB2/ST18"                 = "#BA2D22",
    "5.1-CHRM4/FOXP4"             = "#EDA0E4",
    "5.2-CHRM4/DGKZ"              = "#EC83E0",
    "6.1-FOXP4/PROX1"             = "#0d8e84",
    "6.2-FOXP4/PROX1"             = "#66c5bc",
    "6.3-FOXP4/FOSB"              = "#b2dfdb",
    "7.1-ZFHX4/KCNH4"             = "#DA81A8",
    "7.2.1-ZFHX4/FOSB"            = "#C94683",
    "7.2.2-ZFHX4/ERBB4"           = "#B80464",
    "8.1-LMX1A"                   = "#FABFB7",
    "8.2-LMX1A/vglut2"            = "salmon4",
    "8.3-LMX1A/vgat/vglut"        = "salmon1",
    "8.4-LMX1A/nos1"              = "salmon3",
    "9.1.1-Oligo/NFASC"           = "#87ab69",
    "9.1.2-MF-Oligo-MPZ"          = "#a3c585",
    "9.1-Oligo/APOE"              = "#4b6043",
    "9.2-Oligo/CAMK2G"            = "#658354",
    "9.3-Oligo-pnOls-PTPRD"       = "#c7ddb5",
    "10.1-PNN"                    = "#ED823E",
    "10.2-PNN"                    = "#F19C7C",
    "11.1.1-IRX1/CALN1"           = "#8b6c5c",
    "11.1.2-IRX1/CALN1"           = "#a08679",
    "11.2-IRX1/PLA2G4A"           = "#bca89f",
    "12.1-Magno-OXT"              = "#F51B1B",
    "12.2.1-Magno-AVP"            = "#F76060",
    "12.2.2-Magno-AVP"            = "#F76060",
    "13.1.1-Astro-Proto"          = "lightgoldenrod3",
    "13.1.2-Astro-Trans"          = "lightgoldenrod1",
    "13.2.1-a-Tany"               = "goldenrod3",
    "13.2.2-b-Tany"               = "#F2D67F",
    "13.3-Pituicyte-Pit"          = "#F2D61F",
    "14.1-MG"                     = "#20b3c4",
    "14.2-MG"                     = "cadetblue2",
    "15.1-CALB2/TSHZ2"            = "#E4928E",
    "15.2-NXPH1/BCL11B"           = "#FF9F9F",
    "16.1-SF1/PACAP"              = "#4444C4",
    "16.2-NPY/CBLN1"              = "#000080",
    "17.1.1-PACAP/ZNF536"         = "#B3544D",
    "17.1.2-PACAP/ZNF536"         = "#F7574A",
    "17.2-PRDM8/EBF1"             = "#E37E76",
    "18-Astro-Fibro"              = "#F2D67F",
    "19.1-SATB1/PROX1"            = "#ade5e2",
    "19.2-ESR1"                   = "#1763AB",
    "19.3-LHX6/ARX"               = "#1C7ED5",
    "20.1.1-GHRH"                 = "#A4D8FF",
    "20.1.2-FOSB"                 = "#339AF0",
    "20.2-KNDy"                   = "#1871C2",
    "20.3-FOXP1/ISL1"             = "#238BE6",
    "20.4-LepR"                   = "#74C0FC",
    "21.1.1-RORA/LEF1"            = "#C3DAFF",
    "21.1.2-QKI/FAT2"             = "#92b4fe",
    "21.1.3-RORA/LEF1"            = "#79B8F4",
    "21.2-QKI/KCNC3"              = "#A7C8FF",
    "21.3-GATA2"                  = "#729efd",
    "22-TIDA"                     = "steelblue3",
    "23.1-Gonado-Pit"             = "#9f87ca",
    "23.2-Folliculo-Pit"          = "#d9b9e2",
    "24-Thermo"                   = "lightsteelblue3",
    "25-OPC"                      = "#aedcae",
    "26-FOXP4/RP3V"               = "indianred",
    "27.1-Somato-Pit"             = "#8953a9",
    "27.2.1-ptThyro-Pit"          = "#bba8d9",
    "27.2.2-pdThyro-Pit"          = "#734f96",
    "28-Lacto-Pit"                = "#bb8fce",
    "29.1-SHOX2/TCF7L2"           = "#EB4F16",
    "29.2-SHOX2/TCF7L2"           = "#EB6E40",
    "30-q-hmRG"                   = "goldenrod1",
    "31.1-CP"                     = "#DD8C00",
    "31.2-BMEC"                   = "#EB9F62",
    "31.3-PVF"                    = "#D49200",
    "32-ParvoPre-TRH"             = "#E66055",
    "33-Cortico-Pit"              = "#7c6192",
    "34-CARTPT/CRHR1"             = "#D631C7",
    "35-MEIS2"                    = "#7865EF",
    "36.1-Pit-SOX2/PROP1"         = "#854abf",
    "36.2-Melano-Pit"             = "#c098f2",
    "37-TH"                       = "#4125F0",
    "38.1-NFOL"                   = "#ACCA82",
    "38.2-NFOL"                   = "#75975e",
    "39-Chol"                     = "#816FF0"
  )
  
  color_key_child_df <- data.frame(
    celltype = names(color_key_child),
    color    = unname(color_key_child),
    stringsAsFactors = FALSE
  )
  
  pdf(file.path(outdir, paste0(grp, "_netCellCellInteractions_", level, '.pdf')), width = 16, height = 18)
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  color.use.plot <- color_key_child[all_nodes]
  
  net.count.plot  <- pad_net_matrix(cellchat.obj_subset@net$count,  all_nodes)
  net.weight.plot <- pad_net_matrix(cellchat.obj_subset@net$weight, all_nodes)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.count.plot, 
                       vertex.weight = groupSize, vertex.size.max	= 5,
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Number of Interactions: Top 10%", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.1,
                       angle_text = TRUE,
                       remove.isolate = F,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.weight.plot, 
                       vertex.weight = groupSize, vertex.size.max	= 5, 
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Interaction Weights/Strength: Top 10%", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.1,
                       angle_text = TRUE,
                       remove.isolate = F,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  pit_celltypes <- c(
    "13.3-Pituicyte-Pit",
    "23.1-Gonado-Pit",
    "23.2-Folliculo-Pit",
    "27.1-Somato-Pit",
    "27.2.1-ptThyro-Pit",
    "27.2.2-pdThyro-Pit",
    "28-Lacto-Pit",
    "33-Cortico-Pit",
    "36.1-Pit-SOX2/PROP1",
    "36.2-Melano-Pit"
  )
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.weight.plot, 
                       vertex.weight = groupSize, 
                       vertex.size.max	= 5,
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Targeting Pituitary CellTypes (interaction strength - Top 50%)", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.5, angle_text = TRUE,
                       targets.use = pit_celltypes,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.weight.plot, 
                       vertex.weight = groupSize, 
                       weight.scale = T, 
                       vertex.size.max	= 5,
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Originating from Pituitary CellTypes (interaction strength - Top 50%", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.5, angle_text = TRUE,
                       sources.use = pit_celltypes,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  dev.off()
  
  pdf(file.path(outdir, paste0(grp, "_netInteractionWeights_per_", level, '.pdf')), width = 25, height = 28)
  mat <- net.weight.plot
  par(mfrow = c(2,2), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle_mod(mat2, 
                         vertex.weight = groupSize, 
                         vertex.size.max	= 5, 
                         weight.scale = T, 
                         edge.weight.max = max(mat), 
                         angle_text = TRUE,
                         color.use = color.use.plot,
                         title.name = paste0(rownames(mat)[i], '\n\n\n\n'))
  }
  dev.off()
  
  saveRDS(
    cellchat.obj,
    file.path(outdir, paste0(level, "_cellchat_tel_", grp, "_.rds"))
  )
  
  ############################################
  ## Export robust LR-level table for ggplot
  ############################################
  
  df.net <- subsetCommunication(cellchat.obj, thresh = 0.05)
  write_csv(df.net, file.path(outdir, paste0(level, "_cellchat.", grp, "_netCellCommunications.csv")))
  
  df.netP <- subsetCommunication(cellchat.obj, thresh = 0.05, slot.name = 'netP')
  write_csv(df.netP, file.path(outdir, paste0(level, "_cellchat.", grp, "_netSignalingPathways.csv")))
  
  gc()
  
  library(pheatmap)
  
  # apply safe ordering (prevents NA)
  df.net$source <- factor(as.character(df.net$source), levels = all_nodes)
  df.net$target <- factor(as.character(df.net$target), levels = all_nodes)
  
  # Quick checks
  levels(cellchat.obj@idents)
  all(levels(cellchat.obj@idents) == all_nodes)
  
  ################################
  ## 4. RIDGE PLOTS (FIXED)
  ################################
  
  library(ggridges)
  
  library(ggplot2)
  library(ggridges)
  
  # enforce your desired celltype order ONCE
  df.net$source <- factor(as.character(df.net$source), levels = all_nodes)
  df.net$target <- factor(as.character(df.net$target), levels = all_nodes)
  
  df.net$source_color <- as.character(df.net$source)
  df.net$source_color <- color_key_child_df$color[match(df.net$source, color_key_child_df$celltype)]
  
  df.net$target_color <- as.character(df.net$target)
  df.net$target_color <- color_key_child_df$color[match(df.net$target, color_key_child_df$celltype)]
  
  pdf(file.path(outdir,paste0(grp, "_netComminicationDistribution_", level, '.pdf')), width = 10, height = 20)
  
  ggplot(
    df.net,
    aes(
      x = prob,
      y = source,
      fill = source_color
    )
  ) +
    geom_density_ridges(
      scale = 2.4,
      alpha = 0.85,
      color = "gray",
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    scale_x_log10() +
    scale_y_discrete(
      limits = rev(all_nodes),
      drop = FALSE
    ) +
    scale_fill_identity(guide = "none") +
    theme_ridges(grid = FALSE) +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      text = element_text(size = 13),
      axis.text = element_text(size = 11),
      axis.title.x = element_text(size = 13, hjust = 0.5),
      axis.title.y = element_text(size = 13, hjust = 0.5),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = paste0(grp, ' - ', level, ": CellChat communication probability distributions by sender"),
      x = "Communication probability (log10 scale)",
      y = "Sender cell type"
    )
  
  ggplot(
    df.net,
    aes(
      x = prob,
      y = target,
      fill = target_color
    )
  ) +
    geom_density_ridges(
      scale = 2.4,
      alpha = 0.85,
      color = "gray",
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    scale_x_log10() +
    scale_y_discrete(
      limits = rev(all_nodes),
      drop = FALSE
    ) +
    scale_fill_identity(guide = "none") +
    theme_ridges(grid = FALSE) +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      text = element_text(size = 13),
      axis.text = element_text(size = 11),
      axis.title.x = element_text(size = 13, hjust = 0.5),
      axis.title.y = element_text(size = 13, hjust = 0.5),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = paste0(grp, ' - ', level, ": CellChat communication probability distributions by receiver"),
      x = "Communication probability (log10 scale)",
      y = "Receiver cell type"
    )
  
  ggplot(
    df.net,
    aes(
      x = prob,
      y = reorder(pathway_name, prob, FUN = median),
      fill = pathway_name
    )
  ) +
    geom_density_ridges(
      scale = 2.6,
      alpha = 0.85,
      color = "gray",
      linewidth = 0.3
    ) +
    scale_x_log10() +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    theme_ridges(grid = FALSE) +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 11),
      axis.title.x = element_text(size = 14, hjust = 0.5),
      axis.title.y = element_text(size = 14, hjust = 0.5),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5)
    ) +
    labs(
      title = paste0(grp, ' - ', level, ": Distribution of CellChat communication probabilities by pathway"),
      subtitle = paste0(grp, ' - ', level, ": Pathways ordered by median probability"),
      x = "Communication probability (log10 scale)",
      y = "Signaling pathway"
    )
  
  df.net %>%
    group_by(source, source_color) %>%
    summarise(total_outgoing = sum(prob), .groups = "drop") %>%
    ggplot(
      aes(
        x = reorder(source, total_outgoing),
        y = total_outgoing,
        fill = source_color
      )
    ) +
    geom_col() +
    scale_fill_identity() +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.x  = element_text(color = "black"),
      axis.text.y  = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title   = element_text(color = "black")
    ) +
    labs(
      title = "Total outgoing cell–cell signaling by source cell type",
      x = "Source (signaling) cell type",
      y = "Summed communication probability"
    )
  
  df.net %>%
    group_by(target, target_color) %>%
    summarise(total_incoming = sum(prob), .groups = "drop") %>%
    ggplot(
      aes(
        x = reorder(target, total_incoming),
        y = total_incoming,
        fill = target_color
      )
    ) +
    geom_col() +
    scale_fill_identity() +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.x  = element_text(color = "black"),
      axis.text.y  = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title   = element_text(color = "black")
    ) +
    labs(
      title = "Total incoming cell–cell signaling by target cell type",
      x = "Target (receiving) cell type",
      y = "Summed communication probability"
    )
  
  dev.off()
  
  pdf(file.path(outdir,paste0(grp, "_netComminicationStrength_heatmap_", level, '.pdf')), width = 18, height = 15)
  mat <- df.net %>%
    group_by(source, target) %>%
    summarise(weight = sum(prob), .groups = "drop") %>%
    pivot_wider(names_from = target, values_from = weight, values_fill = 0) %>%
    column_to_rownames("source") %>%
    as.matrix()
  
  pheatmap(
    mat,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#F7FBFF", 'grey90', 'lightblue1', "lightblue", "#41B6C4", 'deepskyblue3', "#225EA8", "#253494"))(100),
    main = "Cell–cell communication strength (summed probabilities)\nRows = source (senders), columns = target (receivers)",
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45
  )
  dev.off()
  
  library(igraph)
  library(ggraph)
  
  pdf(file.path(outdir,paste0(grp, "_LigandReceptor_interactions groupedBy_signaling.pathway", level, '.pdf')), width = 35, height = 90)
  
  colors <- c("lightblue", "#41B6C4", 'deepskyblue3', "#225EA8", "#253494")
  
  df.net %>%
    group_by(interaction_name, pathway_name) %>%
    summarise(
      mean_prob = mean(prob),
      n_pairs = n(),
      .groups = "drop"
    ) %>%
    ggplot(aes(
      x = pathway_name,
      y = interaction_name,
      size = mean_prob,
      color = n_pairs
    )) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(trans = "log10") +
    scale_color_gradientn(
      colors = colors
    ) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_line(color = "gray"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.text.y   = element_text(size = 15, color = "black"),
      axis.text.x   = element_text(size = 15,color = "black"),
      axis.title.x  = element_text(size = 15,color = "black"),
      axis.title.y  = element_text(size = 15,color = "black"),
      legend.text   = element_text(size = 15,color = "black"),
      legend.title  = element_text(size = 15,color = "black"),
      plot.title    = element_text(size = 18, color = "black")
    ) +
    labs(
      title = "Ligand–receptor interactions grouped by signaling pathway",
      x = "Signaling pathway",
      y = "Ligand–receptor pair",
      size = "Mean communication probability",
      color = "# interactions"
    ) +
    RotatedAxis()
  
  dev.off()
  
  saveRDS(
    cellchat.obj,
    file.path(outdir, paste0("cellchat_tel_", grp, "_L9_Child_CellType_netCentrality.rds"))
  )
  
  # add to environment
  assign(paste0("cellchat_tel_", grp, "_L9_Child_CellType"), cellchat.obj, envir = .GlobalEnv)
  
}








levels_L8 <- levels(factor(as.character(tel$L8_Parent_CellType)))

for (grp in group_list) {
  obj <- tel
  # obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay") 
  
  level <- 'L8_Parent_CellType'
  Idents(obj) <- 'L8_Parent_CellType'
  levels(obj)
  
  data.input <- GetAssayData(obj, assay = 'SCT', layer = 'data') # normalized data matrix
  
  data.input.tel <- as.matrix(data.input[homologs$mzebra_seurat_gene_id[which(!is.na(homologs$one_to_one_human))],])
  rownames(data.input.tel) <- homologs$one_to_one_human[which(!is.na(homologs$one_to_one_human))]
  head(rownames(data.input.tel))
  
  obj@assays$SCT@data <- data.input.tel
  head(rownames(obj@assays$SCT@data))
  
  gc()
  
  meta <- obj@meta.data
  cell.use = rownames(meta)[meta$group == grp] # extract the cell names from disease data
  
  # Prepare input data for CelChat analysis
  data.input.tel = data.input.tel[, cell.use]
  meta = meta[cell.use, ]
  unique(meta$L8_Parent_CellType)
  
  meta$samples <- as.factor(meta$individual)
  
  #################################################
  ## IMPORTANT: DROP unused levels BEFORE computeCommunProb()
  ## (computeCommunProb fails if idents contain zero-cell levels)
  #################################################
  
  meta$L8_Parent_CellType <- droplevels(factor(as.character(meta$L8_Parent_CellType)))
  
  cellchat.obj <- createCellChat(object = data.input.tel, group.by = "L8_Parent_CellType", assay = "SCT", meta = meta)
  
  cellchat.obj <- addMeta(cellchat.obj, meta = meta)
  cellchat.obj <- setIdent(cellchat.obj, ident.use = "L8_Parent_CellType") # set default cell identity
  
  gc()
  
  # Database (human)
  CellChatDB <- CellChatDB.human
  cellchat.obj@DB <- CellChatDB
  
  # Subset signaling genes to reduce compute
  cellchat.obj <- subsetData(cellchat.obj)
  
  ############################################
  ## Preprocessing + overexpression
  ############################################
  
  # Conservative overexpression settings to reduce run/noise effects
  cellchat.obj_varFeats <- identifyOverExpressedGenes(
    cellchat.obj,
    thresh.p   = 0.05,
    min.cells  = 5,
    return.object = FALSE
  )
  length(unique(cellchat.obj_varFeats$features))
  
  outdir <- file.path("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/tel_ob/cellchat", grp, "parent_clusters")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  write_csv(cellchat.obj_varFeats, file.path(outdir, paste0(level, "_cellchat.", grp, "_OverExpressedGenes.rds")))
  
  cellchat.obj <- identifyOverExpressedGenes(
    cellchat.obj,
    thresh.p   = 0.05,
    min.cells  = 5
  )
  
  cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
  
  ############################################
  ## Optional smoothing (CellChat-documented)
  ############################################
  
  data(PPI.human)
  PPI.human <- as(PPI.human, "sparseMatrix")
  cellchat.obj <- smoothData(cellchat.obj, adj = PPI.human)
  
  ############################################
  ## Compute communication probability
  ############################################
  
  cellchat.obj <- computeCommunProb(
    cellchat.obj,
    population.size = FALSE,
    raw.use = FALSE
  )
  
  # Aggregate to pathway level (more robust biologically)
  cellchat.obj <- computeCommunProbPathway(cellchat.obj)
  cellchat.obj <- aggregateNet(cellchat.obj)
  
  ############################################
  ## FILTERING INSIDE THE CELLCHAT OBJECT
  ############################################
  
  cellchat.obj_subset <- filterCommunication(
    cellchat.obj,
    min.cells     = 2,
    min.samples   = 8,
    rare.keep     = TRUE,
    nonFilter.keep = FALSE
  )
  
  #################################################
  ## PLOTTING FIX:
  ## Re-introduce missing master celltypes ONLY FOR PLOTTING
  ## by padding the net matrices + groupSize + color vector.
  #################################################
  
  pad_net_matrix <- function(mat, all_nodes) {
    out <- matrix(
      0,
      nrow = length(all_nodes),
      ncol = length(all_nodes),
      dimnames = list(all_nodes, all_nodes)
    )
    common <- intersect(rownames(mat), all_nodes)
    out[common, common] <- mat[common, common, drop = FALSE]
    out
  }
  
  idents_ordered <- c(
    "1-PV1",
    "2-LYPD1",
    "3-NPY",
    "4-EBF1",
    "5-CHRM4",
    "6-GRIK2",
    "7-ZFHX4",
    "8-LMX1A",
    "9-Oligo",
    "10-PNN",
    "11-LHX9",
    "12-Magno",
    "13-Astroglia",
    "14-MG",
    "15-PITX2",
    "16-SOX5",
    "17-ZEB2",
    "18-Astro",
    "19-VIP",
    "20-FOXP2",
    "21-KCNC3",
    "22-TIDA",
    "23-Pit-PITX1",
    "24-Thermo",
    "25-OPC",
    "26-FOXP4",
    "27-Pit-POU1F1",
    "28-Pit-PRL",
    "29-SHOX2",
    "30-RG",
    "31-Vascular",
    "32-ParvoPre",
    "33-Pit-Cortico",
    "34-CARTPT",
    "35-MEIS2",
    "36-Pit-SOX2",
    "37-TH",
    "38-NFOL",
    "39-Chol"
  )
  
  # Apply to CellChat object
  cellchat.obj@idents <- factor(cellchat.obj@idents, levels = idents_ordered)
  
  all_nodes <-idents_ordered
  
  color_key_parent <- c(
    "1-PV1"            = "#E993CB",
    "2-LYPD1"          = "#EE88BA",
    "3-NPY"            = "#BBC7FF",
    "4-EBF1"           = "#BA2D22",
    "5-CHRM4"          = "#EDA0E4",
    "6-GRIK2"          = "#66c5bc",
    "7-ZFHX4"          = "#DA81A8",
    "8-LMX1A"          = "#FABFB7",
    "9-Oligo"          = "#87ab69",
    "10-PNN"           = "#F19C7C",
    "11-LHX9"          = "#8b6c5c",
    "12-Magno"         = "#F51B1B",
    "13-Astroglia"     = "goldenrod3",
    "14-MG"            = "#20b3c4",
    "15-PITX2"         = "#E4928E",
    "16-SOX5"          = "#4444C4",
    "17-ZEB2"          = "#B3544D",
    "18-Astro"         = "#F2D67F",
    "19-VIP"           = "#1763AB",
    "20-FOXP2"         = "#A4D8FF",
    "21-KCNC3"         = "#729efd",
    "22-TIDA"          = "steelblue3",
    "23-Pit-PITX1"     = "#9f87ca",
    "24-Thermo"        = "lightsteelblue3",
    "25-OPC"           = "#658354",
    "26-FOXP4"         = "indianred",
    "27-Pit-POU1F1"    = "#734f96",
    "28-Pit-PRL"       = "#bb8fce",
    "29-SHOX2"         = "#EB6E40",
    "30-RG"            = "goldenrod1",
    "31-Vascular"      = "#DD8C00",
    "32-ParvoPre"      = "#E66055",
    "33-Pit-Cortico"   = "#7c6192",
    "34-CARTPT"        = "#D631C7",
    "35-MEIS2"         = "#7865EF",
    "36-Pit-SOX2"      = "#854abf",
    "37-TH"            = "#4125F0",
    "38-NFOL"          = "#ACCA82",
    "39-Chol"          = "#816FF0"
  )
  
  color_key_parent_df <- data.frame(
    celltype = names(color_key_parent),
    color    = unname(color_key_parent),
    stringsAsFactors = FALSE
  )
  
  pdf(file.path(outdir, paste0(grp, "_netCellCellInteractions_", level, '.pdf')), width = 16, height = 18)
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  color.use.plot <- color_key_parent[all_nodes]
  
  net.count.plot  <- pad_net_matrix(cellchat.obj_subset@net$count,  all_nodes)
  net.weight.plot <- pad_net_matrix(cellchat.obj_subset@net$weight, all_nodes)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.count.plot, 
                       vertex.weight = groupSize, vertex.size.max	= 5,
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Number of Interactions: Top 10%", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.1,
                       angle_text = TRUE,
                       remove.isolate = F,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.weight.plot, 
                       vertex.weight = groupSize, vertex.size.max	= 5, 
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Interaction Weights/Strength: Top 10%", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.1,
                       angle_text = TRUE,
                       remove.isolate = F,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  pit_celltypes <- c(
    "13-Astroglia",
    "23-Pit-PITX1",
    "27-Pit-POU1F1",
    "28-Pit-PRL",
    "33-Pit-Cortico",
    "36-Pit-SOX2"
  )
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.weight.plot, 
                       vertex.weight = groupSize, 
                       vertex.size.max	= 5,
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Targeting Pituitary CellTypes (interaction strength - Top 50%)", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.5, angle_text = TRUE,
                       targets.use = pit_celltypes,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.weight.plot, 
                       vertex.weight = groupSize, 
                       weight.scale = T, 
                       vertex.size.max	= 5,
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Originating from Pituitary CellTypes (interaction strength - Top 50%", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.5, angle_text = TRUE,
                       sources.use = pit_celltypes,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  dev.off()
  
  pdf(file.path(outdir, paste0(grp, "_netInteractionWeights_per_", level, '.pdf')), width = 25, height = 28)
  mat <- net.weight.plot
  par(mfrow = c(2,2), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle_mod(mat2, 
                         vertex.weight = groupSize, 
                         vertex.size.max	= 5, 
                         weight.scale = T, 
                         edge.weight.max = max(mat), 
                         angle_text = TRUE,
                         color.use = color.use.plot,
                         title.name = paste0(rownames(mat)[i], '\n\n\n\n'))
  }
  dev.off()
  
  saveRDS(
    cellchat.obj,
    file.path(outdir, paste0(level, "_cellchat_tel_", grp, "_.rds"))
  )
  
  ############################################
  ## Compute centrality on FINAL filtered network
  ############################################
  
  # cellchat.obj <- netAnalysis_computeCentrality(
  #   cellchat.obj,
  #   slot.name = "netP"
  # )
  
  ############################################
  ## Export robust LR-level table for ggplot
  ############################################
  
  df.net <- subsetCommunication(cellchat.obj, thresh = 0.05)
  write_csv(df.net, file.path(outdir, paste0(level, "_cellchat.", grp, "_netCellCommunications.csv")))
  
  df.netP <- subsetCommunication(cellchat.obj, thresh = 0.05, slot.name = 'netP')
  write_csv(df.netP, file.path(outdir, paste0(level, "_cellchat.", grp, "_netSignalingPathways.csv")))
  
  gc()
  
  library(pheatmap)
  
  # (Recommended) apply to df.net for consistent ggplot facet/axis ordering
  df.net$source <- factor(df.net$source, levels = idents_ordered)
  df.net$target <- factor(df.net$target, levels = idents_ordered)
  
  # Quick checks
  levels(cellchat.obj@idents)
  all(levels(cellchat.obj@idents) == idents_ordered)
  
  ################################
  ## 4. RIDGE PLOTS (FIXED)
  ################################
  
  library(ggridges)
  
  library(ggplot2)
  library(ggridges)
  
  # enforce your desired celltype order ONCE
  df.net$source <- factor(df.net$source, levels = idents_ordered)
  df.net$target <- factor(df.net$target, levels = idents_ordered)
  
  df.net$source_color <- as.character(df.net$source)
  df.net$source_color <- color_key_parent_df$color[match(df.net$source, color_key_parent_df$celltype)]
  
  df.net$target_color <- as.character(df.net$target)
  df.net$target_color <- color_key_parent_df$color[match(df.net$target, color_key_parent_df$celltype)]
  
  pdf(file.path(outdir,paste0(grp, "_netComminicationDistribution_", level, '.pdf')), width = 10, height = 20)
  
  ggplot(
    df.net,
    aes(
      x = prob,
      y = source,
      fill = source_color
    )
  ) +
    geom_density_ridges(
      scale = 2.4,
      alpha = 0.85,
      color = "gray",
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    scale_x_log10() +
    scale_y_discrete(
      limits = rev(idents_ordered),
      drop = FALSE
    ) +
    scale_fill_identity(guide = "none") +
    theme_ridges(grid = FALSE) +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      text = element_text(size = 13),
      axis.text = element_text(size = 11),
      axis.title.x = element_text(size = 13, hjust = 0.5),
      axis.title.y = element_text(size = 13, hjust = 0.5),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = paste0(grp, ' - ', level, ": CellChat communication probability distributions by sender"),
      x = "Communication probability (log10 scale)",
      y = "Sender cell type"
    )
  
  ggplot(
    df.net,
    aes(
      x = prob,
      y = target,
      fill = target_color
    )
  ) +
    geom_density_ridges(
      scale = 2.4,
      alpha = 0.85,
      color = "gray",
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    scale_x_log10() +
    scale_y_discrete(
      limits = rev(idents_ordered),
      drop = FALSE
    ) +
    scale_fill_identity(guide = "none") +
    theme_ridges(grid = FALSE) +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      text = element_text(size = 13),
      axis.text = element_text(size = 11),
      axis.title.x = element_text(size = 13, hjust = 0.5),
      axis.title.y = element_text(size = 13, hjust = 0.5),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = paste0(grp, ' - ', level, ": CellChat communication probability distributions by receiver"),
      x = "Communication probability (log10 scale)",
      y = "Receiver cell type"
    )
  
  ggplot(
    df.net,
    aes(
      x = prob,
      y = reorder(pathway_name, prob, FUN = median),
      fill = pathway_name
    )
  ) +
    geom_density_ridges(
      scale = 2.6,
      alpha = 0.85,
      color = "gray",
      linewidth = 0.3
    ) +
    scale_x_log10() +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    theme_ridges(grid = FALSE) +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 11),
      axis.title.x = element_text(size = 14, hjust = 0.5),
      axis.title.y = element_text(size = 14, hjust = 0.5),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5)
    ) +
    labs(
      title = paste0(grp, ' - ', level, ": Distribution of CellChat communication probabilities by pathway"),
      subtitle = paste0(grp, ' - ', level, ": Pathways ordered by median probability"),
      x = "Communication probability (log10 scale)",
      y = "Signaling pathway"
    )
  
  df.net %>%
    group_by(source, source_color) %>%
    summarise(total_outgoing = sum(prob), .groups = "drop") %>%
    ggplot(
      aes(
        x = reorder(source, total_outgoing),
        y = total_outgoing,
        fill = source_color
      )
    ) +
    geom_col() +
    scale_fill_identity() +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.x  = element_text(color = "black"),
      axis.text.y  = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title   = element_text(color = "black")
    ) +
    labs(
      title = "Total outgoing cell–cell signaling by source cell type",
      x = "Source (signaling) cell type",
      y = "Summed communication probability"
    )
  
  df.net %>%
    group_by(target, target_color) %>%
    summarise(total_incoming = sum(prob), .groups = "drop") %>%
    ggplot(
      aes(
        x = reorder(target, total_incoming),
        y = total_incoming,
        fill = target_color
      )
    ) +
    geom_col() +
    scale_fill_identity() +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.x  = element_text(color = "black"),
      axis.text.y  = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title   = element_text(color = "black")
    ) +
    labs(
      title = "Total incoming cell–cell signaling by target cell type",
      x = "Target (receiving) cell type",
      y = "Summed communication probability"
    )
  
  dev.off()
  
  pdf(file.path(outdir,paste0(grp, "_netComminicationStrength_heatmap_", level, '.pdf')), width = 18, height = 15)
  mat <- df.net %>%
    group_by(source, target) %>%
    summarise(weight = sum(prob), .groups = "drop") %>%
    pivot_wider(names_from = target, values_from = weight, values_fill = 0) %>%
    column_to_rownames("source") %>%
    as.matrix()
  
  pheatmap(
    mat,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#F7FBFF", 'grey90', 'lightblue1', "lightblue", "#41B6C4", 'deepskyblue3', "#225EA8", "#253494"))(100),
    main = "Cell–cell communication strength (summed probabilities)\nRows = source (senders), columns = target (receivers)",
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45
  )
  dev.off()
  
  library(igraph)
  library(ggraph)
  
  pdf(file.path(outdir,paste0(grp, "_LigandReceptor_interactions groupedBy_signaling.pathway", level, '.pdf')), width = 35, height = 90)
  
  colors <- c("lightblue", "#41B6C4", 'deepskyblue3', "#225EA8", "#253494")
  
  df.net %>%
    group_by(interaction_name, pathway_name) %>%
    summarise(
      mean_prob = mean(prob),
      n_pairs = n(),
      .groups = "drop"
    ) %>%
    ggplot(aes(
      x = pathway_name,
      y = interaction_name,
      size = mean_prob,
      color = n_pairs
    )) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(trans = "log10") +
    scale_color_gradientn(
      colors = colors
    ) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_line(color = "gray"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.text.y   = element_text(size = 15, color = "black"),
      axis.text.x   = element_text(size = 15,color = "black"),
      axis.title.x  = element_text(size = 15,color = "black"),
      axis.title.y  = element_text(size = 15,color = "black"),
      legend.text   = element_text(size = 15,color = "black"),
      legend.title  = element_text(size = 15,color = "black"),
      plot.title    = element_text(size = 18, color = "black")
    ) +
    labs(
      title = "Ligand–receptor interactions grouped by signaling pathway",
      x = "Signaling pathway",
      y = "Ligand–receptor pair",
      size = "Mean communication probability",
      color = "# interactions"
    ) +
    RotatedAxis()
  
  dev.off()
  
  saveRDS(
    cellchat.obj,
    file.path(outdir, paste0("cellchat_tel_", grp, "_L8_Parent_CellType_netCentrality.rds"))
  )
  
  # add to environment
  assign(paste0("cellchat_tel_", grp, "_L8_Parent_CellType"), cellchat.obj, envir = .GlobalEnv)
  
}





levels_L2 <- levels(factor(as.character(tel$L2_Predicted_Region)))

for (grp in group_list) {
  #################################################
  ## L2_Predicted_Region
  #################################################
  
  obj <- tel
  obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay") 
  
  level <- 'L2_Predicted_Region'
  Idents(obj) <- 'L2_Predicted_Region'
  levels(obj)
  
  data.input <- GetAssayData(obj, assay = 'RNA', layer = 'data') # normalized data matrix
  
  data.input.tel <- as.matrix(data.input[homologs$mzebra_seurat_gene_id[which(!is.na(homologs$one_to_one_human))],])
  rownames(data.input.tel) <- homologs$one_to_one_human[which(!is.na(homologs$one_to_one_human))]
  head(rownames(data.input.tel))
  
  obj@assays$RNA@data <- data.input.tel
  head(rownames(obj@assays$RNA@data))
  
  gc()
  
  meta <- obj@meta.data
  cell.use = rownames(meta)[meta$group == grp]
  
  data.input.tel = data.input.tel[, cell.use]
  meta = meta[cell.use, ]
  unique(meta$L2_Predicted_Region)
  
  meta$samples <- as.factor(meta$individual)
  
  #################################################
  ## IMPORTANT: DROP unused levels BEFORE computeCommunProb()
  ## (computeCommunProb fails if idents contain zero-cell levels)
  #################################################
  
  meta$L2_Predicted_Region <- droplevels(factor(as.character(meta$L2_Predicted_Region)))
  
  cellchat.obj <- createCellChat(object = data.input.tel, group.by = "L2_Predicted_Region", assay = "RNA", meta = meta)
  
  cellchat.obj <- addMeta(cellchat.obj, meta = meta)
  cellchat.obj <- setIdent(cellchat.obj, ident.use = "L2_Predicted_Region")
  
  
  levels(cellchat.obj@idents)
  groupSize <- as.numeric(table(cellchat.obj@idents))
  groupSize <- groupSize[which(groupSize != 0)]
  
  gc()
  
  CellChatDB <- CellChatDB.human
  cellchat.obj@DB <- CellChatDB
  
  cellchat.obj <- subsetData(cellchat.obj)
  
  cellchat.obj_varFeats <- identifyOverExpressedGenes(
    cellchat.obj,
    thresh.p   = 0.01,
    min.cells  = 10,
    return.object = FALSE
  )
  length(unique(cellchat.obj_varFeats$features))
  
  outdir <- file.path("/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/mbuna_multiome_project/seurat_data/tel_ob/cellchat", grp, "predicted_region")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  write_csv(cellchat.obj_varFeats, file.path(outdir, paste0(level, "_cellchat.", grp, "_OverExpressedGenes.rds")))
  
  cellchat.obj <- identifyOverExpressedGenes(
    cellchat.obj,
    thresh.p   = 0.01,
    min.cells  = 10
  )
  
  cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
  
  data(PPI.human)
  PPI.human <- as(PPI.human, "sparseMatrix")
  cellchat.obj <- smoothData(cellchat.obj, adj = PPI.human)
  
  cellchat.obj <- computeCommunProb(
    cellchat.obj,
    population.size = FALSE,
    raw.use = FALSE
  )
  
  cellchat.obj <- computeCommunProbPathway(cellchat.obj)
  cellchat.obj <- aggregateNet(cellchat.obj)
  
  cellchat.obj_subset <- filterCommunication(
    cellchat.obj,
    min.cells     = 5,
    min.samples   = 3,
    rare.keep     = TRUE,
    nonFilter.keep = FALSE
  )
  
  #################################################
  ## PLOTTING FIX:
  ## Re-introduce missing regions ONLY FOR PLOTTING
  ## by padding the net matrices + groupSize.
  #################################################
  
  pad_net_matrix <- function(mat, all_nodes) {
    out <- matrix(
      0,
      nrow = length(all_nodes),
      ncol = length(all_nodes),
      dimnames = list(all_nodes, all_nodes)
    )
    common <- intersect(rownames(mat), all_nodes)
    out[common, common] <- mat[common, common, drop = FALSE]
    out
  }
  
  all_nodes <- levels_L2
  
  net.count.plot  <- pad_net_matrix(cellchat.obj_subset@net$count,  all_nodes)
  net.weight.plot <- pad_net_matrix(cellchat.obj_subset@net$weight, all_nodes)
  
  # Full-length groupSize INCLUDING 0-count regions (do NOT filter zeros!)
  groupSize_full <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize_full <- as.numeric(groupSize_full)
  
  #################################################
  ## (Optional but consistent) re-expand idents/meta labels
  ## for downstream ordering / ridge plots.
  #################################################
  
  # cellchat.obj@idents <- factor(as.character(cellchat.obj@idents), levels = levels_L2)
  # if (!is.null(cellchat.obj@meta[[level]])) {
  #   cellchat.obj@meta[[level]] <- factor(as.character(cellchat.obj@meta[[level]]), levels = levels_L2)
  # }
  # 
  # cellchat.obj_subset@idents <- factor(as.character(cellchat.obj_subset@idents), levels = levels_L2)
  # if (!is.null(cellchat.obj_subset@meta[[level]])) {
  #   cellchat.obj_subset@meta[[level]] <- factor(as.character(cellchat.obj_subset@meta[[level]]), levels = levels_L2)
  # }
  
  idents_ordered <- c(
    "3V",
    "AC/PF",
    "AH",
    "AH/POA",
    "AntPit",
    "ARC",
    "BL",
    "DMH",
    "ECM/PNN",
    "tel",
    "ID",
    "LH",
    "LH/PZ",
    "LH/VMH",
    "MBH",
    "ME",
    "ME/ARC",
    "ME/POA",
    "MMN",
    "MN",
    "MZ",
    "PF/MMT",
    "PF/MT",
    "PM",
    "POA",
    "PostPit",
    "PSTN",
    "PVN",
    "PVS",
    "RM",
    "SCN",
    "SON/PVN",
    "STN",
    "SuM",
    "VMH",
    "ZI"
  )
  
  # Apply to CellChat object
  cellchat.obj@idents <- factor(cellchat.obj@idents, levels = idents_ordered)
  
  color_key_region <- c(
    "3V"      = "#D49200",
    "AC/PF"   = "#DA81A8",
    "AH"      = "#1C7ED5",
    "AH/POA"  = "#91A7FF",
    "AntPit"  = "#8953a9",
    "ARC"     = "#4125F0",
    "BL"      = "#D631C7",
    "DMH"     = "#816FF0",
    "ECM/PNN" = "#F19C7C",
    "tel"    = "#20b3c4",
    "ID"      = "#8b6c5c",
    "LH"      = "#ff9acd",
    "LH/PZ"   = "#66c5bc",
    "LH/VMH"  = "#E6589F",
    "MBH"     = "#ACCA82",
    "ME"      = "#75975e",
    "ME/ARC"  = "goldenrod3",
    "ME/POA"  = "#0d8e84",
    "MMN"     = "#BA2D22",
    "MN"      = "salmon3",
    "MZ"      = "goldenrod1",
    "PF/MMT"  = "lightgoldenrod3",
    "PF/MT"   = "#F2D67F",
    "PM"      = "#729efd",
    "POA"     = "#EB6E40",
    "PostPit" = "#c098f2",
    "PSTN"    = "#E4928E",
    "PVN"     = "#E66055",
    "PVS"     = "#EC83E0",
    "RM"      = "#74C0FC",
    "SCN"     = "#B3544D",
    "SON/PVN" = "#F76060",
    "STN"     = "lightsteelblue3",
    "SuM"     = "indianred",
    "VMH"     = "#4444C4",
    "ZI"      = "#7865EF"
  )
  color_key_region_df <- data.frame(
    celltype = names(color_key_region),
    color    = unname(color_key_region),
    stringsAsFactors = FALSE
  )
  
  color.use.plot <- color_key_region[all_nodes]
  
  # (Recommended) apply to df.net for consistent ggplot facet/axis ordering
  df.net$source <- factor(df.net$source, levels = idents_ordered)
  df.net$target <- factor(df.net$target, levels = idents_ordered)
  
  # Quick checks
  levels(cellchat.obj@idents)
  all(levels(cellchat.obj@idents) == idents_ordered)
  
  pdf(file.path(outdir, paste0(grp, "_netCellCellInteractions_", level, '.pdf')), width = 16, height = 18)
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.count.plot, 
                       vertex.weight = groupSize, vertex.size.max	= 5,
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Number of Interactions: Top 10%", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.1,
                       angle_text = TRUE,
                       remove.isolate = F,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.weight.plot, 
                       vertex.weight = groupSize, vertex.size.max	= 5, 
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Interaction Weights/Strength: Top 10%", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.1,
                       angle_text = TRUE,
                       remove.isolate = F,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  pit_celltypes <- c(
    "PostPit",
    "AntPit"
  )
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.weight.plot, 
                       vertex.weight = groupSize, 
                       vertex.size.max	= 5,
                       weight.scale = T, 
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Targeting Pituitary Regions (interaction strength - Top 50%)", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.5, angle_text = TRUE,
                       targets.use = pit_celltypes,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  groupSize <- table(factor(as.character(cellchat.obj_subset@idents), levels = all_nodes))
  groupSize <- as.numeric(groupSize)
  
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_circle_mod(net.weight.plot, 
                       vertex.weight = groupSize, 
                       weight.scale = T, 
                       vertex.size.max	= 5,
                       label.edge= F, 
                       title.name = paste0('\n', grp, " Cell–Cell Signaling Originating from Pituitary Regions (interaction strength - Top 50%", '\n LEVEL: ', level, '\n\n\n'),
                       top = 0.5, angle_text = TRUE,
                       sources.use = pit_celltypes,
                       label.offset.dist = 0.76,
                       color.use = color.use.plot)
  
  dev.off()
  
  pdf(file.path(outdir, paste0(grp, "_netInteractionWeights_per_", level, '.pdf')), width = 25, height = 28)
  mat <- net.weight.plot
  par(mfrow = c(2,2), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle_mod(mat2, 
                         vertex.weight = groupSize, 
                         vertex.size.max	= 5, 
                         weight.scale = T, 
                         edge.weight.max = max(mat), 
                         angle_text = TRUE,
                         color.use = color.use.plot,
                         title.name = paste0(rownames(mat)[i], '\n\n\n\n'))
  }
  dev.off()
  
  saveRDS(
    cellchat.obj,
    file.path(outdir, paste0("cellchat_tel_", grp, "_L2_Predicted_Region.rds"))
  )
  
  # cellchat.obj <- netAnalysis_computeCentrality(
  #   cellchat.obj,
  #   slot.name = "netP"
  # )
  
  df.net <- subsetCommunication(cellchat.obj, thresh = 0.05)
  write_csv(df.net, file.path(outdir, paste0(level, "_cellchat.", grp, "_netCellCommunications.csv")))
  
  df.netP <- subsetCommunication(cellchat.obj, thresh = 0.05, slot.name = 'netP')
  write_csv(df.netP, file.path(outdir, paste0(level, "_cellchat.", grp, "_netSignalingPathways.csv")))
  
  
  gc()
  
  library(pheatmap)
  
  ################################
  ## 4. RIDGE PLOTS (FIXED)
  ################################
  
  library(ggridges)
  
  library(ggplot2)
  library(ggridges)
  
  # enforce your desired celltype order ONCE
  df.net$source <- factor(df.net$source, levels = idents_ordered)
  df.net$target <- factor(df.net$target, levels = idents_ordered)
  
  df.net$source_color <- as.character(df.net$source)
  df.net$source_color <- color_key_region_df$color[match(df.net$source, color_key_region_df$region)]
  
  df.net$target_color <- as.character(df.net$target)
  df.net$target_color <- color_key_region_df$color[match(df.net$target, color_key_region_df$celltype)]
  
  pdf(file.path(outdir,paste0(grp, "_netComminicationDistribution_", level, '.pdf')), width = 10, height = 20)
  
  ################################
  ## 4A. Ridge plot by sender (grid lines + grey ridge outlines)
  ################################
  
  ggplot(
    df.net,
    aes(
      x = prob,
      y = source,
      fill = source_color
    )
  ) +
    geom_density_ridges(
      scale = 2.4,
      alpha = 0.85,
      color = "gray",
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    scale_x_log10() +
    scale_y_discrete(
      limits = rev(idents_ordered),
      drop = FALSE
    ) +
    scale_fill_identity(guide = "none") +
    theme_ridges(grid = FALSE) +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      text = element_text(size = 13),
      axis.text = element_text(size = 11),
      axis.title.x = element_text(size = 13, hjust = 0.5),
      axis.title.y = element_text(size = 13, hjust = 0.5),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = paste0(grp, ' - ', level, ": CellChat communication probability distributions by sender"),
      x = "Communication probability (log10 scale)",
      y = "Sender cell type"
    )
  
  ################################
  ## 4B. Ridge plot by receiver (grid lines + grey ridge outlines)
  ################################
  
  ggplot(
    df.net,
    aes(
      x = prob,
      y = target,
      fill = target_color
    )
  ) +
    geom_density_ridges(
      scale = 2.4,
      alpha = 0.85,
      color = "gray",
      linewidth = 0.3,
      na.rm = TRUE
    ) +
    scale_x_log10() +
    scale_y_discrete(
      limits = rev(idents_ordered),
      drop = FALSE
    ) +
    scale_fill_identity(guide = "none") +
    theme_ridges(grid = FALSE) +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      text = element_text(size = 13),
      axis.text = element_text(size = 11),
      axis.title.x = element_text(size = 13, hjust = 0.5),
      axis.title.y = element_text(size = 13, hjust = 0.5),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = paste0(grp, ' - ', level, ": CellChat communication probability distributions by receiver"),
      x = "Communication probability (log10 scale)",
      y = "Receiver cell type"
    )
  
  ################################
  ## 4C. Ridge plot by signaling pathway (grid lines + grey ridge outlines)
  ################################
  
  ggplot(
    df.net,
    aes(
      x = prob,
      y = reorder(pathway_name, prob, FUN = median),
      fill = pathway_name
    )
  ) +
    geom_density_ridges(
      scale = 2.6,
      alpha = 0.85,
      color = "gray",
      linewidth = 0.3
    ) +
    scale_x_log10() +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    theme_ridges(grid = FALSE) +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 11),
      axis.title.x = element_text(size = 14, hjust = 0.5),
      axis.title.y = element_text(size = 14, hjust = 0.5),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5)
    ) +
    labs(
      title = paste0(grp, ' - ', level, ": Distribution of CellChat communication probabilities by pathway"),
      subtitle = paste0(grp, ' - ', level, ": Pathways ordered by median probability"),
      x = "Communication probability (log10 scale)",
      y = "Signaling pathway"
    )
  
  df.net %>%
    group_by(source, source_color) %>%
    summarise(total_outgoing = sum(prob), .groups = "drop") %>%
    ggplot(
      aes(
        x = reorder(source, total_outgoing),
        y = total_outgoing,
        fill = source_color
      )
    ) +
    geom_col() +
    scale_fill_identity() +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.x  = element_text(color = "black"),
      axis.text.y  = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title   = element_text(color = "black")
    ) +
    labs(
      title = "Total outgoing cell–cell signaling by source cell type",
      x = "Source (signaling) cell type",
      y = "Summed communication probability"
    )
  
  df.net %>%
    group_by(target, target_color) %>%
    summarise(total_incoming = sum(prob), .groups = "drop") %>%
    ggplot(
      aes(
        x = reorder(target, total_incoming),
        y = total_incoming,
        fill = target_color
      )
    ) +
    geom_col() +
    scale_fill_identity() +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.x  = element_text(color = "black"),
      axis.text.y  = element_text(color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      plot.title   = element_text(color = "black")
    ) +
    labs(
      title = "Total incoming cell–cell signaling by target cell type",
      x = "Target (receiving) cell type",
      y = "Summed communication probability"
    )
  
  dev.off()
  
  
  pdf(file.path(outdir,paste0(grp, "_netComminicationStrength_heatmap_", level, '.pdf')), width = 18, height = 15)
  mat <- df.net %>%
    group_by(source, target) %>%
    summarise(weight = sum(prob), .groups = "drop") %>%
    pivot_wider(names_from = target, values_from = weight, values_fill = 0) %>%
    column_to_rownames("source") %>%
    as.matrix()
  
  pheatmap(
    mat,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#F7FBFF", 'grey90', 'lightblue1', "lightblue", "#41B6C4", 'deepskyblue3', "#225EA8", "#253494"))(100),
    main = "Cell–cell communication strength (summed probabilities)\nRows = source (senders), columns = target (receivers)",
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45
  )
  dev.off()
  
  library(igraph)
  library(ggraph)
  
  pdf(file.path(outdir,paste0(grp, "_LigandReceptor_interactions groupedBy_signaling.pathway", level, '.pdf')), width = 35, height = 90)
  
  colors <- c("lightblue", "#41B6C4", 'deepskyblue3', "#225EA8", "#253494")
  
  df.net %>%
    group_by(interaction_name, pathway_name) %>%
    summarise(
      mean_prob = mean(prob),
      n_pairs = n(),
      .groups = "drop"
    ) %>%
    ggplot(aes(
      x = pathway_name,
      y = interaction_name,
      size = mean_prob,
      color = n_pairs
    )) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(trans = "log10") +
    scale_color_gradientn(
      colors = colors
    ) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_line(color = "gray"),
      panel.grid.major.x = element_line(color = "gray"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.text.y   = element_text(size = 15, color = "black"),
      axis.text.x   = element_text(size = 15,color = "black"),
      axis.title.x  = element_text(size = 15,color = "black"),
      axis.title.y  = element_text(size = 15,color = "black"),
      legend.text   = element_text(size = 15,color = "black"),
      legend.title  = element_text(size = 15,color = "black"),
      plot.title    = element_text(size = 18, color = "black")
    ) +
    labs(
      title = "Ligand–receptor interactions grouped by signaling pathway",
      x = "Signaling pathway",
      y = "Ligand–receptor pair",
      size = "Mean communication probability",
      color = "# interactions"
    ) +
    RotatedAxis()
  
  
  dev.off()
  
  saveRDS(
    cellchat.obj,
    file.path(outdir, paste0("cellchat_tel_", grp, "_L2_Predicted_Region_netCentrality.rds"))
  )
  
  assign(paste0("cellchat_tel_", grp, "_L2_Predicted_Region"), cellchat.obj, envir = .GlobalEnv)
}
