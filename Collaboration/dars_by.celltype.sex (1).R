
peak.feats <- head(VariableFeatures(tel_final[['peaks']]), 20000)


gc()
gc()

DefaultAssay(tel_final) <- 'peaks'
Idents(tel_final) <- "celltype.sex"
table(Idents(tel_final))

BgrB_male_peak.markers_1_Astro.MG <- FindMarkers(
  object = tel_final,
  ident.1 = "1_Astro/MG_MALE",
  ident.2 = "1_Astro/MG_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_1_Astro.MG$pct.diff <- ((abs(BgrB_male_peak.markers_1_Astro.MG$pct.1 - BgrB_male_peak.markers_1_Astro.MG$pct.2)) / ((BgrB_male_peak.markers_1_Astro.MG$pct.1 + BgrB_male_peak.markers_1_Astro.MG$pct.2)/2)) * 100
BgrB_male_peak.markers_1_Astro.MG <- BgrB_male_peak.markers_1_Astro.MG[order(BgrB_male_peak.markers_1_Astro.MG$pct.diff, decreasing = T),]


BgrB_male_peak.markers_1_Astro.MG$CellType <- "1_Astro/MG"
BgrB_male_peak.markers_1_Astro.MG$sex <- "MALE"
BgrB_male_peak.markers_1_Astro.MG$Gene <- rownames(BgrB_male_peak.markers_1_Astro.MG)
BgrB_male_peak.markers_1_Astro.MG <- BgrB_male_peak.markers_1_Astro.MG %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_1_Astro.MG <- FindMarkers(
  object = tel_final,
  ident.1 = "1_Astro/MG_FEMALE",
  ident.2 = "1_Astro/MG_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_1_Astro.MG$pct.diff <- ((abs(BgrB_female_peak.markers_1_Astro.MG$pct.1 - BgrB_female_peak.markers_1_Astro.MG$pct.2)) / ((BgrB_female_peak.markers_1_Astro.MG$pct.1 + BgrB_female_peak.markers_1_Astro.MG$pct.2)/2)) * 100
BgrB_female_peak.markers_1_Astro.MG <- BgrB_female_peak.markers_1_Astro.MG[order(BgrB_female_peak.markers_1_Astro.MG$pct.diff, decreasing = T),]


BgrB_female_peak.markers_1_Astro.MG$CellType <- "1_Astro/MG"
BgrB_female_peak.markers_1_Astro.MG$sex <- "FEMALE"
BgrB_female_peak.markers_1_Astro.MG$Gene <- rownames(BgrB_female_peak.markers_1_Astro.MG)
BgrB_female_peak.markers_1_Astro.MG <- BgrB_female_peak.markers_1_Astro.MG %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_1_Astro.MG <- rbind(BgrB_male_peak.markers_1_Astro.MG, BgrB_female_peak.markers_1_Astro.MG)

saveRDS(BgrB_sex_peak.markers_1_Astro.MG, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_1_Astro.MG.rds')


BgrB_male_peak.markers_2_OPC_Oligo <- FindMarkers(
  object = tel_final,
  ident.1 = "2_OPC/Oligo_MALE",
  ident.2 = "2_OPC/Oligo_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_2_OPC_Oligo$pct.diff <- ((abs(BgrB_male_peak.markers_2_OPC_Oligo$pct.1 - BgrB_male_peak.markers_2_OPC_Oligo$pct.2)) / ((BgrB_male_peak.markers_2_OPC_Oligo$pct.1 + BgrB_male_peak.markers_2_OPC_Oligo$pct.2)/2)) * 100
BgrB_male_peak.markers_2_OPC_Oligo <- BgrB_male_peak.markers_2_OPC_Oligo[order(BgrB_male_peak.markers_2_OPC_Oligo$pct.diff, decreasing = T),]


BgrB_male_peak.markers_2_OPC_Oligo$CellType <- "2_OPC/Oligo"
BgrB_male_peak.markers_2_OPC_Oligo$sex <- "MALE"
BgrB_male_peak.markers_2_OPC_Oligo$Gene <- rownames(BgrB_male_peak.markers_2_OPC_Oligo)
BgrB_male_peak.markers_2_OPC_Oligo <- BgrB_male_peak.markers_2_OPC_Oligo %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_2_OPC_Oligo <- FindMarkers(
  object = tel_final,
  ident.1 = "2_OPC/Oligo_FEMALE",
  ident.2 = "2_OPC/Oligo_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_2_OPC_Oligo$pct.diff <- ((abs(BgrB_female_peak.markers_2_OPC_Oligo$pct.1 - BgrB_female_peak.markers_2_OPC_Oligo$pct.2)) / ((BgrB_female_peak.markers_2_OPC_Oligo$pct.1 + BgrB_female_peak.markers_2_OPC_Oligo$pct.2)/2)) * 100
BgrB_female_peak.markers_2_OPC_Oligo <- BgrB_female_peak.markers_2_OPC_Oligo[order(BgrB_female_peak.markers_2_OPC_Oligo$pct.diff, decreasing = T),]


BgrB_female_peak.markers_2_OPC_Oligo$CellType <- "2_OPC/Oligo"
BgrB_female_peak.markers_2_OPC_Oligo$sex <- "FEMALE"
BgrB_female_peak.markers_2_OPC_Oligo$Gene <- rownames(BgrB_female_peak.markers_2_OPC_Oligo)
BgrB_female_peak.markers_2_OPC_Oligo <- BgrB_female_peak.markers_2_OPC_Oligo %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_2_OPC_Oligo <- rbind(BgrB_male_peak.markers_2_OPC_Oligo, BgrB_female_peak.markers_2_OPC_Oligo)

saveRDS(BgrB_sex_peak.markers_2_OPC_Oligo, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_2_OPC_Oligo.rds')


BgrB_male_peak.markers_3_Peri <- FindMarkers(
  object = tel_final,
  ident.1 = "3_Peri_MALE",
  ident.2 = "3_Peri_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_3_Peri$pct.diff <- ((abs(BgrB_male_peak.markers_3_Peri$pct.1 - BgrB_male_peak.markers_3_Peri$pct.2)) / ((BgrB_male_peak.markers_3_Peri$pct.1 + BgrB_male_peak.markers_3_Peri$pct.2)/2)) * 100
BgrB_male_peak.markers_3_Peri <- BgrB_male_peak.markers_3_Peri[order(BgrB_male_peak.markers_3_Peri$pct.diff, decreasing = T),]


BgrB_male_peak.markers_3_Peri$CellType <- "3_Peri"
BgrB_male_peak.markers_3_Peri$sex <- "MALE"
BgrB_male_peak.markers_3_Peri$Gene <- rownames(BgrB_male_peak.markers_3_Peri)
BgrB_male_peak.markers_3_Peri <- BgrB_male_peak.markers_3_Peri %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_3_Peri <- FindMarkers(
  object = tel_final,
  ident.1 = "3_Peri_FEMALE",
  ident.2 = "3_Peri_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_3_Peri$pct.diff <- ((abs(BgrB_female_peak.markers_3_Peri$pct.1 - BgrB_female_peak.markers_3_Peri$pct.2)) / ((BgrB_female_peak.markers_3_Peri$pct.1 + BgrB_female_peak.markers_3_Peri$pct.2)/2)) * 100
BgrB_female_peak.markers_3_Peri <- BgrB_female_peak.markers_3_Peri[order(BgrB_female_peak.markers_3_Peri$pct.diff, decreasing = T),]


BgrB_female_peak.markers_3_Peri$CellType <- "3_Peri"
BgrB_female_peak.markers_3_Peri$sex <- "FEMALE"
BgrB_female_peak.markers_3_Peri$Gene <- rownames(BgrB_female_peak.markers_3_Peri)
BgrB_female_peak.markers_3_Peri <- BgrB_female_peak.markers_3_Peri %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_3_Peri <- rbind(BgrB_male_peak.markers_3_Peri, BgrB_female_peak.markers_3_Peri)

saveRDS(BgrB_sex_peak.markers_3_Peri, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_3_Peri.rds')


BgrB_male_peak.markers_4_In <- FindMarkers(
  object = tel_final,
  ident.1 = "4_In_MALE",
  ident.2 = "4_In_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_4_In$pct.diff <- ((abs(BgrB_male_peak.markers_4_In$pct.1 - BgrB_male_peak.markers_4_In$pct.2)) / ((BgrB_male_peak.markers_4_In$pct.1 + BgrB_male_peak.markers_4_In$pct.2)/2)) * 100
BgrB_male_peak.markers_4_In <- BgrB_male_peak.markers_4_In[order(BgrB_male_peak.markers_4_In$pct.diff, decreasing = T),]


BgrB_male_peak.markers_4_In$CellType <- "4_In"
BgrB_male_peak.markers_4_In$sex <- "MALE"
BgrB_male_peak.markers_4_In$Gene <- rownames(BgrB_male_peak.markers_4_In)
BgrB_male_peak.markers_4_In <- BgrB_male_peak.markers_4_In %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_4_In <- FindMarkers(
  object = tel_final,
  ident.1 = "4_In_FEMALE",
  ident.2 = "4_In_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_4_In$pct.diff <- ((abs(BgrB_female_peak.markers_4_In$pct.1 - BgrB_female_peak.markers_4_In$pct.2)) / ((BgrB_female_peak.markers_4_In$pct.1 + BgrB_female_peak.markers_4_In$pct.2)/2)) * 100
BgrB_female_peak.markers_4_In <- BgrB_female_peak.markers_4_In[order(BgrB_female_peak.markers_4_In$pct.diff, decreasing = T),]


BgrB_female_peak.markers_4_In$CellType <- "4_In"
BgrB_female_peak.markers_4_In$sex <- "FEMALE"
BgrB_female_peak.markers_4_In$Gene <- rownames(BgrB_female_peak.markers_4_In)
BgrB_female_peak.markers_4_In <- BgrB_female_peak.markers_4_In %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_4_In <- rbind(BgrB_male_peak.markers_4_In, BgrB_female_peak.markers_4_In)

saveRDS(BgrB_sex_peak.markers_4_In, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_4_In.rds')


BgrB_male_peak.markers_5_In <- FindMarkers(
  object = tel_final,
  ident.1 = "5_In_MALE",
  ident.2 = "5_In_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_5_In$pct.diff <- ((abs(BgrB_male_peak.markers_5_In$pct.1 - BgrB_male_peak.markers_5_In$pct.2)) / ((BgrB_male_peak.markers_5_In$pct.1 + BgrB_male_peak.markers_5_In$pct.2)/2)) * 100
BgrB_male_peak.markers_5_In <- BgrB_male_peak.markers_5_In[order(BgrB_male_peak.markers_5_In$pct.diff, decreasing = T),]


BgrB_male_peak.markers_5_In$CellType <- "5_In"
BgrB_male_peak.markers_5_In$sex <- "MALE"
BgrB_male_peak.markers_5_In$Gene <- rownames(BgrB_male_peak.markers_5_In)
BgrB_male_peak.markers_5_In <- BgrB_male_peak.markers_5_In %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_5_In <- FindMarkers(
  object = tel_final,
  ident.1 = "5_In_FEMALE",
  ident.2 = "5_In_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_5_In$pct.diff <- ((abs(BgrB_female_peak.markers_5_In$pct.1 - BgrB_female_peak.markers_5_In$pct.2)) / ((BgrB_female_peak.markers_5_In$pct.1 + BgrB_female_peak.markers_5_In$pct.2)/2)) * 100
BgrB_female_peak.markers_5_In <- BgrB_female_peak.markers_5_In[order(BgrB_female_peak.markers_5_In$pct.diff, decreasing = T),]


BgrB_female_peak.markers_5_In$CellType <- "5_In"
BgrB_female_peak.markers_5_In$sex <- "FEMALE"
BgrB_female_peak.markers_5_In$Gene <- rownames(BgrB_female_peak.markers_5_In)
BgrB_female_peak.markers_5_In <- BgrB_female_peak.markers_5_In %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_5_In <- rbind(BgrB_male_peak.markers_5_In, BgrB_female_peak.markers_5_In)

saveRDS(BgrB_sex_peak.markers_5_In, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_5_In.rds')

BgrB_male_peak.markers_6_In <- FindMarkers(
  object = tel_final,
  ident.1 = "6_In_MALE",
  ident.2 = "6_In_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_6_In$pct.diff <- ((abs(BgrB_male_peak.markers_6_In$pct.1 - BgrB_male_peak.markers_6_In$pct.2)) / ((BgrB_male_peak.markers_6_In$pct.1 + BgrB_male_peak.markers_6_In$pct.2)/2)) * 100
BgrB_male_peak.markers_6_In <- BgrB_male_peak.markers_6_In[order(BgrB_male_peak.markers_6_In$pct.diff, decreasing = T),]


BgrB_male_peak.markers_6_In$CellType <- "6_In"
BgrB_male_peak.markers_6_In$sex <- "MALE"
BgrB_male_peak.markers_6_In$Gene <- rownames(BgrB_male_peak.markers_6_In)
BgrB_male_peak.markers_6_In <- BgrB_male_peak.markers_6_In %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_6_In <- FindMarkers(
  object = tel_final,
  ident.1 = "6_In_FEMALE",
  ident.2 = "6_In_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_6_In$pct.diff <- ((abs(BgrB_female_peak.markers_6_In$pct.1 - BgrB_female_peak.markers_6_In$pct.2)) / ((BgrB_female_peak.markers_6_In$pct.1 + BgrB_female_peak.markers_6_In$pct.2)/2)) * 100
BgrB_female_peak.markers_6_In <- BgrB_female_peak.markers_6_In[order(BgrB_female_peak.markers_6_In$pct.diff, decreasing = T),]


BgrB_female_peak.markers_6_In$CellType <- "6_In"
BgrB_female_peak.markers_6_In$sex <- "FEMALE"
BgrB_female_peak.markers_6_In$Gene <- rownames(BgrB_female_peak.markers_6_In)
BgrB_female_peak.markers_6_In <- BgrB_female_peak.markers_6_In %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_6_In <- rbind(BgrB_male_peak.markers_6_In, BgrB_female_peak.markers_6_In)

saveRDS(BgrB_sex_peak.markers_6_In, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_6_In.rds')

BgrB_male_peak.markers_7_In <- FindMarkers(
  object = tel_final,
  ident.1 = "7_In_MALE",
  ident.2 = "7_In_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_7_In$pct.diff <- ((abs(BgrB_male_peak.markers_7_In$pct.1 - BgrB_male_peak.markers_7_In$pct.2)) / ((BgrB_male_peak.markers_7_In$pct.1 + BgrB_male_peak.markers_7_In$pct.2)/2)) * 100
BgrB_male_peak.markers_7_In <- BgrB_male_peak.markers_7_In[order(BgrB_male_peak.markers_7_In$pct.diff, decreasing = T),]


BgrB_male_peak.markers_7_In$CellType <- "7_In"
BgrB_male_peak.markers_7_In$sex <- "MALE"
BgrB_male_peak.markers_7_In$Gene <- rownames(BgrB_male_peak.markers_7_In)
BgrB_male_peak.markers_7_In <- BgrB_male_peak.markers_7_In %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_7_In <- FindMarkers(
  object = tel_final,
  ident.1 = "7_In_FEMALE",
  ident.2 = "7_In_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_7_In$pct.diff <- ((abs(BgrB_female_peak.markers_7_In$pct.1 - BgrB_female_peak.markers_7_In$pct.2)) / ((BgrB_female_peak.markers_7_In$pct.1 + BgrB_female_peak.markers_7_In$pct.2)/2)) * 100
BgrB_female_peak.markers_7_In <- BgrB_female_peak.markers_7_In[order(BgrB_female_peak.markers_7_In$pct.diff, decreasing = T),]


BgrB_female_peak.markers_7_In$CellType <- "7_In"
BgrB_female_peak.markers_7_In$sex <- "FEMALE"
BgrB_female_peak.markers_7_In$Gene <- rownames(BgrB_female_peak.markers_7_In)
BgrB_female_peak.markers_7_In <- BgrB_female_peak.markers_7_In %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_7_In <- rbind(BgrB_male_peak.markers_7_In, BgrB_female_peak.markers_7_In)

saveRDS(BgrB_sex_peak.markers_7_In, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_7_In.rds')

BgrB_male_peak.markers_8_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "8_Ex_MALE",
  ident.2 = "8_Ex_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_8_Ex$pct.diff <- ((abs(BgrB_male_peak.markers_8_Ex$pct.1 - BgrB_male_peak.markers_8_Ex$pct.2)) / ((BgrB_male_peak.markers_8_Ex$pct.1 + BgrB_male_peak.markers_8_Ex$pct.2)/2)) * 100
BgrB_male_peak.markers_8_Ex <- BgrB_male_peak.markers_8_Ex[order(BgrB_male_peak.markers_8_Ex$pct.diff, decreasing = T),]


BgrB_male_peak.markers_8_Ex$CellType <- "8_Ex"
BgrB_male_peak.markers_8_Ex$sex <- "MALE"
BgrB_male_peak.markers_8_Ex$Gene <- rownames(BgrB_male_peak.markers_8_Ex)
BgrB_male_peak.markers_8_Ex <- BgrB_male_peak.markers_8_Ex %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_8_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "8_Ex_FEMALE",
  ident.2 = "8_Ex_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_8_Ex$pct.diff <- ((abs(BgrB_female_peak.markers_8_Ex$pct.1 - BgrB_female_peak.markers_8_Ex$pct.2)) / ((BgrB_female_peak.markers_8_Ex$pct.1 + BgrB_female_peak.markers_8_Ex$pct.2)/2)) * 100
BgrB_female_peak.markers_8_Ex <- BgrB_female_peak.markers_8_Ex[order(BgrB_female_peak.markers_8_Ex$pct.diff, decreasing = T),]


BgrB_female_peak.markers_8_Ex$CellType <- "8_Ex"
BgrB_female_peak.markers_8_Ex$sex <- "FEMALE"
BgrB_female_peak.markers_8_Ex$Gene <- rownames(BgrB_female_peak.markers_8_Ex)
BgrB_female_peak.markers_8_Ex <- BgrB_female_peak.markers_8_Ex %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_8_Ex <- rbind(BgrB_male_peak.markers_8_Ex, BgrB_female_peak.markers_8_Ex)

saveRDS(BgrB_sex_peak.markers_8_Ex, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_8_Ex.rds')

BgrB_male_peak.markers_9_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "9_Ex_MALE",
  ident.2 = "9_Ex_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_9_Ex$pct.diff <- ((abs(BgrB_male_peak.markers_9_Ex$pct.1 - BgrB_male_peak.markers_9_Ex$pct.2)) / ((BgrB_male_peak.markers_9_Ex$pct.1 + BgrB_male_peak.markers_9_Ex$pct.2)/2)) * 100
BgrB_male_peak.markers_9_Ex <- BgrB_male_peak.markers_9_Ex[order(BgrB_male_peak.markers_9_Ex$pct.diff, decreasing = T),]


BgrB_male_peak.markers_9_Ex$CellType <- "9_Ex"
BgrB_male_peak.markers_9_Ex$sex <- "MALE"
BgrB_male_peak.markers_9_Ex$Gene <- rownames(BgrB_male_peak.markers_9_Ex)
BgrB_male_peak.markers_9_Ex <- BgrB_male_peak.markers_9_Ex %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_9_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "9_Ex_FEMALE",
  ident.2 = "9_Ex_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_9_Ex$pct.diff <- ((abs(BgrB_female_peak.markers_9_Ex$pct.1 - BgrB_female_peak.markers_9_Ex$pct.2)) / ((BgrB_female_peak.markers_9_Ex$pct.1 + BgrB_female_peak.markers_9_Ex$pct.2)/2)) * 100
BgrB_female_peak.markers_9_Ex <- BgrB_female_peak.markers_9_Ex[order(BgrB_female_peak.markers_9_Ex$pct.diff, decreasing = T),]


BgrB_female_peak.markers_9_Ex$CellType <- "9_Ex"
BgrB_female_peak.markers_9_Ex$sex <- "FEMALE"
BgrB_female_peak.markers_9_Ex$Gene <- rownames(BgrB_female_peak.markers_9_Ex)
BgrB_female_peak.markers_9_Ex <- BgrB_female_peak.markers_9_Ex %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_9_Ex <- rbind(BgrB_male_peak.markers_9_Ex, BgrB_female_peak.markers_9_Ex)

saveRDS(BgrB_sex_peak.markers_9_Ex, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_9_Ex.rds')

BgrB_male_peak.markers_10_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "10_Ex_MALE",
  ident.2 = "10_Ex_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_10_Ex$pct.diff <- ((abs(BgrB_male_peak.markers_10_Ex$pct.1 - BgrB_male_peak.markers_10_Ex$pct.2)) / ((BgrB_male_peak.markers_10_Ex$pct.1 + BgrB_male_peak.markers_10_Ex$pct.2)/2)) * 100
BgrB_male_peak.markers_10_Ex <- BgrB_male_peak.markers_10_Ex[order(BgrB_male_peak.markers_10_Ex$pct.diff, decreasing = T),]


BgrB_male_peak.markers_10_Ex$CellType <- "10_Ex"
BgrB_male_peak.markers_10_Ex$sex <- "MALE"
BgrB_male_peak.markers_10_Ex$Gene <- rownames(BgrB_male_peak.markers_10_Ex)
BgrB_male_peak.markers_10_Ex <- BgrB_male_peak.markers_10_Ex %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_10_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "10_Ex_FEMALE",
  ident.2 = "10_Ex_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_10_Ex$pct.diff <- ((abs(BgrB_female_peak.markers_10_Ex$pct.1 - BgrB_female_peak.markers_10_Ex$pct.2)) / ((BgrB_female_peak.markers_10_Ex$pct.1 + BgrB_female_peak.markers_10_Ex$pct.2)/2)) * 100
BgrB_female_peak.markers_10_Ex <- BgrB_female_peak.markers_10_Ex[order(BgrB_female_peak.markers_10_Ex$pct.diff, decreasing = T),]


BgrB_female_peak.markers_10_Ex$CellType <- "10_Ex"
BgrB_female_peak.markers_10_Ex$sex <- "FEMALE"
BgrB_female_peak.markers_10_Ex$Gene <- rownames(BgrB_female_peak.markers_10_Ex)
BgrB_female_peak.markers_10_Ex <- BgrB_female_peak.markers_10_Ex %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_10_Ex <- rbind(BgrB_male_peak.markers_10_Ex, BgrB_female_peak.markers_10_Ex)

saveRDS(BgrB_sex_peak.markers_10_Ex, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_10_Ex.rds')

BgrB_male_peak.markers_11_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "11_Ex_MALE",
  ident.2 = "11_Ex_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_11_Ex$pct.diff <- ((abs(BgrB_male_peak.markers_11_Ex$pct.1 - BgrB_male_peak.markers_11_Ex$pct.2)) / ((BgrB_male_peak.markers_11_Ex$pct.1 + BgrB_male_peak.markers_11_Ex$pct.2)/2)) * 100
BgrB_male_peak.markers_11_Ex <- BgrB_male_peak.markers_11_Ex[order(BgrB_male_peak.markers_11_Ex$pct.diff, decreasing = T),]


BgrB_male_peak.markers_11_Ex$CellType <- "11_Ex"
BgrB_male_peak.markers_11_Ex$sex <- "MALE"
BgrB_male_peak.markers_11_Ex$Gene <- rownames(BgrB_male_peak.markers_11_Ex)
BgrB_male_peak.markers_11_Ex <- BgrB_male_peak.markers_11_Ex %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_11_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "11_Ex_FEMALE",
  ident.2 = "11_Ex_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_11_Ex$pct.diff <- ((abs(BgrB_female_peak.markers_11_Ex$pct.1 - BgrB_female_peak.markers_11_Ex$pct.2)) / ((BgrB_female_peak.markers_11_Ex$pct.1 + BgrB_female_peak.markers_11_Ex$pct.2)/2)) * 100
BgrB_female_peak.markers_11_Ex <- BgrB_female_peak.markers_11_Ex[order(BgrB_female_peak.markers_11_Ex$pct.diff, decreasing = T),]


BgrB_female_peak.markers_11_Ex$CellType <- "11_Ex"
BgrB_female_peak.markers_11_Ex$sex <- "FEMALE"
BgrB_female_peak.markers_11_Ex$Gene <- rownames(BgrB_female_peak.markers_11_Ex)
BgrB_female_peak.markers_11_Ex <- BgrB_female_peak.markers_11_Ex %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_11_Ex <- rbind(BgrB_male_peak.markers_11_Ex, BgrB_female_peak.markers_11_Ex)

saveRDS(BgrB_sex_peak.markers_11_Ex, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_11_Ex.rds')

BgrB_male_peak.markers_12_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "12_Ex_MALE",
  ident.2 = "12_Ex_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_12_Ex$pct.diff <- ((abs(BgrB_male_peak.markers_12_Ex$pct.1 - BgrB_male_peak.markers_12_Ex$pct.2)) / ((BgrB_male_peak.markers_12_Ex$pct.1 + BgrB_male_peak.markers_12_Ex$pct.2)/2)) * 100
BgrB_male_peak.markers_12_Ex <- BgrB_male_peak.markers_12_Ex[order(BgrB_male_peak.markers_12_Ex$pct.diff, decreasing = T),]


BgrB_male_peak.markers_12_Ex$CellType <- "12_Ex"
BgrB_male_peak.markers_12_Ex$sex <- "MALE"
BgrB_male_peak.markers_12_Ex$Gene <- rownames(BgrB_male_peak.markers_12_Ex)
BgrB_male_peak.markers_12_Ex <- BgrB_male_peak.markers_12_Ex %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_12_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "12_Ex_FEMALE",
  ident.2 = "12_Ex_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_12_Ex$pct.diff <- ((abs(BgrB_female_peak.markers_12_Ex$pct.1 - BgrB_female_peak.markers_12_Ex$pct.2)) / ((BgrB_female_peak.markers_12_Ex$pct.1 + BgrB_female_peak.markers_12_Ex$pct.2)/2)) * 100
BgrB_female_peak.markers_12_Ex <- BgrB_female_peak.markers_12_Ex[order(BgrB_female_peak.markers_12_Ex$pct.diff, decreasing = T),]


BgrB_female_peak.markers_12_Ex$CellType <- "12_Ex"
BgrB_female_peak.markers_12_Ex$sex <- "FEMALE"
BgrB_female_peak.markers_12_Ex$Gene <- rownames(BgrB_female_peak.markers_12_Ex)
BgrB_female_peak.markers_12_Ex <- BgrB_female_peak.markers_12_Ex %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_12_Ex <- rbind(BgrB_male_peak.markers_12_Ex, BgrB_female_peak.markers_12_Ex)

saveRDS(BgrB_sex_peak.markers_12_Ex, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_12_Ex.rds')

BgrB_male_peak.markers_13_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "13_Ex_MALE",
  ident.2 = "13_Ex_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_13_Ex$pct.diff <- ((abs(BgrB_male_peak.markers_13_Ex$pct.1 - BgrB_male_peak.markers_13_Ex$pct.2)) / ((BgrB_male_peak.markers_13_Ex$pct.1 + BgrB_male_peak.markers_13_Ex$pct.2)/2)) * 100
BgrB_male_peak.markers_13_Ex <- BgrB_male_peak.markers_13_Ex[order(BgrB_male_peak.markers_13_Ex$pct.diff, decreasing = T),]


BgrB_male_peak.markers_13_Ex$CellType <- "13_Ex"
BgrB_male_peak.markers_13_Ex$sex <- "MALE"
BgrB_male_peak.markers_13_Ex$Gene <- rownames(BgrB_male_peak.markers_13_Ex)
BgrB_male_peak.markers_13_Ex <- BgrB_male_peak.markers_13_Ex %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_13_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "13_Ex_FEMALE",
  ident.2 = "13_Ex_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_13_Ex$pct.diff <- ((abs(BgrB_female_peak.markers_13_Ex$pct.1 - BgrB_female_peak.markers_13_Ex$pct.2)) / ((BgrB_female_peak.markers_13_Ex$pct.1 + BgrB_female_peak.markers_13_Ex$pct.2)/2)) * 100
BgrB_female_peak.markers_13_Ex <- BgrB_female_peak.markers_13_Ex[order(BgrB_female_peak.markers_13_Ex$pct.diff, decreasing = T),]


BgrB_female_peak.markers_13_Ex$CellType <- "13_Ex"
BgrB_female_peak.markers_13_Ex$sex <- "FEMALE"
BgrB_female_peak.markers_13_Ex$Gene <- rownames(BgrB_female_peak.markers_13_Ex)
BgrB_female_peak.markers_13_Ex <- BgrB_female_peak.markers_13_Ex %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_13_Ex <- rbind(BgrB_male_peak.markers_13_Ex, BgrB_female_peak.markers_13_Ex)

saveRDS(BgrB_sex_peak.markers_13_Ex, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_13_Ex.rds')

BgrB_male_peak.markers_14_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "14_Ex_MALE",
  ident.2 = "14_Ex_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_14_Ex$pct.diff <- ((abs(BgrB_male_peak.markers_14_Ex$pct.1 - BgrB_male_peak.markers_14_Ex$pct.2)) / ((BgrB_male_peak.markers_14_Ex$pct.1 + BgrB_male_peak.markers_14_Ex$pct.2)/2)) * 100
BgrB_male_peak.markers_14_Ex <- BgrB_male_peak.markers_14_Ex[order(BgrB_male_peak.markers_14_Ex$pct.diff, decreasing = T),]


BgrB_male_peak.markers_14_Ex$CellType <- "14_Ex"
BgrB_male_peak.markers_14_Ex$sex <- "MALE"
BgrB_male_peak.markers_14_Ex$Gene <- rownames(BgrB_male_peak.markers_14_Ex)
BgrB_male_peak.markers_14_Ex <- BgrB_male_peak.markers_14_Ex %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_14_Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "14_Ex_FEMALE",
  ident.2 = "14_Ex_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_14_Ex$pct.diff <- ((abs(BgrB_female_peak.markers_14_Ex$pct.1 - BgrB_female_peak.markers_14_Ex$pct.2)) / ((BgrB_female_peak.markers_14_Ex$pct.1 + BgrB_female_peak.markers_14_Ex$pct.2)/2)) * 100
BgrB_female_peak.markers_14_Ex <- BgrB_female_peak.markers_14_Ex[order(BgrB_female_peak.markers_14_Ex$pct.diff, decreasing = T),]


BgrB_female_peak.markers_14_Ex$CellType <- "14_Ex"
BgrB_female_peak.markers_14_Ex$sex <- "FEMALE"
BgrB_female_peak.markers_14_Ex$Gene <- rownames(BgrB_female_peak.markers_14_Ex)
BgrB_female_peak.markers_14_Ex <- BgrB_female_peak.markers_14_Ex %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_14_Ex <- rbind(BgrB_male_peak.markers_14_Ex, BgrB_female_peak.markers_14_Ex)

saveRDS(BgrB_sex_peak.markers_14_Ex, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_14_Ex.rds')

BgrB_male_peak.markers_15_In.Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "15_In/Ex_MALE",
  ident.2 = "15_In/Ex_FEMALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_male_peak.markers_15_In.Ex$pct.diff <- ((abs(BgrB_male_peak.markers_15_In.Ex$pct.1 - BgrB_male_peak.markers_15_In.Ex$pct.2)) / ((BgrB_male_peak.markers_15_In.Ex$pct.1 + BgrB_male_peak.markers_15_In.Ex$pct.2)/2)) * 100
BgrB_male_peak.markers_15_In.Ex <- BgrB_male_peak.markers_15_In.Ex[order(BgrB_male_peak.markers_15_In.Ex$pct.diff, decreasing = T),]


BgrB_male_peak.markers_15_In.Ex$CellType <- "15_In/Ex"
BgrB_male_peak.markers_15_In.Ex$sex <- "MALE"
BgrB_male_peak.markers_15_In.Ex$Gene <- rownames(BgrB_male_peak.markers_15_In.Ex)
BgrB_male_peak.markers_15_In.Ex <- BgrB_male_peak.markers_15_In.Ex %>% dplyr::filter(p_val <= 0.05)


BgrB_female_peak.markers_15_In.Ex <- FindMarkers(
  object = tel_final,
  ident.1 = "15_In/Ex_FEMALE",
  ident.2 = "15_In/Ex_MALE",
  test.use = 'LR',
  latent.vars = c('nCount_peaks'),
  features = peak.feats, 
  assay = 'peaks', only.pos = T
)

BgrB_female_peak.markers_15_In.Ex$pct.diff <- ((abs(BgrB_female_peak.markers_15_In.Ex$pct.1 - BgrB_female_peak.markers_15_In.Ex$pct.2)) / ((BgrB_female_peak.markers_15_In.Ex$pct.1 + BgrB_female_peak.markers_15_In.Ex$pct.2)/2)) * 100
BgrB_female_peak.markers_15_In.Ex <- BgrB_female_peak.markers_15_In.Ex[order(BgrB_female_peak.markers_15_In.Ex$pct.diff, decreasing = T),]


BgrB_female_peak.markers_15_In.Ex$CellType <- "15_In/Ex"
BgrB_female_peak.markers_15_In.Ex$sex <- "FEMALE"
BgrB_female_peak.markers_15_In.Ex$Gene <- rownames(BgrB_female_peak.markers_15_In.Ex)
BgrB_female_peak.markers_15_In.Ex <- BgrB_female_peak.markers_15_In.Ex %>% dplyr::filter(p_val <= 0.05)

BgrB_sex_peak.markers_15_In.Ex <- rbind(BgrB_male_peak.markers_15_In.Ex, BgrB_female_peak.markers_15_In.Ex)

saveRDS(BgrB_sex_peak.markers_15_In.Ex, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_15_In.Ex.rds')



BgrB_sex_peak.markers_primary_celltype <- rbind(BgrB_sex_peak.markers_1_Astro.MG, BgrB_sex_peak.markers_2_OPC_Oligo, BgrB_sex_peak.markers_3_Peri, BgrB_sex_peak.markers_4_In, BgrB_sex_peak.markers_5_In, BgrB_sex_peak.markers_6_In, BgrB_sex_peak.markers_7_In, BgrB_sex_peak.markers_8_Ex, BgrB_sex_peak.markers_9_Ex, BgrB_sex_peak.markers_10_Ex, BgrB_sex_peak.markers_11_Ex, BgrB_sex_peak.markers_12_Ex, BgrB_sex_peak.markers_13_Ex, BgrB_sex_peak.markers_14_Ex, BgrB_sex_peak.markers_15_In.Ex)

saveRDS(BgrB_sex_peak.markers_primary_celltype, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_primary_celltype.rds')

BgrB_sex_peak.markers_primary_celltype_v2 <- BgrB_sex_peak.markers_primary_celltype %>% dplyr::filter(p_val_adj <= 0.05)

saveRDS(BgrB_sex_peak.markers_primary_celltype_v2, '/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/BgrB_sex_peak.markers_primary.celltype_significant.rds')



##############

BgrB_sex_peak.markers_primary.celltype_significant <- readRDS('/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/tests_kt/BgrB_sex_peak.markers_primary.celltype_significant.rds')
BgrB_sex_peak.markers_primary.celltype <- readRDS('/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/tests_kt/BgrB_sex_peak.markers_primary_celltype.rds')
tel_final <- readRDS('/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/FINAL_multiome_rds/Tel_OB/tel.ob_integrated_atac.gex_FINAL_04.17.24_KL.rds')
# tel <- readRDS('/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/MultiomeProjects/BgrB_multiome/tests_kt/tel.ob_corrected_integration_4.23.24.rds')
degs <- readRDS('/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/Katie/BgrB_atac_gex_multiome/FINAL_multiome_files/Final_objects/BgrB_Tel_OB/degs_dars/degs_by.secondary.celltype_per.male.female&bhve.ctrl.rds')
mz_gene_list <- read_csv('/Users/kathrynleatherbury/Dropbox (GaTech)/BioSci-Streelman/Katie/BgrB_atac_gex_multiome/gene_info_3.csv')
tel_final$celltype.sex <- paste(tel_final$primary_celltypes, tel_final$SEX, sep = "_")
tel_final$celltype.cond <- paste(tel_final$celltype.sex, tel_final$COND, sep = "_")
Idents(tel_final) <- 'celltype.cond'
tel_final$celltype.sample <- paste(tel_final$primary_celltypes, tel_final$SAMPLE, sep = "_")
Idents(tel_final) <- 'celltype.sample'

tel.ob <- LinkPeaks(object = tel.ob, peak.assay = "peaks", expression.assay = "SCT", genes.use = c('fabp7', 'fosb'))

CoveragePlot(
  object = tel_final,
  region = "fosb",
  features = "fosb",
  expression.assay = "SCT",
  extend.upstream = 5000,
  extend.downstream = 5000,
  annotation=T, 
  assay = 'peaks', group.by = 'celltype.sample', sep = c('-', '-'), 
  scale.factor = 1e6, 
  ymax = 75, idents = c('7_In_FEMALE_BHVE', '7_In_MALE_BHVE', '13_Ex_FEMALE_BHVE', '13_Ex_MALE_BHVE', '2_OPC/Oligo_FEMALE_BHVE', '2_OPC/Oligo_MALE_BHVE', '3_Peri_FEMALE_BHVE', '3_Peri_MALE_BHVE')
)

CoveragePlot(
  object = tel_final,
  region = "ptgfrn",
  features = "ptgfrn",
  expression.assay = "SCT",
  extend.upstream = 5000,
  extend.downstream = 1000,
  annotation=T, 
  assay = 'peaks', group.by = 'celltype.sample', sep = c('-', '-'), 
  scale.factor = 1e6, heights = c(7,1),
  ymax = 75, idents = c('1_Astro/MG_MALE_BHVE','1_Astro/MG_MALE_CTRL', '2_OPC/Oligo_MALE_BHVE', '2_OPC/Oligo_MALE_CTRL', '3_Peri_MALE_BHVE', '3_Peri_MALE_CTRL', '4_In_MALE_BHVE','4_In_MALE_CTRL','5_In_MALE_BHVE','5_In_MALE_CTRL', '6_In_MALE_BHVE', '6_In_MALE_CTRL', '7_In_MALE_BHVE','7_In_MALE_CTRL','8_Ex_MALE_BHVE','8_Ex_MALE_CTRL','9_Ex_MALE_BHVE', '9_Ex_MALE_CTRL','10_Ex_MALE_BHVE','10_Ex_MALE_CTRL', '11_Ex_MALE_BHVE','11_Ex_MALE_CTRL', '12_Ex_MALE_BHVE', '12_Ex_MALE_CTRL','13_Ex_MALE_BHVE','13_Ex_MALE_CTRL','14_Ex_MALE_BHVE','14_Ex_MALE_CTRL','15_In/Ex_MALE_BHVE','15_In/Ex_MALE_CTRL')
)

CoveragePlot(
  object = tel_final,
  region = "ptgfr",
  features = "ptgfr",
  expression.assay = "SCT",
  extend.upstream = 1000,
  extend.downstream = 1000,
  annotation=T, 
  assay = 'peaks', group.by = 'celltype.sample', sep = c('-', '-'), 
  scale.factor = 1e6, heights = c(7,1),
  ymax = 25, idents = c('1_Astro/MG_FEMALE_BHVE','1_Astro/MG_FEMALE_CTRL', '2_OPC/Oligo_FEMALE_BHVE', '2_OPC/Oligo_FEMALE_CTRL', '3_Peri_FEMALE_BHVE', '3_Peri_FEMALE_CTRL', '4_In_FEMALE_BHVE','4_In_FEMALE_CTRL','5_In_FEMALE_BHVE','5_In_FEMALE_CTRL', '6_In_FEMALE_BHVE', '6_In_FEMALE_CTRL', '7_In_FEMALE_BHVE','7_In_FEMALE_CTRL','8_Ex_FEMALE_BHVE','8_Ex_FEMALE_CTRL','9_Ex_FEMALE_BHVE', '9_Ex_FEMALE_CTRL','10_Ex_FEMALE_BHVE','10_Ex_FEMALE_CTRL', '11_Ex_FEMALE_BHVE','11_Ex_FEMALE_CTRL', '12_Ex_FEMALE_BHVE', '12_Ex_FEMALE_CTRL','13_Ex_FEMALE_BHVE','13_Ex_FEMALE_CTRL','14_Ex_FEMALE_BHVE','14_Ex_FEMALE_CTRL','15_In/Ex_FEMALE_BHVE','15_In/Ex_FEMALE_CTRL')
)


CoveragePlot(
  object = tel_final,
  region = "fosb",
  features = "fosb",
  expression.assay = "SCT",
  extend.upstream = 5000,
  extend.downstream = 5000,
  annotation=T, 
  assay = 'peaks', group.by = 'celltype.sample', split.by = 'primary_celltypes'
)
c('1_Astro/MG_MALE_BHVE ','1_Astro/MG_MALE_CTRL', '2_OPC/Oligo_MALE_BHVE', '2_OPC/Oligo_MALE_CTRL', '3_Peri_MALE_BHVE', '3_Peri_MALE_CTRL', '4_In_MALE_BHVE','4_In_MALE_CTRL','5_In_MALE_BHVE','5_In_MALE_CTRL', '6_In_MALE_BHVE', '6_In_MALE_CTRL', '7_In_MALE_BHVE','7_In_MALE_CTRL','8_Ex_MALE_BHVE','8_Ex_MALE_CTRL','9_Ex_MALE_BHVE', '9_Ex_MALE_CTRL','10_Ex_MALE_BHVE','10_Ex_MALE_CTRL', '11_Ex_MALE_BHVE','11_Ex_MALE_CTRL', '12_Ex_MALE_BHVE', '12_Ex_MALE_CTRL','13_Ex_MALE_BHVE','13_Ex_MALE_CTRL','14_Ex_MALE_BHVE','14_Ex_MALE_CTRL','15_In/Ex_MALE_BHVE','15_In/Ex_MALE_CTRL')


p1 <- CoveragePlot(
  object = tel_final,
  region = "NC_036797.1-11998585-11999537",
  extend.upstream = 40000,
  extend.downstream = 40000, assay = 'peaks', 
  annotation = F, 
  peaks = F,
  group.by = 'celltype.sex', 
  sep = c('-', '-'), 
  scale.factor = 1e7, 
  ymax = 100,
  heights = c(7,1), split.by = 'SEX') + scale_fill_manual(values = c('seagreen', 'seagreen2', 'seagreen', 'seagreen2', 'seagreen2', 'seagreen2', 'seagreen2', 'seagreen', 'seagreen', 'seagreen', 'seagreen', 'seagreen', 'seagreen', 'seagreen', 'seagreen2'))

p2 <- CoveragePlot(
  object = tel_final,
  region = "NC_036797.1-11998585-11999537",
  extend.upstream = 40000,
  extend.downstream = 40000, assay = 'peaks', 
  annotation = F, 
  peaks = F,
  group.by = 'SEX', 
  sep = c('-', '-'), 
  scale.factor = 1e7, 
  ymax = 100, split.by = 'primary_celltypes') + scale_fill_manual(values = c('violetred1', 'turquoise3'))

ggarrange(p2, p1, heights = c(1.3, 5),
           ncol = 1, nrow = 2, align = "v")


wrap_plots(p1, p2, ncol = 1)

p1 <- CoveragePlot(
  object = tel_final,
  region = "NW_020192755.1-106158-107010",
  extend.upstream = 40000,
  extend.downstream = 40000, assay = 'peaks', 
  annotation = F, 
  peaks = F,
  group.by = 'primary_celltypes', 
  sep = c('-', '-'), 
  scale.factor = 1e7, 
  ymax = 125,
  heights = c(7,1)) + scale_fill_manual(values = c('red', 'orange1', 'gold1', 'green2', 'green2', 'green2', 'green2', 'turquoise3', 'turquoise3', 'turquoise3', 'turquoise3', 'turquoise3', 'turquoise3', 'turquoise3', 'blue'))

p2 <- CoveragePlot(
  object = tel_final,
  region = "NW_020192755.1-106158-107010",
  extend.upstream = 40000,
  extend.downstream = 40000, assay = 'peaks', 
  annotation = F, 
  peaks = F,
  group.by = 'SEX', 
  sep = c('-', '-'), 
  scale.factor = 1e7, 
  ymax = 125) + scale_fill_manual(values = c('gray30', 'gray70'))

ggarrange(p2, p1, heights = c(1.5, 5),
          ncol = 1, nrow = 2, align = "v")


# spidr gene is important for ovarian function

FeaturePlot(tel_final, features = c('spidr'), reduction = 'sctUMAP', order = T, split.by = 'SAMPLE')
FeaturePlot_scCustom(tel_final, features = c('spidr'), reduction = 'sctUMAP', order = T, label = F, split.by = 'SAMPLE', combine = F, colors_use = viridis(option = 'C', n = 8, direction = -1))

p1 <- DotPlot_scCustom(tel_final, group.by = 'celltype.cond', features = 'spidr', scale = T, flip_axes = F, assay = 'SCT', idents = c('9_Ex_MALE_BHVE', '9_Ex_FEMALE_BHVE', '6_In_MALE_BHVE', '6_In_FEMALE_BHVE', '4_In_MALE_BHVE', '4_In_FEMALE_BHVE'), colors_use = wes_palette("Zissou1", 100, type = "continuous"))

p1$data$id <- factor(p1$data$id, levels = c('9_Ex_MALE_BHVE', '9_Ex_FEMALE_BHVE', '6_In_MALE_BHVE', '6_In_FEMALE_BHVE', '4_In_MALE_BHVE', '4_In_FEMALE_BHVE'))
p1


CoveragePlot(
  object = tel_final,
  region = "NW_020192755.1-106158-107010", ymax = 'q98', 
  extend.upstream = 40000,
  extend.downstream = 40000, assay = 'peaks', annotation = T, group.by = 'primary_celltypes', sep = c('-', '-'), scale.factor = 1e7, ymax = 362.6)

CoveragePlot(
  object = tel_final,
  region = "NC_036793.1-26139025-261401647", 
  extend.upstream = 40000,
  extend.downstream = 40000, assay = 'peaks', annotation = T, group.by = 'SAMPLE', sep = c('-', '-'))



p2 <- CoveragePlot(
  object = tel_final,
  region = "NC_036797.1-11998585-11999537",
  extend.upstream = 40000,
  extend.downstream = 40000, assay = 'peaks', 
  annotation = F, 
  peaks = F,
  group.by = 'celltype.cond', 
  sep = c('-', '-'), 
  scale.factor = 1e7, 
  ymax = 125, idents = c('4_In_FEMALE_BHVE', '4_In_MALE_BHVE', '6_In_FEMALE_BHVE', '6_In_MALE_BHVE', '9_Ex_FEMALE_BHVE', '9_Ex_MALE_BHVE', '8_Ex_FEMALE_BHVE', '8_Ex_MALE_BHVE')) + scale_fill_manual(values = c('violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise', 'violetred1', 'darkturquoise'))

p2$data$group <- factor(p2$data$group, levels = c(c('4_In_FEMALE_BHVE', '4_In_MALE_BHVE', '6_In_FEMALE_BHVE', '6_In_MALE_BHVE', '9_Ex_FEMALE_BHVE', '9_Ex_MALE_BHVE', '8_Ex_FEMALE_BHVE', '8_Ex_MALE_BHVE')))
p2

CoveragePlot(
  object = tel_final,
  region = "NC_036793.1-26139025-261401647",
  extend.upstream = 40000,
  extend.downstream = 40000, assay = 'peaks', 
  annotation = F, 
  peaks = F,
  group.by = 'SEX', 
  sep = c('-', '-'))




##################



p1 <- CoveragePlot(
  object = tel_final,
  region = "NC_036797.1-11998585-11999537",
  extend.upstream = 20000,
  extend.downstream = 20000, assay = 'peaks', 
  annotation = T, 
  peaks = F,
  group.by = 'primary_celltypes', 
  scale.factor = 1e7, 
  ymax = 200,
  heights = c(7,1)) + scale_fill_manual(values = c('red', 'orange1', 'gold1', 'green2', 'green2', 'green2', 'green2', 'turquoise3', 'turquoise3', 'turquoise3', 'turquoise3', 'turquoise3', 'turquoise3', 'turquoise3', 'blue'))

p2 <- CoveragePlot(
  object = tel_final,
  region = "NC_036797.1-11998585-11999537",
  extend.upstream = 20000,
  extend.downstream = 20000, assay = 'peaks', 
  annotation = F, 
  peaks = F,
  group.by = 'SEX', 
  scale.factor = 1e7, 
  ymax = 200) + scale_fill_manual(values = c('magenta', 'turquoise3'))

ggarrange(p2, p1, heights = c(1.5, 5),
          ncol = 1, nrow = 2, align = "v")

DimPlot(rock.sand, reduction = 'sctUMAP', group.by = 'primary_celltypes', split.by = 'GROUP', label = T)
FeaturePlot(rock.sand, reduction = 'sctUMAP', features = 'slc17a6', order = T)


p3 <- FeaturePlot_scCustom(rock.sand, features = c('gad2'), reduction = 'sctUMAP', order = T, label = F, split.by = 'GROUP', combine = F, colors_use = viridis(option = 'D', n = 8, direction = -1))

LabelClusters(p3[[1]], id = "ident", clusters = celltypes, fontface = "bold", color = "black", position = 'median', repel = F) + NoAxes()
LabelClusters(p3[[2]], id = "ident", clusters = celltypes, fontface = "bold", color = "black", position = 'median', repel = F) + NoAxes()


