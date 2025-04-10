{
  library(parallel)
  library(clusterProfiler)
  library(blme)
  library(Seurat)
  library(tidyverse)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(SeuratObject)
  library(Signac)
  #library(CytoTRACE)
  # SCTRransform_mean_plot <- readRDS("R/Gabe/SCTRransform_mean_plot.rds")
  #mac.neg.bin <- readRDS(file = 'R/Gabe/mac.neg.bin.rds')
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  #mean_cell <- readRDS('R/Gabe/mean_cell.rds')
  library(openxlsx)
  #clown_go <- readRDS('R/Gabe/clown_go.rds')
  
}
multiome_object <-readRDS("~/Desktop/snRNA-seq R Files 122524/RNA Object.rds")
data <- read_excel("Reference/Complete Data Frame (Hormones, Behavior, Size, Gonads) GG.xlsx")

ggplot(data, aes(x = Log10_Volume, y = length_final_cm, color = Status))+
  geom_point()+
  geom_smooth(method = 'lm')

ggplot(data, aes(x = Log10_Volume, y = mass_final_cm, color = Status))+
  geom_point()+
  geom_smooth(method = 'lm')

mass_model <- lm(Log10_Volume~Status+mass_final_cm, data = data)


length_model <- lm(Log10_Volume~Status+length_final_cm, data = data)

both_model <- lm(Log10_Volume~Status+length_final_cm +mass_final_cm, data = data)

anova(mass_model,
      length_model,
      both_model)

summary(length_model) #adj r2 = 0.689 ### best barely
summary(mass_model) #adj r = 0.6717
summary(both_model) #adj r = 0.6718

library(ggeffects)

length_model_corrected <- ggeffect(length_model, c('Status'))%>%
  as.data.frame()
length_model_corrected$x = factor(length_model_corrected$x, levels = c('S','EP','NM','M','D','E','NF','F'))
 
 ggplot(length_model_corrected,aes(x = x, y = predicted))+
  geom_pointrange(aes(x = x, y = predicted, ymin = predicted -std.error, ymax = predicted + std.error))+
  labs(x = 'Sex',y = 'Predicted Gonad Volume Corrected for Length +/- SE')

 library(emmeans)
pairs(emmeans(length_model, 'Status'), adjust = 'none')
