library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(lme4)
library(car)
library(emmeans)
library(CytoTRACE)

# ---- data prep ----
obj   = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
sub_6 = subset(obj, res0.8_50nn_40PC_45LSI == 6)
sub_6$Status = factor(sub_6$Status, levels = c('NRM', 'M', 'D', 'E', 'NF', 'F'))

sub_6@meta.data$Phase = as.character(sub_6@meta.data$Status)
sub_6@meta.data$Phase = ifelse(sub_6@meta.data$Phase == 'D',  'I',  sub_6@meta.data$Phase)
sub_6@meta.data$Phase = ifelse(sub_6@meta.data$Phase == 'E',  'LI', sub_6@meta.data$Phase)
sub_6@meta.data$Phase = factor(sub_6@meta.data$Phase, levels = c('NRM', 'M', 'I', 'LI', 'NF', 'F'))

cyto      = CytoTRACE(sub_6@assays$RNA$data %>% as.matrix())
sub_6$cyto = cyto$CytoTRACE


# ---- model ----
temp = sub_6@meta.data%>%
  subset( Status != 'NRM')%>%
  group_by(individual,Phase)%>%
  summarize(mean_cyto = mean(cyto))

model  = lm(mean_cyto ~ Phase, data = temp)
av_    = anova(model, test ='Chisq')

pairs(emmeans(model, 'Phase'), adjust ='none')


