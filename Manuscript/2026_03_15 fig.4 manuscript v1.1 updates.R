library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)
library(CytoTRACE)
library(lme4)

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_6 = FindSubCluster(obj, 6, graph.name='harmony.wsnn')
Idents(sub_6) <- 'sub.cluster'
sub_6 = subset(sub_6, final_clusters ==6)
sub_6$sub.cluster = factor(sub_6$sub.cluster, levels = c(paste0('6_',0:3)))
sub_6$Status = factor(sub_6$Status, levels = c('NRM','M',"D",'E','NF','F'))

### proportion of 6 #####
total_cells = sub_6@meta.data%>%
  group_by(Status, individual)%>%
  subset(Status !='NRM')%>%
  summarize(total_cells = n())

######degs
degs = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_6 = subset(degs, cluster == 6)

p_value_annotation <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return(paste0("p = ", round(p, 3)))
  }
}

obj@meta.data$Phase = obj@meta.data$Status
obj@meta.data$Phase = ifelse(obj@meta.data$Phase == 'D', 'I', obj@meta.data$Phase)
obj@meta.data$Phase = ifelse(obj@meta.data$Phase == 'E', 'LI', obj@meta.data$Phase)
obj@meta.data$Phase = factor(obj@meta.data$Phase, levels = c('NRM','M','I','LI','NF','F'))

make_summary = function(dat, numerator_col, denominator_col) {
  dat %>%
    filter(Phase != "NF") %>%
    group_by(Phase) %>%
    summarise(
      mean_prop = mean(.data[[numerator_col]] / .data[[denominator_col]], na.rm = TRUE),
      se = sd(.data[[numerator_col]] / .data[[denominator_col]], na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
}

#### GnRH + proportion ####
sub_6$GnRH = ifelse(sub_6@assays$RNA$data['LOC111571064',] > 0, T, F)

gnrh = sub_6@meta.data %>%
  group_by(individual, Status) %>%
  summarize(n_GnRH = sum(GnRH == T))

tot = sub_6@meta.data %>%
  group_by(individual, Status) %>%
  summarize(total_cells = n())

tog = gnrh %>%
  right_join(tot, by = c('individual'))

tog_6_gnrh = subset(tog, Status.x != 'NRM')

gnrh_glm = glm(cbind(tog_6_gnrh$n_GnRH, tog_6_gnrh$total_cells - tog_6_gnrh$n_GnRH) ~ Status.x,
               data = tog_6_gnrh,
               family = binomial('logit'))

anova(gnrh_glm, test = 'Chisq')

tog_6_gnrh$Phase = as.character(tog_6_gnrh$Status.x)
tog_6_gnrh$Phase = ifelse(tog_6_gnrh$Phase == 'D', 'I', tog_6_gnrh$Phase)
tog_6_gnrh$Phase = ifelse(tog_6_gnrh$Phase == 'E', 'LI', tog_6_gnrh$Phase)
tog_6_gnrh$Phase = factor(tog_6_gnrh$Phase, levels = c('M','I','LI','NF','F'))

gnrh_pairs = pairs(emmeans(gnrh_glm, 'Status.x'), adjust = 'none')

tog_6_gnrh$prop = tog_6_gnrh$n_GnRH / tog_6_gnrh$total_cells

summary_gnrh = make_summary(tog_6_gnrh, 'n_GnRH', 'total_cells')

gnrh_prop = ggplot(tog_6_gnrh, aes(x = Phase, y = n_GnRH / total_cells)) +
  geom_bar(data = summary_gnrh,
           aes(x = Phase, y = mean_prop, fill = Phase),
           stat = 'identity', inherit.aes = FALSE) +
  geom_errorbar(data = summary_gnrh,
                aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.4, inherit.aes = FALSE) +
  geom_point(size = 1) +
  theme_classic() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = 'Phase', y = '% GnRH of 6') +
  theme(legend.position = 'none') +
  geom_signif(xmin = c(1), xmax = c(1.9),
              y_position = c(max(tog_6_gnrh$prop * 1.1)),
              annotation = c('p = 0.514'), color = "black",
              tip_length = c(0, 0), textsize = 3) +
  geom_signif(xmin = c(1), xmax = c(5),
              y_position = c(max(tog_6_gnrh$prop * 1.3)),
              annotation = c('p = 0.127'), color = "black",
              tip_length = c(0, 0), textsize = 3) +
  geom_signif(xmin = c(2.1), xmax = c(5),
              y_position = c(max(tog_6_gnrh$prop * 1.1)),
              annotation = c('p = 0.372'), color = "black",
              tip_length = c(0, 0), textsize = 3)

 # ggsave(plot = gnrh_prop,
  #     file = 'gnrh_prop_6.svg',
   #    device = "svg",
    #   units = "in",
     #  width = 2,
      # height = 2,
       #path = "Manuscript/Plots/Manuscript v1.1/Fig.4")


### cckb proportion #####
sub_6$cckb = ifelse(sub_6@assays$RNA$data['cckb',] > 0, T, F)

cckb = sub_6@meta.data %>%
  group_by(individual, Status) %>%
  summarize(n_cckb = sum(cckb == T))

tot = sub_6@meta.data %>%
  group_by(individual, Status) %>%
  summarize(total_cells = n())

tog = cckb %>%
  right_join(tot, by = c('individual'))

tog_6_cckb = subset(tog, Status.x != c('NRM'))

cckb_glm = glm(cbind(tog_6_cckb$n_cckb, tog_6_cckb$total_cells - tog_6_cckb$n_cckb) ~ Status.x,
               data = tog_6_cckb,
               family = binomial('logit'))

anova(cckb_glm, test = 'Chisq')

tog_6_cckb$Phase = as.character(tog_6_cckb$Status.x)
tog_6_cckb$Phase = ifelse(tog_6_cckb$Phase == 'D', 'I', tog_6_cckb$Phase)
tog_6_cckb$Phase = ifelse(tog_6_cckb$Phase == 'E', 'LI', tog_6_cckb$Phase)
tog_6_cckb$Phase = factor(tog_6_cckb$Phase, levels = c('M','I','LI','NF','F'))

cckb_pairs = pairs(emmeans(cckb_glm, 'Status.x'), adjust = 'none')

tog_6_cckb$prop = tog_6_cckb$n_cckb / tog_6_cckb$total_cells

summary_cckb = make_summary(tog_6_cckb, 'n_cckb', 'total_cells')

cckb_prop = ggplot(tog_6_cckb, aes(x = Phase, y = n_cckb / total_cells)) +
  geom_bar(data = summary_cckb,
           aes(x = Phase, y = mean_prop, fill = Phase),
           stat = 'identity', inherit.aes = FALSE) +
  geom_errorbar(data = summary_cckb,
                aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.4, inherit.aes = FALSE) +
  geom_point(size = 1) +
  theme_classic() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = 'Phase', y = '% cckb of 6') +
  theme(legend.position = 'none') +
  geom_signif(xmin = c(1), xmax = c(1.9),
              y_position = c(max(tog_6_cckb$prop * 1.1)),
              annotation = c('p = 0.629'), color = "black",
              tip_length = c(0, 0), textsize = 3) +
  geom_signif(xmin = c(1), xmax = c(5),
              y_position = c(max(tog_6_cckb$prop * 1.3)),
              annotation = c('*'), color = "black",
              tip_length = c(0, 0), textsize = 3) +
  geom_signif(xmin = c(2.1), xmax = c(5),
              y_position = c(max(tog_6_cckb$prop * 1.1)),
              annotation = c('*'), color = "black",
              tip_length = c(0, 0), textsize = 3)
cckb_prop

  #ggsave(plot = cckb_prop,
   #    file = 'gnrh_prop_6.svg',
    #   device = "svg",
     #  units = "in",
      # width = 2,
       #height = 2,
       #path = "Manuscript/Plots/Manuscript v1.1/Fig.4")

### cyto by degs ####
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')

degs_plasticity = c(
  'LOC111588913', # nxph1
  'cntn4',
  'LOC111567620', #cntnapl4
  'pcdh10b',
  'sdc2',
  'LOC111585095', # dsel
  'bcan'
)

degs_neurotransmission = c(
  'gabra3',
  'nmbr',
  'npy7r',
  'pgr',
  'LOC111577100', # cd87
  'grm5a',
  'LOC111568069', #pgrl
  'tacr3a',
  'drd3'
)


sub_6_plasticity = sub_6@assays$RNA$data[degs_plasticity,]
sub_6_nt = sub_6@assays$RNA$data[degs_neurotransmission,]

plasticity_positive_cells = colSums(sub_6_plasticity)>0
nt_positive_cells = colSums(sub_6_nt)>0


sub_6$plasticity = ifelse(plasticity_positive_cells, T, F)
sub_6$nt = ifelse(nt_positive_cells, T, F)

FeaturePlot(sub_6, 'plasticity')
FeaturePlot(sub_6, 'nt')

DimPlot(sub_6)

DotPlot(sub_6, 'plasticity')+
  coord_flip()

DotPlot(sub_6, 'nt')+ # strong 6_0
  coord_flip()

cyto = CytoTRACE(sub_6@assays$RNA$data%>%as.matrix())

sub_6$cyto = cyto$CytoTRACE
FeaturePlot(sub_6, 'cyto')

sub_6@meta.data$plasticity = factor(sub_6@meta.data$plasticity, levels = c(T,F))

cyto_ecm_dat = sub_6@meta.data %>%
  group_by(individual, Status, plasticity) %>%
  summarize(mean_cyto = mean(cyto)) %>%
  subset(Status != 'NRM')
cyto_ecm_dat$Status =as.character(cyto_ecm_dat$Status)
cyto_ecm_dat$Status = factor(cyto_ecm_dat$Status, levels = c('M', 'D','E','NF','F'))

cyto_ecm = ggplot(cyto_ecm_dat, aes(x = Status, y = mean_cyto, fill = Status)) +
  geom_boxplot(data = subset(cyto_ecm_dat, Status != 'NF'),
               aes(x = Status, y = mean_cyto, fill = Status),
               outlier.shape = NA) +
  geom_point(position = position_dodge(.75)) +
  theme_minimal() +
  facet_wrap(~plasticity) +
  theme(legend.position = 'none') +
  labs(y = 'Mean CytoTRACE', x = 'Phase') +
  scale_x_discrete(drop = FALSE) +
  ylim(0.2, .9)
cyto_ecm

ggsave(plot = cyto_ecm,
       file =  'cyto_ecm.svg',
       device = "svg",
       units = "in",
       width = 3,
       height = 3,
       path = "Manuscript/Plots/Manuscript v1.1/Fig.4")

plas_cyto_mod = lmer(cyto ~ (Status * plasticity) + nCount_RNA + (1|individual),
                     data = subset(sub_6@meta.data, Status != 'NRM'))
car::Anova(plas_cyto_mod, 3)
pairs(emmeans(plas_cyto_mod, 'Status', by = 'plasticity'), adjust = 'none')

nt_cyto_mod = lmer(cyto ~ nCount_RNA+(Status * nt) + (1|individual),
                     data = subset(sub_6@meta.data, Status != 'NRM'))
car::Anova(nt_cyto_mod, 3)
pairs(emmeans(nt_cyto_mod, 'Status', by = 'nt'), adjust = 'none')

sub_6@meta.data$nt = factor(sub_6@meta.data$nt, levels = c(T,F))

nt_cyto_dat = sub_6@meta.data %>%
  group_by(individual, Status, nt) %>%
  summarize(mean_cyto = mean(cyto)) %>%
  subset(Status != 'NRM')

nt_cyto_dat$Status =as.character(nt_cyto_dat$Status)
nt_cyto_dat$Status = factor(nt_cyto_dat$Status, levels = c('M', 'D','E','NF','F'))

nt_cyto = ggplot(nt_cyto_dat, aes(x = Status, y = mean_cyto, fill = Status)) +
  geom_boxplot(data = subset(nt_cyto_dat, Status != 'NF'),
               aes(x = Status, y = mean_cyto, fill = Status),
               outlier.shape = NA) +
  geom_point(position = position_dodge(.75)) +
  theme_minimal() +
  facet_wrap(~nt) +
  theme(legend.position = 'none') +
  labs(y = 'Mean CytoTRACE', x = 'Phase') +
  scale_x_discrete(drop = FALSE)+
  ylim(0.3, 0.7)
nt_cyto

#ggsave(plot = nt_cyto,
 #      file =  'nt_cyto.svg',
  #     device = "svg",
   #    units = "in",
    #   width = 3,
     #  height = 3,
      # path = "Manuscript/Plots/Manuscript v1.1/Fig.4")
######cckb cyto #####

sub_6$cckb = ifelse(sub_6@assays$RNA$data['cckb',]>0, T, F)
sub_6@meta.data$cckb = factor(sub_6@meta.data$cckb, levels = c(T,F))

cckb_cyto_dat = sub_6@meta.data %>%
  group_by(individual, Status, cckb) %>%
  summarize(mean_cyto = mean(cyto)) %>%
  subset(Status != 'NRM')

cckb_cyto_dat$Status =as.character(cckb_cyto_dat$Status)
cckb_cyto_dat$Status = factor(cckb_cyto_dat$Status, levels = c('M', 'D','E','NF','F'))

cckb_cyto = ggplot(cckb_cyto_dat, aes(x = Status, y = mean_cyto, fill = Status)) +
  geom_boxplot(data = subset(cckb_cyto_dat, Status != 'NF'),
               aes(x = Status, y = mean_cyto, fill = Status),
               outlier.shape = NA) +
  geom_point(position = position_dodge(.75)) +
  theme_minimal() +
  facet_wrap(~cckb) +
  theme(legend.position = 'none') +
  labs(y = 'Mean CytoTRACE', x = 'Phase') +
  scale_x_discrete(drop = FALSE)+
  ylim(0.3, 0.8)
cckb_cyto

#ggsave(plot = cckb_cyto,
 #      file =  'cckb_cyto.svg',
  #     device = "svg",
   #    units = "in",
    #   width = 3,
     #  height = 3,
      # path = "Manuscript/Plots/Manuscript v1.1/Fig.4")

cck_cyto_mod = lmer(cyto ~ (Status * cckb) + nCount_RNA + (1|individual),
                     data = subset(sub_6@meta.data, Status != 'NRM'))

car::Anova(cck_cyto_mod, 3) ##### lets fucking gooooooo
pairs(emmeans(cck_cyto_mod, 'Status', by ='cckb'), adjust = 'none')

#again mmodel_ecm_1#again more immature

sub_6$gnrh1 = ifelse(sub_6@assays$RNA$data['LOC111571064',]>0, T, F)
ggplot(sub_6@meta.data%>%
         group_by(individual, Status, gnrh1)%>%
         summarize(mean_cyto = mean(cyto))%>%
         subset(Status!='NRM'), aes(x = Status, y = mean_cyto, fill= Status))+
  geom_boxplot()+
  geom_point(position = position_dodge(.75))+
  theme_minimal()+
  facet_wrap(~gnrh1)+
  theme(legend.position = 'none')+
  labs(y = 'Mean CytoTRACE',x = 'Phase' )

# in this case it is significantly lower, though with the high scores of some Is you coulr
# argue that there are also new cells with them
mode = lmer(cyto~nCount_RNA+(gnrh1*Status)+(1|individual), data = sub_6@meta.data)
summary(mode)
car::Anova(mode,3)
#AHAH GNRH1 cells are NOT new

#### proportions ####
proportion_test = function(gene){
  
  expr = sub_6@assays$RNA$data[gene,] > 0
  
  dat = sub_6@meta.data
  dat$expr = expr
  
  agg = dat %>%
    dplyr::filter(Status != "NRM") %>%
    dplyr::group_by(individual, Status) %>%
    dplyr::summarize(
      success = sum(expr),
      total = dplyr::n(),
      .groups = "drop"
    )
  
  mod = glm(
    cbind(success, total - success) ~ Status,
    data = agg,
    family = binomial("logit")
  )
  
  return(mod)
}

proportion_plot = function(gene){
  
  expr = sub_6@assays$RNA$data[gene,] > 0
  
  dat = sub_6@meta.data
  dat$expr = expr
  
  agg = dat %>%
    dplyr::filter(Status != "NRM") %>%
    dplyr::group_by(individual, Status) %>%
    dplyr::summarize(
      success = sum(expr),
      total = dplyr::n(),
      .groups = "drop"
    )
  
  agg$Phase = as.character(agg$Status)
  agg$Phase = ifelse(agg$Phase == 'D', 'I', agg$Phase)
  agg$Phase = ifelse(agg$Phase == 'E', 'LI', agg$Phase)
  agg$Phase = factor(agg$Phase, levels = c('M', 'I', 'LI', 'NF', 'F'))
  
  agg_summary = agg %>%
    dplyr::filter(Phase != "NF") %>%
    dplyr::group_by(Phase) %>%
    dplyr::summarize(
      mean_prop = mean(success / total, na.rm = TRUE),
      se = sd(success / total, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )
  
  p = ggplot(agg, aes(x = Phase, y = success / total)) +
    geom_bar(data = agg_summary,
             aes(x = Phase, y = mean_prop, fill = Phase),
             stat = 'identity', inherit.aes = FALSE) +
    geom_errorbar(data = agg_summary,
                  aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                  width = 0.4, inherit.aes = FALSE) +
    geom_point(size = 1) +
    theme_classic() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = 'Phase', y = paste0('% ', gene, '+ of 6')) +
    theme(legend.position = 'none')
  
  return(p)
}

proportion_test('tacr3a')
proportion_plot('tacr3a')

proportion_test('cckb')
proportion_plot('cckb')


cckb_prop_test =proportion_test('cckb')
anova(cckb_prop_test, test= 'Chisq')
pairs(emmeans(cckb_prop_test, 'Status'), adjust = 'none')

cckb_prop = proportion_plot('cckb')+
    geom_signif(xmin = c(1), xmax = c(1.9),
              y_position = c(0.3),
              annotation = c('p = 0.629'), color = "black",
              tip_length = c(0, 0), textsize = 3) +
      geom_signif(xmin = c(2.1), xmax = c(5),
              y_position = c(0.3),
              annotation = c('*'), color = "black",
              tip_length = c(0, 0), textsize = 3) +
      geom_signif(xmin = c(1), xmax = c(5),
              y_position = c(0.35),
              annotation = c('*'), color = "black",
              tip_length = c(0, 0), textsize = 3) 



#ggsave(plot = cckb_prop,
 #      file =  'cckb_prop.svg',
  #     device = "svg",
   #    units = "in",
    #   width = 3,
     #  height = 3,
      # path = "Manuscript/Plots/Manuscript v1.1/Fig.4")

gnrh1_test =proportion_test('LOC111571064')
anova(gnrh1_test, test ='Chisq')
pairs(emmeans(gnrh1_test, 'Status'), adjust ='none')

gnrh1_prop = proportion_plot('LOC111571064')+
      geom_signif(xmin = c(1), xmax = c(1.9),
              y_position = c(0.1),
              annotation = c('p = 0.514'), color = "black",
              tip_length = c(0, 0), textsize = 3) +
      geom_signif(xmin = c(2.1), xmax = c(5),
              y_position = c(0.1),
              annotation = c('p = 0.372'), color = "black",
              tip_length = c(0, 0), textsize = 3) +
      geom_signif(xmin = c(1), xmax = c(5),
              y_position = c(0.12),
              annotation = c('p = 0.127'), color = "black",
              tip_length = c(0, 0), textsize = 3) 


ggsave(plot = gnrh1_prop,
       file =  'gnrh1_prop.svg',
       device = "svg",
       units = "in",
       width = 3,
       height = 3,
       path = "Manuscript/Plots/Manuscript v1.1/Fig.4")






