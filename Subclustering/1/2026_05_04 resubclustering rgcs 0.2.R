library(Seurat)
library(ggplot2)
library(tidyverse)
library(emmeans)
library(ggsignif)

clown_go = readRDS('Functions/clown_go2')

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
obj = FindSubCluster(obj,
                     1, 
                     'harmony.wsnn',
                     resolution = 0.5,
                     subcluster.name = 'sub_res0.5')

obj = FindSubCluster(obj,
                     1, 
                     'harmony.wsnn',
                     resolution = 0.2,
                     subcluster.name = 'sub_res0.2')

obj = FindSubCluster(obj,
                     1, 
                     'harmony.wsnn',
                     resolution = 0.1,
                     subcluster.name = 'sub_res0.1')

DimPlot(obj, reduction = 'harmony_wnn.umap', group.by = 'sub_res0.2')

DimPlot(obj, reduction = 'harmony_wnn.umap', group.by = 'sub_res0.1')

sub_1 = subset(obj, final_clusters == 1)

marks_0.2 = FindAllMarkers(sub_1, group.by = 'sub_res0.2', only.pos = T)
marks_0.1 = FindAllMarkers(sub_1, group.by = 'sub_res0.1', only.pos = T)


# gonna try with res 0.2

Idents(sub_1) = 'sub_res0.2'

DotPlot(sub_1, 
        c('gfap',
          'sox2',
          'shha',
          'LOC111577263',
          'fgfr3',
          'elavl3',
          'LOC111576017',
          'sox4b',
          'pax6b'
          ))+
  coord_flip()

# 1_0 expresses a lot of txn factors
#1_1 solute carriers
# 1_2 different txn factors and a lot of aroB
# 1_3 gfap


#readin visium data
markers = read.csv("~/Desktop/Visium/2026-Mar-Visium/6P17_D1/outs/nppa_RGC_markers.csv")

sig_markers = marks_0.2$gene[marks_0.2$p_val_adj<0.05]

marks_sub = subset(markers, FeatureID %in% sig_markers)

nPPa_markers = subset(marks_sub,
                      nppa_RGCs.Log2.Fold.Change > 1 & 
                        nppa_RGCs.P.Value  <0.05
                        )
sub_1 = AddModuleScore_UCell(sub_1, 
                       features = list(nppa= c(nPPa_markers$FeatureID)), 
                       name = 'nPPA_RGC')
DotPlot(sub_1, 'nppanPPA_RGC', scale = F)
# so all types are represented except 1_0, 1_1 is most common
library(CytoTRACE)
sub_1$CytoTRACE = CytoTRACE(sub_1@assays$RNA$data%>%as.matrix())$CytoTRACE

DotPlot(sub_1, 'CytoTRACE', scale = F)
# 1_1 has the highest cytoTRACE score

clown_go(nPPa_markers$FeatureID)%>%dotplot()
clown_go(marks_0.2$gene[marks_0.2$cluster=='1_1' & marks_0.2$p_val_adj<0.05])%>%dotplot()

clown_go(marks_0.2$gene[marks_0.2$cluster=='1_0' & marks_0.2$p_val_adj<0.05])%>%dotplot()
# 1_0 still showing these markers 
clown_go(marks_0.2$gene[marks_0.2$cluster=='1_3' & marks_0.2$p_val_adj<0.05])%>%dotplot()
# gfap positive cluster

clown_go(marks_0.2$gene[marks_0.2$cluster=='1_2' & marks_0.2$p_val_adj<0.05])%>%dotplot()
# arob + 


DimPlot(sub_1)
# this seems more right to me?


## props ###
cells_subcluster = sub_1@meta.data%>%
  group_by(individual, Status, sub_res0.2)%>%
  summarize(n_cells = n())

# summarize total cells
cells_total = sub_1@meta.data%>%
  group_by(individual)%>%
  summarize(total_cells = n())

# join and get proportions
joined = cells_subcluster%>%
  right_join(cells_total, by = 'individual')%>%
  mutate(proportion = n_cells/total_cells)%>%
  subset(individual != 'GH')

joined$Status = factor(joined$Status, levels =c('M',
                                 'D',
                                 "E",
                                 'NF',
                                 'F'))
ggplot(joined, aes(x = Status, y = proportion))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~sub_res0.2, scales = 'free')
# 1_0 spikes and 1_1 dips

# this might be a better way to think about it, 1_1 is the stem and 
# 1_0 is the progenitor?, they seem to change complementarily


clown_go(marks_0.2$gene[marks_0.2$cluster=='1_1' & marks_0.2$p_val_adj<0.05])%>%dotplot()

clown_go(marks_0.2$gene[marks_0.2$cluster=='1_0' & marks_0.2$p_val_adj<0.05])%>%dotplot()
# 1_0 still showing these markers 

# glm ####
out_dat = data.frame()
for(i in 0:3){
  dat_for_mod = subset(joined,sub_res0.2 == paste0('1_',i) )
  model_mat = cbind(dat_for_mod$n_cells, dat_for_mod$total_cells- dat_for_mod$n_cells)
  
  mod = glm(model_mat~Status,
            data = dat_for_mod, 
            family = binomial('logit'))  
  
  av_ = anova(mod, test ='Chisq')
  av_p = av_$`Pr(>Chi)`[2]
  
  pair_ = pairs(emmeans(mod, 'Status'), adjust = 'none')%>%as.data.frame()
  m_d_p = pair_$p.value[pair_$contrast=='M - D']
  m_f_p = pair_$p.value[pair_$contrast=='M - F']
  d_f_p = pair_$p.value[pair_$contrast=='D - F']
  
  newd = data.frame(sub.cluster =i, 
                    av_p = round(av_p, 5),
                    m_d_p = m_d_p, 
                    m_f_p = m_f_p, 
                    d_f_p = d_f_p)

  out_dat = rbind(out_dat, newd)
  
}

### have claude plot it for me ####
colors = c('#1965B0', '#4EB265', '#F7F056', '#7BAFDE', '#DC050C')

status_to_phase = c(
  "D"='I',
  'M' = 'M',
  'F' ='F',
  'NF' ='NF',
  'NRM' = 'NRM',
  'E' = 'LI'
)

p_value_annotation <- function(p) {
  if (p < 0.001) {
    return(list(label = "***", size = 7, adjust = 0.6))
  } else if (p < 0.01) {
    return(list(label = "**", size = 7, adjust = 0.6))
  } else if (p < 0.05) {
    return(list(label = "*", size = 7, adjust = 0.6))
  } else {
    return(list(label = paste0("p = ", round(p, 3)), size = 2, adjust = 0))
  }
}


# --- Prep proportions (mirrors your existing pattern) ---
total_cells_sub1 <- sub_1@meta.data %>%
  group_by(individual, Status) %>%
  summarize(total_cells = n(), .groups = 'drop')

n_cells_sub1 <- sub_1@meta.data %>%
  group_by(individual, sub_res0.2) %>%
  summarize(ncells = n(), .groups = 'drop')

joint_sub1 <- total_cells_sub1 %>%
  left_join(n_cells_sub1, by = 'individual')

joint_sub1$Status <- as.character(joint_sub1$Status)

# --- GLM loop across subclusters ---
dat_sub1 <- data.frame()

for(i in unique(joint_sub1$sub_res0.2)){
  sub <- subset(joint_sub1, sub_res0.2 == i & Status %in% c('M', 'D', 'F'))
  mat <- cbind(sub$ncells, sub$total_cells - sub$ncells)
  
  mod  <- glm(mat ~ Status, data = sub, family = 'binomial')
  p    <- anova(mod, test = 'Chisq')
  p.value <- p$`Pr(>Chi)`[2]
  
  pair <- pairs(emmeans(mod, 'Status'), adjust = 'none') %>% as.data.frame()
  
  newd <- data.frame(
    cluster       = i,
    p             = p.value,
    d_m_p.value   = pair$p.value[pair$contrast == 'D - M'],
    f_m_p.value   = pair$p.value[pair$contrast == 'F - M'],
    d_f_p.value   = pair$p.value[pair$contrast == 'D - F'],
    d_m_estimate  = pair$estimate[pair$contrast == 'D - M'],
    f_m_estimate  = pair$estimate[pair$contrast == 'F - M'],
    d_f_estimate  = pair$estimate[pair$contrast == 'D - F']
  )
  dat_sub1 <- rbind(dat_sub1, newd)
}

dat_sub1$av_q.value <- p.adjust(dat_sub1$p, 'fdr', nrow(dat_sub1))

# --- Adapted proportion_plotter for subclusters ---
# Only change: cluster_col points to sub_res0.2 and NF is kept
proportion_plotter_sub <- function(data_frame, cluster,
                                   cluster_col = 'sub_res0.2',
                                   d_m_p = NULL, f_m_p = NULL, d_f_p = NULL) {
  
  colors = c('#1965B0', '#4EB265', '#F7F056', '#7BAFDE', '#DC050C')
  
  dag <- subset(data_frame, data_frame[[cluster_col]] == cluster & 
                  Status != 'NRM')  # keep NF, only drop NRM
  dag$Phase <- status_to_phase[dag$Status]
  dag$Phase <- factor(dag$Phase, levels = c('M', 'I', 'LI', 'NF', 'F'))
  dag$prop  <- dag$ncells / dag$total_cells
  
  dag_summary <- dag %>%
    group_by(Phase) %>%
    summarise(
      mean_prop = mean(prop, na.rm = TRUE),
      se        = sd(prop, na.rm = TRUE) / sqrt(n()),
      .groups   = 'drop'
    )
  
  y_max      <- max(dag$prop, na.rm = TRUE)
  y_lower    <- y_max * 1.10
  y_upper    <- y_max * 1.35
  y_axis_max <- y_max * 1.50
  
  d_m_fmt <- if (!is.null(d_m_p)) p_value_annotation(d_m_p) else list(label = '', size = 3, adjust = 0)
  f_m_fmt <- if (!is.null(f_m_p)) p_value_annotation(f_m_p) else list(label = '', size = 3, adjust = 0)
  d_f_fmt <- if (!is.null(d_f_p)) p_value_annotation(d_f_p) else list(label = '', size = 3, adjust = 0)
  
  p <- ggplot(dag, aes(x = Phase, y = prop)) +
    geom_bar(data = dag_summary,
             aes(x = Phase, y = mean_prop, fill = Phase),
             stat = 'identity', width = 0.6, inherit.aes = FALSE) +
    geom_errorbar(data = dag_summary,
                  aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                  width = 0.4, linewidth = 0.5, inherit.aes = FALSE) +
    geom_point(size = 1, position = position_jitter(width = 0.1, seed = 42)) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, y_axis_max),
                       expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    theme_classic() +
    labs(x = 'Phase', y = '% of 1_RGC',
         title = paste0('Subcluster ', cluster)) +
    theme(legend.position = 'none')
  
  if (!is.null(d_m_p)) {
    p <- p + geom_signif(xmin = 1, xmax = 1.9, y_position = y_lower,
                         annotation = d_m_fmt$label, color = 'black',
                         tip_length = c(0, 0), textsize = d_m_fmt$size,
                         vjust = d_m_fmt$adjust)
  }
  if (!is.null(d_f_p)) {
    p <- p + geom_signif(xmin = 2.1, xmax = 5, y_position = y_lower,
                         annotation = d_f_fmt$label, color = 'black',
                         tip_length = c(0, 0), textsize = d_f_fmt$size,
                         vjust = d_f_fmt$adjust)
  }
  if (!is.null(f_m_p)) {
    p <- p + geom_signif(xmin = 1, xmax = 5, y_position = y_upper,
                         annotation = f_m_fmt$label, color = 'black',
                         tip_length = c(0, 0), textsize = f_m_fmt$size,
                         vjust = f_m_fmt$adjust)
  }
  
  return(p)
}

# --- Generate and save plots ---
plots_sub1 <- list()

for(i in unique(dat_sub1$cluster)){
  plots_sub1[[i]] <- proportion_plotter_sub(
    joint_sub1, i,
    d_m_p = dat_sub1$d_m_p.value[dat_sub1$cluster == i],
    f_m_p = dat_sub1$f_m_p.value[dat_sub1$cluster == i],
    d_f_p = dat_sub1$d_f_p.value[dat_sub1$cluster == i]
  )
}

# View all together
wrap_plots(plots_sub1)

# Save individually
#for(i in names(plots_sub1)){
 # ggsave(
#    plot     = plots_sub1[[i]],
#    file     = paste0("subclust_", i, "_prop.svg"),
#    device   = "svg",
#    units    = "in",
#    width    = 2, height = 2,
#    path     = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2.1/new_rgsc_subs_prop_cells'
#  )
#}

#### plot UMAP ####
dimplot = DimPlot(sub_1, label = T, raster =T, 
                  pt.size = 4,
                  raster.dpi= c(1024, 1024))+
  theme_void()+
  theme(legend.position = "none")
dimplot

#ggsave(plot = dimplot,
#       file = "UMAP_sub1.svg",
#       device = "svg",
#       units = "in",
#       width = 4,
#       height = 4,
#       path = "Manuscript/Plots/Manuscript v1.2.1/")

#### plot GO ####

go_1_1= clown_go(marks_0.2$gene[marks_0.2$cluster=='1_1' & marks_0.2$p_val_adj<0.05])%>%dotplot()

go_1_0=clown_go(marks_0.2$gene[marks_0.2$cluster=='1_0' & marks_0.2$p_val_adj<0.05])%>%dotplot()
# 1_0 still showing these markers 
go_1_3 =clown_go(marks_0.2$gene[marks_0.2$cluster=='1_3' & marks_0.2$p_val_adj<0.05])%>%dotplot()
# gfap positive cluster

go_1_2 = clown_go(marks_0.2$gene[marks_0.2$cluster=='1_2' & marks_0.2$p_val_adj<0.05])%>%
  dotplot()

for(i in 0:3){
  
  p = get(paste0('go_1_',i))+
    labs(title =paste0('1_',i))
  
    #ggsave(plot = p,
    #     file = paste0("go_1_", i, ".svg"),
    #     device = 'svg',
    #     units = "in",
    #     width = 5,
    #     height = 3,
    #     path = "Manuscript/Plots/Manuscript v1.2.1/RGC_sub_marker_go")

}


#### markers #####
mark_plot = DotPlot(sub_1, 
        c('gfap',
          'gja1b',
          'LOC111577263',
          'angpt1',
          'slc6a11b',
          'slc6a1b',
          'slc4a4a',
          'sema4gb',
           'pax6b',
          'sox2',
          'foxg1a',
          'LOC111582483',
          'shha',
          'fgfr2',
           'fgfr3',
          'vegfab'
          ))+
  coord_flip()

#    ggsave(plot = mark_plot,
 #        file = "RGC_sub_markers_dot.svg",
  #       device = 'svg',
   #      units = "in",
    #     width = 4.75,
     #    height = 5,
      #   path = "Manuscript/Plots/Manuscript v1.2.1/")

#### plot cyto ####
library(multcomp)

cyt_dat = sub_1@meta.data%>%
  group_by(sub_res0.2, individual, Status)%>%
  summarize(mean_cyto = mean(CytoTRACE))%>%
  subset(Status!= 'NRM')

#run stats #
cyt_model = lm(mean_cyto~sub_res0.2, data = cyt_dat)
anova(cyt_model, test = 'Chisq')

p1 = pairs(emmeans(cyt_model, 'sub_res0.2'), adjust = 'none')

cld_cyt = cld(emmeans(cyt_model, 'sub_res0.2'), Letters = letters, adjust = "none", alpha = 0.05)
cld_cyt_df = as.data.frame(cld_cyt)

cyt_dat$sub_res0.2 = factor(cyt_dat$sub_res0.2, levels = c('1_1',
                                                           '1_0',
                                                           '1_2',
                                                           '1_3'))
cyt_plot =ggplot(cyt_dat, aes(x = (sub_res0.2), y = mean_cyto))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(0.75), shape = 1)+
  theme_minimal()+
  labs(y = 'CytoTRACE', x = 'Subcluster')+
    geom_text(data = cld_cyt_df, aes(x = sub_res0.2, y = 1, label = .group), 
          size = 3, inherit.aes = FALSE)

    ggsave(plot = cyt_plot,
         file = "RGC_sub_cyto.svg",
         device = 'svg',
         units = "in",
         width = 2,
         height = 2,
         path = "Manuscript/Plots/Manuscript v1.2.1/")
    

    

