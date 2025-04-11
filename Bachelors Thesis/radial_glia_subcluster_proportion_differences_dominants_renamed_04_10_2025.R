#radial glial differences in subcluster proportion
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
  library(glmGamPoi)
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)
  library(emmeans)
  library(CytoTRACE)
  library(ggrepel)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(Polychrome)
  library(scCustomize)
  library(hdWGCNA)
  
  P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
  swatch(P40)
  names(P40) <- NULL
  
  mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
  prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
  define_degs_prop<- readRDS('Functions/define_degs_prop.rds')
  mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
  prop_deg_function.rds<- readRDS('Functions/DEG_functions/prop_deg_function.rds')
  define_behavior_degs<- readRDS('Functions/define_behavior_degs')
  clown_go<- readRDS('Functions/clown_go')
  define_degs<- readRDS('Functions/define_degs')
  
}

radial_glia <- readRDS(file='A:/anemonefish_multiome_radial_glia_03_20_2025/radial_glia_object.rds')

radial_glia_count_data_all <- radial_glia@meta.data%>%
  dplyr::select(-c('brown', 'yellow','blue','turquoise','grey'))%>%
  group_by(individual)%>%
  summarize(n_cells_total = n())
radial_glia_count_data_all

radial_glia_count_data_sub <- radial_glia@meta.data%>%
  dplyr::select(-c('brown', 'yellow','blue','turquoise','grey'))%>%
  group_by(individual, Status, sub)%>%
  summarize(n_cells_sub = n())
radial_glia_count_data_sub

radial_glia_counts = radial_glia_count_data_sub%>%
  right_join(radial_glia_count_data_all, by = 'individual')%>%
  filter(Status %in% c('M','D','F'))
radial_glia_counts

radial_glia_counts$n_cells_not_in_cluster = radial_glia_counts$n_cells_total-radial_glia_counts$n_cells_sub

out <- data.frame()
for(i in unique(radial_glia$sub)){
  data_subset = subset(radial_glia_counts, sub == i)
  
  model_matrix = cbind(data_subset$n_cells_sub, data_subset$n_cells_not_in_cluster)
  
  model = glmer(model_matrix ~ Status +(1|individual), family = binomial('logit'), data =data_subset)
  
  av = car::Anova(model, type = 'III')%>%as.data.frame()
  
  newd = data.frame(sub =i,
                    av = av$`Pr(>Chisq)`[2],
                    singular = isSingular(model))
  
  out <- rbind(out, newd)
}

#sex differences in 2_3 and 2_4
data_subset = subset(radial_glia_counts, sub == '2_4')

model_matrix = cbind(data_subset$n_cells_sub, data_subset$n_cells_not_in_cluster)

model = glmer(model_matrix ~ Status +(1|individual), family = binomial('logit'), data =data_subset)
car::Anova(model, type = 'III')%>%as.data.frame()

pairs(emmeans(model, 'Status'), adjust ='none')
#df different but all close tbh

plot_2_4_data = radial_glia_count_data_sub%>%
  right_join(radial_glia_count_data_all, by = 'individual')%>%
  filter(Status !='NRM' & sub =='2_4')%>%
  mutate(prop = n_cells_sub/n_cells_total)
plot_2_4_data

plot_2_4_data$Status <- ifelse(plot_2_4_data$Status=='D','I',plot_2_4_data$Status)
plot_2_4_data$Status <- ifelse(plot_2_4_data$Status=='E','LI',plot_2_4_data$Status)

plot_2_4_data$Status <- factor(plot_2_4_data$Status, levels = c('M','I','LI','NF','F'))

prop_plot_2_4 <- ggplot(plot_2_4_data, aes(x = Status, y = prop))+
  geom_boxplot(aes(color = Status, fill = Status), alpha = 0.25, outlier.shape=NA)+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'Proportion of Cells', title = '2_4')+
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(.2), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)+
  ylim(min(plot_2_4_data$prop),.23)
prop_plot_2_4

ggsave(plot = prop_plot_2_4,
       file = "prop_plot_2_4.svg",
       device = "svg",
       units = "in",
       width = 1.8,
       height = 2,
       path = "Bachelors Thesis/Plots_dominants_renamed/Figure 5")

### ok fine maybe I do care about 2_3

data_subset_23 = subset(radial_glia_counts, sub == '2_3')

model_matrix_23 = cbind(data_subset_23$n_cells_sub, data_subset_23$n_cells_not_in_cluster)

model_23 = glmer(model_matrix_23 ~ Status +(1|individual), family = binomial('logit'), data =data_subset_23)
car::Anova(model_23, type = 'III')%>%as.data.frame()

pairs(emmeans(model, 'Status'), adjust ='none')
#df is the major difference

plot_2_3_data = radial_glia_count_data_sub%>%
  right_join(radial_glia_count_data_all, by = 'individual')%>%
  filter(Status !='NRM' & sub =='2_3')%>%
  mutate(prop = n_cells_sub/n_cells_total)
plot_2_3_data

plot_2_3_data$Status <- ifelse(plot_2_3_data$Status=='D','I',plot_2_3_data$Status)
plot_2_3_data$Status <- ifelse(plot_2_3_data$Status=='E','LI',plot_2_3_data$Status)

plot_2_3_data$Status <- factor(plot_2_3_data$Status, levels = c('M','I','LI','NF','F'))

prop_plot_2_3 <- ggplot(plot_2_3_data, aes(x = Status, y = prop))+
  geom_boxplot(aes(color = Status, fill = Status), alpha = 0.25, outlier.shape=NA)+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'Proportion of Cells', title = '2_3')+
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(.2), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)+
  ylim(min(plot_2_3_data$prop),.23)
prop_plot_2_3

ggsave(plot = prop_plot_2_3,
       file = "prop_plot_2_3.svg",
       device = "svg",
       units = "in",
       width = 1.25,
       height = 2,
       path = "Bachelors Thesis/Plots_dominants_renamed/Figure 5")



