### New prop DEGs analysis

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
  
  library(Polychrome)
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
names(P40) <- NULL

mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
define_degs_prop<- readRDS('Functions/define_degs_prop.rds')
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
prop_deg_function<- readRDS('Functions/DEG_functions/prop_deg_function.rds')

}

###Analysis ####
for(i in c(0:31)){
  if(i == 15 | i == 30){next}
  print(i)
  data <- read.csv(
              paste0('~/Desktop/multiome_poa/DEG Outputs/Prop DEGs 011525/prop_degs_cluster_', i, '.csv'))
    
  data <-define_degs_prop(data) 
   assign(paste0('prop_degs_cluster_', i), data, envir = .GlobalEnv)
  }

together_data <- data.frame()
for(i in 0:31){
    if(i == 15 | i == 30){next}
  data <- get(paste0('prop_degs_cluster_',i))
  together_data <- rbind(together_data, data)
}

together_data_summed <- together_data%>%
    subset(!is.na(class)& is.na(warning))%>%
  group_by(cluster, class)%>%
  summarize(class_count = n())

together_data_chisq <-  together_data%>%
  subset(!is.na(class)& is.na(warning))%>%
  group_by(cluster)%>%
  summarize(class_count = n())

total_sum <- sum(together_data_chisq$class_count)
mean_sum <- total_sum / nrow(together_data_chisq)
chisq_test <- chisq.test(together_data_chisq$class_count, p = rep(1/nrow(together_data_chisq), nrow(together_data_chisq)))
expected <- chisq_test$expected
residuals <- (together_data_chisq$class_count - expected) / sqrt(expected)
together_data_chisq$residuals <- residuals
together_data_chisq$significant <- together_data_chisq$residuals > 2
together_data_chisq$issignif <- ifelse(together_data_chisq$significant==T, '*',NA)

together_data_defined_summed_plot <- together_data_summed%>%
  right_join(together_data_chisq, by = 'cluster')

`%notin%`<- Negate(`%in%`) 

dpg_bar_plot <- ggplot(subset(together_data_defined_summed_plot, cluster %notin% c(15,30)), aes(x = as.factor(cluster), y = class_count.x, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DPGs", fill = "Expression Pattern") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
  scale_fill_manual(values = P40)+
geom_text(aes(label = issignif, y = class_count.y), size = 8)+
    theme_minimal()+
    theme(axis.text.x = element_text(size = 8, vjust =0.4, angle = -90),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=10),
        legend.background  = element_rect(color = 'white'),
        legend.position = c(.7,.6),
        legend.direction = 'vertical',
        legend.title.position = 'top')+
  ylim(0,550)
dpg_bar_plot

ggsave(plot = dpg_bar_plot,
       file = "dpg_bar_plot.svg",
       device = "svg",
       units = "in",
       width = 3.5,
       height = 3,
       path = "Bachelors Thesis/Plots/Figure 2")




  

