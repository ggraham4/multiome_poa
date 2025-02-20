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
clown_go<- readRDS('Functions/clown_go')
`%notin%`<- Negate(`%in%`) 


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


together_data_defined_summed_plot$colors <- NA
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Early Upregulated', 'red', NA)
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Early Downregulated', 'green', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Late Downregulated', 'blue', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Late Upregulated', 'pink', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Progressively Upregulated', 'cyan', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Progressively Downregulated', 'hotpink', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Terminally Downregulated', 'orange', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Terminally Upregulated', 'maroon', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Transiently Downregulated', 'yellow', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Transiently Upregulated', 'darkgreen', together_data_defined_summed_plot$colors )

together_data_defined_summed_plot$class <- factor(together_data_defined_summed_plot$class, 
                                                  levels = rev(c('Transiently Upregulated',
                                                             'Transiently Downregulated',
                                                             'Terminally Upregulated',
                                                             "Terminally Downregulated",
                                                             'Progressively Upregulated',
                                                             'Progressively Downregulated',
                                                             'Late Upregulated',
                                                             'Late Downregulated',
                                                             'Early Upregulated',
                                                             'Early Downregulated')))


dpg_bar_plot <- ggplot(subset(together_data_defined_summed_plot, cluster %notin% c(15,30)), aes(x = as.factor(cluster), y = class_count.x, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DPGs", fill = "Expression Pattern") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
    scale_fill_manual(values = unique(together_data_defined_summed_plot$colors))+
geom_text(aes(label = issignif, y = class_count.y), size = 8)+
    theme_minimal()+
    theme(axis.text.x = element_text(size = 10, vjust =0.4, angle = -90),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size=8),
        legend.title  = element_text(size=12),
        legend.background  = element_rect(color = 'white'),
        #legend.position = c(.7,.6),
        legend.direction = 'vertical',
        legend.title.position = 'top',
        legend.position = 'none'
        )+
  ylim(0,550)
dpg_bar_plot

ggsave(plot = dpg_bar_plot,
       file = "dpg_bar_plot.svg",
       device = "svg",
       units = "in",
       width = 4,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 2")




dpg_bar_plot2 <- ggplot(subset(together_data_defined_summed_plot, cluster %notin% c(15,30)), aes(x = as.factor(cluster), y = class_count.x, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DPGs", fill = "Expression Pattern") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
    scale_fill_manual(values = unique(together_data_defined_summed_plot$colors))+
  theme(legend.direction = 'horizontal',
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24))
dpg_bar_plot2

#ggsave(plot = dpg_bar_plot2,
#       file = "legend.svg",
#       device = "svg",
#       units = "in",
#       width = 20,
#       height = 10,
#       path = "Bachelors Thesis/Plots/Figure 2")



#### i have no idea why I cant replicate the go terms the data are identical (I checked) and the definition functions are identical (I checked),
### I hope its not something to do with how I subset
together_data_test <- together_data[together_data$issignif=='*' & is.na(together_data$warning) & !is.na(together_data$issignif),]

#### something mustc change here
prop_data_genes <- together_data$gene[together_data$issignif=='*' & is.na(together_data$warning) & !is.na(together_data$issignif)]
unique_prop_data_genes <-unique(prop_data_genes)

exp_data_genes <- together_data_exp$gene[together_data_exp$issignif=='*' & is.na(together_data_exp$warning) & !is.na(together_data_exp$issignif)]
unique_exp_data_genes <-unique(exp_data_genes)

all_genes <- unique(append(unique_exp_data_genes,unique_prop_data_genes))

together_genes <- data.frame(genes = all_genes)

together_genes$type <- ifelse(together_genes$genes%in%prop_data_genes & !(together_genes$genes%in%exp_data_genes),
                                           'prop',
                                           NA)
######## ^^^ referencing investigating prop and expression DEGs



go_all <- clown_go(together_data_test$gene[!is.na(together_data_test$issignif) & together_data_test$cluster %notin% c(15,30) & is.na(together_data_test$warning)])

go_all_plot <- dotplot(go_all)+ 
  #coord_flip()+
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = -45),
        legend.position = 'none')+
  scale_size_continuous(range = c(1,3))
go_all_plot

  

