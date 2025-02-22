#Negative binomial lower stringency 
{
  library(parallel)

  library(Seurat)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(emmeans)
  library(openxlsx)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
library(biomaRt)
  library(Polychrome)
  P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
names(P40) <- NULL

`%notin%` <- Negate(`%in%`)

  mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
clown_go<- readRDS('Functions/clown_go')
define_degs<- readRDS('Functions/define_degs')

}

### ANALYSIS ####
together_data <- data.frame()
for(i in 0:31){
  if(i == 30 | i == 15){next}
  data <- read.csv(paste0('DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_', i, '.csv'))
  #data <- get(paste0('results_cluster',i))
  data$cluster <- i
  together_data <- rbind(together_data, data)
}

together_data_defined <- define_degs(together_data)
together_data_defined_only_signif <- together_data_defined[!is.na(together_data_defined$class),]

together_data_defined_summed <- together_data_defined%>%
  subset(!is.na(class))%>%
  group_by(class, cluster)%>%
  summarize(class_count = n())

ggplot(together_data_defined_summed, aes(x = cluster, y = class_count, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DEGs") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
    scale_x_continuous(breaks = c(0:31))+
  scale_y_continuous()+
  scale_fill_manual(values = P40)

#### Chisq test for enrichment ####
together_data_chisq <-  together_data_defined%>%
  subset(!is.na(class))%>%
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

together_data_defined_summed_plot <- together_data_defined_summed%>%
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

deg_bar_plot <- ggplot(together_data_defined_summed_plot, aes(x = as.factor(cluster), y = class_count.x, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DEGs", fill = 'Expression Pattern') +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
geom_text(aes(label = issignif, y = class_count.y), size = 8)+
    scale_fill_manual(values = unique(together_data_defined_summed_plot$colors))+
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
  ylim(0,60)
deg_bar_plot

#ggsave(plot = deg_bar_plot,
 #      file = "deg_bar_plot.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 4,
     #  height = 2,
      # path = "Bachelors Thesis/Plots/Figure 2")


go_all <- clown_go(together_data_defined$gene[!is.na(together_data_defined$class) & together_data_defined$cluster %notin%c(15,30)])

go_all_plot <- dotplot(go_all)+ 
  #coord_flip()+
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = -45),
        legend.position = 'none')+
  scale_size_continuous(range = c(1,3))
go_all_plot

ggsave(plot = go_all_plot,
       file = "go_all_plot.svg",
       device = "svg",
       units = "in",
       width = 5.5,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 2")


                        