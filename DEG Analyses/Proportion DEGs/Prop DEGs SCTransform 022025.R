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
prop_deg_function_sctransform<- readRDS('Functions/DEG_functions/prop_deg_function_sctransform.rds')
clown_go<- readRDS('Functions/clown_go')

}

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

options(future.globals.maxSize = 23 * 1024^3) 

#obj <- SCTransform(obj)
#saveRDS(obj, '/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object with SCT.rds')

obj <- readRDS( '/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object with SCT.rds')

for (i in c(0:31)) {
  if(i == 15 | i == 30){next}
  print(i)
  start <-Sys.time()
   prop_deg_data<-  prop_deg_function_sctransform(object = obj, cluster = i)

  assign(paste0('prop_degs_cluster_', i), prop_deg_data, envir = .GlobalEnv)
  
    write.csv(prop_deg_data, 
              paste0('DEG Outputs/022025 Prop DEGs SCTransform Matrix/prop_degs_cluster_', i, '.csv'))
end <- Sys.time()
print(end-start)
    }

###Analysis ####
together_data <- data.frame()
for(i in c(0:31)){
    if(i == 15 | i == 30){next}
  print(i)
  data <- read.csv(
              paste0('DEG Outputs/022025 Prop DEGs SCTransform Matrix/prop_degs_cluster_', i, '.csv'))
  data <-define_degs_prop(data) 
     together_data <- rbind(together_data, data)
  }


together_data_summed <- together_data%>%
    subset(!is.na(class)& is.na(warning) & singular==F & prop_expressing_D >0.05 & prop_expressing_F>0.05 & prop_expressing_M>0.05)%>%
  group_by(cluster, class)%>%
  summarize(class_count = n())

together_data_chisq <-  together_data%>%
    subset(!is.na(class)& is.na(warning) & singular==F & prop_expressing_D >0.05 & prop_expressing_F>0.05 & prop_expressing_M>0.05)%>%
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
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Early Upregulated', 'blue', NA)
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Early Downregulated', 'red', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Late Downregulated', 'pink', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Late Upregulated', 'orange', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Progressively Upregulated', 'cyan', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Progressively Downregulated', 'hotpink', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Terminally Downregulated', 'maroon', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Terminally Upregulated', 'darkgreen', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Transiently Downregulated', 'yellow', together_data_defined_summed_plot$colors )
together_data_defined_summed_plot$colors <- ifelse(together_data_defined_summed_plot$class == 'Transiently Upregulated', 'green', together_data_defined_summed_plot$colors )

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


dpg_bar_plot <- ggplot(together_data_defined_summed_plot, aes(x = as.factor(cluster), y = class_count.x, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DPGs", fill = 'Expression Pattern') +
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
  ylim(0,130)
dpg_bar_plot

tab <-  as.data.frame(table(together_data_subset$class))
sum(tab$Freq)

length(together_data$gene[!is.na(together_data$class)])

ggsave(plot = dpg_bar_plot,
       file = "dpg_bar_plot.svg",
       device = "svg",
       units = "in",
       width = 4,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 2")


together_data_subset <- together_data%>%
    subset(!is.na(class)& is.na(warning) & singular==F & prop_expressing_D >0.05 & prop_expressing_F>0.05 & prop_expressing_M>0.05)

clown_go(together_data_subset$gene[together_data_subset$class=='Late Downregulated'])%>%dotplot()
clown_go(together_data_subset$gene[together_data_subset$class=='Early Upregulated'])%>%dotplot()
clown_go(together_data_subset$gene[together_data_subset$class=='Transiently Upregulated'])%>%dotplot()

go_early_up <- clown_go(together_data_subset$gene[together_data_subset$class=='Early Upregulated'])%>%clusterProfiler::dotplot()


go_early_up_plot <-go_early_up+ 
  #coord_flip()+
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = -45),
        legend.position = 'none')+
  geom_point(fill = '#7c5c6c', color = 'black')+
  scale_size_continuous(range = c(1,3))
go_early_up_plot

ggsave(plot = go_early_up_plot,
       file = "go_early_up_plot.svg",
       device = "svg",
       units = "in",
       width = 5.5,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 2")


  

