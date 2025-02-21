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

}

#obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

#options(future.globals.maxSize = 23 * 1024^3) 

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
              paste0('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/Prop DEGs 011525/prop_degs_cluster_', i, '.csv'))
  data <-define_degs_prop(data) 
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


ggplot(together_data_defined_summed_plot, aes(x = cluster, y = class_count.x, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DEGs") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
    scale_x_continuous(breaks = c(0:31))+
  scale_y_continuous()+
  scale_fill_manual(values = P40)+
geom_text(aes(label = issignif, y = class_count.y), size = 10)





  

