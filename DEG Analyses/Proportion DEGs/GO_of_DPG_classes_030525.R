#go of all prop deg classes
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


options(future.globals.maxSize = 23 * 1024^3) 

#obj <- readRDS( '/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object with SCT.rds')

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

together_data_defined_only_signif <- together_data[!is.na(together_data$class),]

table(together_data_defined_only_signif$class)

##### GENES THAT CHANGE BEFORE GONADS CHANGE -- early DEGs
early_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Early Upregulated', "Early Downregulated")]

clown_go(early_degs)%>%dotplot()

### genes that change after gonads change -- late DEGs
late_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Late Upregulated', "Late Downregulated")]

clown_go(late_degs)%>%dotplot()

### genes that change as gonads change -- transient DEGs
transient_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Transiently Upregulated', "Transiently Downregulated")]

clown_go(transient_degs)%>%dotplot()

## genes that change at the ends of sex change -- terminal DEGs
terminal_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Terminally Upregulated', "Terminally Downregulated")]
clown_go(terminal_degs)%>%dotplot()

# progressive degs -- are there even any 
progressive_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Progressively Upregulated', "Progressively Downregulated")]

clown_go(progressive_degs)%>%dotplot()


### early up 
early_up_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Early Upregulated')]
clown_go(early_up_degs)%>%dotplot()

### early down
early_down_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Early Downregulated')]
clown_go(early_down_degs)%>%dotplot()

###late up 
late_up_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Late Upregulated')]
clown_go(late_up_degs)%>%dotplot()

### late down
late_down_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Late Downregulated')]
clown_go(late_down_degs)%>%dotplot()

### term up 
term_up_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Terminally Upregulated')]
clown_go(term_up_degs)%>%dotplot()

## term down
term_down_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Terminally Downregulated')]
clown_go(term_down_degs)%>%dotplot()

all_dpgs <- together_data_defined_only_signif$gene
clown_go(all_dpgs)%>%dotplot()









