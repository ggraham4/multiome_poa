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
define_degs<- readRDS('Functions/define_degs')

mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
prop_deg_function<- readRDS('Functions/DEG_functions/prop_deg_function.rds')
clown_go<- readRDS('Functions/clown_go')
`%notin%`<- Negate(`%in%`) 


}

###Analysis ####
together_data <- data.frame()

for(i in c(0:31)){
  if(i == 15 | i == 30){next}
  print(i)
  data <- read.csv(
              paste0('~/Desktop/multiome_poa/DEG Outputs/Prop DEGs 011525/prop_degs_cluster_', i, '.csv'))
    
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

#ggsave(plot = dpg_bar_plot,
#       file = "dpg_bar_plot.svg",
#       device = "svg",
#       units = "in",
#       width = 4,
#       height = 2,
#       path = "Bachelors Thesis/Plots/Figure 2")




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
#### Trying to replicate #####

### Read in exp degs ####
together_data_exp <- data.frame()
for(i in 0:31){
  if(i == 30 | i == 15){next}
  data <- read.csv(paste0('DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_', i, '.csv'))
  #data <- get(paste0('results_cluster',i))
  data <- define_degs(data)
  data$cluster <- i
  data$type = 'exp'
  together_data_exp <- rbind(together_data_exp, data)
}

together_data$type = 'prop'

colnames(together_data)
colnames(together_data_exp)
together_data_exp$issignif <- ifelse(!is.na(together_data_exp$class), '*', NA)

together_data_prop = dplyr::select(together_data, c('gene','type','issignif', 'warning', 'singular'))
together_data_exp = dplyr::select(together_data_exp, c('gene','type','issignif', 'warning', 'singular'))

together_data_full = rbind(together_data_prop,together_data_exp)
together_data_full_signif <- subset(together_data_full, issignif =='*' & is.na(warning))


#### something mustc change here
prop_data_genes <- together_data_full_signif$gene[together_data$type == 'prop']
unique_prop_data_genes <-unique(prop_data_genes)

exp_data_genes <- together_data_full_signif$gene[together_data_full_signif$type =='exp']
unique_exp_data_genes <-unique(exp_data_genes)

all_genes <- unique(append(unique_exp_data_genes,unique_prop_data_genes))

go_prop <- clown_go(unique_prop_data_genes)%>%
  dotplot()
go_prop
###THIS IS YET ANOTHER DIFFERENT RESULT WHY

go_exp <- clown_go(unique_exp_data_genes)%>%
  dotplot
go_exp ### ok this is correct

go_both <- clown_go(together_data_full_signif$gene)%>%
  dotplot()
go_both

###was it in non overlapping genes
prop_no_overlap <- unique_prop_data_genes[unique_prop_data_genes %notin% unique_exp_data_genes]
go_prop_no_overlap <- clown_go(prop_no_overlap)%>%
  dotplot()
go_prop_no_overlap
#ok thats where it came from but its still not right

######## ^^^ referencing investigating prop and expression DEGs #### i cant fully replicate it but clearly I did something wrong #####
###### moving on #####

#start over
together_data <- data.frame()

for(i in c(0:31)){
  if(i == 15 | i == 30){next}
  print(i)
  data <- read.csv(
              paste0('~/Desktop/multiome_poa/DEG Outputs/Prop DEGs 011525/prop_degs_cluster_', i, '.csv'))
    
  data <-define_degs_prop(data) 
  together_data <- rbind(together_data, data)
  }

together_data_subset <- subset(together_data, !is.na(class)&
                                 is.na(warning) & 
                                 cluster %notin% c(15,30))


go_all <- clown_go(together_data_subset$gene)

go_dpg_plot <- dotplot(go_all)+ 
  #coord_flip()+
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = -45),
        legend.position = 'none')+
  scale_size_continuous(range = c(1,3))
go_dpg_plot

#ggsave(plot = go_dpg_plot,
#       file = "go_dpg_plot.svg",
#       device = "svg",
#       units = "in",
#       width = 4,
#       height = 2,
#       path = "Bachelors Thesis/Plots/Figure 2")




enriched_degs <-   unlist(strsplit(go_all$geneID[go_all$Description =='regulation of small GTPase mediated signal transduction'], "/"))

together_data_subset$cluster[together_data_subset$gene %in% enriched_degs
                             ]%>% table()%>% plot()

together_data_subset$cluster[together_data_subset$gene %in% 'sipa1l1'
                             ]%>% table()

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

prop_cluster_plot(obj,
                  'sipa1l1',
                  7)
prop_cluster_plot(obj,
                  'sipa1l1',
                  8) # late downreg

together_data_subset$class[together_data_subset$gene %in%enriched_degs
                             ]%>% table()
### all late downregulated


##### explore the together data a little ######
together_data_filtered <- subset(together_data,
                                 is.na(warning) & 
                                 cluster %notin% c(15,30))

##negative binomially distributed
hist(together_data_filtered$cells_expressing_D)
hist(together_data_filtered$cells_expressing_F)
hist(together_data_filtered$cells_expressing_M)

#Kind of like a normal dist idrk
hist(log(together_data_filtered$cells_expressing_M))
hist(log(together_data_filtered$cells_expressing_D))
hist(log(together_data_filtered$cells_expressing_F))

## majority of genes are singular
table(together_data_filtered$singular)

table(together_data_filtered$issignif)

#negative binomially distributed
hist(together_data_filtered$prop_expressing_D)
hist(together_data_filtered$prop_expressing_F)
hist(together_data_filtered$prop_expressing_M)

#basically normal
hist(log(together_data_filtered$prop_expressing_D))
hist(log(together_data_filtered$prop_expressing_F))
hist(log(together_data_filtered$prop_expressing_M))

hist(together_data_filtered$d_m_q.value)
hist(together_data_filtered$d_f_q.value)
hist(together_data_filtered$f_m_q.value)

#0.44
min(together_data_filtered$prop_expressing_D[!is.na(together_data_filtered$issignif)])

#0.099
min(together_data_filtered$prop_expressing_F[!is.na(together_data_filtered$issignif)])

#0.38
min(together_data_filtered$prop_expressing_M[!is.na(together_data_filtered$issignif)])

#not that many < 1%
length(together_data_filtered$gene[together_data_filtered$prop_expressing_M<0.01 & !is.na(together_data_filtered$issignif)])
length(together_data_filtered$gene[together_data_filtered$prop_expressing_F<0.01 & !is.na(together_data_filtered$issignif)])
length(together_data_filtered$gene[together_data_filtered$prop_expressing_D<0.01 & !is.na(together_data_filtered$issignif)])

#Fair amount less than 5%
length(together_data_filtered$gene[together_data_filtered$prop_expressing_M<0.05 & !is.na(together_data_filtered$issignif)])
length(together_data_filtered$gene[together_data_filtered$prop_expressing_F<0.05 & !is.na(together_data_filtered$issignif)])
length(together_data_filtered$gene[together_data_filtered$prop_expressing_D<0.05 & !is.na(together_data_filtered$issignif)])

together_data_filtered$prop_expressing_D[(together_data_filtered$cluster==19) & (together_data_filtered$gene == 'tacr3a')]
together_data_filtered$prop_expressing_F[(together_data_filtered$cluster==19) & (together_data_filtered$gene == 'tacr3a')]
together_data_filtered$prop_expressing_M[(together_data_filtered$cluster==19) & (together_data_filtered$gene == 'tacr3a')]
together_data_filtered$class[(together_data_filtered$cluster==19) & (together_data_filtered$gene == 'tacr3a')]

together_data_filtered$prop_expressing_D[(together_data_filtered$cluster==19) & (together_data_filtered$gene == 'pgr')]
together_data_filtered$prop_expressing_F[(together_data_filtered$cluster==19) & (together_data_filtered$gene == 'pgr')]
together_data_filtered$prop_expressing_M[(together_data_filtered$cluster==19) & (together_data_filtered$gene == 'pgr')]


#i think 5% is a good cutoff
together_data_filtered_cutoff_5 <- together_data_filtered[
                                         (!is.na(together_data_filtered$issignif))&
                                            (together_data_filtered$prop_expressing_M>0.05)&
                                            (together_data_filtered$prop_expressing_F>0.05) &
                                           (together_data_filtered$prop_expressing_D >0.05)
                                           ,]

together_data_filtered_exclude <- together_data_filtered[
                                         (is.na(together_data_filtered$issignif))|
                                            (together_data_filtered$prop_expressing_M<0.05)|
                                            (together_data_filtered$prop_expressing_F<0.05) |
                                           (together_data_filtered$prop_expressing_D <0.05)
                                           ,]



together_data_filtered_cutoff_5_summed<- as.data.frame(table(together_data_filtered_cutoff_5$class, together_data_filtered_cutoff_5$cluster))

ggplot(together_data_filtered_cutoff_5_summed, aes(x = as.factor(Var2), y = Freq, fill = Var1), color = 'black')+
  geom_bar(stat = 'identity', position = 'stack')
### still the same bias though, late up or early down
### rerun model with offset of size factors
