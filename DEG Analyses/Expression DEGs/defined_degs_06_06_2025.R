library(dplyr)
define_degs2<- readRDS('Functions/define_degs2')

degs = "DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/"
degs_dir = dir("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering")

deg_data = data.frame()
for(i in degs_dir){
  dat = read.csv(file = paste0(degs, i))
  dat_filtered = dat%>%
    dplyr::select(-c(d_m_q.value,
             f_m_q.value,
             d_f_q.value))
  data_defined = define_degs2(dat_filtered)
  deg_data = rbind(deg_data, data_defined)
}

### test to confirm no classifications failed
deg_data_test <- subset(deg_data, av_q.value<0.05)
deg_data_test$wrong <- ifelse((deg_data_test$d_m_p.value > 0.05 & 
                                 deg_data_test$f_m_p.value > 0.05 & 
                                 deg_data_test$d_f_p.value > 0.05), 'fail',NA )

## ok so just to be clear, the pairs(emmeans()) function I used in the DEG function
## calls tukey by default, I then subset to genes where a type III anova has
## a q value less than 0.05 and apply my pairwise classification 

## save DEGs only, 
deg_data_0.05 <- subset(deg_data, av_q.value<0.05)
#write.csv(deg_data_0.05, 'DEG Outputs/defined_degs_anova_0.05_06_06_2025.csv')
deg_data_0.05 = read.csv('DEG Outputs/defined_degs_anova_0.05_06_06_2025.csv')
deg_data_0.05 = subset(deg_data_0.05, singular == F)

#### DEG classes ###
deg_classes <- table(deg_data_0.05$cluster, deg_data_0.05$class)%>%as.data.frame.matrix()

## let me formally ask the question I want to ask
#> is there a statistically signifiacnt difference in a given cluster i in its
#> composition of DEG classes 
#> or, is it better to ask if a cluster is made up 

# first, I just want to see the sums
class_sums = colSums(deg_classes)
cluster_sums = rowSums(deg_classes)

## i think chisq ##
chi_data = data.frame()
for(i in rownames(deg_classes)){
  idx = which(rownames(deg_classes) == i)
  
  cluster_of_interest = deg_classes[i, ]
  other_clusters = colSums(deg_classes[-idx, ])
  
  test_matrix = rbind(cluster_of_interest,
                      other_clusters)
  
  test = chisq.test(test_matrix)
  
  chi_results = data.frame(cluster = i,
                           p_value = test$p.value)
  
  chi_data = rbind(chi_results, chi_data)
  
  
}

chi_data$issignif = ifelse(chi_data$p_value<0.05, '*', NA)

evs = class_sums/27

### plot
library(ggplot2)
plot_data = deg_data_0.05%>%
  group_by(cluster, class)%>%
  summarize(n_degs = n())

ggplot(plot_data, aes(x = as.factor(cluster), y = n_degs, fill = class))+
  geom_bar(stat = 'identity')

