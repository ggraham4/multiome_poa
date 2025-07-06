library(ggplot2)
library(tidyverse)
library(ggsignif)
library(forcats)

data = read_csv("/Users/ggraham/Desktop/multiome_poa/Sex Classifier and Linearity/progress_classifier/individual_predictions_06_30_2025.csv")

### first, lets group by sex and cluster and gather means

data_grouped_cluster_sex = data%>%
  group_by(cluster, status)%>%
  summarize(mean_pred = mean(prediction),
            se_pred = sd(prediction)/sqrt(n()))

data_grouped_cluster = data%>%
  group_by(cluster)%>%
  summarize(mean_pred = mean(prediction),
            se_pred = sd(prediction)/sqrt(n()))

plot_1 = ggplot(data_grouped_cluster_sex, aes(x = as.factor(cluster), mean_pred, y = mean_pred, shape = status, color =status))+
   #geom_errorbar(aes(x = as.factor(cluster), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred))+
   geom_point(size = 3)+
  scale_shape_manual(values = c(1,2,3))+
  labs(x ='Cluster', y= 'Mean +/- SE Prediction')
  
#the nf seems wrong, I want to look at the raw data 
for(i in 0:25){
  print(i)
  if(i == 21){next}
  data <- read.csv(
    paste0('~/Desktop/multiome_poa/Sex Classifier and Linearity/progress_classifier/degs_0.1_f_m/cluster_',i,'.csv'))
  assign(paste0('cluster_',i,'_degs'), data)
}

## 15 and it only has 2 degs
ggplot(cluster_15_degs, aes(x = slc6a8, y = faua, color = Status, shape = Status))+
  geom_point(size = 4)#+
  geom_smooth(method = 'lm')
## you would definitely say that the NF is closer to the male expression than the female

#what about 7
cluster_7_degs$Status = factor(cluster_7_degs$Status, levels =c('M','D','E','NF','F') )
ggplot(cluster_7_degs, aes(x = Status, y = ppiab))+
  geom_point()
# again, the NF is definitely closer to the male than the female 
# BUT, this is also a clear example of nonlinearity, dominants are higher than females and 
#males

#2 
cluster_2_prcomp = prcomp(cluster_2_degs[,4:ncol(cluster_2_degs)], scale = T)
fviz_pca_ind(cluster_2_prcomp, habillage =cluster_2_degs$Status, addEllipses = T)

#3 
cluster_3_prcomp = prcomp(cluster_3_degs[,4:ncol(cluster_3_degs)], scale = T)
fviz_pca_ind(cluster_3_prcomp, habillage =cluster_3_degs$Status, addEllipses = T)

ggplot(cluster_12_degs, aes(x = hectd2, y = nfixb, color = Status, shape = Status))+
  geom_point(size = 4)#+

#what about 13
cluster_13_degs$Status = factor(cluster_13_degs$Status, levels =c('M','D','E','NF','F') )
ggplot(cluster_13_degs, aes(x = Status, y = tcf7l2))+
  geom_point()



#0 
cluster_0_prcomp = prcomp(cluster_0_degs[,4:ncol(cluster_0_degs)], scale = T)
fviz_pca_ind(cluster_0_prcomp, habillage =cluster_0_degs$Status, addEllipses = T)





