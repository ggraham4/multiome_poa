#    install.packages("arrow")
library(arrow)
library(tidyverse)
library(ggplot2)

observed_100 = read_parquet("Classifying Newborn Neurons/output_22/2025_10_23_100_restarts_observed.parquet")
#“For each bootstrap iteration, what fraction of newborn neurons were classified as each mature neuronal subtype (observed vs shuffled labels).”

shuffled_100 =read_parquet("Classifying Newborn Neurons/output_22/2025_10_23_100_restarts_shuffled.parquet")
#Same as above, but with shuffled labels (null expectation)

again_observed_100 = read_parquet("Classifying Newborn Neurons/output_22/2025_10_23_again_100_restarts_observed.parquet")
# “Average 5-fold cross-validated classification performance on mature neurons themselves (observed vs shuffled), expressed as % of cells of each class predicted as each class.”


again_shuffled_100 = read_parquet("Classifying Newborn Neurons/output_22/2025_10_23_again_100_restarts_shuffled.parquet")
# Same as above but with shuffled labels (null baseline for training performance)


colMeans(observed_100)%>%sort()
#this is why I shouldve used the subset...

colMeans(again_shuffled_100[,-1])%>%sort()

hist(again_shuffled_100[,-1]%>%as.matrix())
hist(observed_100[,-1]%>%as.matrix())

data_for_plot <- data.frame(percentage = again_shuffled_100[,-1]%>%
                              as.matrix()%>%
                              as.numeric(),
                            label = 'null', 
                            cluster = 'na')
data_for_plot2 = data.frame(percentage = colMeans(observed_100[,-1])%>%
                              as.matrix()%>%
                              as.numeric(),
                            label = 'empirical',
                            cluster = colnames(observed_100[,-1]))

data_for_plot_joined = rbind(data_for_plot, data_for_plot2)
data_for_plot_joined$percentage = as.numeric(data_for_plot_joined$percentage)
data_for_plot_joined$cluster = as.factor(data_for_plot_joined$cluster)
data_for_plot_joined$label = as.factor(data_for_plot_joined$label)

ggplot(data_for_plot, aes(x = percentage))+
  geom_histogram()+
  geom_histogram(data = data_for_plot2, aes(x = percentage, fill = cluster))+
  geom_text(data = data_for_plot2, aes(x = percentage, y =1, label = cluster))

sum(again_shuffled_100[,-1]%>%as.matrix()%>%as.numeric()%>%na.omit()<0.651558074)/length(again_shuffled_100[,-1]%>%as.matrix())
again_shuffled_100[,-1]%>%as.matrix()%>%as.numeric()%>%na.omit()%>%range()

# empirical p values
means = colMeans(observed_100)%>%sort()

# two sided p-value 
(1-sum(again_shuffled_100[,-1]%>%as.matrix()%>%as.numeric()%>%na.omit()< means['3'])/length(again_shuffled_100[,-1]%>%as.matrix()))/2
#0

(1-sum(again_shuffled_100[,-1]%>%as.matrix()%>%as.numeric()%>%na.omit()< means['6'])/length(again_shuffled_100[,-1]%>%as.matrix()))/2
#0

(1-sum(again_shuffled_100[,-1]%>%as.matrix()%>%as.numeric()%>%na.omit()< means['5'])/length(again_shuffled_100[,-1]%>%as.matrix()))/2

(1-sum(again_shuffled_100[,-1]%>%as.matrix()%>%as.numeric()%>%na.omit()< means['21'])/length(again_shuffled_100[,-1]%>%as.matrix()))/2


### These are the only two significant p values
ggplot(data_for_plot, aes(x = percentage))+
  geom_histogram()+
  geom_histogram(data = data_for_plot2, aes(x = percentage, fill = cluster))+
  geom_text(data = data_for_plot2, aes(x = percentage, y =1, label = cluster))+
  labs(x = 'Percentage of Cells in 22 Predicted to Cluster', y = 'Count')



DotPlot(obj, 'neurod1')
