library(ggplot2)
library(ggsignif)
library(tidyverse)

data <- read_csv("Measures/all_data.csv")

data$Status_dummy = data$Status
data$Status_dummy <- factor(data$Status_dummy, levels = c('M',
                                                          'F',
                                                          'D',
                                                          'S',
                                                          'E',
                                                          'EP',
                                                          'NM',
                                                          'NF'
))

data$Status <- ifelse(data$Status %in% c('NM','NF'), NA, data$Status)

### Behaviors day 2 ####
data$Behaviors_Day_2 <- as.numeric(data$Behaviors_Day_2)


beh_model <- lm(Behaviors_Day_2~Status_dummy, data = data)
beh_pairs <- as.data.frame(pairs(emmeans(beh_model, 'Status_dummy'), adjust = 'none'))
beh_pairs$issignif <- ifelse(beh_pairs$p.value<0.05, '*', NA)


behavior_plot <- ggplot(data, aes(x = Status_dummy, y = Behaviors_Day_2))+
  geom_boxplot(data = subset(data, !is.na(Status)),aes(group = Status, fill = Status, , color = Status), alpha = 0.25, outlier.shape = NA)+
  geom_point(size = 2, color = 'black', shape = 1, position = position_dodge2(0.5))+
  theme_classic()+
  geom_signif(xmin = c(1), xmax = c(4), y_position = c(max(data$Behaviors_Day_2)*1.2), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(3), xmax = c(4), y_position = c(max(data$Behaviors_Day_2)*1.05), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(5), xmax = c(6), y_position = c(max(data$Behaviors_Day_2)*1.05), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  labs(x  ='Sex', y = 'Parental Care Behaviors')+
  ylim(c(0, max(data$Behaviors_Day_2)*1.25))+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0), legend.position = 'none')
behavior_plot

#ggsave(plot = behavior_plot,
#       file = "behavior_plot.tiff",
#       device = "tiff",
#       units = "in",
#       width = 2.5,
#       height = 2.5,
#       path = "Bachelors Thesis/Plots/Figure 1")

### Time day 2 #### 
data$Time_Day_2 <- as.numeric(data$Time_Day_2)

time_model <- lm(Time_Day_2~Status_dummy, data = data)
time_pairs <- as.data.frame(pairs(emmeans(time_model, 'Status_dummy'), adjust = 'none'))
time_pairs$issignif <- ifelse(time_pairs$p.value<0.05, '*', NA)
time_pairs

time_plot <- ggplot(data, aes(x = Status_dummy, y = Time_Day_2))+
  geom_boxplot(data = subset(data, !is.na(Status)),aes(group = Status, fill = Status, , color = Status), alpha = 0.25, outlier.shape = NA)+
  geom_point(size = 2, color = 'black', shape = 1, position = position_dodge2(0.5))+
  theme_classic()+
  labs(x  ='Sex', y = 'Time in Nest (s)')+
  geom_signif(xmin = c(1), xmax = c(2), y_position = c(max(data$Time_Day_2)*1.05), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(1), xmax = c(4), y_position = c(max(data$Time_Day_2)*1.2), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(3), xmax = c(4), y_position = c(max(data$Time_Day_2)*1.05), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=5)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0), legend.position = 'none')+
  ylim(c(0, max(data$Time_Day_2)*1.25))
time_plot

#ggsave(plot = time_plot,
#       file = "time_plot.tiff",
#       device = "tiff",
#       units = "in",
#       width = 2.5,
#       height = 2.5,
#       path = "Bachelors Thesis/Plots/Figure 1")

#### Testosterone ####
data$Log_11KT <- as.numeric(data$Log_11KT)

kt_model <- lm(Log_11KT~Status_dummy, data = data)
kt_pairs <- as.data.frame(pairs(emmeans(kt_model, 'Status_dummy'), adjust = 'none'))
kt_pairs$issignif <- ifelse(kt_pairs$p.value<0.05, '*', NA)
kt_pairs

kt_data <- subset(data, !is.na(Log_11KT))
kt_plot <- ggplot(kt_data, aes(x = Status_dummy, y = Log_11KT))+
  geom_boxplot(data = subset(data, !is.na(Status)),aes(group = Status, fill = Status, , color = Status), alpha = 0.25, outlier.shape = NA)+
  geom_point(size = 2, color = 'black', shape = 1, position = position_dodge2(0.5))+
  theme_classic()+
  labs(x  ='Sex', y = 'Log 11KT (pg/ml)')+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(max(kt_data$Log_11KT)*1.05), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(1), xmax = c(4), y_position = c(max(kt_data$Log_11KT)*1.5), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(1), xmax = c(6), y_position = c(max(kt_data$Log_11KT)*1.7), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
  geom_signif(xmin = c(2), xmax = c(4), y_position = c(max(kt_data$Log_11KT)*1.3), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
   geom_signif(xmin = c(2.1), xmax = c(3), y_position = c(max(kt_data$Log_11KT)*1.05), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=5)+
   theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0), legend.position = 'none')+
  ylim(c(0, max(kt_data$Log_11KT)*1.8))
kt_plot

ggsave(plot = kt_plot,
       file = "kt_plot.tiff",
       device = "tiff",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Bachelors Thesis/Plots/Figure 1")


