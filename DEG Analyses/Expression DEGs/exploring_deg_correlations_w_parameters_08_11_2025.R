{
  library(ggridges)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(CytoTRACE)
  library(BiocManager)
  library(lme4)
  library(enrichR)
  `%notin%` = Negate(`%in%`)
  library(car)
  library(emmeans)
  
  clown_go = readRDS("Functions/clown_go2")  
  mecd = readRDS("Functions/mean_expression_cluster_data.rds")
  mecp = readRDS("Functions/mean_expression_cluster_plot.rds")
  }

#define functions
individual_mean_expression = function(object, gene, clustering= 'sub.cluster'){
  library(stringr)    
  options(dplyr.summarise.inform = FALSE)
  
  counts <- t(as.matrix(object@assays$RNA$data))
  Counts_of_interest <- as.data.frame(counts[,gene])
  Counts_of_interest$expression <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual
  Counts_of_interest$Status <- object@meta.data$Status
  Counts_of_interest$cluster <- object@meta.data[[clustering]]
  
  results <-Counts_of_interest%>%
    group_by(individual, Status, cluster)%>%
    summarize(mean = mean(expression),
              se = sd(expression)/sqrt(n()))
  results$Sex <- results$Status
  
  
  results$Sex <- str_sub(results$individual, -1)
  results$Sex[results$individual == 'T17D'] = 'NF'
  results$Sex[results$individual == 'A12D'] = 'E'
  results$Sex[results$individual == 'T11D'] = 'E'
  results$Sex[results$individual == 'GH'] = 'NRM'
  return(results)
  
  
}

plot_mean_expression <- function(data_frame){
  data_frame$Sex = factor(data_frame$Sex, levels = c('NRM','M', 'D','E','NF',"F"))
  plot <- ggplot(data_frame, aes(x = cluster, y = mean, color = Sex))+
    geom_point(position = position_dodge(0.75))+
    geom_boxplot(alpha = 0)+
    theme_minimal()+
    labs(y = 'Mean Normalized Expression', x ='Sub Cluster')
  
  return(plot)
  
  
}

#read in object
obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")


mecp(obj, 
     'esr2b',
     9,
     'res0.8_50nn_40PC_45LSI')

mecp(obj, 
     'LOC111577263',
     1,
     'res0.8_50nn_40PC_45LSI')

mecp(obj, 
     'pgr',
     6,
     'res0.8_50nn_40PC_45LSI')


mecp(obj, 
     'LOC111568069',
     6,
     'res0.8_50nn_40PC_45LSI')

esr_data = mecd(obj, 
     'esr2b',
     9,
     'res0.8_50nn_40PC_45LSI')

kt_data = obj@meta.data%>%
  group_by(individual)%>%
  summarize(mean = mean(Log_11KT))

esr_data$test = kt_data$mean

t=cor.test(esr_data$mean, esr_data$test)# significant
plot(esr_data$mean, esr_data$test) # meh r2

cor_and_plot = function(gene, cluster, parameter){
  # Get the parameter name as a string
  param_name <- deparse(substitute(parameter))
  
  gene_mean = mecd(obj, 
                   gene,
                   cluster,
                   'res0.8_50nn_40PC_45LSI')
  
  parameter_mean = obj@meta.data %>%
    group_by(individual) %>%
    summarize(mean_ = mean({{parameter}}))
  
  gene_mean$param = parameter_mean$mean_
  
  test = cor.test(gene_mean$mean, gene_mean$param)
  cat('\n P.value = ',test$p.value, '\n')
  cat('\n R = ',test$estimate, '\n')
  
  plot <- ggplot(gene_mean, aes(x = mean, y = param, color = Status))+
    geom_point()+
    geom_smooth(method = 'lm')+
    labs(x = paste0(gene,'_',cluster), y = param_name)
  print(plot)
}
cor_and_plot('LOC111577263', 1, Percent_Testicular)

#the problem with this is I am trying to do comparisons across sex when maybe I shouldnt
cover = function(gene, cluster, parameter){
  gene_mean = mecd(obj, 
                   gene,
                   cluster,
                   'res0.8_50nn_40PC_45LSI')
  
  parameter_mean = obj@meta.data%>%
    group_by(individual)%>%
    summarize(mean_ = mean({{parameter}}))
  
  gene_mean$param = parameter_mean$mean_
  
  test = lm(mean~param+Status, data = gene_mean)
return(anova(test, test = 'Chisq'))
  
}


### ok so I want to discuss the two kinds of sex change I think are happening
#> 1) correlated with gonads
#> 2) independent of gonads
#> 
#> 1 I think is pgr-like in 6 and esr2b in 9
#> 2 I think is brain aromatase in 1 and pgr in 6
#> 
#> Gonadal correlated
#esr2b, 9
cover('esr2b', 9, Percent_Testicular) #*
cor_and_plot('esr2b', 9, Percent_Testicular) #seems negative
cover('esr2b', 9, Log_11KT) #*
cor_and_plot('esr2b', 9, Log_11KT) #certainly no convincing correaltion, dominant and female on top of eachother

#pgrlike 6
cover('LOC111568069', 6, Percent_Testicular) #** 
cor_and_plot('LOC111568069', 6, Percent_Testicular) #divergent, interesting
cover('LOC111568069', 6, Log_11KT) #***
cor_and_plot('LOC111568069', 6, Log_11KT) #not divergent, d and f basically on top of eachother suggesting this change occurs prior to gonads?

# Gonadal independent
cover('pgr', 6, Percent_Testicular) #na
cor_and_plot('pgr', 6, Percent_Testicular) 
cover('pgr', 6, Log_11KT) #*
cor_and_plot('pgr', 6, Log_11KT) #positive, same as below

cover('LOC111577263', 1, Percent_Testicular) #*
cor_and_plot('LOC111577263', 1, Percent_Testicular) 
cover('LOC111577263', 1, Log_11KT) #*
cor_and_plot('LOC111577263', 1, Log_11KT) # positive, dominant and MALE more or less continuous suggesting this occurs after gonads?

#> ok this is interesting, so I think I can say that the ones i am calling gonadal correlated, 
#> the adoption of the female phenotype occurs during gonadal transformation or prior as dominants and females form a continuum distinct from males, thus changes 
#> occur prior to  6mos
#> 
#> the ones that I am calling gonadal independent seem to only change after gonads/hormones have changed as dominants still look male\
#> thus adoption of female phenotype occurs after 6 most
#> 
#raw params
param_data = read.csv("Measures/all_data.csv")

#to test my belief, I am going to do the same analysis but with behavior, which I know changes last
# my predictio n is that males and dominants will form a continuum and females will be distinct
ggplot(param_data, aes(x = Behaviors_Day_2, y = Log_11KT, color = Status))+
  geom_point()+
  geom_smooth(method = 'lm')
# and my prediction was right,


cor.test(param_data$Behaviors_Day_2, param_data$Log_11KT) #highly correlated
cor.test(param_data$Behaviors_Day_2, param_data$Percent_Testicular) # correlated
#so behavior (in my opinion) comes from the hormones thus I dont need to talk about how esr2b is correlated with
#> behavior because that is autocorrelated with the physiology
#> 


