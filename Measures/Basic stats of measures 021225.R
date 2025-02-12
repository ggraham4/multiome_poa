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
  #library(CytoTRACE)
  # SCTRransform_mean_plot <- readRDS("R/Gabe/SCTRransform_mean_plot.rds")
  #mac.neg.bin <- readRDS(file = 'R/Gabe/mac.neg.bin.rds')
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  #mean_cell <- readRDS('R/Gabe/mean_cell.rds')
  library(openxlsx)
  #clown_go <- readRDS('R/Gabe/clown_go.rds')
  
}
multiome_object <-readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

 data.frame(time = multiome_object$Time_Day_2,
                       individual = multiome_object$individual,
                       status = multiome_object$Status)%>%
  group_by(individual, status)%>%
  summarize(time = mean(time))%>%
  ggplot(aes(x = status, y = time))+
  geom_point()+
  geom_boxplot()
 
 data.frame(time = multiome_object$Vo,
            individual = multiome_object$individual,
            status = multiome_object$Status)%>%
   group_by(individual, status)%>%
   summarize(time = mean(time))%>%
   ggplot(aes(x = status, y = time))+
   geom_point()+
   geom_boxplot()
 
 
 data.frame(time = multiome_object$Behaviors_Day_2,
            individual = multiome_object$individual,
            status = multiome_object$Status)%>%
   group_by(individual, status)%>%
   summarize(time = mean(time))%>%
   ggplot(aes(x = status, y = time))+
   geom_point()+
   geom_boxplot()
 
 data.frame(time = multiome_object$Log_11KT,
            individual = multiome_object$individual,
            status = multiome_object$Status)%>%
   group_by(individual, status)%>%
   summarize(time = mean(time))%>%
   ggplot(aes(x = status, y = time))+
   geom_point()+
   geom_boxplot()
 
 gonad_data <- data.frame(Test = multiome_object$Percent_Testicular,
            individual = multiome_object$individual,
            status = multiome_object$Status)%>%
   group_by(individual, status)%>%
   summarize(Test = mean(Test))
 gonad_data_no_e <- subset(gonad_data, status !='E')

 t.test(gonad_data_no_e$Test~gonad_data_no_e$status)
 
 DotPlot(multiome_object, 'vegfd')+
   coord_flip()
 