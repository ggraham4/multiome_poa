library(ggplot2)
library(ggsignif)
library(tidyverse)
library(emmeans)
library(Seurat)

### DIY WGCNA using genes of interest

### read in object
obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

###extract data from a few genes of interest

#hsd3b1 catalyzes androstenediol to testosterone and dhea to androstenedione
hsd3b1_expression = obj@assays$RNA$data['hsd3b1',]

#hsd17b14 catalyzes androstenediol to DHEA and estradiol to estrone
hsd17b14_expression = obj@assays$RNA$data['hsd17b14',]

#hsd17b1 catalyzes DHEA to androstenediol and estrone to estradione and DHT to 3 beta diol
hsd17b1_expression = obj@assays$RNA$data['hsd17b1',]

#hsd17b2 catalyzes test to androstenedione and estrone to estradiol, and androstenediol to DHEA 
#hsd17b2_expression = obj@assays$RNA$data['hsd17b2',]
#i cant find this one, moving on

#cyp19a1b catalyzes androgens to estrogens
cyp19a1b_expression = obj@assays$RNA$data['LOC111577263',]

#hsd17b5_expression = obj@assays$RNA$data['hsd17b5',]

### these data are not normally distributed, so we use spearman correlation 

#first step is HSD3B1 vs HSD17B1, lets pick based on which has the higher expresison 
mean(hsd17b1_expression)# 1.518279e-05
mean(hsd3b1_expression)#0.01512033 #ok 3b1 wins, next

#second step is EITHER HSD17B5 or cyp19a1b
# i cant find hsd17b5, move on 
#getting confused







