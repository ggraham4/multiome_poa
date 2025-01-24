###Testing Permutation analysis with lasso
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
  library('glmGamPoi')
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
  library(glmnet)
  mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')


}

multiome_obj <- readRDS("~/Desktop/multiome_poa/multiome_poa/RNA Object.rds")

cluster_19_data <- read.csv('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/122324 Neg Bin with Doms/cluster_19.csv')

lasso_scorer_permute <- function(deg_data= cluster_19_data,
                         cluster=19,
                         n_iter = 1:1000
                         ){
  #shutup dplyr
    options(dplyr.summarise.inform = FALSE)

  #subset data to cluster and no nas
  data<- subset(deg_data, cluster == cluster & !is.na(gene))

  #ignore clusters w/o degs
  if(nrow(data)<1){return(NULL)}

  #condirm numeric
  data$d_f_q.value <- as.numeric(data$d_f_q.value)
  data$d_m_q.value <- as.numeric(data$d_m_q.value)
  data$f_m_q.value <- as.numeric(data$f_m_q.value)

  #extract DEGs
  degs <- data$gene[data$f_m_q.value<0.05|
                      data$d_f_q.value<0.05|
                      data$d_m_q.value<0.05]
  degs <- degs[!is.na(degs)]
  
  #if there are none, ignore
  if(length(degs) ==0){return(NULL)}

  #summarize their expression 
degs_expression <- data.frame()
for(i in degs){
  gene_expression <-mean_expression_cluster_data(
    multiome_obj,
    i,
    cluster
    )

  new_data <-data.frame(
    cluster= cluster,
    gene = i,
    mean = gene_expression$mean,
    sex = gene_expression$Sex,
    individual = gene_expression$individual
    )
  degs_expression <- rbind(new_data, degs_expression)
}

### this is where the permutations are going to be ###
#Define sexes as binomial  
degs_expression$fish.class <- ifelse(degs_expression$sex == 'F',1,NA)
degs_expression$fish.class <- ifelse(degs_expression$sex == 'M',0,degs_expression$fish.class)
#remove the other weird sexes
degs_expression <- subset(degs_expression, sex == 'F'|
                            sex=='M'|
                            sex == 'D')

pivoted_data<- degs_expression%>%
  pivot_wider(names_from = gene, 
              values_from = mean)

training_data <- pivoted_data[!is.na(pivoted_data$fish.class),]
training_data <- training_data[complete.cases(training_data[, 5:ncol(training_data)]), ]

x.training <- as.matrix(na.omit(training_data[,5:ncol(training_data)]))
  if(ncol(x.training)<2){return(NULL)}

y.training <- training_data$fish.class

lasso <-suppressWarnings(cv.glmnet(y = y.training, x = x.training, family = "binomial", alpha = 1, lambda = NULL))

min <- lasso$lambda.min

lasso.final <- suppressWarnings(glmnet(x=x.training, y=y.training, alpha = 1, family = "binomial",
                      lambda = min))

test_data <- pivoted_data[is.na(pivoted_data$fish.class),]
test_data <- test_data[complete.cases(test_data[, 5:ncol(test_data)]), ]

x.test <- as.matrix(na.omit(test_data[,5:ncol(test_data)]))
  if(ncol(x.test)<2){return(NULL)}

probabilities <- lasso.final %>% predict(newx = x.test, type = 'response')

predicted_data <- as.data.frame(unname(probabilities))
made.data <- as.data.frame(unname(probabilities))
made.data$fish <- test_data$individual
made.data$probabilities <- probabilities
made.data$sexes = 'real'

real.data <- data.frame(probabilities = mean(probabilities),
                        sexes = 'real')

pivoted_data2 <- pivoted_data
perm_function <- function(arg, 
                          pivoted_data=pivoted_data2){
  message(arg)
  pivoted_data$fake_sex <- sample(pivoted_data$sex, nrow(pivoted_data), replace = F)
 
  pivoted_data$fish.class <- ifelse(pivoted_data$fake_sex == 'F',1,NA)
pivoted_data$fish.class <- ifelse(pivoted_data$fake_sex == 'M',0,pivoted_data$fish.class)
  training_data <- pivoted_data[!is.na(pivoted_data$fish.class),]
  
training_data <- training_data[complete.cases(training_data[, 5:ncol(training_data)-1]), ]

x.training <- as.matrix(na.omit(training_data[,5:ncol(training_data)]))
  if(ncol(x.training)<2){return(NULL)}

y.training <- training_data$fish.class

lasso <-suppressWarnings(cv.glmnet(y = y.training, x = x.training, family = "binomial", alpha = 1, lambda = NULL))

min <- lasso$lambda.min

lasso.final <- suppressWarnings(glmnet(x=x.training, y=y.training, alpha = 1, family = "binomial",
                      lambda = min))

test_data <- pivoted_data[is.na(pivoted_data$fish.class),]
test_data <- test_data[complete.cases(test_data[, 5:ncol(test_data)]), ]

x.test <- as.matrix(na.omit(test_data[,5:ncol(test_data)]))
  if(ncol(x.test)<2){return(NULL)}

probabilities <- lasso.final %>% predict(newx = x.test, type = 'response')

predicted_data <- as.data.frame(unname(probabilities))
made.data <- as.data.frame(mean(unname(probabilities)))
made.data$sexes = 'fake'


return(made.data)
  }

perm_test <- sapply(n_iter, perm_function)
perm_test_p <- as.data.frame(t(perm_test))

colnames(perm_test_p) <- c('probabilities', 'sexes')
perm_test_p$probabilities <- as.numeric(perm_test_p$probabilities)

final_test <- rbind(perm_test_p, real.data)

print(ggplot(final_test, aes(x = probabilities))+
        geom_histogram()+
        geom_vline(data = subset(final_test, sexes == 'real'), aes(xintercept = probabilities), linetype = 2, color = 'red')
)

p.value <- 'none calculated'

if(final_test$probabilities[final_test$sexes=='real']>0.5){
p.value <- sum(final_test$probabilities>final_test$probabilities[final_test$sexes=='real'])/nrow(final_test)
}

if(final_test$probabilities[final_test$sexes=='real']<0.5){
p.value <- sum(final_test$probabilities<final_test$probabilities[final_test$sexes=='real'])/nrow(final_test)
}
cat('\np.value = ',p.value)

return(final_test)
}

clust_19_perms <- lasso_scorer_permute(deg_data = cluster_19_data, 
                                       cluster = 19,
                                       n_iter = 1:5000)


