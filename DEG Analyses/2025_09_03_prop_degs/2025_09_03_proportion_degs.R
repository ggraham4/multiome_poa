library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(lme4)
library(emmeans)
library(parallel)

clown_go = readRDS('Functions/clown_go2')

object = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

proportionDEGs <- function(object, 
                           cluster=6,
                           clustering ='res0.8_50nn_40PC_45LSI',
                           alpha=0.05){
  
  dataMatrix = as.matrix(object@assays$RNA$data[,object@meta.data[[clustering]]==cluster])%>%t()
  binaryMatrix = (dataMatrix>0)%>%as.data.frame()
  binaryMatrix$individual = object$individual[object@meta.data[[clustering]]==cluster]
  binaryMatrix$Status = object$Status[object@meta.data[[clustering]]==cluster]

  glmerDataMaker <- function(gene){
   subsetData = binaryMatrix[,c(gene, 'Status', 'individual')] 
   subsetData$gene = subsetData[,1]
   
  options(dplyr.summarise.inform = FALSE)
  
   tableData = subsetData%>%
     group_by(individual, Status)%>%
     subset(Status %in% c('M','D','F'))%>%
     summarize(cellsExpressing = sum(gene ==T),
               cellsFailing = sum(gene ==F))
   return(tableData)
  }
  
  glmerAnalyzer <- function(gene){
    message(paste0('Gene ', which(genesList==gene), ' of ', length(genesList)))
    dataNow <- glmerDataMaker(gene)
    glmerMat <- as.matrix(cbind(dataNow[,3], dataNow[,4]))
    
    if(sum(glmerMat)==0){return(NULL)}
    if(sum(glmerMat[,1])==0){return(NULL)}
    
    model = tryCatch({
      glmer(glmerMat~Status + (1|individual), 
            family = binomial(),
            dataNow,
            control = glmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4)))
    }, error = function(e) {
      return(NULL)
    })
    
    if(is.null(model)){return(NULL)}
    
    type_iii = car::Anova(model, type ='III')
    warning = ifelse(length(model@optinfo$conv$lme4$code) != 0, substr(model@optinfo$conv$lme4$messages, 1, 50), NA)

    newd = data.frame(gene =gene,
                      av_p.value = type_iii$`Pr(>Chisq)`[2],
                      singular = isSingular(model),
                      warning = warning)
    
    return(newd)
  }
  
  max = ncol(object@assays$RNA$data[,object@meta.data[[clustering]]==cluster])
  genesList = rownames(object@assays$RNA$data[,object@meta.data[[clustering]]==cluster])
  
  anovas = mclapply(genesList, glmerAnalyzer, mc.cores = detectCores()-1)
  anovasAnalyzed = do.call(rbind, anovas)
  anovasAnalyzed$av_q.value = p.adjust(anovasAnalyzed$av_p.value, 'fdr', nrow(anovasAnalyzed))
  
  anovaSignificant=NULL
  
  anovaSignificant = subset(anovasAnalyzed, av_q.value<alpha & is.na(warning))
  
  if(nrow(anovaSignificant)<1){return(NULL)}
  if(is.null(anovaSignificant) ){return(NULL)}
  
  significantGenes = anovaSignificant$gene
  
  anovaSignificant$d_m_p.value = NA
  anovaSignificant$d_f_p.value = NA
  anovaSignificant$f_m_p.value = NA
  
  anovaSignificant$d_m_estimate = NA
  anovaSignificant$d_f_estimate = NA
  anovaSignificant$f_m_estimate = NA


  pairwiseAnalyzer <- function(gene){
    emptyDF = data.frame()
    for(gene in significantGenes){
    message(paste0('Gene ', which(significantGenes==gene), ' of ', length(significantGenes)))
    dataNow <- glmerDataMaker(gene)
    glmerMat <- as.matrix(cbind(dataNow[,3], dataNow[,4]))
    
    if(sum(glmerMat)==0){return(NULL)}
    if(sum(glmerMat[,1])==0){return(NULL)}
    
    model = tryCatch({
      glmer(glmerMat~Status + (1|individual), 
            family = binomial(),
            dataNow,
            control = glmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4)))
    }, error = function(e) {
      return(NULL)
    })
    
    if(is.null(model)){return(NULL)}
    
    pairwiseComparisons = pairs(emmeans(model,'Status'), adjust = 'none' )%>%as.data.frame()
    
    pairsDF = data.frame(gene = gene,
                        d_m_estimate= pairwiseComparisons$estimate[pairwiseComparisons$contrast == 'D - M'],
                         d_m_p.value = pairwiseComparisons$p.value[pairwiseComparisons$contrast == 'D - M'],
                        f_m_estimate= pairwiseComparisons$estimate[pairwiseComparisons$contrast == 'F - M'],
                         f_m_p.value = pairwiseComparisons$p.value[pairwiseComparisons$contrast == 'F - M'],
                         d_f_estimate= pairwiseComparisons$estimate[pairwiseComparisons$contrast == 'D - F'],
                         d_f_p.value = pairwiseComparisons$p.value[pairwiseComparisons$contrast == 'D - F']
                         )
    emptyDF = rbind(emptyDF, pairsDF)
    }
    return(emptyDF)
    
  }
  if(length(significantGenes) < 1 | all(is.na(significantGenes))){return(NULL)}
  
  pairwiseComparisons = pairwiseAnalyzer(significantGenes)
  
  final_list = anovasAnalyzed%>%
    subset(is.na(warning))%>%
    full_join(pairwiseComparisons, by = 'gene')
    
  df_sorted <- final_list[order(final_list$av_q.value), ]  

  return(df_sorted)
}

degs_6 = proportionDEGs(object)

clown_go(degs_6$gene[degs_6$av_q.value<0.05 & degs_6$singular==F])%>%dotplot()


for(i in 16:26){
  start = Sys.time()
  degs =  proportionDEGs(object,
                         i)
  
  write.csv(degs, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2025_09_04_proportion_DEGs/cluster_',i,'.csv'))
  end = Sys.time()
  print(end-start)
}
