### cell chat analysis

library(Seurat)
library(CellChat)
library(tidyverse)
library(lme4)
library(patchwork)
library(parallel)
library(future)
options(stringsAsFactors = FALSE)
 
#read in data       
human_named = readRDS("A:/optimal_clustering_05_06_2025/RNA_object_human_names.rds")

# change cluster names cause IG cluster 0 isnt allowed
meta = human_named@meta.data
human_named$renamed = as.character(human_named$res0.8_50nn_40PC_45LSI)
human_named$renamed = ifelse(human_named$renamed ==0, 'o',human_named$renamed)

#use all cells, make meta data ====
cell.use = rownames(meta)
meta = meta[cell.use, ]

# set up for cell chat ====
cellchat= createCellChat(object = human_named, meta = human_named@meta.data, group.by = 'renamed')

#add meta data ====
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "renamed") 
groupSize <- as.numeric(table(cellchat@idents))

# set interaction db ====
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
"Users can update CellChatDB by adding their own curated ligand-receptor pairs.Please check our tutorial on how to do it.
"
#may be worth doing to add DHEA if it is not included
# a lot from cellphone is not included that I would like to add

#use all interaction ----
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use

# preprocessing ----
# raw use = true, my data is of good quality
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multicore", workers = detectCores()-1) 
options(future.globals.maxSize = 1000 * 1024^10)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# inference of comm network ----
#for now, I am going to use raw use, though can be more conservative later, though O
cellchat <- computeCommunProb(cellchat, raw.use =T) 
#examine
df.net <- subsetCommunication(cellchat)
#man this misses a lot of interactions

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#saveRDS(cellchat,'A:/2025_09_26_cellchat_allclusters_raw_use.rds')

### analyze ----
# there is no way to do cellchat by sex after the fact, so I will have to manually do it ig
rm(cellChat)
big_cellchat_list= list()
individuals = unique(human_named$individual)
for(ind in individuals){
  print(ind)
  
  human_named_temp = subset(human_named, individual== ind)
  s = unique(human_named_temp$Status)
  
  cellchat <- createCellChat(object = human_named_temp, group.by = "renamed", assay = "RNA")
  
  rm(human_named_temp)
  CellChatDB <- CellChatDB.human 

  CellChatDB.use <- CellChatDB 
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multicore", workers = detectCores()-1) 
  options(future.globals.maxSize = 1000 * 1024^10)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat, raw.use =T,
                                type = "truncatedMean", trim = 0.05# most liberal
                                ) 
  df.net <- subsetCommunication(cellchat)

  cellchat <- aggregateNet(cellchat)
    
  big_cellchat_list[[ind]]= cellchat
  
}


cellchat_merged <- mergeCellChat(big_cellchat_list, add.names = names(big_cellchat_list), cell.prefix = TRUE)

saveRDS(cellchat_merged, 'A:/2025_09_26_cellchat_all/cellchat_by_ind.rds')
cellchat_merged = readRDS( 'A:/2025_09_26_cellchat_all/cellchat_by_ind.rds')
## do statistics

compareInteractions(cellchat_merged, 'weight')
compareInteractions(cellchat_merged, 'count')

# Function to extract pathway data for a single pathway
# Function to extract pathway data for a single pathway (modified to use @net instead of @netP)
extract_pathway_data <- function(interaction, cellchat_merged, clusters) {
  
  pathway_data = data.frame()
  
  for(source_cluster in clusters){
    for(target_cluster in clusters){
      
      interaction_data = data.frame()
      
      for(ind in names(cellchat_merged@net)){
        
        # Check if this individual has any interaction data
        if(is.null(cellchat_merged@net[[ind]]$prob) || length(dim(cellchat_merged@net[[ind]]$prob)) < 3) {
          prob_value = 0
        } else {
          # Extract probability for this specific interaction
          prob_matrix = cellchat_merged@net[[ind]]$prob
          pathway_names = dimnames(prob_matrix)[[3]]
          
          # Check if interaction exists in this individual's data
          if(is.null(pathway_names) || !interaction %in% pathway_names || 
             !source_cluster %in% dimnames(prob_matrix)[[1]] || 
             !target_cluster %in% dimnames(prob_matrix)[[2]]){
            
            prob_value = 0
            
          } else {
            prob_value = prob_matrix[source_cluster, target_cluster, interaction]
          }
        }
        
        temp_row = data.frame(
          individual = ind,
          pathway = interaction,
          source_cluster = source_cluster,
          target_cluster = target_cluster,
          probability = prob_value,
          stringsAsFactors = FALSE
        )
        
        interaction_data = rbind(interaction_data, temp_row)
      }
      
      pathway_data = rbind(pathway_data, interaction_data)
    }
  }
  
  return(pathway_data)
}

# Run in parallel (modified to use @net instead of @netP)
# Extract the probability matrices first to avoid serialization issues
ind = 'T11D'
pathways = cellchat_merged@net[[ind]]$prob[1,,] %>% colnames()

# Extract all probability matrices and metadata upfront
prob_matrices = list()
individual_names = names(cellchat_merged@net)

for(ind_name in individual_names){
  if(!is.null(cellchat_merged@net[[ind_name]]$prob) && length(dim(cellchat_merged@net[[ind_name]]$prob)) >= 3){
    prob_matrices[[ind_name]] = cellchat_merged@net[[ind_name]]$prob
  }
}

# Get clusters from the first available matrix
clusters <- rownames(prob_matrices[[1]])

# Modified function that works with pre-extracted matrices
extract_pathway_data_fixed <- function(interaction, prob_matrices, clusters) {
  
  pathway_data = data.frame()
  
  for(source_cluster in clusters){
    for(target_cluster in clusters){
      
      for(ind in names(prob_matrices)){
        
        prob_matrix = prob_matrices[[ind]]
        pathway_names = dimnames(prob_matrix)[[3]]
        
        # Check if interaction exists in this individual's data
        if(is.null(pathway_names) || !interaction %in% pathway_names || 
           !source_cluster %in% dimnames(prob_matrix)[[1]] || 
           !target_cluster %in% dimnames(prob_matrix)[[2]]){
          
          prob_value = 0
          
        } else {
          prob_value = prob_matrix[source_cluster, target_cluster, interaction]
        }
        
        temp_row = data.frame(
          individual = ind,
          pathway = interaction,
          source_cluster = source_cluster,
          target_cluster = target_cluster,
          probability = prob_value,
          stringsAsFactors = FALSE
        )
        
        pathway_data = rbind(pathway_data, temp_row)
      }
    }
  }
  
  return(pathway_data)
}

# Run in parallel with extracted data
library(parallelsugar)
pathway_comparison_data_list <- mclapply(pathways, 
                                         extract_pathway_data_fixed, 
                                         prob_matrices = prob_matrices,
                                         clusters = clusters,
                                         mc.cores = detectCores()-1)
pathway_comparison_data <- do.call(rbind, pathway_comparison_data_list)

#stats time
library(readxl)
measures = read_excel("Reference/Complete Data Frame (Hormones, Behavior, Size, Gonads) GG.xlsx")
measures_selected = measures %>%
  select(Fish, Status)

pathway_comparison_data = pathway_comparison_data%>%
  right_join(measures_selected, by= join_by('individual'=='Fish'))

pathway_comparison_data = pathway_comparison_data%>%
  subset(Status %in% c('M','D','F'))

pathways = (unique(pathway_comparison_data$pathway))
pathways_len = length(unique(pathway_comparison_data$pathway))
library(emmeans)
out_data = data.frame()
  for(pathway in unique(pathway_comparison_data$pathway)){
    message('Pathway ', which(pathways == pathway), ' of ',pathways_len)
    pathway_data = data.frame()
    for(sender in unique(pathway_comparison_data$source_cluster)){
      for(receiver in unique(pathway_comparison_data$target_cluster)){
        dat = pathway_comparison_data[pathway_comparison_data$pathway  == pathway  &
                                               pathway_comparison_data$ source_cluster == sender&
                                               pathway_comparison_data$target_cluster == receiver,]
        
        if(length(unique(dat$probability))==1){next}
        test = lm(probability ~Status, data = dat)
        av_ = anova(test, test = 'Chisq')
        av_p.value = av_$`Pr(>F)`[1]
        
        pair= pairs(emmeans(test, "Status"), adjust = 'none')%>%as.data.frame()
        d_m_estimate = pair$estimate[pair$contrast== 'D - M']
        d_f_estimate = pair$estimate[pair$contrast== 'D - F']
        f_m_estimate = pair$estimate[pair$contrast== 'F - M']
        
        d_m_p.value = pair$p.value[pair$contrast== 'D - M']
        d_f_p.value = pair$p.value[pair$contrast== 'D - F']
        f_m_p.value = pair$p.value[pair$contrast== 'F - M']
        
        newd = data.frame(pathway = pathway,
                          sender = sender,
                          receiver = receiver,
                          av_p.value = av_p.value,
                          d_m_estimate= d_m_estimate,
                          d_f_estimate = d_f_estimate,
                          f_m_estimate = f_m_estimate,
                          d_m_p.value= d_m_p.value,
                          d_f_p.value = d_f_p.value,
                          f_m_p.value = f_m_p.value
                          )
        pathway_data = rbind(pathway_data, newd)
        
      }
    }
    pathway_data$av_q.value = p.adjust(pathway_data$av_p.value, 'fdr',nrow(pathway_data))
    out_data = rbind(pathway_data, out_data)
  }

out_dat2 = data.frame()
for(sender_val in unique(pathway_comparison_data$source_cluster)){
  for(receiver_val in unique(pathway_comparison_data$target_cluster)){
    temp = subset(out_data, sender == sender_val & receiver == receiver_val)
    temp$cluster_pair_q.value <- p.adjust(temp$av_p.value, 'fdr', nrow(temp))
    out_dat2 = rbind(temp, out_dat2)
  }
}

write.csv(out_data, 'A:/2025_09_26_cellchat_all/pathway_data.csv')
write.csv(out_dat2, 'A:/2025_09_26_cellchat_all/pathway_data_cluster_q.value.csv')

