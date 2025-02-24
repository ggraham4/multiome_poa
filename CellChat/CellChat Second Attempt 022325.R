#Trying out cellchat
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
  library(glmGamPoi)
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
  library(CytoTRACE)
  library(ggrepel)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(Polychrome)
  library(scCustomize)
  library(CellChat)
  
  
  P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
  swatch(P40)
  names(P40) <- NULL
  
  mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
  prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
  define_degs_prop<- readRDS('Functions/define_degs_prop.rds')
  mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
  prop_deg_function.rds<- readRDS('Functions/DEG_functions/prop_deg_function.rds')
  define_behavior_degs<- readRDS('Functions/define_behavior_degs')
  clown_go<- readRDS('Functions/clown_go')
  define_degs<- readRDS('Functions/define_degs')
  
}

#read in seurat object
obj <- readRDS('C:/Users/Gabe/Desktop/RNA object mouse names.rds')

#for some reason I need to do this because it doesnt tolerate 0s in cluster labels even if they are characters ??
obj$cluster_labels <-as.character(obj$harmony.wnn_res0.4_clusters)
obj$cluster_labels <- gsub("\\b0\\b", "o", obj$cluster_labels)
obj$cluster_labels <- paste0('Cluster ',obj$cluster_labels)


cellchat <- createCellChat(object = obj, 
                           #meta = meta, 
                           group.by = "cluster_labels")

cellchat <- setIdent(cellchat, ident.use = "cluster_labels") # set "labels" as default cell identity


# call up cellchat database
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling

# set the used database in the object
cellchat@DB <- CellChatDB.use


# need to run this subsetData command then preprocessing continues
cellchat <- subsetData(cellchat) # not sure what this does

# identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat) # this runs once for each cluster
cellchat <- identifyOverExpressedInteractions(cellchat) # this runs quickly

# compute communication probability - slow step, about ten minutes
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# the subsetCommunication() function could be useful at this point
# it is used to access inferred cell-cell communications of interest
# like subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
# which returns the inferred signaling from 1&2 to 4&5

# compute communication at a pathway level - much faster step
cellchat <- computeCommunProbPathway(cellchat)

# compute aggregated network 
cellchat <- aggregateNet(cellchat)

##I am going to save this and pick it up tomorrow 

saveRDS(cellchat, 'C:/Users/Gabe/Desktop/Local R Files/cell chat second attempt 022325.rds')

cellchat <- readRDS('C:/Users/Gabe/Desktop/Local R Files/cell chat second attempt 022325.rds')

####################
cellchat@net$pval

df.net <- subsetCommunication(cellchat)
df.net$adjusted_p_value = p.adjust(df.net$pval, 'fdr', nrow(df.net))
df.net_sig <- subset(df.net, adjusted_p_value<0.05 & receptor != 'Insr')

### which cluster has the most enriched inputs
table(df.net_sig$target)%>%sort()
#cluster 27 is easily the most targeted

table(df.net_sig$source)%>%sort()
#and the biggest sender
#What is cluster 27 in the old dataset?

###also I want to plot this
##fix labels
idents = as.array(as.character(cellchat@idents))
idents2 = str_remove(celchat@idents, 'Cluster')

cellchat <- updateClusterLabels(
  cellchat,
  old.cluster.name = cellchat@idents,
  new.cluster.name = idents2,
  new.cluster.metaname = "new.labels"
)

unique(df.net$pathway_name)

netVisual_aggregate(cellchat, "CCK", layout = "circle")
title('CCK')

netVisual_aggregate(cellchat, "WNT", layout = "circle")
title('WNT')

netVisual_aggregate(cellchat, "IGF", layout = "circle")
title('WNT')

netVisual_aggregate(cellchat, "COMPLEMENT", layout = "circle")
title('COMPLETMENT')

netVisual_aggregate(cellchat, "ncWNT", layout = "circle")
title('ncWNT')

netVisual_aggregate(cellchat, "EDN", layout = "circle")
title('EDN')

netVisual_aggregate(cellchat, "OPIOID", layout = "circle")
titsle('OPIOID')

netVisual_aggregate(cellchat, "INSULIN", layout = "circle")
title('INSULIN')

netVisual_aggregate(cellchat, "SOMATOSTATIN", layout = "circle")
title('SOMATOSTATIN')

netVisual_aggregate(cellchat, "SEMA3", layout = "circle")
title('SEMA3')

netVisual_aggregate(cellchat, "GCG", layout = "circle")
title('GCG')

netVisual_aggregate(cellchat, "OSTN", layout = "circle")
title('OSTN')

cluster_27_cellchat <- df.net_sig[df.net_sig$target=="Cluster 27" & df.net_sig$pathway_name=='WNT',]


