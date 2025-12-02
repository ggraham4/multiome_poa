##### hodgkin huxley predictor #####
#read in stuff
library(Seurat)
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

obj_neurons = subset(obj, final_clusters %notin% c(1,
                                                   2,
                                                   20,
                                                   15,
                                                   26,
                                                   11,
                                                   13
                                                   ))

reference =read.csv('Reference/genes updated.csv')
expr_z = obj_neurons@assays$RNA$data%>%as.matrix()%>%t()
extract_pc1 = function(family_genes){
  
  pca <- prcomp(expr_z[,family_genes], center=FALSE, scale.=FALSE)
  E_pc1 <- pca$x[,1]
  return(E_pc1)
}

### Step 1, find the mean expression of the components of a hodgkin huxley neuron ####

#Find genes encoding voltage gated  channels
voltage_gated_k_genes = reference$NIH_accession[intersect(grep('potassium',reference$NIH_description),
                                                  grep('voltage',reference$NIH_description) )]

cell_wise_voltage_gated_k = extract_pc1(voltage_gated_k_genes)

voltage_gated_na_genes = reference$NIH_accession[intersect(grep('sodium',reference$NIH_description),
                                                  grep('voltage',reference$NIH_description) )]

cell_wise_voltage_gated_na = extract_pc1(voltage_gated_na_genes)


voltage_gated_cl_genes = reference$NIH_accession[intersect(grep('chloride',reference$NIH_description),
                                                  grep('voltage',reference$NIH_description) )]
#cell wise
cell_wise_voltage_gated_cl = extract_pc1(voltage_gated_cl_genes)

k_leak_channels <- reference$NIH_accession[
  grepl("kcnk", reference$NIH_accession, ignore.case = TRUE) |
  grepl("two pore", reference$NIH_description, ignore.case = TRUE) |
  grepl("k2p", reference$NIH_description, ignore.case = TRUE) 
]

na_leak_channels =   reference$NIH_accession[grepl("nalcn", reference$NIH_accession, ignore.case = TRUE)]

cell_wise_k_leak = extract_pc1(k_leak_channels)
cell_wise_na_leak =  extract_pc1(na_leak_channels)
# convert to z score
z_k   <- scale(cell_wise_voltage_gated_k)[,1]
z_na  <- scale(cell_wise_voltage_gated_na)[,1]
z_cl  <- scale(cell_wise_voltage_gated_cl)[,1]
z_k_leak<- scale(cell_wise_k_leak)[,1]
z_na_leak =  scale(cell_wise_na_leak)[,1]
# map to multiplier
map_to_multiplier <- function(z, beta = 0.25, m_min = 0.2, m_max = 5){
  m <- 10^(beta * z)
  m[m < m_min] <- m_min
  m[m > m_max] <- m_max
  return(m)
}

g_k_mult    <- map_to_multiplier(z_k)
g_na_mult   <- map_to_multiplier(z_na)
g_cl_mult   <- map_to_multiplier(z_cl)
g_leak_k_mult <- map_to_multiplier(z_k_leak)
g_leak_na_mult =  map_to_multiplier(z_na_leak)

obj_neurons$g_k_mult    <- g_k_mult
obj_neurons$g_na_mult   <- g_na_mult
obj_neurons$g_cl_mult   <- g_cl_mult
obj_neurons$g_leak_k_mult <- g_leak_k_mult
obj_neurons$g_leak_na_mult <- g_leak_na_mult

# solving time - chatgpt code #####
library(deSolve)
hh_model <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # classic gates
    am <- 0.1 * (V+40) / (1 - exp(-(V+40)/10))
    bm <- 4 * exp(-(V+65)/18)
    ah <- 0.07 * exp(-(V+65)/20)
    bh <- 1 / (1 + exp(-(V+35)/10))
    an <- 0.01*(V+55)/(1-exp(-(V+55)/10))
    bn <- 0.125*exp(-(V+65)/80)

    dmdt <- am*(1-m) - bm*m
    dhdt <- ah*(1-h) - bh*h
    dndt <- an*(1-n) - bn*n

    # mult-adjusted conductances
    gNa_eff   <- gNa_base  * g_na_mult
    gK_eff    <- gK_base   * g_k_mult
    gCl_eff   <- gCl_base  * g_cl_mult
    gK_leak   <- gL_K_base * g_leak_k_mult
    gNa_leak  <- gL_Na_base * g_leak_na_mult

    INa      <- gNa_eff  * m^3 * h * (V - ENa)
    IK       <- gK_eff   * n^4 * (V - EK)
    ICl      <- gCl_eff  * (V - ECl)
    IK_leak  <- gK_leak  * (V - EK)
    INa_leak <- gNa_leak * (V - ENa)

    dVdt <- (Iext - INa - IK - ICl - IK_leak - INa_leak) / Cm

    list(c(dVdt, dmdt, dhdt, dndt))
  })
}

simulate_cell <- function(i, meta) {
  parms <- list(
    g_na_mult   = meta$g_na_mult[i],
    g_k_mult    = meta$g_k_mult[i],
    g_cl_mult   = meta$g_cl_mult[i],
    g_leak_k_mult = meta$g_leak_k_mult[i],
    g_leak_na_mult = meta$g_leak_na_mult[i],
    gNa_base=120, gK_base=36, gCl_base=0.3,
    gL_K_base=0.1, gL_Na_base=0.05,
    ENa=50, EK=-77, ECl=-54.4,
    Iext = 10, Cm=1
  )

  init <- c(V=-65, m=0.05, h=0.6, n=0.32)

  times <- seq(0, 50, by=0.01)

  out <- ode(y=init, times=times, func=hh_model, parms=parms)

  # Return AP metrics (peak, width, #spikes, etc)
  V <- out[,"V"]
  n_spikes <- sum(diff(V > 0) == 1)

  return(n_spikes)
}

library(future.apply)
plan(multisession, workers = 6)

results <- future_sapply(seq_len(nrow(obj_neurons)), simulate_cell, meta=obj_neurons@meta.data)
obj_neurons$predicted_spikes <- results

DotPlot(obj_neurons, 'predicted_spikes')


sex_spike_plotter <- function(cluster){

    plot_data = obj_neurons@meta.data %>%
    group_by(individual, Status, final_clusters) %>%
      subset(final_clusters == cluster)%>%
    summarize(mean_spikes = mean(predicted_spikes), .groups = 'drop') 
  
    plot_data$Status = factor(  plot_data$Status, levels = c('NRM', 'M', 'D',"E",'NF','F'))
  # Create plot
  p = ggplot(plot_data, aes(x = Status, y = mean_spikes)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(.5)) +
    theme_minimal()
  
  return(p)
  
  
  
}
sex_spike_plotter(10)+
  labs(title = "Predicted spikes cluster 10")
sex_spike_plotter(6)

sex_spike_plotter(22)+
    labs(title = "Predicted spikes cluster 22")

sex_spike_plotter(9)
sex_spike_plotter(3)
sex_spike_plotter(4)

sex_spike_plotter(16)

sex_spike_plotter(17)


