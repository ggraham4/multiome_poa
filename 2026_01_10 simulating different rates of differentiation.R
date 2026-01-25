simulator_1_1 = function(stem_cells=10,
                         probability_of_symmetric_division=0.5,
                         symmetric_divisions=5,
                         male_exit = 0.25,
                         intermediate_exit = 0.75, 
                         female_exit = 0.25){


sim_data = data.frame(Sex = rep(c('M','I','F'), 20),
                      stem_cells = stem_cells)

# Fixed division calculator - now uses proper exponential growth
div_calculator = function(n_cells){
  rng = rbinom(1, symmetric_divisions, probability_of_symmetric_division)
  return(n_cells * (2^rng))  # proper symmetric division: cells double
}

# Apply division calculator to each row independently
sim_data$divided_stem_cells = sapply(sim_data$stem_cells, div_calculator)

# Add randomization to differentiation probabilities
sim_data$diff_prob = sapply(sim_data$Sex, function(sex){
  base_prob = ifelse(sim_data$Sex == 'M', male_exit, 
                     ifelse(sim_data$Sex == 'I', intermediate_exit, female_exit))
  
  rbeta(1, base_prob * 10, (1 - base_prob) * 10)
})

sim_data$neurons = sim_data$divided_stem_cells * sim_data$diff_prob
sim_data$stem_cells_remaining = sim_data$divided_stem_cells - sim_data$neurons
sim_data$Sex = factor(sim_data$Sex, levels = c('M', "I",'F'))

p=ggplot(sim_data, aes(x = Sex, y = (divided_stem_cells)))+
    geom_boxplot()+
    geom_jitter()
return(p)

}

simulator_summary_1_1 = function(params, n_reps = 20){
  stem_cells = params[1]
  probability_of_symmetric_division = params[2]
  symmetric_divisions = round(params[3])
  
  sim_data = data.frame(Sex = rep(c('M','D','F'), n_reps),
                        stem_cells = stem_cells)
  
  div_calculator = function(n_cells){
    rng = rbinom(1, symmetric_divisions, probability_of_symmetric_division)
    return(n_cells * (2^rng))
  }
  
  sim_data$divided_stem_cells = sapply(sim_data$stem_cells, div_calculator)
  
  # Calculate summary stats by group
  summary_stats = sim_data %>%
    group_by(Sex) %>%
    summarize(mean_cells = mean(divided_stem_cells),
              sd_cells = sd(divided_stem_cells))
  
  return(summary_stats)
}

objective_function = function(params, n_sims = 10){
  # Run multiple simulations to get stable estimates
  all_means = matrix(NA, nrow = n_sims, ncol = 3)
  all_sds = matrix(NA, nrow = n_sims, ncol = 3)
  
  for(i in 1:n_sims){
    sim_result = simulator_summary_1_1(params)
    all_means[i,] = sim_result$mean_cells
    all_sds[i,] = sim_result$sd_cells
  }
  
  avg_means = colMeans(all_means)
  avg_sds = colMeans(all_sds)
  
  # Calculate sum of squared errors for means
  mean_error = sum((avg_means - target_means)^2)
  
  # Calculate sum of squared errors for SDs (weighted less)
  sd_error = sum((avg_sds - target_sds)^2)
  
  # Combined objective (you can adjust weights)
  total_error = mean_error + 0.5 * sd_error
  
  return(total_error)
}



stem_cells = 10
probability_of_symmetric_division = 0.5 
symmetric_divisions = 5

sim_data = data.frame(Sex = rep(c('M','I','F'), 20),
                      stem_cells = stem_cells)

# Fixed division calculator - now uses proper exponential growth
div_calculator = function(n_cells){
  rng = rbinom(1, symmetric_divisions, probability_of_symmetric_division)
  return(n_cells * (2^rng))  # proper symmetric division: cells double
}

# Apply division calculator to each row independently
sim_data$divided_stem_cells = sapply(sim_data$stem_cells, div_calculator)

# Add randomization to differentiation probabilities
sim_data$diff_prob = sapply(sim_data$Sex, function(sex){
  base_prob = ifelse(sex == 'M', 0.25, 
                     ifelse(sex == 'I', 0.75, 0.25))
  rbeta(1, base_prob * 10, (1 - base_prob) * 10)
})

sim_data$neurons = sim_data$divided_stem_cells * sim_data$diff_prob
sim_data$stem_cells_remaining = sim_data$divided_stem_cells - sim_data$neurons
sim_data$Sex = factor(sim_data$Sex, levels = c('M', "I",'F'))

ggplot(sim_data, aes(x = Sex, y = (stem_cells_remaining)))+
    geom_boxplot()+
    geom_jitter()

#compare it to the real data
library(Seurat)
obj =  readRDS('~/Desktop/optimal_clustering_rna_only.rds')
obj = FindSubCluster(obj, 1, 'harmony.wsnn')
sub_1 = subset(obj, final_clusters==1)
Idents(sub_1) = 'sub.cluster'
DimPlot(sub_1)

n_cells = sub_1@meta.data%>%
  group_by(individual, Status)%>%
  summarize(n_cells = n())
n_sub = sub_1@meta.data%>%
  group_by(individual, sub.cluster)%>%
  summarize(n_sub = n())%>%
  left_join(n_cells, by = 'individual')%>%
  mutate(prop_cells = n_sub/n_cells)

n_sub$Status = factor(n_sub$Status, levels = c('NRM',
                                               'M',
                                               'D',
                                               'E',
                                               'NF',
                                               'F'))
ggplot(subset(n_sub, Status %in% c('M','D','F')), aes(x = Status, y = (prop_cells)))+
    geom_boxplot()+
    geom_jitter()+
  facet_wrap(~sub.cluster, scales = 'free')

 ggplot(subset(n_sub, Status %in% c('M','D','F')), aes(x = Status, y = (n_sub)))+
    geom_boxplot()+
    geom_jitter()+
  facet_wrap(~sub.cluster, scales = 'free')
 
 plot_2 =ggplot(subset(n_sub, Status %in% c('M','D','F') & sub.cluster =='1_1'), aes(x = Status, y = (n_sub)))+
    geom_boxplot()+
    geom_jitter()


# I believe 1_1 is my number of cells after symmetrical divisions, and 1_3 is my number of cells 
# that are differentiating into some other cell type (the remnants)


simulator_1_1(stem_cells=8,
                         probability_of_symmetric_division=.6,
                         symmetric_divisions=3,
                         male_exit = 0.25,
                         intermediate_exit = 0.75, 
                         female_exit = 0.25)

plot_2

data_1_1 = n_sub[n_sub$sub.cluster == '1_1' & n_sub$Status %in% c('M','D','F'),]
data_1_1%>%
  group_by(Status)%>%
  summarize(mean_cells = mean(n_sub),
            sd_cells = sd(n_sub))
#M            30.8     9.04
#D            29.6    15.6 
#F            35.7    17.7 

## optimize
initial_params = c(8, 0.6, 3)
target_means = c(M = 30.8, D = 29.6, F = 35.7)
target_sds = c(M = 9.04, D = 15.6, F = 17.7)

# Use optim with bounds
result = optim(par = initial_params,
               fn = objective_function,
               method = "L-BFGS-B",
               lower = c(1, 0.01, 1),
               upper = c(40, 0.99, 20),
               control = list(trace = 1, maxit = 50))
cat("\nOptimal parameters:\n")
cat("stem_cells:", result$par[1], "\n")
cat("probability_of_symmetric_division:", result$par[2], "\n")
cat("symmetric_divisions:", round(result$par[3]), "\n")
set.seed(1)
plot_2_op = simulator_1_1(stem_cells=8,
                         probability_of_symmetric_division=0.5927467,
                         symmetric_divisions=3,
                         male_exit = 0.25,
                         intermediate_exit = 0.75, 
                         female_exit = 0.25)
library(patchwork)
(plot_2 +ylim(0,100)+plot_2_op+ylim(0,100))
# decently close


### 1_3

one_three =data_1_3%>%
  group_by(individual)%>%
  summarize(cells1_3 = mean(n_sub)
  )

rat = data_1_1%>%
  group_by(individual, Status)%>%
  summarize(cells_1_1 = mean(n_sub))%>%
  left_join(one_three, by = 'individual')%>%
  mutate(ratio_3_1 = cells1_3/cells_1_1)

hist(rat$ratio_3_1)
rat%>%
  group_by(Status)%>%
  summarize(mean = mean(ratio_3_1),
            sd = sd(ratio_3_1))

# this is unexpected 

rat2 = rat%>%
  pivot_longer(cols = c(cells_1_1, cells1_3))

ggplot(rat2, aes(x = name, y = value, color = Status, group =individual ))+
  geom_point()+
  geom_line()
