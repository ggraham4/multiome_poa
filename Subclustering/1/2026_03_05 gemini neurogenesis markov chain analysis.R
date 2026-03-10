# ... [Assuming your 'joined' and 'proportion_matrix_data' steps are the same] ...

# We need to define the lineage order clearly
lineage_order <- c('1_1', '1_3', '1_2', '1_0', '9', 'other_cells')

# Updated Transition Function
# Logic: Probability of moving from A to B is based on the 'share' 
# of the remaining downstream population.
estimate_transitions = function(data_wide) {
  
  # Ensure we only have numeric proportion columns for the math
  mat <- as.matrix(data_wide %>% ungroup() %>% select(all_of(lineage_order)))
  
  # Initialize a dataframe for results
  res <- data.frame(individual = data_wide$individual, Status = data_wide$Status)
  
  # Calculate Forward Probabilities P(Forward)
  # P(A -> B) = Prop(B) / (Prop(A) + Prop(B))
  # This models the 'pull' into the next state.
  res$p.1_1_to_1_3 = mat[,'1_3'] / (mat[,'1_1'] + mat[,'1_3'])
  res$p.1_3_to_1_2 = mat[,'1_2'] / (mat[,'1_3'] + mat[,'1_2'])
  res$p.1_2_to_1_0 = mat[,'1_0'] / (mat[,'1_2'] + mat[,'1_0'])
  res$p.1_0_to_9   = mat[,'9']   / (mat[,'1_0'] + mat[,'9'])
  res$p.9_to_other   = mat[,'other_cells']   / (mat[,'9'] + mat[,'other_cells'])

  # Calculate Stay Probabilities P(Stay)
  # P(Stay) = 1 - P(Forward)
  res$pStay.1_1 = 1 - res$p.1_1_to_1_3
  res$pStay.1_3 = 1 - res$p.1_3_to_1_2
  res$pStay.1_2 = 1 - res$p.1_2_to_1_0
  res$pStay.1_0 = 1 - res$p.1_0_to_9
  res$pStay.9 =1 -res$p.9_to_other
  return(res)
}

# Execute
transition_df <- estimate_transitions(proportion_matrix_data)

# Visualize the 'Kinetics'
transitions_long <- transition_df %>%
  pivot_longer(cols = starts_with("p"), names_to = "Type", values_to = "Prob")

# Filter for just the forward movements to see where the lineage 'speeds up'
forward_plot <- transitions_long %>% filter(!str_detect(Type, "pStay"))

forward_plot$Status = factor(forward_plot$Status, levels = c('M',
                                                                     'D',
                                                                     'E',
                                                                     'NF',
                                                                     'F'))

ggplot(forward_plot, aes(x = Status, y = Prob, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2) +
  theme_minimal() +
  labs(title = "Estimated Transition Probabilities by Group",
       subtitle = "Higher values indicate faster progression to the next stage",
       y = "Probability of Transition")

ggplot(forward_plot, aes(x = Status, y = Prob, fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2) +
  facet_wrap(~Type, scales = "free_y") +
  theme_minimal() +
  labs(title = "Estimated Transition Probabilities by Group",
       subtitle = "Higher values indicate faster progression to the next stage",
       y = "Probability of Transition")




model_outs = data.frame()
for(i in unique(forward_plot$Type)){
  model = lm(Prob~Status, data = subset(forward_plot, Type ==i ))
  chisq = anova(model, test = 'Chisq')
  av_p = chisq$`Pr(>F)`[1]
  newd = data.frame(transition= i,
                    p = av_p)
  model_outs = rbind(model_outs, newd)

}

# no significant differences, but even still I think worth reporting, or maybe 
# this is good reason to do it at the trajectory level
