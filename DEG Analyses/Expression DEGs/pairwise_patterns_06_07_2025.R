### working on pairwise comparison analysis

library(gtools)

possible_relationships = c('<','>','NS')
pairs = c('d_m', 'd_f','f_m')

all_permutations = permutations(n = length(pairs), r = length(possible_relationships), v = possible_relationships, repeats.allowed = T)

colnames(all_permutations) = pairs

all_permutations = as.data.frame(all_permutations)
### here, < means the label on the left is smaller and vice versa
### and NS means no significant difference

#i am going to do this in excel
write.csv(all_permutations, '/Users/ggraham/Desktop/multiome_poa/DEG Analyses/Expression DEGs/pairwise_patterns.csv')
