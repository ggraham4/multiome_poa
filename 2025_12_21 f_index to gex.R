# spearmann rank correlation to test for similarities to theory curves
library(Seurat)
library(tidyverse)
library(lme4)
library(ggplot2)
library(hdWGCNA)
mecp = readRDS('Functions/mean_expression_cluster_plot.rds')
clown_go<- readRDS('Functions/clown_go2')

scale_2 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

test_formula=readRDS( 'Functions/Theory/f_index_testicular.rds')

f_indexer =function(data.frame){
  vars= c(
  'Log_11KT',
  'mass_final'
)

  newd = data.frame[complete.cases(data.frame[,vars]),]
newd = newd%>%
  subset(!Phase%in% c('S'))

X <- scale(newd[, vars])
pca <- prcomp(X, center = TRUE, scale. = TRUE)


ifelse(pca$rotation['mass_final','PC1']<0, newd$f_index <- pca$x[,1]*-1, newd$f_index <- pca$x[,1])

return(newd)
}

# as an example, lets use gonad data
measures = read.csv('Measures/all_data.csv')
measures$mass_final = measures$mass_final_cm
status_to_phase = list('D'='I',
                       'E' = 'LI',
                       'EP' = 'LIP',
                       'S' = 'IP',
                       'M'='M',
                       'F'='F',
                       'NF'='NF',
                       'NM'='NM')
measures$Phase = unlist(status_to_phase[measures$Status])
measures=f_indexer(measures)

ggplot(measures, aes(x = Phase, y = f_index))+
  geom_boxplot()
#good

measures$testis_curve = sapply( measures$f_index,test_formula)

ggplot(measures, aes(x = f_index, y = testis_curve))+
  geom_point()
#beautiful

#seurat time
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
obj@meta.data$Status = factor(obj@meta.data$Status, levels = c('NRM',
                                                               'M',
                                                               'D',
                                                               'E',
                                                               'NF',
                                                               'F'))



obj@meta.data= obj@meta.data %>%
   left_join(measures, by = join_by('individual'=='Fish'))

#lets try with 6
meta_6 = obj@meta.data%>%
  subset(final_clusters==6)
rna_6 = obj@assays$RNA$data[,obj@meta.data$final_clusters == 6]%>%as.matrix()
rna_6 = t(rna_6)

genes = colnames(rna_6)

individual_to_cell = data.frame(individual = paste0(meta_6$individual),
                                cell = paste0(rownames(rna_6)))%>%
  distinct()

n_cells = nrow(rna_6)
# cycle through genes

ind <- factor(meta_6$individual)

sum_expr <- rowsum(rna_6, ind)          
n_cells  <- table(ind)

mean_expr <- sweep(sum_expr, # arrays to be swept
                   1, #row or column
                   n_cells,  #by
                   "/") # operation

# order f_index by mean_expr
f_df = obj@meta.data%>%
  group_by(individual)%>%
  summarize(testis_curve = mean(testis_curve, na.rm =T))


f_ordered <- f_df[match(rownames(mean_expr), f_df$individual), ]%>%
  subset(!is.na(f_df$testis_curve))

mean_expr_no_na = mean_expr[!is.na(f_df$testis_curve),]

test = cor.test(mean_expr_no_na[,1], 
         f_ordered$testis_curve,
         method = 'spearman')

rownames(mean_expr_no_na) ==f_ordered$individual

cor_df = data.frame()
for(i in 1:length(genes)){
  print(i)
  
  test = cor.test(mean_expr_no_na[,i], 
         f_ordered$testis_curve,
         method = 'spearman')
  
  newd = data.frame(gene = genes[i],
                    S = test$statistic,
                    rho = test$estimate,
                    p.value = test$p.value
                    )
cor_df = rbind(newd, cor_df)
  }

cor_df_successful = cor_df%>%
  subset(!is.na(p.value))

cor_df_successful$q.value = p.adjust(cor_df_successful$p.value, 'fdr', nrow(cor_df_successful))

mecp(obj, 'espn',6)
# thats exactly backward lol
mecp(obj, 'ddx56',6)
# yearh I see it now 

mecp(obj, 'LOC129348661',6)

ggplot(f_df[!is.nan(f_df$testis_curve),], aes(x = mean_expr_no_na[,'LOC129348661'],y= testis_curve))+
  geom_point()

# no dice, lets try a different way

library(mgcv)

### ------------------------------------------------------------
### 1. Build per-individual f_index + testis curve table
### ------------------------------------------------------------

f_df <- obj@meta.data %>%
  group_by(individual) %>%
  summarize(
    f_index      = mean(f_index, na.rm = TRUE),
    testis_curve = mean(testis_curve, na.rm = TRUE),
    .groups = "drop"
  )

# explicitly drop individuals with missing f_index or testis curve
f_df <- f_df %>%
  filter(!is.na(f_index), !is.na(testis_curve))

### ------------------------------------------------------------
### 2. Pseudobulk mean expression (you already did this correctly)
### ------------------------------------------------------------

meta_6 <- obj@meta.data %>%
  subset(final_clusters == 6)

rna_6 <- obj@assays$RNA$data[, meta_6$final_clusters == 6] %>%
  as.matrix() %>%
  t()

genes <- colnames(rna_6)

ind <- factor(meta_6$individual)

sum_expr <- rowsum(rna_6, ind)
n_cells  <- table(ind)

mean_expr <- sweep(sum_expr, 1, n_cells, "/")

# keep only individuals that survived f_df filtering
mean_expr <- mean_expr[rownames(mean_expr) %in% f_df$individual, ]

# reorder f_df to match mean_expr rows
f_df <- f_df[match(rownames(mean_expr), f_df$individual), ]

stopifnot(all(rownames(mean_expr) == f_df$individual))

### ------------------------------------------------------------
### 3. Define common f_index grid + testis trajectory
### ------------------------------------------------------------

f_grid <- seq(
  min(f_df$f_index, na.rm = TRUE),
  max(f_df$f_index, na.rm = TRUE),
  length.out = 100
)

testis_grid <- test_formula(f_grid)

### ------------------------------------------------------------
### 4. Loop over genes: GAM → trajectory correlation
### ------------------------------------------------------------

results <- vector("list", length(genes))

for (i in seq_along(genes)) {

  gene <- genes[i]
  message(i, " / ", length(genes), ": ", gene)

  df_gene <- data.frame(
    expr    = mean_expr[, gene],
    f_index = f_df$f_index
  )

  # strict NA handling
  df_gene <- df_gene[complete.cases(df_gene), ]

  # need enough points to fit a smooth
  if (nrow(df_gene) < 6) next

  # fit GAM
  gam_fit <- try(
    gam(expr ~ s(f_index, k = 5), data = df_gene),
    silent = TRUE
  )
  if (inherits(gam_fit, "try-error")) next

  # check if smooth is significant (real structure vs noise)
  s_tab <- summary(gam_fit)$s.table
  p_smooth <- s_tab[1, "p-value"]
  if (is.null(s_tab) || is.na(p_smooth) || p_smooth > 0.05) next

  # evaluate gene curve on common grid
  gene_grid <- predict(
    gam_fit,
    newdata = data.frame(f_index = f_grid)
  )

  # trajectory correlation
  traj_cor <- suppressWarnings(
    cor.test(gene_grid, testis_grid, method = "spearman")
  )

  # derivative (change-alignment) correlation
  d_gene   <- diff(gene_grid)
  d_testis <- diff(testis_grid)

  deriv_cor <- suppressWarnings(
    cor.test(d_gene, d_testis, method = "spearman")
  )

  results[[i]] <- data.frame(
    gene              = gene,
    rho_traj          = traj_cor$estimate,
    p_traj            = traj_cor$p.value,
    rho_deriv         = deriv_cor$estimate,
    p_deriv           = deriv_cor$p.value,
    edf               = s_tab[1, "edf"]
  )
}

### ------------------------------------------------------------
### 5. Combine + FDR correction
### ------------------------------------------------------------

cor_df <- bind_rows(results)

cor_df_mut <- cor_df %>%
    subset(p_traj !=0 )%>%
  mutate(
    q_traj  = p.adjust(p_traj,  method = "fdr"),
    q_deriv = p.adjust(p_deriv, method = "fdr")
  )%>%
  subset(q_deriv!=0)


### ------------------------------------------------------------
### 6. Example sanity check plot for one gene
### ------------------------------------------------------------

mecp(obj, 'igfbp3',6)

gene_example <- "igfbp3"

df_plot <- data.frame(
  f_index = f_df$f_index,
  expr    = mean_expr[, gene_example]
)

gam_ex <- gam(expr ~ s(f_index, k = 5), data = df_plot)

plot(f_grid, testis_grid, type = "l", lwd = 2,
     xlab = "f_index", ylab = "Scaled value",
     main = gene_example)
lines(
  f_grid,
  scale_2(predict(gam_ex, newdata = data.frame(f_index = f_grid))),
  col = "red", lwd = 2
)
legend("topright", c("Testis curve", "Gene curve"),
       col = c("black", "red"), lwd = 2)


clown_go2(cor_df_mut$gene[cor_df_mut$q_traj<0.05])
