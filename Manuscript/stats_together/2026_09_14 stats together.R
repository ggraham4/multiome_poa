
library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
  clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)
  library(AUCell)

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
colors = c('#1965B0', '#4EB265', '#F7F056', '#7BAFDE', '#DC050C')
pairwise_names <- c(
  "m_i",
  "m_f",
  "i_f"
)

level_map <- c(
  "m" = "M",
  "i" = "I",
  "f" = "F"
)


populate_statistics <- function(...) {

  models <- list(...)

  results <- lapply(
    names(models),
    function(model_name) {

      fit <- models[[model_name]]

      # ------------------------------------------------------
      # Model information
      # ------------------------------------------------------

      model_formula <- formula(fit)

      model_terms <- attr(
        terms(model_formula),
        "term.labels"
      )


      # ------------------------------------------------------
      # Identify grouping variable
      # ------------------------------------------------------

      if("Phase" %in% model_terms) {

        grouping_variable <- "Phase"

      } else if("Status" %in% model_terms) {

        grouping_variable <- "Status"

      } else {

        stop(
          paste0(
            "Model '",
            model_name,
            "' does not contain either Phase or Status."
          )
        )
      }


      # ------------------------------------------------------
      # Identify covariates
      # ------------------------------------------------------

      covariates <- setdiff(
        model_terms,
        grouping_variable
      )


      # ------------------------------------------------------
      # Determine whether this is a binomial GLM
      # ------------------------------------------------------

      is_binomial <- (
        inherits(fit, "glm") &&
        family(fit)$family == "binomial"
      )


      # ======================================================
      # BINOMIAL GLM
      # ======================================================

      if(is_binomial) {

        # ----------------------------------------------------
        # Test grouping variable with likelihood-ratio test
        # ----------------------------------------------------

        null_formula <- update(
          model_formula,
          paste0(". ~ . - ", grouping_variable)
        )

        null_fit <- update(
          fit,
          formula = null_formula
        )

        lr <- anova(
          null_fit,
          fit,
          test = "Chisq"
        )

        p_value <- lr$`Pr(>Chi)`[2]

        deviance <- lr$Deviance[2]

        df <- lr$Df[2]


        # ----------------------------------------------------
        # Create output
        # ----------------------------------------------------

        out <- tibble(

          model = paste0(
            "glm(",
            deparse(model_formula),
            ", family = binomial)"
          ),

          figure = NA_character_,

          anova_p.value = p_value,
          sum_of_squares = deviance,
          Df = df,
          f_value = NA_real_,

          covariate = if(length(covariates) > 0) {
            covariates[1]
          } else {
            NA_character_
          },

          covariate_p.value = NA_real_,
          covariate_sum_of_squares = NA_real_,
          covariate_Df = NA_real_,
          covariate_f_value = NA_real_,

          p.value_adjustment = "none",
          p.value_adjusted = p_value,

          pairwise_test_statistic = "t.ratio"
        )


        # ----------------------------------------------------
        # Initialize ONLY M-I, M-F, I-F
        # ----------------------------------------------------

        for(pair in pairwise_names) {

          out[[paste0(pair, "_p.value")]] <- NA_real_

          out[[paste0(pair, "_statistic")]] <- NA_real_
        }


        # ----------------------------------------------------
        # EMMEANS
        # ----------------------------------------------------

        emm <- emmeans(
          fit,
          specs = as.formula(
            paste0("~ ", grouping_variable)
          )
        )


        pw <- pairs(
          emm,
          adjust = "none"
        ) %>%
          as.data.frame()


        # ----------------------------------------------------
        # Normalize contrast names
        # ----------------------------------------------------

        pw$contrast_norm <- pw$contrast %>%
          str_replace_all("\\s+", "") %>%
          str_replace_all("−", "-")


        # ----------------------------------------------------
        # If Status is used:
        #
        # D corresponds to I.
        #
        # Therefore:
        # M-D = M-I
        # M-F = M-F
        # D-F = I-F
        # ----------------------------------------------------

        if(grouping_variable == "Status") {

          pw$contrast_norm <- pw$contrast_norm %>%
            str_replace_all("^D-", "I-") %>%
            str_replace_all("-D$", "-I")
        }


        # ----------------------------------------------------
        # Fill M-I, M-F, I-F
        # ----------------------------------------------------

        for(pair in pairwise_names) {

          parts <- str_split(
            pair,
            "_",
            simplify = TRUE
          )

          group1 <- level_map[parts[1]]
          group2 <- level_map[parts[2]]

          target <- paste0(
            group1,
            "-",
            group2
          )

          match_row <- pw %>%
            filter(
              contrast_norm == target
            )


          if(nrow(match_row) == 0) {
            next
          }


          out[[paste0(pair, "_p.value")]] <-
            match_row$p.value[1]

          out[[paste0(pair, "_statistic")]] <-
            match_row$t.ratio[1]
        }


        return(out)
      }


      # ======================================================
      # LINEAR MODEL
      # ======================================================

      if(inherits(fit, "lm")) {

        if(grouping_variable != "Phase") {

          stop(
            paste0(
              "Linear model '",
              model_name,
              "' uses '",
              grouping_variable,
              "'. Only Phase is currently supported for lm models."
            )
          )
        }


        # ----------------------------------------------------
        # Partial ANOVA tests
        # ----------------------------------------------------

        drop_tab <- drop1(
          fit,
          test = "F"
        )


        # ----------------------------------------------------
        # Phase statistics
        # ----------------------------------------------------

        phase_row <- which(
          rownames(drop_tab) == "Phase"
        )

        phase_p <-
          drop_tab$`Pr(>F)`[phase_row]

        phase_ss <-
          drop_tab$`Sum of Sq`[phase_row]

        phase_df <-
          drop_tab$Df[phase_row]

        phase_f <-
          drop_tab$`F value`[phase_row]


        # ----------------------------------------------------
        # Covariate statistics
        # ----------------------------------------------------

        if(length(covariates) > 0) {

          covariate <- covariates[1]

          cov_row <- which(
            rownames(drop_tab) == covariate
          )

          if(length(cov_row) > 0) {

            cov_p <-
              drop_tab$`Pr(>F)`[cov_row]

            cov_ss <-
              drop_tab$`Sum of Sq`[cov_row]

            cov_df <-
              drop_tab$Df[cov_row]

            cov_f <-
              drop_tab$`F value`[cov_row]

          } else {

            cov_p <- NA_real_
            cov_ss <- NA_real_
            cov_df <- NA_real_
            cov_f <- NA_real_
          }

        } else {

          covariate <- NA_character_
          cov_p <- NA_real_
          cov_ss <- NA_real_
          cov_df <- NA_real_
          cov_f <- NA_real_
        }


        # ----------------------------------------------------
        # Create output
        # ----------------------------------------------------

        out <- tibble(

          model = paste0(
            "lm(",
            deparse(model_formula),
            ")"
          ),

          figure = NA_character_,

          anova_p.value = phase_p,
          sum_of_squares = phase_ss,
          Df = phase_df,
          f_value = phase_f,

          covariate = covariate,
          covariate_p.value = cov_p,
          covariate_sum_of_squares = cov_ss,
          covariate_Df = cov_df,
          covariate_f_value = cov_f,

          p.value_adjustment = "none",
          p.value_adjusted = phase_p,

          pairwise_test_statistic = "t.ratio"
        )


        # ----------------------------------------------------
        # Initialize ONLY M-I, M-F, I-F
        # ----------------------------------------------------

        for(pair in pairwise_names) {

          out[[paste0(pair, "_p.value")]] <- NA_real_

          out[[paste0(pair, "_statistic")]] <- NA_real_
        }


        # ----------------------------------------------------
        # EMMEANS for Phase
        # ----------------------------------------------------

        emm <- emmeans(
          fit,
          specs = ~ Phase
        )


        pw <- pairs(
          emm,
          adjust = "none"
        ) %>%
          as.data.frame()


        # ----------------------------------------------------
        # Normalize contrasts
        # ----------------------------------------------------

        pw$contrast_norm <- pw$contrast %>%
          str_replace_all("\\s+", "") %>%
          str_replace_all("−", "-")


        # ----------------------------------------------------
        # Fill M-I, M-F, I-F
        # ----------------------------------------------------

        for(pair in pairwise_names) {

          parts <- str_split(
            pair,
            "_",
            simplify = TRUE
          )

          group1 <- level_map[parts[1]]
          group2 <- level_map[parts[2]]

          target <- paste0(
            group1,
            "-",
            group2
          )

          match_row <- pw %>%
            filter(
              contrast_norm == target
            )


          if(nrow(match_row) == 0) {
            next
          }


          out[[paste0(pair, "_p.value")]] <-
            match_row$p.value[1]

          out[[paste0(pair, "_statistic")]] <-
            match_row$t.ratio[1]
        }


        return(out)
      }


      # ------------------------------------------------------
      # Unsupported model
      # ------------------------------------------------------

      stop(
        paste0(
          "Model '",
          model_name,
          "' has unsupported class: ",
          paste(class(fit), collapse = ", ")
        )
      )
    }
  )

  bind_rows(results)
}

# ####fig 3 a, b#####
status_to_phase = c(
  "D"='I',
  'M' = 'M',
  'F' ='F',
  'NF' ='NF',
  'NRM' = 'NRM',
  'E' = 'LI'
)

neurons_only <- subset(obj, 
                     #oligos
                     res0.8_50nn_40PC_45LSI!=2&
                     #microglia
                    res0.8_50nn_40PC_45LSI!=11&
                    #opcs
                    res0.8_50nn_40PC_45LSI!=13&
                    #dividing glia
                    res0.8_50nn_40PC_45LSI!=26&
                    #leuko
                    res0.8_50nn_40PC_45LSI!=20&
                    #ependymal
                    res0.8_50nn_40PC_45LSI!=15
                    &  res0.8_50nn_40PC_45LSI!=1
                    )

total_cells_neuron = neurons_only@meta.data%>%
  group_by(individual, Status)%>%
  summarize(total_cells = n())

n_cells_neuron=neurons_only@meta.data%>%
  group_by(individual, res0.8_50nn_40PC_45LSI)%>%
  summarize(ncells = n())

joint_neuron = total_cells_neuron%>%
  left_join(n_cells_neuron, by = 'individual')

joint_neuron$Status = as.character(joint_neuron$Status)
joint_neuron$Status = factor(joint_neuron$Status, levels = c('M','D','F'))

  sub= subset(joint_neuron, res0.8_50nn_40PC_45LSI ==6 & Status %in% c('M','D','F'))
  mat = cbind(sub$ncells, sub$total_cells-sub$ncells)
  

total_cells = obj@meta.data%>%
  group_by(individual, Status)%>%
  summarize(total_cells = n())

n_cells=obj@meta.data%>%
  group_by(individual, res0.8_50nn_40PC_45LSI)%>%
  summarize(ncells = n())

joint = total_cells%>%
  left_join(n_cells, by = 'individual')

joint$Status = as.character(joint$Status)
joint$Status = factor(joint$Status, levels = c('M','D','F'))

  sub_1= subset(joint, res0.8_50nn_40PC_45LSI ==1 & Status %in% c('M','D','F'))
  mat_1 = cbind(sub_1$ncells, sub_1$total_cells-sub_1$ncells)
  
  
mod_6 = glm(mat ~ Status, data = sub, family = 'binomial')
car::Anova(mod_6, type = 'III')
pairs(emmeans::emmeans(mod_6, 'Status'), adjust = 'none')

mod_1 = glm(mat_1 ~ Status, data = sub_1, family = 'binomial')
car::Anova(mod_1, type = 'III')
pairs(emmeans::emmeans(mod_1, 'Status'), adjust = 'none')
# Fig 4 gh ####

obj = FindSubCluster(obj,
                     1, 'harmony.wsnn', resolution = 0.2, subcluster.name = 'sub_res0.2')
sub_1 = subset(obj, final_clusters ==1)
Idents(sub_1) <- 'sub_res0.2'
sub_1 = subset(sub_1, final_clusters ==1)
sub_1$Status = factor(sub_1$Status, levels = c('NRM','M',"D",'E','NF','F'))

DimPlot(sub_1)

# 11 and 10
cells_ind = sub_1@meta.data%>%
  group_by(individual)%>%
  summarize(n_cells = n())

cells_sub_ind = sub_1@meta.data%>%
  group_by(individual, Status, sub_res0.2)%>%
  summarize(n_cells_in = n())%>%
    subset(Status%in%c('M','D','F'))


cells_total = cells_ind%>%
  right_join(cells_sub_ind, by = 'individual')
cells_total$prop = cells_total$n_cells_in/cells_total$n_cells

cells10 = subset(cells_total, sub_res0.2 == paste0('1_',0) & Status !='NRM')
matrix_10 = cbind(cells10$n_cells_in,cells10$n_cells-cells10$n_cells_in)

mod_10 = glm(matrix_10~Status, family = binomial(), data = cells10)
 car::Anova(mod_10, 3)
pairs(emmeans(mod_10, 'Status'), adjust ='none')%>%as.data.frame()

cells11 = subset(cells_total, sub_res0.2 == paste0('1_',1) & Status !='NRM')
matrix_11 = cbind(cells11$n_cells_in,cells11$n_cells-cells11$n_cells_in)

mod_11 = glm(matrix_11~Status, family = binomial(), data = cells11)
car::Anova(mod_11, 3)
pairs(emmeans(mod_11, 'Status'), adjust ='none')%>%as.data.frame()


#### fig 4 i ####
.subset_if_needed <- function(obj, subcluster) {
if (!is.null(subcluster)) {
obj <- subset(obj, sub_res0.2 == subcluster)
}
obj
}

go_module_aucell <- function(term, obj, subcluster = NULL) {
set.seed(1)
obj <- .subset_if_needed(obj, subcluster)

term2gene <- readRDS("Function Scripts/Dependencies/Term2gene_clown_go2.rds")
term2name <- readRDS('/Users/ggraham/Desktop/multiome_poa/Function Scripts/Dependencies/Term2name.rds')

go_terms <- term2gene %>% left_join(term2name, by = 'go_id')

term_genes        <- go_terms$aocellaris_name[go_terms$go_id == term]
term_name         <- unique(go_terms$go_name[go_terms$go_id == term])
term_genes_in_obj <- unique(term_genes[term_genes %in% rownames(obj)])

print(paste0(length(term_genes_in_obj), ' genes found for: ', term_name))

if (length(term_genes_in_obj) < 1) { return(NULL) }

expr_matrix  <- obj@assays$RNA$data
cell_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE, verbose = FALSE)

gene_set  <- setNames(list(term_genes_in_obj), term)
auc_scores <- AUCell_calcAUC(gene_set, cell_rankings, verbose = FALSE)
scores     <- as.numeric(getAUC(auc_scores)[1, ])

return(list(scores = scores, term_name = term_name))
}

model_go_aucell <- function(term, obj, subcluster = NULL) {
set.seed(1)
obj <- .subset_if_needed(obj, subcluster)

result           <- go_module_aucell(term, obj)
obj$aucell_score <- result$scores

mod_df <- obj@meta.data %>%
group_by(individual, Status) %>%
summarize(mean_mod = mean(aucell_score), .groups = "drop") %>%
mutate(Status = factor(Status, levels = c('NRM', 'M', 'D', 'E', 'NF', 'F'))) %>%
subset(Status %in% c('M','D','F'))

lm(mean_mod ~ Status, data = mod_df)
}

neuron_diff = model_go_aucell('GO:0030182', sub_1, subcluster = '1_0')
neuron_diff%>%anova(test ='Chisq')
pairs(emmeans::emmeans(neuron_diff, 'Status'), adjust = 'none')


#fig 5 de
## cyto ecm 5d ####
library(CytoTRACE)

sub_6 = FindSubCluster(obj, 6, graph.name = "harmony.wsnn")
Idents(sub_6) <- "sub.cluster"
sub_6 = subset(sub_6, final_clusters == 6)

sub_6$Status = factor(
  sub_6$Status,
  levels = c("NRM", "M", "D", "E", "NF", "F")
)

cyto = CytoTRACE(sub_6@assays$RNA$data %>% as.matrix())
sub_6$cyto = cyto$CytoTRACE

degs_plasticity = c(
  "LOC111588913",
  "cntn4",
  "LOC111567620",
  "pcdh10b",
  "sdc2",
  "LOC111585095",
  "bcan",
  "LOC111568896"
)

sub_6$gene_pos = colSums(
  sub_6@assays$RNA$data[degs_plasticity, ]
) > 0

temp = sub_6@meta.data %>%
  filter(Status %in% c('M','D','F'))

model = lmer(
  cyto ~ nCount_RNA + Status * gene_pos + (1 | individual),
  data = temp
)

anova_type3 = car::Anova(model, type = 3)

type3_results = anova_type3 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  select(term, everything())

print(type3_results)

emm = emmeans(model, ~ Status | gene_pos)

pairwise_results = pairs(
  emm,
  adjust = "none"
) %>%
  as.data.frame()

print(pairwise_results)

### fig 5e ####
model_5e = lmer(
  cyto~nCount_RNA+Status+(1|individual),
  data = temp
)
model_5e_type3 = car::Anova(model_5e, type = 3)

emm = emmeans(model_5e, ~Status)

pairwise_results = pairs(
  emm,
  adjust = "none"
) %>%
  as.data.frame()

### 6e- d, s10d-e ####
genes_interest = c(
  "drd3",
  "tacr3a",
  "cckb",
  "pgr",
  "LOC111568069"
)

for (gene in genes_interest) {

  temp = sub_6@meta.data
  temp$gene_pos = sub_6@assays$RNA$data[gene, ] > 0
  temp = temp %>%
    filter(Status %in% c("M","D","F"), gene_pos == TRUE)

  if (nrow(temp) == 0) {
    next
  }

  model = lmer(
    cyto ~ nCount_RNA + Status + (1 | individual),
    data = temp
  )

  anova_table = car::Anova(model, type = 3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term")

  assign(
    paste0("anova_", gene),
    anova_table,
    envir = .GlobalEnv
  )

  pairwise = pairs(
    emmeans(model, "Status"),
    adjust = "none"
  ) %>%
    as.data.frame() %>%
    mutate(gene = gene, .before = 1)

  assign(
    paste0("pairwise_", gene),
    pairwise,
    envir = .GlobalEnv
  )
}

###5fg####

# fig 6 a-c

#fig s6 a,b

#fig s7 b c,

#fig s8 a-d

#fig s9

#fig s10 a-e

