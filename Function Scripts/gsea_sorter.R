sorting_function = function(df, sorting_parameter = 'avg_log2FC', decreasing = TRUE){
    ranks = df[[sorting_parameter]]
    names(ranks) <- rownames(df)
    ranks_sorted <- sort(ranks, decreasing = TRUE)
return(ranks_sorted)

}

saveRDS(sorting_function,'Functions/gsea_sorter.rds')

