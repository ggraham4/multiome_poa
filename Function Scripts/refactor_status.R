refactor =function(dataframe){
  dataframe$Status =factor(dataframe$Status,
                               levels =c(
                                 'NRM',
                                 'M',
                                 'D',
                                 'E',
                                 'NF',
                                 'F'
                               ))
   dataframe
}

saveRDS(refactor, 'Functions/refactor_status.rds')
