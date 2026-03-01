### adding cellphone data ato the github


direct = dir('A:/cellphonedb_05_12_2025', full.names=T,recursive = FALSE
)
folders = list.dirs(direct ,recursive = T)

for(i in folders){
  print(i)
  
  if(i %in% c('A:/cellphonedb_05_12_2025/v5.0.0',
              "A:/cellphonedb_05_12_2025/v5.0.0/sources",
              "A:/cellphonedb_05_12_2025/GH" )){
    next
  }
  else{
    subdir = dir(i)
    analysis_means_path = list.files(path = i, pattern = 'analysis_means')
    
    analysis_means = read_tsv(paste0(i,'/',analysis_means_path ))
    
    write_tsv(analysis_means, paste0("Python/CellPhoneDb/Analysis Means/",analysis_means_path))
    
  }
}