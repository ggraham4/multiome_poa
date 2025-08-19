### final defined DEGs fr fr this time
setwd("C:/Users/Gabe/Desktop/multiome_poa")
define_degs2 = readRDS('Functions/define_degs2')
path_prefix = "DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/"
directory = dir("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering")
coalesced_data = data.frame()

for(i in directory){
  print(i)
  data = read.csv(paste0(path_prefix, i))
    data = subset(data, av_q.value<0.05 & singular == F)
  data2 = define_degs2(data)
  coalesced_data= rbind(coalesced_data, data2)
  
}

write.csv(coalesced_data, 'DEG Outputs/FINAL_degs_classified_08_11_2025.csv')
