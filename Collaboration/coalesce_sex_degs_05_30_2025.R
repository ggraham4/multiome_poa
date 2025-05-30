directory = dir("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering")

define_degs<- readRDS('Functions/define_degs')

full_data <- data.frame()
for(i in directory){
  path <- "DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering"
  data <- read.csv(paste0(path,'/', i))
  data <- define_degs(data)
  full_data <- rbind(data, full_data)
}

subset_data <- subset(full_data, !is.na(class))

write.csv(subset_data, "Collaboration/sex_degs_05_30_2025.csv")
