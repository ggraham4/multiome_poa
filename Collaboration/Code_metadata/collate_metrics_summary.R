#Code to summarize metrics_summary.csv from cellranger count runs

#run on worker node and in directory containing all cellranger count output folders (one per sample):
# $ srun --pty bash
# $ module load R/4.1.2-IGB-gcc-8.2.0
# $ Rscript pathToFile/collate_metrics_summary.R

#output is 



args <- commandArgs(TRUE)

# find all metrics_summary.csv filenames 

filenames <- dir(pattern = "metrics_summary.csv$", recursive = T, full.names = T)

# get short sample names
temp <- strsplit(filenames, "/")
names(filenames) <- sapply(temp, function(x) {
  Lpos <- grep("outs", x)
  x[Lpos-1]
})


# read in the first file

colmetrics <- read.csv(filenames[1])



#read in rest and merge

for (i in 2:length(filenames)) {
  colmetrics <- rbind(colmetrics, read.csv(filenames[i]))
} 

#Add in short sample names

colmetrics <- data.frame(Sample = names(filenames), colmetrics)


write.table(colmetrics, file = "collated_metrics_summaries.txt",
             row.names = FALSE, sep = "\t")
write.csv(colmetrics, file = "collated_metrics_summaries.csv",
             row.names = FALSE)


