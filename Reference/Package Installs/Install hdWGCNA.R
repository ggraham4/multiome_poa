library(devtools)
BiocManager::install('WGCNA', force = T)
library(WGCNA)

BiocManager::install()
BiocManager::install(c('UCell',
                       'GeneOverlap'))
devtools::install_github('wjawaid/enrichR', force = T)
devtools::install_github('smorabit/hdWGCNA', ref='dev')

library(hdWGCNA)
