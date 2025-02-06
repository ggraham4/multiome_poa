#Install Cell Chat
library(devtools)
library(BiocManager)
#install dependencies: ComplexHeatmap, BiocGenerics
devtools::install_github("jokergoo/ComplexHeatmap")
BiocManager::install("BiocGenerics")

# next try installing CellChat
devtools::install_github("sqjin/CellChat")

library(CellChat)
