#Install Cell Chat
library(devtools)
library(BiocManager)
#install dependencies: ComplexHeatmap, BiocGenerics
devtools::install_github("jokergoo/ComplexHeatmap")
BiocManager::install("BiocGenerics")

# next try installing NeuronChat
devtools::install_github("Wei-BioMath/NeuronChat")

library(NeuronChat)
