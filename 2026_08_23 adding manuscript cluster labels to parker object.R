library(Seurat)
library(tidyverse)

# ── Load object ──────────────────────────────────────────────────────────────
parker_path <- "A:/Anemonefish POA Legacy R Objects/parker_2024_realigned_ocellaris.rds"
parker_2024 <- readRDS(parker_path)

# ── clusters20 -> parent_clusters mapping ───────────────────────────────────
# names  = OLD clusters20 IDs
# values = NEW parent cluster labels
# old cluster 13 intentionally excluded (never had a valid new label)
old_to_new_20 <- c('1'='1',
                   '2'='2',
                   '6'='3',
                   '9'='4',
                   '10'='5',
                   '12'='6',
                   '15'='7',
                   '16'='8',
                   '0'='9',
                   '7'='10',
                   '14'='11',
                   '18'='12',
                   '3'='13',
                   '4'='14',
                   '5'='15',
                   '11'='16',
                   '8'='17',
                   '17'='18',
                   '19'='19')

# ── clusters49 -> child_clusters mapping ────────────────────────────────────
# This vector as originally written is new -> old:
#   names  = NEW child cluster labels (e.g. '1a','1b',...,'18')
#   values = OLD clusters49 IDs
new_to_old_49 <- c('1a'='0',
                   '1b'='7',
                   '1c'='11',
                   '1d'='13',
                   '1e'='27',
                   '1f'='29',
                   '1g'='47',
                   '2a'='2',
                   '2b'='6',
                   '2c'='17',
                   '2d'='20',
                   '2e'='23',
                   '2f'='30',
                   '2g'='31',
                   '2h'='41',
                   '2i'='43',
                   '3'='9',
                   '4'='24',
                   '5'='26',
                   '6a'='34',
                   '6b'='38',
                   '7'='37',
                   '8'='39',
                   '9a'= '1',
                   '9b'='4',
                   '9c'='5',
                   '9d'='12',
                   '9e'='16',
                   '9f'='18',
                   '10'='21',
                   '11'='36',
                   '12'='42',
                   '13a'='8',
                   '13b'='10',
                   '13c'='19',
                   '13d'='25',
                   '13e'='45',
                   '13f'='48',
                   '14a'='3',
                   '14b'='15',
                   '14c'='32',
                   '14d'='44',
                   '15a'='14',
                   '15b'='33',
                   '16'='28',
                   '17'='22',
                   '18'='40',
                   '19'='46'
)

# Invert to old -> new: names = OLD clusters49 IDs, values = NEW child labels
# old cluster 35 intentionally excluded (never had a valid new label)
old_to_new_49 <- setNames(names(new_to_old_49), new_to_old_49)

# Sanity check inversion is 1:1
stopifnot(!anyDuplicated(new_to_old_49))

# ── Exclude cells belonging to unmapped old clusters ────────────────────────
# clusters20 == '13' and clusters49 == '35' were never assigned new labels
parker_2024 <- subset(parker_2024, clusters20 != '13')
parker_2024 <- subset(parker_2024, clusters49 != '35')

# ── Apply relabeling ─────────────────────────────────────────────────────────
parker_2024$parent_clusters <- unname(old_to_new_20[as.character(parker_2024$clusters20)])
parker_2024$child_clusters  <- unname(old_to_new_49[as.character(parker_2024$clusters49)])

DimPlot(parker_2024, group.by = 'clusters49', label =T)
DimPlot(parker_2024, group.by = 'child_clusters', label = T)
DimPlot(parker_2024, group.by = 'parent_clusters', label =T)

# ── Overwrite the saved object ──────────────────────────────────────────────
saveRDS(parker_2024, parker_path)
