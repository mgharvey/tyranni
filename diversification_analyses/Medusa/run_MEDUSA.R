setwd("tyranni/diversification_analyses/Medusa")
getwd()

library(ape)
library(MEDUSA)

tree1 <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')
res1 <- MEDUSA(tree1)
summ1 <- medusaSummary(res1, show.tip.label=FALSE, main="AOS")
save(tree1, res1, summ1, file="medusa_T400F_AOSHM_out.RData")
