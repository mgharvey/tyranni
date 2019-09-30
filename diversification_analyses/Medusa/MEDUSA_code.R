setwd('/Users/derryberry/Documents/medusa')
getwd()

library(ape)
library(MEDUSA)

tree1 <- read.tree('/Users/derryberry/Documents/tyranni_trees/timetrees/T400F_AOS_HowardMoore.tre')
res1 <- MEDUSA(tree1)
summ1 <- medusaSummary(res1, show.tip.label=FALSE, main="AOS")
save(tree1, res1, summ1, file="/Users/derryberry/Documents/medusa/medusa_T400F_AOSHM_out.RData")
