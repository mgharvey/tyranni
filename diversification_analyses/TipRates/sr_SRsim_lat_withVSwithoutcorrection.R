setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/sr')
getwd()

library(ape)
library(phangorn)
library(phytools)
library(diversitree)
source("../tbsim.R")

results <- vector()

# without node correction

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/sr_AOSHM.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$lats.mid)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait <- abs(trait)
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$sr)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

e.res <- tbsim(subtree, trait, nsim=1000, dx)
results <- rbind(results, c("AOSHM", e.res))

# with node correction

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/sr_AOSHM.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$lats.mid)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait <- abs(trait)
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$sr_nd_resid)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

e.res <- tbsim(subtree, trait, nsim=1000, dx)
results <- rbind(results, c("AOSHM", e.res))


write.table(results, file = "sr_tbsim_lat_withVSwithoutcorrection.txt")
