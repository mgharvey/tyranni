setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/tb')
getwd()

library(ape)
library(phangorn)
library(phytools)
library(diversitree)
source("../tbsim.R")

results <- vector()

# AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/tb_AOSHM.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$lats.mid)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait <- abs(trait)
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

e.res <- tbsim(subtree, trait, nsim=1000, dx)
results <- rbind(results, c("AOSHM", e.res))


# AOSC

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/tb_AOSC.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSC.txt")

trait <- as.numeric(trait.table$lats.mid)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait <- abs(trait)
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

e.res <- tbsim(subtree, trait, nsim=1000, dx)
results <- rbind(results, c("AOSC", e.res))

# IUCN

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/tb_IUCN.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_IUCN.txt")

trait <- as.numeric(trait.table$lats.mid)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait <- abs(trait)
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

e.res <- tbsim(subtree, trait, nsim=1000, dx)
results <- rbind(results, c("IUCN", e.res))

# HGAPF HM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/tb_HGAPF_AOSHM.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$lats.mid)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait <- abs(trait)
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

e.res <- tbsim(subtree, trait, nsim=1000, dx)
results <- rbind(results, c("HGAPF_HM", e.res))

# Astral HM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/tb_Astral_AOSHM.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$lats.mid)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait <- abs(trait)
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

e.res <- tbsim(subtree, trait, nsim=1000, dx)
results <- rbind(results, c("Astral", e.res))

# Astral missingdropped HM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/tb_Astral_AOSHM_missingdropped.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$lats.mid)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait <- abs(trait)
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

e.res <- tbsim(subtree, trait, nsim=1000, dx)
results <- rbind(results, c("Astral_missingdropped", e.res))

write.table(results, file = "tb_tbsim_lat.txt")
