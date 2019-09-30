setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/sr')
getwd()

library(ape)
library(phangorn)
library(phytools)
library(diversitree)
source("../FISSE.R")

results <- vector()

# AOSHM

dx.table <- read.table("./text_files/sr_AOSHM.txt")
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(dx.table)))
tree <- force.ultrametric(tree)
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$regions)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 0
trait[trait == 2] <- 1
trait[trait == 1] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$sr)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("AOSHM", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

# AOSHC

dx.table <- read.table("./text_files/sr_AOSC.txt")
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(dx.table)))
tree <- force.ultrametric(tree)
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSC.txt")

trait <- as.numeric(trait.table$regions)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 0
trait[trait == 2] <- 1
trait[trait == 1] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$sr)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("AOSC", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

# IUCN

dx.table <- read.table("./text_files/sr_IUCN.txt")
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(dx.table)))
tree <- force.ultrametric(tree)
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_IUCN.txt")

trait <- as.numeric(trait.table$regions)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 0
trait[trait == 2] <- 1
trait[trait == 1] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$sr)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("IUCN", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

# HGAPF HM

dx.table <- read.table("./text_files/sr_HGAPF_AOSHM.txt")
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(dx.table)))
tree <- force.ultrametric(tree)
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$regions)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 0
trait[trait == 2] <- 1
trait[trait == 1] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$sr)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("HGAPF_HM", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

# Astral HM

dx.table <- read.table("./text_files/sr_Astral_AOSHM.txt")
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(dx.table)))
tree <- force.ultrametric(tree)
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$regions)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 0
trait[trait == 2] <- 1
trait[trait == 1] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$sr)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("Astral", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

# Astral missingdropped HM

dx.table <- read.table("./text_files/sr_Astral_AOSHM_missingdropped.txt")
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(dx.table)))
tree <- force.ultrametric(tree)
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$regions)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 0
trait[trait == 2] <- 1
trait[trait == 1] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$sr)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("Astral_missingdropped", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

colnames(results) <- c("test", "lambda0", "lambda1", "null_mean_diff", "pval")
write.table(results, file = "sr_FISSE_OldWorld_NewWorld.txt")
