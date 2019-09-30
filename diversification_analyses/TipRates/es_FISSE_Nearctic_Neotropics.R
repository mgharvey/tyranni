setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/es')
getwd()

library(ape)
library(phangorn)
library(phytools)
library(diversitree)
source("../FISSE.R")

results <- vector()

# AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/es_AOSHM.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

sub.trait.table <- trait.table[trait.table$region %in% c("Nearctic", "Neotropics"),]
trait <- as.numeric(sub.trait.table$regions)
names(trait) <- as.character(sub.trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("AOSHM", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

# AOSHC

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/es_AOSC.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSC.txt")

sub.trait.table <- trait.table[trait.table$region %in% c("Nearctic", "Neotropics"),]
trait <- as.numeric(sub.trait.table$regions)
names(trait) <- as.character(sub.trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("AOSC", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

# IUCN

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/es_IUCN.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_IUCN.txt")

sub.trait.table <- trait.table[trait.table$region %in% c("Nearctic", "Neotropics"),]
trait <- as.numeric(sub.trait.table$regions)
names(trait) <- as.character(sub.trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("IUCN", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

# HGAPF HM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/es_HGAPF_AOSHM.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

sub.trait.table <- trait.table[trait.table$region %in% c("Nearctic", "Neotropics"),]
trait <- as.numeric(sub.trait.table$regions)
names(trait) <- as.character(sub.trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("HGAPF_HM", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

# Astral HM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/es_Astral_AOSHM.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

sub.trait.table <- trait.table[trait.table$region %in% c("Nearctic", "Neotropics"),]
trait <- as.numeric(sub.trait.table$regions)
names(trait) <- as.character(sub.trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("Astral", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

# Astral missingdropped HM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/es_Astral_AOSHM_missingdropped.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

sub.trait.table <- trait.table[trait.table$region %in% c("Nearctic", "Neotropics"),]
trait <- as.numeric(sub.trait.table$regions)
names(trait) <- as.character(sub.trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0
trait[trait == 3] <- 1
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
f.res.NW
results <- rbind(results,  c("Astral_missingdropped", f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff, f.res.NW$pval))

colnames(results) <- c("test", "lambda0", "lambda1", "null_mean_diff", "pval")
write.table(results, file = "es_FISSE_Nearctic_Neotropics.txt")
