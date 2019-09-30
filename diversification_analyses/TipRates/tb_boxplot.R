setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/tb')
getwd()

library(ape)
library(phangorn)
library(phytools)
library(diversitree)

results <- vector()

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/tb_AOSHM.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt", stringsAsFactors=FALSE)

trait <- trait.table$regions
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == "Neotropics"] <- "NT" # Neotropics
trait[trait == "Nearctic"] <- "NA" # Nearctic
trait[trait == "Africa"] <- "OW" # Nearctic
trait[trait == "AsiaAustralia"] <- "OW" # Nearctic
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

# Colors from map 
library(RColorBrewer)
Pal <- rev(brewer.pal(n=11, c("RdYlBu")))
map.colors <- c("#74ADD1", "#E0F3F8", "#FEE090")

pdf("boxplot.pdf", width=5.5, height=4)
bp <- boxplot(dx~trait, col=map.colors, ylab="log Species Age", outline=FALSE)
dev.off()
