setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/richness')
getwd()

library(ape)
library(phangorn)
library(phytools)
library(diversitree)

results <- vector()

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("../es/text_files/es_AOSHM.txt")
dx.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt", stringsAsFactors=FALSE)
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_overlap_data.txt", sep=" ", header=TRUE, row.names=1)
trait <- rowSums(trait.table, na.rm=TRUE)

# Need to convert trait names from species names to tipnamecodes
name.map.file <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')
name.map <- name.map.file$tipnamecodes
names(name.map) <- name.map.file$aos.howardmoore.species
name.map.used <- name.map[name.map %in% tree$tip.label]
name.map.used.unique <- name.map.used[!duplicated(name.map.used)]
length(name.map.used.unique)
newnames <- as.character(name.map.used.unique[names(trait)])
names(trait) <- newnames

regions <- dx.table$regions
names(regions) <- as.character(dx.table$V1)
regions <- regions[tree$tip.label]
regions <- regions[!is.na(trait)]
regions[regions == "Neotropics"] <- "NT" # Neotropics
regions[regions == "Nearctic"] <- "NA" # Nearctic
regions[regions == "Africa"] <- "OW" # Nearctic
regions[regions == "AsiaAustralia"] <- "OW" # Nearctic
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(regions)))
regions <- regions[subtree$tip.label]
trait <- trait[subtree$tip.label]

# Colors from map 
library(RColorBrewer)
Pal <- rev(brewer.pal(n=11, c("RdYlBu")))
map.colors <- c("#FEE090", "#A50026", "#E0F3F8")

pdf("boxplot.pdf", width=5.5, height=4)
bp <- boxplot(log(trait)~regions, col=map.colors, ylab="log Species Richness", outline=FALSE)
dev.off()
