setwd("tyranni/diversification_analyses/PGLS")
getwd()

library(ape)
library(caper)
library(phytools)
library(MEDUSA)
library(BAMMtools)
library(scales)
library(phylotate)

# Get ES
es.table <- read.table("../SummaryStats/es/text_files/es_AOSHM.txt")
es <- es.table$x
names(es) <- rownames(es.table)

# Get tree/diversity data

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
#overlap.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_overlap_data_allbirds.txt", sep=" ", header=TRUE, row.names=1)
#overlap.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_overlap_data_passerines.txt", sep=" ", header=TRUE, row.names=1)
overlap.table <- read.table("tyranni/other_data/Range_overlap_data.txt", sep=" ", header=TRUE, row.names=1)
overlap <- rowSums(overlap.table, na.rm=TRUE)

# Need to convert overlap names from species names to tipnamecodes
name.map.file <- read.csv('../../Species_name_map_uids.csv')
name.map <- name.map.file$tipnamecodes
names(name.map) <- name.map.file$aos.howardmoore.species
name.map.used <- name.map[name.map %in% tree$tip.label]
name.map.used.unique <- name.map.used[!duplicated(name.map.used)]
length(name.map.used.unique)
newnames <- as.character(name.map.used.unique[names(overlap)])
names(overlap) <- newnames

# Diversity vs. speciation rate
region.data <- read.table("tyranni/other_data/Range_data_AOSHM_Olson_broad.txt")
areas <- cbind(region.data[,3:11])
rownames(areas) <- as.character(region.data$V1)
region.names <- c("WI", "AM", "AN", "NA", "OW", "PA", "DT", "AF", "CA")
colnames(areas) <- region.names
NewWorld.only <- areas[!colnames(areas)[apply(areas,1,which.max)] == "OW",]
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(NewWorld.only))) # Drop Old World
overlap <- overlap[names(overlap) %in% tree$tip.label]

overlap <- overlap[tree$tip.label]
overlap <- overlap[!is.na(overlap)]
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(overlap)))
overlap <- overlap[subtree$tip.label]

# Range size
range.data <- read.table("tyranni/other_data/Range_data_AOSHM.txt", header=1)
areas <- range.data$areas
names(areas) <- range.data$V1

# ES vs richness + range area

dframe <- data.frame(subtree$tip.label, es[subtree$tip.label], overlap, areas[subtree$tip.label])
colnames(dframe) <- c("Species", "Speciation", "Species.Richness", "Range.Area")
data <- comparative.data(data=dframe, phy=subtree, names.col="Species")
full <- pgls(Speciation ~ Species.Richness + Range.Area, data=data)
sum <- summary(full)
sum
full.log <- pgls(Speciation ~ Species.Richness + log(Range.Area), data=data)
sum <- summary(full.log)
sum


# Model comparison (leave-one-out approach)
drop.area <- pgls(Speciation ~ Species.Richness, data=data)
drop.richness <- pgls(Speciation ~ Range.Area, data=data)
drop.area$aicc-full$aicc
drop.richness$aicc-full$aicc

# Akaike weights
full.wgt <- exp(-0.5 * full$aicc)
drop.area.wgt <- exp(-0.5 * drop.area$aicc)
drop.richness.wgt <- exp(-0.5 * drop.richness$aicc)
full.wgt/(full.wgt+drop.area.wgt+drop.richness.wgt)
drop.area.wgt/(full.wgt+drop.area.wgt+drop.richness.wgt)
drop.richness.wgt/(full.wgt+drop.area.wgt+drop.richness.wgt)

# 
richness.area <- pgls(Species.Richness ~ Range.Area, data=data)
summary(richness.area)
