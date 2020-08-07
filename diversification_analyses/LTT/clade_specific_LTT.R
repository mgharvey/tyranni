setwd("tyranni/diversification_analyses/LTT")
getwd()

library(ape)
library(diversitree)
library(phytools)
library(DDD)
library(RColorBrewer)

# Get mean crown age of families
sample.data <- read.csv("../../Species_name_map_uids.csv", header=TRUE)
family.map <- sample.data$howardmoore.family
names(family.map) <- sample.data$tipnamecodes

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')

# Which families in tree?
families.present <- unique(family.map[tree$tip.label])

crown.ages <- vector()
for(i in 1:length(families.present)) {
	family.tips <- tree$tip.label[which(family.map[tree$tip.label] == families.present[i])]
	if(length(family.tips) > 1){
		subclade <- extract.clade(tree, findMRCA(tree, tip=family.tips))
		crown.ages <- c(crown.ages, max(nodeHeights(subclade)))
		#plot(extract.clade(tree, findMRCA(tree, tip=family.tips)))
	} else {
		crown.ages <- c(crown.ages, NA)
	}
}
names(crown.ages) <- families.present
crown.ages
median(crown.ages, na.rm=TRUE)

# Trim using time threshold
x <- median(crown.ages, na.rm=TRUE) # Mean age of a family
h <- nodeHeights(tree)
t <- max(h)-x # time from the root
h1 <- which(h[,1] < t) # identify all edges crossing time x
h2 <- which(h[,2] > t)
ii <- intersect(h1, h2)
nodes <- tree$edge[ii,2] # all daughter nodes of those edges

# LTT plot for each subclade
nodes <- nodes[nodes > length(tree$tip.label)]

# Get list of lineages with >20 species
ltt.nodes <- vector()
n.cols <- 0
for(i in 1:length(nodes)) {
	subclade <- extract.clade(tree, nodes[i])
	if(length(subclade$tip.label) > 20) {
		n.cols <- n.cols+1
		ltt.nodes <- c(ltt.nodes, nodes[i])
	}
}
ltt.nodes

# Make LTT plots
ltt.families <- vector()
ltt.cols <- colorRampPalette(brewer.pal(8, "Set1"))(n.cols) 
for(i in 1:length(ltt.nodes)) {
	subclade <- extract.clade(tree, ltt.nodes[i])
	ltt.families <- c(ltt.families, as.character(family.map[subclade$tip.label[1]]))
	if(i == 1) {
		ltt.plot(subclade, ylim=c(1,500), xlim=c(-x,0), log="y", col=ltt.cols[i])		
	} else {
		ltt.lines(subclade, col=ltt.cols[i])				
	}
	par(new=TRUE)
}
legend("topleft", ltt.families, lty=1, col=ltt.cols, bty="n")



