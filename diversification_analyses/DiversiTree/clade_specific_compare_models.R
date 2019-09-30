setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/ltt_plots/clade_specific_plots')
getwd()

library(ape)
library(phytools)

# Get mean crown age of families
sample.data <- read.csv("/Users/michaelharvey/Documents/research/Tyranni/Species_name_map_uids.csv", header=TRUE)
family.map <- sample.data$howardmoore.family
names(family.map) <- sample.data$tipnamecodes

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')

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
pdf("ltt_clades.pdf", width=5, height=5)
ltt.families <- vector()
ltt.sizes <- vector()
nodes <- nodes[nodes > length(tree$tip.label)]
for(i in 1:length(nodes)) {
	subclade <- extract.clade(tree, nodes[i])
	ltt.families <- c(ltt.families, as.character(family.map[subclade$tip.label[1]]))
	ltt.sizes <- c(ltt.sizes, length(subclade$tip.label))
	if(i == 1) {
		ltt.plot(subclade, ylim=c(1,500), xlim=c(-x,0), log="y", col="gray")		
	} else {
		ltt.lines(subclade, col="gray")				
	}
	par(new=TRUE)
}
dev.off()
names(ltt.sizes) <- ltt.families
sort(ltt.sizes, decreasing=TRUE)
