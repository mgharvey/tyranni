setwd("tyranni/diversification_analyses/Diversitree_DDD")
getwd()

library(ape)
library(diversitree)
library(phytools)
library(DDD)
library(paleotree)

# Get mean crown age of families
sample.data <- read.csv("../../Species_name_map_uids_Jetz.csv", header=TRUE)
family.map <- sample.data$howardmoore.family
names(family.map) <- sample.data$tipnamecodes

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')

results <- c("Clade", "Gamma", "P-Value", "Trimmed.Gamma", "Trimmed.P-Value")

# Measure gamma for each family

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

# Model comparison and comparison LTT plots
ltt.families <- vector()
ltt.sizes <- vector()
nodes <- nodes[nodes > length(tree$tip.label)]
for(i in 1:length(nodes)) {
	subclade <- extract.clade(tree, nodes[i])
	subclade <- force.ultrametric(subclade)
	if(length(subclade$tip.label) > 20) {
		family <- as.character(family.map[subclade$tip.label[1]])
		print(subclade$tip.label[1])
		gamma <- ltt(subclade, plot=FALSE, gamma=FALSE)
		print(gammatest(gamma))
		
		# Gamme after trimming most recent 2.5 My
		tr.subclade <- timeSliceTree(subclade, sliceTime=2.5)
		tr.gamma <- ltt(tr.subclade, plot=FALSE, gamma=FALSE)
		
		results <- rbind(results, c(family, gammatest(gamma)$gamma, gammatest(gamma)$p, gammatest(tr.gamma)$gamma, gammatest(tr.gamma)$p))

	}
}

write.csv(results, "Gamma_results2.csv")


