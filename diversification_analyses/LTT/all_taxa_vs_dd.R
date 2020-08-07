setwd("tyranni/diversification_analyses/LTT")
getwd()

library(ape)
library(phytools)
library(diversitree)
library(DDD)

# AOSHM bd

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
# Fit birth-death model
dd.fit <- dd_ML(branching.times(tree), initparsopt=c(0.2, 0.01, 1500), idparsopt=1:3)

pdf("AOSHM_dd.pdf")
ltt.plot(tree, log="y", ylim=c(1,7500))
for(i in 1:1000) {
	sim <- dd_sim(unlist(dd.fit[1:3]), age=max(nodeHeights(tree)))
	ltt.lines(sim[[1]], col="gray")
}
ltt.lines(tree)
dev.off()
