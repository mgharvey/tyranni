setwd("tyranni/diversification_analyses/LTT")
getwd()

library(ape)
library(phytools)
library(diversitree)

# AOSHM pb

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
pdf("AOSHM.pdf")
ltt.plot(tree, log="y", ylim=c(1,7500))
# Fit pure birth model
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
sims <- pbtree(b=pb.fit$par, d=0, t=max(nodeHeights(tree)), nsim=1000)
for(i in 1:length(sims)) {
	ltt.lines(sims[[i]], col="gray")
}
ltt.lines(tree)
dev.off()

# Dups Dropped pb

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_dupsdropped.tre')
tree <- force.ultrametric(tree)
pdf("dupdropped.pdf")
ltt.plot(tree, log="y", ylim=c(1,7500))
# Fit pure birth model
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
sims <- pbtree(b=pb.fit$par, d=0, t=max(nodeHeights(tree)), nsim=1000)
for(i in 1:length(sims)) {
	ltt.lines(sims[[i]], col="gray")
}
ltt.lines(tree)
dev.off()

# AOSC pb

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_Clements.tre')
tree <- force.ultrametric(tree)
pdf("AOSC.pdf")
ltt.plot(tree, log="y", ylim=c(1,7500))
# Fit pure birth model
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
sims <- pbtree(b=pb.fit$par, d=0, t=max(nodeHeights(tree)), nsim=1000)
for(i in 1:length(sims)) {
	ltt.lines(sims[[i]], col="gray")
}
ltt.lines(tree)
dev.off()

# IUCN pb

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_IUCN.tre')
tree <- force.ultrametric(tree)
pdf("IUCN.pdf")
ltt.plot(tree, log="y", ylim=c(1,7500))
# Fit pure birth model
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
sims <- pbtree(b=pb.fit$par, d=0, t=max(nodeHeights(tree)), nsim=1000)
for(i in 1:length(sims)) {
	ltt.lines(sims[[i]], col="gray")
}
ltt.lines(tree)
dev.off()

# 1 My pb

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_1My_collapsed.tre')
tree <- force.ultrametric(tree)
pdf("1My.pdf")
ltt.plot(tree, log="y", ylim=c(1,7500))
# Fit pure birth model
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
sims <- pbtree(b=pb.fit$par, d=0, t=max(nodeHeights(tree)), nsim=1000)
for(i in 1:length(sims)) {
	ltt.lines(sims[[i]], col="gray")
}
ltt.lines(tree)
dev.off()

# 2 My pb

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_2My_collapsed.tre')
tree <- force.ultrametric(tree)
pdf("2My.pdf")
ltt.plot(tree, log="y", ylim=c(1,7500))
# Fit pure birth model
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
sims <- pbtree(b=pb.fit$par, d=0, t=max(nodeHeights(tree)), nsim=1000)
for(i in 1:length(sims)) {
	ltt.lines(sims[[i]], col="gray")
}
ltt.lines(tree)
dev.off()

# HGAPF AOSHM pb

tree <- read.tree('tyranni/species_trees/final_timetrees/HGAPF_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
pdf("HGAPF_AOSHM.pdf")
ltt.plot(tree, log="y", ylim=c(1,7500))
# Fit pure birth model
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
sims <- pbtree(b=pb.fit$par, d=0, t=max(nodeHeights(tree)), nsim=1000)
for(i in 1:length(sims)) {
	ltt.lines(sims[[i]], col="gray")
}
ltt.lines(tree)
dev.off()

# Astral AOSHM pb

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
pdf("Astral_AOSHM.pdf")
ltt.plot(tree, log="y", ylim=c(1,7500))
# Fit pure birth model
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
sims <- pbtree(b=pb.fit$par, d=0, t=max(nodeHeights(tree)), nsim=1000)
for(i in 1:length(sims)) {
	ltt.lines(sims[[i]], col="gray")
}
ltt.lines(tree)
dev.off()

# Astral AOSHM missingdropped pb

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
tree <- force.ultrametric(tree)
pdf("Astral_AOSHM_missingdropped.pdf")
ltt.plot(tree, log="y", ylim=c(1,7500))
# Fit pure birth model
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
sims <- pbtree(b=pb.fit$par, d=0, t=max(nodeHeights(tree)), nsim=1000)
for(i in 1:length(sims)) {
	ltt.lines(sims[[i]], col="gray")
}
ltt.lines(tree)
dev.off()


