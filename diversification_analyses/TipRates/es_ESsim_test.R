setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/es')
getwd()

library(ape)
library(phangorn)
library(phytools)
library(diversitree)
source("../essim.R")

results <- vector()

# AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/es_AOSHM.txt")
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

trait <- as.numeric(trait.table$lats.mid)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait <- abs(trait)
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]
dx <- as.numeric(dx.table$x)
names(dx) <- rownames(dx.table)
dx <- dx[subtree$tip.label]

e.res <- essim(subtree, trait, nsim=1000, dx)
e.res

phy <- subtree
nsim = 1000
	
es <- dx

es <- es[phy$tip.label] # log transform
trait <- trait[phy$tip.label]
	
# Pearson's correlation between log inverse equal splits statistic and trait
res <- cor.test(es, trait, method="pearson")
plot(es~trait)
# Fit Brownian motion model to get diffusion rate and root state estimates
vv <- vcv.phylo(as.phylo(phy))
onev <- matrix(rep(1, length(trait)), nrow=length(trait), ncol=1)
root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% trait))
rate <- as.vector((t(trait-root) %*% solve(vv) %*% (trait-root))/length(trait))
root
rate

# Brownian simulations 
sims <- t(rmvnorm(nsim, sigma=rate*vv))
rownames(sims) <- rownames(vv)
		
# Pearson's correlations of simulated datasets
sim.r <- sapply(1:nsim, function(x) cor.test(es[as.vector(rownames(sims))], sims[,x], method="pearson")$estimate)
	
# Calculate the two-tailed p value
corr <- res$estimate
upper <- (length(sim.r[sim.r >= corr])+1)/(nsim+1)
lower <- (length(sim.r[sim.r <= corr])+1)/(nsim+1)
pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed

if(missing(return.es)) { # output just rho and p value
	result <- as.vector(c(corr, pval))
	names(result) <- c("rho", "P Value")
	return(result)
} else { # output rho, p value, and list of es values
	result <- as.vector(c(corr, pval, list(es)))
	names(result) <- c("rho", "P Value", "es")
	return(result)		
}
