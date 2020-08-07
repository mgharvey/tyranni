setwd("tyranni/diversification_analyses/ClaDS")
getwd()

options(scipen=999)

library(ape)
library(RPANDA)

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')

### Or use simulated test data ###
#obj <- sim_ClaDS(lambda_0=0.1, mu_0=0.5, sigma_lamb=0.7, alpha_lamb=0.90, condition="taxa", taxa_stop = 20, prune_extinct = TRUE)
#tree <- obj$tree

#sampler <- fit_ClaDS0(tree=tree, name="ClaDS0_1.Rdata", nCPU=1, pamhLocalName="local1_", iteration=1000000, thin=20000, update=1000, adaptation=5)

load("ClaDS0_1.Rdata")

plot_ClaDS0_chains(Cl0_chains, burn=0, param=1:4)

MAPS <- getMAPS_ClaDS0(tree, Cl0_chains, thin=10)
plot_ClaDS_phylo(tree, MAPS[-(1:3)])

# branch-specific speciation rates in the order of phylo$edges
edge.rates <- MAPS[-(1:3)] 
# get the node numbers of the tips
nodes <- sapply(tree$tip.label, function(x,y) which(y==x), y=tree$tip.label)
# get the edge rates for those nodes
tip.rates <- setNames(edge.rates[sapply(nodes, function(x,y) which(y==x),y=tree$edge[,2])], names(nodes))
write.csv(tip.rates, "ClaDS0_1_tiprates.csv")
plot(tip.rates~es[names(tip.rates)])



Cl0_chains
1000000/1000
length(Cl0_chains[[1]][,1])
