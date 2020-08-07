setwd("tyranni/diversification_analyses/Diversitree_DDD")
getwd()

library(ape)
library(diversitree)
library(phytools)
library(DDD)

# Get mean crown age of families
sample.data <- read.csv("../../Species_name_map_uids_Jetz.csv", header=TRUE)
family.map <- sample.data$howardmoore.family
names(family.map) <- sample.data$tipnamecodes

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')

# Log-likelihood of diversity dependent model
dd.res <- dd_ML(branching.times(tree), initparsopt=c(0.2, 0.01, 1500), idparsopt=1:3)
print(dd.res$loglik)
print(dd.res$lambda)
print(dd.res$mu)
print(dd.res$K)
ddLL <- dd.res$loglik

# Log-likelihood of birth-death model
bd.res <- bd_ML(branching.times(tree), initparsopt=c(0.2, 0.01), idparsopt=1:2)
print(bd.res$loglik)
print(bd.res$lambda0)
print(bd.res$mu0)
bdLL <- bd.res$loglik

# AIC calculation
n <- length(branching.times(tree))
K <- 3
dd.AICc <- -2*ddLL+2*K*(n/(n-K-1))
K <- 2
bd.AICc <- -2*bdLL+2*K*(n/(n-K-1))

print(dd.AICc)
print(bd.AICc)
print(dd.AICc-bd.AICc)
