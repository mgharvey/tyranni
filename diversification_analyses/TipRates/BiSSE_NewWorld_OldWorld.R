setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/es')
getwd()

library(ape)
library(phangorn)
library(phytools)
library(diversitree)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
trait.table <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM.txt")

# BiSSE Old World vs. New World

trait <- as.numeric(trait.table$regions)
names(trait) <- as.character(trait.table$V1)
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
trait[trait == 4] <- 0 # Neotropics
trait[trait == 3] <- 0 # Nearctic
trait[trait == 2] <- 1 # Australasia
trait[trait == 1] <- 1 # Africa
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]

lik <- make.bisse(subtree, trait)
p <- starting.point.bisse(subtree)
fit <- find.mle(lik, p)

lik.l <- constrain(lik, lambda1 ~ lambda0)
fit.l <- find.mle(lik.l, p[argnames(lik.l)])
res <- anova(fit, equal.l=fit.l)

print(res)
print(coef(fit))