setwd('~/Desktop/QuaSSE/')
getwd()

library(ape)
library(phytools)
library(diversitree)

tree <- read.tree('~/Desktop/QuaSSE/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
trait.table <- read.table("~/Desktop/QuaSSE/Range_overlap_data.txt", sep=" ", header=TRUE, row.names=1)
trait <- rowSums(trait.table, na.rm=TRUE)

# Need to convert trait names from species names to tipnamecodes
name.map.file <- read.csv('~/Desktop/QuaSSE/Species_name_map_uids.csv')
name.map <- name.map.file$tipnamecodes
names(name.map) <- name.map.file$aos.howardmoore.species
name.map.used <- name.map[name.map %in% tree$tip.label]
name.map.used.unique <- name.map.used[!duplicated(name.map.used)]
newnames <- as.character(name.map.used.unique[names(trait)])
names(trait) <- newnames
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
trait <- trait[subtree$tip.label]

# QuaSSE

# constant
p.constant <- starting.point.quasse(subtree, trait)
xr <- range(trait) + c(-1,1) * 20 * p.constant["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])
nodrift <- function(f) constrain(f, drift ~ 0)
lik.e <- make.quasse(subtree, trait, 1, constant.x, linear.x)
p.e <- c(p.constant[1:2], m.m=0, p.constant[3])
fit.e <- find.mle(nodrift(lik.e), p.e, verbose=0)			
print("extinction variable complete")
print(coef(fit.e))
print(logLik(fit.e))
