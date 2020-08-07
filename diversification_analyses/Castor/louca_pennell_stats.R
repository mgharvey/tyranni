setwd("tyranni/diversification_analyses/Castor")
getwd()

library(ape)
library(phytools)
library(castor)

# Get mean crown age of families
sample.data <- read.csv("../../Species_name_map_uids_Jetz.csv", header=TRUE)
family.map <- sample.data$howardmoore.family
names(family.map) <- sample.data$tipnamecodes

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')

# Fit pulled diversification rate to whole tree (do 100x for final run)
#pdr <- fit_hbd_pdr_on_grid(tree, Ntrials=1)
#pdr

# PDR
#pdr$fitted_PDR

# Present-day speciation rate
#pdr$fitted_rholambda0

# Fit PDR using a grid of times
tree.height <- max(branching.times(tree))
ages <- seq(0, tree.height, length.out=25)
pdr.grid <- fit_hbd_pdr_on_grid(tree, oldest_age=tree.height, age_grid=ages, Ntrials=1)
x.target <- seq(0, tree.height, length.out=1000)
spline_vals <- evaluate_spline(ages, pdr.grid$fitted_PDR, splines_degree=3, Xtarget=x.target)

pdf("PDR_through_time.pdf")
plot(x=x.target, y=spline_vals, type='l', ylab="Pulled Diversification Rate", xlab="Time (My)")
# "A relatively constant pdr over time would be indicative of constant - or only slowly changing - speciation and extinction rates" 
dev.off()