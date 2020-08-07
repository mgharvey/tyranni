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
#ages <- seq(0, tree.height, length.out=25)
ages <- seq(0, 35, length.out=25)
#pdr.grid <- fit_hbd_pdr_on_grid(tree, oldest_age=tree.height, age_grid=ages, Ntrials=1)
pdr.grid <- fit_hbd_pdr_on_grid(tree, oldest_age=35, age_grid=ages, Ntrials=10, Nbootstraps=100, min_PDR=-100, max_PDR=+100, Nthreads=2, max_model_runtime=1)
save(pdr.grid, file="pdr_fit_truncated.Rdata")
x.target <- seq(0, 35, length.out=1000)
spline_vals <- evaluate_spline(ages, pdr.grid$fitted_PDR, splines_degree=3, Xtarget=x.target)
spline_vals.lower <- evaluate_spline(ages, pdr.grid$CI50lower$PDR, splines_degree=3, Xtarget=x.target)
spline_vals.upper <- evaluate_spline(ages, pdr.grid$CI50upper$PDR, splines_degree=3, Xtarget=x.target)

pdf("PDR_through_time_truncated_SE.pdf")
plot(1, type="n", xlim=c(0,30), ylim=c(-20,20), ylab="Pulled Diversification Rate", xlab="Time (My)")
polygon(c(x.target, rev(x.target)), c(spline_vals.upper, rev(spline_vals.lower)), col="lightgray", border=NA)
lines(x=x.target, y=spline_vals, xlim=c(0,30))
# "A relatively constant pdr over time would be indicative of constant - or only slowly changing - speciation and extinction rates" 
dev.off()


#load("pdr_fit_truncated.Rdata")
