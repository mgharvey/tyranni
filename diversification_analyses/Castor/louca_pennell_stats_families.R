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
		print(family)
		
		# Fit PDR using a grid of times
		tree.height <- max(branching.times(subclade))
		ages <- seq(0, tree.height, length.out=15)
		pdr.grid <- fit_hbd_pdr_on_grid(subclade, oldest_age=tree.height, age_grid=ages, Ntrials=10, Nbootstraps=100, min_PDR=-100, max_PDR=+100, Nthreads=2, max_model_runtime=1)
		save(pdr.grid, file=paste0("./pdr_fits_families/pdr_fit_", family, ".Rdata"))
		x.target <- seq(0, tree.height, length.out=1000)
		spline_vals <- evaluate_spline(ages, pdr.grid$fitted_PDR, splines_degree=3, Xtarget=x.target)
		spline_vals.lower <- evaluate_spline(ages, pdr.grid$CI50lower$PDR, splines_degree=3, Xtarget=x.target)
		spline_vals.upper <- evaluate_spline(ages, pdr.grid$CI50upper$PDR, splines_degree=3, Xtarget=x.target)

		pdf(paste0("./pdr_fits_families/PDR_through_time_", family, ".pdf"))
		plot(1, type="n", xlim=c(0,tree.height), ylim=c(-20,20), ylab="Pulled Diversification Rate", xlab="Time (My)")
		polygon(c(x.target, rev(x.target)), c(spline_vals.upper, rev(spline_vals.lower)), col="lightgray", border=NA)
		lines(x=x.target, y=spline_vals, xlim=c(0,tree.height))
		# "A relatively constant pdr over time would be indicative of constant - or only slowly changing - speciation and extinction rates" 
		dev.off()

	}
}


