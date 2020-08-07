setwd("tyranni/diversification_analyses/BAMM")
getwd()

require(BAMMtools)
require(caper)
require(scales)

##################################################

# Get all the data:

##################################################

# Get the tree
tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- ladderize(tree)
#plot(tree, show.tip.label=FALSE)
#write.tree(tree, "ladderized_HGAP_spp.tre")
#write.tree(tree, "ladderized_HGAP_1My.tre")

# Get event data from a BAMM run on the tree above
#ed <- getEventData(tree, './bamm_100prior_AOSHM/100prior_AOSHM_event_data.txt', burnin=0.1, nsamples=1250) 
#save(ed, file="./bamm_100prior_AOSHM/100prior_AOSHM_eventsample.rda")
load("./bamm_100prior_AOSHM/100prior_AOSHM_eventsample.rda")

##################################################

# BAMM Plots:

##################################################

#plot(ed)

# Credible shift set
#css <- credibleShiftSet(ed, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
#css$number.distinct
#summary(css)
#plot.credibleshiftset(css)

# Best shift set
best <- getBestShiftConfiguration(ed, expectedNumberOfShifts=100)
plot.ed <- plot(ed, pal="RdYlBu", logcolor=TRUE)
addBAMMshifts(best, cex=1, bg="black")
addBAMMlegend(plot.ed)

# Identify number of tips subtending those shifts

subtree.sizes <- vector()
shift.nodes <- best$eventData[[1]]$node
for(i in 1:length(shift.nodes)) {
	subtree <- extract.clade(tree, shift.nodes[i])
	subtree.size <- length(subtree$tip.label)
	subtree.sizes <- c(subtree.sizes, subtree.size)
}
sum(subtree.sizes[2:length(subtree.sizes)])
sum(subtree.sizes[2:length(subtree.sizes)])/length(ed$tip.label)

st <- max(branching.times(tree))
plotRateThroughTime(ed, intervalCol="black", avgCol="black", start.time=st, ylim=c(0,1), cex.axis=1)
