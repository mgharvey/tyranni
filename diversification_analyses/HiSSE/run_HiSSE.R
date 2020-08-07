setwd("tyranni/diversification_analyses/HiSSE")
getwd()

library(ape)
library(phytools)
library(diversitree)
library(hisse)

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
trait.table <- read.table("tyranni/other_data/Range_overlap_data.txt", sep=" ", header=TRUE, row.names=1)
trait <- rowSums(trait.table, na.rm=TRUE)

# Need to convert trait names from species names to tipnamecodes
name.map.file <- read.csv('../../Species_name_map_uids.csv')
name.map <- name.map.file$tipnamecodes
names(name.map) <- name.map.file$aos.howardmoore.species
name.map.used <- name.map[name.map %in% tree$tip.label]
name.map.used.unique <- name.map.used[!duplicated(name.map.used)]
newnames <- as.character(name.map.used.unique[names(trait)])
names(trait) <- newnames
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
#trait <- trait[1:30]
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))

# Remove Old World and/or Nearctic species
region.data <- read.table("tyranni/other_data/Range_data_AOSHM_Olson_broad.txt")
areas <- cbind(region.data[,3:11])
rownames(areas) <- as.character(region.data$V1)
region.names <- c("WI", "AM", "AN", "NA", "OW", "PA", "DT", "AF", "CA")
colnames(areas) <- region.names
NewWorld <- areas[!colnames(areas)[apply(areas,1,which.max)] == "OW",]
subtree <- drop.tip(subtree, setdiff(subtree$tip.label, rownames(NewWorld)))
trait <- trait[subtree$tip.label]

# Convert trait to binary
hist(trait)
mean(trait)
median(trait)
max(trait)
# Suite of thresholds between species-rich/-poor
#thresholds <- c(100,200,300,400,500,600) 
thresholds <-  seq(max(trait)/10, max(trait)-(max(trait)/10), by=max(trait)/10) 
thresholds 

sink("HiSSE_NewWorld_parameter_estimates.txt")
all.res <- list()

for(i in 1:length(thresholds)) {
	
	binary.trait <- trait
	binary.trait[trait > thresholds[i]] <- 1
	binary.trait[trait <= thresholds[i]] <- 0

	# HiSSE expects a 2-column data frame with names in col1 and trait in col2
	trait.table <- data.frame(cbind(names(binary.trait), binary.trait))

	# HiSSE, first on a null model where rates vary independently of the trait

	trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
	trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))
	trans.rates.hisse[!is.na(trans.rates.hisse) & !trans.rates.hisse == 0] = 1
	null.2.hisse <- hisse(subtree, trait.table, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.hisse)

	null.logL <- null.2.hisse$loglik
	null.AIC <- null.2.hisse$AIC

	# HiSSE on a model where rates vary with the trait

	trans.rates.bisse <- TransMatMaker(hidden.states=FALSE)

	pp.bisse.no.hidden <- hisse(subtree, trait.table, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse, output.type="raw")
	pp.bisse.no.hidden # This is the TDD result summary
	print(paste0("Threshold = ", thresholds[i]))
	print(pp.bisse.no.hidden)

	# Compare this model to the null model implemented in HiSSE

	alt.logL <- pp.bisse.no.hidden$loglik
	alt.AIC <- pp.bisse.no.hidden$AIC

	logL <- c(null.logL, alt.logL)
	AIC <- c(null.AIC, alt.AIC)

	res <- as.data.frame(logL, row.names = c("null", "BiSSE"))
	res$AIC <- AIC

	# For kicks, can compare to the BiSSE null model in which rates do not vary
	trans.rates.bisse <- TransMatMaker(hidden.states=FALSE)
	bisse.null <- hisse(subtree, trait.table, hidden.states=FALSE, turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse, output.type="raw")

	bisse.null.logL <- bisse.null$loglik
	bisse.null.AIC <- bisse.null$AIC

	res <- rbind(data.frame(logL = bisse.null.logL, AIC = bisse.null.AIC), res)

	row.names(res) <- c("bisse null", "hisse null", "richness dependent")

	# But what about HiSSE on an additional, unobserved state?

	trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
	trans.rates.hisse <- ParDrop(trans.rates.hisse, c(2,3,5,7,8,9,10,12))

	hisse.tyranni <- hisse(subtree, trait.table, hidden.states=TRUE, turnover.anc=c(1,2,0,3), eps.anc=c(1,2,0,3), trans.rate=trans.rates.hisse, output.type="raw")

	hisse.logL <- hisse.tyranni$loglik
	hisse.AIC <- hisse.tyranni$AIC
	res <- rbind(res, data.frame(logL = hisse.logL, AIC = hisse.AIC, row.names = "richness with hidden state"))
	all.res[[i]] <- res # Add to list of results

}

sink()

threshold.list <- rep(thresholds, rep(4, length(thresholds)))
write.csv(cbind(threshold.list, do.call(rbind, all.res)), file = "HiSSE_NewWorld_model_comparison_results.txt")