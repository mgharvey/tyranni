setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/es')
getwd()

library(ape)
library(phangorn)
library(phytools)
library(diversitree)
source("../FISSE.R")

# AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
dx.table <- read.table("./text_files/es_AOSHM.txt")
trait.data <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM_Olson_broad.txt")
biogeobears <- cbind(as.character(trait.data$V1), sprintf("%09d", trait.data$regions.vectors)) # Add leading zeros back in, massage into matrix
areas <- cbind(trait.data[,3:11])
rownames(areas) <- as.character(trait.data$V1)

region.names <- c("WI", "AM", "AN", "NA", "OW", "PA", "DT", "AF", "CA")
colnames(areas) <- region.names

res <- vector()
for(i in 1:length(region.names)) {
	for(j in 1:length(region.names)) {
		if(i < j) {			
			
			target.areas <- c(region.names[i], region.names[j]) # Choose areas
			subset.a <- areas[colnames(areas)[apply(areas,1,which.max)] == target.areas[1],]
			subset.b <- areas[colnames(areas)[apply(areas,1,which.max)] == target.areas[2],]

			trait <- c(rep(0, nrow(subset.a)), rep(1, nrow(subset.b)))
			names(trait) <- c(rownames(subset.a), rownames(subset.b))
			subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))
			trait <- trait[subtree$tip.label]
			dx <- as.numeric(dx.table$x)
			names(dx) <- rownames(dx.table)
			dx <- dx[subtree$tip.label]

			f.res.NW <- FISSE.binary.mod(tree, trait, dx, incomplete=TRUE)
			stats <- c(region.names[i], region.names[j], f.res.NW$pval, f.res.NW$lambda0, f.res.NW$lambda1, f.res.NW$null_mean_diff)
			print(stats)
			res <- rbind(res, stats)
			
		}
	}
}

colnames(res) <- c("Region1", "Region2", "P-Value", "Lambda0", "Lambda1", "Null_Mean_Diff")
write.table(res, "Olson_broad_es_comparison.txt")