setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/summary_stats/nh/text_files')
getwd()

library(ape)
library(phytools)

# AOSHM (maybe just do this, not other taxonomies)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)

# Calculate node heights
rootnode <- length(tree$tip.label) + 1
ht.mean <- numeric(length(tree$tip.label))
ht.median <- numeric(length(tree$tip.label))
for (i in 1:length(ht.mean)){
	node <- i
	qx <- 0
	hts <- vector()
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el	
		hts <- c(hts, qx)	
	}
	ht.mean[i] <- mean(as.numeric(hts))
	ht.median[i] <- median(as.numeric(hts))
}		
names(ht.mean) <- tree$tip.label
names(ht.median) <- tree$tip.label

write.table(cbind(ht.mean, ht.median), "nh_mean_median_AOSHM.txt")

# Dups Dropped

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')
tree <- force.ultrametric(tree)

# Calculate node heights
rootnode <- length(tree$tip.label) + 1
ht.mean <- numeric(length(tree$tip.label))
ht.median <- numeric(length(tree$tip.label))
for (i in 1:length(ht.mean)){
	node <- i
	qx <- 0
	hts <- vector()
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el	
		hts <- c(hts, qx)	
	}
	ht.mean[i] <- mean(as.numeric(hts))
	ht.median[i] <- median(as.numeric(hts))
}		
names(ht.mean) <- tree$tip.label
names(ht.median) <- tree$tip.label

write.table(cbind(ht.mean, ht.median), "nh_mean_median_dupsdropped.txt")

# AOSC

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')
tree <- force.ultrametric(tree)

# Calculate node heights
rootnode <- length(tree$tip.label) + 1
ht.mean <- numeric(length(tree$tip.label))
ht.median <- numeric(length(tree$tip.label))
for (i in 1:length(ht.mean)){
	node <- i
	qx <- 0
	hts <- vector()
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el	
		hts <- c(hts, qx)	
	}
	ht.mean[i] <- mean(as.numeric(hts))
	ht.median[i] <- median(as.numeric(hts))
}		
names(ht.mean) <- tree$tip.label
names(ht.median) <- tree$tip.label

write.table(cbind(ht.mean, ht.median), "nh_mean_median_AOSC.txt")

# IUCN

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')
tree <- force.ultrametric(tree)

# Calculate node heights
rootnode <- length(tree$tip.label) + 1
ht.mean <- numeric(length(tree$tip.label))
ht.median <- numeric(length(tree$tip.label))
for (i in 1:length(ht.mean)){
	node <- i
	qx <- 0
	hts <- vector()
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el	
		hts <- c(hts, qx)	
	}
	ht.mean[i] <- mean(as.numeric(hts))
	ht.median[i] <- median(as.numeric(hts))
}		
names(ht.mean) <- tree$tip.label
names(ht.median) <- tree$tip.label

write.table(cbind(ht.mean, ht.median), "nh_mean_median_IUCN.txt")

# 1 My

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_1My_collapsed.tre')
tree <- force.ultrametric(tree)

# Calculate node heights
rootnode <- length(tree$tip.label) + 1
ht.mean <- numeric(length(tree$tip.label))
ht.median <- numeric(length(tree$tip.label))
for (i in 1:length(ht.mean)){
	node <- i
	qx <- 0
	hts <- vector()
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el	
		hts <- c(hts, qx)	
	}
	ht.mean[i] <- mean(as.numeric(hts))
	ht.median[i] <- median(as.numeric(hts))
}		
names(ht.mean) <- tree$tip.label
names(ht.median) <- tree$tip.label

write.table(cbind(ht.mean, ht.median), "nh_mean_median_1My.txt")

# 2 My

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_2My_collapsed.tre')
tree <- force.ultrametric(tree)

# Calculate node heights
rootnode <- length(tree$tip.label) + 1
ht.mean <- numeric(length(tree$tip.label))
ht.median <- numeric(length(tree$tip.label))
for (i in 1:length(ht.mean)){
	node <- i
	qx <- 0
	hts <- vector()
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el	
		hts <- c(hts, qx)	
	}
	ht.mean[i] <- mean(as.numeric(hts))
	ht.median[i] <- median(as.numeric(hts))
}		
names(ht.mean) <- tree$tip.label
names(ht.median) <- tree$tip.label

write.table(cbind(ht.mean, ht.median), "nh_mean_median_2My.txt")

# HGAPF AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)

# Calculate node heights
rootnode <- length(tree$tip.label) + 1
ht.mean <- numeric(length(tree$tip.label))
ht.median <- numeric(length(tree$tip.label))
for (i in 1:length(ht.mean)){
	node <- i
	qx <- 0
	hts <- vector()
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el	
		hts <- c(hts, qx)	
	}
	ht.mean[i] <- mean(as.numeric(hts))
	ht.median[i] <- median(as.numeric(hts))
}		
names(ht.mean) <- tree$tip.label
names(ht.median) <- tree$tip.label

write.table(cbind(ht.mean, ht.median), "nh_mean_median_HGAPF_AOSHM.txt")

# Astral AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)

# Calculate node heights
rootnode <- length(tree$tip.label) + 1
ht.mean <- numeric(length(tree$tip.label))
ht.median <- numeric(length(tree$tip.label))
for (i in 1:length(ht.mean)){
	node <- i
	qx <- 0
	hts <- vector()
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el	
		hts <- c(hts, qx)	
	}
	ht.mean[i] <- mean(as.numeric(hts))
	ht.median[i] <- median(as.numeric(hts))
}		
names(ht.mean) <- tree$tip.label
names(ht.median) <- tree$tip.label

write.table(cbind(ht.mean, ht.median), "nh_mean_median_Astral_AOSHM.txt")

# Astral AOSHM missingdropped

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
tree <- force.ultrametric(tree)

# Calculate node heights
rootnode <- length(tree$tip.label) + 1
ht.mean <- numeric(length(tree$tip.label))
ht.median <- numeric(length(tree$tip.label))
for (i in 1:length(ht.mean)){
	node <- i
	qx <- 0
	hts <- vector()
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el	
		hts <- c(hts, qx)	
	}
	ht.mean[i] <- mean(as.numeric(hts))
	ht.median[i] <- median(as.numeric(hts))
}		
names(ht.mean) <- tree$tip.label
names(ht.median) <- tree$tip.label

write.table(cbind(ht.mean, ht.median), "nh_mean_median_Astral_AOSHM_missingdropped.txt")

