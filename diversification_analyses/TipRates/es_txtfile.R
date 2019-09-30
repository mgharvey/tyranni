setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/summary_stats/es/text_files')
getwd()

library(ape)
library(phytools)

# AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)

# Get ES values
rootnode <- length(tree$tip.label) + 1
es <- numeric(length(tree$tip.label))
for (i in 1:length(es)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es[i] <- log(1/qx)
}		
names(es) <- tree$tip.label

write.table(es, "es_AOSHM.txt")

# Dups Dropped

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')
tree <- force.ultrametric(tree)

# Get ES values
rootnode <- length(tree$tip.label) + 1
es <- numeric(length(tree$tip.label))
for (i in 1:length(es)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es[i] <- log(1/qx)
}		
names(es) <- tree$tip.label

write.table(es, "es_dupsdropped.txt")

# AOSC

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')
tree <- force.ultrametric(tree)

# Get ES values
rootnode <- length(tree$tip.label) + 1
es <- numeric(length(tree$tip.label))
for (i in 1:length(es)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es[i] <- log(1/qx)
}		
names(es) <- tree$tip.label

write.table(es, "es_AOSC.txt")

# IUCN

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')
tree <- force.ultrametric(tree)

# Get ES values
rootnode <- length(tree$tip.label) + 1
es <- numeric(length(tree$tip.label))
for (i in 1:length(es)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es[i] <- log(1/qx)
}		
names(es) <- tree$tip.label

write.table(es, "es_IUCN.txt")

# 1 My

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_1My_collapsed.tre')
tree <- force.ultrametric(tree)

# Get ES values
rootnode <- length(tree$tip.label) + 1
es <- numeric(length(tree$tip.label))
for (i in 1:length(es)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es[i] <- log(1/qx)
}		
names(es) <- tree$tip.label

write.table(es, "es_1My.txt")

# 2 My

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_2My_collapsed.tre')
tree <- force.ultrametric(tree)

# Get ES values
rootnode <- length(tree$tip.label) + 1
es <- numeric(length(tree$tip.label))
for (i in 1:length(es)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es[i] <- log(1/qx)
}		
names(es) <- tree$tip.label

write.table(es, "es_2My.txt")

# HGAPF AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)

# Get ES values
rootnode <- length(tree$tip.label) + 1
es <- numeric(length(tree$tip.label))
for (i in 1:length(es)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es[i] <- log(1/qx)
}		
names(es) <- tree$tip.label

write.table(es, "es_HGAPF_AOSHM.txt")

# Astral AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)

# Get IS values
rootnode <- length(tree$tip.label) + 1
es <- numeric(length(tree$tip.label))
for (i in 1:length(es)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es[i] <- log(1/qx)
}		
names(es) <- tree$tip.label

write.table(es, "es_Astral_AOSHM.txt")

# Astral AOSHM missingdropped

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
tree <- force.ultrametric(tree)

# Get IS values
rootnode <- length(tree$tip.label) + 1
es <- numeric(length(tree$tip.label))
for (i in 1:length(es)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es[i] <- log(1/qx)
}		
names(es) <- tree$tip.label

write.table(es, "es_Astral_AOSHM_missingdropped.txt")

