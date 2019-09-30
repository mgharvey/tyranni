setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/summary_stats/es/tree_plot')
getwd()

library(ape)
library(caper)
library(BAMMtools)
library(monogeneaGM)
library(TeachingDemos)
source("../../colorMap.R")

##################################################

# Color as mean of daughter branches:

##################################################


tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

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

# follow Jetz et al. and color internal branches using mean of descendant branches
rootnode <- length(tree$tip.label)+1
es.jetz <- numeric(nrow(tree$edge))
for (i in 1:nrow(tree$edge)){
	daughters <- clade.members(tree$edge[i,2], tree, tip.labels=TRUE)
	if (length(daughters) == 1) { # If terminal branch
		node <- tree$edge[i,2]
		index <- 1
		qx <- 0
		while (node != rootnode){
			el <- tree$edge.length[tree$edge[,2] == node]
			node <- tree$edge[,1][tree$edge[,2] == node]			
			qx <- qx + el* (1 / 2^(index-1))			
			index <- index + 1
		}
		es.jetz[i] <- log(1/qx)		
	} else { # If internal branch
		es.jetz[i] <- mean(es[names(es) %in% daughters])
	}
}		

# color internal branches using their es values
rootnode <- length(tree$tip.label) + 1
es.all <- numeric(nrow(tree$edge))
for (i in 1:length(es.all)){
	node <- tree$edge[i,2]
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es.all[i] <- log(1/qx)
}	
which(es.all=="Inf")
plot(es.all~es.jetz)

# color internal branches using the inverse of their branch lengths
rootnode <- length(tree$tip.label) + 1
ibl.all <- numeric(nrow(tree$edge))
for (i in 1:length(ibl.all)){
	ibl.all[i] <- log(1/tree$edge.length[i])
}	

RdYlBu =  rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695"))
class(RdYlBu)

# Using ES for each branch
pdf("AOSHM_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(es.all, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(es.all, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu)(length(es.all)), nticks=5, min(es.all, na.rm=TRUE), max(es.all, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# Using Jetz method
#colorbreaks <- assignColorBreaks(es.jetz, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
#colorobj <- colorMap(es.jetz, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)

# Using inverse of branch length
#colorbreaks <- assignColorBreaks(ibl.all, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
#colorobj <- colorMap(ibl.all, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)

# AOSC:

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# color internal branches using their es values
rootnode <- length(tree$tip.label) + 1
es.all <- numeric(nrow(tree$edge))
for (i in 1:length(es.all)){
	node <- tree$edge[i,2]
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es.all[i] <- log(1/qx)
}	

RdYlBu =  rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695"))
class(RdYlBu)

# Using ES for each branch
pdf("AOSC_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(es.all, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(es.all, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu)(length(es.all)), nticks=5, min(es.all, na.rm=TRUE), max(es.all, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# IUCN:

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# color internal branches using their es values
rootnode <- length(tree$tip.label) + 1
es.all <- numeric(nrow(tree$edge))
for (i in 1:length(es.all)){
	node <- tree$edge[i,2]
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es.all[i] <- log(1/qx)
}	

RdYlBu =  rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695"))
class(RdYlBu)

# Using ES for each branch
pdf("IUCN_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(es.all, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(es.all, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu)(length(es.all)), nticks=5, min(es.all, na.rm=TRUE), max(es.all, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# 1My:

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_1My_collapsed.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# color internal branches using their es values
rootnode <- length(tree$tip.label) + 1
es.all <- numeric(nrow(tree$edge))
for (i in 1:length(es.all)){
	node <- tree$edge[i,2]
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es.all[i] <- log(1/qx)
}	

RdYlBu =  rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695"))
class(RdYlBu)

# Using ES for each branch
pdf("1My_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(es.all, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(es.all, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu)(length(es.all)), nticks=5, min(es.all, na.rm=TRUE), max(es.all, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# 2My:

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_2My_collapsed.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# color internal branches using their es values
rootnode <- length(tree$tip.label) + 1
es.all <- numeric(nrow(tree$edge))
for (i in 1:length(es.all)){
	node <- tree$edge[i,2]
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es.all[i] <- log(1/qx)
}	

RdYlBu =  rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695"))
class(RdYlBu)

# Using ES for each branch
pdf("2My_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(es.all, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(es.all, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu)(length(es.all)), nticks=5, min(es.all, na.rm=TRUE), max(es.all, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# dupsdropped:

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# color internal branches using their es values
rootnode <- length(tree$tip.label) + 1
es.all <- numeric(nrow(tree$edge))
for (i in 1:length(es.all)){
	node <- tree$edge[i,2]
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es.all[i] <- log(1/qx)
}	

RdYlBu =  rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695"))
class(RdYlBu)

# Using ES for each branch
pdf("dupsdropped_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(es.all, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(es.all, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu)(length(es.all)), nticks=5, min(es.all, na.rm=TRUE), max(es.all, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# HGAPF:

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# color internal branches using their es values
rootnode <- length(tree$tip.label) + 1
es.all <- numeric(nrow(tree$edge))
for (i in 1:length(es.all)){
	node <- tree$edge[i,2]
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es.all[i] <- log(1/qx)
}	

RdYlBu =  rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695"))
class(RdYlBu)

# Using ES for each branch
pdf("HGAPF_AOSHM_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(es.all, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(es.all, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu)(length(es.all)), nticks=5, min(es.all, na.rm=TRUE), max(es.all, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# Astral:

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# color internal branches using their es values
rootnode <- length(tree$tip.label) + 1
es.all <- numeric(nrow(tree$edge))
for (i in 1:length(es.all)){
	node <- tree$edge[i,2]
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es.all[i] <- log(1/qx)
}	

RdYlBu =  rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695"))
class(RdYlBu)

# Using ES for each branch
pdf("Astral_AOSHM_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(es.all, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(es.all, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu)(length(es.all)), nticks=5, min(es.all, na.rm=TRUE), max(es.all, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# Astral missingdropped:

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# color internal branches using their es values
rootnode <- length(tree$tip.label) + 1
es.all <- numeric(nrow(tree$edge))
for (i in 1:length(es.all)){
	node <- tree$edge[i,2]
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- tree$edge.length[tree$edge[,2] == node]
		node <- tree$edge[,1][tree$edge[,2] == node]			
		qx <- qx + el* (1 / 2^(index-1))			
		index <- index + 1
	}
	es.all[i] <- log(1/qx)
}	

RdYlBu =  rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695"))
class(RdYlBu)

# Using ES for each branch
pdf("Astral_AOSHM_missingdropped_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(es.all, 64, "es", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(es.all, "RdYlBu", colorbreaks, logcolor, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(tree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu)(length(es.all)), nticks=5, min(es.all, na.rm=TRUE), max(es.all, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

