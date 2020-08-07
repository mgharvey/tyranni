setwd('/Users/eebuser/Documents/Harvey/research/Tyranni/tyranni_revision/BEAST')
getwd()

library(ape)
library(phytools)
library(phylotate)

# Get family names
sample.data <- read.csv("../../Species_name_map_uids_Jetz.csv", header=TRUE)
family.map <- sample.data$howardmoore.family
names(family.map) <- sample.data$tipnamecodes

# Get primary tree
primary <- read_annotated('TyranniTreewithBars.tre')
#primary$node.comment # CI

# Get beast tree
beast <- read_annotated('TyranniFamilies50kMCC.tre')
beast <- ladderize(beast)
#beast$node.comment # CI

# Trim primary to tips on BEAST tree
primary.sub <- drop.tip(primary, setdiff(primary$tip.label, beast$tip.label))
primary.sub <- ladderize(primary.sub)

# Cophylogeny plot
assoc <- cbind(primary.sub$tip.label, primary.sub$tip.label)
obj <- cophylo(primary.sub, beast, rotate=TRUE, assoc=assoc, lwd=0.1) # Rotating is slow

# This is some code added just to plot the two trees with family names, and output those names for Illustrator
examl.fam <- obj$tree[[1]]
examl.fam$tip.label <- as.character(family.map[examl.fam$tip.label])
examl.fam$tip.label[is.na(examl.fam$tip.label)] <- "Psittacidae"
beast.fam <- obj$tree[[2]]
beast.fam$tip.label <- as.character(family.map[beast.fam$tip.label])
beast.fam$tip.label[is.na(beast.fam$tip.label)] <- "Psittacidae"
#par(mfrow=c(1,2))
#plot(examl.fam, label.offset=1)
#plot(beast.fam, label.offset=1)
#par(mfrow=c(1,1))
#pdf("family_list.pdf")
#plot(examl.fam, label.offset=1)
#dev.off()

#pdf("cophylogeny_primary_beast.pdf", width=6, height=6)
#plot(obj, tip.lty=9, link.lty=1, tip.lwd=0.1, link.lwd = 0.1, pts=FALSE, ftype="off")
#dev.off()

# Get all hpd info for primary tree
split.hpd.strings <- strsplit(gsub("\\}","",gsub("\\&age_interval=\\{","", primary$node.comment)),",")
all.mins <- unlist(lapply(split.hpd.strings, "[", 1))
all.maxs <- unlist(lapply(split.hpd.strings, "[", 2))
all.avgs <- c(rep(0, Ntip(primary)), branching.times(primary))

# Subset to nodes present in primary subtree and get min/max values for those
nodes.retained <- matchNodes(primary.sub, primary, method="distances")
nodes.retained[,2]
primary.min <- as.numeric(all.mins[nodes.retained[,2]])
primary.max <- as.numeric(all.maxs[nodes.retained[,2]])
primary.avg <- as.numeric(all.avgs[nodes.retained[,2]])

par(mfrow=c(1,1))
plotTree.errorbars(primary.sub, cbind(primary.min, primary.max))

# CIs for beast tree
pat <- "height_95%_HPD=\\{.+,.+\\},height_median"
hpd.strings <- regmatches(beast$node.comment, regexpr(pat, beast$node.comment))
split.hpd.strings <- strsplit(gsub("\\},height_median","",gsub("height_95%_HPD=\\{","", hpd.strings)),",")
beast.min <- as.numeric(unlist(lapply(split.hpd.strings, "[", 1)))[(Ntip(beast)+1):length(split.hpd.strings)]
beast.max <- as.numeric(unlist(lapply(split.hpd.strings, "[", 2)))[(Ntip(beast)+1):length(split.hpd.strings)]
beast.avg <- branching.times(beast)

plotTree.errorbars(beast, cbind(beast.min, beast.max))

# Match nodes as much as possible (given topological difference)
node.map <- matchNodes(primary.sub, beast, method="descendants")

# Reorganize min/max values for plotting
times.table <- data.frame(cbind(node.map, primary.min, beast.min[node.map[,2]-Ntip(beast)], primary.max, beast.max[node.map[,2]-Ntip(beast)], primary.avg, beast.avg[node.map[,2]-Ntip(beast)]))
times.table <- times.table[-1,] # Get rid of root divergence
times.table <- cbind(times.table, 1:nrow(times.table))
missing.nodes <- which(is.na(times.table[,2]))
times.table <- times.table[-missing.nodes,] # Get rid of columns where nodes don't align
colnames(times.table) <- c("tr1","tr2","primary.min","beast.min","primary.max","beast.max","primary.avg","beast.avg","node.order")

times.table.reorg <- data.frame(cbind(c(times.table$primary.min, times.table$beast.min), c(times.table$primary.max, times.table$beast.max), c(times.table$primary.avg, times.table$beast.avg), c(times.table$node.order, times.table$node.order), c(rep(1, nrow(times.table)), rep(2, nrow(times.table))))) # 1=primary, 2=beast
colnames(times.table.reorg) <- c("min", "max", "avg", "node", "tree")
times.table.reorg <- times.table.reorg[order(times.table.reorg$node),]

#pdf("relative_divtimes_pointvalues.pdf", width=9.5, height=6)
layout(matrix(c(1,2,1,2), ncol=2, byrow=TRUE), widths=c(2,1.4))
x.pos <- 1:(nrow(times.table.reorg)*1.5)
plot(x.pos[x.pos%%3!=0], times.table.reorg$avg, ylim=c(0, max(times.table.reorg$max)), col=NULL, xaxt="n", ylab="Divergence Time (My)", xlab="Node")
segments(x.pos[x.pos%%3!=0], times.table.reorg$min, x.pos[x.pos%%3!=0], times.table.reorg$max, col=c("blue","orange"), lwd=4)
points(x.pos[x.pos%%3!=0], times.table.reorg$avg, ylim=c(0, max(times.table.reorg$max)), pch=16, col="black")
axis(1, at=seq(1.5, nrow(times.table.reorg)*1.5, by=3), labels=times.table$node.order, cex.axis=.6)
legend("topright", c("ExaML+TreePL", "BEAST"), col=c("blue", "orange"), lty=1, lwd=4, bty="n")

# Plot a tree showing the node map
primary.sub.node.lab <- makeNodeLabel(primary.sub, prefix="")
primary.sub.node.lab <- drop.tip(primary.sub.node.lab, "Aratinga_L38084")
primary.sub.node.lab$node.label <- as.numeric(primary.sub.node.lab$node.label)-1
primary.sub.node.lab$tip.label <- as.character(family.map[primary.sub.node.lab$tip.label])
plot(primary.sub.node.lab, label.offset=1, show.node.label=TRUE)
#dev.off()


