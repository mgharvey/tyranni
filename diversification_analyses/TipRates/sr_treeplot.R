setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/sr/tree_plot')
getwd()

library(ape)
library(phytools)


# AOSHM:

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# Get es values
sr.table <- read.table("../text_files/sr_AOSHM.txt")
sr <- sr.table$sr
names(sr) <- rownames(sr.table)

# Plot (based on code from Rafael Maia)
bardata <- sr[tree$tip.label]

# plot tree
pdf("AOSHM.pdf", width=8, height=8)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(tree, type='fan', show.tip.label=F, edge.width=0.5, x.lim=c(-70,70), y.lim=c(-70,70), edge.color="gray")
# set angles, convert to radians
angle <- seq(0,360, length.out=Ntip(tree)+1)  # to rotate so it starts counting on right
angle <- angle[-length(angle)] # remove last one so it doesn't overlap with first one
angle <- angle*pi/180 # convert to radians
start <- 45 # set starting point for bars where tips end
scaledat <- 80 # scale to apply so its not disproportional to plotting area
radius <- (bardata*scaledat) + start # "radius" determines size of bars
segments(cos(angle)*start, sin(angle)*start, cos(angle)*radius, sin(angle)*radius, lwd=0.5, col="black", lend=1)
dev.off()

# AOSC:

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# Get es values
sr.table <- read.table("../text_files/sr_AOSC.txt")
sr <- sr.table$sr
names(sr) <- rownames(sr.table)

# Plot (based on code from Rafael Maia)
bardata <- sr[tree$tip.label]

# plot tree
pdf("AOSC.pdf", width=8, height=8)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(tree, type='fan', show.tip.label=F, edge.width=0.5, x.lim=c(-70,70), y.lim=c(-70,70), edge.color="gray")
# set angles, convert to radians
angle <- seq(0,360, length.out=Ntip(tree)+1)  # to rotate so it starts counting on right
angle <- angle[-length(angle)] # remove last one so it doesn't overlap with first one
angle <- angle*pi/180 # convert to radians
start <- 45 # set starting point for bars where tips end
scaledat <- 80 # scale to apply so its not disproportional to plotting area
radius <- (bardata*scaledat) + start # "radius" determines size of bars
segments(cos(angle)*start, sin(angle)*start, cos(angle)*radius, sin(angle)*radius, lwd=0.5, col="black", lend=1)
dev.off()

# IUCN:

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# Get es values
sr.table <- read.table("../text_files/sr_IUCN.txt")
sr <- sr.table$sr
names(sr) <- rownames(sr.table)

# Plot (based on code from Rafael Maia)
bardata <- sr[tree$tip.label]

# plot tree
pdf("IUCN.pdf", width=8, height=8)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(tree, type='fan', show.tip.label=F, edge.width=0.5, x.lim=c(-70,70), y.lim=c(-70,70), edge.color="gray")
# set angles, convert to radians
angle <- seq(0,360, length.out=Ntip(tree)+1)  # to rotate so it starts counting on right
angle <- angle[-length(angle)] # remove last one so it doesn't overlap with first one
angle <- angle*pi/180 # convert to radians
start <- 45 # set starting point for bars where tips end
scaledat <- 80 # scale to apply so its not disproportional to plotting area
radius <- (bardata*scaledat) + start # "radius" determines size of bars
segments(cos(angle)*start, sin(angle)*start, cos(angle)*radius, sin(angle)*radius, lwd=0.5, col="black", lend=1)
dev.off()

# 1My:

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_1My_collapsed.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# Get es values
sr.table <- read.table("../text_files/sr_1My.txt")
sr <- sr.table$sr
names(sr) <- rownames(sr.table)

# Plot (based on code from Rafael Maia)
bardata <- sr[tree$tip.label]

# plot tree
pdf("1My.pdf", width=8, height=8)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(tree, type='fan', show.tip.label=F, edge.width=0.5, x.lim=c(-70,70), y.lim=c(-70,70), edge.color="gray")
# set angles, convert to radians
angle <- seq(0,360, length.out=Ntip(tree)+1)  # to rotate so it starts counting on right
angle <- angle[-length(angle)] # remove last one so it doesn't overlap with first one
angle <- angle*pi/180 # convert to radians
start <- 45 # set starting point for bars where tips end
scaledat <- 80 # scale to apply so its not disproportional to plotting area
radius <- (bardata*scaledat) + start # "radius" determines size of bars
segments(cos(angle)*start, sin(angle)*start, cos(angle)*radius, sin(angle)*radius, lwd=0.5, col="black", lend=1)
dev.off()

# 2My:

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_2My_collapsed.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# Get es values
sr.table <- read.table("../text_files/sr_2My.txt")
sr <- sr.table$sr
names(sr) <- rownames(sr.table)

# Plot (based on code from Rafael Maia)
bardata <- sr[tree$tip.label]

# plot tree
pdf("2My.pdf", width=8, height=8)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(tree, type='fan', show.tip.label=F, edge.width=0.5, x.lim=c(-70,70), y.lim=c(-70,70), edge.color="gray")
# set angles, convert to radians
angle <- seq(0,360, length.out=Ntip(tree)+1)  # to rotate so it starts counting on right
angle <- angle[-length(angle)] # remove last one so it doesn't overlap with first one
angle <- angle*pi/180 # convert to radians
start <- 45 # set starting point for bars where tips end
scaledat <- 80 # scale to apply so its not disproportional to plotting area
radius <- (bardata*scaledat) + start # "radius" determines size of bars
segments(cos(angle)*start, sin(angle)*start, cos(angle)*radius, sin(angle)*radius, lwd=0.5, col="black", lend=1)
dev.off()

# dupsdropped:

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# Get es values
sr.table <- read.table("../text_files/sr_dupsdropped.txt")
sr <- sr.table$sr
names(sr) <- rownames(sr.table)

# Plot (based on code from Rafael Maia)
bardata <- sr[tree$tip.label]

# plot tree
pdf("dupsdropped.pdf", width=8, height=8)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(tree, type='fan', show.tip.label=F, edge.width=0.5, x.lim=c(-70,70), y.lim=c(-70,70), edge.color="gray")
# set angles, convert to radians
angle <- seq(0,360, length.out=Ntip(tree)+1)  # to rotate so it starts counting on right
angle <- angle[-length(angle)] # remove last one so it doesn't overlap with first one
angle <- angle*pi/180 # convert to radians
start <- 45 # set starting point for bars where tips end
scaledat <- 80 # scale to apply so its not disproportional to plotting area
radius <- (bardata*scaledat) + start # "radius" determines size of bars
segments(cos(angle)*start, sin(angle)*start, cos(angle)*radius, sin(angle)*radius, lwd=0.5, col="black", lend=1)
dev.off()

# HGAPF:

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# Get es values
sr.table <- read.table("../text_files/sr_HGAPF_AOSHM.txt")
sr <- sr.table$sr
names(sr) <- rownames(sr.table)

# Plot (based on code from Rafael Maia)
bardata <- sr[tree$tip.label]

# plot tree
pdf("HGAPF_AOSHM.pdf", width=8, height=8)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(tree, type='fan', show.tip.label=F, edge.width=0.5, x.lim=c(-70,70), y.lim=c(-70,70), edge.color="gray")
# set angles, convert to radians
angle <- seq(0,360, length.out=Ntip(tree)+1)  # to rotate so it starts counting on right
angle <- angle[-length(angle)] # remove last one so it doesn't overlap with first one
angle <- angle*pi/180 # convert to radians
start <- 45 # set starting point for bars where tips end
scaledat <- 350 # scale to apply so its not disproportional to plotting area
radius <- (bardata*scaledat) + start # "radius" determines size of bars
segments(cos(angle)*start, sin(angle)*start, cos(angle)*radius, sin(angle)*radius, lwd=0.5, col="black", lend=1)
dev.off()

# Astral:

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# Get es values
sr.table <- read.table("../text_files/sr_Astral_AOSHM.txt")
sr <- sr.table$sr
names(sr) <- rownames(sr.table)

# Plot (based on code from Rafael Maia)
bardata <- sr[tree$tip.label]

# plot tree
pdf("Astral_AOSHM.pdf", width=8, height=8)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(tree, type='fan', show.tip.label=F, edge.width=0.5, x.lim=c(-70,70), y.lim=c(-70,70), edge.color="gray")
# set angles, convert to radians
angle <- seq(0,360, length.out=Ntip(tree)+1)  # to rotate so it starts counting on right
angle <- angle[-length(angle)] # remove last one so it doesn't overlap with first one
angle <- angle*pi/180 # convert to radians
start <- 45 # set starting point for bars where tips end
scaledat <- 80 # scale to apply so its not disproportional to plotting area
radius <- (bardata*scaledat) + start # "radius" determines size of bars
segments(cos(angle)*start, sin(angle)*start, cos(angle)*radius, sin(angle)*radius, lwd=0.5, col="black", lend=1)
dev.off()

# Astral missingdropped:

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
tree <- force.ultrametric(tree)
tree <- ladderize(tree)

# Get es values
sr.table <- read.table("../text_files/sr_Astral_AOSHM_missingdropped.txt")
sr <- sr.table$sr
names(sr) <- rownames(sr.table)

# Plot (based on code from Rafael Maia)
bardata <- sr[tree$tip.label]

# plot tree
pdf("Astral_AOSHM_missingdropped.pdf", width=8, height=8)
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(tree, type='fan', show.tip.label=F, edge.width=0.5, x.lim=c(-70,70), y.lim=c(-70,70), edge.color="gray")
# set angles, convert to radians
angle <- seq(0,360, length.out=Ntip(tree)+1)  # to rotate so it starts counting on right
angle <- angle[-length(angle)] # remove last one so it doesn't overlap with first one
angle <- angle*pi/180 # convert to radians
start <- 45 # set starting point for bars where tips end
scaledat <- 80 # scale to apply so its not disproportional to plotting area
radius <- (bardata*scaledat) + start # "radius" determines size of bars
segments(cos(angle)*start, sin(angle)*start, cos(angle)*radius, sin(angle)*radius, lwd=0.5, col="black", lend=1)
dev.off()



