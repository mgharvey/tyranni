setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/sr/tree_plot')
getwd()

library(ape)
library(caper)
library(BAMMtools)
library(monogeneaGM)
library(TeachingDemos)
source("../../colorMap.R")

# AOSHM

timetree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')

# Process tree with raw branch lengths (non-ultrametric)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_result.CML_T400FcompleteT0010.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
# Drop low quality samples
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/TreePL/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual)
tree <- ladderize(tree)

# Remove tips not in timetree
c.tree <- drop.tip(tree, setdiff(tree$tip.label, timetree$tip.label))
timetree <- drop.tip(timetree, setdiff(timetree$tip.label, tree$tip.label))

#par(mfrow=c(1,2))
#plot(timetree, show.tip.label=FALSE)
#plot(c.tree, show.tip.label=FALSE)
#c.tree$tip.label[!(c.tree$tip.label %in% timetree$tip.label)]
#all.equal.phylo(c.tree, timetree, use.edge.length=TRUE)
#all.equal.phylo(c.tree, timetree, use.edge.length=FALSE)
#all.equal.phylo(c.tree, timetree, use.edge.length=FALSE, index.return=TRUE)

comp <- all.equal.phylo(c.tree, timetree, use.edge.length=FALSE, index.return=TRUE)
#comp[,1] # timetree
#comp[,2] # c.tree
node.map <- comp[,2]
names(node.map) <- as.character(comp[,1])

c.tree.lengths <- vector()
timetree.lengths <- vector()
for(i in 1:length(timetree$edge.length)) { # for each branch in ref tree
	c.edge.node1 <- as.numeric(node.map[as.character(timetree$edge[i,][1])]) # find matching node for node1 in target tree
	c.edge.node2 <- as.numeric(node.map[as.character(timetree$edge[i,][2])]) # find matching node for node2 in target tree
	c.tree.index <- intersect(which(c.tree$edge[,1] == c.edge.node1), which(c.tree$edge[,2] == c.edge.node2)) # get that branch in target
	c.tree.edge.length <- c.tree$edge.length[c.tree.index] # get it's length
	c.tree.lengths <- c(c.tree.lengths, c.tree.edge.length)
	timetree.lengths <- c(timetree.lengths, timetree$edge.length[i])
}
ratio <- log(c.tree.lengths/timetree.lengths)

RdYlBu =  rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695"))
class(RdYlBu)
RdYlBu.extended =  rev(c("#a50026","#a50026","#a50026","#a50026",
                 "#a50026","#d73027","#f46d43","#fdae61","#fee090",
                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",
                 "#313695","#313695","#313695","#313695","#313695"))
class(RdYlBu.extended)

# Using ES for each branch
pdf("AOSHM_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(ratio, 64, "sr", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(ratio, "RdYlBu.extended", colorbreaks, logcolor=FALSE, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(timetree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu.extended)(length(ratio)), nticks=5, min(ratio, na.rm=TRUE), max(ratio, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

#nrow(c.tree$edge)
#length(timetree.lengths)
#length(c.tree.lengths)
#plot(c.tree.lengths~timetree.lengths)
#plot(timetree, show.tip.label=FALSE, edge.color=)


# AOSC:

timetree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')

# Process tree with raw branch lengths (non-ultrametric)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_result.CML_T400FcompleteT0010.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual) # Drop low quality samples
tree <- ladderize(tree)

# Remove tips not in timetree
c.tree <- drop.tip(tree, setdiff(tree$tip.label, timetree$tip.label))
timetree <- drop.tip(timetree, setdiff(timetree$tip.label, tree$tip.label))

comp <- all.equal.phylo(c.tree, timetree, use.edge.length=FALSE, index.return=TRUE)
node.map <- comp[,2]
names(node.map) <- as.character(comp[,1])

c.tree.lengths <- vector()
timetree.lengths <- vector()
for(i in 1:length(timetree$edge.length)) { # for each branch in ref tree
	c.edge.node1 <- as.numeric(node.map[as.character(timetree$edge[i,][1])]) # find matching node for node1 in target tree
	c.edge.node2 <- as.numeric(node.map[as.character(timetree$edge[i,][2])]) # find matching node for node2 in target tree
	c.tree.index <- intersect(which(c.tree$edge[,1] == c.edge.node1), which(c.tree$edge[,2] == c.edge.node2)) # get that branch in target
	c.tree.edge.length <- c.tree$edge.length[c.tree.index] # get it's length
	c.tree.lengths <- c(c.tree.lengths, c.tree.edge.length)
	timetree.lengths <- c(timetree.lengths, timetree$edge.length[i])
}
ratio <- log(c.tree.lengths/timetree.lengths)

# Using ES for each branch
pdf("AOSC_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(ratio, 64, "sr", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(ratio, "RdYlBu.extended", colorbreaks, logcolor=FALSE, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(timetree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu.extended)(length(ratio)), nticks=5, min(ratio, na.rm=TRUE), max(ratio, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# IUCN:

timetree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')

# Process tree with raw branch lengths (non-ultrametric)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_result.CML_T400FcompleteT0010.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual) # Drop low quality samples
tree <- ladderize(tree)

# Remove tips not in timetree
c.tree <- drop.tip(tree, setdiff(tree$tip.label, timetree$tip.label))
timetree <- drop.tip(timetree, setdiff(timetree$tip.label, tree$tip.label))

comp <- all.equal.phylo(c.tree, timetree, use.edge.length=FALSE, index.return=TRUE)
node.map <- comp[,2]
names(node.map) <- as.character(comp[,1])

c.tree.lengths <- vector()
timetree.lengths <- vector()
for(i in 1:length(timetree$edge.length)) { # for each branch in ref tree
	c.edge.node1 <- as.numeric(node.map[as.character(timetree$edge[i,][1])]) # find matching node for node1 in target tree
	c.edge.node2 <- as.numeric(node.map[as.character(timetree$edge[i,][2])]) # find matching node for node2 in target tree
	c.tree.index <- intersect(which(c.tree$edge[,1] == c.edge.node1), which(c.tree$edge[,2] == c.edge.node2)) # get that branch in target
	c.tree.edge.length <- c.tree$edge.length[c.tree.index] # get it's length
	c.tree.lengths <- c(c.tree.lengths, c.tree.edge.length)
	timetree.lengths <- c(timetree.lengths, timetree$edge.length[i])
}
ratio <- log(c.tree.lengths/timetree.lengths)

# Using ES for each branch
pdf("IUCN_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(ratio, 64, "sr", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(ratio, "RdYlBu.extended", colorbreaks, logcolor=FALSE, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(timetree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu.extended)(length(ratio)), nticks=5, min(ratio, na.rm=TRUE), max(ratio, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# 1My:

timetree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_1My_collapsed.tre')

# Process tree with raw branch lengths (non-ultrametric)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_result.CML_T400FcompleteT0010.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual) # Drop low quality samples
tree <- ladderize(tree)

# Remove tips not in timetree
c.tree <- drop.tip(tree, setdiff(tree$tip.label, timetree$tip.label))
timetree <- drop.tip(timetree, setdiff(timetree$tip.label, tree$tip.label))

comp <- all.equal.phylo(c.tree, timetree, use.edge.length=FALSE, index.return=TRUE)
node.map <- comp[,2]
names(node.map) <- as.character(comp[,1])

c.tree.lengths <- vector()
timetree.lengths <- vector()
for(i in 1:length(timetree$edge.length)) { # for each branch in ref tree
	c.edge.node1 <- as.numeric(node.map[as.character(timetree$edge[i,][1])]) # find matching node for node1 in target tree
	c.edge.node2 <- as.numeric(node.map[as.character(timetree$edge[i,][2])]) # find matching node for node2 in target tree
	c.tree.index <- intersect(which(c.tree$edge[,1] == c.edge.node1), which(c.tree$edge[,2] == c.edge.node2)) # get that branch in target
	c.tree.edge.length <- c.tree$edge.length[c.tree.index] # get it's length
	c.tree.lengths <- c(c.tree.lengths, c.tree.edge.length)
	timetree.lengths <- c(timetree.lengths, timetree$edge.length[i])
}
ratio <- log(c.tree.lengths/timetree.lengths)

# Using ES for each branch
pdf("1My_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(ratio, 64, "sr", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(ratio, "RdYlBu.extended", colorbreaks, logcolor=FALSE, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(timetree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu.extended)(length(ratio)), nticks=5, min(ratio, na.rm=TRUE), max(ratio, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# 2My:

timetree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_2My_collapsed.tre')

# Process tree with raw branch lengths (non-ultrametric)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_result.CML_T400FcompleteT0010.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual) # Drop low quality samples
tree <- ladderize(tree)

# Remove tips not in timetree
c.tree <- drop.tip(tree, setdiff(tree$tip.label, timetree$tip.label))
timetree <- drop.tip(timetree, setdiff(timetree$tip.label, tree$tip.label))

comp <- all.equal.phylo(c.tree, timetree, use.edge.length=FALSE, index.return=TRUE)
node.map <- comp[,2]
names(node.map) <- as.character(comp[,1])

c.tree.lengths <- vector()
timetree.lengths <- vector()
for(i in 1:length(timetree$edge.length)) { # for each branch in ref tree
	c.edge.node1 <- as.numeric(node.map[as.character(timetree$edge[i,][1])]) # find matching node for node1 in target tree
	c.edge.node2 <- as.numeric(node.map[as.character(timetree$edge[i,][2])]) # find matching node for node2 in target tree
	c.tree.index <- intersect(which(c.tree$edge[,1] == c.edge.node1), which(c.tree$edge[,2] == c.edge.node2)) # get that branch in target
	c.tree.edge.length <- c.tree$edge.length[c.tree.index] # get it's length
	c.tree.lengths <- c(c.tree.lengths, c.tree.edge.length)
	timetree.lengths <- c(timetree.lengths, timetree$edge.length[i])
}
ratio <- log(c.tree.lengths/timetree.lengths)

# Using ES for each branch
pdf("2My_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(ratio, 64, "sr", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(ratio, "RdYlBu.extended", colorbreaks, logcolor=FALSE, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(timetree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu.extended)(length(ratio)), nticks=5, min(ratio, na.rm=TRUE), max(ratio, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# dupsdropped:

timetree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')

# Process tree with raw branch lengths (non-ultrametric)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_result.CML_T400FcompleteT0010.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual) # Drop low quality samples
tree <- ladderize(tree)

# Remove tips not in timetree
c.tree <- drop.tip(tree, setdiff(tree$tip.label, timetree$tip.label))
timetree <- drop.tip(timetree, setdiff(timetree$tip.label, tree$tip.label))

comp <- all.equal.phylo(c.tree, timetree, use.edge.length=FALSE, index.return=TRUE)
node.map <- comp[,2]
names(node.map) <- as.character(comp[,1])

c.tree.lengths <- vector()
timetree.lengths <- vector()
for(i in 1:length(timetree$edge.length)) { # for each branch in ref tree
	c.edge.node1 <- as.numeric(node.map[as.character(timetree$edge[i,][1])]) # find matching node for node1 in target tree
	c.edge.node2 <- as.numeric(node.map[as.character(timetree$edge[i,][2])]) # find matching node for node2 in target tree
	c.tree.index <- intersect(which(c.tree$edge[,1] == c.edge.node1), which(c.tree$edge[,2] == c.edge.node2)) # get that branch in target
	c.tree.edge.length <- c.tree$edge.length[c.tree.index] # get it's length
	c.tree.lengths <- c(c.tree.lengths, c.tree.edge.length)
	timetree.lengths <- c(timetree.lengths, timetree$edge.length[i])
}
ratio <- log(c.tree.lengths/timetree.lengths)

# Using ES for each branch
pdf("dupsdropped_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(ratio, 64, "sr", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(ratio, "RdYlBu.extended", colorbreaks, logcolor=FALSE, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(timetree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu.extended)(length(ratio)), nticks=5, min(ratio, na.rm=TRUE), max(ratio, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# HGAPF:

timetree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')

# Process tree with raw branch lengths (non-ultrametric)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_result.CML_HGAPFcompleteT0013.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual) # Drop low quality samples
tree <- ladderize(tree)

# Remove tips not in timetree
c.tree <- drop.tip(tree, setdiff(tree$tip.label, timetree$tip.label))
timetree <- drop.tip(timetree, setdiff(timetree$tip.label, tree$tip.label))

comp <- all.equal.phylo(c.tree, timetree, use.edge.length=FALSE, index.return=TRUE)
node.map <- comp[,2]
names(node.map) <- as.character(comp[,1])

c.tree.lengths <- vector()
timetree.lengths <- vector()
for(i in 1:length(timetree$edge.length)) { # for each branch in ref tree
	c.edge.node1 <- as.numeric(node.map[as.character(timetree$edge[i,][1])]) # find matching node for node1 in target tree
	c.edge.node2 <- as.numeric(node.map[as.character(timetree$edge[i,][2])]) # find matching node for node2 in target tree
	c.tree.index <- intersect(which(c.tree$edge[,1] == c.edge.node1), which(c.tree$edge[,2] == c.edge.node2)) # get that branch in target
	c.tree.edge.length <- c.tree$edge.length[c.tree.index] # get it's length
	c.tree.lengths <- c(c.tree.lengths, c.tree.edge.length)
	timetree.lengths <- c(timetree.lengths, timetree$edge.length[i])
}
ratio <- log(c.tree.lengths/timetree.lengths)

# Using ES for each branch
pdf("HGAPF_AOSHM_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(ratio, 64, "sr", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(ratio, "RdYlBu.extended", colorbreaks, logcolor=FALSE, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(timetree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu.extended)(length(ratio)), nticks=5, min(ratio, na.rm=TRUE), max(ratio, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# Astral:

timetree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')

# Process tree with raw branch lengths (non-ultrametric)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_TreeFile.BML_T400Fcomplete_20180325PT400F.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual) # Drop low quality samples
tree <- ladderize(tree)

# Remove tips not in timetree
c.tree <- drop.tip(tree, setdiff(tree$tip.label, timetree$tip.label))
timetree <- drop.tip(timetree, setdiff(timetree$tip.label, tree$tip.label))

comp <- all.equal.phylo(c.tree, timetree, use.edge.length=FALSE, index.return=TRUE)
node.map <- comp[,2]
names(node.map) <- as.character(comp[,1])

c.tree.lengths <- vector()
timetree.lengths <- vector()
for(i in 1:length(timetree$edge.length)) { # for each branch in ref tree
	c.edge.node1 <- as.numeric(node.map[as.character(timetree$edge[i,][1])]) # find matching node for node1 in target tree
	c.edge.node2 <- as.numeric(node.map[as.character(timetree$edge[i,][2])]) # find matching node for node2 in target tree
	c.tree.index <- intersect(which(c.tree$edge[,1] == c.edge.node1), which(c.tree$edge[,2] == c.edge.node2)) # get that branch in target
	c.tree.edge.length <- c.tree$edge.length[c.tree.index] # get it's length
	c.tree.lengths <- c(c.tree.lengths, c.tree.edge.length)
	timetree.lengths <- c(timetree.lengths, timetree$edge.length[i])
}
ratio <- log(c.tree.lengths/timetree.lengths)

# Using ES for each branch
pdf("Astral_AOSHM_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(ratio, 64, "sr", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(ratio, "RdYlBu.extended", colorbreaks, logcolor=FALSE, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(timetree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu.extended)(length(ratio)), nticks=5, min(ratio, na.rm=TRUE), max(ratio, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()

# Astral missingdropped:

timetree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')

# Process tree with raw branch lengths (non-ultrametric)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_TreeFile.BML_T400Fcomplete_20180325PT400F.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual) # Drop low quality samples
tree <- ladderize(tree)

# Remove tips not in timetree
c.tree <- drop.tip(tree, setdiff(tree$tip.label, timetree$tip.label))
timetree <- drop.tip(timetree, setdiff(timetree$tip.label, tree$tip.label))

comp <- all.equal.phylo(c.tree, timetree, use.edge.length=FALSE, index.return=TRUE)
node.map <- comp[,2]
names(node.map) <- as.character(comp[,1])

c.tree.lengths <- vector()
timetree.lengths <- vector()
for(i in 1:length(timetree$edge.length)) { # for each branch in ref tree
	c.edge.node1 <- as.numeric(node.map[as.character(timetree$edge[i,][1])]) # find matching node for node1 in target tree
	c.edge.node2 <- as.numeric(node.map[as.character(timetree$edge[i,][2])]) # find matching node for node2 in target tree
	c.tree.index <- intersect(which(c.tree$edge[,1] == c.edge.node1), which(c.tree$edge[,2] == c.edge.node2)) # get that branch in target
	c.tree.edge.length <- c.tree$edge.length[c.tree.index] # get it's length
	c.tree.lengths <- c(c.tree.lengths, c.tree.edge.length)
	timetree.lengths <- c(timetree.lengths, timetree$edge.length[i])
}
ratio <- log(c.tree.lengths/timetree.lengths)

# Using ES for each branch
pdf("Astral_AOSHM_missingdropped_allBranches.pdf", width=8, height=8)
colorbreaks <- assignColorBreaks(ratio, 64, "sr", logcolor=FALSE, method="linear", JenksSubset=20000)
colorobj <- colorMap(ratio, "RdYlBu.extended", colorbreaks, logcolor=FALSE, color.interval=NULL)
#plot(tree, edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
plot(timetree, type="fan", edge.color = colorobj$cols, cex=0.07, show.tip.label=FALSE)
subplot(colorBar(colorRampPalette(RdYlBu.extended)(length(ratio)), nticks=5, min(ratio, na.rm=TRUE), max(ratio, na.rm=TRUE)),
	x=grconvertX(c(0.05,0.12), from='npc'), y=grconvertY(c(0.75,1), from='npc'), type='fig', pars=list( mar=c(1.5,1.5,0,0)+0.1))
dev.off()


