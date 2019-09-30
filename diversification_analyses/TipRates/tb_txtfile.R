setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/summary_stats/tb/text_files')
getwd()

library(ape)
library(phytools)

# AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)

# Calculate terminal edge lengths
n <- length(tree$tip.label)
# based on post on Liam Revell's blog:
tb <- setNames(tree$edge.length[sapply(1:n, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)	
tb <- log(tb[tree$tip.label]) # log transform	

write.table(tb, "tb_AOSHM.txt")

# Dups Dropped

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')
tree <- force.ultrametric(tree)

# Calculate terminal edge lengths
n <- length(tree$tip.label)
# based on post on Liam Revell's blog:
tb <- setNames(tree$edge.length[sapply(1:n, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)	
tb <- log(tb[tree$tip.label]) # log transform	

write.table(tb, "tb_dupsdropped.txt")

# AOSC

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')
tree <- force.ultrametric(tree)

# Calculate terminal edge lengths
n <- length(tree$tip.label)
# based on post on Liam Revell's blog:
tb <- setNames(tree$edge.length[sapply(1:n, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)	
tb <- log(tb[tree$tip.label]) # log transform	

write.table(tb, "tb_AOSC.txt")

# IUCN

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')
tree <- force.ultrametric(tree)

# Calculate terminal edge lengths
n <- length(tree$tip.label)
# based on post on Liam Revell's blog:
tb <- setNames(tree$edge.length[sapply(1:n, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)	
tb <- log(tb[tree$tip.label]) # log transform	

write.table(tb, "tb_IUCN.txt")

# 1 My

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_1My_collapsed.tre')
tree <- force.ultrametric(tree)

# Calculate terminal edge lengths
n <- length(tree$tip.label)
# based on post on Liam Revell's blog:
tb <- setNames(tree$edge.length[sapply(1:n, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)	
tb <- log(tb[tree$tip.label]) # log transform	

write.table(tb, "tb_1My.txt")

# 2 My

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_2My_collapsed.tre')
tree <- force.ultrametric(tree)

# Calculate terminal edge lengths
n <- length(tree$tip.label)
# based on post on Liam Revell's blog:
tb <- setNames(tree$edge.length[sapply(1:n, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)	
tb <- log(tb[tree$tip.label]) # log transform	

write.table(tb, "tb_2My.txt")

# HGAPF AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)

# Calculate terminal edge lengths
n <- length(tree$tip.label)
# based on post on Liam Revell's blog:
tb <- setNames(tree$edge.length[sapply(1:n, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)	
tb <- log(tb[tree$tip.label]) # log transform	

write.table(tb, "tb_HGAPF_AOSHM.txt")

# Astral AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)

# Calculate terminal edge lengths
n <- length(tree$tip.label)
# based on post on Liam Revell's blog:
tb <- setNames(tree$edge.length[sapply(1:n, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)	
tb <- log(tb[tree$tip.label]) # log transform	

write.table(tb, "tb_Astral_AOSHM.txt")

# Astral AOSHM missingdropped

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
tree <- force.ultrametric(tree)

# Calculate terminal edge lengths
n <- length(tree$tip.label)
# based on post on Liam Revell's blog:
tb <- setNames(tree$edge.length[sapply(1:n, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)	
tb <- log(tb[tree$tip.label]) # log transform	

write.table(tb, "tb_Astral_AOSHM_missingdropped.txt")

