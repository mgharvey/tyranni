setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/summary_stats/sr/text_files')
getwd()

library(ape)
library(phytools)

# Process tree with raw branch lengths (non-ultrametric)

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_result.CML_T400FcompleteT0010.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual)
tree <- ladderize(tree)

# AOSHM

# Trim to this classification
name.map <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')
aos.howardmoore.species <- name.map$aos.howardmoore.species
names(aos.howardmoore.species) <- name.map$tipnamecodes
tip.species <- aos.howardmoore.species[tree$tip.label]
c.tree <- drop.tip(tree, tree$tip.label[duplicated(tip.species)])

# Get subst rate (root to tip distance) and node density for each tip
rootnode <- length(c.tree$tip.label) + 1
sr <- numeric(length(c.tree$tip.label))
nd <- numeric(length(c.tree$tip.label))
for (i in 1:length(sr)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- c.tree$edge.length[c.tree$edge[,2] == node]
		node <- c.tree$edge[,1][c.tree$edge[,2] == node]			
		qx <- qx + el		
		index <- index + 1
	}
	sr[i] <- qx
	nd[i] <- index
}		
names(sr) <- c.tree$tip.label
names(nd) <- c.tree$tip.label
res <- lm(sr~nd)
sr_nd_resid <- residuals(res)[names(sr)]

write.table(cbind(sr, sr_nd_resid), "sr_AOSHM.txt")

# Dups Dropped

# Trim to this classification
name.map <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')
aos.howardmoore.species <- name.map$aos.howardmoore.species
names(aos.howardmoore.species) <- name.map$tipnamecodes
tips.to.drop <- c("Thamphil_caecen_MZU81155", "Thamphil_caecen_MZU91197", "Thamphil_caecen_MZUSP91661", "Thamphil_caecen_MZUSP91196", "Thamphil_caecen_MZUSP82637", "Thamphil_caecen_MZUSP93176", "Thamphil_caecen_MZUSPML2030", "Thamphil_caecen_LSU34687", "Thamphil_caecen_UWBM77213", "Thamphil_caecen_LSU38284", "Thamphil_caecen_LSU50102", "Thamphil_caecen_LSU1684", "Thamphil_caecen_MZUSPALG085", "Thamphil_caecen_MZUSPALG086", "Thamphil_caecen_MZUSPALG080", "Thamphil_caecen_MZUSPALG020", "Thamphil_caecen_MPEGCPEII037", "Dysith_mental_LSU69289", "Dysith_mental_LSU69439", "Dysith_mental_LSU26463", "Dysith_mental_IAvHBT4874", "Dysith_mental_ANSP19107", "Dysith_mental_L33085", "Dysith_mental_LSU22639", "Dysith_mental_MZUSPALG046", "Dysith_mental_MZUSP85741", "Dysith_mental_MZUSPALG087", "Dysith_mental_MZUSP85740", "Dysith_mental_MZUSP79196", "Dysith_mental_MZUSP79194", "Dysith_mental_MZUSP91760", "Dysith_mental_MZUSP93175", "Dysith_mental_MZUSP79195", "Dysith_mental_MZU94637", "Myiorn_aurris_MZU79943", "Myiorn_aurris_MZU80464", "Myiorn_aurris_MZU79945", "Myiorn_aurris_MZU79944", "Myiorn_aurris_MZUMY16", "Myiorn_aurris_MZUMY02", "Myiorn_aurris_MZUMY05", "Myiorn_aurris_MZUMY12", "Myiorn_aurris_MZUMY18", "Myiorn_aurris_MZUMY17", "Myiorn_aurris_MZUMY03", "Myiorn_aurris_MZUMY04", "Myiorn_aurris_MZUMY13", "Myiorn_albtri_46229", "Myiorn_albtri_46059", "Myiorn_albtri_40513", "Myiorn_albtri_LSU40633", "Myiorn_albtri_58394", "Myiorn_spnov_MZUCAX01", "Myiorn_spnov_MZUMAPI19", "Myiorn_spnov_MZUMAPI11", "Myiorn_spnov_MZUMAPI10", "Myiorn_aurris_MZUBOR1430", "Myiorn_aurris_MZUBOR1429", "Myiorn_aurris_25832", "Myiorn_aurris_MZUCMAL53", "Myiorn_aurris_MZUMY37", "Myiorn_aurris_MZU94449", "Myiorn_aurris_MZUFRAPQP04", "Myiorn_aurris_MZUBOR1431", "Myiorn_aurris_MZUFRAPQP05", "Myiorn_aurris_MZUMY38", "Myiorn_aurris_MZUMY36", "Myiorn_aurris_MZU86003", "Myiorn_aurris_MZUCMAL34", "Myiorn_aurris_MZUMY26", "Myiorn_aurris_MZUCMAL50", "Myiorn_aurris_MZUCMAL46", "Myiorn_aurris_MZUMY23", "Myiorn_aurris_MZUMY24", "Platyr_mysceu_MZUSP91268", "Platyr_mysceu_MZUSP93234", "Platyr_mysceu_MZUSP92485", "Platyr_mysceu_MZUSP92487", "Platyr_mysceu_MZUSP92484", "Platyr_mysceu_MZUSP92167", "Platyr_mysceu_MZUSP93536", "Platyr_mysceu_MZUSP91466", "Platyr_mysceu_MZUSP90974", "Platyr_mysceu_MZUSP90978", "Platyr_mysceu_MZUSP90975", "Platyr_mysceu_MZUSP91465", "Platyr_mysceu_MZUSPALG39", "Platyr_mysceu_MZUSP85805", "Platyr_mysceu_MZUSP85833", "Platyr_mysceu_MZUSPALG159", "Platyr_mysceu_MZUSP85760", "Platyr_mysceu_L28375")
c.tree <- drop.tip(tree, tips.to.drop)

# Get subst rate (root to tip distance) and node density for each tip
rootnode <- length(c.tree$tip.label) + 1
sr <- numeric(length(c.tree$tip.label))
nd <- numeric(length(c.tree$tip.label))
for (i in 1:length(sr)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- c.tree$edge.length[c.tree$edge[,2] == node]
		node <- c.tree$edge[,1][c.tree$edge[,2] == node]			
		qx <- qx + el		
		index <- index + 1
	}
	sr[i] <- qx
	nd[i] <- index
}		
names(sr) <- c.tree$tip.label
names(nd) <- c.tree$tip.label
res <- lm(sr~nd)
sr_nd_resid <- residuals(res)[names(sr)]

write.table(cbind(sr, sr_nd_resid), "sr_dupsdropped.txt")

# AOSC

# Trim to this classification
name.map <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')
aos.clements.species <- name.map$aos.clements.species
names(aos.clements.species) <- name.map$tipnamecodes
tip.species <- aos.clements.species[tree$tip.label]
c.tree <- drop.tip(tree, tree$tip.label[duplicated(tip.species)])

# Get subst rate (root to tip distance) and node density for each tip
rootnode <- length(c.tree$tip.label) + 1
sr <- numeric(length(c.tree$tip.label))
nd <- numeric(length(c.tree$tip.label))
for (i in 1:length(sr)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- c.tree$edge.length[c.tree$edge[,2] == node]
		node <- c.tree$edge[,1][c.tree$edge[,2] == node]			
		qx <- qx + el		
		index <- index + 1
	}
	sr[i] <- qx
	nd[i] <- index
}		
names(sr) <- c.tree$tip.label
names(nd) <- c.tree$tip.label
res <- lm(sr~nd)
sr_nd_resid <- residuals(res)[names(sr)]

write.table(cbind(sr, sr_nd_resid), "sr_AOSC.txt")

# IUCN

# Trim to this classification
name.map <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')
iucn.species <- name.map$iucn.species
names(iucn.species) <- name.map$tipnamecodes
tip.species <- iucn.species[tree$tip.label]
c.tree <- drop.tip(tree, tree$tip.label[duplicated(tip.species)])

# Get subst rate (root to tip distance) and node density for each tip
rootnode <- length(c.tree$tip.label) + 1
sr <- numeric(length(c.tree$tip.label))
nd <- numeric(length(c.tree$tip.label))
for (i in 1:length(sr)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- c.tree$edge.length[c.tree$edge[,2] == node]
		node <- c.tree$edge[,1][c.tree$edge[,2] == node]			
		qx <- qx + el		
		index <- index + 1
	}
	sr[i] <- qx
	nd[i] <- index
}		
names(sr) <- c.tree$tip.label
names(nd) <- c.tree$tip.label
res <- lm(sr~nd)
sr_nd_resid <- residuals(res)[names(sr)]

write.table(cbind(sr, sr_nd_resid), "sr_IUCN.txt")

# 1 My

# Trim time-calibrated tree to tips in our current subst. rate tree
tc.tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')
tc.tree <- force.ultrametric(tc.tree)
tc.tree <- drop.tip(tc.tree, setdiff(tc.tree$tip.label, tree$tip.label))

age <- 1 # Age below which to collapse nodes (My)
h <- nodeHeights(tc.tree)
t <- max(h)-age # time from the root
h1 <- which(h[,1] < t) # identify all edges crossing time x
h2 <- which(h[,2] > t)
ii <- intersect(h1, h2)
nodes <- tc.tree$edge[ii,2] # all daughter nodes of those edges
getDescendants <- phytools:::getDescendants
tips <- lapply(nodes, getDescendants, tree=tc.tree) # find all descendants from each edge
keep <- tc.tree$tip.label[sapply(tips,function(x,y) x[x<=Ntip(y)][1],y=tc.tree)]
c.tree <- drop.tip(tree, setdiff(tree$tip.label, keep))

# Get subst rate (root to tip distance) and node density for each tip
rootnode <- length(c.tree$tip.label) + 1
sr <- numeric(length(c.tree$tip.label))
nd <- numeric(length(c.tree$tip.label))
for (i in 1:length(sr)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- c.tree$edge.length[c.tree$edge[,2] == node]
		node <- c.tree$edge[,1][c.tree$edge[,2] == node]			
		qx <- qx + el		
		index <- index + 1
	}
	sr[i] <- qx
	nd[i] <- index
}		
names(sr) <- c.tree$tip.label
names(nd) <- c.tree$tip.label
res <- lm(sr~nd)
sr_nd_resid <- residuals(res)[names(sr)]

write.table(cbind(sr, sr_nd_resid), "sr_1My.txt")

# 2 My

# Trim time-calibrated tree to tips in our current subst. rate tree
tc.tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')
tc.tree <- force.ultrametric(tc.tree)
tc.tree <- drop.tip(tc.tree, setdiff(tc.tree$tip.label, tree$tip.label))

age <- 2 # Age below which to collapse nodes (My)
h <- nodeHeights(tc.tree)
t <- max(h)-age # time from the root
h1 <- which(h[,1] < t) # identify all edges crossing time x
h2 <- which(h[,2] > t)
ii <- intersect(h1, h2)
nodes <- tc.tree$edge[ii,2] # all daughter nodes of those edges
getDescendants <- phytools:::getDescendants
tips <- lapply(nodes, getDescendants, tree=tc.tree) # find all descendants from each edge
keep <- tc.tree$tip.label[sapply(tips,function(x,y) x[x<=Ntip(y)][1],y=tc.tree)]
c.tree <- drop.tip(tree, setdiff(tree$tip.label, keep))

# Get subst rate (root to tip distance) and node density for each tip
rootnode <- length(c.tree$tip.label) + 1
sr <- numeric(length(c.tree$tip.label))
nd <- numeric(length(c.tree$tip.label))
for (i in 1:length(sr)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- c.tree$edge.length[c.tree$edge[,2] == node]
		node <- c.tree$edge[,1][c.tree$edge[,2] == node]			
		qx <- qx + el		
		index <- index + 1
	}
	sr[i] <- qx
	nd[i] <- index
}		
names(sr) <- c.tree$tip.label
names(nd) <- c.tree$tip.label
res <- lm(sr~nd)
sr_nd_resid <- residuals(res)[names(sr)]

write.table(cbind(sr, sr_nd_resid), "sr_2My.txt")

# HGAPF AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_result.CML_HGAPFcompleteT0013.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495", "Acropt_ortnyx_PhAMC1246")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual)
tree <- ladderize(tree)

# Trim to this classification
name.map <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')
aos.howardmoore.species <- name.map$aos.howardmoore.species
names(aos.howardmoore.species) <- name.map$tipnamecodes
tip.species <- aos.howardmoore.species[tree$tip.label]
c.tree <- drop.tip(tree, tree$tip.label[duplicated(tip.species)])

# Get subst rate (root to tip distance) and node density for each tip
rootnode <- length(c.tree$tip.label) + 1
sr <- numeric(length(c.tree$tip.label))
nd <- numeric(length(c.tree$tip.label))
for (i in 1:length(sr)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- c.tree$edge.length[c.tree$edge[,2] == node]
		node <- c.tree$edge[,1][c.tree$edge[,2] == node]			
		qx <- qx + el		
		index <- index + 1
	}
	sr[i] <- qx
	nd[i] <- index
}		
names(sr) <- c.tree$tip.label
names(nd) <- c.tree$tip.label
res <- lm(sr~nd)
sr_nd_resid <- residuals(res)[names(sr)]

write.table(cbind(sr, sr_nd_resid), "sr_HGAPF_AOSHM.txt")

# Astral AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_TreeFile.BML_T400Fcomplete_20180325PT400F.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895", "Empinax_vircen_L64495", "Acropt_ortnyx_PhAMC1246")
tree <- drop.tip(tree, tips.to.drop)
tree <- root(tree, "Gallus_L36208")
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
tree <- drop.tip(tree, outgroups)
low.qual <- as.character(read.table("~/Documents/research/Tyranni/v2/tree_building/trim_long_terminal_branches/Terminal_branches_trimmed.txt")$low.qual)
tree <- drop.tip(tree, low.qual)
tree <- ladderize(tree)

# Trim to this classification
name.map <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')
aos.howardmoore.species <- name.map$aos.howardmoore.species
names(aos.howardmoore.species) <- name.map$tipnamecodes
tip.species <- aos.howardmoore.species[tree$tip.label]
c.tree <- drop.tip(tree, tree$tip.label[duplicated(tip.species)])

# Get subst rate (root to tip distance) and node density for each tip
rootnode <- length(c.tree$tip.label) + 1
sr <- numeric(length(c.tree$tip.label))
nd <- numeric(length(c.tree$tip.label))
for (i in 1:length(sr)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- c.tree$edge.length[c.tree$edge[,2] == node]
		node <- c.tree$edge[,1][c.tree$edge[,2] == node]			
		qx <- qx + el		
		index <- index + 1
	}
	sr[i] <- qx
	nd[i] <- index
}		
names(sr) <- c.tree$tip.label
names(nd) <- c.tree$tip.label
res <- lm(sr~nd)
sr_nd_resid <- residuals(res)[names(sr)]

write.table(cbind(sr, sr_nd_resid), "sr_Astral_AOSHM.txt")

# Astral AOSHM missingdropped

# Missing information
missing <- read.table("~/Documents/research/Tyranni/v2/tree_building/alignments/T400Fcomplete_missingloci.txt")
tips.to.drop <- as.character(missing[missing$V2 > 250,]$V1)
tree <- drop.tip(tree, tips.to.drop)
tree <- ladderize(tree)

# Trim to this classification
name.map <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')
aos.howardmoore.species <- name.map$aos.howardmoore.species
names(aos.howardmoore.species) <- name.map$tipnamecodes
tip.species <- aos.howardmoore.species[tree$tip.label]
c.tree <- drop.tip(tree, tree$tip.label[duplicated(tip.species)])

# Get subst rate (root to tip distance) and node density for each tip
rootnode <- length(c.tree$tip.label) + 1
sr <- numeric(length(c.tree$tip.label))
nd <- numeric(length(c.tree$tip.label))
for (i in 1:length(sr)){
	node <- i
	index <- 1
	qx <- 0
	while (node != rootnode){
		el <- c.tree$edge.length[c.tree$edge[,2] == node]
		node <- c.tree$edge[,1][c.tree$edge[,2] == node]			
		qx <- qx + el		
		index <- index + 1
	}
	sr[i] <- qx
	nd[i] <- index
}		
names(sr) <- c.tree$tip.label
names(nd) <- c.tree$tip.label
res <- lm(sr~nd)
sr_nd_resid <- residuals(res)[names(sr)]

write.table(cbind(sr, sr_nd_resid), "sr_Astral_AOSHM_missingdropped.txt")

