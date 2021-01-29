setwd("/Users/michaelharvey/Documents/research/Tyranni/v2/map_diversity")
getwd()

library(ape)
library(rgdal)
library(rgeos)
library(raster)
library(sf)


#################
### Load data ###
#################

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/examl/ExaML_result.CML_T400FcompleteT0010.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895") #T400F
tree <- drop.tip(tree, tips.to.drop)
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
intree <- drop.tip(tree, outgroups)
intree <- ladderize(intree)
tree.species <- intree$tip.label

# Range data
sf <- st_read("~/Documents/research/range_evolution/birdlife2017_moll.gpkg")

# Name map
name.map.c <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')

# Set up output directory
fixedDir <- '~/Documents/research/Tyranni/v2/map_diversity/range_maps_shape_AOSHM/'

# Subset table name.map.c to unique tips present in tree
name.map.c.unique <- subset(name.map.c[name.map.c$tipnamecodes %in% intree$tip.label,], !duplicated(tipnamecodes))
#unique.range.species <- unique(name.map.c.unique$iucn.species[!is.na(name.map.c.unique$iucn.species)])

# Also drop tips without range data
no.ranges <- as.character(unique(name.map.c[is.na(name.map.c$iucn.species),]$tipnamecodes))
no.ranges <- no.ranges[no.ranges %in% intree$tip.label]
range.tree <- drop.tip(intree, no.ranges)
#write.tree(range.tree, "~/Documents/research/Tyranni/v2/map_diversity/range_tree.tre")

###############################################################
### Extract portions of range of interest and save to files ###
###############################################################

# Lookup table of names
lookup.map <-  name.map.c.unique$aos.howardmoore.species
names(lookup.map) <- name.map.c.unique$tipnamecodes

# Loop over all species present in tree
for(i in 1:length(range.tree$tip.label)) {

	# How many range.species are there in the phy.species for that tip
	phy.species.i <- unique(name.map.c[name.map.c$tipnamecodes == range.tree$tip.label[i],]$aos.howardmoore.species)
	range.species.i <- as.character(unique(name.map.c[name.map.c$aos.howardmoore.species == phy.species.i,]$iucn.species))
	range.species.i <- range.species.i[!is.na(range.species.i)]
	
	# Are any of those other range.species in the tree? # Not using this unless we revert to IUCN taxonomy
	#if(length(range.species.i) > 1) {
	#	tips.range.species.i <- as.character(name.map.c.unique[name.map.c.unique$range.species %in% range.species.i,]$tipnamecodes)
	#	if(length(tips.range.species.i[!is.na(tips.range.species.i)] %in% range.tree$tip.label) == 1) {
	#		# If not, then combine their ranges for calculations
	#		range.species.i <- range.species.i
	#	} else {
	#		# If so, then use only the range.species for that tip
	#		range.species.i <- as.character(lookup.map[names(lookup.map) == range.tree$tip.label[i]]) # Seems to be error here
	#	}
	#} 
	
	# Get range	
	sp.data.i <- sf[sf$sciname %in% range.species.i,]
	sp.data.i <- sp.data.i[sp.data.i$origin == 1,] # Native only
	sp.data.i <- sp.data.i[sp.data.i$seasonal <= 2,] # Breeding and resident ranges only
	sp.data.i <- st_union(sp.data.i)
	
	out.dir <- paste0(fixedDir, gsub(' ', '_', range.tree$tip.label[i]))
	dir.create(out.dir)
	st_write(sp.data.i, dsn=out.dir, driver="ESRI Shapefile", update=TRUE)

	print(paste0(i, ": ", range.tree$tip.label[i]))
		
}

