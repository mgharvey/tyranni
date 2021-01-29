setwd("/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/range_maps")
getwd()

require(ape)
require(rgdal)
require(rgeos)
require(raster)
require(geosphere)
require(UScensus2010)
 
# The input file geodatabase
fgdb = "/Users/michaelharvey/Documents/research/Tyranni/v1/ILS_vs_latitude/ranges/BOTW.gdb" # Download from IUCN website
 
# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
fc_list = ogrListLayers(fgdb)
 
# Read the feature class
fc = readOGR(dsn=fgdb,layer="All_Species") # This takes ~1 hour
 
# Name map
name.map.c <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/ExaML_result.CML_T400FcompleteT0010.tre')
tips.to.drop <- c("Lepcol_alblin_LSUMZ7508", "Siptor_strcol_L6202", "Myiopag_cotta_LSUMZ143820", "Rhegma_berxh_LSU81109", "Rhegma_berxh_MZUSP76895") #T400F
tree <- drop.tip(tree, tips.to.drop)
outgroups <- c("Gallus_L36208", "Cariama_Y101002", "Micrastur_L11298", "Nestor_merlis_A11054", "Nestor_notab_A13087", "Aratinga_L38084", "Acanth_chlris_O3RIFL", "Certhia_L64261", "Sitta_L30550", "Lipau_vocran_GUNA27", "Passarella_L30382", "Passer_L52749", "Corvus_brac_LSUMZ53041", "Corvus_coride_A17778", "Artamus_L73536", "Pomast_temlis_A19627", "Chlamy_nuclis_A19590", "Menura_A2360", "Stiltia_isalla_ANWC28635")
intree <- drop.tip(tree, outgroups)
intree <- ladderize(intree)
tree.species <- intree$tip.label

# regions
region.data <- readOGR(dsn="./Olson_regions_broad/", layer = "Olson_regions_broad")
regions <- spTransform(region.data, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#plot(regions[9,])
region.names <- c("WI", "AM", "AN", "NA", "OW", "PA", "DT", "AF", "CA")
# West Indies, Patagonia, Amazon, Atlantic Forest, Dry Triangle, Andes, C. America, N. America, Old World
regions$OBJECTID <- region.names
	
# Subset table name.map.c to unique tips present in tree
name.map.c.unique <- subset(name.map.c[name.map.c$tipnamecodes %in% intree$tip.label,], !duplicated(tipnamecodes))
#unique.range.species <- unique(name.map.c.unique$iucn.species[!is.na(name.map.c.unique$iucn.species)])

# Also drop tips without range data
no.ranges <- as.character(unique(name.map.c[is.na(name.map.c$iucn.species),]$tipnamecodes))
no.ranges <- no.ranges[no.ranges %in% intree$tip.label]
range.tree <- drop.tip(intree, no.ranges)
#write.tree(range.tree, "~/Documents/research/Tyranni/v2/map_diversity/range_tree.tre")

# AOSHM

# Lookup table of names
lookup.map <-  name.map.c.unique$aos.howardmoore.species
names(lookup.map) <- name.map.c.unique$iucn.species

# Loop over all species present in tree
regions.vectors <- vector()
region.areas.vectors <- vector()
for(i in 1:length(range.tree$tip.label)) {

	# How many range.species are there in the phy.species for that tip
	phy.species.i <- unique(name.map.c[name.map.c$tipnamecodes == range.tree$tip.label[i],]$aos.howardmoore.species)
	range.species.i <- as.character(unique(name.map.c[name.map.c$aos.howardmoore.species == phy.species.i,]$iucn.species))
	range.species.i <- range.species.i[!is.na(range.species.i)]
	#range.species.i 
	
	# Get range	
	sp.data <- fc[fc$SCINAME %in% range.species.i,]
	#sp.data@data

	sp.data <- sp.data[sp.data$ORIGIN == 1,] # Native only
			
	# Reduce to breeding and resident distributions
	sp.data <- sp.data[sp.data$SEASONAL <= 2,] # Breeding and resident ranges only
	sp.data <- aggregate(sp.data, by = "SEASONAL")
	
	regions.vector <- rep(0, length(region.names))
	region.areas.vector <- rep(0, length(region.names))
	overlaps <- intersect(sp.data, regions)
	overlaps <- aggregate(overlaps, by = "OBJECTID")

	if(!is.null(overlaps)) { # If found in area polygons
		# Get list of regions in which it occurs
		regions.vector[region.names %in% overlaps$OBJECTID[which(areaPoly(overlaps) > (0.01*sum(areaPoly(sp.data))))]] <- 1  # Areas containing > 1% of sp. range
		# Get areas of occurrence in each region divided by total area of that region
		region.areas.vector[region.names %in% overlaps$OBJECTID] <- areaPoly(overlaps)/areaPoly(regions[which(regions$OBJECTID %in% overlaps$OBJECTID),])
	} else { # Rare (if any) cases where not found within polygons, assign to nearest
		closest.region <- region.names[which.min(gDistance(sp.data, regions, byid=TRUE))]
		regions.vector[region.names %in% closest.region] <- 1
		# Keep region areas all at 0?
	}
	
	regions.vectors <- c(regions.vectors, paste(regions.vector, collapse=""))
	region.areas.vectors <- rbind(region.areas.vectors, region.areas.vector)
	print(paste0(i, ". ", range.tree$tip.label[i], ": ", paste(region.names[which(regions.vector == 1)], collapse=","), " [", paste(region.areas.vector, collapse=","),"]"))
	
}

# Write statistics to file
data <- data.frame(cbind(range.tree$tip.label, regions.vectors), region.areas.vectors)
write.table(data, "Range_data_AOSHM_Olson_broad.txt")

