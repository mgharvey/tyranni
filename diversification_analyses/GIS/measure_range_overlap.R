setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/range_maps/')
getwd()

library(ape)
library(phytools)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(raster)

# Name map
name.map.c <- read.csv('~/Documents/research/Tyranni/Species_name_map_uids.csv')

# Tree
tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)

# Also drop tips without range data
no.ranges <- as.character(unique(name.map.c[is.na(name.map.c$iucn.species),]$tipnamecodes))
no.ranges <- no.ranges[no.ranges %in% tree$tip.label]
range.tree <- drop.tip(tree, no.ranges)
#write.tree(range.tree, "~/Documents/research/Tyranni/v2/map_diversity/range_tree.tre")

# Subset table name.map.c to unique tips present in tree
name.map <- name.map.c$aos.howardmoore.species
names(name.map) <- name.map.c$tipnamecodes

reverse.map <- name.map.c$tipnamecodes
names(reverse.map) <- name.map.c$aos.howardmoore.species

range.files <- list.files('./range_maps_shape_AOSHM/')

species <- unique(as.character(name.map[range.tree$tip.label]))
ranges <- vector()
for(i in 1:length(species)) {
	tips <- name.map.c[name.map.c$aos.howardmoore.species == species[i],]$tipnamecodes
	tips <- tips[!is.na(tips)]
	tip <- as.character(tips[tips %in% range.files][1])
	if(is.null(tip)) {
		print(paste("No range data for: ", species[i]))		
		ranges <- c(ranges, NA)
	} else {
		print(paste0("Getting range data for ", i, ": ", species[i]))
		range.data.i <- readOGR(dsn= paste0('./range_maps_shape_AOSHM/', tip, '/'), layer=tip, verbose=FALSE)	
		range.i.raw <- spTransform(range.data.i, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
		range.i.raw <- gBuffer(range.i.raw, byid=TRUE, width=0)
		range.i <- try(unionSpatialPolygons(range.i.raw, rep(as.character(i), length(range.i.raw@polygons))))
		if(class(range.i) == "try-error") { # If still an error, simplify shapes
			print("Correcting range error")
			range.i.raw <- gSimplify(range.i.raw, tol=0.00001)
			range.i <- try(unionSpatialPolygons(range.i.raw, rep(as.character(i), length(range.i.raw@polygons))))
		} 
		if(class(range.i) == "try-error") { # If still an error, simplify shapes
			print("Range error not corrected, NA inserted")
			range.i <- NA
		} else {
			range.i <- gBuffer(range.i, byid=TRUE, width=0)
		}
		ranges <- c(ranges, range.i)		
	}
}

results <- vector()
for(i in 1:length(species)) {
	print(paste0("Measuring overlap for ", i, ": ", species[i]))			
	j.results <- vector()
	for(j in 1:length(species)) {
		print(j)
		if(i == j) {
			j.result <- NA
		} else {
			i.j.overlap <- try(raster::intersect(ranges[i][[1]], ranges[j][[1]]))
			if(class(i.j.overlap) == "try-error") { # If still an error, simplify shapes
				i.j.overlap <- gIntersection(ranges[i][[1]], ranges[j][[1]])
				if(sum(area(i.j.overlap)) > (0.01*sum(area(ranges[i][[1]])))) {
					j.result <- 1
				} else {
					j.result <- 0	
				}						
			} else {
				if(is.null(i.j.overlap)) {
					j.result <- 0
				} else if(length(i.j.overlap) == 0) { # This is for some weird polygons that fail with intersect()
					i.j.overlap <- gIntersection(ranges[i][[1]], ranges[j][[1]])
					if(class(i.j.overlap) == "SpatialCollections") {
						i.j.overlap <- i.j.overlap@polyobj
					}
					if(sum(area(i.j.overlap)) > (0.01*sum(area(ranges[i][[1]])))) {
						j.result <- 1
					} else {
						j.result <- 0				
					}				
				} else { # For the non-weird ones
					if(sum(area(i.j.overlap)) > (0.01*sum(area(ranges[i][[1]])))) {
						j.result <- 1
					} else {
						j.result <- 0				
					}
				}											
			}
		}
		j.results <- c(j.results, j.result)	
	}
	results <- rbind(results, j.results)
}

rownames(results) <- species
colnames(results) <- species
results
# Write statistics to file
write.table(results, "Range_overlap_data.txt")

