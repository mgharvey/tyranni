setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/climate/')
getwd()

options(scipen=999)

library(ape)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(raster)

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

# Also drop tips without range data
no.ranges <- as.character(unique(name.map.c[is.na(name.map.c$iucn.species),]$tipnamecodes))
no.ranges <- no.ranges[no.ranges %in% intree$tip.label]
range.tree <- drop.tip(intree, no.ranges)
#write.tree(range.tree, "~/Documents/research/Tyranni/v2/map_diversity/range_tree.tre")


# Subset table name.map.c to unique tips present in tree
name.map <- name.map.c$aos.howardmoore.species
names(name.map) <- name.map.c$tipnamecodes

# Bio1 is temperature, Bio12 is annual precip (see: http://chelsa-climate.org/bioclim/)
# Bio6 min temp coldest month, bio17 precip of driest quarter
# Bio5 max temp warmest month, bio16 precip wettest quarter
# Bio4 is temperature seasonality, Bio15 precip seasonality

var.names <- vector()
# Contemporary Bioclim data (1979-2013)
# LGM (ca. 21 Ky)
# MIS19 (ca. 785 Ky)
# mid-Pliocene warm period (3.264-3.025 Ma)
# Pliocene M2 (ca. 3.3 Ma)

for(i in c(1,4,5,6,12,15,16,17)) {
	assign(paste0("cont", i), raster(paste0("./bio2_5m/CHELSA_cur_V1_2/bio_", i, ".tif")))
	assign(paste0("lgm", i), raster(paste0("./bio2_5m/chelsa_LGM_v1_2/bio_", i, ".tif")))
	var.names <- c(var.names, paste0("cont", i, "mean"))
	var.names <- c(var.names, paste0("cont", i, "median"))
	var.names <- c(var.names, paste0("cont", i, "var"))
	var.names <- c(var.names, paste0(i, ".lgm.cont.change.mean"))
	var.names <- c(var.names, paste0(i, ".lgm.cont.change.median"))
	if(i %in% c(1,4,12,15,16,17)) {	
		assign(paste0("mis", i), raster(paste0("./bio2_5m/MIS19_v1/bio_", i, ".tif")))
		assign(paste0("plw", i), raster(paste0("./bio2_5m/mPWP_v1/bio_", i, ".tif")))
		assign(paste0("plm", i), raster(paste0("./bio2_5m/M2_v1/bio_", i, ".tif")))
		var.names <- c(var.names, paste0(i, ".total.change.avg.mean"))
		var.names <- c(var.names, paste0(i, ".total.change.avg.median"))
	}
}
var.names

alt <- getData("worldclim", var="alt", res=2.5)
proj4string(alt) <- proj4string(cont1)

elev.table <- read.csv("../elevation/nature25794-s3b.csv", header=TRUE)
elev.min <- vector()
elev.max <- vector()
uniq.species <- unique(elev.table$Species)
for(i in 1:length(uniq.species)) {
	sub.table <- elev.table[elev.table$Species == uniq.species[i],]
	overall.min <- min(sub.table$MinElev)
	overall.max <- max(sub.table$MaxElev)
	elev.max <- c(elev.max, overall.max)
	elev.min <- c(elev.min, overall.min)
}
names(elev.min) <- as.character(uniq.species)
names(elev.max) <- as.character(uniq.species)

# Create a function for trimming one raster to another
trim.to.range <- function(x,y) {
	x2 <- resample(x, y)
	cont1.mod <- overlay(x2, y, fun = function(x3, y2) {
		x3[is.na(y2[])] <- NA
		return(x3)
	})
}

results <- vector()
species <- vector()
for(i in 1:length(range.tree$tip.label)) {
		
	print(paste0(i, ": ", range.tree$tip.label[i]))
	if(as.character(name.map[range.tree$tip.label[i]]) %in% species) {
		
		sp.results <- results[which(species == as.character(name.map[range.tree$tip.label[i]]))[1],]
		
	} else {

		# Get range
		range.data <- readOGR(dsn= paste0('../range_maps/range_maps_shape_AOSHM/', range.tree$tip.label[i], '/'), layer=range.tree$tip.label[i])	
		range.raw <- spTransform(range.data, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
		range.raw <- gBuffer(range.raw, byid=TRUE, width=0)
		range <- try(unionSpatialPolygons(range.raw, rep(as.character(i), length(range.raw@polygons))))
		if(class(range) == "try-error") { # If still an error, simplify shapes
			print(paste("Range error for: ", range.tree$tip.label[i]))
			range.raw <- gSimplify(range.raw, tol=0.00001)
			range <- try(unionSpatialPolygons(range.raw, rep(as.character(i), length(range.raw@polygons))))
		} 
		range <- gBuffer(range, byid=TRUE, width=0)
		
		# Create a raster containing cells within the range and elevation range of this species
		this.species <- as.character(name.map[range.tree$tip.label[i]])
		if(this.species %in% names(elev.min)) {			
			sub.alt <- crop(alt, range)
			sub.alt.trimmed <- mask(sub.alt, range)
			sub.alt.trimmed <- calc(sub.alt.trimmed, fun=function(x){ x[x > elev.max[this.species]] <- NA; return(x)} )
			sub.alt.trimmed <- calc(sub.alt.trimmed, fun=function(x){ x[x < elev.min[this.species]] <- NA; return(x)} )				
		} else {
			sub.alt <- crop(alt, range)
			sub.alt.trimmed <- mask(sub.alt, range)
		}
		
		# Get Bio data for species
		sp.results <- vector()
		for(j in c(1,4,5,6,12,15,16,17)) {
			cont <- trim.to.range(get(paste0("cont",j)), sub.alt.trimmed)
			vals1 <- getValues(cont)
			vals1 <- vals1[!is.na(vals1)]
			sp.results <- c(sp.results, mean(vals1), median(vals1), var(vals1))
			lgm <- trim.to.range(get(paste0("lgm",j)), sub.alt.trimmed)
			cont.lgm.diff <- abs(lgm-cont)
			cont.lgm.diff.vals <- getValues(cont.lgm.diff)
			sp.results <- c(sp.results, mean(cont.lgm.diff.vals, na.rm=TRUE), median(cont.lgm.diff.vals, na.rm=TRUE))
		
			if(j %in% c(1,4,12,15,16,17)) {	
				mis <- trim.to.range(get(paste0("mis",j)), sub.alt.trimmed)
				plw <- trim.to.range(get(paste0("plw",j)), sub.alt.trimmed)
				plm <- trim.to.range(get(paste0("plm",j)), sub.alt.trimmed)

				lgm.mis.diff <- abs(mis-lgm)
				mis.plw.diff <- abs(plw-mis)
				plw.plm.diff <- abs(plm-plw)
		
				#plot(unlist(extract(resampled2, range))~unlist(extract(get(paste0("cont",j)), range)))
		
				all.diffs <- mean(cont.lgm.diff, lgm.mis.diff, mis.plw.diff, plw.plm.diff, na.rm=TRUE)
				all.diffs.vals <- unlist(extract(all.diffs, range))
		
				sp.results <- c(sp.results, mean(all.diffs.vals, na.rm=TRUE), median(all.diffs.vals, na.rm=TRUE))
		}
					
			#par(mfrow=c(1,5))
			#plot(cont)
			#plot(lgm)
			#plot(mis)
			#plot(plw)
			#plot(plm)
			
		}
	}
	
	species <- c(species, as.character(name.map[range.tree$tip.label[i]]))
	print(sp.results)
	results <- rbind(results, sp.results)

}
colnames(results) <- var.names
rownames(results) <- range.tree$tip.label

# Write statistics to file
write.table(results, "Climate_data_3.txt")

