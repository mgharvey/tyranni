setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/summary_stats/nh/maps')
getwd()

library(raster)
library(maptools)
library(maps)
library(rgdal)
library(rworldmap)
library(sp)

# I was too lazy to replace all instances of "es" with "nh" in script

# AOSHM

# Get es values
es.table <- read.table("../text_files/nh_mean_median_AOSHM.txt")
es <- log(es.table$ht.mean)
names(es) <- rownames(es.table)

# Loop to get range maps
file.NAs <- vector()
polygons <- vector()
my.es <- vector()
for(i in 1:length(es)) {
	print(names(es[i]))
	my.file <- paste0('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/range_maps/range_maps_shape_AOSHM/', names(es[i]), '/', names(es[i]), '.shp')
	if(file.exists(my.file)) {
		range <- readShapeSpatial(my.file)	
		range.dissolved <- try(unionSpatialPolygons(range, rep(as.character(i), length(range@polygons))))
		if(class(range.dissolved) == "try-error") {
			print(paste("Range error for: ", names(es[i])))
		} else {
			polygons <- c(polygons, range.dissolved)			
			my.es <- c(my.es, es[i])
		}
	} else {
		print(paste("No range data for: ", names(es[i])))
		file.NAs <- c(file.NAs, names(es[i]))
	}
}
print(file.NAs)
joined = SpatialPolygons(lapply(polygons, function(x){x@polygons[[1]]}))

# Make an empty raster with 200 km resolution
r <- raster(xmn = -12000000, xmx = 16000000, ymn = -7000000, ymx = 8000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # whole world
#r <- raster(xmn = -12000000, xmx = -2000000, ymn = -3000000, ymx = 3000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # Neotropics only

# For each species, make raster for es and another for weights
es.rasters <- stack()
count.rasters <- stack()
for(i in 1:length(joined)) {
	print(i)
	cell.count <- lengths(extract(r, joined[i]))
	es.raster <- r
	count.raster <- r
	es.raster[joined[i],] <- my.es[i]
	count.raster[joined[i],] <- 1/cell.count
	es.rasters <- stack(es.rasters, es.raster)
	count.rasters <- stack(count.rasters, count.raster)
}
wm <- weighted.mean(x=es.rasters, w=count.rasters, na.rm=TRUE)
library(RColorBrewer)
Pal <- rev(brewer.pal(n=11, c("RdYlBu")))

# Using initial projection of range maps
pdf(file="nh_AOSHM_log.pdf", height=5, width=9)
par(mar=c(1,1,1,5))
proj4string(wm) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline <- readOGR('/Users/michaelharvey/Documents/research/seabirds/sampling_labwork/map/ne_50m_coastline/ne_50m_coastline.shp')
coastline2 <- spTransform(coastline, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
sPDF <- getMap()
sPDF <- spTransform(sPDF, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="")
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black", col="gray")
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="", add=TRUE)
plot(wm, col=Pal, add=TRUE)
plot(coastline2, lwd=0.5, col="black", add=TRUE)
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black")
dev.off()

# AOSC

# Get es values
es.table <- read.table("../text_files/nh_mean_median_AOSC.txt")
es <- log(es.table$ht.mean)
names(es) <- rownames(es.table)

# Loop to get range maps
file.NAs <- vector()
polygons <- vector()
my.es <- vector()
for(i in 1:length(es)) {
	print(names(es[i]))
	my.file <- paste0('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/range_maps/range_maps_shape_AOSC/', names(es[i]), '/', names(es[i]), '.shp')
	if(file.exists(my.file)) {
		range <- readShapeSpatial(my.file)	
		range.dissolved <- try(unionSpatialPolygons(range, rep(as.character(i), length(range@polygons))))
		if(class(range.dissolved) == "try-error") {
			print(paste("Range error for: ", names(es[i])))
		} else {
			polygons <- c(polygons, range.dissolved)			
			my.es <- c(my.es, es[i])
		}
	} else {
		print(paste("No range data for: ", names(es[i])))
		file.NAs <- c(file.NAs, names(es[i]))
	}
}
print(file.NAs)
joined = SpatialPolygons(lapply(polygons, function(x){x@polygons[[1]]}))

# Make an empty raster with 200 km resolution
r <- raster(xmn = -12000000, xmx = 16000000, ymn = -7000000, ymx = 8000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # whole world
#r <- raster(xmn = -12000000, xmx = -2000000, ymn = -3000000, ymx = 3000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # Neotropics only

# For each species, make raster for es and another for weights
es.rasters <- stack()
count.rasters <- stack()
for(i in 1:length(joined)) {
	print(i)
	cell.count <- lengths(extract(r, joined[i]))
	es.raster <- r
	count.raster <- r
	es.raster[joined[i],] <- my.es[i]
	count.raster[joined[i],] <- 1/cell.count
	es.rasters <- stack(es.rasters, es.raster)
	count.rasters <- stack(count.rasters, count.raster)
}
wm <- weighted.mean(x=es.rasters, w=count.rasters, na.rm=TRUE)
library(RColorBrewer)
Pal <- rev(brewer.pal(n=11, c("RdYlBu")))

# Using initial projection of range maps
pdf(file="nh_AOSC_log.pdf", height=5, width=9)
par(mar=c(1,1,1,5))
proj4string(wm) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline <- readOGR('/Users/michaelharvey/Documents/research/seabirds/sampling_labwork/map/ne_50m_coastline/ne_50m_coastline.shp')
coastline2 <- spTransform(coastline, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
sPDF <- getMap()
sPDF <- spTransform(sPDF, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="")
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black", col="gray")
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="", add=TRUE)
plot(wm, col=Pal, add=TRUE)
plot(coastline2, lwd=0.5, col="black", add=TRUE)
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black")
dev.off()

# IUCN

# Get es values
es.table <- read.table("../text_files/nh_mean_median_IUCN.txt")
es <- log(es.table$ht.mean)
names(es) <- rownames(es.table)

# Loop to get range maps
file.NAs <- vector()
polygons <- vector()
my.es <- vector()
for(i in 1:length(es)) {
	print(names(es[i]))
	my.file <- paste0('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/range_maps/range_maps_shape_IUCN/', names(es[i]), '/', names(es[i]), '.shp')
	if(file.exists(my.file)) {
		range <- readShapeSpatial(my.file)	
		range.dissolved <- try(unionSpatialPolygons(range, rep(as.character(i), length(range@polygons))))
		if(class(range.dissolved) == "try-error") {
			print(paste("Range error for: ", names(es[i])))
		} else {
			polygons <- c(polygons, range.dissolved)			
			my.es <- c(my.es, es[i])
		}
	} else {
		print(paste("No range data for: ", names(es[i])))
		file.NAs <- c(file.NAs, names(es[i]))
	}
}
print(file.NAs)
joined = SpatialPolygons(lapply(polygons, function(x){x@polygons[[1]]}))

# Make an empty raster with 200 km resolution
r <- raster(xmn = -12000000, xmx = 16000000, ymn = -7000000, ymx = 8000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # whole world
#r <- raster(xmn = -12000000, xmx = -2000000, ymn = -3000000, ymx = 3000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # Neotropics only

# For each species, make raster for es and another for weights
es.rasters <- stack()
count.rasters <- stack()
for(i in 1:length(joined)) {
	print(i)
	cell.count <- lengths(extract(r, joined[i]))
	es.raster <- r
	count.raster <- r
	es.raster[joined[i],] <- my.es[i]
	count.raster[joined[i],] <- 1/cell.count
	es.rasters <- stack(es.rasters, es.raster)
	count.rasters <- stack(count.rasters, count.raster)
}
wm <- weighted.mean(x=es.rasters, w=count.rasters, na.rm=TRUE)
library(RColorBrewer)
Pal <- rev(brewer.pal(n=11, c("RdYlBu")))

# Using initial projection of range maps
pdf(file="nh_IUCN_log.pdf", height=5, width=9)
par(mar=c(1,1,1,5))
proj4string(wm) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline <- readOGR('/Users/michaelharvey/Documents/research/seabirds/sampling_labwork/map/ne_50m_coastline/ne_50m_coastline.shp')
coastline2 <- spTransform(coastline, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
sPDF <- getMap()
sPDF <- spTransform(sPDF, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="")
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black", col="gray")
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="", add=TRUE)
plot(wm, col=Pal, add=TRUE)
plot(coastline2, lwd=0.5, col="black", add=TRUE)
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black")
dev.off()

# HGAPF AOSHM

# Get es values
es.table <- read.table("../text_files/nh_mean_median_HGAPF_AOSHM.txt")
es <- log(es.table$ht.mean)
names(es) <- rownames(es.table)

# Loop to get range maps
file.NAs <- vector()
polygons <- vector()
my.es <- vector()
for(i in 1:length(es)) {
	print(names(es[i]))
	my.file <- paste0('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/range_maps/range_maps_shape_AOSHM/', names(es[i]), '/', names(es[i]), '.shp')
	if(file.exists(my.file)) {
		range <- readShapeSpatial(my.file)	
		range.dissolved <- try(unionSpatialPolygons(range, rep(as.character(i), length(range@polygons))))
		if(class(range.dissolved) == "try-error") {
			print(paste("Range error for: ", names(es[i])))
		} else {
			polygons <- c(polygons, range.dissolved)			
			my.es <- c(my.es, es[i])
		}
	} else {
		print(paste("No range data for: ", names(es[i])))
		file.NAs <- c(file.NAs, names(es[i]))
	}
}
print(file.NAs)
joined = SpatialPolygons(lapply(polygons, function(x){x@polygons[[1]]}))

# Make an empty raster with 200 km resolution
r <- raster(xmn = -12000000, xmx = 16000000, ymn = -7000000, ymx = 8000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # whole world
#r <- raster(xmn = -12000000, xmx = -2000000, ymn = -3000000, ymx = 3000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # Neotropics only

# For each species, make raster for es and another for weights
es.rasters <- stack()
count.rasters <- stack()
for(i in 1:length(joined)) {
	print(i)
	cell.count <- lengths(extract(r, joined[i]))
	es.raster <- r
	count.raster <- r
	es.raster[joined[i],] <- my.es[i]
	count.raster[joined[i],] <- 1/cell.count
	es.rasters <- stack(es.rasters, es.raster)
	count.rasters <- stack(count.rasters, count.raster)
}
wm <- weighted.mean(x=es.rasters, w=count.rasters, na.rm=TRUE)
library(RColorBrewer)
Pal <- rev(brewer.pal(n=11, c("RdYlBu")))

# Using initial projection of range maps
pdf(file="nh_HGAPF_AOSHM_log.pdf", height=5, width=9)
par(mar=c(1,1,1,5))
proj4string(wm) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline <- readOGR('/Users/michaelharvey/Documents/research/seabirds/sampling_labwork/map/ne_50m_coastline/ne_50m_coastline.shp')
coastline2 <- spTransform(coastline, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
sPDF <- getMap()
sPDF <- spTransform(sPDF, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="")
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black", col="gray")
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="", add=TRUE)
plot(wm, col=Pal, add=TRUE)
plot(coastline2, lwd=0.5, col="black", add=TRUE)
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black")
dev.off()

# Astral AOSHM

# Get es values
es.table <- read.table("../text_files/nh_mean_median_Astral_AOSHM.txt")
es <- log(es.table$ht.mean)
names(es) <- rownames(es.table)

# Loop to get range maps
file.NAs <- vector()
polygons <- vector()
my.es <- vector()
for(i in 1:length(es)) {
	print(names(es[i]))
	my.file <- paste0('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/range_maps/range_maps_shape_AOSHM/', names(es[i]), '/', names(es[i]), '.shp')
	if(file.exists(my.file)) {
		range <- readShapeSpatial(my.file)	
		range.dissolved <- try(unionSpatialPolygons(range, rep(as.character(i), length(range@polygons))))
		if(class(range.dissolved) == "try-error") {
			print(paste("Range error for: ", names(es[i])))
		} else {
			polygons <- c(polygons, range.dissolved)			
			my.es <- c(my.es, es[i])
		}
	} else {
		print(paste("No range data for: ", names(es[i])))
		file.NAs <- c(file.NAs, names(es[i]))
	}
}
print(file.NAs)
joined = SpatialPolygons(lapply(polygons, function(x){x@polygons[[1]]}))

# Make an empty raster with 200 km resolution
r <- raster(xmn = -12000000, xmx = 16000000, ymn = -7000000, ymx = 8000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # whole world
#r <- raster(xmn = -12000000, xmx = -2000000, ymn = -3000000, ymx = 3000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # Neotropics only

# For each species, make raster for es and another for weights
es.rasters <- stack()
count.rasters <- stack()
for(i in 1:length(joined)) {
	print(i)
	cell.count <- lengths(extract(r, joined[i]))
	es.raster <- r
	count.raster <- r
	es.raster[joined[i],] <- my.es[i]
	count.raster[joined[i],] <- 1/cell.count
	es.rasters <- stack(es.rasters, es.raster)
	count.rasters <- stack(count.rasters, count.raster)
}
wm <- weighted.mean(x=es.rasters, w=count.rasters, na.rm=TRUE)
library(RColorBrewer)
Pal <- rev(brewer.pal(n=11, c("RdYlBu")))

# Using initial projection of range maps
pdf(file="nh_Astral_AOSHM_log.pdf", height=5, width=9)
par(mar=c(1,1,1,5))
proj4string(wm) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline <- readOGR('/Users/michaelharvey/Documents/research/seabirds/sampling_labwork/map/ne_50m_coastline/ne_50m_coastline.shp')
coastline2 <- spTransform(coastline, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
sPDF <- getMap()
sPDF <- spTransform(sPDF, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="")
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black", col="gray")
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="", add=TRUE)
plot(wm, col=Pal, add=TRUE)
plot(coastline2, lwd=0.5, col="black", add=TRUE)
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black")
dev.off()

# Astral AOSHM missing dropped

# Get es values
es.table <- read.table("../text_files/nh_mean_median_Astral_AOSHM_missingdropped.txt")
es <- log(es.table$ht.mean)
names(es) <- rownames(es.table)

# Loop to get range maps
file.NAs <- vector()
polygons <- vector()
my.es <- vector()
for(i in 1:length(es)) {
	print(names(es[i]))
	my.file <- paste0('/Users/michaelharvey/Documents/research/Tyranni/v2/tree_analyses/range_maps/range_maps_shape_AOSHM/', names(es[i]), '/', names(es[i]), '.shp')
	if(file.exists(my.file)) {
		range <- readShapeSpatial(my.file)	
		range.dissolved <- try(unionSpatialPolygons(range, rep(as.character(i), length(range@polygons))))
		if(class(range.dissolved) == "try-error") {
			print(paste("Range error for: ", names(es[i])))
		} else {
			polygons <- c(polygons, range.dissolved)			
			my.es <- c(my.es, es[i])
		}
	} else {
		print(paste("No range data for: ", names(es[i])))
		file.NAs <- c(file.NAs, names(es[i]))
	}
}
print(file.NAs)
joined = SpatialPolygons(lapply(polygons, function(x){x@polygons[[1]]}))

# Make an empty raster with 200 km resolution
r <- raster(xmn = -12000000, xmx = 16000000, ymn = -7000000, ymx = 8000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # whole world
#r <- raster(xmn = -12000000, xmx = -2000000, ymn = -3000000, ymx = 3000000, vals=0, crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") # Neotropics only

# For each species, make raster for es and another for weights
es.rasters <- stack()
count.rasters <- stack()
for(i in 1:length(joined)) {
	print(i)
	cell.count <- lengths(extract(r, joined[i]))
	es.raster <- r
	count.raster <- r
	es.raster[joined[i],] <- my.es[i]
	count.raster[joined[i],] <- 1/cell.count
	es.rasters <- stack(es.rasters, es.raster)
	count.rasters <- stack(count.rasters, count.raster)
}
wm <- weighted.mean(x=es.rasters, w=count.rasters, na.rm=TRUE)
library(RColorBrewer)
Pal <- rev(brewer.pal(n=11, c("RdYlBu")))

# Using initial projection of range maps
pdf(file="nh_Astral_AOSHM_missingdropped_log.pdf", height=5, width=9)
par(mar=c(1,1,1,5))
proj4string(wm) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline <- readOGR('/Users/michaelharvey/Documents/research/seabirds/sampling_labwork/map/ne_50m_coastline/ne_50m_coastline.shp')
coastline2 <- spTransform(coastline, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
sPDF <- getMap()
sPDF <- spTransform(sPDF, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="")
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black", col="gray")
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="", add=TRUE)
plot(wm, col=Pal, add=TRUE)
plot(coastline2, lwd=0.5, col="black", add=TRUE)
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black")
dev.off()

