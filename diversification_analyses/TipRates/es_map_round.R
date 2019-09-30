setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/summary_stats/es/maps')
getwd()

library(raster)
library(maptools)
library(maps)
library(rgdal)
library(rworldmap)
library(sp)
library(ggplot2)

# Get es values
es.table <- read.table("../text_files/es_AOSHM.txt")
es <- es.table$x
names(es) <- rownames(es.table)

# Loop to get range maps
file.NAs <- vector()
polygons <- vector()
my.es <- vector()
for(i in 1:length(es)) {
	print(names(es[i]))
	my.file <- paste0('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/range_maps/range_maps_shape_AOSHM/', names(es[i]), '/', names(es[i]), '.shp')
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

saveRDS(wm, "wm_raster.rds")

# Using initial projection of range maps
#pdf(file="es_AOSHM.pdf", height=5, width=9)
par(mar=c(1,1,1,5))
proj4string(wm) <- "+proj=ortho +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
coastline <- readOGR('/Users/michaelharvey/Documents/research/seabirds/sampling_labwork/map/ne_50m_coastline/ne_50m_coastline.shp')
coastline2 <- spTransform(coastline, CRS("+proj=ortho +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

## start plot & extract coordinates from orthographic map
o <- c(0,-60,0) # orientation

## overlay world map
xy <- map("world", proj="orthographic", orientation=o, fill=TRUE, col="white", add=TRUE)
xy <- na.omit(data.frame(do.call(cbind, xy[c("x","y")])))
## draw a circle around the points for coloring the ocean 
polygon(max(xy$x)*sin(seq(0,2*pi,length.out=100)),max(xy$y)*cos(seq(0,2*pi,length.out=100)), 
        col="gray", border="black", lwd=2)
map("world", proj="orthographic", orientation=o, fill=TRUE, col="white", add=TRUE)
map("world", proj="orthographic", orientation=o, fill=TRUE, col="white", lwd=0.000001, add=TRUE)
map(wm, col=Pal, add=TRUE)

sPDF <- getMap()
sPDF <- spTransform(sPDF, CRS("+proj=ortho +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="")
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black", col="gray")
mapCountryData(sPDF, nameColumnToPlot="continent", colourPalette=c("white", "white"), addLegend=FALSE, borderCol="white", mapTitle="", add=TRUE)
plot(wm, col=Pal, add=TRUE)
plot(coastline2, lwd=0.5, col="black", add=TRUE)
polygon((bbox(coastline2)[1,2]*1.024)*sin(seq(0,2*pi,length.out=100)),(bbox(coastline2)[2,2]*1.024)*cos(seq(0,2*pi,length.out=100)), lwd=1, border="black")
#dev.off()


# World map
worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id
world.df <- world.points[,c("long","lat","group", "region")]

wm_df <- as.data.frame(wm, xy=TRUE)
str(wm_df)

map <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group)) +
  geom_tile(data=wm_df, aes(x=x, y=y, fill=layer), alpha=0.8) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  coord_map("ortho", orientation=c(0, 280, 0))
map








