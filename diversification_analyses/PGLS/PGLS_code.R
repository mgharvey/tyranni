setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/pgls')
getwd()

library(ape)
library(caper)
library(phangorn)
library(phytools)
library(diversitree)

results <- vector()

# AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
all.dframe <- read.csv('pgls_data.csv', row.names=1)

# Remove Old World and/or Nearctic species

region.data <- read.table("~/Documents/research/Tyranni/v2/div_analyses/range_maps/Range_data_AOSHM_Olson_broad.txt")
areas <- cbind(region.data[,3:11])
rownames(areas) <- as.character(region.data$V1)
region.names <- c("WI", "AM", "AN", "NA", "OW", "PA", "DT", "AF", "CA")
colnames(areas) <- region.names
NewWorld.only <- areas[!colnames(areas)[apply(areas,1,which.max)] == "OW",]
Neotropics.only <- areas[!colnames(areas)[apply(areas,1,which.max)] == "NA",]
#tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(NewWorld.only)))
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(Neotropics.only)))
dframe <- all.dframe[all.dframe$species %in% tree$tip.label,]
nrow(dframe)

data <- comparative.data(data=dframe, phy=tree, names.col="species")

# PGLS model with a large set of variables (additive)
test1 <- pgls(nh.statistic ~ elev.max + elev.range + lat.max + lat.range + temp + temp.stability + min.temp.cold.q + precip + precip.stability + min.precip.dry.q, data=data)
summary(test1)




# Leave-one-out model-testing approach
test2 <- pgls(dr.statistic ~ elev.range + lat.max + lat.range + temp + temp.stability + min.temp.cold.q + precip + precip.stability + min.precip.dry.q, data=data)
test3 <- pgls(dr.statistic ~ elev.max + lat.max + lat.range + temp + temp.stability + min.temp.cold.q + precip + precip.stability + min.precip.dry.q, data=data)
test4 <- pgls(dr.statistic ~ elev.max + elev.range + lat.range + temp + temp.stability + min.temp.cold.q + precip + precip.stability + min.precip.dry.q, data=data)
test5 <- pgls(dr.statistic ~ elev.max + elev.range + lat.max + temp + temp.stability + min.temp.cold.q + precip + precip.stability + min.precip.dry.q, data=data)
test6 <- pgls(dr.statistic ~ elev.max + elev.range + lat.max + lat.range + temp.stability + min.temp.cold.q + precip + precip.stability + min.precip.dry.q, data=data)
test7 <- pgls(dr.statistic ~ elev.max + elev.range + lat.max + lat.range + temp + min.temp.cold.q + precip + precip.stability + min.precip.dry.q, data=data)
test8 <- pgls(dr.statistic ~ elev.max + elev.range + lat.max + lat.range + temp + temp.stability + precip + precip.stability + min.precip.dry.q, data=data)
test9 <- pgls(dr.statistic ~ elev.max + elev.range + lat.max + lat.range + temp + temp.stability + min.temp.cold.q + precip.stability + min.precip.dry.q, data=data)
test10 <- pgls(dr.statistic ~ elev.max + elev.range + lat.max + lat.range + temp + temp.stability + min.temp.cold.q + precip + min.precip.dry.q, data=data)
test11 <- pgls(dr.statistic ~ elev.max + elev.range + lat.max + lat.range + temp + temp.stability + min.temp.cold.q + precip + precip.stability, data=data)
AIC(test1, test2, test3, test4, test5, test6, test7, test8, test9, test10, test11)
test11$aicc-test1$aicc
test10$aicc-test1$aicc
test9$aicc-test1$aicc
test8$aicc-test1$aicc
test7$aicc-test1$aicc
test6$aicc-test1$aicc
test5$aicc-test1$aicc
test4$aicc-test1$aicc
test3$aicc-test1$aicc
test2$aicc-test1$aicc

# Make additive models for all combinations of a set of variables (based on Liam Revell code)
# First remove rows missing climate/geographic data
dframe.nomissing <- dframe[!is.na(dframe$temp),]
dframe.nomissing <- dframe.nomissing[!is.na(dframe.nomissing$lat.max),]
tree.nomissing <- drop.tip(tree, setdiff(tree$tip.label, dframe.nomissing$species))
sapply(dframe.nomissing, function(x) sum(is.na(x)))

y <- "nh.statistic"
x <- c("elev.max", "lat.max", "temp.stability", "min.temp.cold.q", "precip.stability", "min.precip.dry.q")
dredge <- function(y,x,fn=pgls.SEy,...) {
	n <- length(x)
	ii <- list()
	for(i in 1:n) ii <- c(ii, combn(x,i,simplify=FALSE))
	fits <- list()
	fits[[1]] <- fn(as.formula(paste(y,"~1")),...)
	for(i in 1:length(ii)) {
		print(paste0("Fitting model ", i, " of ", length(ii)))
		print(ii[i][[1]])
		model <- as.formula(paste(y, "~", paste(ii[[i]], collapse="+")))
		fits[[i+1]] <- fn(model,...)
	}
	fits
}

fits.bm <- dredge(y, x, data=dframe.nomissing, tree=tree.nomissing, method="ML")
aic <- sapply(fits.bm, AIC)
fits.bm[aic==min(aic)][[1]]

names <- unlist(lapply(ii, function(x) paste(x, collapse='_')))
names <- c("none", names)
names(aic) <- names
aic
which(aic==min(aic))

# Make models for all possible combinations of two variables with interaction ("*")
y <- "dr.statistic"
x <- c("elev.max", "elev.range", "lat.max", "lat.range", "temp", "temp.stability", "min.temp.cold.q", "precip", "precip.stability", "min.precip.dry.q")
dredge <- function(y,x,fn=pgls.SEy,...) {
	n <- length(x)
	ii <- list()
	for(i in 1:1) ii <- c(ii, combn(x,2,simplify=FALSE))
	fits <- list()
	fits[[1]] <- fn(as.formula(paste(y,"~1")),...)
	for(i in 1:length(ii)) {
		print(paste0("Fitting model ", i, " of ", length(ii)))
		print(ii[i][[1]])
		model <- as.formula(paste(y, "~", paste(ii[[i]], collapse="*")))
		fits[[i+1]] <- fn(model,...)
	}
	fits
}

fits.inter.bm <- dredge(y, x, data=dframe.nomissing, tree=tree.nomissing, method="ML")
aic.inter <- sapply(fits.inter.bm, AIC)
fits.inter.bm[aic.inter==min(aic.inter)][[1]]

names <- unlist(lapply(ii, function(x) paste(x, collapse='_')))
names <- c("none", names)
names(aic.inter) <- names
aic.inter
which(aic.inter==min(aic.inter))
