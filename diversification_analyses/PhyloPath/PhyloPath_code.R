setwd("tyranni/diversification_analyses/PhyloPath")
getwd()

library(ape)
library(phytools)
library(phylopath)

# Get data

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
all.dframe <- read.csv('pgls_data2.csv', row.names=1)

# Remove Old World and/or Nearctic species

region.data <- read.table("tyranni/other_data/Range_data_AOSHM_Olson_broad.txt")
areas <- cbind(region.data[,3:11])
rownames(areas) <- as.character(region.data$V1)
region.names <- c("WI", "AM", "AN", "NA", "OW", "PA", "DT", "AF", "CA")
colnames(areas) <- region.names
NewWorld.only <- areas[!colnames(areas)[apply(areas,1,which.max)] == "OW",]
#Neotropics.only <- areas[!colnames(areas)[apply(areas,1,which.max)] == "NA",]
tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(NewWorld.only)))
#tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(Neotropics.only)))
dframe <- all.dframe[all.dframe$species %in% tree$tip.label,]
nrow(dframe)

dframe[1,]

# Set up tests
m <- define_model_set(
	null = c(),
	direct = c(dr.statistic~elev.max + temp + temp.stability + temp.seasonality + precip + precip.stability + precip.seasonality),
	indirect = c(dr.statistic~div, div~elev.max + temp + temp.stability + temp.seasonality + precip + precip.stability + precip.seasonality),
	both = c(dr.statistic~elev.max + temp + temp.stability + temp.seasonality + precip + precip.stability + precip.seasonality, dr.statistic~div, div~elev.max + temp + temp.stability + temp.seasonality + precip + precip.stability + precip.seasonality)
)

# Set up tests only strong variables
# temp.mean, temp.stability, temp.seasonality, precip.mean, precip.seasonality
m <- define_model_set(
	null = c(),
	direct = c(dr.statistic~elev.max + temp + temp.stability + temp.seasonality + precip + precip.seasonality),
	indirect = c(dr.statistic~div, div~elev.max + temp + temp.stability + temp.seasonality + precip + precip.seasonality),
	both = c(dr.statistic~elev.max + temp + temp.stability + temp.seasonality + precip + precip.seasonality, dr.statistic~div, div~elev.max + temp + temp.stability + temp.seasonality + precip + precip.stability + precip.seasonality)
)


plot_model_set(m, text_size=1.5, edge_width=0.5, arrow = grid::arrow(type = "closed", 10, grid::unit(15, "points")))
par(mfrow=c(1,3))
pdf(file = "DAGs.pdf", height=3, width=4)
plot(m$direct, text_size=1.5, edge_width=0.5, arrow = grid::arrow(type = "closed", 8, grid::unit(6, "points")))
plot(m$indirect, text_size=1.5, edge_width=0.5, arrow = grid::arrow(type = "closed", 8, grid::unit(6, "points")), main="Indirect Model")
plot(m$both, text_size=1.5, edge_width=0.5, arrow = grid::arrow(type = "closed", 8, grid::unit(6, "points")), main="Both Model")
dev.off()

p <- phylo_path(m, dframe, tree)
p
s <- summary(p)
s
pdf(file = "PPA_results.pdf", width=6, height=6)
plot(s)
dev.off()

b <- best(p)
plot(b, text_size=3.5)
coef_plot(b, error_bar = "se", order_by="strength", to="dr.statistic")


