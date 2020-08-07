setwd("tyranni/diversification_analyses/QuaSSE")
getwd()

library(ape)
library(phytools)
library(diversitree)

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
trait.table <- read.table("tyranni/other_data/Range_overlap_data.txt", sep=" ", header=TRUE, row.names=1)
trait <- log1p(rowSums(trait.table, na.rm=TRUE))

# Need to convert trait names from species names to tipnamecodes
name.map.file <- read.csv('../../Species_name_map_uids.csv')
name.map <- name.map.file$tipnamecodes
names(name.map) <- name.map.file$aos.howardmoore.species
name.map.used <- name.map[name.map %in% tree$tip.label]
name.map.used.unique <- name.map.used[!duplicated(name.map.used)]
newnames <- as.character(name.map.used.unique[names(trait)])
names(trait) <- newnames
trait <- trait[tree$tip.label]
trait <- trait[!is.na(trait)]
#trait <- trait[1:30]
subtree <- drop.tip(tree, setdiff(tree$tip.label, names(trait)))

# Remove Old World and/or Nearctic species
region.data <- read.table("tyranni/other_data/Range_data_AOSHM_Olson_broad.txt")
areas <- cbind(region.data[,3:11])
rownames(areas) <- as.character(region.data$V1)
region.names <- c("WI", "AM", "AN", "NA", "OW", "PA", "DT", "AF", "CA")
colnames(areas) <- region.names
NewWorld <- areas[!colnames(areas)[apply(areas,1,which.max)] == "OW",]
subtree <- drop.tip(subtree, setdiff(subtree$tip.label, rownames(NewWorld)))
trait <- trait[subtree$tip.label]
save(list=ls(), file="log.rda")


# QuaSSE

# constant

p.constant <- starting.point.quasse(subtree, trait)
xr <- range(trait) + c(-1,1) * 20 * p.constant["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])
nodrift <- function(f) constrain(f, drift ~ 0)
lik.constant <- make.quasse(subtree, trait, 1, constant.x, constant.x)
fit.constant <- find.mle(nodrift(lik.constant), p.constant, verbose=7)
print("constant complete")
print(fit.constant$par.full)
print(coef(fit.constant))
print(logLik(fit.constant))
cat("constant\n", file="QuaSSE_NewWorld_parameter_estimates_LOG.txt")
write.table(cbind(names(fit.constant$par.full), as.character(fit.constant$par.full)), file="QuaSSE_NewWorld_parameter_estimates_LOG.txt", append=TRUE, row.names=FALSE, col.names=FALSE)
save(list=ls(), file="log.rda")


# speciation rate varies
p.s <- c(fit.constant$par[1], l.m=0, fit.constant$par[2:3])
lik.s <- make.quasse(subtree, trait, 1, linear.x, constant.x)
fit.s <- find.mle(nodrift(lik.s), p.s, verbose=7)
print("speciation variable complete")
print(fit.s$par.full)
print(coef(fit.s))
print(logLik(fit.s))
cat("variable.speciation\n", file="QuaSSE_NewWorld_parameter_estimates_LOG.txt", append=TRUE)
write.table(cbind(names(fit.s$par.full), as.character(fit.s$par.full)), file="QuaSSE_NewWorld_parameter_estimates_LOG.txt", append=TRUE, row.names=FALSE, col.names=FALSE)
save(list=ls(), file="log.rda")


# extinction rate varies
p.e <- c(fit.constant$par[1:2], m.m=0, fit.constant$par[3])
lik.e <- make.quasse(subtree, trait, 1, constant.x, linear.x)
fit.e <- find.mle(nodrift(lik.e), p.e, verbose=7)
print("extinction variable complete")
print(fit.e$par.full)
print(coef(fit.e))
print(logLik(fit.e))
cat("variable.extinction\n", file="QuaSSE_NewWorld_parameter_estimates_LOG.txt", append=TRUE)
write.table(cbind(names(fit.e$par.full), as.character(fit.e$par.full)), file="QuaSSE_NewWorld_parameter_estimates_LOG.txt", append=TRUE, row.names=FALSE, col.names=FALSE)
save(list=ls(), file="log.rda")


# speciation and extinction rates vary
p.se <- c(fit.constant$par[1], l.m=0, fit.constant$par[2], m.m=0, fit.constant$par[3])
lik.se <- make.quasse(subtree, trait, 1, linear.x, linear.x)
fit.se <- find.mle(nodrift(lik.se), p.se, verbose=7)
print("speciation and extinction variable complete")
print(fit.se$par.full)
print(coef(fit.se))
print(logLik(fit.se))
cat("variable.speciation.extinction\n", file="QuaSSE_NewWorld_parameter_estimates_LOG.txt", append=TRUE)
write.table(cbind(names(fit.se$par.full), as.character(fit.se$par.full)), file="QuaSSE_NewWorld_parameter_estimates_LOG.txt", append=TRUE, row.names=FALSE, col.names=FALSE)
save(list=ls(), file="log.rda")
# comparisons
res.all <- anova(fit.constant, spec = fit.s, ext = fit.e, spec.ext = fit.se)
print(res.all)
res.s.c <- anova(fit.s, constant=fit.constant)
print(res.s.c)
res.e.c <- anova(fit.e, constant=fit.constant)
print(res.e.c)
res.se.c <- anova(fit.se, constant=fit.constant)
print(res.se.c)

cat("all\n", file="QuaSSE_NewWorld_model_comparisons.txt")
capture.output(res.all, file="QuaSSE_NewWorld_model_comparisons.txt", append=TRUE)
cat("\nspeciationVSconstant\n", file="QuaSSE_NewWorld_model_comparisons.txt", append=TRUE)
capture.output(res.s.c, file="QuaSSE_NewWorld_model_comparisons.txt", append=TRUE)
cat("\nextinctionVSconstant\n", file="QuaSSE_NewWorld_model_comparisons.txt", append=TRUE)
capture.output(res.e.c, file="QuaSSE_NewWorld_model_comparisons.txt", append=TRUE)
cat("\nspeciationextinctionVSconstant\n", file="QuaSSE_NewWorld_model_comparisons.txt", append=TRUE)
capture.output(res.se.c, file="QuaSSE_NewWorld_model_comparisons.txt", append=TRUE)
