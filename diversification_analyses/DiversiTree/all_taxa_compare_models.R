setwd('/Users/michaelharvey/Documents/research/Tyranni/v2/div_analyses/ltt_plots')
getwd()

library(ape)
library(phytools)
library(diversitree)

results <- vector()

# AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
bd.mod <- make.bd(tree)
bd.fit <- find.mle(bd.mod, c(0.1,0.5), method="optim", lower=0)
results <- rbind(results, c("AOSHM", pb.fit$lnLik, bd.fit$lnLik, bd.fit$par[2], AIC(bd.fit)-AIC(pb.fit), anova(bd.fit, pb.fit)[2,5]))

# Dups Dropped

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_dupsdropped.tre')
tree <- force.ultrametric(tree)
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
bd.mod <- make.bd(tree)
bd.fit <- find.mle(bd.mod, c(0.1,0.5), method="optim", lower=0)
results <- rbind(results, c("Dups dropped", pb.fit$lnLik, bd.fit$lnLik, bd.fit$par[2], AIC(bd.fit)-AIC(pb.fit), anova(bd.fit, pb.fit)[2,5]))

# AOSC

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_AOS_Clements.tre')
tree <- force.ultrametric(tree)
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
bd.mod <- make.bd(tree)
bd.fit <- find.mle(bd.mod, c(0.1,0.5), method="optim", lower=0)
results <- rbind(results, c("AOSC", pb.fit$lnLik, bd.fit$lnLik, bd.fit$par[2], AIC(bd.fit)-AIC(pb.fit), anova(bd.fit, pb.fit)[2,5]))

# IUCN

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_IUCN.tre')
tree <- force.ultrametric(tree)
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
bd.mod <- make.bd(tree)
bd.fit <- find.mle(bd.mod, c(0.1,0.5), method="optim", lower=0)
results <- rbind(results, c("IUCN", pb.fit$lnLik, bd.fit$lnLik, bd.fit$par[2], AIC(bd.fit)-AIC(pb.fit), anova(bd.fit, pb.fit)[2,5]))

# 1 My

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_1My_collapsed.tre')
tree <- force.ultrametric(tree)
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
bd.mod <- make.bd(tree)
bd.fit <- find.mle(bd.mod, c(0.1,0.5), method="optim", lower=0)
results <- rbind(results, c("1 My threshold", pb.fit$lnLik, bd.fit$lnLik, bd.fit$par[2], AIC(bd.fit)-AIC(pb.fit), anova(bd.fit, pb.fit)[2,5]))

# 2 My

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_2My_collapsed.tre')
tree <- force.ultrametric(tree)
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
bd.mod <- make.bd(tree)
bd.fit <- find.mle(bd.mod, c(0.1,0.5), method="optim", lower=0)
results <- rbind(results, c("2 My threshold", pb.fit$lnLik, bd.fit$lnLik, bd.fit$par[2], AIC(bd.fit)-AIC(pb.fit), anova(bd.fit, pb.fit)[2,5]))

# HGAPF AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/HGAPF_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
bd.mod <- make.bd(tree)
bd.fit <- find.mle(bd.mod, c(0.1,0.5), method="optim", lower=0)
results <- rbind(results, c("HGAPF AOSHM", pb.fit$lnLik, bd.fit$lnLik, bd.fit$par[2], AIC(bd.fit)-AIC(pb.fit), anova(bd.fit, pb.fit)[2,5]))

# Astral AOSHM

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore.tre')
tree <- force.ultrametric(tree)
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
bd.mod <- make.bd(tree)
bd.fit <- find.mle(bd.mod, c(0.1,0.5), method="optim", lower=0)
results <- rbind(results, c("Astral AOSHM", pb.fit$lnLik, bd.fit$lnLik, bd.fit$par[2], AIC(bd.fit)-AIC(pb.fit), anova(bd.fit, pb.fit)[2,5]))

# Astral AOSHM missingdropped

tree <- read.tree('~/Documents/research/Tyranni/v2/tree_building/examl/final_timetrees/T400F_Astral_AOS_HowardMoore_highmissingdropped.tre')
tree <- force.ultrametric(tree)
pb.mod <- make.yule(tree)
pb.fit <- find.mle(pb.mod, 0.1)
bd.mod <- make.bd(tree)
bd.fit <- find.mle(bd.mod, c(0.1,0.5), method="optim", lower=0)
results <- rbind(results, c("Astral AOSHM missing dropped", pb.fit$lnLik, bd.fit$lnLik, bd.fit$par[2], AIC(bd.fit)-AIC(pb.fit), anova(bd.fit, pb.fit)[2,5]))

write.table(results, "model_comparison_results2.txt")


