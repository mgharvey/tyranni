setwd("tyranni/diversification_analyses/TESS")
getwd()

library(TESS)
library(phytools)

tree <- read.tree('tyranni/species_trees/final_timetrees/T400F_AOS_HowardMoore.tre')
#tree <- read.tree('T400F_AOS_HowardMoore_Neotropics.tre')
tree <- force.ultrametric(tree)

# Set priors
samplingFraction <- 1
numExpectedRateChanges <- 2 # 5
numExpectedMassExtinctions <- 2
pMassExtinctionPriorShape2 <- 100
expectedSurvivalProbability <- 0.5
pMassExtinctionPriorShape1 <- - pMassExtinctionPriorShape2 * expectedSurvivalProbability / (expectedSurvivalProbability - 1)

times <- as.numeric(branching.times(tree))

set.seed(12345)
tess.analysis(tree,
	empiricalHyperPriors = TRUE,
	samplingProbability = samplingFraction,
	numExpectedRateChanges = numExpectedRateChanges,
	numExpectedMassExtinctions = numExpectedMassExtinctions,
	pMassExtinctionPriorShape1 = pMassExtinctionPriorShape1,
	pMassExtinctionPriorShape2 = pMassExtinctionPriorShape2,
	dir = "comet_hyperpriors")

output <- tess.process.output("comet_hyperpriors",
	numExpectedRateChanges = numExpectedRateChanges,
	numExpectedMassExtinctions = numExpectedMassExtinctions)

pdf("TESS_speciation_plot.pdf", width=6, height=6)

par(mfrow=c(1,1))
tess.plot.output(output, fig.types = c("speciation rates"), ylim=c(0,0.3))

dev.off()
	
pdf("all_TESS_plots.pdf", width=6, height=7)

layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)

tess.plot.output(output,
	fig.types = c("speciation rates",
		"speciation shift times",
		"extinction rates",
		"extinction shift times",
		"mass extinction Bayes factors",
		"mass extinction times"),
	las=2)
	
dev.off()