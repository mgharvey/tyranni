nhsim <- function(phy, trait, nsim = 1000, nh, return.nh=FALSE) {
	
	require(ape)
	require(mvtnorm)
	
	if(missing(nh)) { # If inverse equal splits statistics not provided, calculate it
		rootnode <- length(tree$tip.label) + 1
		ht.mean <- numeric(length(tree$tip.label))
		ht.median <- numeric(length(tree$tip.label))
		for (i in 1:length(ht.mean)){
			node <- i
			qx <- 0
			hts <- vector()
			while (node != rootnode){
				el <- tree$edge.length[tree$edge[,2] == node]
				node <- tree$edge[,1][tree$edge[,2] == node]			
				qx <- qx + el	
				hts <- c(hts, qx)	
			}
			ht.mean[i] <- mean(as.numeric(hts))
			ht.median[i] <- median(as.numeric(hts))
		}		
		names(ht.mean) <- tree$tip.label
		names(ht.median) <- tree$tip.label
		nh <- ht.mean
	}
	
	nh <- nh[phy$tip.label] # log transform
	trait <- trait[phy$tip.label]
	
	# Pearson's correlation between log inverse equal splits statistic and trait
	res <- cor.test(nh, trait, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates
	vv <- vcv.phylo(as.phylo(phy))
	onev <- matrix(rep(1, length(trait)), nrow=length(trait), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% trait))
	rate <- as.vector((t(trait-root) %*% solve(vv) %*% (trait-root))/length(trait))
	
	# Brownian simulations 
	sims <- t(rmvnorm(nsim, sigma=rate*vv))
	rownames(sims) <- rownames(vv)
		
	# Pearson's correlations of simulated datasets
	sim.r <- sapply(1:nsim, function(x) cor.test(nh[as.vector(rownames(sims))], sims[,x], method="pearson")$estimate)
	
	# Calculate the two-tailed p value
	corr <- res$estimate
	upper <- (length(sim.r[sim.r >= corr])+1)/(nsim+1)
	lower <- (length(sim.r[sim.r <= corr])+1)/(nsim+1)
	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed

	if(missing(return.nh)) { # output just rho and p value
		result <- as.vector(c(corr, pval))
		names(result) <- c("rho", "P Value")
		return(result)
	} else { # output rho, p value, and list of nh values
		result <- as.vector(c(corr, pval, list(nh)))
		names(result) <- c("rho", "P Value", "nh")
		return(result)		
	}

}
