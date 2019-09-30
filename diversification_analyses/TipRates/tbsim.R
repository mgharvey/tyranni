tbsim <- function(phy, trait, nsim = 1000, tb, return.tb=FALSE) {
	
	require(ape)
	require(mvtnorm)
	
	if(missing(tb)) { # If inverse equal splits statistics not provided, calculate it
		# Calculate terminal edge lengths
		n <- length(tree$tip.label)
		# based on post on Liam Revell's blog:
		tb <- setNames(tree$edge.length[sapply(1:n, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)	
	}
	
	tb <- tb[phy$tip.label] # log transform
	trait <- trait[phy$tip.label]
	
	# Pearson's correlation between log inverse equal splits statistic and trait
	res <- cor.test(tb, trait, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates
	vv <- vcv.phylo(as.phylo(phy))
	onev <- matrix(rep(1, length(trait)), nrow=length(trait), ncol=1)
	root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% trait))
	rate <- as.vector((t(trait-root) %*% solve(vv) %*% (trait-root))/length(trait))
	
	# Brownian simulations 
	sims <- t(rmvnorm(nsim, sigma=rate*vv))
	rownames(sims) <- rownames(vv)
		
	# Pearson's correlations of simulated datasets
	sim.r <- sapply(1:nsim, function(x) cor.test(tb[as.vector(rownames(sims))], sims[,x], method="pearson")$estimate)
	
	# Calculate the two-tailed p value
	corr <- res$estimate
	upper <- (length(sim.r[sim.r >= corr])+1)/(nsim+1)
	lower <- (length(sim.r[sim.r <= corr])+1)/(nsim+1)
	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed

	if(missing(return.tb)) { # output just rho and p value
		result <- as.vector(c(corr, pval))
		names(result) <- c("rho", "P Value")
		return(result)
	} else { # output rho, p value, and list of tb values
		result <- as.vector(c(corr, pval, list(tb)))
		names(result) <- c("rho", "P Value", "tb")
		return(result)		
	}

}
