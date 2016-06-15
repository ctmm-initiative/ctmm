fitalignedUDs <- function(data.list,
						  ctmm.list,
					      res    = 10,
					      debias = TRUE,
					      error  = 0.001,
					      silent = FALSE) {
	n.instances <- length(data.list)

	#initialize per-instance vars that need to be tracked
	HP.list   <- vector("list", n.instances)
	H.list    <- vector("list", n.instances)
	UD.list   <- vector("list", n.instances)
	n.samples.total <- 0

	#step 1: compute the optimal bandwidths for each instance
	if (!silent)
		cat("Computing optimal bandwidths\n")
	for (i.instance in 1:n.instances) {
		#alias for readability
		instance <- data.list[[i.instance]]
		n.samples <- length(instance$x)
		n.samples.total <- n.samples.total + n.samples

		#bandwidth calculations for this instance
		HP.list[[i.instance]] <- akde.bandwidth(data = instance, CTMM = ctmm.list[[i.instance]], verbose = TRUE)
		H.list[[i.instance]]  <- ctmm:::prepare.H(instance, HP.list[[i.instance]]$H)
		#reshape H
		H.list[[i.instance]]  <- array(H.list[[i.instance]], c(n.samples, 2 * 2)) #TODO: magic numbers
	}

	#step 2: compute a grid for all UDs
	if (!silent)
		cat("Computing the grid\n")
	H.all <- do.call("rbind", H.list)
	H.all <- array(H.all, c(n.samples.total, 2, 2))
	data.all <- do.call("rbind", data.list)
	dr <- vector("numeric", 2)
	dr[1] <- sqrt(min(sapply(HP.list, function(i) {
		                 	                         i$H[1,1]
		                 	                      })))
	dr[2] <- sqrt(min(sapply(HP.list, function(i) {
		                 	                         i$H[2,2]
		                 	                      })))
	dr    <- dr / res

	grid    <- ctmm:::kde.grid(data.all, H.all, alpha = error, dr = dr)
	DX.orig <- grid$DX
	DY.orig <- grid$DY

	#now align all UDs to the grid
	sample.position <- 1 #track position of next sample to read out
	for (i.instance in 1:n.instances) {
		if (!silent)
			cat(sprintf("Aligning UD for instance %d/%d\n", i.instance, n.instances))

		#repeated work, but cheap and more readable than storing from step 1
		instance <- data.list[[i.instance]]
		n.samples <- length(instance$x)

		#remember when we reshaped H? yeah undo that
		H.list[[i.instance]] <- array(H.list[[i.instance]], c(n.samples, 2, 2))

		grid$DX <- DX.orig[sample.position:(sample.position + n.samples - 1)]
		grid$DY <- DY.orig[sample.position:(sample.position + n.samples - 1)]

		if (debias) {
			debias <- HP.list[[i.instance]]$bias
		}
		UD.list[[i.instance]] <- ctmm:::kde(data.list[[i.instance]], 
			                                H     = H.list[[i.instance]], 
			                                grid  = grid,
			                                alpha = error,
			                                bias  = debias)

		sample.position <- sample.position + n.samples
	}

	return(UD.list)
}