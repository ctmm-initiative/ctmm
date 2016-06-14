OI.response <- function(UD.list,
	                    trunc.vec = seq(0.05, 0.95, 0.05),
	                    OI.method = "bc",
	                    do.plot   = TRUE) {
	n.instances <- length(UD.list)
	n.trunc     <- length(trunc.vec)

	OI <- array(data = 0, dim = c(n.instances, n.instances, n.trunc, n.trunc))
	for (i.instance in 1:n.instances) {
		for (i.trunc in 1:n.trunc) {
			#truncate the PDF of the ith dist'n
			UD.list[[i.instance]]$PDF.trunc <- UD.list[[i.instance]]$PDF * 
			                                   (UD.list[[i.instance]]$CDF < trunc.vec[i.trunc])
			for (j.instance in 1:n.instances) {
				for (j.trunc in 1:n.trunc) {
					#truncate the PDF of the jth dist'n
					UD.list[[j.instance]]$PDF.trunc <- UD.list[[j.instance]]$PDF * 
			                                           (UD.list[[j.instance]]$CDF < trunc.vec[j.trunc])

			        OI[i.instance, j.instance, i.trunc, j.trunc] <- 
			           OverlapDiscrete(UD.list[[i.instance]]$PDF.trunc,
			                           UD.list[[j.instance]]$PDF.trunc,
			                           method = OI.method,
			                           dA = UD.list[[i.instance]]$dA) #all UDs should have the same dA
				}
			}
		}
	}

	if (do.plot) {
		#plot equal-truncation OI response curves
		pdf("OI_response.pdf", width = 14, height = 14)
		par(mfrow = c(n.instances - 1, n.instances - 1))
		for (i.instance in 1:(n.instances - 1)) {
			for (j.instance in 1:(n.instances - 1)) {
				if (j.instance < i.instance) {
					plot.new()
				} else {
					OI.vec = diag(OI[i.instance, j.instance + 1, , ])

					if (!any(is.finite(OI.vec))) {
						plot.new()
					} else {
						plot(trunc.vec, OI.vec,
							 xlab = "trunc", ylab = OI.method,
							 main = sprintf("Overlap of Animal #%d and #%d", i.instance, j.instance + 1),
							 ylim = c(min(OI.vec[is.finite(OI.vec)]), max(OI.vec[is.finite(OI.vec)])))
					}
				}
			}
		}
		dev.off()
	}

	return(OI)
}