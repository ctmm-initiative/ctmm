overlay.pdf.grid <- function(PDF.list,
	                         silent = T) {
	#plot width in inches
	#plot height must respect aspect ratio of PDF
	im.width  <- 3

	pdf.width  <- dim(PDF.list[[1]])[1]
	pdf.height <- dim(PDF.list[[1]])[2]
	im.height  <- im.width * pdf.height / pdf.width

	n.instances <- length(PDF.list)

	#draw grid
	i.grob <- 1
	n.grobs <- (n.instances - 1) * n.instances / 2
	grobs.list <- vector("list", n.grobs)
	for (i.instance in 1:(n.instances - 1)) {
		for (j.instance in (i.instance + 1):n.instances) {
			if (!silent)
				cat(sprintf("Overlapping animals %d and %d\n", i.instance, j.instance))
			PDF.overlay.png <- array(0, c(dim(PDF.list[[i.instance]]), 3))
			PDF.overlay.png[,,1] <- PDF.list[[i.instance]] / max(PDF.list[[i.instance]])
			PDF.overlay.png[,,3] <- PDF.list[[j.instance]] / max(PDF.list[[j.instance]])

			#render the image (to a grob)
			grobs.list[[i.grob]] <- rasterGrob(PDF.overlay.png, 
				                               vp     = viewport(angle=90), 
				                               interpolate = T,
				                               width  = unit(im.height, "in"), #width and height get swapped here because the 90deg
				                               height = unit(im.width,  "in"))   #viewport rotation is applied after the height/width
			i.grob <- i.grob + 1
		}
	}

	#compute the layout matrix
	#basically enumerate the upper triangular part of a matrix by row
	#R doesn't like this, so instead we enumerate the lower tri part by col and then take the transpose
	layout.matrix <- matrix(NA, n.instances - 1, n.instances - 1)
	lt.mask <- lower.tri(layout.matrix, diag = T)
	layout.matrix[lt.mask] <- 1:n.grobs
	layout.matrix <- t(layout.matrix)

	pdf("overlay.pdf.grid.pdf", 
		width  = unit((n.instances - 1) * im.width,  "in"), 
		height = unit((n.instances - 1) * im.height, "in"))
	grid.arrange(grobs = grobs.list, layout_matrix = layout.matrix)
	invisible(dev.off())
}