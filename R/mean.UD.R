###################
# average aligned UDs
mean.UD <- function(x,weights=NULL,sample=TRUE,...)
{
  n <- length(x)
  axes <- x[[1]]$axes

  if(is.null(weights))
  {
    if(x[[1]]@type=="occurrence") # time weighted by default
    { weights <- sapply(x,function(y){y$W}) }
    else
    { weights <- rep(1,length(x)) }
  }

  weights <- weights/max(weights)
  names(weights) <- names(x)
  WEIGHT <- sum(weights)

  # list of individual models
  CTMM <- lapply(x,function(y){y@CTMM})
  # population model
  CTMM <- mean.ctmm(CTMM,weights=weights,sample=sample)
  # population stationary distribution
  if(sample) { CTMM <- mean_pop(CTMM) }

  # harmonic mean bandwidth matrix
  H <- 0
  for(i in 1:n) { H <- H + weights[i] * PDsolve(x[[i]]$H) }
  H <- H/WEIGHT
  H <- PDsolve(H)

  info <- mean_info(x)
  type <- unique(sapply(x,function(y){attr(y,"type")}))
  if(length(type)>1) { stop("Distribution types ",type," differ.") }
  dV <- prod(x[[1]]$dr)

  GRID <- grid.union(x) # r,dr of grid union
  DIM <- c(length(GRID$r$x),length(GRID$r$y))
  PDF <- matrix(0,DIM[1],DIM[2]) # initialize Joint PDF

  for(i in 1:n)
  {
    SUB <- grid.intersection(list(GRID,x[[i]]))
    PDF[SUB[[1]]$x,SUB[[1]]$y] <- PDF[SUB[[1]]$x,SUB[[1]]$y] + weights[i] * x[[i]]$PDF[SUB[[2]]$x,SUB[[2]]$y]
  }
  PDF <- PDF / WEIGHT

  x <- GRID
  x$weights <- weights
  x$axes <- axes
  x$PDF <- PDF
  x$CDF <- pmf2cdf(PDF*dV)
  if(type!="occurrence") { x$DOF.area <- DOF.area(CTMM) }
  # x$H <- H
  x$H <- NULL

  x <- new.UD(x,info=info,type=type,CTMM=CTMM)

  return(x)
}
