# overlap <- function(object,...) UseMethod("overlap") #S3 generic

# forwarding function for list of a particular datatype
overlap <- function(object,level=0.95,debias=TRUE,...)
{
  CLASS <- class(object[[1]])

  # check for consistent projections
  PROJ <- projection(object)
  if(length(PROJ)==0) { stop("Missing projection.") }
  if(length(PROJ)>1) { stop("Inconsistent projections.") }

  if(CLASS=="ctmm")
  { OverlapFun <- overlap.ctmm }
  # else if(CLASS=="telemetry")
  # {
  #   # Generate aligned UDs
  #   object <- akde(object,...)
  #   OverlapFun <- overlap.UD
  # }
  else if(CLASS=="UD")
  { OverlapFun <- overlap.UD }
  else { stop(CLASS," object class not supported by overlap.") }

  n <- length(object)
  OVER <- array(0,c(n,n,3))
  # tabulate overlaps
  for(i in 1:n)
  {
    for(j in (i+1):n)
    { if(j<=n) { OVER[i,j,] <- OverlapFun(object[c(i,j)],level=level,debias=debias,...) } }
  }

  # symmetrize matrix
  OVER <- OVER + aperm(OVER,c(2,1,3))

  # fix diagonals
  diag(OVER[,,1]) <- diag(OVER[,,2]) <- diag(OVER[,,3]) <- 1

  dimnames(OVER) <- list(names(object),names(object),NAMES.CI)

  return(OVER)
  # utils::getS3method("overlap",CLASS)(object,...)
}


#####################
overlap.ctmm <- function(object,level=0.95,debias=TRUE,COV=TRUE,...)
{
  CTMM1 <- object[[1]]
  CTMM2 <- object[[2]]

  # ML parameters
  par <- NULL
  parscale <- NULL
  lower <- NULL
  upper <- NULL

  # first mean
  par <- c(par, CTMM1$mu[1,] )
  parscale <- c(parscale, sqrt(diag(CTMM1$sigma)) )
  lower <- c(lower, c(-Inf,-Inf) )
  upper <- c(upper, c(Inf,Inf) )

  # first covariance
  par <- c(par,CTMM1$sigma@par)
  parscale <- c(parscale, c(CTMM1$sigma@par[1:2],pi/2) )
  lower <- c(lower, c(0,0,-Inf) )
  upper <- c(upper, c(Inf,Inf,Inf) )

  # second mean
  par <- c(par,CTMM2$mu[1,])
  parscale <- c(parscale, sqrt(diag(CTMM2$sigma)) )
  lower <- c(lower, c(-Inf,-Inf) )
  upper <- c(upper, c(Inf,Inf) )

  # second covariance
  par <- c(par,CTMM2$sigma@par)
  parscale <- c(parscale, c(CTMM2$sigma@par[1:2],pi/2) )
  lower <- c(lower, c(0,0,-Inf) )
  upper <- c(upper, c(Inf,Inf,Inf) )

  D <- function(p)
  {
    CTMM1$mu[1,] <- p[1:2]
    CTMM1$sigma <- covm(p[3:5],isotropic=CTMM1$isotropic)

    CTMM2$mu[1,] <- p[6:7]
    CTMM2$sigma <- covm(p[8:10],isotropic=CTMM2$isotropic)

    return(BhattacharyyaD(CTMM1,CTMM2))
  }

  VAR <- 0
  # propagate uncertainty - slightly wasteful if isotropic
  if(COV)
  {
    # grad <- numDeriv::grad(D,par)
    grad <- genD(par,D,lower=lower,upper=upper,parscale=parscale,mc.cores=1,order=1,...)$gradient

    VAR <- VAR + abs(grad[1:2] %*% CTMM1$COV.mu %*% grad[1:2])
    if(CTMM1$isotropic)
    { VAR <- VAR + abs(grad[3] * CTMM1$COV[1,1] * grad[3]) }
    else
    { VAR <- VAR + abs(grad[3:5] %*% CTMM1$COV[1:3,1:3] %*% grad[3:5]) }

    VAR <- VAR + abs(grad[6:7] %*% CTMM2$COV.mu %*% grad[6:7])
    if(CTMM2$isotropic)
    { VAR <- VAR + abs(grad[8] * CTMM2$COV[1,1] * grad[8]) }
    else
    { VAR <- VAR + abs(grad[8:10] %*% CTMM2$COV[1:3,1:3] %*% grad[8:10]) }

    VAR <- as.numeric(VAR)
  }

  # this quantity is roughly chi-square
  MLE <- BhattacharyyaD(CTMM1,CTMM2)
  DOF <- 2*MLE^2/VAR


  # approximate debiasing, correct for IID, equal covariance, REML
  ########################
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]
  COV.mu <- CTMM1$COV.mu + CTMM2$COV.mu

  sigma <- (CTMM1$sigma + CTMM2$sigma)/2
  DIM <- nrow(sigma)

  s0 <- mean(diag(sigma)^2)
  s1 <- mean(diag(CTMM1$sigma)^2)
  s2 <- mean(diag(CTMM2$sigma)^2)

  n1 <- DOF.area(CTMM1)
  n2 <- DOF.area(CTMM2)

  # approximate average Wishart DOF
  # using mean variance - additive & rotationally invariant
  n0 <- 4 * s0 / (s1/n1 + s2/n2)

  # clamp the DOF not to diverge
  n0 <- clamp(n0,DIM+2,Inf)
  n1 <- clamp(n1,DIM+2,Inf)
  n2 <- clamp(n2,DIM+2,Inf)

  # expectation value of log det Wishart
  ElogW <- function(s,n) { log(det(s)) + mpsigamma(n/2,dim=DIM) - DIM*log(n/2) }

  # inverse Wishart expectation value pre-factor
  BIAS <- n0/(n0-DIM-1)
  # mean terms
  BIAS <- sum(diag((BIAS*outer(mu) + COV.mu) %*% PDsolve(sigma)))/8
  # AMGM covariance terms
  BIAS <- BIAS + max(ElogW(sigma,n0)/2 - ElogW(CTMM1$sigma,n1)/4 - ElogW(CTMM2$sigma,n2)/4 , 0)
  # this is actually the expectation value?

  # relative bias instead of absolute bias
  BIAS <- BIAS/MLE
  # would subtract off estimte to get absolute bias

  # error corrections
  BIAS <- as.numeric(BIAS)
  if(MLE==0) { BIAS <- 1 }
  #####################

  if(level)
  {
    if(debias) { MLE <- MLE/BIAS }

    CI <- chisq.ci(MLE,COV=VAR,alpha=1-level)

    # transform from (square) distance to overlap measure
    CI <- exp(-rev(CI))
    names(CI) <- NAMES.CI

    return(CI)
  }
  else # return BD ingredients
  { return(list(MLE=MLE,VAR=VAR,DOF=DOF,BIAS=BIAS)) }
}

####################
# square distance between stationary Gaussian distributions
BhattacharyyaD <- function(CTMM1,CTMM2)
{
  sigma <- (CTMM1$sigma + CTMM2$sigma)/2
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]

  D <- as.numeric(mu %*% PDsolve(sigma) %*% mu)/8 + log(det(sigma)/sqrt(det(CTMM1$sigma)*det(CTMM2$sigma)))/2

  return(D)
}

#####################
#overlap density function
overlap.UD <- function(object,level=0.95,debias=TRUE,...)
{
  CTMM <- list(attr(object[[1]],"CTMM"),attr(object[[2]],"CTMM"))
  type <- c(attr(object[[1]],"type"),attr(object[[2]],"type"))
  type <- type[type!="range"]
  if(length(type)) { stop(type," overlap is not generally meaningful, biologically.") }

  dr <- object[[1]]$dr
  dA <- prod(dr)

  OVER <- sqrt(object[[1]]$PDF * object[[2]]$PDF)

  # overlap point estimate
  OVER <- sum(OVER)*dA

  if(!is.null(CTMM))
  {
    # calculate Gaussian overlap distance^2 variance, bias, etc.
    CI <- overlap.ctmm(CTMM,level=FALSE)

    # Bhattacharyya distances
    D <- -log(OVER)

    # relative debias
    if(debias){ D <- D/CI$BIAS }

    # calculate new distance^2 with KDE point estimate
    CI <- chisq.ci(D,COV=CI$VAR,alpha=1-level)

    # transform from (square) distance to overlap measure
    OVER <- exp(-rev(CI))

  }
  else
  { OVER <- c(NA,OVER,NA) }

  names(OVER) <- NAMES.CI
  return(OVER)
}


###################
# average aligned UDs
mean.UD <- function(x,...)
{
  info <- mean.info(x)
  dV <- prod(x[[1]]$dr)
  n <- length(x)
  N <- sum(sapply(x,function(y){ y$DOF.area }))

  M <- lapply(x,function(ud){ ud$PDF })
  M <- Reduce('+',M)
  M <- M/n # now the average PDF

  x <- x[[1]]
  x$PDF <- M
  attr(x,"info") <- info
  x$DOF.area <- 0*x$DOF.area + N

  x$H <- NULL
  x$h <- NULL
  x$DOF.H <- NA
  x$bias <- NULL
  x$weights <- NULL
  x$MISE <- NULL

  M <- M*dV # now the average cell PMF
  x$CDF <- pmf2cdf(M)

  # to be continued....
  attr(x,"CTMM") <- ctmm()

  return(x)
}
