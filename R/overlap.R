# overlap <- function(object,...) UseMethod("overlap") #S3 generic

# forwarding function for list of a particular datatype
overlap <- function(object,level=0.95,debias=TRUE,...)
{
  check.projections(object)

  CLASS <- class(object[[1]])[1]

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


# approximate-exact Wishart DOFs
DOF.wishart <- function(CTMM)
{
  AXES <- length(CTMM$axes)
  # extract covariance matrix SIGMA and covariance of covariance matrix COV
  SIGMA <- CTMM$sigma # Wishart matrix / n
  if(CTMM$isotropic) # chi^2
  { PAR <- 'major' }
  else # approximate Wishart DOF (exact if Wishart)
  { PAR <- c('major','minor') }
  EST <- SIGMA@par[PAR]
  DOF <- CTMM[['COV']][PAR,PAR]
  if(length(DOF)==length(PAR)^2)
  { DOF <- (2/AXES) * c(EST %*% PDsolve(DOF) %*% EST) } # average multiple DOFs if not Wishart
  if(length(DOF)==0) { DOF <- 0 }
  return(DOF)
}


# n/(n-dim-1) for n>=dim+2
# leading order in 1/n for n<=dim+2
# matched to first derivative
soft.clamp <- function(n,DIM)
{
  if(n >= DIM+2) { return(n) }

  A <- -2*DIM - 5*DIM^2 - 4*DIM^3 - DIM^4
  B <- 4 + 16*DIM + 25*DIM^2 + 19*DIM^3 + 7*DIM^4 + DIM^5

  # same to leading order, with matching first derivative
  BIAS <- 1 + (DIM+1)/n + A/n^2 + B/n^3

  # n that would give the above bias with n/(n-dim-1) formula
  (DIM+1)*BIAS/(BIAS-1)
}


#####################
overlap.ctmm <- function(object,level=0.95,debias=TRUE,COV=TRUE,method="Bhattacharyya",distance=FALSE,...)
{
  CTMM1 <- object[[1]]
  CTMM2 <- object[[2]]
  DIM <- length(CTMM1$axes)

  if(method=="Bhattacharyya") { Dfunc <- BhattacharyyaD }
  else if(method=="Mahalanobis") { Dfunc <- MahalanobisD }
  else if(method=="Euclidean") { Dfunc <- EuclideanD }
  STUFF <- gauss.comp(Dfunc,object,COV=COV)
  MLE <- c(STUFF$MLE)
  VAR <- c(STUFF$COV)
  # this quantity is roughly chi-square
  DOF <- 2*MLE^2/VAR

  # approximate debiasing, correct for IID, equal covariance, REML
  ########################
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]
  COV.mu <- CTMM1$COV.mu + CTMM2$COV.mu

  if(method=="Euclidean") { sigma <- diag(1,DIM) }
  else { sigma <- (CTMM1$sigma + CTMM2$sigma)/2 }

  # trace variances
  s0 <- mean(diag(sigma))
  s1 <- mean(diag(CTMM1$sigma))
  s2 <- mean(diag(CTMM2$sigma))

  # approximate average Wishart DOFs
  # n1 <- DOF.wishart(CTMM1)
  # n2 <- DOF.wishart(CTMM2)
  # the above can be flaky
  n1 <- DOF.area(CTMM1)
  n2 <- DOF.area(CTMM2)
  # using mean variance - additive & rotationally invariant
  n0 <- 4 * s0^2 / (s1^2/n1 + s2^2/n2)
  # dim cancels out

  # hard clamp before soft clamp
  n1 <- clamp(n1,1,Inf)
  n2 <- clamp(n2,1,Inf)
  n0 <- clamp(n0,2,Inf)

  # clamp the DOF not to diverge <=DIM+1
  n0 <- soft.clamp(n0,DIM)
  n1 <- soft.clamp(n1,DIM)
  n2 <- soft.clamp(n2,DIM)

  # expectation value of log det Wishart
  ElogW <- function(s,n) { log(det(s)) + mpsigamma(n/2,dim=DIM) - DIM*log(n/2) }

  # inverse Wishart expectation value pre-factor
  BIAS <- n0/(n0-DIM-1)
  if(method=="Euclidean") { BIAS <- 0 } # don't include this term
  # BIAS <- clamped.bias(n0,DIM)
  # mean terms
  BIAS <- sum(diag((BIAS*outer(mu) + COV.mu) %*% PDsolve(sigma)))
  if(method=="Bhattacharyya")
  {
    BIAS <- BIAS/8
    # AMGM covariance terms
    BIAS <- BIAS + max( ElogW(sigma,n0)/2 - ElogW(CTMM1$sigma,n1)/4 - ElogW(CTMM2$sigma,n2)/4 , 0)
    # this is actually the expectation value?
  }

  # relative bias instead of absolute bias
  BIAS <- BIAS/MLE
  # would subtract off estimate to get absolute bias

  # error corrections
  BIAS <- as.numeric(BIAS)
  if(MLE==0) { BIAS <- 1 }
  #####################

  if(level)
  {
    if(debias) { MLE <- MLE/BIAS }

    CI <- chisq.ci(MLE,VAR=VAR,alpha=1-level)
    if(distance) { return(CI) } # return distance

    # transform from (square) distance to overlap measure
    CI <- exp(-rev(CI))
    names(CI) <- NAMES.CI

    return(CI)
  }
  else # return BD ingredients
  { return(list(MLE=MLE,VAR=VAR,DOF=DOF,BIAS=BIAS)) }
}


#####################
#overlap density function
overlap.UD <- function(object,level=0.95,debias=TRUE,...)
{
  CTMM <- list(attr(object[[1]],"CTMM"),attr(object[[2]],"CTMM"))
  type <- c(attr(object[[1]],"type"),attr(object[[2]],"type"))
  type <- type[type!="range"]
  if(length(type)) { stop(type," overlap is not generally meaningful, biologically.") }

  # check resolution and subset to overlapping grid
  object <- same.grids(object)
  # can now be null mass

  dr <- object[[1]]$dr
  dA <- prod(dr)

  OVER <- object[[1]]$PDF * object[[2]]$PDF
  if(!is.null(OVER)) { OVER <- sqrt(OVER) }

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
    CI <- chisq.ci(D,VAR=CI$VAR,alpha=1-level)

    # transform from (square) distance to overlap measure
    OVER <- exp(-rev(CI))

  }
  else
  { OVER <- c(NA,OVER,NA) }

  names(OVER) <- NAMES.CI
  return(OVER)
}
