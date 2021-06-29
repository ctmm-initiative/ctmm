# TODO link function adjustment
# TODO zero before sum of log.det, quadratic
# TODO multivariate case
ctpm.loglike <- function(data,CTMM,link=identity,REML=FALSE,profile=TRUE,zero=0,verbose=FALSE)
{
  lag <- attr(data,"lag")
  trait <- c( get.telemetry(data,axes=CTMM$axes) )
  n <- length(trait)
  K <- length(CTMM$tau)

  range <- CTMM$range
  if(!range)
  { REML <- TRUE }
  else # catch bad parameters
  {
    if(length(CTMM$tau) && CTMM$tau[1]==Inf)
    { return(-Inf) }
  }

  VAR.MULT <- n/(n-REML)

  # default identity link function for traits
  if(class(link)=="function")
  {
    if(identical(link,identity))
    { link <- list(fn=identity,grad=function(x){rep(1,length(x))}) }
    else if(identical(link,log))
    { link <- list(fn=log,grad=function(x){1/x}) }
  }
  # full transform is used so that link function can be fitted / selected
  grad <- link$grad(trait)
  trait <- link$fn(trait)

  TEST <- CTMM
  TEST$sigma <- covm(1,axes=CTMM$axes,CTMM$isotropic) # unit variance for SVF/ACF below
  COR <- svf.func(TEST)
  if(range) # stationary likelihood
  { COR <- COR$ACF }
  else # likelihood that avoids explicit conditioning
  { COR <- COR$SVF }
  COR <- Vectorize(COR)
  COR <- COR(lag)
  dim(COR) <- dim(lag)
  if(!range) { COR <- -COR } # ACF = COV - SVF

  EIGEN <- eigen(COR,symmetric=TRUE)
  # SVF case needs to drop the large negative eigenvalue, which pairs to infinite variance term dropped
  if(!range)
  {
    EIGEN$values <- EIGEN$values[-n]
    EIGEN$vectors <- EIGEN$vectors[,-n]
  }
  N <- length(EIGEN$values)

  if(tail(EIGEN$values,1)<=0) { return(-Inf) }

  if(range) # profile mean
  {
    # this is O(n^3)
    # W <- colSums(solve(COR))
    # this is obtuse but O(n^2)
    W <- c( (colSums(EIGEN$vectors)/EIGEN$values) %*% t(EIGEN$vectors) )
    COV.mu <- 1/sum(W)
    W <- W * COV.mu
    mu <- c(W %*% trait)
  }
  else # ignored projection (from infinite variance limit)
  {
    mu <- mean(trait)
    COV.mu <- Inf
  }

  # detrend mean
  trait <- trait - mu

  # slow calculation
  # Q <- c(trait %*% iCOR %*% trait)
  # obtuse but fast
  Q <- (trait %*% EIGEN$vectors)^2/EIGEN$values
  # not yet summed
  log.det <- log(EIGEN$values)
  # not yet summed
  log.det <- mean(log.det) # per N

  # profile the variance / diffusion rate
  if(profile)
  {
    sigma <- sum(Q) / (n-REML)
    Q <- 0 # - MLE per N
  }
  else
  {
    sigma <- attr(CTMM$sigma,"par")[1]
    Q <- mean(Q)/sigma - 1/VAR.MULT # - MLE per N
  }
  COV.mu <- sigma * COV.mu

  # unchanging constants
  LL.CONST <- - 1/2/VAR.MULT - 1/2*log(2*pi)

  loglike <- -Q/2 - log(sigma)/2 - log.det/2 # per N
  loglike <- N*(loglike + (LL.CONST-zero/N))

  if(range && REML) { loglike <- loglike + log(COV.mu)/2 }

  if(verbose)
  {
    # assign variables
    if(profile)
    {
      sigma <- covm(sigma,isotropic=TRUE,axes=CTMM$axes)
      CTMM$sigma <- sigma
    }

    CTMM <- ctmm.repair(CTMM,K=K)

    # mu <- drift@shift(mu,mu.center) # translate back to origin from center
    CTMM$mu <- mu
    CTMM$COV.mu <- COV.mu

    CTMM$loglike <- loglike
    #attr(CTMM,"info") <- attr(data,"info")

    CTMM <- ctmm.ctmm(CTMM)
    return(CTMM)
  }
  else
  {
    return(loglike)
  }
}
