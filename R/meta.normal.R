# TODO make sure this code works when there is no natural variance !!!!!!!!!!!!

# population-level parameter estimates for normally distributed parameters and parameter uncertainties
# VARS boolean denotes whether or not there is natural variance
# MEANS boolean denotes whether or not there is a non-zero mean
meta.normal <- function(MU,SIGMA,MEANS=TRUE,VARS=TRUE,isotropic=FALSE,debias=TRUE,precision=1/2)
{
  if(length(dim(MU))<2)
  {
    NAMES <- "x"
    N <- length(MU)
    DIM <- 1
    MU <- array(MU,c(N,1))
    SIGMA <- array(SIGMA,c(N,1,1))
  }
  else
  {
    NAMES <- colnames(MU)
    N <- dim(MU)[1]
    DIM <- dim(MU)[2]
  }

  # do we give variance to each dimension
  VARS <- array(VARS,DIM)
  ZEROV <- !VARS

  # do we give variance to each dimension
  MEANS <- array(MEANS,DIM)
  ZEROM <- !MEANS

  # observations
  OBS <- array(TRUE,c(N,DIM))
  for(i in 1:N) { OBS[i,] <- diag(cbind(SIGMA[i,,]))<Inf }

  tol <- .Machine$double.eps^precision
  REML <- debias

  # initial guesses
  mu <- colMeans(MU)
  if(any(ZEROM)) { mu[ZEROM] <- 0 }

  sigma <- 0
  for(i in 1:N) { sigma <- sigma + outer(MU[i,]-mu) }
  sigma <- sigma/(N-REML)
  if(isotropic) { sigma <- diag( mean(diag(sigma)), DIM ) }
  if(any(ZEROV))
  {
    sigma[ZEROV,] <- 0
    sigma[,ZEROV] <- 0
  }

  ERROR <- Inf
  loglike <- loglike.old <- -Inf
  while(ERROR>tol && loglike>=loglike.old)
  {
    loglike.old <- loglike
    sigma.old <- sigma

    # estimate mu exactly
    P <- array(0,c(N,DIM,DIM))
    mu <- P.mu <- 0
    for(i in 1:N)
    {
      P[i,,] <- PDsolve(sigma + SIGMA[i,,])
      P.mu <- P.mu + P[i,,]
      mu <- mu + c(P[i,,] %*% MU[i,])
    }
    COV.mu <- PDsolve(P.mu)
    mu <- c(COV.mu %*% mu)
    # we are solving the unconstrained mu above and then projecting back to the constrained mu below
    if(any(ZEROM)) { mu[ZEROM] <- 0 }

    loglike <- REML/2*log(abs(det(COV.mu))) - DIM*(N-REML)/2*log(2*pi)
    # gradient with respect to sigma
    RHS <- 0
    LHS <- P.mu
    for(i in 1:N)
    {
      D <- mu - MU[i,]
      RHS <- RHS + (P[i,,] %*% outer(D) %*% P[i,,])
      if(debias) { LHS <- LHS - (P[i,,] %*% COV.mu %*% P[i,,]) }
      # set aside infinite uncertainty measurements
      loglike <- loglike + 1/2*log(abs(det(P[i,OBS[i,],OBS[i,]]))) - 1/2*c(D %*% P[i,,] %*% D)
    }

    K <- sqrtm(sigma,force=TRUE) %*% PDfunc(LHS,function(m){1/sqrt(m)},pseudo=TRUE)
    sigma <- K %*% RHS %*% t(K)
    # we are solving the unconstrained sigma above and then projecting back to the constrained sigma below
    if(isotropic) { sigma <- diag( mean(diag(sigma)), DIM ) }
    if(any(ZEROV))
    {
      sigma[ZEROV,] <- 0
      sigma[,ZEROV] <- 0
    }

    # Standardized error
    ERROR <- sigma - sigma.old # absolute error
    ERROR <- ERROR %*% ERROR # square to make positive
    K <- PDsolve(sigma)
    ERROR <- K %*% ERROR %*% K # standardize to make ~1 unitless
    ERROR <- sum(abs(diag(ERROR)),na.rm=TRUE)
  }

  # non-zero unique sigma parameters
  DUP <- upper.tri(sigma,diag=TRUE)
  DUP[ZEROV,] <- FALSE
  DUP[,ZEROV] <- FALSE

  # negative log-likelihood
  nlog.like <- function(par,zero=0,REML=debias)
  {
    if(isotropic)
    { sigma <- diag(par,DIM) }
    else
    {
      sigma <- array(0,c(DIM,DIM))
      sigma[DUP] <- par
      sigma <- t(sigma)
      sigma[DUP] <- par
    }

    #if(any(eigen(sigma,only.values=TRUE)$values<0)) { return(-Inf) }

    # estimate mu exactly | sigma
    P <- array(0,c(N,DIM,DIM))
    mu <- P.mu <- 0
    for(i in 1:N)
    {
      P[i,,] <- PDsolve(sigma + SIGMA[i,,])
      P.mu <- P.mu + P[i,,]
      mu <- mu + c(P[i,,] %*% MU[i,])
    }
    COV.mu <- PDsolve(P.mu)
    mu <- c(COV.mu %*% mu)
    if(any(ZEROM)) { mu[ZEROM] <- 0 }

    # sum up log-likelihood
    loglike <- REML/2*log(abs(det(COV.mu))) - DIM*(N-REML)/2*log(2*pi)
    for(i in 1:N)
    {
      D <- mu - MU[i,]
      loglike <- loglike + 1/2*log(abs(det(P[i,OBS[i,],OBS[i,]]))) - 1/2*c(D %*% P[i,,] %*% D)
    }
    loglike <- nant(loglike,-Inf)
    return(-loglike)
  }

  if(isotropic)
  {
    par <- mean(diag(sigma))

    parscale <- par
    lower <- 0
  }
  else
  {
    par <- sigma[DUP]

    parscale <- sqrt( diag(sigma) )
    parscale <- parscale %o% parscale
    parscale <- parscale[DUP]

    lower <- array(-Inf,dim(sigma))
    diag(lower) <- 0
    lower <- lower[DUP]
  }

  # implement zero if this is needed
  #par <- optimizer(par,nlog.like,parscale=parscale,lower=lower,upper=Inf)

  DIFF <- genD(par,nlog.like,parscale=parscale,lower=lower,upper=Inf)
  COV.sigma <- cov.loglike(DIFF$hessian,DIFF$gradient)

  loglike <- -nlog.like(par,REML=FALSE) # not REML for AIC/BIC

  # AIC
  n <- N
  q <- DIM
  qk <- sum(MEANS)
  if(isotropic)
  { nu <- 1 }
  else
  { nu <- sum(DUP) }
  K <- qk + nu

  AIC <- 2*K - 2*loglike
  AICc <- (q*n-qk)*2*K/max(q*n-K-nu,0) - 2*loglike
  BIC <- K*log(N) - 2*loglike

  names(mu) <- NAMES
  dimnames(COV.mu) <- list(NAMES,NAMES)
  dimnames(sigma) <- list(NAMES,NAMES)
  # sigma <- sigma[VARS,VARS]
  if(isotropic)
  { NAMES <- "sigma-sigma" }
  else
  {
    NAMES <- outer(NAMES,NAMES,function(x,y){paste0(x,"-",y)})
    NAMES <- NAMES[DUP]
  }
  dimnames(COV.sigma) <- list(NAMES,NAMES)

  return(list(mu=mu,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike,AIC=AIC,AICc=AICc,BIC=BIC,isotropic=isotropic))
}
