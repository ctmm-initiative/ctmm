meta.normal <- function(MU,SIGMA,VARS=TRUE,debias=TRUE,isotropic=FALSE,precision=1/2)
{
  if(length(dim(MU))<2)
  {
    N <- length(MU)
    DIM <- 1
    MU <- array(MU,c(N,1))
    SIGMA <- array(SIGMA,c(N,1,1))
  }
  else
  {
    N <- dim(MU)[1]
    DIM <- dim(MU)[2]
  }

  # do we give variance to each dimension
  VARS <- array(VARS,DIM)
  ZERO <- !VARS

  if(all(VARS))
  { return(meta.normal.all(MU,SIGMA,debias=debias,isotropic=isotropic,precision=precision)) }
  else if (all(ZERO))
  { return(meta.normal.null(MU,SIGMA,debias=debias,isotropic=isotropic,precision=precision)) }
  # else some of both

  STUFF.V <- meta.normal.all(MU[,VARS],SIGMA[,VARS,VARS],debias=debias,isotropic=isotropic,precision=precision)
  STUFF.Z <- meta.normal.all(MU[,ZERO],SIGMA[,ZERO,ZERO],debias=debias,isotropic=isotropic,precision=precision)

  mu <- array(0,DIM)
  mu[VARS] <- STUFF.V$MU
  mu[ZERO] <- STUFF.Z$MU

  sigma <- COV.mu <- array(0,c(DIM,DIM))
  sigma[VARS,VARS] <- STUFF.V$sigma
  sigma[ZERO,ZERO] <- STUFF.Z$sigma
  COV.mu[VARS,VARS] <- STUFF.V$COV.mu
  COV.mu[ZERO,ZERO] <- STUFF.Z$COV.mu

  DUP <- upper.tri(array(0,c(DIM,DIM)),diag=TRUE)
  COV.sigma <- array(0,c(1,1)*sum(DUP))
  DUP.V <- DUP # start with canonical DUP pairs
  DUP.V[ZERO,] <- DUP.V[,ZERO] <- 0 # only consider non-zero pairs
  DUP.V <- DUP.V[DUP] # subset that have variance
  DUP.Z <- DUP # start with canonical DUP pairs
  DUP.Z[VARS,] <- DUP.Z[,VARS] <- 0 # only consider zero pairs
  DUP.Z <- DUP.Z[DUP] # subset that don't have variance
  COV.sigma[VARS,VARS] <- STUFF.V$COV.sigma
  COV.sigma[ZERO,ZERO] <- STUFF.Z$COV.sigma

  loglike <- STUFF.V$loglike + STUFF.Z$loglike
  AIC <- STUFF.V$AIC + STUFF.Z$AIC
  AICc <- STUFF.V$AICc + STUFF.Z$AICc
  BIC <- STUFF.V$BIC + STUFF.Z$BIC

  return(list(mu=mu,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike,AIC=AIC,AICc=AICc,BIC=BIC,isotropic=isotropic))
}


meta.normal.null <- function(MU,SIGMA,debias=TRUE,isotropic=FALSE,precision=1/2)
{
  N <- dim(MU)[1]
  DIM <- dim(MU)[2]

  # estimate mu exactly
  P <- array(0,c(N,DIM,DIM))
  mu <- P.mu <- 0
  for(i in 1:N)
  {
    P[i,,] <- PDsolve(SIGMA[i,,])
    P.mu <- P.mu + P[i,,]
    mu <- mu + c(P[i,,] %*% MU[i,])
  }
  COV.mu <- PDsolve(P.mu)
  mu <- c(COV.mu %*% mu)

  # sum up log-likelihood
  loglike <- -DIM*(N)/2*log(2*pi)
  for(i in 1:N)
  {
    D <- mu - MU[i,]
    sigma <- array(SIGMA[i,,],c(DIM,DIM)) # R drops dim=1 and det() is not defined for scalars ... :(
    loglike <- loglike - 1/2*log(det(sigma)) - 1/2*c(D %*% P[i,,] %*% D)
  }

  # AIC
  n <- N
  q <- DIM
  k <- 1
  nu <- 0
  K <- q*k + nu

  AIC <- 2*K - 2*loglike
  AICc <- q*(n-k)*2*K/max(q*n-K-nu,0) - 2*loglike
  BIC <- K*log(N) - 2*loglike

  sigma <- array(0,c(DIM,DIM))
  DUP <- upper.tri(sigma,diag=TRUE)
  DUP <- sum(DUP)
  COV.sigma <- array(0,c(DUP,DUP))

  return(list(mu=mu,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike,AIC=AIC,AICc=AICc,BIC=BIC,isotropic=isotropic))
}


# population-level parameter estimates for normally distributed parameters and parameter uncertainties
meta.normal.all <- function(MU,SIGMA,debias=TRUE,isotropic=FALSE,precision=1/2)
{
  N <- dim(MU)[1]
  DIM <- dim(MU)[2]

  tol <- .Machine$double.eps^precision
  REML <- debias

  # initial guesses
  mu <- colMeans(MU)
  sigma <- 0
  for(i in 1:N) { sigma <- sigma + outer(MU[i,]-mu) }
  sigma <- sigma/(N-REML)

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

    loglike <- REML/2*log(det(COV.mu)) - DIM*(N-REML)/2*log(2*pi)
    # gradient with respect to sigma
    RHS <- 0
    LHS <- P.mu
    for(i in 1:N)
    {
      D <- mu - MU[i,]
      RHS <- RHS + (P[i,,] %*% outer(D) %*% P[i,,])
      if(debias) { LHS <- LHS - (P[i,,] %*% COV.mu %*% P[i,,]) }
      loglike <- loglike - 1/2*log(det(sigma + SIGMA[i,,])) - 1/2*c(D %*% P[i,,] %*% D)
    }

    K <- PDfunc(sigma,sqrt,force=TRUE) %*% PDfunc(LHS,function(m){1/sqrt(m)},pseudo=TRUE)
    sigma <- K %*% RHS %*% t(K)

    # Standardized error
    ERROR <- sigma - sigma.old # absolute error
    K <- PDfunc(sigma,function(m){1/sqrt(m)},pseudo=TRUE)
    ERROR <- K %*% ERROR %*% K # standardize to make ~1
    ERROR <- ERROR %*% ERROR # square to make positive
    ERROR <- sqrt(mean(diag(ERROR %*% ERROR))) # error in standard deviations
  }

  DUP <- upper.tri(sigma,diag=TRUE)

  # negative log-likelihood
  nlog.like <- function(par,zero=0,REML=debias)
  {
    sigma <- array(0,c(DIM,DIM))
    sigma[DUP] <- par
    sigma <- t(sigma)
    sigma[DUP] <- par

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

    # sum up log-likelihood
    loglike <- REML/2*log(det(COV.mu)) - DIM*(N-REML)/2*log(2*pi)
    for(i in 1:N)
    {
      D <- mu - MU[i,]
      loglike <- loglike - 1/2*log(det(sigma + SIGMA[i,,])) - 1/2*c(D %*% P[i,,] %*% D)
    }
    return(-loglike)
  }

  par <- sigma[DUP]

  parscale <- sqrt( diag(sigma) )
  parscale <- parscale %o% parscale
  parscale <- parscale[DUP]

  lower <- array(-Inf,c(DIM,DIM))
  diag(lower) <- 0
  lower <- lower[DUP]

  # implement zero if this is needed
  #par <- optimizer(par,nlog.like,parscale=parscale,lower=lower,upper=Inf)

  DIFF <- genD(par,nlog.like,parscale=parscale,lower=lower,upper=Inf)
  COV.sigma <- cov.loglike(DIFF$hessian,DIFF$gradient)

  loglike <- -nlog.like(par,REML=FALSE) # not REML for AIC/BIC

  # AIC
  n <- N
  q <- DIM
  k <- 1
  nu <- (q^2+q)/2
  K <- q*k + nu

  AIC <- 2*K - 2*loglike
  AICc <- q*(n-k)*2*K/max(q*n-K-nu,0) - 2*loglike
  BIC <- K*log(N) - 2*loglike

  return(list(mu=mu,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike,AIC=AIC,AICc=AICc,BIC=BIC,isotropic=isotropic))
}
