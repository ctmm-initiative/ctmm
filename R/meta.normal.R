# TODO make sure this code works when there is no natural variance !!!!!!!!!!!!

# population-level parameter estimates for normally distributed parameters and parameter uncertainties
# VARS Boolean denotes whether or not there is natural variance
# MEANS Boolean denotes whether or not there is a non-zero mean
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

  # zero out Inf-VAR observations, in case of extreme point estimates
  MU[!OBS] <- 0

  tol <- .Machine$double.eps^precision
  REML <- debias

  ####################
  # initial non-zero guesses
  mu <- colMeans(MU)
  if(any(ZEROM)) { mu[ZEROM] <- 0 }

  sigma <- 0
  for(i in 1:N) { sigma <- sigma + outer(MU[i,]-mu) }
  sigma <- sigma/(N-REML)
  COV.mu <- sigma/N
  if(any(ZEROV)) { sigma[ZEROV,] <- sigma[,ZEROV] <- 0 }
  if(any(ZEROM)) { COV.mu[ZEROM,] <- COV.mu[,ZEROM] <- 0 }
  if(isotropic) { sigma <- diag( mean(diag(sigma)), DIM ) }

  # non-zero unique sigma parameters
  DUP <- upper.tri(sigma,diag=TRUE)
  DUP[ZEROV,] <- FALSE
  DUP[,ZEROV] <- FALSE

  # negative log-likelihood
  CONSTRAIN <- TRUE
  nloglike <- function(par,REML=debias,verbose=FALSE)
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

    # check for bad sigma matrices
    V <- abs(diag(sigma[VARS,VARS]))
    V <- sqrt(V)
    # don't divide by zero
    TEST <- V<=.Machine$double.eps
    if(any(TEST)) { V[TEST] <- 1 }
    V <- V %o% V
    S <- sigma[VARS,VARS]/V
    S <- eigen(S)
    if(CONSTRAIN)
    { if(any(S$values<0)) { return(Inf) } }
    else
    {
      S$values <- pmax(S$values,.Machine$double.eps)
      S <- S$vectors %*% diag(S$values) %*% t(S$vectors)
      sigma[VARS,VARS] <- S * V
    }

    # estimate mu exactly | sigma
    P <- array(0,c(N,DIM,DIM))
    mu <- P.mu <- 0
    for(i in 1:N)
    {
      P[i,,] <- PDsolve(sigma + SIGMA[i,,])
      P.mu <- P.mu + P[i,,]
      mu <- mu + c(P[i,,] %*% MU[i,])
    }
    COV.mu <- array(0,c(DIM,DIM))
    COV.mu[MEANS,MEANS] <- PDsolve(P.mu[MEANS,MEANS])
    mu <- c(COV.mu %*% mu)
    # if(any(ZEROM)) { mu[ZEROM] <- 0 } # should be okay

    # sum up log-likelihood
    loglike <- REML/2*log(abs(det(COV.mu[MEANS,MEANS]))) - DIM*(N-REML)/2*log(2*pi)
    RHS <- 0
    LHS <- P.mu
    for(i in 1:N)
    {
      D <- mu - MU[i,]
      # set aside infinite uncertainty measurements
      loglike <- loglike + 1/2*log(abs(det(P[i,OBS[i,],OBS[i,]]))) - 1/2*c(D %*% P[i,,] %*% D)

      # gradient with respect to sigma, under sum and trace
      if(verbose)
      {
        RHS <- RHS + (P[i,,] %*% outer(D) %*% P[i,,])
        if(debias) { LHS <- LHS - (P[i,,] %*% COV.mu %*% P[i,,]) }
      }
    }
    loglike <- nant(loglike,-Inf)
    if(!verbose) { return(-loglike) }

    K <- sqrtm(sigma[VARS,VARS],pseudo=TRUE)
    K <- K %*% PDfunc(LHS[VARS,VARS],function(m){1/sqrt(m)},pseudo=TRUE)
    sigma.old <- sigma
    sigma[VARS,VARS] <- K %*% RHS[VARS,VARS] %*% t(K)
    # we are solving the unconstrained sigma above and then projecting back to the constrained sigma below
    # not sure how approximate this is, but will do numerical optimization afterwards
    if(isotropic) { sigma <- diag( mean(diag(sigma)), DIM ) }

    R <- list(loglike=loglike,mu=mu,COV.mu=COV.mu,sigma=sigma,sigma.old=sigma.old)
    return(R)
  }

  # extract parameters from sigma matrix
  sigma2par <- function(sigma)
  {
    if(isotropic)
    { par <- mean(diag(sigma)) }
    else
    { par <- sigma[DUP] }

    return(par)
  }

  ##############
  # zero sigma solution
  ZSOL <- nloglike(0,verbose=TRUE)

  #############
  # non zero sigma iterative solution
  count.0 <- 0
  ERROR <- Inf
  SOL <- list(loglike=-Inf,sigma=sigma)
  while(ERROR>tol)
  {
    par <- sigma2par( SOL$sigma )
    NSOL <- nloglike(par,verbose=TRUE)

    # did likelihood fail to increase?
    if(NSOL$loglike<=SOL$loglike)
    {
      SOL$sigma <- SOL$sigma.old
      break
    }

    # is sigma shrinking to zero?
    if(ZSOL$loglike>NSOL$loglike)
    {
      if(count.0>10) # will converge to zero slowly
      {
        SOL <- ZSOL
        break
      }
      else # accelerate to zero
      {
        NSOL$sigma <- NSOL$sigma/2
        count.0 <- count.0 + 1
      }
    }
    else # proceed as usual
    {
      # Standardized error
      ERROR <- (NSOL$sigma - SOL$sigma)[VARS,VARS] # absolute error
      ERROR <- ERROR %*% ERROR # square to make positive
      K <- PDsolve(NSOL$sigma[VARS,VARS])
      ERROR <- K %*% ERROR %*% K # standardize to make ~1 unitless
      ERROR <- sum(abs(diag(ERROR)),na.rm=TRUE)
      count.0 <- 0
    }

    SOL <- NSOL
  } # end while
  # pull out results
  loglike <- SOL$loglike
  mu <- SOL$mu
  COV.mu <- SOL$COV.mu
  sigma <- SOL$sigma
  par <- sigma2par( SOL$sigma )

  if(isotropic)
  {
    # in case sigma is zero
    MIN <- mean(diag(COV.mu))

    parscale <- pmax(par,MIN)
    lower <- 0
  }
  else
  {
    parscale <- sqrt( diag(sigma) )
    parscale <- parscale %o% parscale
    parscale <- parscale[DUP]

    # in case sigma is zero
    MIN <- sqrt( abs( diag(COV.mu) ) )
    MIN <- MIN %o% MIN
    MIN <- MIN[DUP]

    parscale <- pmax(parscale,MIN)

    lower <- array(-Inf,dim(sigma))
    diag(lower) <- 0
    lower <- lower[DUP]
  }

  # not sure if the iterative solution always works in constrained problems
  # if(any(ZEROV) || any(ZEROM))
  if(TRUE)
  {
    SOL <- optimizer(par,nloglike,parscale=parscale,lower=lower,upper=Inf)
    par <- SOL$par
    loglike <- -SOL$value

    # needed for optimization and backtracking
    SOL <- nloglike(par,verbose=TRUE)
    loglike <- SOL$loglike
    mu <- SOL$mu
    COV.mu <- SOL$COV.mu
    sigma <- SOL$sigma.old
    par <- sigma2par(sigma)
  }

  # uncertainty estimates
  CONSTRAIN <- FALSE # numderiv doesn't deal well with boundaries
  DIFF <- genD(par,nloglike,parscale=parscale,lower=lower,upper=Inf)
  COV.sigma <- cov.loglike(DIFF$hessian,DIFF$gradient)

  loglike <- -nloglike(par,REML=FALSE) # non-REML for AIC/BIC

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
