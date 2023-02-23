# TODO make sure this code works when there is no natural variance !!!!!!!!!!!!

# population-level parameter estimates for normally distributed parameters and parameter uncertainties
# VARS Boolean denotes whether or not there is natural variance
# MEANS Boolean denotes whether or not there is a non-zero mean
meta.normal <- function(MU,SIGMA,MEANS=TRUE,VARS=TRUE,diagonal=FALSE,isotropic=FALSE,debias=TRUE,weights=NULL,precision=1/2)
{
  tol <- .Machine$double.eps^precision
  REML <- debias

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

  if(is.null(weights))
  { weights <- rep(1,N) }
  else
  { weights <- weights/mean(weights) }

  # do we give variance to each dimension
  VARS <- array(VARS,DIM)
  # do we give a mean to each dimension
  MEANS <- array(MEANS,DIM)

  SCALE <- apply(MU,2,stats::mad)
  if(isotropic) { SCALE[] <- stats::median(SCALE) }
  ZERO <- SCALE<.Machine$double.eps
  if(any(ZERO)) # do not divide by zero or near zero
  {
    if(any(!ZERO))
    { SCALE[ZERO] <- min(SCALE[!ZERO]) } # do some rescaling... assuming axes are similar
    else
    { SCALE[ZERO] <- 1 } # do no rescaling
  }
  SCALE2 <- SCALE^2

  # observations
  OBS <- array(TRUE,c(N,DIM))
  for(i in 1:N) { OBS[i,] <- diag(cbind(SIGMA[i,,])) < SCALE2/.Machine$double.eps }
  # zero out Inf-VAR observations, in case of extreme point estimates
  MU[!OBS] <- 0
  # zero out Inf-VAR correlations, in case of extreme point estimates
  for(i in 1:N) { for(j in which(!OBS[i,])) { SIGMA[i,j,-j] <- SIGMA[i,-j,j] <- 0 } }

  NOBS <- colSums(OBS) # number of observations per parameter
  SUB <- NOBS<=1 # can't calculate variance for these
  if(any(SUB)) { VARS[SUB] <- FALSE }
  SUB <- NOBS==0 # can't calculate mean for these
  if(any(SUB)) { MEANS[SUB] <- FALSE }

  ZEROV <- !VARS
  ZEROM <- !MEANS

  ####################
  # robust rescaling in case of outliers with ill-conditioned COV matrices
  SHIFT <- apply(MU,2,stats::median)
  if(any(ZEROM)) { SHIFT[ZEROM] <- 0 }
  # apply shift
  MU <- t( t(MU) - SHIFT )
  # initial guess of mu
  mu <- SHIFT
  mu[] <- 0

  # apply scale
  MU <- t( t(MU)/SCALE )
  SIGMA <- aperm(SIGMA,c(2,3,1)) # [x,y,id]
  SIGMA <- SIGMA/SCALE
  SIGMA <- aperm(SIGMA,c(2,3,1)) # [y,id,x]
  SIGMA <- SIGMA/SCALE
  SIGMA <- aperm(SIGMA,c(2,3,1)) # [id,x,y]
  # initial guess of sigma
  sigma <- diag(1,nrow=length(mu))

  COV.mu <- sigma/N
  if(any(ZEROV)) { sigma[ZEROV,] <- sigma[,ZEROV] <- 0 }
  if(any(ZEROM)) { COV.mu[ZEROM,] <- COV.mu[,ZEROM] <- 0 }

  # non-zero unique sigma parameters
  UP <- upper.tri(sigma,diag=FALSE)
  UP[ZEROV,] <- FALSE
  UP[,ZEROV] <- FALSE

  DUP <- upper.tri(sigma,diag=TRUE)
  DUP[ZEROV,] <- FALSE
  DUP[,ZEROV] <- FALSE

  INF <- list(loglike=-Inf,mu=mu,COV.mu=sigma,sigma=sigma,sigma.old=sigma)

  # negative log-likelihood
  CONSTRAIN <- TRUE # constrain likelihood to positive definite
  COR <- TRUE # use correlations rather than covariances
  nloglike <- function(par,REML=debias,verbose=FALSE,zero=0)
  {
    zero <- -zero # now log-likelihood
    sigma <- par2sigma(par)

    # check for bad sigma matrices
    if(any(VARS))
    {
      V <- abs(diag(sigma[VARS,VARS]))
      V <- sqrt(V)
      # don't divide by zero
      TEST <- V<=.Machine$double.eps
      if(any(TEST)) { V[TEST] <- 1 }
      V <- V %o% V
      S <- sigma[VARS,VARS]/V
      S <- eigen(S)
      if(CONSTRAIN)
      {
        if(any(S$values<0))
        {
          if(!verbose)
          { return(Inf) }
          else
          { return(INF) }
        }
      }
      else
      {
        S$values <- pmax(S$values,.Machine$double.eps)
        S <- S$vectors %*% diag(S$values) %*% t(S$vectors)
        sigma[VARS,VARS] <- S * V
      }
    }

    # estimate mu exactly | sigma
    S <- P <- array(0,c(N,DIM,DIM))
    mu <- P.mu <- 0
    for(i in 1:N)
    {
      S[i,,] <- sigma + SIGMA[i,,]
      P[i,,] <- PDsolve(S[i,,],force=TRUE)
      P[i,,] <- PDclamp(P[i,,])
      P.mu <- P.mu + weights[i]*P[i,,]
      mu <- mu + weights[i]*c(P[i,,] %*% MU[i,])
    }
    COV.mu <- array(0,c(DIM,DIM))
    COV.mu[MEANS,MEANS] <- PDsolve(P.mu[MEANS,MEANS],force=TRUE)
    # COV.mu[MEANS,MEANS] <- cov.loglike(P.mu[MEANS,MEANS])
    mu <- c(COV.mu %*% mu)
    # if(any(ZEROM)) { mu[ZEROM] <- 0 } # should be okay

    # sum up log-likelihood
    loglike <- REML/2*log(abs(det(COV.mu[MEANS,MEANS,drop=FALSE]))) # - DIM*(N-REML)/2*log(2*pi)
    zero <- zero + DIM*(N-REML)/2*log(2*pi)
    zero <- 2*zero # for 1/2 below
    RHS <- 0
    LHS <- P.mu
    for(i in 1:N)
    {
      D <- mu - MU[i,]
      # set aside infinite uncertainty measurements
      loglike <- loglike - weights[i]/2*( PDlogdet(cbind(S[i,OBS[i,],OBS[i,]]),force=TRUE) + max(c(D %*% P[i,,] %*% D),0) + zero )

      # gradient with respect to sigma, under sum and trace
      if(verbose)
      {
        RHS <- RHS + weights[i]*(P[i,,] %*% outer(D) %*% P[i,,])
        if(debias) { LHS <- LHS - weights[i]*(P[i,,] %*% COV.mu %*% P[i,,]) }
      }
    }
    LHS <- unnant(LHS)
    loglike <- nant(loglike,-Inf)
    if(!verbose) { return(-loglike) }

    sigma.old <- sigma
    # update sigma
    if(any(VARS))
    {
      K <- sqrtm(sigma[VARS,VARS],pseudo=TRUE)
      K <- K %*% PDfunc(LHS[VARS,VARS],function(m){1/sqrt(m)},pseudo=TRUE)
      sigma[VARS,VARS] <- K %*% RHS[VARS,VARS] %*% t(K)
      # we are solving the unconstrained sigma above and then projecting back to the constrained sigma below
      # not sure how approximate this is, but will do numerical optimization afterwards
      if(isotropic)
      { sigma <- diag( mean(diag(sigma)), DIM ) }
      else if(diagonal)
      { sigma <- diag( diag(sigma) , DIM ) }
    }

    R <- list(loglike=loglike,mu=mu,COV.mu=COV.mu,sigma=sigma,sigma.old=sigma.old)
    return(R)
  }

  # extract parameters from sigma matrix
  sigma2par <- function(sigma)
  {
    if(isotropic)
    { par <- mean(diag(sigma)) }
    else if(diagonal)
    { par <- diag(sigma)[VARS] }
    else
    {
      if(COR)
      {
        D <- diag(sigma)
        d <- sqrt(D)
        d[d<=.Machine$double.eps] <- 1
        d <- d %o% d
        sigma <- sigma/d
        diag(sigma) <- D
      }

      par <- sigma[DUP]
    }

    return(par)
  }

  # construct sigma matrix from parameters
  par2sigma <- function(par)
  {
    if(isotropic)
    { sigma <- diag(par,DIM) }
    else if(diagonal)
    {
      sigma <- diag(0,DIM)
      diag(sigma)[VARS] <- par
    }
    else
    {
      sigma <- array(0,c(DIM,DIM))
      sigma[DUP] <- par
      sigma <- t(sigma)
      sigma[DUP] <- par

      if(COR)
      {
        D <- diag(sigma)
        d <- sqrt(D)
        d[d<=.Machine$double.eps] <- 1
        d <- d %o% d
        sigma <- sigma*d
        diag(sigma) <- D
      }
    }

    return(sigma)
  }

  ##############
  # zero sigma solution
  ZSOL <- nloglike(0,verbose=TRUE)

  #############
  # non zero sigma iterative solution
  if(any(VARS))
  {
    ERROR <- Inf
    SOL <- INF
  }
  else
  {
    ERROR <- 0
    SOL <- ZSOL
  }

  count <- 0
  count.0 <- 0
  while(ERROR>tol && count<100)
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
      K <- PDsolve(NSOL$sigma[VARS,VARS],pseudo=TRUE)
      ERROR <- K %*% ERROR %*% K # standardize to make ~1 unitless
      ERROR <- sum(abs(diag(ERROR)),na.rm=TRUE)
      count.0 <- 0
    }

    SOL <- NSOL
    count <- count + 1
  } # end while

  # end check
  if(ZSOL$loglike>SOL$loglike)
  { SOL <- ZSOL }

  # pull out results
  loglike <- SOL$loglike
  mu <- SOL$mu
  COV.mu <- SOL$COV.mu
  sigma <- SOL$sigma
  par <- sigma2par( SOL$sigma )

  parscale <- par
  lower <- 0
  upper <- Inf

  set.parscale <- function(known=FALSE)
  {
    parscale <<- par
    MIN <- 0

    if(isotropic)
    {
      # in case sigma is zero
      if(!known) { MIN <- mean(diag(COV.mu)) }
      lower <<- 0
      upper <<- Inf
    }
    else if(diagonal)
    {
      # in case sigma is zero
      if(!known) { MIN <- diag(COV.mu)[VARS] }
      lower <<- 0
      upper <<- Inf
    }
    else
    {
      if(COR)
      {
        parscale <<- sigma
        parscale[] <<- 1

        lower <<- array(-1,dim(sigma))
        upper <<- array(+1,dim(sigma))
      }
      else
      {
        parscale <<- sqrt( diag(sigma) )
        parscale <<- parscale %o% parscale

        lower <<- array(-Inf,dim(sigma))
        upper <<- array(+Inf,dim(sigma))
      }

      diag(parscale) <<- diag(sigma)
      parscale <<- parscale[DUP]

      diag(lower) <<- 0
      lower <<- lower[DUP]

      diag(upper) <<- Inf
      upper <<- upper[DUP]

      if(!known)
      {
        # in case sigma is zero
        if(COR)
        {
          MIN <- COV.mu
          MIN[] <- 1
        }
        else
        {
          MIN <- sqrt( abs( diag(COV.mu) ) )
          MIN <- MIN %o% MIN
        }

        diag(MIN) <- diag(COV.mu)
        MIN <- MIN[DUP]
      } # end unknown sigma
    } # end unstructured sigma

    MIN <- pmax(MIN,1)
    parscale <<- pmax(parscale,MIN)
  }
  set.parscale()

  # not sure if the iterative solution always works in constrained problems
  # if(any(ZEROV) || any(ZEROM))
  if(any(VARS))
  {
    # liberal parscale
    COR <- TRUE
    CONSTRAIN <- TRUE
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

    # uncertainty estimates
    COR <- FALSE
    CONSTRAIN <- FALSE # numderiv doesn't deal well with boundaries
    par <- sigma2par(sigma)
    set.parscale(TRUE) # more accurate parscale for numderiv
    DIFF <- genD(par,nloglike,parscale=parscale,lower=lower,upper=Inf)
    COV.sigma <- cov.loglike(DIFF$hessian,DIFF$gradient)
  }
  else
  {
    if(isotropic)
    { P <- 1 }
    else if(diagonal)
    { P <- DIM }
    else
    { P <- DIM*(DIM+1)/2 }
    COV.sigma <- matrix(0,P,P)
  }

  loglike <- -nloglike(par,REML=FALSE) # non-REML for AIC/BIC

  # rescale
  loglike <- loglike - log(prod(SCALE))
  mu <- mu*SCALE + SHIFT
  sigma <- t(sigma*SCALE)*SCALE
  COV.mu <- t(COV.mu*SCALE)*SCALE
  if(any(VARS))
  {
    SCALE <- SCALE %o% SCALE
    if(isotropic)
    { SCALE <- SCALE[1] }
    else if(diagonal)
    { SCALE <- diag(SCALE) }
    else
    { SCALE <- SCALE[DUP] }
    COV.sigma <- t(COV.sigma*SCALE)*SCALE
  }

  # AIC
  n <- N
  q <- DIM
  qk <- sum(MEANS)
  if(isotropic)
  { nu <- min(1,sum(VARS)) }
  else if(diagonal)
  { nu <- sum(VARS) }
  else
  { nu <- sum(DUP) }
  K <- qk + nu

  AIC <- 2*K - 2*loglike
  if(nu==0) # no variance parameters estimated, no bias
  { AICc <- AIC }
  else # some variance parameters estimated
  {
    AICc <- (q*n-qk)*2*K/max(q*n-K-nu,0) - 2*loglike
    AICc <- nant(AICc,Inf)
  }
  BIC <- K*log(N) - 2*loglike

  names(mu) <- NAMES
  dimnames(COV.mu) <- list(NAMES,NAMES)
  dimnames(sigma) <- list(NAMES,NAMES)
  # sigma <- sigma[VARS,VARS]
  if(isotropic)
  {
    NAMES2 <- "sigma-sigma"
    dimnames(COV.sigma) <- list(NAMES2,NAMES2)
  }
  else if(diagonal)
  {
    # basic structure
    NAMES2 <- NAMES[VARS]
    NAMES2 <- paste0(NAMES2,"-",NAMES2)
    dimnames(COV.sigma) <- list(NAMES2,NAMES2)

    # promote to general structure below
    NAMES3 <- outer(NAMES,NAMES,function(x,y){paste0(x,"-",y)})
    if(any(VARS))
    { NAMES3 <- NAMES3[DUP] }
    else
    { NAMES3 <- NAMES3[upper.tri(NAMES3,diag=TRUE)] }
    # zeroed COV matrix
    TEMP <- matrix(0,length(NAMES3),length(NAMES3))
    dimnames(TEMP) <- list(NAMES3,NAMES3)
    # copy over non-zero COV
    TEMP[NAMES2,NAMES2] <- COV.sigma
    COV.sigma <- TEMP
  }
  else
  {
    NAMES2 <- outer(NAMES,NAMES,function(x,y){paste0(x,"-",y)})
    if(any(VARS))
    { NAMES2 <- NAMES2[DUP] }
    else
    { NAMES2 <- NAMES2[upper.tri(NAMES2,diag=TRUE)] }
    dimnames(COV.sigma) <- list(NAMES2,NAMES2)
  }

  R <- list(mu=mu,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike,AIC=AIC,AICc=AICc,BIC=BIC,isotropic=isotropic)
  R$VARS <- VARS # need to pass this if variance was turned off due to lack of data
  R$MEANS <- MEANS
  return(R)
}
