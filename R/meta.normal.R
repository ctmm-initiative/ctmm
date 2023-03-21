# population-level parameter estimates for normally distributed parameters and parameter uncertainties
# MEANS Boolean denotes whether or not there is a non-zero mean
# VARS Boolean denotes whether or not there is natural variance-covariance
meta.normal <- function(MU,SIGMA,MEANS=TRUE,VARS=TRUE,isotropic=FALSE,GUESS=NULL,debias=TRUE,weights=NULL,precision=1/2)
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

  # do we give a mean to each dimension
  MEANS <- array(MEANS,DIM)
  # do we give a variance to each dimension after J-transformation
  if(length(VARS)==1) { VARS <- array(VARS,DIM) }
  if(length(dim(VARS))<2) { VARS <- outer(VARS,FUN="&") }
  if(isotropic) { VARS <- diag(diag(VARS)) }

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

  # observations
  OBS <- array(TRUE,c(N,DIM))
  for(i in 1:N) { OBS[i,] <- diag(cbind(SIGMA[i,,])) < SCALE^2/.Machine$double.eps }
  # zero out Inf-VAR observations, in case of extreme point estimates
  MU[!OBS] <- 0
  # zero out Inf-VAR correlations, in case of extreme point estimates
  for(i in 1:N) { for(j in which(!OBS[i,])) { SIGMA[i,j,-j] <- SIGMA[i,-j,j] <- 0 } }

  NOBS <- colSums(OBS) # number of observations per parameter
  SUB <- NOBS==0 # can't calculate mean for these
  if(any(SUB)) { MEANS[SUB] <- FALSE }
  SUB <- NOBS<=1 # can't calculate variance for these
  if(any(SUB)) { VARS[SUB,] <- VARS[,SUB] <- FALSE }
  vars <- diag(VARS)

  ####################
  # robust rescaling in case of outliers with ill-conditioned COV matrices
  SHIFT <- apply(MU,2,stats::median)
  if(any(!MEANS)) { SHIFT[!MEANS] <- 0 }
  # apply shift
  MU <- t( t(MU) - SHIFT )
  # initial guess of mu
  mu <- array(0,DIM)

  # apply scale
  MU <- t( t(MU)/SCALE )
  SIGMA <- aperm(SIGMA,c(2,3,1)) # [x,y,id]
  SIGMA <- SIGMA/SCALE
  SIGMA <- aperm(SIGMA,c(2,3,1)) # [y,id,x]
  SIGMA <- SIGMA/SCALE
  SIGMA <- aperm(SIGMA,c(2,3,1)) # [id,x,y]
  # initial guess of sigma
  sigma <- diag(diag(VARS))

  COV.mu <- sigma/N
  if(any(!MEANS)) { COV.mu[!MEANS,] <- COV.mu[,!MEANS] <- 0 }

  # potential unique sigma parameters
  TRI <- upper.tri(sigma,diag=TRUE)
  # non-zero unique sigma parameters
  DUP <- TRI & VARS

  # extract parameters from sigma matrix
  COR <- TRUE # use correlations rather than covariances
  sigma2par <- function(sigma)
  {
    if(isotropic)
    {
      par <- mean(sigma[VARS])
      par <- nant(par,0) # no VARS
    }
    else
    {
      if(COR)
      {
        D <- diag(sigma)
        d <- sqrt(D)
        d[!vars | d<=.Machine$double.eps] <- 1
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
    sigma <- diag(0,DIM)

    sigma[DUP] <- par
    sigma <- t(sigma)
    sigma[DUP] <- par

    if(!isotropic && COR)
    {
      D <- diag(sigma)
      d <- sqrt(abs(D))
      # d[d<=.Machine$double.eps] <- 1
      d <- d %o% d
      sigma <- sigma*d
      diag(sigma) <- D
    }

    return(sigma)
  }

  INF <- list(loglike=-Inf,mu=mu,COV.mu=sigma,sigma=sigma,sigma.old=sigma)
  if(!is.null(GUESS))
  {
    SUB <- MEANS & GUESS$mu!=0
    GUESS$mu <- (GUESS$mu-SHIFT)/SCALE
    INF$mu[SUB] <- GUESS$mu[SUB]

    SUB <- VARS & GUESS$sigma!=0
    GUESS$sigma <- t(GUESS$sigma/SCALE)/SCALE
    INF$sigma[SUB] <- GUESS$sigma[SUB]
    sigma <- INF$sigma

    INF$sigma.old <- INF$sigma
    INF$COV.mu <- INF$sigma/N
    if(any(!MEANS)) { INF$COV.mu[!MEANS,] <- INF$COV.mu[,!MEANS] <- 0 }
  }

  par <- sigma2par(sigma)

  # negative log-likelihood
  CONSTRAIN <- TRUE # constrain likelihood to positive definite
  nloglike <- function(par,REML=debias,verbose=FALSE,zero=0)
  {
    if(!verbose) { INF <- Inf }
    zero <- -zero # log-likelihood calculated below
    sigma <- par2sigma(par)

    # check for bad sigma matrices
    if(any(VARS))
    {
      # not sure how this happened once
      if(any(abs(par)==Inf)) { return(INF) }

      V <- abs(diag(sigma)[vars])
      V <- sqrt(V)
      # don't divide by zero
      TEST <- V<=.Machine$double.eps
      if(any(TEST)) { V[TEST] <- 1 }
      V <- V %o% V
      S <- sigma[vars,vars]/V
      S <- eigen(S)
      if(any(S$values<0))
      {
        if(CONSTRAIN)
        { return(INF) }
        else
        {
          S$values <- pmax(S$values,0)
          S <- S$vectors %*% diag(S$values,length(S$values)) %*% t(S$vectors)
          sigma[vars,vars] <- S * V
          sigma[!VARS] <- 0
        }
      } # end if negative variances
    } # end bad sigma check

    # estimate mu exactly | sigma
    S <- P <- P.inf <- array(0,c(N,DIM,DIM))
    mu <- mu.inf <- P.mu <- P.mu.inf <- 0
    for(i in 1:N)
    {
      S[i,,] <- sigma + SIGMA[i,,]
      P[i,,] <- PDsolve(S[i,,],force=TRUE) # finite and infinite weights
    }
    P.inf <- P==Inf # infinite weights
    P[P.inf] <- 0 # all finite weights now
    for(i in 1:N)
    {
      mu <- mu + weights[i]*c(P[i,,] %*% MU[i,])
      mu.inf <- mu.inf + weights[i]*(diag(P.inf[i,,]) * MU[i,])
      P.mu <- P.mu + weights[i]*P[i,,]
      P.mu.inf <- P.mu.inf + weights[i]*diag(P.inf[i,,])
    }
    COV.mu <- COV.mu.inf <- array(0,c(DIM,DIM))
    COV.mu[MEANS,MEANS] <- PDsolve(P.mu[MEANS,MEANS],force=TRUE)
    mu <- c(COV.mu %*% mu)
    # now handle infinite weight part
    SUB <- P.mu.inf>0
    mu[SUB] <- mu.inf[SUB]/P.mu.inf[SUB]
    # now put infinity back for likelihood
    P[P.inf] <- Inf

    # sum up log-likelihood
    loglike <- -REML/2*PDlogdet(P.mu[MEANS,MEANS]) # - DIM*(N-REML)/2*log(2*pi)
    zero <- 2/N*zero + DIM*(1-REML/N)*log(2*pi) # for below
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
      # we are solving the unconstrained sigma above and then projecting back to the constrained sigma below
      K <- sqrtm(sigma[vars,vars],pseudo=TRUE)
      K <- K %*% PDfunc(LHS[vars,vars],function(m){1/sqrt(m)},pseudo=TRUE)
      sigma[vars,vars] <- K %*% RHS[vars,vars] %*% t(K)

      # not sure how approximate this is, but will do numerical optimization afterwards
      sigma[!VARS] <- 0
      if(isotropic) { sigma[VARS] <- mean(sigma[VARS]) }
    }

    R <- list(loglike=loglike,mu=mu,COV.mu=COV.mu,sigma=sigma,sigma.old=sigma.old)
    return(R)
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
      ERROR <- (NSOL$sigma - SOL$sigma)[vars,vars] # absolute error
      ERROR <- ERROR %*% ERROR # square to make positive
      K <- PDsolve(NSOL$sigma[vars,vars],pseudo=TRUE)
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
      lower <<- 0
      upper <<- Inf

      # in case sigma is zero
      if(!known) { MIN <- mean(COV.mu[VARS]) }
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
    # HESS <- DIFF$hessian - outer(DIFF$gradient) # Hessian penalized by non-zero gradient
  }
  else
  { COV.sigma <- diag(0,DIM*(DIM+1)/2) }

  loglike <- -nloglike(par,REML=FALSE) # non-REML for AIC/BIC

  # rescale
  loglike <- loglike - log(prod(SCALE))
  mu <- SCALE*mu + SHIFT
  sigma <- t(sigma*SCALE)*SCALE
  COV.mu <- t(COV.mu*SCALE)*SCALE
  if(any(VARS))
  {
    SCALE <- SCALE %o% SCALE
    if(isotropic)
    { SCALE <- SCALE[1] }
    else
    { SCALE <- SCALE[DUP] }

    COV.sigma <- t(COV.sigma*SCALE)*SCALE
    # HESS <- t(HESS/SCALE)/SCALE
  }

  # AIC
  n <- N
  q <- DIM
  qn <- sum(OBS)
  qk <- sum(MEANS)
  if(isotropic)
  { nu <- max(VARS) }
  else
  { nu <- sum(DUP) }
  K <- qk + nu

  AIC <- 2*K - 2*loglike
  AICc <- (qn-qk)/max(qn-K-nu,0)*2*K - 2*loglike
  AICc <- nant(AICc,Inf)
  BIC <- K*log(N) - 2*loglike

  names(mu) <- NAMES
  dimnames(COV.mu) <- list(NAMES,NAMES)
  dimnames(sigma) <- list(NAMES,NAMES)
  dimnames(VARS) <- list(NAMES,NAMES)

  # full matrix to copy into - all possible matrix elements
  NAMES2 <- outer(NAMES,NAMES,function(x,y){paste0(x,"-",y)})
  COV.SIGMA <- diag(0,DIM^2)
  dimnames(COV.SIGMA) <- list(NAMES2,NAMES2)
  # copy over non-zero elements
  COV.SIGMA[which(DUP),which(DUP)] <- COV.sigma
  # unique matrix elements only
  COV.sigma <- COV.SIGMA[which(TRI),which(TRI)]
  # dimnames(HESS) <- list(NAMES2,NAMES2)

  R <- list(mu=mu,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike,AIC=AIC,AICc=AICc,BIC=BIC,isotropic=isotropic)
  R$VARS <- VARS # need to pass this if variance was turned off due to lack of data
  R$MEANS <- MEANS
  # R$HESS <- HESS
  return(R)
}
