# population-level parameter estimates for normally distributed parameters and parameter uncertainties
# Y: response point estimates
# SY: response uncertainty covariance
# INT Boolean denotes whether or not there is an intercept
# VARS Boolean denotes whether or not there is process variance-covariance in Y (not estimation uncertainty)
# isotropic: VARS is an isotropic covariance matrix
# X: predictor point estimates
# SX: predictor uncertainty variance
# without X, slopes will not be estimated
# D: design submatrix: Y = X %*% DSM %*% beta + ...
meta.normal <- function(Y,SY=FALSE,X=FALSE,SX=FALSE,DSM=NULL,INT=TRUE,VARS=TRUE,isotropic=FALSE,GUESS=NULL,debias=TRUE,weights=NULL,precision=1/2,WARN=TRUE,...)
{
  tol <- .Machine$double.eps^precision
  REML <- debias
  UNBIASED <- TRUE # use BLUE beta estimator, rather than maximum likelihood (biased)

  if(length(dim(Y))<2)
  {
    NAMES <- "x"
    N <- length(Y)
    DIM <- 1
    Y <- array(Y,c(N,1))
    SY <- array(SY,c(N,1,1))
  }
  else
  {
    NAMES <- colnames(Y)
    N <- dim(Y)[1]
    DIM <- dim(Y)[2]
  }

  if(is.null(weights))
  { weights <- rep(1,N) }
  else
  { weights <- weights/mean(weights) }

  # do we give an intercept to each dimension
  INT <- array(INT,DIM)
  names(INT) <- NAMES
  WINT <- which(INT)
  names(WINT) <- NAMES[INT]
  # do we give a variance to each dimension after J-transformation
  if(length(VARS)==1) { VARS <- array(VARS,DIM) }
  if(length(dim(VARS))<2) { VARS <- outer(VARS,FUN="&") }
  if(isotropic) { VARS <- diag(diag(VARS)) }
  dimnames(VARS) <- list(NAMES,NAMES)

  # finite observations
  OBS <- t( apply(SY,1,function(s){diag(cbind(s))<Inf}) ) # [DIM,N]
  dim(OBS) <- c(DIM,N) # R drops dimensions :(
  OBS <- t(OBS)

  # natural scales of the Y-data
  SHIFT <- rep(0,DIM)
  SCALE <- rep(1,DIM)
  for(i in 1:DIM)
  {
    TEMP <- Y[OBS[,i],i]
    if(is.null(TEMP))
    {
      SHIFT[i] <- 0
      SCALE[i] <- 1
    }
    else
    {
      SHIFT[i] <- stats::median(TEMP)
      SCALE[i] <- stats::mad(TEMP)

      if(SCALE[i]<.Machine$double.eps) { SCALE[i] <- sqrt(stats::var(TEMP)) }
      if(is.na(SCALE[i]) || SCALE[i]<.Machine$double.eps) { SCALE[i] <- 1 }
    }
  }
  SHIFT[!INT] <- 0 # don't shift if no intercept
  if(isotropic) { SCALE[] <- stats::median(SCALE) }

  # well-conditioned Y-observations (more stringent)
  OBS <- apply(SY,1,function(s){diag(cbind(s))<SCALE^2/.Machine$double.eps}) # [DIM,N]
  dim(OBS) <- c(DIM,N) # R drops dimensions :(
  OBS <- t(OBS)

  if(any(X!=0))
  {
    if(is.null(dim(X))) { X <- array(X,c(N,1)) }
    XDIM <- ncol(X)

    if(is.null(DSM)) # design submatrix inferred from X-columns
    {
      if(XDIM==DIM)
      {
        DSM <- diag(1,nrow=DIM)
        dimnames(DSM) <- list(NAMES,NAMES)
        XNAMES <- colnames(X) <- NAMES
      }
      else if(XDIM<DIM)
      {
        XNAMES <- colnames(X)
        DSM <- array(0,c(XDIM,DIM))
        dimnames(DSM) <- list(XNAMES,NAMES)
        for(x in XNAMES) { DSM[x,x] <- 1 }
      }
    }

    if(length(dim(SX))==2) # promote variance to covariance
    {
      TEMP <- SX
      SX <- vapply(1:nrow(SX),function(i){diag(SX[i],nrow=XDIM)},diag(0,nrow=XDIM))
      SX <- aperm(SX,c(3,1,2))
      dimnames(SX) <- list(NULL,XNAMES,XNAMES)
    }

    XOBS <- t( apply(SX,1,function(s){diag(cbind(s))<Inf}) ) # [X,N]
    dim(XOBS) <- c(XDIM,N) # R drops dimensions :(
    XOBS <- t(XOBS)

    XSHIFT <- rep(0,XDIM)
    XSCALE <- rep(1,XDIM)
    for(i in 1:XDIM)
    {
      TEMP <- X[XOBS[,i],i]
      if(is.null(TEMP))
      {
        XSHIFT[i] <- 0
        XSCALE[i] <- 1
      }
      else
      {
        XSHIFT[i] <- stats::median(TEMP)
        XSCALE[i] <- stats::mad(TEMP)

        if(XSCALE[i]<.Machine$double.eps) { XSCALE[i] <- sqrt(stats::var(TEMP)) }
        if(XSCALE[i]<.Machine$double.eps) { XSCALE[i] <- 1 }
      }
    }

    # well-conditioned x-observations (more stringent)
    XOBS <- apply(SX,1,function(s){diag(cbind(s))<XSCALE^2/.Machine$double.eps}) # [X,N]
    dim(XOBS) <- c(XDIM,N) # R drops dimensions :(
    XOBS <- t(XOBS)

    # don't shift if no intercept
    SUB <- apply(DSM,1,function(D){sum(D[!INT])}) # [Y,X] # [X]
    SUB <- as.logical(SUB)
    XSHIFT[SUB] <- 0

    X <- t( t(X)-XSHIFT )
    X <- t( t(X)/XSCALE )

    SX <- aperm(SX,c(3,1:2))
    SX <- SX/XSCALE
    SX <- aperm(SX,c(3,1:2))
    SX <- SX/XSCALE
    SX <- aperm(SX,c(3,1:2))

    # zero out Inf-VAR correlations, in case of extreme point estimates
    for(i in 1:N) { for(j in which(!XOBS[i,])) { SX[i,j,-j] <- SX[i,-j,j] <- 0 } }

    # zero out Inf-VAR observations, in case of extreme point estimates
    X[!XOBS] <- 0
  } # end if(any(X))
  else
  { XNAMES <- BSCALE <- NULL }

  bDSM <- as.logical(DSM)
  dim(bDSM) <- dim(DSM) # R drops dimensions :(
  BDIM <- sum(bDSM)
  iDSM <- which(bDSM,arr.ind=TRUE)
  if(BDIM)
  {
    BX <- iDSM[,1]
    BY <- iDSM[,2]

    BSCALE <- SCALE[BY]/XSCALE[BX]

    BETA <- array(0,dim(DSM))

    # can't use y if can't use associated x
    OBS <- OBS & !( (!XOBS) %*% BETA ) # [N,Y]
  }
  else
  { BX <- BY <- NULL }

  # zero out Inf-VAR observations, in case of extreme point estimates
  Y[!OBS] <- 0
  # zero out Inf-VAR correlations, in case of extreme point estimates
  for(i in 1:N) { for(j in which(!OBS[i,])) { SY[i,j,-j] <- SY[i,-j,j] <- 0 } }

  NOBS <- colSums(OBS) # number of observations per parameter
  SUB <- NOBS==0 # can't calculate mean for these
  if(any(SUB)) { INT[SUB] <- FALSE }
  SUB <- NOBS<=1 # can't calculate variance for these
  if(any(SUB)) { VARS[SUB,] <- VARS[,SUB] <- FALSE }
  vars <- diag(VARS)

  ####################
  # robust rescaling in case of outliers with ill-conditioned COV matrices
  Y <- t( t(Y) - SHIFT )
  Y[!OBS] <- 0
  # initial guess of mu
  MDIM <- sum(INT)
  mu <- array(0,DIM)
  beta <- array(0,BDIM)

  # sorting indices for mu = c(beta,mu)
  BI <- 1%:%BDIM
  MI <- BDIM + 1%:%MDIM

  # Jacobian: [(beta,mu),Y]
  JAC <- array(0,c(BDIM+MDIM,DIM))
  for(j in 1%:%MDIM) { JAC[MI[j],WINT[j]] <- 1 }

  # apply scale
  Y <- t( t(Y)/SCALE )
  SY <- aperm(SY,c(2,3,1)) # [x,y,id]
  SY <- SY/SCALE
  SY <- aperm(SY,c(2,3,1)) # [y,id,x]
  SY <- SY/SCALE
  SY <- aperm(SY,c(2,3,1)) # [id,x,y]
  # initial guess of sigma
  sigma <- diag(diag(VARS),nrow=DIM)

  COV.mu <- array(0,c(1,1)*(DIM+MDIM))

  # potential unique sigma parameters
  TRI <- upper.tri(sigma,diag=TRUE)
  # non-zero unique sigma parameters
  DUP <- TRI & VARS

  if(isotropic)
  { VDIM <- max(VARS) }
  else
  { VDIM <- sum(DUP) }

  # sorting indices for par = c(beta,sigma)
  VI <- BDIM + 1%:%VDIM

  # extract parameters from sigma matrix
  COR <- TRUE # use correlations rather than covariances
  sigma2par <- function(BETA,sigma)
  {
    if(length(dim(BETA)))
    { BETA <- BETA[bDSM] }
    beta <- BETA

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

    return(c(beta,par))
  }

  # construct sigma matrix from parameters
  par2sigma <- function(par)
  {
    if(BDIM)
    {
      beta <- par[1:length(BDIM)]
      # names(beta) <- XNAMES

      BETA <- array(0,c(XDIM,DIM))
      rownames(BETA) <- XNAMES
      colnames(BETA) <- NAMES
      BETA[iDSM] <- beta

      par <- par[-(1:length(BDIM))]
    }
    else
    { BETA <- NULL }

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

    dimnames(sigma) <- list(NAMES,NAMES)

    return(list(BETA=BETA,sigma=sigma))
  }

  INF <- list(loglike=-Inf,mu=mu,COV.mu=sigma,sigma=sigma,sigma.old=sigma,beta=beta,beta.old=beta)

  # negative log-likelihood with mu (intercept) exactly profiled
  CONSTRAIN <- TRUE # constrain likelihood to positive definite
  nloglike <- function(par,beta.fixed=FALSE,REML=debias,verbose=FALSE,zero=0)
  {
    if(!verbose) { INF <- Inf }
    zero <- -zero # log-likelihood calculated below
    TEMP <- par2sigma(par)
    sigma <- TEMP$sigma
    BETA <- TEMP$BETA

    if(beta.fixed)
    { Y <- Y - (X %*% BETA) }

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

    S <- P <- P.inf <- array(0,c(N,DIM,DIM))
    dimnames(S) <- dimnames(P) <- dimnames(P.inf) <- list(NULL,NAMES,NAMES)
    mu <- mi <- P.mi <- array(0,BDIM+MDIM) # combine mu <- c(beta,mu)
    names(mu) <- names(mi) <- names(P.mi) <- c(XNAMES,NAMES)
    COV.mu <- P.mu <- array(0,c(BDIM+MDIM,BDIM+MDIM))
    dimnames(P.mu) <- list(c(XNAMES,NAMES),c(XNAMES,NAMES))

    for(i in 1:N)
    {
      S[i,,] <- sigma + SY[i,,]
      if(BDIM) { S[i,,] <- S[i,,] + t(BETA) %*% SX[i,,] %*% BETA } # model-2 regression
      P[i,,] <- PDsolve(S[i,,],force=TRUE) # finite and infinite weights
    }
    P.INF <- P==Inf # infinite weights
    P[P.INF] <- 0 # all finite weights now, will re-introduce Inf after to avoid NaNs

    # mu, beta | beta, sigma
    for(i in 1:N)
    {
      # infinite precision terms (diagonal only
      Pi <- diag(cbind(P.INF[i,,])) # R will drop dimensions here

      if(MDIM)
      {
        mu[MI] <- mu[MI] + weights[i] * c(P[i,INT,] %*% Y[i,])
        mi[MI] <- mi[MI] + weights[i] * Pi[INT] * Y[i,INT] # infinitely weighted contributions from diagonal

        P.mu[MI,MI] <- P.mu[MI,MI] + weights[i] * P[i,INT,INT]
        P.mi[MI] <- P.mi[MI] + weights[i] * Pi[INT]
      }

      if(BDIM)
      {
        mu[BI] <- mu[BI] + weights[i] * c(X[i,BX] %*% P[i,BY,] %*% Y[i,])
        mi[BI] <- mi[BI] + weights[i] * X[i,BX] * Pi[BY] * Y[i,BY]

        P.mu[BI,BI] <- P.mu[BI,BI] + weights[i] * c(X[i,BX] %*% P[i,BY,BY] %*% X[i,BX])
        P.mi[BI] <- P.mi[BI] + weights[i] * X[i,BX] * Pi[BY] * X[i,BX]
      }

      if(BDIM && MDIM)
      {
        P.mu[BI,MI] <- P.mu[BI,MI] + weights[i] * c(X[i,BX] %*% P[i,BY,INT])
        P.mu[MI,BI] <- P.mu[MI,BI] + weights[i] * c(cbind(P[i,INT,BY]) %*% X[i,BX]) # R drops dimensions :(
      }
    } # end for(i in 1:N)

    if(beta.fixed) # mu | beta, var
    { COV.mu[MI,MI] <- PDsolve(P.mu[MI,MI],force=TRUE) }
    else # mu, beta | var
    { COV.mu <- PDsolve(P.mu,force=TRUE) }
    mu <- c(COV.mu %*% mu)
    if(beta.fixed)
    { mu[BI] <- par[BI] }

    # now handle infinite weight part
    SUB <- (P.mi>0)
    mu[SUB] <- mi[SUB]/P.mi[SUB]
    # now put infinity back for likelihood
    P[P.INF] <- Inf

    # sum up log-likelihood
    zero <- 2/N*(zero + (N*DIM-REML*(MDIM+BDIM))/2*log(2*pi) ) # for below
    loglike <- 0
    # REML correction # - DIM*(N-REML)/2*log(2*pi)
    loglike <- loglike - REML/2*PDlogdet(P.mu) # approximate but complete
    # loglike <- loglike - REML/2*PDlogdet(P.mu[MI,MI]) # incomplete but exact

    MU <- numeric(DIM)
    MU[INT] <- mu[MI]
    # updated BETA for residuals
    if(!beta.fixed)
    { BETA[iDSM] <- mu[BI] }

    Y <- t( t(Y) - MU )
    if(BDIM && !beta.fixed)
    { Y <- Y - (X %*% BETA) }

    RHS <- LHS <- array(0,c(DIM,DIM))
    for(i in 1:N)
    {
      # set aside infinite uncertainty measurements
      loglike <- loglike - weights[i]/2*( PDlogdet(cbind(S[i,OBS[i,],OBS[i,]]),force=TRUE) + max(c(Y[i,] %*% P[i,,] %*% Y[i,]),0) + zero )

      # gradient with respect to sigma, under sum and trace
      if(verbose)
      {
        RHS <- RHS + weights[i] * (P[i,,] %*% outer(Y[i,]) %*% P[i,,])
        LHS <- LHS + weights[i] * P[i,,]

        # smaller REML contribution # I think this code is approximate for structured matrices
        if(debias && (BDIM+MDIM))
        {
          # jacobian terms (X)
          for(j in 1%:%BDIM) { JAC[BI[j],BY[j]] <- X[i,BX[j]] }

          # individual contribution
          J <- JAC %*% P[i,,]
          J <- t(J) %*% COV.mu %*% J

          LHS <- LHS - weights[i] * J
        }
      }
    } # end if(i in 1:N)

    loglike <- nant(loglike,-Inf)
    if(!verbose) { return(-loglike) }

    LHS <- fixNaN(LHS)

    sigma.old <- sigma
    # update sigma
    # LHS == RHS # 1/sigma == 1/sigma %*% S %*% 1/sigma
    if(any(VARS))
    {
      LHS <- LHS[vars,vars]
      RHS <- RHS[vars,vars]

      if(isotropic)
      {
        LHS[VARS] <- mean(LHS[VARS])
        RHS[VARS] <- mean(RHS[VARS])
      }

      # method 1 (faster)
      # K <- PDsolve(LHS) # sigma

      # method 2 (slower)
      K <- sqrtm(sigma[vars,vars],pseudo=TRUE) %*% PDfunc(LHS,function(m){1/sqrt(m)},pseudo=TRUE) # sigma

      sigma[vars,vars] <- K %*% RHS %*% t(K) # S
      sigma[!VARS] <- 0

      if(isotropic) { sigma[VARS] <- mean(sigma[VARS]) }
    }

    beta.old <- par[BI]
    beta <- mu[BI]
    mu <- mu[MI]

    R <- list(loglike=loglike,mu=mu,COV.mu=COV.mu,sigma=sigma,sigma.old=sigma.old,beta=beta,beta.old=beta.old)
    return(R)
  }

  # profile beta
  if(BDIM)
  {
    nloglike.profile <- function(par,verbose=FALSE,zero=0,...)
    {
      SOL <- nloglike(par,verbose=TRUE,zero=zero,...)

      ERROR <- Inf
      count <- 0
      while(ERROR>tol && count<100)
      {
        par[BI] <- SOL$beta
        NSOL <- nloglike(par,verbose=TRUE,zero=zero,...)

        E <- (NSOL$beta - SOL$beta)^2 / pmax( NSOL$beta^2 , diag(cbind(NSOL$COV.mu))[BI], 1 )
        ERROR <- sum(E,na.rm=TRUE)

        SOL <- NSOL
        count <- count + 1
      } # end while

      if(!verbose) { SOL <- -SOL$loglike }
      return(SOL)
    }

    # don't supply beta and profile beta
    nloglike.nobeta <- function(par,zero=0,...)
    {
      par <- c(beta,par)
      nloglike.profile(par,zero=zero,...)
    }

    # fix beta, do not profile beta
    nloglike.bfixed <- function(par,zero=0,...) { nloglike(par,beta.fixed=TRUE,zero=zero,...) }
  } # end if(BDIM)
  else
  {
    nloglike.profile <- nloglike
    nloglike.nobeta <- nloglike
    nloglike.bfixed <- nloglike
  }

  par <- sigma2par(beta,sigma)
  SOL <- nloglike.profile(par,verbose=TRUE)
  ZSOL <- nloglike.profile(0*par,verbose=TRUE)
  # starting point from previous model
  if(!is.null(GUESS))
  {
    # GUESS$mu <- (GUESS$mu-SHIFT)/SCALE
    GUESS$sigma <- t(GUESS$sigma/SCALE)/SCALE
    if(BDIM)
    { GUESS$beta <- GUESS$beta/BSCALE }
    GUESS$par <- sigma2par(GUESS$beta,GUESS$sigma)

    GSOL <- nloglike.profile(GUESS$par,verbose=TRUE)
    if(GSOL$loglike>SOL$loglike) { SOL <- GSOL }
  }

  ERROR <- Inf
  count <- 0
  count.0 <- 0
  while(ERROR>tol && count<100)
  {
    par <- sigma2par( SOL$beta, SOL$sigma )
    NSOL <- nloglike.profile(par,verbose=TRUE)

    # did likelihood fail to increase?
    if(NSOL$loglike<=SOL$loglike)
    {
      # SOL$beta <- SOL$beta.old
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
      ERROR <- 0

      # Standardized error^2 from unknown sigma
      if(any(VARS))
      {
        E <- (NSOL$sigma - SOL$sigma)[vars,vars] # absolute error
        E <- E %*% E # square to make positive
        K <- PDsolve(NSOL$sigma[vars,vars],pseudo=TRUE)
        E <- K %*% E %*% K # standardize to make ~1 unitless
        ERROR <- sum(abs(diag(E)),na.rm=TRUE)
      }

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
  beta <- SOL$beta
  sigma <- SOL$sigma
  par <- sigma2par( SOL$beta, SOL$sigma )

  parscale <- pmax( par, 1 )
  lower <- c( rep(-Inf,BDIM) , rep(0,VDIM) )
  upper <- array(Inf,BDIM+VDIM)

  set.parscale <- function(known=FALSE)
  {
    PS <- par[VI]
    MIN <- 0

    if(isotropic)
    {
      lower[VI] <<- 0
      upper[VI] <<- Inf

      # in case sigma is zero
      if(!known) { MIN <- mean(COV.mu[VARS]) }
    }
    else
    {
      if(COR)
      {
        PS <- array(1,dim(sigma))

        LO <- array(-1,dim(sigma))
        UP <- array(+1,dim(sigma))
      }
      else
      {
        PS <- sqrt( diag(sigma) )
        PS <- PS %o% PS

        LO <- array(-Inf,dim(sigma))
        UP <- array(+Inf,dim(sigma))
      }

      diag(PS) <- diag(sigma)
      PS <- PS[DUP]

      diag(LO) <- 0
      diag(UP) <- Inf

      lower[VI] <<- LO[DUP]
      upper[VI] <<- UP[DUP]

      if(!known)
      {
        # in case sigma is zero
        if(COR)
        {
          MIN <- COV.mu[MI,MI]
          MIN[] <- 1
        }
        else
        {
          MIN <- sqrt( abs( diag(COV.mu)[MI] ) )
          MIN <- MIN %o% MIN
        }

        diag(MIN) <- diag(COV.mu)[MI]
        MIN <- MIN[DUP]

        if(isotropic) { MIN <- max(MIN) }
      } # end unknown sigma
    } # end unstructured sigma

    MIN <- pmax(MIN,1)
    parscale[VI] <<- pmax(PS,MIN)
  }
  set.parscale()

  # iterative solution does not always work
  if(VDIM || (BDIM && !UNBIASED)) # VDIM+BDIM if ML beta
  {
    # liberal parscale
    COR <- TRUE
    CONSTRAIN <- TRUE
    set.parscale()
    ## self-consistent - unbiased beta
    if(UNBIASED)
    { SOL <- optimizer(par[VI],nloglike.nobeta,parscale=parscale[VI],lower=lower[VI],upper=upper[VI]) }
    else
    { SOL <- optimizer(par,nloglike.bfixed,parscale=parscale,lower=lower,upper=upper) }
    par <- SOL$par
    loglike <- -SOL$value

    # needed for optimization and backtracking
    if(UNBIASED)
    { SOL <- nloglike.nobeta(par,verbose=TRUE) }
    else
    { SOL <- nloglike.bfixed(par,verbose=TRUE) }
    loglike <- SOL$loglike
    mu <- SOL$mu
    COV.mu <- SOL$COV.mu
    beta <- SOL$beta
    sigma <- SOL$sigma.old
    par <- sigma2par(beta,sigma)

    # uncertainty estimates
    COR <- FALSE
    CONSTRAIN <- FALSE # numderiv doesn't deal well with boundaries
    par <- sigma2par(beta,sigma)
    set.parscale(TRUE) # more accurate parscale for numderiv
    if(UNBIASED)
    { DIFF <- genD(par[VI],nloglike.nobeta,parscale=parscale[VI],lower=lower[VI],upper=upper[VI]) }
    else
    { DIFF <- genD(par,nloglike.bfixed,parscale=parscale,lower=lower,upper=Inf) }
    COV.sigma <- cov.loglike(DIFF$hessian,DIFF$gradient,WARN=WARN)
    # HESS <- DIFF$hessian - outer(DIFF$gradient) # Hessian penalized by non-zero gradient

    # correlation correction
    if(BDIM && !UNBIASED)
    {
      C <- stats::cov2cor( nloglike(par,verbose=TRUE)$COV.mu ) # approximate beta-mu correlations
      VAR <- c( diag(COV.sigma)[BI], diag(COV.mu)[MI] )
      SD <- sqrt(VAR)
      BASE <- C * (SD %o% SD)
      BASE[BI,BI] <- COV.sigma[BI,BI]
      BASE[MI,MI] <- COV.mu[MI,MI]
      COV.mu <- BASE
      COV.sigma <- COV.sigma[-BI,-BI,drop=FALSE]
    }
  }
  else
  { COV.sigma <- diag(0,DIM*(DIM+1)/2) }

  loglike <- -nloglike(par,REML=FALSE) # non-REML for AIC/BIC

  # rescale
  loglike <- loglike - N*DIM*log(prod(SCALE))

  if(BDIM && MDIM) # change to intercept from x-shift (equivalent after scaling)
  { mu[BY] <- mu[BY] - beta * (XSHIFT/XSCALE)[BX] }

  mu <- SCALE[INT]*mu + SHIFT[INT]
  beta <- BSCALE*beta

  CSCALE <- c(BSCALE,SCALE[INT])
  COV.mu <- t(COV.mu*CSCALE)*CSCALE

  sigma <- t(sigma*SCALE)*SCALE

  if(any(VARS))
  {
    CSCALE <- SCALE %o% SCALE
    if(isotropic)
    { CSCALE <- CSCALE[1] }
    else
    { CSCALE <- CSCALE[DUP] }

    COV.sigma <- t(COV.sigma*CSCALE)*CSCALE
    # HESS <- t(HESS/SCALE)/SCALE
  }

  # AIC
  n <- N
  q <- DIM
  qn <- sum(OBS)
  qk <- MDIM + BDIM
  nu <- VDIM
  K <- qk + nu

  AIC <- 2*K - 2*loglike
  AICc <- (qn-qk)/max(qn-K-nu,0)*2*K - 2*loglike
  AICc <- nant(AICc,Inf)
  BIC <- K*log(N) - 2*loglike

  names(mu) <- NAMES[INT] # only non-zero intercepts
  if(BDIM) { names(beta) <- paste0(NAMES[BY],"/",XNAMES[BX]) }
  CNAMES <- c(names(beta),names(mu))
  dimnames(COV.mu) <- list(CNAMES,CNAMES)
  dimnames(sigma) <- list(NAMES,NAMES) # all variances, including zero-variances
  dimnames(VARS) <- list(NAMES,NAMES)

  # full matrix to copy into - all possible matrix elements
  NAMES2 <- outer(NAMES,NAMES,function(x,y){paste0(x,"-",y)})
  COV.SY <- diag(0,DIM^2)
  dimnames(COV.SY) <- list(NAMES2,NAMES2)
  # copy over non-zero elements
  COV.SY[which(DUP),which(DUP)] <- COV.sigma
  # unique matrix elements only
  COV.sigma <- COV.SY[which(TRI),which(TRI)]
  # dimnames(HESS) <- list(NAMES2,NAMES2)

  R <- list(mu=mu,beta=beta,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike,AIC=AIC,AICc=AICc,BIC=BIC,isotropic=isotropic)
  R$VARS <- VARS # need to pass this if variance was turned off due to lack of data
  R$INT <- INT
  # R$HESS <- HESS
  return(R)
}


cross.terms <- function(TERMS,unique=TRUE)
{
  TERMS <- outer(TERMS,TERMS,function(x,y){paste0(x,"-",y)})
  if(unique) { TERMS <- TERMS[upper.tri(TERMS,diag=TRUE)] }
  return(TERMS)
}
