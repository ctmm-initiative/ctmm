# simulate an estimated model from the parameter estimate covariances
emulate.ctmm <- function(object,data=NULL,fast=FALSE,...)
{
  CTMM <- object
  if(!fast) { return( emulate(data,CTMM,fast=fast,...) ) }

  # needs to be updated for periodic stuff
  Sigma <- CTMM$COV.mu
  DIM <- dim(Sigma)
  if(length(DIM)==2)
  {  CTMM$mu <- MASS::mvrnorm(mu=CTMM$mu,Sigma=Sigma) }
  else if(length(DIM)==4)
  {
    # flatten covariance block-matrix
    Sigma <- aperm(Sigma,c(2,1,3,4)) # (x,k,k,x) -> (k,x,k,x)
    Sigma <- array(Sigma,c(prod(DIM[2:1]),prod(DIM[3:4]))) # (k,x,k,x) -> (k*x,k*x)
    # flatten mean block-vector
    mu <- CTMM$mu # k,x
    mu <- array(mu,prod(DIM[1:2])) # k,x -> k*x
    mu <- as.numeric(mu) # aRg!
    mu <- MASS::mvrnorm(mu=mu,Sigma=Sigma) # (k*x,k*x) %*% k*x -> k*x
    mu <- array(mu,DIM[2:1]) # k*x -> k,x
    CTMM$mu <- mu
  }
  else stop("COV.mu has odd dimension.")

  COV <- CTMM$COV
  NAMES <- dimnames(COV)[[1]]
  par <- get.parameters(CTMM,NAMES)

  # positive variables (potential)
  PAR <- c("area","tau position","tau velocity","error")
  # included variables that are positive
  PAR <- PAR[PAR %in% NAMES]
  # log transform positive parameters
  for(P in PAR)
  {
    COV[P,] <- COV[P,]/par[P]
    COV[,P] <- t(t(COV[,P])/par[P])

    par[P] <- log(par[P])
  }

  par <- MASS::mvrnorm(mu=par,Sigma=COV)

  # transform log back to positive parameters
  par[PAR] <- exp(par[PAR])

  CTMM <- set.parameters(CTMM,par)

  return(CTMM)
}


# simulate an estimated model by simulating data and then fitting a model to that data
emulate.telemetry <- function(object,CTMM,fast=FALSE,...)
{
  if(fast) { return( emulate(CTMM,data=object,fast=fast,...) ) }

  # copy over time and error structure
  FRAME <- object
  FRAME[,CTMM$axes] <- NULL # delete location information

  # simulate data of the same model and sampling
  data <- simulate(CTMM,data=FRAME)

  # fit model to simulated data
  CTMM <- ctmm.fit(data,CTMM,method=CTMM$method,...)

  return(CTMM)
}


#####################################
# fast here debiases for fast emulate
ctmm.boot <- function(data,CTMM,fast=FALSE,n=1E4,mc.cores=NULL,...)
{
  if(is.null(mc.cores)) { mc.cores <- detectCores() } # Windows safe wrapper

  NAMES <- dimnames(CTMM$COV)[[1]] # parameters to extract & debias

  Replicate <- function(i)
  {
    FIT <- emulate(CTMM,data=data,fast=FALSE,...)
    FIT <- get.parameters(FIT,NAMES)
    return(FIT)
  }

  PARS <- mclapply(1:n,Replicate,mc.cores=mc.cores)
  PARS <- simplify2array(PARS)
  rownames(PARS) <- NAMES

  sigma <- attr(CTMM$sigma,"par")

  # create MSD / average variance data row
  if(CTMM$isotropic)
  { STAT <- PARS['area',] }
  else
  { STAT <- PARS['area',] * cosh(PARS['eccentricity',]/2) }
  PARS <- rbind(PARS,STAT)
  rownames(PARS)[nrow(PARS)] <- 'var'

  ###################################
  ### debias important quantities ###

  # debias p1 given p1 -> p2
  # p1: null model parameter
  # p2 = p1 + bias
  # lower if p>0
  debias <- function(p1,p2,lower=FALSE)
  {
    if(lower && p1<p2) # relative/scale debias formula (same to first order)
    { p0 <- p1^2 / p2 } # remember sign for magnitudes...
    else # absolute/shift debias formula (exact on first order)
    { p0 <- 2*p1 - p2 }

    return(p0)
  }

  ### (SCALE) standard area debias
  STAT <- mean(PARS['area',])
  sigma['area'] <- debias(sigma['area'],STAT,lower=TRUE)

  if(!CTMM$isotropic)
  {
    ### (SCALE) average variance debias
    STAT <- mean( PARS['var',] )
    # cosh factor
    STAT <- STAT/sigma['area']
    # (SCALE) debiased via cosh factor
    STAT <- debias(cosh(sigma['eccentricity']/2),STAT,lower=TRUE)
    if(STAT>1)
    {
      sigma['eccentricity'] <- 2*acosh(STAT)

      # (SHIFT) debiased angle
      sigma['angle'] <- debias(sigma['angle'],mean(PARS['angle',]))
    }
    else # eccentricity is explained entirely by bias
    {
      sigma['eccentricity'] <- 0
      sigma['angle'] <- 0
      CTMM$isotropic <- FALSE
    }
  }

  # set sigma parameters
  CTMM$sigma <- covm(sigma,isotropic=CTMM$isotropic)
  VAR <- sigma['area'] * cosh(sigma['eccentricity']/2)

  ### (SCALE) average diffusion rate debias
  if('tau position' %in% NAMES)
  {
    # average diffusion rate
    STAT <- mean( PARS['var',] / PARS['tau position',] )
    # tau position
    STAT <- VAR/STAT
    # (SCALE) dibias via tau
    CTMM$tau[1] <- debias(CTMM$tau[1],STAT,lower=TRUE)
  }

  ### (SCALE) MS velocity debias
  if('tau velocity' %in% NAMES)
  {
    if(CTMM$range)
    {
      STAT <- mean( PARS['var',] / PARS['tau position',] / PARS['tau velocity',] )
      STAT <- VAR / CTMM$tau[1] / STAT
    }
    else
    {
      STAT <- mean( PARS['var',] / PARS['tau velocity',] )
      STAT <- VAR / STAT
    }
    # debias via tau velocity
    CTMM$tau[2] <- debias(CTMM$tau[2],STAT,lower=TRUE)
  }

  ### (SCALE) MS velocity debias
  if('circle' %in% NAMES)
  {
    SIGN <- sign(CTMM$circle)
    STAT <- mean( PARS['var',] * PARS['circle',]^2 )
    # abs(circle)
    STAT <- sqrt( STAT / VAR ) * SIGN
    CTMM$circle <- debias(CTMM$circle,STAT)
    if(sign(CTMM$circle)!=SIGN) { CTMM$circle <- FALSE } # circle is entirely explained by bias
  }

  ### (SCALE) error
  if('error' %in% NAMES)
  {
    STAT <- mean( PARS['error',] )
    CTMM$error <- debias(CTMM$error,STAT,lower=TRUE)
  }

  ### parameter estimate covariances ###

  P0 <- get.parameters(CTMM,NAMES)
  PARS <- PARS[NAMES,] # discard derived parameters

  # positive variables (potential)
  POS <- c("area","tau position","tau velocity","error")
  # included variables that are positive
  POS <- POS[POS %in% NAMES]
  # log transform positive parameters
  if(length(POS)) { PARS[POS,] <- log(PARS[POS,]) }

  COV <- stats::cov(t(PARS))

  # log transform back
  if(length(POS))
  {
    #P0 <- mean(t(PARS[POS,]))

    COV[POS,] <- P0[POS] * COV[POS,]
    COV[,POS] <- t(P0[POS] * t(COV[,POS]))
  }

  # recalculate trend terms
  CTMM <- ctmm.loglike(data,CTMM,verbose=TRUE,profile=FALSE)

  CTMM$COV <- COV

  return(CTMM)
}
