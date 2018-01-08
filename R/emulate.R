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
# multi-estimator parametric bootstrap
# + concurrent double-bootstrap AICc
ctmm.boot <- function(data,CTMM,method=c('HREML','pREML','pHREML'),error=0.01,mc.cores=NULL,AICc=FALSE,...)
{
  if(is.null(mc.cores)) { mc.cores <- detectCores() } # Windows safe wrapper

  # toss location information for simulations unconditioned
  FRAME <- data
  FRAME[,CTMM$axes] <- NULL

  # parameters to extract & debias
  NAMES <- dimnames(CTMM$COV)[[1]]
  DIM <- length(NAMES) * length(method)

  # positive variables (potential)
  POS <- c("area","tau position","tau velocity","error")
  # included variables that are positive
  POS <- POS[POS %in% NAMES]

  # log transform positive parameters
  Transform <- function(p,inverse=FALSE)
  {
    if(length(POS))
    {
      if(inverse) { p[POS,] <- exp(p[POS,]) }
      else { p[POS,] <- log(p[POS,]) }
    }

    return(p)
  }

  # simulate with errors & return parameter estimates
  Replicate <- function(i,DATA=simulate(CTMM,data=FRAME))
  {
    # apply all estimators, recycling previous fit as first guess
    FIT <- list(CTMM=CTMM) # first guess (to be removed)
    for(m in method) { FIT <- c(FIT,ctmm.fit(FRAME,FIT[[length(FIT)]],method=m,COV=FALSE,...)) }
    FIT <- FIT[-1] # remove first (fake) guess

    # extract parameters
    FIT <- sapply(FIT,function(f){ get.parameters(f,NAMES) }) # [parameter,method]
    rownames(FIT) <- NAMES
    colnames(FIT) <- method

    # log transform positive parameters
    FIT <- Transform(FIT)

    return(FIT)
  }

  # outer loop: designing the weight matrix until estimate converges
  EST <- Replicate(0,DATA=data) # estimates to weight (ultimately)
  P0 <- Transform(get.parameters(CTMM,NAMES)) # null model for parametric bootstrap
  ERROR <- Inf
  while(ERROR>=error)
  {
    MEAN <- S1 <- array(0,DIM)
    COV <- S2 <- array(0,c(DIM,DIM))
    ERROR <- Inf
    N <- 0

    # inner loop: adding to ensemble until average converges
    while(ERROR>=error)
    {
      # calculate burst of mc.cores estimates
      ADD <- mclapply(1:mc.cores,Replicate,mc.cores=mc.cores)
      ADD <- simplify2array(ADD) # [par,method,run]
      ADD <- array(ADD,c(DIM,mc.cores)) # [par*method,run]

      # rolling mean and variance
      # update sums
      S1 <- S1 + rowSums(ADD)
      S2 <- S2 + (ADD %.% t(ADD))
      # update sample size
      N <- N + mc.cores
      # normalize sums
      MEAN <- S1/N
      COV <- (S2-N^2*outer(MEAN))/max(1,N-1) # !!! subtract mean^2

      # recalculate error
      if(N>10) # may want to increase to 100
      {
        BIAS <- MEAN - P0
        ERROR <- sqrt((BIAS %*% PDsolve(COV) %*% BIAS)/N)
      }
    }

    ## recalculate weighted estimator
    P1 <- EST - BIAS # debiased estimates
    M <- expm::sqrtm(PDsolve(COV))
    P1 <- M %*% P1 # standardized esitmates
    M <- M %*% t(array(diag(length(NAMES)),c(length(NAMES),DIM))) # LSM
    P1 <- t(M) %*% P1
    M <- t(M) %*% M # precision matrix
    COV <- PDsolve(M)
    P1 <- M %*% P1

    ## recalculate error & store best model
    BIAS <- P1 - P0 # now the change
    ERROR <- sqrt(BIAS %*% M %*% BIAS)
    CTMM <- set.parameters(CTMM,Transform(P1,inverse=TRUE))
    P1 -> P0
  }

  # transform COV
  if(length(POS))
  {
    COV[POS,] <- COV[POS,]/P0[POS]
    COV[,POS] <- COV[,POS]/P0[POS]
  }
  CTMM$COV <- COV

  return(CTMM)
}
