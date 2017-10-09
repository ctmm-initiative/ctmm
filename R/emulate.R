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
    COV[,P] <- COV[,P]/par[P]

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


# fast here debiases for fast emulate
ctmm.boot <- function(data,CTMM,fast=FALSE,n=1E4,mc.cores=NULL,...)
{
  if(is.null(mc.cores)) { mc.cores <- detectCores() } # Windows safe wrapper

  NAMES <- dimnames(CTMM$COV)[[1]] # parameters to extract & debias

  Replicate <- function()
  {
    FIT <- emulate(CTMM,data=data,fast=FALSE,...)
    FIT <- get.parameters(FIT,NAMES)
    return(FIT)
  }

  P0 <- get.parameters(CTMM,NAMES)
  PS <- mclapply(1:n,Replicate,mc.cores=mc.cores)



}
