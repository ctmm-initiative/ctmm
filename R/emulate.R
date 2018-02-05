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
ctmm.boot <- function(data,CTMM,method=c('HREML','pREML','pHREML'),error=0.01,mc.cores=1,trace=TRUE,...)
{
  method <- match.arg(method,c('ML','HREML','pREML','pHREML','REML'),several.ok=TRUE)
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
    DIM <- dim(p)
    if(length(DIM)<2)
    {
      NAMES <- names(p)
      p <- cbind(p)
    }

    if(length(POS))
    {
      if(inverse) { p[POS,] <- exp(p[POS,]) }
      else { p[POS,] <- log(p[POS,]) }
    }

    # restore dimensionality and array names
    if(length(DIM)<2)
    {
      dim(p) <- DIM
      names(p) <- NAMES
    }

    return(p)
  }


  # simulate with errors & return parameter estimates
  Replicate <- function(i=0,DATA=simulate(CTMM,data=FRAME,precompute=precompute),...)
  {
    # apply all estimators, recycling previous fit as first guess
    FIT <- list(CTMM=CTMM) # first guess (to be removed)
    for(m in method) { FIT[[length(FIT)+1]] <- ctmm.fit(DATA,FIT[[length(FIT)]],method=m,COV=FALSE,...) }
    FIT <- FIT[-1] # remove first (fake) guess

    # extract parameters
    FIT <- sapply(FIT,function(f){ get.parameters(f,NAMES) }) # [parameter,method]
    rownames(FIT) <- NAMES
    colnames(FIT) <- method

    # log transform positive parameters
    FIT <- Transform(FIT)

    return(FIT)
  }

  # add batch of calculations to intermediate results
  # rolling mean and variance
  rolling_update <- function(ADD) # [par*method,run]
  {
    # update sums
    S1 <<- S1 + rowSums(ADD) # sum over runs
    S2 <<- S2 + (ADD %.% t(ADD)) # inner product over runs
    # update sample size
    N <<- N + mc.cores
    # normalize sums
    MEAN <<- S1/N
    COV <<- (S2-N*outer(MEAN))/max(1,N-1)
  }

  W <- array(1/length(method),c(length(NAMES),length(method))) # weights
  J <- array(1,c(length(method),1)) # de-projection matrix

  ### outer loop: designing the weight matrix until estimate converges
  EST <- Replicate(0,DATA=data) # fundamental estimates to weight (ultimately)
  dim(EST) <- DIM
  P0 <- Transform(get.parameters(CTMM,NAMES)) -> P1 # null model for parametric bootstrap
  if(trace>1) { print(P0) }
  ERROR <- Inf -> ERROR.OLD
  tol <- error*sqrt(length(NAMES)) # eigen-value tolerance for pseudo-inverse
  if(trace) { pb <- utils::txtProgressBar(style=3) }
  while(ERROR>error)
  {
    # initialize results
    MEAN <- S1 <- array(0,DIM)
    COV <- S2 <- array(0,c(DIM,DIM))
    ERROR <- Inf
    N <- 0

    # INITIAL simulation precomputes Kalman filter matrix
    precompute <- TRUE
    ADD <- Replicate() # [par,method]
    dim(ADD) <- c(DIM,1) # [par*method,1]
    rolling_update(ADD)

    # inner loop: adding to ensemble until average converges
    precompute <- -1 # now use precomputed Kalman filter matrices
    MAX <- 1/error^2
    while(N<=MAX && ERROR>=error) # standard error ratio criterias
    {
      # calculate burst of mc.cores estimates
      ADD <- mclapply(1:mc.cores,Replicate,mc.cores=mc.cores)
      ADD <- simplify2array(ADD) # [par,method,run]
      dim(ADD) <- c(DIM,mc.cores) # [par*method,run]
      rolling_update(ADD)

      # alternative test if bias is relatively well resolved
      if(N>20)
      {
        BIAS <- MEAN - EST
        ERROR <- sqrt(sum(BIAS^2/diag(COV))/N/DIM)
      }

      if(trace) { utils::setTxtProgressBar(pb,clamp(max(N/MAX,(error/ERROR)^2))) }
    } # END INNER LOOP

    ## recalculate weighted estimator
    JACOB <- array(0,c(length(NAMES),DIM))
    for(i in 1:length(NAMES))
    {
      SUB <- array(FALSE,c(length(NAMES),length(method)))
      SUB[i,] <- TRUE

      H <- PDsolve(COV[SUB,SUB],pseudo=TRUE,tol=tol)
      VAR <- PDsolve(t(J) %*% H %*% J)

      W[i,] <- VAR %*% t(J) %*% H
      P1[i] <- c(W[i,] %*% (EST[SUB]-BIAS[SUB]))

      JACOB[i,] <- array(W*SUB,DIM)
    }
    # covariance of weighted estimator
    COV.W <- JACOB %*% COV %*% t(JACOB)

    ## recalculate error & store best model
    BIAS <- P1 - P0 # how much does best estimate change?
    # previous error estimate
    ERROR <- sqrt(abs(c(BIAS %*% PDsolve(COV.W) %*% BIAS))/length(NAMES)/2) # 2-point comparison correction

    # check to make sure relative error is decreasing
    if(ERROR>ERROR.OLD)
    {
      if(trace>1) { sprintf("Iterative convergence ceased at error=%f over tolerance of %f",ERROR.OLD,error)  }
      break
    }
    COV.W -> COV.OLD
    ERROR -> ERROR.OLD

    # order of current error, which we do not know
    ERROR <- ERROR^2
    names(P1) <- NAMES
    if(trace>1) { print(P1) }
    P <- Transform(P1,inverse=TRUE)
    CTMM <- set.parameters(CTMM,P)
    P1 -> P0
  } # END OUTER LOOP

  COV.W -> COV
  dimnames(COV) <- list(NAMES,NAMES)

  # transform COV
  if(length(POS))
  {
    COV[POS,] <- COV[POS,] * P[POS]
    COV[,POS] <- t(t(COV[,POS]) * P[POS])
  }
  CTMM$COV <- COV

  # fix COV[mean] with one last likelihood evaluation
  CTMM <- ctmm.loglike(data,CTMM,profile=FALSE,verbose=TRUE)

  if(trace) { close(pb) }
  return(CTMM)
}
