# simulate an estimated model from the parameter estimate covariances
emulate.ctmm <- function(object,data=NULL,fast=FALSE,...)
{
  if(!fast && is.null(data)) { stop("fast=FALSE requires data.") }

  CTMM <- object
  if(!fast) { return( emulate.telemetry(data,CTMM,fast=fast,...) ) }

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
  PAR <- c("area","tau position","tau velocity","tau","error")
  # included variables that are positive
  PAR <- PAR[PAR %in% NAMES]
  # log transform positive parameters
  for(P in PAR)
  {
    COV[P,] <- COV[P,]/par[P]
    COV[,P] <- t(t(COV[,P])/par[P])

    par[P] <- log(par[P])
  }

  # variances should be comparable after transformation --- better condition number
  TEST <- c(diag(COV),eigen(COV,only.values=TRUE)$values) <= .Machine$double.eps
  if(any(TEST)) { stop("Hessian not negative semidefinite. Try fast=FALSE.") }

  par <- MASS::mvrnorm(mu=par,Sigma=COV)

  # transform log back to positive parameters
  par[PAR] <- exp(par[PAR])
  # keep tau sorted
  # PAR <- c('tau position','tau velocity')
  # if(all(PAR %in% NAMES)) { par[PAR] <- sort(par[PAR],decreasing=TRUE) }

  CTMM <- set.parameters(CTMM,par)

  return(CTMM)
}


# simulate an estimated model by simulating data and then fitting a model to that data
emulate.telemetry <- function(object,CTMM,fast=FALSE,...)
{
  if(fast) { return( emulate.ctmm(CTMM,data=object,fast=fast,...) ) }

  # copy over time and error structure
  FRAME <- object
  FRAME[,CTMM$axes] <- NULL # delete location information

  # simulate data of the same model and sampling
  object[,CTMM$axes] <- simulate(CTMM,data=FRAME)[,CTMM$axes]
  # error columns are copied over this way

  # fit model to simulated data
  CTMM <- ctmm.fit(object,CTMM,method=CTMM$method,COV=FALSE,...)
  # would never need COV matrix

  return(CTMM)
}


#####################################
# multi-estimator parametric bootstrap
# + concurrent double-bootstrap AICc
ctmm.boot <- function(data,CTMM,method=CTMM$method,multiplicative=TRUE,robust=FALSE,error=0.01,cores=1,trace=TRUE,...)
{
  method <- match.arg(method,c('ML','HREML','pREML','pHREML','REML'),several.ok=TRUE)
  cores <- resolveCores(cores,fast=FALSE)

  # toss location information for simulations unconditioned
  FRAME <- data
  FRAME[,CTMM$axes] <- NULL

  # parameters to extract & debias
  NAMES <- dimnames(CTMM$COV)[[1]]
  DIM <- length(NAMES) * length(method)

  # positive variables (potential)
  POS <- c("area","eccentricity","tau position","tau velocity","tau","error")
  # included variables that are positive
  POS <- POS[POS %in% NAMES]

  # eccentricity
  if("eccentricity" %in% NAMES)
  {
    # ensure eccentricity is positive
    if(attr(CTMM$sigma,'par')['eccentricity']<0)
    {
      attr(CTMM$sigma,'par')['eccentricity'] <- abs( attr(CTMM$sigma,'par')['eccentricity'] )
      attr(CTMM$sigma,'par')['angle'] <- ( (attr(CTMM$sigma,'par')['angle'] + pi/2) %% pi ) - pi/2
    }
  }

  # log transform positive parameters
  Transform <- function(p,inverse=FALSE)
  {
    DIM <- dim(p)
    if(length(DIM)<2)
    {
      NAMES <- names(p)
      p <- cbind(p)
    }

    # which slot will contain the variance
    if("eccentricity" %in% NAMES) { VAR <- "eccentricity" } else { VAR <- "area" }

    # transform to parameters of interest
    if(!inverse)
    {
      # transform eccentricity to variance (not area)
      if("eccentricity" %in% NAMES) { p["eccentricity",] <- p['area',]*cosh(p["eccentricity",]/2) }

      # if(CTMM$omega) # underdamped
      # {
      #   # mean square velocity and average diffusion
      #   p["omega",] <- p[VAR,]*(1/p["tau",]^2+p["omega",]^2) # MSV
      #   p["tau",] <- p[VAR,]/p["tau",] # DIF
      #   # refer back to get.taus() for how this was parameterized
      # }
      # else if("tau" %in% NAMES) # critically damped
      # {
      #   p["tau",] <- p[VAR,]/p["tau",]^2 # mean square Velocity
      # }
      # else # overdamped
      # {
      #   # transform tau position to diffusion rate
      #   if("tau position" %in% NAMES)
      #   { DIF <- p[VAR,]/p["tau position",] }
      #
      #   # mean square velocity
      #   if("tau velocity" %in% NAMES)
      #   {
      #     if(!CTMM$range) { MSV <- p[VAR,]/(p["tau velocity",]) }
      #     else { MSV <- p[VAR,]/(p["tau position",]*p["tau velocity",]) }
      #   }
      #
      #   if("tau velocity" %in% NAMES) { p["tau velocity",] <- MSV }
      #   if("tau position" %in% NAMES) { p["tau position",] <- DIF }
      # }

      # log transform
      if(multiplicative && length(POS)) { p[POS,] <- log(p[POS,]) }
    }

    # transform back from parameters of interest
    if(inverse)
    {
      if(multiplicative && length(POS)) { p[POS,] <- exp(p[POS,]) }

      # if(CTMM$omega) # underdamped
      # {
      #   p["tau",] <- p[VAR,]/p["tau",]
      #   p["omega",] <- p["omega",]/p[VAR,] - 1/p["tau",]^2 # omega^2
      #   p["omega",] <- pmax(p["omega",],0) # bias should be other way
      #   p["omega",] <- sqrt(p["omega",])
      # }
      # else if("tau" %in% NAMES) # critically damped
      # {
      #   p["tau",] <- sqrt(p[VAR,]/p["tau",])
      # }
      # else # overdamped
      # {
      #   if("tau velocity" %in% NAMES) { p["tau velocity",] -> MSV }
      #   if("tau position" %in% NAMES) { p["tau position",] -> DIF }
      #
      #   if("tau position" %in% NAMES) { p["tau position",] <- p[VAR,]/DIF }
      #
      #   if("tau velocity" %in% NAMES)
      #   {
      #     if(CTMM$range) { p["tau velocity",] <- p[VAR,]/(MSV*p["tau position",]) }
      #     else { p["tau velocity",] <- p[VAR,]/MSV }
      #   }
      # }

      if("eccentricity" %in% NAMES) { p['eccentricity',] <- 2*acosh(max(p['eccentricity',]/p['area',],1)) }
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
  precompute <- FALSE
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

    FIT <- Transform(FIT)

    return(FIT)
  }

  # add batch of calculations to intermediate results
  # rolling mean and variance
  rolling_update <- function(ADD) # [par*method,run]
  {
    N <<- N + ncol(ADD) # update sample size
    SAMPLE <<- cbind(SAMPLE,ADD) # update sample (transformed)
    if(!robust)
    {
      # update sums
      S1 <<- S1 + rowSums(ADD) # sum over runs
      S2 <<- S2 + (ADD %.% t(ADD)) # inner product over runs

      # normalize sums
      AVE <<- S1/N
      COV <<- (S2-N*outer(AVE))/max(1,N-1)

      if(any(abs(AVE)==Inf)) { stop('Some scale parameter estimates on boundary. Requires multiplicative=FALSE or robust=TRUE.') }
    }
    else if(N>=DIM+1)
    {
      STUFF <- rcov(SAMPLE)
      AVE <<- STUFF$median
      COV <<- STUFF$COV

      if(any(abs(AVE)==Inf) || any(abs(diag(COV))==Inf)) { stop('Too many scale parameter estimates on boundary. Requires multiplicative=FALSE.') }
    }
  }

  W <- array(1/length(method),c(length(NAMES),length(method))) # weights
  J <- array(1,c(length(method),1)) # de-projection matrix

  ### outer loop: designing the weight matrix until estimate converges
  EST <- Replicate(0,DATA=data) # original estimates that we are working with to optimally weight (transformed)
  dim(EST) <- DIM # (NAMES*methods)
  P0 <- Transform(get.parameters(CTMM,NAMES)) -> P1 # null model for parametric bootstrap
  if(trace>1) { print(P0) }
  ERROR <- Inf -> ERROR.OLD
  tol <- error*sqrt(length(NAMES)) # eigen-value tolerance for pseudo-inverse
  MAX <- max(1/error^2,DIM+1) # maximum sample size
  if(trace) { pb <- utils::txtProgressBar(style=3) }
  while(ERROR>error)
  {
    # initialize results
    SAMPLE <- NULL # sampling distribution sample (transformed)
    AVE <- array(0,DIM) # average of EST
    COV <-  array(0,c(DIM,DIM)) # covariation of EST
    if(!robust) # more efficient computation with moment accumulation
    {
      S1 <- AVE
      S2 <- COV
    }
    ERROR <- Inf
    N <- 0

    # INITIAL simulation precomputes Kalman filter matrix
    precompute <- TRUE
    ADD <- Replicate() # [par,method] (transformed)
    dim(ADD) <- c(DIM,1) # [par*method,1]
    rolling_update(ADD)

    # inner loop: adding to ensemble until average converges
    precompute <- -1 # now use precomputed Kalman filter matrices
    while(N<=MAX && ERROR>=error) # standard error ratio criterias
    {
      # calculate burst of cores estimates
      ADD <- plapply(1:cores,Replicate,cores=cores,fast=FALSE)
      ADD <- simplify2array(ADD) # [par,method,run]
      dim(ADD) <- c(DIM,cores) # [par*method,run]
      rolling_update(ADD)

      # alternative test if bias is relatively well resolved
      if(N>=DIM+1)
      {
        BIAS <- AVE - EST
        ERROR <- sqrt(sum(BIAS^2/diag(COV))/N/DIM)
        if(is.nan(ERROR)) { ERROR <- Inf }
      }

      if(trace) { utils::setTxtProgressBar(pb,clamp(max(N/MAX,(error/ERROR)^2))) }
    } # END INNER LOOP

    ## recalculate weighted estimator
    ## !!! add optimal method to use cross-parameter estimate covariance
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
    ERROR <- P1 - P0 # how much does best estimate change?
    ERROR <- sqrt(abs(c(ERROR %*% PDsolve(COV.W) %*% ERROR))/length(NAMES)/2) # 2-point comparison correction
    if(is.nan(ERROR)) { ERROR <- Inf }

    # check to make sure relative error is decreasing
    if(ERROR>ERROR.OLD)
    {
      # maybe allow this once or twice with a count???
      if(trace>1) { warning("Iterative convergence ceased at error=",ERROR.OLD," over tolerance of ",error)  }
      break
    }
    ERROR -> ERROR.OLD

    # order of current error, which we do not know
    ERROR <- ERROR^2
    names(P1) <- NAMES
    if(trace>1) { print(P1) }
    P <- Transform(P1,inverse=TRUE)
    CTMM <- set.parameters(CTMM,P)
    P1 -> P0
  } # END OUTER LOOP

  ## calculate COV (somewhat robustly in untransformed basis)
  # calculate all weighted esitmates in sample
  SAMPLE <- SAMPLE - c(BIAS) # debiased estimates
  dim(SAMPLE) <- c(length(NAMES),length(method),N)
  SAMPLE <- sapply(1:length(NAMES),function(p) W[p,] %*% array(SAMPLE[p,,],c(length(method),N)) ) # R arrays are insane
  dim(SAMPLE) <- c(N,length(NAMES)) # R arrays are insane, part 2
  SAMPLE <- t(SAMPLE)
  # transform weighted estimates to canonical coordinates
  rownames(SAMPLE) <- NAMES
  SAMPLE <- Transform(SAMPLE,inverse=TRUE)

  # robust covariance estimate
  # if(!robust) { COV <- cov(t(SAMPLE)) } else { COV <- rcov(SAMPLE) }
  COV <- stats::cov(t(SAMPLE))
  dimnames(COV) <- list(NAMES,NAMES)
  CTMM$COV <- COV

  # we lost eccentricity
  if("eccentricity" %in% NAMES && P["eccentricity"]==0) { CTMM <- simplify.ctmm(CTMM,"eccentricity") }
  # what other parameters could have been lost?

  # fix COV[mean] with a verbose likelihood evaluation
  CTMM <- ctmm.loglike(data,CTMM,profile=FALSE,verbose=TRUE)

  if(trace) { close(pb) }
  return(CTMM)
}
