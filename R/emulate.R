# simulate an estimated model from the parameter estimate covariances
emulate.ctmm <- function(object,data=NULL,fast=FALSE,...)
{
  if(!fast && is.null(data)) { stop("fast=FALSE requires data.") }

  CTMM <- object
  if(!fast) { return( emulate.telemetry(data,CTMM,fast=fast,...) ) }

  AXES <- length(CTMM$axes)

  Sigma <- CTMM$COV.mu
  DIM <- dim(Sigma)
  if(length(DIM)==AXES) # stationary mean
  { if(CTMM$range) { CTMM$mu <- MASS::mvrnorm(mu=CTMM$mu,Sigma=Sigma) } }
  else if(length(DIM)==4)
  {
    mu <- CTMM$mu # (k,x)

    # dont sample stationary mean for BM/IOU
    if(!CTMM$range)
    {
      mu <- mu[-1,,drop=FALSE] # (k-1,x)
      Sigma <- Sigma[,-1,-1,,drop=FALSE] # (x,k,k,x) -> (x,k-1,k-1,x)
      DIM <- dim(Sigma)
    }

    if(DIM[3]) # make sure there's something left to sample
    {
      # flatten covariance block-matrix
      Sigma <- aperm(Sigma,c(2,1,3,4)) # (x,k,k,x) -> (k,x,k,x)
      Sigma <- array(Sigma,c(prod(DIM[2:1]),prod(DIM[3:4]))) # (k,x,k,x) -> (k*x,k*x)
      # flatten mean block-vector
      mu <- array(mu,prod(DIM[1:2])) # k,x -> k*x
      mu <- as.numeric(mu) # aRg!
      mu <- MASS::mvrnorm(mu=mu,Sigma=Sigma) # (k*x,k*x) %*% k*x -> k*x
      mu <- array(mu,DIM[2:1]) # k*x -> k,x

      if(CTMM$range)
      { CTMM$mu <- mu }
      else
      { CTMM$mu[-1] <- mu }
    }
  }
  else stop("COV.mu has odd dimension.")

  COV <- CTMM$COV
  NAMES <- dimnames(COV)[[1]]

  par <- get.parameters(CTMM,NAMES)
  # positive variables (potential)
  PAR <- POSITIVE.PARAMETERS
  # included variables that are positive
  PAR <- PAR[PAR %in% NAMES]
  # log transform positive parameters
  for(P in PAR)
  {
    COV[P,] <- COV[P,]/par[P]
    COV[,P] <- t(t(COV[,P])/par[P])

    if(par[P]==0) { stop("fast=TRUE (CLT) not possible when ",P," = ",par[P]) }
    par[P] <- log(par[P])
  }

  # variances should be comparable after transformation --- better condition number
  TEST <- any( diag(COV) <= .Machine$double.eps ) # first test
  TEST <- TEST || any( eigen(COV)$values <= .Machine$double.eps ) # only run eigen() if first passes
  if(any(TEST)) { stop("Hessian not negative semidefinite. Try fast=FALSE.") }

  par <- MASS::mvrnorm(mu=par,Sigma=COV)

  # transform log back to positive parameters
  par[PAR] <- exp(par[PAR])
  # keep tau sorted
  PAR <- c('tau position','tau velocity')
  if(all(PAR %in% NAMES) && par['tau velocity']>=par['tau position'])
  { par['tau velocity'] <- par['tau position']*(1-.Machine$double.eps) }

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


# remove zero parameters
ctmm.reduce <- function(CTMM)
{
  # parameters to extract & debias
  NAMES <- dimnames(CTMM$COV)[[1]]

  # delete zero parameters
  if(length(NAMES))
  {
    P <- get.parameters(CTMM,NAMES)
    ZERO <- !P & (NAMES %in% POSITIVE.PARAMETERS)
    names(ZERO) <- NAMES

    # remove angle if removed major and minor axes (!isotropic)
    if("minor" %in% NAMES && ZERO['major'] && ZERO['minor']) { ZERO['angle'] <- TRUE }
    # restructure tau if contains zeros
    P <- "tau velocity"; if(P %in% NAMES && ZERO[P]) { CTMM$tau <- CTMM$tau[1] }
    P <- "tau position"; if(P %in% NAMES && ZERO[P]) { CTMM$tau <- NULL }
    P <- "tau"; if(P %in% NAMES && ZERO[P]) { CTMM$tau <- NULL }

    if(any(ZERO)) { CTMM$features <- CTMM$features[ CTMM$features %nin% NAMES[ZERO] ] }
    NAMES <- NAMES[!ZERO]
    CTMM$COV <- CTMM$COV[!ZERO,!ZERO]
  }

  return(CTMM)
}


#####################################
# multi-estimator parametric bootstrap
# + concurrent double-bootstrap AICc
ctmm.boot <- function(data,CTMM,method=CTMM$method,AICc=FALSE,iterate=FALSE,robust=FALSE,error=0.01,clamp=0.001,cores=1,trace=TRUE,...)
{
  if("COV" %nin% names(CTMM)) { stop("CTMM needs to be output from ctmm.select or ctmm.fit.") }

  CTMM <- ctmm.reduce(CTMM) # remove zero parameters
  if(!length(CTMM$COV)) { return(CTMM) } # nothing to do

  # parameters to extract & debias
  NAMES <- dimnames(CTMM$COV)[[1]]

  # positive variables (potential)
  POS <- POSITIVE.PARAMETERS
  # included variables that are positive
  POS <- POS[POS %in% NAMES]

  prior <- clamp
  if(prior)
  {
    DOF <- CTMM$COV
    P <- get.parameters(CTMM,NAMES)
    if(length(NAMES)>length(POS)) { P[NAMES %nin% POS] <- 1 }
    DOF <- t(DOF/P)/P
  }

  method <- match.arg(method,c('ML','HREML','pREML','pHREML','REML'),several.ok=FALSE)
  cores <- resolveCores(cores,fast=FALSE)

  # toss location information for simulations unconditioned
  FRAME <- data
  FRAME[,CTMM$axes] <- NULL

  Transform <- function(p,inverse=FALSE)
  {
    DIM <- dim(p)
    if(length(DIM)<2)
    {
      NAMES <- names(p)
      p <- cbind(p)
    }

    # transform to parameters of interest (to debias)
    if("minor" %in% NAMES)
    {
      # linearize covariance matrix for debiasing
      P <- c("major","minor","angle")
      if(!inverse)
      { p[P,] <- sapply(1:ncol(p),function(i){sigma.construct(p[P,i])[c(1,2,4)]}) }
      else
      { p[P,] <- sapply(1:ncol(p),function(i){sigma.destruct(matrix(p[P,i][c(1,2,2,3)],2,2))}) }
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
    FIT <- ctmm.fit(DATA,CTMM,method=method,COV=FALSE,...)
    if(AICc)
    {
      # double expectation value
      DATA <- simulate(CTMM,data=FRAME,precompute=precompute)
      # update KLD sum
      AIC <- -2*ctmm.loglike(DATA,FIT,profile=FALSE)
    }
    else
    { AIC <- 0 }
    FIT <- get.parameters(FIT,NAMES)
    names(FIT) <- NAMES
    FIT <- Transform(FIT)
    attr(FIT,"AIC") <- AIC # attach to slot
    return(FIT)
  }

  # add batch of calculations to intermediate results
  # rolling mean and variance
  rolling_update <- function(ADD) # [par,run]
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

      if(any(abs(AVE)==Inf)) { stop('Some scale parameter estimates on boundary. Requires robust=TRUE.') }
    }
    else if(N>=length(NAMES)+1)
    {
      STUFF <- rcov(SAMPLE)
      AVE <<- STUFF$median
      COV <<- STUFF$COV

      if(any(abs(AVE)==Inf) || any(abs(diag(COV))==Inf)) { stop('Too many scale parameter estimates on boundary.') }
    }
  }

  # estimate truth from original estimate, correcting for estimated bias
  estimator <- function(STEP=1)
  {
    P <- EST - BIAS
    names(P) <- NAMES
    # fixes boundary crossings, etc.
    P <- Transform(P,inverse=TRUE)
    P <- clean.parameters(P)
    P <- Transform(P)

    return(P)
  }

  ### outer loop: designing the weight matrix until estimate converges
  P0 <- P1 <- EST <- Transform(get.parameters(CTMM,NAMES)) # null model for parametric bootstrap
  #if(trace>1) { print(Transform(P0,inverse=TRUE)) }
  PROGRESS <- TRUE
  ERROR <- Inf -> ERROR.OLD
  BEST <- list()
  tol <- error*sqrt(length(NAMES)) # eigen-value tolerance for pseudo-inverse
  MAX <- max(1/error^2,length(NAMES)+1) # maximum sample size
  if(trace) { pb <- utils::txtProgressBar(style=3) }
  while(ERROR>error && PROGRESS)
  {
    # base progress
    PG <- sqrt(error/ERROR.OLD)

    # initialize results
    SAMPLE <- NULL # sampling distribution sample (transformed)
    AVE <- array(0,length(NAMES)) # average of EST
    COV <-  array(0,c(length(NAMES),length(NAMES))) # covariation of EST
    if(!robust) # more efficient computation with moment accumulation
    {
      S1 <- AVE
      S2 <- COV
    }
    ERROR <- Inf
    N <- 0

    # INITIAL simulation precomputes Kalman filter matrix
    precompute <- TRUE
    ADD <- Replicate() # [par] (transformed)
    AIC <- attr(ADD,"AIC") # to store ensemble of -2*KLD
    dim(ADD) <- c(length(NAMES),1) # [par,1]
    rolling_update(ADD)

    # inner loop: adding to ensemble until average converges
    precompute <- -1 # now use precomputed Kalman filter matrices
    ERROR <- Inf
    while(N<=MAX && ERROR>=error) # standard error ratio criterias
    {
      # calculate burst of cores estimates
      ADD <- plapply(1:cores,Replicate,cores=cores,fast=FALSE)
      AIC <- c(AIC, sapply(ADD, function(add){attr(add,"AIC")} ) )
      ADD <- simplify2array(ADD) # [par,run]
      dim(ADD) <- c(length(NAMES),cores) # [par,run]
      rolling_update(ADD)

      # alternative test if bias is relatively well resolved
      if(N>=length(NAMES)+1)
      {
        BIAS <- AVE - P0 # fixed for multiple iterations
        # ERROR <- sqrt(sum(BIAS^2/diag(COV))/length(NAMES)/N)
        ERROR <- sqrt(c(BIAS %*% PDsolve(COV) %*% BIAS)/length(NAMES)/N)
        if(is.nan(ERROR)) { ERROR <- Inf }
      }

      if(trace) { utils::setTxtProgressBar(pb,PG+(1-PG)*max(N/MAX,sqrt(error/ERROR))) }
    } # END INNER LOOP
    VAR.AIC <- stats::var(AIC)
    AIC <- mean(AIC)

    ## recalculate error & store best model
    ERROR <- AVE - EST # how far off were we?
    ERROR <- sqrt(abs(c(ERROR %*% PDsolve(COV) %*% ERROR))/length(NAMES))
    if(is.nan(ERROR)) { ERROR <- Inf }

    # check to make sure relative error is decreasing
    names(BIAS) <- NAMES
    PROGRESS <- ERROR<ERROR.OLD && EST['major']>BIAS['major']
    # also make sure variance doesn't collapse
    if(PROGRESS) # made progress
    {
      ERROR -> ERROR.OLD

      BEST$COV <- COV
      BEST$BIAS <- BIAS
      BEST$AIC <- AIC
      BEST$VAR.AIC <- VAR.AIC

      # final iteration --- calculate final parameter estimate COV
      if(ERROR<=error || !iterate)
      {
        ## calculate COV (somewhat robustly in untransformed basis)
        # calculate all weighted estimates in sample
        SAMPLE <- SAMPLE - c(BIAS) # debiased estimates
        dim(SAMPLE) <- c(length(NAMES),N)
        rownames(SAMPLE) <- NAMES
        SAMPLE <- Transform(SAMPLE,inverse=TRUE)

        # robust covariance estimate
        # if(!robust) { COV <- cov(t(SAMPLE)) } else { COV <- rcov(SAMPLE) }
        COV <- stats::cov(t(SAMPLE))
        dimnames(COV) <- list(NAMES,NAMES)
        CTMM$COV <- COV
      }
    }
    else # back up to previous and quit
    {
      # failed on first round
      if(ERROR.OLD==Inf) { return(CTMM) }

      # go back to previous (best) result
      BEST$COV -> COV
      BEST$BIAS -> BIAS
      BEST$AIC -> AIC
      BEST$VAR.AIC -> VAR.AIC
    }

    ## calculate best estimate
    P1 <- estimator()

    P <- Transform(P1,inverse=TRUE)
    CTMM <- set.parameters(CTMM,P)
    #if(trace) { message("\nRelative error = ",ERROR) }
    P1 -> P0

    #if(trace>1) { print(P) }
  } # END OUTER LOOP

  # fix COV[mean] with a verbose likelihood evaluation
  CTMM <- ctmm.loglike(data,CTMM,profile=FALSE,verbose=TRUE)

  # fix AIC/BIC/MSPE
  CTMM <- ic.ctmm(CTMM,nrow(data))

  # store simulated AICc
  if(AICc) { CTMM$AICc <- AIC; CTMM$VAR.AICc <- VAR.AIC }

  # crude fix to CIs at small N
  if(prior)
  {
    DOF2 <- CTMM$COV
    P <- get.parameters(CTMM,NAMES)
    if(length(NAMES)>length(POS)) { P[NAMES %nin% POS] <- 1 } #!!! invalid argument type error
    DOF2 <- t(DOF2/P)/P

    R <- clamp*sqrt(N) # <= clamp/error

    DOF <- ( DOF + R*DOF2 )/( 1 + R )
    CTMM$COV <- t(DOF*P)*P
  }

  if(trace) { close(pb) }
  return(CTMM)
}
