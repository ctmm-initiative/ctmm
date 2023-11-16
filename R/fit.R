# smallest resolutions in data (for soft bounding parameter guesstimates)
telemetry.mins <- function(data,axes=c('x','y'))
{
  if(class(data)[1]=="phylometry")
  { dt <- c(attr(data,"lag")) }
  else
  { dt <- diff(data$t) }

  dt <- dt[dt>0]

  if(length(dt))
  {
    # smallest frequency
    if(class(data)[1]=="phylometry")
    { df <- 2*pi/max(dt) }
    else
    { df <- 2*pi/(last(data$t)-data$t[1]) }

    dt <- stats::median(dt) # median time difference
  }
  else
  {
    dt <- df <- 1 # no time differences in data
  }

  dz <- get.telemetry(data,axes)
  dz <- apply(dz,2,diff)
  dim(dz) <- c(nrow(data)-1,length(axes)) # R drops length-1 dimensions
  dz <- rowSums( dz^2 )
  dz <- dz[dz>0]
  if(length(dz))
  { dz <- sqrt(min(dz)) } # smallest nonzero distance
  else
  { dz <- 1 } # no distances in data

  return(list(dt=dt,df=df,dz=dz))
}


# return relevant log-likelihood function - avoids extra function calling
get.loglike <- function(data)
{
  CLASS <- class(data)[1]
  if(CLASS=="phylometry") { return(ctpm.loglike) }
  else { return(ctmm.loglike) }
}

###########################################################
# FIT MODEL WITH LIKELIHOOD FUNCTION (convenience wrapper to optim)
ctmm.fit <- function(data,CTMM=ctmm(),method="pHREML",COV=TRUE,control=list(),trace=FALSE)
{
  # check.class(data)

  loglike.fn <- get.loglike(data)

  if(!is.null(control$message)) { message <- control$message }
  # pass trace argument (demoted)
  if(is.null(control$trace) && trace) { control$trace <- trace-1 }

  method <- match.arg(method,c("ML","pREML","pHREML","HREML","REML"))
  axes <- CTMM$axes
  AXES <- length(axes)
  CTMM <- get.mle(CTMM) # if has better start for pREML/HREML/pHREML

  if(is.null(CTMM$sigma))
  {
    K <- length(CTMM$tau)
    if(K==0 && CTMM$range)
    { CTMM$sigma <- covm(stats::cov(get.telemetry(data,axes)),isotropic=CTMM$isotropic,axes=axes) }
    else # above fails for IOU/BM
    {
      CTMM <- ctmm.guess(data,CTMM=CTMM,interactive=FALSE)
      # preserve continuity
      if(K==0 && !CTMM$range) # assume BM
      { CTMM$tau <- Inf }
      else if(K==1) # OU/BM
      { CTMM$tau <- CTMM$tau[1] }
    }
  }

  ### standardize data for numerical stability ###
  DATA <- data[,axes,drop=FALSE]
  # but apply link first
  link <- get.link(CTMM)
  log.like.link <- sum(log(abs(link$grad(DATA))))
  DATA <- link$fn(DATA)
  CTMM$link <- "identity" # does not need to be applied again in loglike function
  # pre-centering the data reduces later numerical error across models (tested)
  SHIFT <- colMeans(DATA)
  data[,axes] <- t( t(DATA) - SHIFT )
  rm(DATA)
  # add mu.center back to the mean value after kalman filter / mean profiling
  # pre-standardizing the data should also help
  SCALE <- sqrt(mean(get.telemetry(data,axes)^2))
  if(SCALE==0) { SCALE <- 1 }
  # standardize time by median diff time
  if(class(data)[1]=="phylometry")
  { DT <- c(attr(data,"lag")) }
  else
  { DT <- diff(data$t) }
  DT <- DT[DT>0]
  if(length(DT)) { TSCALE <- stats::median(DT) } else { TSCALE <- 1 }
  data <- unit.telemetry(data,length=SCALE,time=TSCALE,axes=axes)
  CTMM <- unit.ctmm(CTMM,length=SCALE,time=TSCALE)
  # TSCALE is effectively 1 now, consistent for profiling in ctmm.fit()

  # used for minimum scale of parameter inspection
  n <- nrow(data)
  MINS <- telemetry.mins(data,axes)
  dt <- MINS$dt
  df <- MINS$df
  dz <- MINS$dz

  # unstandardize (includes likelihood adjustment)
  unscale.ctmm <- function(CTMM)
  {
    CTMM <- unit.ctmm(CTMM,length=1/SCALE,time=1/TSCALE)

    # translate back to origin from center
    CTMM <- drift.shift(CTMM,SHIFT)

    # log-likelihood adjustment
    CTMM$loglike <- CTMM$loglike - length(axes)*n*log(SCALE)

    if(any(UERE.FIT)) { CTMM$loglike <- CTMM$loglike - length(axes)*sum(UERE.DOF[UERE.FIT])*log(SCALE) }

    # fix numeric error (from scaling) when it should be logical
    if(any(!UERE.FIT)) { CTMM$error[!UERE.FIT] <- as.logical(CTMM$error[!UERE.FIT]) }

    return(CTMM)
  }

  default <- list(method="pNewton",precision=1/2,maxit=.Machine$integer.max)
  control <- replace(default,names(control),control)
  precision <- control$precision
  op.method <- control$method
  control$method <- NULL

  # clean/validate
  CTMM <- ctmm.ctmm(CTMM)
  CTMM$mu <- NULL # can always profile mu analytically
  range <- CTMM$range

  # save for fitting
  COV.init <- CTMM$COV
  # make sure we can start from previous failed fit
  if(any(is.nan(COV.init) | COV.init==Inf)) { COV.init <- NULL }
  if(!is.null(COV.init)) { TEST <- eigen(COV.init)$values } else { TEST <- FALSE }
  if(any(TEST<=.Machine$double.eps | TEST==Inf)) { COV.init <- NULL }
  # erase previous fitting info if present
  CTMM$COV <- NULL
  CTMM$COV.mu <- NULL
  CTMM$DOF.mu <- NULL

  # evaluate mean function and error matrices for this data once upfront
  CTMM <- ctmm.prepare(data,CTMM,tau=FALSE,calibrate=FALSE) # don't muck with taus, don't calibrate
  TYPE <- DOP.match(axes)
  if(TYPE!="unknown")
  {
    UERE.DOF <- attr(data,"UERE")$DOF[,TYPE]
    names(UERE.DOF) <- rownames(attr(data,"UERE")$DOF)
  }
  else
  {
    UERE.DOF <- 0
    names(UERE.DOF) <- "all"
  }

  UERE.FIT <- CTMM$error>0 & !is.na(UERE.DOF) & UERE.DOF<Inf # will we be fitting any error parameters?
  UERE.FIX <- CTMM$error>0 & (is.na(UERE.DOF) | UERE.DOF==Inf) # are there any fixed error parameters?
  UERE.PAR <- names(UERE.FIT)[UERE.FIT>0] # names of fitted UERE parameters
  # fix numeric error (from scaling) when it should be logical
  if(any(!UERE.FIT)) { CTMM$error[!UERE.FIT] <- as.logical(CTMM$error[!UERE.FIT]) }

  # don't try to fit error class parameters absent from data
  if(any(CTMM$error>0) && "class" %in% names(data) && TYPE!="unknown")
  {
    LEVELS <- levels(data$class)
    UERE.DOF <- UERE.DOF[LEVELS]
    UERE.FIT <- UERE.FIT[LEVELS]
    UERE.FIX <- UERE.FIX[LEVELS]
    UERE.PAR <- UERE.PAR[UERE.PAR %in% LEVELS]
    CTMM$error <- CTMM$error[LEVELS]
  }

  # make sure to include calibration in log-likelihood even if error==FALSE
  CTMM$errors <- any(CTMM$error>0)

  # can we profile any variances?
  profile <- !any(UERE.FIX)
  # doesn't always work to profile with error fitting
  profile <- profile & !any(UERE.FIT)

  ### id and characterize parameters for profiling ###
  pars <- NAMES <- parscale <- lower <- upper <- period <- NULL
  ORIGINAL <- CTMM # original structure of model before fitting
  linear.cov <- FALSE # represent sigma linearly (for perturbation) versus * (for optimization)
  setup.parameters <- function(CTMM,profile=profile,linear=FALSE)
  {
    STUFF <- id.parameters(CTMM,profile=profile,linear=linear,linear.cov=linear.cov,UERE.FIT=UERE.FIT,dt=dt,df=df,dz=dz,STRUCT=ORIGINAL,fit=TRUE)
    NAMES <<- STUFF$NAMES
    parscale <<- STUFF$parscale
    lower <<- STUFF$lower
    upper <<- STUFF$upper
    period <<- STUFF$period
    # reflect <<- STUFF$reflect
    # initial guess for optimization
    pars <<- get.parameters(CTMM,NAMES,profile=profile,linear.cov=linear.cov)
  }
  setup.parameters(CTMM,profile=profile)
  # fix numeric error (from mucking) when it should be logical
  if(any(!UERE.FIT)) { CTMM$error[!UERE.FIT] <- as.logical(CTMM$error[!UERE.FIT]) }

  # degrees of freedom, including the mean, variance/covariance, tau, and error model
  k.mean <- ncol(CTMM$mean.vec)

  # no reason to do pREML with infinite-variance processes
  if(!range && k.mean==1 && method %in% c('pREML','pHREML','HREML')) { method <- "REML" }

  if(method=="REML") { REML <- TRUE }
  else { REML <- FALSE }

  # OPTIMIZATION FUNCTION (fn)
  # optional argument lengths: TAU, TAU+1, TAU+SIGMA
  fn <- function(p,zero=0)
  {
    names(p) <- NAMES

    # catch for zero error
    if(any(CTMM$error[UERE.FIT]>0 & UERE.DOF[UERE.FIT]>0 & p[paste("error",names(CTMM$error[UERE.FIT]))]==0))
    { return(Inf) }
    # otherwise, UERE.DOF is not counted in likelihood

    p <- clean.parameters(p,profile=profile,linear.cov=linear.cov,timelink=CTMM$timelink)
    CTMM <- set.parameters(CTMM,p,profile=profile,linear.cov=linear.cov)

    # negative log likelihood
    return(-loglike.fn(data,CTMM,REML=REML,zero=-zero,profile=profile))
  }

  # shift parameters back near origin to prevent overflow in optimizer
  reset <- function(p) { clean.parameters(p,profile=profile,linear.cov=FALSE,timelink=CTMM$timelink) }

  # construct covariance matrix guess - must account for differences between storage and optimization representations
  # covariance <- function()
  # {
  #   NAMES.init <- rownames(COV.init)
  #
  #   COV <- diag(parscale^2,nrow=length(parscale))
  #   dimnames(COV) <- list(NAMES,NAMES)
  #   COPY <- NAMES.init %in% NAMES
  #   if(any(COPY))
  #   {
  #     COPY <- NAMES.init[COPY]
  #     COV[COPY,COPY] <- COV.init[COPY,COPY]
  #   }
  #
  #   # relative variance transformation
  #   if(profile)
  #   {
  #     SP <- c('major','minor') # major should never be here
  #     SP <- SP[SP %in% COPY]
  #
  #     EP <- paste('error',names(UERE.FIT)[UERE.FIT])
  #     EP <- EP[EP %in% COPY]
  #
  #     MAX <- profiled.var(CTMM,UERE.RMS=CTMM$error[UERE.FIT],AVE=FALSE)
  #     PARS <- c(SP,EP)
  #     COV[PARS,] <- COV[PARS,]/MAX
  #     COV[,PARS] <- COV[,PARS]/MAX
  #   }
  #
  #   return(COV)
  # }

  # is the model under constrained?
  # UNDER <- AXES*k.mean + length(id.parameters(CTMM,profile=FALSE,UERE.FIT=UERE.FIT,STRUCT=ORIGINAL)$NAMES) > AXES*nrow(data)

  ### NOW OPTIMIZE ###
  # profile <- !any(UERE.FIX)
  # if(UNDER) # under constrained model
  # {
  #   CTMM$loglike <- -Inf
  #   CTMM$AIC <- Inf
  #   CTMM$AICc <- Inf
  #   CTMM$BIC <- Inf
  #   CTMM$MSPE <- c(Inf,Inf)
  #   names(CTMM$MSPE) <- c("position","velocity")
  # }
  if(length(NAMES)==0) # EXACT
  {
    if(method %in% c("pHREML","HREML")) { REML <- TRUE } # IID limit pHREML/HREML -> REML

    # Bi-variate Gaussian with zero error
    CTMM <- loglike.fn(data,CTMM=CTMM,REML=REML,verbose=TRUE)

    # pREML perturbation adjustment
    if(method=="pREML")
    {
      VAR.MULT <- (1+k.mean/n)
      CTMM$sigma <- scale.covm(CTMM$sigma,VAR.MULT)
      CTMM$COV.mu <- VAR.MULT * CTMM$COV.mu
    }

    if(method=="pREML") { REML <- TRUE } # uses REML COV formula

    # fast calculation of sigma covariance
    COVSTUFF <- COV.covm(CTMM$sigma,n=n,k=k.mean,REML=REML)
    CTMM$COV <- COVSTUFF$COV

    MLE <- NULL
  }
  else ### all further cases require optimization ###
  {
    if(trace) { message("Maximizing likelihood.") }
    # fix for bad starting conditions
    if(class(data)[1]=="phylometry" && length(CTMM$tau)>1)
    {
      TEST <- fn(pars)
      while(TEST>=Inf)
      {
        if(CTMM$tau[1]==CTMM$tau[2])
        { pars['tau'] <- pars['tau']/2 }
        else
        { pars['tau velocity'] <- pars['tau velocity']/2 }
        TEST <- fn(pars)
      }
    }
    # control$covariance <- covariance() - parameter storage and optimization representations can differ
    control$parscale <- parscale
    control$zero <- TRUE
    RESULT <- optimizer(par=pars,fn=fn,method=op.method,lower=lower,upper=upper,period=period,reset=reset,control=control)
    pars <- clean.parameters(RESULT$par,profile=profile,timelink=CTMM$timelink)
    # copy over hessian from fit to COV.init ?

    # write best estimates over initial guess
    store.pars <- function(pars,profile=profile,finish=TRUE)
    {
      names(pars) <- NAMES
      pars <- clean.parameters(pars,profile=profile,linear.cov=linear.cov,timelink=CTMM$timelink)

      CTMM <- set.parameters(CTMM,pars,profile=profile,linear.cov=linear.cov)

      # this is a wasted evaluation !!! store verbose glob in environment?
      if(finish) { CTMM <- loglike.fn(data,CTMM,REML=REML,verbose=TRUE,profile=TRUE) }

      return(CTMM)
    }
    CTMM <- store.pars(pars,profile=profile,finish=TRUE)

    profile <- FALSE # no longer solving covariance analytically, no matter what
    setup.parameters(CTMM,profile=FALSE)
    ### COV CALCULATION #############
    if(COV || method %in% c("pREML","pHREML"))
    {
      # if pREML, calculate Hessian in safe parameterization and then transform afterwards
      if(trace) { message("Calculating Hessian.") }
      DIFF <- genD(par=pars,fn=fn,zero=-CTMM$loglike,lower=lower,upper=upper,parscale=parscale,Richardson=2,mc.cores=1,control=control)
      hess <- DIFF$hessian
      grad <- DIFF$gradient

      # more robust covariance calculation than straight inverse
      CTMM$COV <- cov.loglike(hess,grad)
      dimnames(CTMM$COV) <- list(NAMES,NAMES)
    }

    # store MLE for faster re-optimization #
    MLE <- unscale.ctmm(CTMM)

    ### pREML correction ########## only do pREML if sufficiently away from boundaries
    PREML.FAIL <- method %in% c("pREML","pHREML") && mat.min(hess) <= .Machine$double.eps*length(NAMES)
    if(method %in% c("pREML","pHREML") && !PREML.FAIL)
    {
      # parameter correction
      REML <- TRUE
      #ML.grad <- grad # save old ML gradient
      if(trace) { message("Calculating REML gradient.") }
      DIFF <- genD(par=pars,fn=fn,lower=lower,upper=upper,parscale=parscale,Richardson=2,order=1,mc.cores=1,control=control)

      J <- diag(length(pars))
      dimnames(J) <- list(NAMES,NAMES)
      # invoke Jacobian (major,minor,angle) -> linear parameterization
      if(!CTMM$isotropic)
      {
        # Jacobian matrix d sigma / d par
        SUB <- names(CTMM$sigma@par)
        J[SUB,SUB] <- J.sigma.par(CTMM$sigma@par)
      }
      if(any(UERE.FIT)) # RMS UERE -> MS UERE (linear parameterization)
      {
        UERE.EPAR <- paste("error",UERE.PAR)
        if(length(UERE.PAR)==1)
        { J[UERE.EPAR,UERE.EPAR] <- 2*CTMM$error[UERE.FIT] }
        else # diag() is annoying
        { diag(J[UERE.EPAR,UERE.EPAR]) <- 2*CTMM$error[UERE.FIT] }
      }

      # apply linear parameter correction
      linear.cov <- TRUE
      setup.parameters(CTMM,profile=FALSE)

      # calculate linear parameter correction
      d.pars <- -c(J %*% CTMM$COV %*% DIFF$gradient)
      # Jacobians cancel out between inverse Hessian and gradient

      # gradient or something else was bad
      PREML.FAIL <- any(is.na(d.pars)) || any(abs(d.pars)==Inf)

      # increment transformed parameters
      # pars <- pars + d.pars
      if(!PREML.FAIL)
      {
        # safety catch for bad models near boundaries
        pars <- line.boxer(d.pars,pars,lower=lower,upper=upper,period=period)
        names(pars) <- NAMES

        # store parameter correction
        profile <- FALSE
        # this can still fail if parameters are crazy
        TEST <- store.pars(pars,profile=FALSE,finish=TRUE)
        if(class(TEST)[1]=="ctmm" && !is.na(TEST$loglike) && TEST$loglike>-Inf)
        {
          CTMM <- TEST
          linear.cov <- FALSE
          setup.parameters(CTMM,profile=FALSE)
        }
        else # parameters crashed likelihood function---this is rare
        { PREML.FAIL <- TRUE }
      }
      linear.cov <- FALSE # 2 fail cases need this
    }
    # revert to ML/HREML if pREML step failed
    if(method %in% c("pREML","pHREML") && PREML.FAIL)
    {
      # shouldn't need to warn if using ctmm.select
      WARN <- TRUE
      N <- sys.nframe()
      if(N>=2)
      {
        for(i in 2:N)
        {
          CALL <- deparse(sys.call(-i))[1]
          CALL <- grepl("ctmm.select",CALL) || grepl("cv.like",CALL) || grepl("ctmm.boot",CALL)
          if(CALL)
          {
            WARN <- FALSE
            break
          }
        }
      }
      # warn if weren't using ctmm.select
      if(WARN) { warning("pREML failure: indefinite ML Hessian or divergent REML gradient.") }

      if(method=='pREML') { method <- 'ML' }
      else if(method=='pHREML') { method <- 'HREML' }
    }
    ### end pREML correction ###

    ### HREML STEP: profile linear REML parameters ###
    if(method %in% c('pHREML','HREML'))
    {
      REML <- TRUE
      profile <- !any(UERE.FIX)

      # profile REML linear parameters numerically if necessary (error || circle)
      setup.parameters(CTMM,profile=profile,linear=TRUE)
      if(length(NAMES))
      {
        REML <- TRUE
        if(trace) { message("Profiling REML likelihood.") }
        # control$covariance <- covariance() parameter storage and optimization representations can differ
        control$parscale <- parscale
        control$zero <- TRUE
        RESULT <- optimizer(par=pars,fn=fn,method=op.method,lower=lower,upper=upper,period=period,reset=reset,control=control)
        pars <- clean.parameters(RESULT$par,timelink=CTMM$timelink)
      }

      # includes free profile
      TEST <- store.pars(pars,profile=profile,finish=TRUE)
      if(class(TEST)[1]=="ctmm")
      { CTMM <- TEST }
      else # parameters crashed likelihood function---this is extremely rare
      {
        if(method=="pHREML") { method <- "pREML" }
        else if(method=="HREML") { method <- "ML" }
      }
    }

    ### FINAL COVARIANCE ESTIMATE ###
    if(COV && method %in% c('pREML','pHREML','HREML')) ### CALCULATE COVARIANCE MATRIX ###
    {
      profile <- FALSE
      setup.parameters(CTMM,profile=FALSE)

      if(trace) { message("Calculating REML Hessian.") }
      # calcualte REML Hessian at pREML parameters
      DIFF <- genD(par=pars,fn=fn,lower=lower,upper=upper,parscale=parscale,Richardson=2,mc.cores=1,control=control)
      # Using MLE gradient, which should be zero off boundary
      CTMM$COV <- cov.loglike(DIFF$hessian,grad)
    }

    if(COV) { dimnames(CTMM$COV) <- list(NAMES,NAMES) }
  } # end optimized estimates

  # model likelihood (not REML for AIC)
  if(method!='ML') { CTMM$loglike <- loglike.fn(data,CTMM=CTMM,REML=FALSE,profile=FALSE) }
  CTMM$method <- method

  # covariance parameters only
  setup.parameters(CTMM,profile=FALSE,linear=FALSE)

  # unstandardize (includes likelihood adjustment)
  CTMM <- unscale.ctmm(CTMM)
  CTMM$features <- NAMES # store all auto-covariance features

  # calculate AIC,AICc,BIC,MSPE,...
  CTMM$link <- link$name
  CTMM$loglike <- CTMM$loglike + log.like.link
  CTMM <- ic.ctmm(CTMM,n)

  # would be temporary ML COV for pREML/pHREML
  if(!COV) { CTMM$COV <- NULL }

  if(method %in% c('pREML','pHREML','HREML') && !is.null(MLE))
  {
    MLE$link <- link$name
    MLE$loglike <- MLE$loglike + log.like.link
    # calculate checksum
    MLE$checksum <- digest::digest(CTMM,algo="md5")
    CTMM$MLE <- MLE
    # now if anyone modifies CTMM, then MLE will not be used
  }

  return(CTMM)
}


#################
# calculate AIC/BIC/AICc/MSPE/...
#################
ic.ctmm <- function(CTMM,n)
{
  axes <- CTMM$axes
  range <- CTMM$range
  k.mean <- length(CTMM$mu)
  method <- CTMM$method

  # all parameters
  q <- length(axes)
  if(is.null(CTMM$formula)) # autocorrelation model
  {
    NAMES <- CTMM$features
    nu <- length(NAMES)

    if(!range)
    {
      k.mean <- k.mean - q
      n <- n - 1
    }
    if(!length(k.mean)) { k.mean <- 0 } # failed fit (bad data or bad parameters)
  }
  else # RSF model
  {
    k.beta <- length(all.vars(CTMM$formula)) # linear terms
    nu.beta <- length(CTMM$beta) - k.beta # quadratic terms (and more)
    k.mean <- 2 + k.beta
    nu <- 1 + nu.beta
  }
  k <- nu + k.mean

  CTMM$AIC <- 2*k - 2*CTMM$loglike
  CTMM$BIC <- log(n)*k - 2*CTMM$loglike

  # IID AICc values
  if(method=='ML')
  { CTMM$AICc <- -2*CTMM$loglike + q*n * 2*k/(q*n-k-nu) }
  else if(method=='pREML')
  { CTMM$AICc <- -2*CTMM$loglike + (q*n)^2/(q*n+k.mean) * 2*k/(q*n-k-nu) }
  else if(method %in% c('pHREML','HREML','REML'))
  { CTMM$AICc <- -2*CTMM$loglike + (q*n-k.mean) * 2*k/(q*n-k-nu) }

  # fix divergence
  if(q*n<=k+nu) { CTMM$AICc <- Inf }

  # Mean square prediction error
  mspe <- function(K=1)
  {
    if(!CTMM$range && K==1) { return(Inf) }

    # velocity MSPE
    if(K>=2)
    {
      if(length(CTMM$tau)<2 || any(CTMM$tau<=0)) { return(Inf) }
      UU <- VV
    }

    MSPE <- CTMM$COV.mu
    if(length(axes)==1)
    {
      if(length(MSPE)==1) # one spatial dimension and one trend component
      { MSPE <- MSPE * UU }
      else # one spatial dimension and many trend components
      { MSPE <- diag(MSPE %*% UU) }
    }
    else
    {
      if(length(UU)==1) # multiple spatial dimensions and one trend component
      { MSPE <- sum(diag(MSPE)) * c(UU) }
      # else if(length(dim(MSPE))==2 && length(UU)>1)
      else if(length(dim(MSPE))==4) # k trend components in multiple dimensions
      {
        MSPE <- lapply(1:length(axes),function(i) MSPE[i,,,i])
        MSPE <- Reduce("+",MSPE)
        MSPE <- diag(MSPE %*% UU)
      }
      else
      { stop("Inconsistent dimensions around COV[mu].") }
    }
    # 0/0 -> 0 (IOU)
    MSPE <- sum(nant(MSPE,0))

    VAR <- sum(diag(CTMM$sigma))
    if(K==2)
    {
      STUFF <- get.taus(CTMM)
      VAR <- VAR * (STUFF$Omega2 + CTMM$circle^2)
    }

    MSPE <- MSPE + VAR

    return(MSPE)
  }
  STUFF <- drift.energy(CTMM)
  UU <- STUFF$UU
  VV <- STUFF$VV
  # Mean square prediction error in locations & velocities
  CTMM$MSPE <- c( mspe(K=1) , mspe(K=2) )
  names(CTMM$MSPE) <- c("position","velocity")

  return(CTMM)
}


###################
# general parameter guessing function
###################
ctmm.guess <- function(data,CTMM=ctmm(),variogram=NULL,name="GUESS",interactive=TRUE)
{
  check.class(data)

  # use intended axes
  if(is.null(variogram)) { variogram = variogram(data,axes=CTMM$axes) }
  else { CTMM$axes <- attr(variogram,"info")$axes }

  n <- length(data$t)
  if(n<=2) { CTMM$isotropic = TRUE }

  CTMM <- ctmm.prepare(data,CTMM,precompute=FALSE,tau=FALSE)

  # mean specific guesswork/preparation
  CTMM <- drift.init(CTMM,data)
  mu <- CTMM$mu
  u <- drift.mean(CTMM,data$t)

  # estimate circulation period if circle=TRUE
  if(CTMM$circle && class(CTMM$circle)[1]=="logical")
  {
    # residuals
    z <- get.telemetry(data,CTMM$axes)
    z <- z - (u %*% CTMM$mu)

    dt <- diff(data$t)
    SUB <- dt>0

    # velocities !!! update this to minimally filtered estimate
    v <- cbind(diff(z[,1]),diff(z[,2])) / dt
    # midpoint locations during velocity v
    z <- cbind(z[-1,1]+z[-n,1],z[-1,2]+z[-n,2])/2

    # average angular momentum
    L <- c(z[SUB,1]%*%v[SUB,2] - z[SUB,2]%*%v[SUB,1]) / (n-1)

    circle <- L / mean(diag(CTMM$sigma))
    # circle <- 2*pi/circle

    CTMM$circle <- circle
  }

  # other stuff...

  variogram.fit(variogram,CTMM=CTMM,name=name,interactive=interactive)
}
