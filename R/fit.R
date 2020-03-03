# smallest resolutions in data (for soft bounding parameter guesstimates)
telemetry.mins <- function(data,axes=c('x','y'))
{
  dt <- diff(data$t)
  dt <- stats::median(dt[dt>0]) # median time difference

  df <- 2*pi/(last(data$t)-data$t[1]) # smallest frequency

  dz <- get.telemetry(data,axes)
  dz <- apply(dz,2,diff)
  dim(dz) <- c(nrow(data)-1,length(axes)) # R drops length-1 dimensions
  dz <- rowSums( dz^2 )
  dz <- sqrt(min(dz[dz>0])) # smallest nonzero distance

  return(list(dt=dt,df=df,dz=dz))
}


###########################################################
# FIT MODEL WITH LIKELIHOOD FUNCTION (convenience wrapper to optim)
ctmm.fit <- function(data,CTMM=ctmm(),method="pHREML",COV=TRUE,control=list(),trace=FALSE)
{
  # DEBUG.FIT <<- list(data=data,CTMM=CTMM,method=method,COV=COV,control=control,trace=trace,SEED=.Random.seed)
  if(!is.null(control$message)) { message <- control$message }

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

  # standardize data for numerical stability
  # pre-centering the data reduces later numerical error across models (tested)
  SHIFT <- colMeans(get.telemetry(data,axes))
  data[,axes] <- t( t(data[,axes,drop=FALSE]) - SHIFT )
  # add mu.center back to the mean value after kalman filter / mean profiling
  # pre-standardizing the data should also help
  SCALE <- sqrt(mean(get.telemetry(data,axes)^2))
  # standardize time by median diff time
  DT <- diff(data$t)
  TSCALE <- stats::median(DT[DT>0])
  data <- unit.telemetry(data,length=SCALE,time=TSCALE)
  CTMM <- unit.ctmm(CTMM,length=SCALE,time=TSCALE)

  # used for minimum scale of parameter inspection
  n <- length(data$t)
  MINS <- telemetry.mins(data,axes)
  dt <- MINS$dt
  df <- MINS$df
  dz <- MINS$dz

  # unstandardize (includes likelihood adjustment)
  unscale.ctmm <- function(CTMM)
  {
    CTMM <- unit.ctmm(CTMM,length=1/SCALE,time=1/TSCALE)
    # log-likelihood adjustment
    CTMM$loglike <- CTMM$loglike - length(axes)*n*log(SCALE)
    # translate back to origin from center
    CTMM$mu <- drift@shift(CTMM$mu,SHIFT)

    return(CTMM)
  }

  default <- list(method="pNewton",precision=1/2,maxit=.Machine$integer.max)
  control <- replace(default,names(control),control)
  precision <- control$precision
  op.method <- control$method
  control$method <- NULL

  # clean/validate
  CTMM <- ctmm.ctmm(CTMM)
  drift <- get(CTMM$mean)
  CTMM$mu <- NULL # can always profile mu analytically
  range <- CTMM$range

  # save for fitting
  COV.init <- CTMM$COV
  # make sure we can start from previous failed fit
  if(any(is.nan(COV.init) | COV.init==Inf)) { COV.init <- NULL }
  if(!is.null(COV.init)) { TEST <- eigen(COV.init,only.values=TRUE)$values } else { TEST <- FALSE }
  if(any(TEST<=.Machine$double.eps | TEST==Inf)) { COV.init <- NULL }
  # erase previous fitting info if present
  CTMM$COV <- NULL
  CTMM$COV.mu <- NULL
  CTMM$DOF.mu <- NULL

  # evaluate mean function and error matrices for this data once upfront
  CTMM <- ctmm.prepare(data,CTMM,tau=FALSE) # don't muck with taus
  UERE <- attr(CTMM$error.mat,"flag") # do we fit the error? Need to know for optimization

  ### id and characterize parameters for profiling ###
  pars <- NAMES <- parscale <- lower <- upper <- period <- NULL
  ORIGINAL <- CTMM # original structure of model before fitting
  linear.cov <- FALSE # represent sigma linearly (for perturbation) versus * (for optimization)
  setup.parameters <- function(CTMM,profile=TRUE,linear=FALSE)
  {
    STUFF <- id.parameters(CTMM,profile=profile,linear=linear,linear.cov=linear.cov,UERE=UERE,dt=dt,df=df,dz=dz,STRUCT=ORIGINAL)
    NAMES <<- STUFF$NAMES
    parscale <<- STUFF$parscale
    lower <<- STUFF$lower
    upper <<- STUFF$upper
    period <<- STUFF$period
    # reflect <<- STUFF$reflect
    # initial guess for optimization
    pars <<- get.parameters(CTMM,NAMES,linear.cov=linear.cov)
  }
  setup.parameters(CTMM,profile=TRUE)
  if("error" %nin% NAMES) { CTMM$error <- as.logical(CTMM$error) } # fix numeric error when it should be logical

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
    p <- clean.parameters(p,linear.cov=linear.cov)
    CTMM <- set.parameters(CTMM,p,linear.cov=linear.cov)

    # negative log likelihood
    return(-ctmm.loglike(data,CTMM,REML=REML,zero=-zero,profile=profile))
  }

  # shift parameters back near origin to prevent overflow in optimizer
  reset <- function(p)
  {
    names(p) <- NAMES
    if(!linear.cov)
    {
      # swap minor and major axes if minor is getting too big
      if("minor" %in% NAMES && "major" %nin% NAMES) # major axis is being profiled
      {
        major <- CTMM$sigma@par['major']
        ERROR <- "error" %in% NAMES
        TEST <- (p['minor']>major)
        if(ERROR) { TEST <- TEST && p['minor']>p['error']^2 }
        if(TEST) # swap major-minor
        {
          # rotate to true minor axis
          p["angle"] <- p["angle"] + pi/2

          # existing variance ratios to maintain
          RATIO <- major/p['minor'] # <1

          p['minor'] <- RATIO * major # old major is the new minor (proportionally)
          if(ERROR) { p['error'] <- sqrt(RATIO) * p['error'] }
        }
      } # end major-minor swap

      # shift back to (-pi/2,pi/2)
      # if("angle" %in% NAMES) { p["angle"] <- (((p["angle"]/pi+1/2) %% 1) - 1/2)*pi }
    } # end !linear.cov

    return(p)
  } # end reset

  # construct covoariance matrix guess
  covariance <- function()
  {
    COV <- diag(parscale^2,nrow=length(parscale))
    dimnames(COV) <- list(NAMES,NAMES)
    COPY <- rownames(COV.init) %in% NAMES
    if(any(COPY))
    {
      COPY <- rownames(COV.init)[COPY]
      COV[COPY,COPY] <- COV.init[COPY,COPY]
    }
    return(COV)
  }

  ### NOW OPTIMIZE ###
  profile <- TRUE
  if(length(NAMES)==0) # EXACT
  {
    if(method %in% c("pHREML","HREML")) { REML <- TRUE } # IID limit pHREML/HREML -> REML

    # Bi-variate Gaussian with zero error
    CTMM <- ctmm.loglike(data,CTMM=CTMM,REML=REML,verbose=TRUE)

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
    control$covariance <- covariance()
    control$parscale <- parscale
    control$zero <- TRUE
    RESULT <- optimizer(par=pars,fn=fn,method=op.method,lower=lower,upper=upper,period=period,reset=reset,control=control)
    pars <- clean.parameters(RESULT$par)
    # copy over hessian from fit to COV.init ?

    # write best estimates over initial guess
    store.pars <- function(pars,profile=TRUE,finish=TRUE)
    {
      names(pars) <- NAMES
      pars <- clean.parameters(pars,linear.cov=linear.cov)

      CTMM <- set.parameters(CTMM,pars,linear.cov=linear.cov)

      # this is a wasted evaluation !!! store verbose glob in environment?
      if(finish) { CTMM <- ctmm.loglike(data,CTMM,REML=REML,verbose=TRUE,profile=profile) }

      return(CTMM)
    }
    CTMM <- store.pars(pars,finish=TRUE)

    profile <- FALSE # no longer solving covariance analytically
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
      # invoke Jacobian (major,minor,angle) -> linear parameterization !!!
      if(!CTMM$isotropic)
      {
        # Jacobian matrix d sigma / d par
        SUB <- names(CTMM$sigma@par)
        J[SUB,SUB] <- J.sigma.par(CTMM$sigma@par)
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
        if(class(TEST)[1]=="ctmm")
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
      warning("pREML failure: indefinite ML Hessian or divergent REML gradient.")
      if(method=='pREML') { method <- 'ML' }
      else if(method=='pHREML') { method <- 'HREML' }
    }
    ### end pREML correction ###

    ### HREML STEP: profile linear REML parameters ###
    if(method %in% c('pHREML','HREML'))
    {
      REML <- TRUE
      profile <- TRUE

      # profile REML linear parameters numerically if necessary (error || circle)
      setup.parameters(CTMM,profile=TRUE,linear=TRUE)
      if(length(NAMES))
      {
        REML <- TRUE
        if(trace) { message("Profiling REML likelihood.") }
        control$covariance <- covariance()
        control$parscale <- parscale
        control$zero <- TRUE
        RESULT <- optimizer(par=pars,fn=fn,method=op.method,lower=lower,upper=upper,period=period,reset=reset,control=control)
        pars <- clean.parameters(RESULT$par)
      }

      # includes free profile
      TEST <- store.pars(pars,profile=TRUE,finish=TRUE)
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
  if(method!='ML') { CTMM$loglike <- ctmm.loglike(data,CTMM=CTMM,REML=FALSE,profile=FALSE) }
  CTMM$method <- method

  # covariance parameters only
  setup.parameters(CTMM,profile=FALSE,linear=FALSE)

  # unstandardize (includes likelihood adjustment)
  CTMM <- unscale.ctmm(CTMM)
  CTMM$features <- NAMES # store all auto-covariance features

  # calculate AIC,AICc,BIC,MSPE,...
  CTMM <- ic.ctmm(CTMM,n)

  # would be temporary ML COV for pREML/pHREML
  if(!COV) { CTMM$COV <- NULL }

  if(method %in% c('pREML','pHREML','HREML') && !is.null(MLE))
  {
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
  NAMES <- CTMM$features
  axes <- CTMM$axes
  range <- CTMM$range
  k.mean <- nrow(CTMM$mu)
  method <- CTMM$method

  nu <- length(NAMES)
  # all parameters
  q <- length(axes)
  if(!range)
  {
    k.mean <- k.mean - 1
    n <- n - 1
  }
  k <- nu + q*k.mean

  CTMM$AIC <- 2*k - 2*CTMM$loglike
  CTMM$BIC <- log(n)*k - 2*CTMM$loglike

  # IID AICc values
  if(method=='ML')
  { CTMM$AICc <- -2*CTMM$loglike + q*n * 2*k/(q*n-k-nu) }
  else if(method=='pREML')
  { CTMM$AICc <- -2*CTMM$loglike + (q*n)^2/(q*n+q*k.mean) * 2*k/(q*n-k-nu) }
  else if(method %in% c('pHREML','HREML','REML'))
  { CTMM$AICc <- -2*CTMM$loglike + (q*n-q*k.mean) * 2*k/(q*n-k-nu) }

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
  drift <- get(CTMM$mean)
  STUFF <- drift@energy(CTMM)
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
  #

  # use intended axes
  if(is.null(variogram)) { variogram = variogram(data,axes=CTMM$axes) }
  else { CTMM$axes <- attr(variogram,"info")$axes }

  n <- length(data$t)
  if(n==2) { CTMM$isotropic = TRUE }

  # mean specific guesswork/preparation
  drift <- get(CTMM$mean)
  CTMM <- drift@init(data,CTMM)
  mu <- CTMM$mu
  u <- drift(data$t,CTMM)

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

  variogram.fit(variogram,CTMM=CTMM,name=name,interactive=interactive)
}
