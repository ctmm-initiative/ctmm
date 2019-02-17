####################################
# log likelihood function
####################################
ctmm.loglike <- function(data,CTMM=ctmm(),REML=FALSE,profile=TRUE,zero=0,verbose=FALSE)
{
  n <- length(data$t)
  AXES <- length(CTMM$axes)

  # save original tau length
  K <- length(CTMM$tau)
  # prepare model for numerics
  CTMM <- ctmm.prepare(data,CTMM)

  drift <- get(CTMM$mean)

  range <- CTMM$range
  isotropic <- CTMM$isotropic

  COVM <- function(...) { covm(...,isotropic=CTMM$isotropic,axes=CTMM$axes) } # enforce model structure
  if(!is.null(CTMM$sigma))
  {
    CTMM$sigma <- M.sigma <- COVM(CTMM$sigma) # model-only covariance
    sigma <- CTMM$sigma # model->profile covariance

    area <- sigma@par[1]
    ecc <- sigma@par[2]
    theta <- sigma@par[3]
  }
  else
  {
    M.sigma <- COVM(diag(AXES))

    area <- 1 # only makes sense with profile
    ecc <- 0
    theta <- 0
  }

  circle <- CTMM$circle

  n <- length(data$t)

  t <- data$t
  # time lags
  dt <- c(Inf,diff(t))

  # data z and mean vector u
  z <- get.telemetry(data,CTMM$axes)
  u <- CTMM$mean.vec
  M <- ncol(u) # number of linear parameters per spatial dimension

  # pre-centering the data reduces later numerical error across models (tested)
  # mu.center <- colMeans(z)
  # z <- t( t(z) - mu.center )
  # add mu.center back to the mean value after kalman filter / mean profiling
  # pre-standardizing the data would also help

  # REML variance debias factor # ML constant
  if(REML) { VAR.MULT <- n/(n-M) } else  { VAR.MULT <- 1 }

  # get the error information
  error <- CTMM$error.mat # note for fitted errors, this is error matrix @ UERE=1 (CTMM$error)
  # are we fitting the error, then the above is not yet normalized
  UERE <- attr(error,"flag")

  # check for bad time intervals
  ZERO <- which(dt==0)
  if(length(ZERO) && K && CTMM$tau[1])
  {
    if(CTMM$error==FALSE) { warning("Duplicate timestamps require an error model.") ; return(-Inf) }
    # check for HDOP==0 just in case
    ZERO <- error[ZERO,,,drop=FALSE]
    ZERO <- apply(ZERO,1,det)
    if(any(ZERO<=0)) { warning("Duplicate timestamps require an error model.") ; return(-Inf) }
  }

  # check for bad variances
  if(area==0)
  {
    ZERO <- rep(FALSE,n)
    for(i in 1:dim(error)[2]) { ZERO <- ZERO | (error[,i,i]==0) }
    if(any(ZERO)) { return(-Inf) }
  }

  ### what kind of profiling is possible
  if((!UERE && !(circle && !isotropic)) || (UERE<3 && isotropic)) # can profile full covariance matrix with unit-variance filter
  {
    PROFILE <- 2
    est.cov <- function(COV) { COVM(VAR.MULT*COV) }
  }
  else if(UERE<3) # can profile variance with fixed-eccentricity 2D or 2x1D filters
  {
    PROFILE <- TRUE
    est.cov <- function(COV)
    {
      if(length(COV)>1) { COV <- mean(diag(COV)) }
      COV <- VAR.MULT * COV
      if(ROTATE) { theta <- 0 } # already along axis
      if(SQUEEZE) { ecc <- 0 } # in circular coordinate system
      COV <- COVM(c(COV,ecc,theta))
      return(COV)
    }
  }
  else # no profiling possible - need exact-covariance filter
  {
    PROFILE <- FALSE
    est.cov <- function(COV) { sigma }
  }

  ### 2D or 1D Kalman filters necessary?
  if(!circle && !isotropic && UERE && UERE<4) # can run 2x1D Kalman filters
  { DIM <- 1/2 }
  else if(UERE==4 || (circle && !isotropic && UERE)) # full 2D filter necessary
  { DIM <- 2 }
  else # can run 1x1D Kalman filter
  { DIM <- 1 }

  ### calibrate unknown errors given PROFILE state
  if(UERE && UERE<3) # calibrate errors
  {
    if(PROFILE) { CTMM$error <- CTMM$error / sqrt(area) } # fix error/variance ratio
    error <- CTMM$error^2 * error
  }

  # orient the data along the major & minor axes of sigma to run 2x1D filters - or to do basic circulation transformation after squeezing
  ROTATE <- !isotropic && (circle || (UERE && UERE<4))
  if(ROTATE)
  {
    R <- rotate(-theta)
    z <- z %*% t(R)

    # minimize numerical errors from transforming back and forth
    M.sigma <- attr(M.sigma,"par")
    M.sigma["angle"] <- 0
    M.sigma <- COVM(M.sigma)

    if(UERE>=4) { error <- rotate.mat(error,-theta) } # rotate error ellipses
  }

  # squeeze from ellipse to circle for circulation model transformation
  SQUEEZE <- (circle && !isotropic)
  if(SQUEEZE)
  {
    z <- squeeze(z,ecc)

    # minimize numerical errors from transforming back and forth
    M.sigma <- attr(M.sigma,"par")
    M.sigma["eccentricity"] <- 0
    M.sigma <- COVM(M.sigma)

    # squeeze error circles into ellipses
    if(UERE) { error <- squeeze.mat(error,ecc) }

    # mean functions can requires squeezing (below)
    # how does !error && circle && !isotropic avoid this?
  }

  ## some cases need separate x & y mean functions
  if(DIM==2) # u -> cbind( (u,0) , (0,u) )
  {
    if(!SQUEEZE)
    { u <- vapply(1:ncol(u),function(i){cbind(u[,i],rep(0,n),rep(0,n),u[,i])},array(0,c(n,4))) }
    else
    { u <- vapply(1:ncol(u),function(i){cbind(u[,i]*exp(-ecc/4),rep(0,n),rep(0,n),u[,i]*exp(+ecc/4))},array(0,c(n,4))) }
    dim(u) <- c(n,4*M) # <- (n,2,2*M)
  }
  else if(circle) # 1D circle can use symmetry with just u -> (u,0) mean function for x
  {
    # this relies on trickery later on, because of no SQUEEZE
    u <- vapply(1:ncol(u),function(i){cbind(u[,i],rep(0,n))},array(0,c(n,2))) # (n,2,M)
    dim(u) <- c(n,2*M)
  }

  if(circle) ## COROTATING FRAME FOR circle=TRUE ##
  {
    R <- rotates(-circle*(t-t[1])) # rotation matrices
    u <- rotates.vec(cbind(z,u),R)
    z <- u[,1:2]
    u <- u[,-(1:2)]
    if(UERE>=4 || (UERE && SQUEEZE)) { error <- rotates.mat(error,R) } # rotate error ellipses
    rm(R)
  }

  #### RUN KALMAN FILTERS ###
  if(DIM==1/2) ### 2 separate 1D Kalman filters instead of 1 2D Kalman filter ###
  {
    # scale eigen variances proportionally
    SIGMA <- exp(c(+1,-1)*ecc/2)
    if(!PROFILE) { SIGMA <- SIGMA * area }

    # major axis likelihood
    CTMM$sigma <- SIGMA[1]
    KALMAN1 <- kalman(cbind(z[,1]),u,dt=dt,CTMM=CTMM,error=error)

    # minor axis likelihood
    CTMM$sigma <- SIGMA[2]
    KALMAN2 <- kalman(cbind(z[,2]),u,dt=dt,CTMM=CTMM,error=error)

    mu <- cbind(KALMAN1$mu,KALMAN2$mu)

    R.sigma <- c(KALMAN1$sigma + KALMAN2$sigma)/2 # residual variance
    if(PROFILE && !profile) { RATIO <- R.sigma/area } # ratio of residual variance to model variance
    if(profile)
    {
      sigma <- est.cov(R.sigma)
      area <- attr(sigma,'par')['area']
    }

    # combine uncorrelated estimates
    COV.mu <- array(c(KALMAN1$iW,diag(0,M),diag(0,M),KALMAN2$iW),c(M,M,2,2)) # -1/Hessian
    if(PROFILE) { COV.mu <- COV.mu * area }
    DOF.mu <- array(c(SIGMA[1]*KALMAN1$W,diag(0,M),diag(0,M),SIGMA[2]*KALMAN2$W),c(M,M,2,2))

    # put in first canonical form (2*M,2*M)
    COV.mu <- aperm(COV.mu,c(3,1,4,2)) ; dim(COV.mu) <- c(2*M,2*M)
    DOF.mu <- aperm(DOF.mu,c(3,1,4,2)) ; dim(DOF.mu) <- c(2*M,2*M)

    logdetCOV <- KALMAN1$logdet + KALMAN2$logdet
    logdetcov <- -log(det(KALMAN1$W)) - log(det(KALMAN2$W))
  }
  else ### 1x 1D or 2D Kalman filter ###
  {
    # prepare variance/covariance of 1D/2D Kalman filters
    if(PROFILE==2 || (PROFILE && DIM==1)) { CTMM$sigma <- 1 } # could be rotated & squeezed
    else if(DIM==1) { CTMM$sigma <- area }
    else if(DIM==2) { CTMM$sigma <- COVM(c(if(PROFILE){1}else{area},if(SQUEEZE){0}else{ecc},if(ROTATE){0}else{theta})) } # circle, !isotropic, UERE=1,2
    # else sigma is full covariance matrix

    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)

    mu <- KALMAN$mu

    R.sigma <- KALMAN$sigma # residual covariance
    if(PROFILE==2 && !profile) { RATIO <- mean(diag( COVM(R.sigma) %*% PDsolve(M.sigma) )) }
    else if(PROFILE && !profile) { RATIO <- mean(diag(R.sigma))/area } # filter already partially standardized via eccentricity
    if(profile) # variance/covariance update based on residual covariance
    {
      sigma <- est.cov(R.sigma)
      area <- attr(sigma,'par')['area']
    }

    logdetCOV <- (AXES/DIM)*KALMAN$logdet # log autocovariance / n

    ### mu covariance terms
    if(DIM==2 || circle) # cases where we have a u(t) for each dimension
    {
      logdetcov <- -log(det(KALMAN$W)) # log cov[beta] absolute

      COV.mu <- KALMAN$iW # (2*M,2*M) from riffle
      # sqrt(sigma)
      SQRTM <- sqrtm(sigma) %o% diag(1,M) # (2,2,M,M)
      SQRTM <- aperm(SQRTM,c(1,3,2,4)) # (2,M,2,M)
      dim(SQRTM) <- c(2*M,2*M)
      DOF.mu <- SQRTM %*% KALMAN$W %*% SQRTM # (2*M,2*M)

      if(PROFILE) # variance was not included
      {
        COV.mu <- COV.mu * area # iW=COV/area
        DOF.mu <- DOF.mu / area # W=area/COV
      }

      logdetcov <- log(det(COV.mu)) # (2*M,2*M)
    }
    else # lower or 1D filter return - cases where we have one u(t) for all dimensions
    {
      logdetcov <- -AXES*log(det(KALMAN$W)) # log cov[beta] absolute

      if(PROFILE)
      {
        COV.mu <- KALMAN$iW %o% sigma # (M,M,2,2)
        DOF.mu <- KALMAN$W %o% diag(AXES)
      }
      else if(!PROFILE)
      {
        COV.mu <- KALMAN$iW %o% diag(AXES)
        DOF.mu <- (area * KALMAN$W) %o% diag(AXES)
      }
      # put either in first canonical form (AXES*M,AXES*M)
      COV.mu <- aperm(COV.mu,c(3,1,4,2)) ; dim(COV.mu) <- c(AXES*M,AXES*M)
      DOF.mu <- aperm(DOF.mu,c(3,1,4,2)) ; dim(DOF.mu) <- c(AXES*M,AXES*M)
    }
  }   ### END KALMAN FILTER RUNS ###

  if(PROFILE) # missing variances from profiling
  {
    logdetCOV <- logdetCOV + AXES*log(area) # per n
    logdetcov <- logdetcov + (AXES*M)*log(area) # absolute
  }

  if(SQUEEZE) # de-squeeze
  {
    sigma <- attr(sigma,"par")
    sigma["eccentricity"] <- ecc
    sigma <- COVM(sigma)

    if(DIM==1) # DIM==2 squeezed u(t), so beta and COV[beta] have normal scale already
    {
      R <- exp(c(+1,-1)*ecc/4)
      mu <- t(R * t(mu)) # (M,2)

      dim(COV.mu) <- c(2,M*2*M)
      COV.mu <- R * COV.mu
      dim(COV.mu) <- c(2,M,2,M)
      COV.mu <- aperm(COV.mu,c(3,4,1,2))
      dim(COV.mu) <- c(2,M*2*M)
      COV.mu <- R * COV.mu
      dim(COV.mu) <- c(2*M,2*M)
    }
  }

  if(ROTATE) # transform results back
  {
    R <- rotate(+theta)

    mu <- mu %*% t(R)

    sigma <- attr(sigma,"par")
    sigma["angle"] <- theta
    sigma <- COVM(sigma)

    dim(COV.mu) <- c(2,M*2*M)
    dim(DOF.mu) <- c(2,M*2*M)

    COV.mu <- R %*% COV.mu
    DOF.mu <- R %*% DOF.mu

    dim(COV.mu) <- c(2,M,2,M)
    dim(DOF.mu) <- c(2,M,2,M)

    COV.mu <- aperm(COV.mu,c(3,4,1,2))
    DOF.mu <- aperm(DOF.mu,c(3,4,1,2))

    dim(COV.mu) <- c(2,M*2*M)
    dim(DOF.mu) <- c(2,M*2*M)

    COV.mu <- R %*% COV.mu
    DOF.mu <- R %*% DOF.mu

    dim(COV.mu) <- c(2*M,2*M)
    dim(DOF.mu) <- c(2*M,2*M)
  }

  # restructure indices from x*n,y*m to x,m,n,y !! CONSIDER DROPPING THIS? !!
  dim(COV.mu) <- c(AXES,M,AXES,M)
  dim(DOF.mu) <- c(AXES,M,AXES,M)
  COV.mu <- aperm(COV.mu,c(1,2,4,3))
  DOF.mu <- aperm(DOF.mu,c(1,2,4,3))

  # should I drop the indices in COV.mu and DOF.mu if possible ?
  COV.mu <- drop(COV.mu)
  DOF.mu <- drop(DOF.mu)

  # isotropic reduction if possible
  if(length(dim(DOF.mu))==2 && AXES==2 && DOF.mu[1,1]==DOF.mu[2,2] && DOF.mu[1,2]==0) { DOF.mu <- mean(diag(DOF.mu)) }
  # is this always correct?

  # re-write all of this to calculate constant, divide constant by n, and then subtract off from sum term by term?
  # likelihood constant/n: 2pi from det second term from variance-profiled quadratic term (which we will subtract if variance is not profiled)
  if(REML)
  { LL.CONST <- -(n-M)/n*AXES/2*log(2*pi) - AXES/2 } # not fixing this second term for REML yet... not wrong, but suboptimal maybe }
  else # ML constant
  { LL.CONST <- -AXES/2*log(2*pi) - AXES/2 }

  ### loglike: ( quadratic term of loglikelihood first )
  if(PROFILE && profile) { RATIO <- 1/VAR.MULT } # everything could be profiled - or transformed to variances whose mean could be profiled
  else if(!PROFILE) { RATIO <- mean(diag(cbind(R.sigma))) } # could profile anything and didn't # residuals fully standardized by model
  # PROFILE && !profile was filter specific !
  # finish off the loglikelihood calculation # RATIO is variance ratio from above
  loglike <- -AXES/2*(RATIO-1) - logdetCOV/2 #
  loglike <- n * (loglike + (LL.CONST-zero/n)) # I expect the last part (all constants) to mostly cancel out
  logdetCOV <- n * logdetCOV

  # mean structure terms
  if(REML) { loglike <- loglike + CTMM$REML.loglike + logdetcov/2 }

  if(verbose)
  {
    # assign variables
    CTMM$sigma <- sigma
    CTMM <- ctmm.repair(CTMM,K=K)

    if(UERE>=3) { CTMM$error <- as.logical(CTMM$error) } # old profile code broke logical error flag - not sure if still necessary
    else if(UERE && UERE<3 && PROFILE) { CTMM$error <- CTMM$error * sqrt(attr(sigma,'par')['area']) } # restore error ratio

    # mu <- drift@shift(mu,mu.center) # translate back to origin from center
    CTMM$mu <- mu
    CTMM$COV.mu <- COV.mu
    CTMM$DOF.mu <- DOF.mu

    CTMM$loglike <- loglike
    attr(CTMM,"info") <- attr(data,"info")

    CTMM <- ctmm.ctmm(CTMM)
    return(CTMM)
  }
  else
  { return(loglike) }
}


# environment for storing MLE when using pREML/pHREML/HREML
MLE.env <- new.env()
empty.env(MLE.env) # default to empty


# smallest resolutions in data (for soft bounding parameter guesstimates)
telemetry.mins <- function(data,axes=c('x','y'))
{
  dt <- stats::median(diff(data$t)) # median time difference

  df <- 2*pi/(last(data$t)-data$t[1]) # smallest frequency

  dz <- get.telemetry(data,axes)
  dz <- apply(dz,2,diff)
  dz <- rowSums( dz^2 )
  dz <- sqrt(min(dz[dz>0])) # smallest nonzero distance

  return(list(dt=dt,df=df,dz=dz))
}

###########################################################
# FIT MODEL WITH LIKELIHOOD FUNCTION (convenience wrapper to optim)
ctmm.fit <- function(data,CTMM=ctmm(),method="ML",COV=TRUE,control=list(),trace=FALSE)
{
  axes <- CTMM$axes
  # standardize data for numerical stability
  # pre-centering the data reduces later numerical error across models (tested)
  SHIFT <- colMeans(get.telemetry(data,axes))
  data[,axes] <- t( t(data[,axes,drop=FALSE]) - SHIFT )
  # add mu.center back to the mean value after kalman filter / mean profiling
  # pre-standardizing the data should also help
  SCALE <- sqrt(mean(get.telemetry(data,axes)^2))
  data <- unit.telemetry(data,length=SCALE)
  CTMM <- unit.ctmm(CTMM,length=SCALE)

  # used for minimum scale of parameter inspection
  n <- length(data$t)
  MINS <- telemetry.mins(data,axes)
  dt <- MINS$dt
  df <- MINS$df
  dz <- MINS$dz

  # unstandardize (includes likelihood adjustment)
  unscale.ctmm <- function(CTMM)
  {
    CTMM <- unit.ctmm(CTMM,length=1/SCALE)
    # log-likelihood adjustment
    CTMM$loglike <- CTMM$loglike - length(axes)*n*log(SCALE)
    # translate back to origin from center
    CTMM$mu <- drift@shift(CTMM$mu,SHIFT)

    return(CTMM)
  }

  method <- match.arg(method,c("ML","pREML","pHREML","HREML","REML"))

  default <- list(method="Nelder-Mead",precision=1/2,maxit=.Machine$integer.max)
  control <- replace(default,names(control),control)
  precision <- control$precision
  optimizer <- control$method
  control$method <- NULL

  if(method=="REML") { REML <- TRUE }
  else { REML <- FALSE }

  # clean/validate
  CTMM <- ctmm.ctmm(CTMM)
  drift <- get(CTMM$mean)
  CTMM$mu <- NULL # can always profile mu analytically
  range <- CTMM$range
  if(is.null(CTMM$sigma)) { CTMM$sigma <- covm(stats::cov(get.telemetry(data,axes)),isotropic=CTMM$isotropic,axes=axes) }

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

  # id and characterize parameters for profiling
  pars <- NAMES <- parscale <- lower <- upper <- period <- NULL
  setup.parameters <- function(CTMM,profile=TRUE,linear=FALSE)
  {
    STUFF <- id.parameters(CTMM,profile=profile,linear=linear,UERE=UERE,dt=dt,df=df,dz=dz)
    NAMES <<- STUFF$NAMES
    parscale <<- STUFF$parscale
    lower <<- STUFF$lower
    upper <<- STUFF$upper
    period <<- STUFF$period
    # initial guess for optimization
    pars <<- get.parameters(CTMM,NAMES)
  }
  setup.parameters(CTMM)
  # fix numeric error when it should be logical
  if(!("error" %in% NAMES)) { CTMM$error <- as.logical(CTMM$error) }

  # degrees of freedom, including the mean, variance/covariance, tau, and error model
  k.mean <- ncol(CTMM$mean.vec)

  # OPTIMIZATION FUNCTION (fn)
  # optional argument lengths: TAU, TAU+1, TAU+SIGMA
  fn <- function(p,zero=0)
  {
    names(p) <- NAMES
    p <- clean.parameters(p)
    CTMM <- set.parameters(CTMM,p)

    # negative log likelihood
    return(-ctmm.loglike(data,CTMM,REML=REML,zero=-zero,profile=profile))
  }

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

  # NOW OPTIMIZE
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
      CTMM$sigma <- VAR.MULT * CTMM$sigma
      CTMM$sigma@par["area"] <- VAR.MULT * CTMM$sigma@par["area"]
      CTMM$COV.mu <- VAR.MULT * CTMM$COV.mu
    }

    if(method=="pREML") { REML <- TRUE } # uses REML COV formula

    # fast calculation of sigma covariance
    COVSTUFF <- COV.covm(CTMM$sigma,n=n,k=k.mean,REML=REML)
    CTMM$COV <- COVSTUFF$COV
  }
  else # all further cases require optimization
  {
    if(trace) { message("Maximizing likelihood.") }
    control$covariance <- covariance()
    control$parscale <- parscale
    control$zero <- TRUE
    RESULT <- Optimizer(par=pars,fn=fn,method=optimizer,lower=lower,upper=upper,period=period,control=control)
    pars <- clean.parameters(RESULT$par)
    # copy over hessian from fit to COV.init ?

    # write best estimates over initial guess
    store.pars <- function(pars,profile=TRUE,finish=TRUE)
    {
      names(pars) <- NAMES
      pars <- clean.parameters(pars)

      CTMM <<- set.parameters(CTMM,pars)

      # this is a wasted evaluation !!! store verbose glob in environment?
      if(finish) { CTMM <<- ctmm.loglike(data,CTMM,REML=REML,verbose=TRUE,profile=profile) }
    }
    store.pars(pars,finish=TRUE)

    profile <- FALSE # no longer solving covariance analytically
    setup.parameters(CTMM,profile=FALSE)

    ### COV CALCULATION #############
    if(COV || method %in% c("pREML","pHREML","HREML"))
    {
      if(trace) { message("Calculating Hessian.") }
      DIFF <- genD(par=pars,fn=fn,zero=-CTMM$loglike,lower=lower,upper=upper,parscale=parscale,Richardson=2,mc.cores=1)
      hess <- DIFF$hessian
      grad <- DIFF$gradient

      # more robust covariance calculation than straight inverse
      CTMM$COV <- cov.loglike(hess,grad)
      dimnames(CTMM$COV) <- list(NAMES,NAMES)
    }

    # store MLE for faster model selection (ML is what is optimized, not pREML or HREML)
    if(method %in% c("pREML","pHREML","HREML"))
    {
      assign("EMPTY",FALSE,pos=MLE.env)
      assign("MLE",unscale.ctmm(CTMM),pos=MLE.env) # convert units back and store
    }
    else
    { empty.env(MLE.env) }

    # pREML correction ############################### only do pREML if sufficiently away from boundaries
    if(method %in% c("pREML","pHREML") && mat.min(hess) > .Machine$double.eps*length(NAMES))
    {
      # parameter correction
      REML <- TRUE
      #ML.grad <- grad # save old ML gradient
      if(trace) { message("Calculating REML gradient.") }
      DIFF <- genD(par=pars,fn=fn,zero=-CTMM$loglike,lower=lower,upper=upper,parscale=parscale,Richardson=2,order=1,mc.cores=1)

      # trying to make this robust here
      # COV is -1/Hessian, grad is of -loglike
      # using least-squares solution in the case of degeneracy (under parscale natural units)
      # hess <- t(hess*parscale)*parscale
      # d.pars <- -c(PDsolve(t(hess)%*%hess,pseudo=TRUE) %*% t(hess)%*%(DIFF$gradient*parscale)) * parscale
      # OK, that didn't work so well... sticking with this
      d.pars <- -c(CTMM$COV %*% DIFF$gradient)

      # increment transformed parameters
      # pars <- pars + d.pars
      # safety catch for bad models near boundaries
      pars <- line.boxer(d.pars,pars,lower=lower,upper=upper,period=period)
      names(pars) <- NAMES

      # store parameter correction only if a correction was made
      if(method=="pREML")
      {
        profile <- FALSE
        store.pars(pars,profile=FALSE,finish=TRUE)
      }
    }
    else if(method %in% c("pREML",'pHREML'))
    {
      warning("pREML failure: indefinite ML Hessian.")
      if(method=='pREML') { method <- 'ML' }
      else if(method=='pHREML') { method <- 'HREML' }
    }

    # profile REML parameters
    if(method %in% c('pHREML','HREML'))
    {
      REML <- TRUE
      profile <- TRUE
      store.pars(pars,profile=TRUE,finish=TRUE)

      # profile REML linear parameters numerically if necessary (error || circle)
      setup.parameters(CTMM,profile=TRUE,linear=TRUE)
      if(length(NAMES))
      {
        REML <- TRUE
        if(trace) { message("Profiling REML likelihood.") }
        control$covariance <- covariance()
        control$parscale <- parscale
        control$zero <- TRUE
        RESULT <- Optimizer(par=pars,fn=fn,method=optimizer,lower=lower,upper=upper,period=period,control=control)
        pars <- clean.parameters(RESULT$par)

        store.pars(pars,profile=TRUE,finish=TRUE)
      }
    }

    # FINAL COVARIANCE ESTIMATE
    TEST <- method %in% c('pREML','pHREML','HREML')
    if(TEST && COV) ### CALCULATE COVARIANCE MATRIX ###
    {
      profile <- FALSE
      setup.parameters(CTMM,profile=FALSE)

      if(trace) { message("Calculating REML Hessian.") }
      # calcualte REML Hessian at pREML parameters
      DIFF <- genD(par=pars,fn=fn,zero=-CTMM$loglike,lower=lower,upper=upper,parscale=parscale,Richardson=2,mc.cores=1)
      # Using MLE gradient, which should be zero off boundary
      CTMM$COV <- cov.loglike(DIFF$hessian,grad)
    }
    else if(TEST) ### don't confuse the ML COV with pREML COV
    { CTMM$COV <- NULL }

    if(COV) { dimnames(CTMM$COV) <- list(NAMES,NAMES) }
  } # end optimized estimates

  # model likelihood
  if(method!='ML') { CTMM$loglike <- ctmm.loglike(data,CTMM=CTMM,REML=FALSE,profile=FALSE) }
  CTMM$method <- method

  # covariance parameters only
  setup.parameters(CTMM,profile=FALSE,linear=FALSE)

  # unstandardize (includes likelihood adjustment)
  CTMM <- unscale.ctmm(CTMM)

  nu <- length(NAMES)
  # all parameters
  q <- length(axes)
  if(!range) { k.mean <- k.mean - 1 }
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

  # Mean square prediction error
  mspe <- function(K=1)
  {
    if(!CTMM$range && K==1) { return(Inf) }

    # velocity MSPE
    if(K==2)
    {
      if(length(CTMM$tau)<2 || any(CTMM$tau<=0)) { return(Inf) }
      UU <- VV
    }

    MSPE <- CTMM$COV.mu
    if(length(dim(MSPE))==2 && length(UU)==1) # multiple spatial dimensions and one trend component
    { MSPE <- sum(diag(MSPE)) * c(UU) }
    else if(length(dim(MSPE))==2 && length(UU)>1) # one spatial dimension and many trend components
    { MSPE <- sum(diag(MSPE %*% UU)) }
    else if(length(dim(MSPE))==4) # k trend components in multiple dimensions
    {
      MSPE <- lapply(1:length(axes),function(i) MSPE[i,,,i])
      MSPE <- Reduce("+",MSPE)
      MSPE <- sum(diag(MSPE %*% UU))
    }

    VAR <- sum(diag(CTMM$sigma))
    if(K==2)
    {
      STUFF <- get.taus(CTMM)
      VAR <- VAR * (STUFF$Omega2 + CTMM$circle^2)
    }

    MSPE <- MSPE + VAR

    return(MSPE)
  }
  STUFF <- drift@energy(CTMM)
  UU <- STUFF$UU
  VV <- STUFF$VV
  # Mean square prediction error in locations & velocities
  CTMM$MSPE <- c( mspe(K=1) , mspe(K=2) )
  names(CTMM$MSPE) <- c("position","velocity")

  return(CTMM)
}


####### calculate variance and variance-covaraince from area and eccentricity information
area2var <- function(CTMM,MEAN=TRUE)
{
  VAR <- mean(diag(CTMM$sigma))
  COV <- CTMM$COV

  NAMES <- rownames(COV)
  if(!CTMM$isotropic)
  {
    NAMES <- c("variance",NAMES[-(1:3)])
    area <- CTMM$sigma@par["area"]
    ecc <- CTMM$sigma@par["eccentricity"]

    # convert area, eccentricity uncertainty into mean variance uncertainty
    grad <- rbind( c(cosh(ecc/2),area*sinh(ecc/2)/2,0) )
    if(!MEAN) { grad <- 2*grad } # total x-y variance or average x-y variance

    P <- nrow(COV)
    if(P>3)
    {
      grad <- rbind( grad , array(0,c(P-3,3)) )
      grad <- cbind( grad , rbind( rep(0,P-3) , diag(1,P-3) ) )
    }

    COV <- grad %*% COV %*% t(grad)
    # backup for infinite covariances
    for(i in 1:nrow(COV))
    {
      if(any( is.nan(COV[i,]) | is.nan(COV[,i]) ))
      {
        COV[i,] <- COV[,i] <- 0
        COV[i,i] <- Inf
      }
    }
  }
  else
  { NAMES <- c("variance",NAMES[-1]) }
  dimnames(COV) <- list(NAMES,NAMES)

  return(COV)
}


###################
# general parameter guessing function
###################
ctmm.guess <- function(data,CTMM=ctmm(),variogram=NULL,name="GUESS",interactive=TRUE)
{
  # use intended axes
  if(is.null(variogram)) { variogram = variogram(data,axes=CTMM$axes) }
  else { CTMM$axes <- attr(variogram,"info")$axes }

  # mean specific guesswork/preparation
  drift <- get(CTMM$mean)
  CTMM <- drift@init(data,CTMM)
  mu <- CTMM$mu
  u <- drift(data$t,CTMM)

  # estimate circulation period if circle=TRUE
  if(CTMM$circle && class(CTMM$circle)=="logical")
  {
    n <- length(data$t)

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
