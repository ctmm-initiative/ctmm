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

  # drift <- get(CTMM$mean)

  range <- CTMM$range
  isotropic <- CTMM$isotropic

  COVM <- function(...) { covm(...,isotropic=CTMM$isotropic,axes=CTMM$axes) } # enforce model structure
  # sigma is current estimate and M.sigma is what we compute the Kalman filter with
  if(!is.null(CTMM$sigma))
  { sigma <- CTMM$sigma <- M.sigma <- COVM(CTMM$sigma) }
  else # better be profiling
  { sigma <- M.sigma <- COVM(diag(AXES)) }
  theta <- sigma@par["angle"] # NA in 1D

  STUFF <- squeezable.covm(CTMM)
  smgm <- STUFF$fact  # ratio of major axis to geometric mean axis
  ECC.EXT <- !STUFF$able # extreme eccentricity --- cannot squeeze data to match variances
  if(AXES==1) { ECC.EXT <- FALSE }

  # simplify filter for no movement process
  if(all(eigenvalues.covm(sigma)<=0))
  {
    CTMM$tau <- NULL
    CTMM$omega <- FALSE
    CTMM$circle <- FALSE
    isotropic <- TRUE
    CTMM$sigma <- COVM(0)
  }

  circle <- CTMM$circle
  if(circle && ECC.EXT) { return(-Inf) } # can't squeeze !!! need 2D Langevin code

  n <- length(data$t)

  t <- data$t
  dt <- c(Inf,diff(t)) # time lags

  # data z and mean vector u
  z <- get.telemetry(data,CTMM$axes)
  u <- CTMM$mean.vec
  M <- ncol(u) # number of linear parameters per spatial dimension

  if(range) { N <- n } else { N <- n-1 } # condition off first point
  # degrees of freedom
  if(REML) { DOF <- n-M } else { DOF <- N }
  # REML variance debias factor # ML constant
  VAR.MULT <- N/DOF

  # pre-centering the data reduces later numerical error across models (tested)
  # mu.center <- colMeans(z)
  # z <- t( t(z) - mu.center )
  # add mu.center back to the mean value after kalman filter / mean profiling
  # pre-standardizing the data would also help

  # get the error information
  error <- CTMM$error.mat # note for fitted errors, this is error matrix @ UERE=1 (CTMM$error)
  # are we fitting the error?, then the above is not yet normalized
  UERE <- attr(error,"flag")

  # check for bad time intervals
  ZERO <- which(dt==0)
  if(length(ZERO) && length(CTMM$tau) && CTMM$tau[1])
  {
    if(CTMM$error==FALSE) { warning("Duplicate timestamps require an error model.") ; return(-Inf) }
    # check for HDOP==0 just in case
    ZERO <- error[ZERO,,,drop=FALSE]
    ZERO <- apply(ZERO,1,det)
    if(any(ZERO<=0)) { warning("Duplicate timestamps require an error model.") ; return(-Inf) }
  }

  # check for bad variances
  if(min(eigenvalues.covm(sigma))<=.Machine$double.eps)
  {
    ZERO <- rep(FALSE,n)
    for(i in 1:dim(error)[2]) { ZERO <- ZERO | (error[,i,i]<=.Machine$double.eps) }
    if(any(ZERO)) { return(-Inf) }
  }

  ### what kind of profiling is possible
  if((!UERE && !(circle && !isotropic)) || (UERE<3 && isotropic)) # can profile full covariance matrix all at once
  { PROFILE <- 2 }
  else if(UERE<3) # can profile (max) variance with fixed-eccentricity 2D or 2x1D filters
  { PROFILE <- TRUE }
  else # no profiling possible - need exact-covariance filter (e.g. UERE>=3)
  { PROFILE <- FALSE }

  ### 2D or 1D Kalman filters necessary?
  if(!circle && !isotropic && UERE && UERE<4 && !ECC.EXT) # can run 2x1D Kalman filters
  { DIM <- 1/2 }
  else if(UERE==4 || (circle && !isotropic && UERE) || ECC.EXT) # full 2D filter necessary
  { DIM <- 2 }
  else # can run 1x1D Kalman filter
  { DIM <- 1 }

  # orient the data along the major & minor axes of sigma to run 2x1D filters - or to do basic circulation transformation after squeezing
  ROTATE <- !isotropic && (circle || (UERE && UERE<4) && !ECC.EXT)
  if(ROTATE)
  {
    R <- rotate(-theta)
    z <- z %*% t(R)
    M.sigma <- rotate.covm(M.sigma,-theta)

    if(UERE>=4) { error <- rotate.mat(error,-theta) } # rotate error ellipses
  }

  # squeeze from ellipse to circle for circulation model transformation
  # don't squeeze if some variances are zero
  SQUEEZE <- circle && !isotropic && !ECC.EXT
  if(SQUEEZE)
  {
    z <- squeeze(z,smgm)
    M.sigma <- squeeze.covm(M.sigma,circle=TRUE)

    # squeeze error circles into ellipses
    if(UERE) { error <- squeeze.mat(error,smgm) }

    # mean functions can requires squeezing (below)
    # how does !error && circle && !isotropic avoid this?
  }

  ## some cases need separate x & y mean functions
  if(DIM==2) # u -> cbind( (u,0) , (0,u) )
  {
    if(!SQUEEZE)
    { u <- vapply(1:ncol(u),function(i){cbind(u[,i],rep(0,n),rep(0,n),u[,i])},array(0,c(n,4))) }
    else
    { u <- vapply(1:ncol(u),function(i){cbind(u[,i]/smgm,rep(0,n),rep(0,n),u[,i]*smgm)},array(0,c(n,4))) }
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

  # largest variance (to profile)
  max.var <- 0
  update.max.var <- function(sigma)
  {
    max.var <<- mean(eigenvalues.covm(sigma))
    if(UERE<3) { max.var <<- max(max.var,CTMM$error^2) }
  }
  update.max.var(M.sigma)

  ### calibrate unknown errors given PROFILE state
  if(UERE && UERE<3) # calibrate errors
  {
    if(PROFILE) { CTMM$error <- CTMM$error / sqrt(max.var) } # fix error/variance ratio
    error <- CTMM$error^2 * error
  }

  # Normalization of filters
  K.sigma <- M.sigma # default Kalman filter
  if(PROFILE==2 && !UERE) # unit-covariance filter
  {
    UNIT <- 2
    K.sigma <- COVM( diag(AXES) )

    # arg COV is ML COV matrix
    update.vars <- function(COV)
    {
      COV <- VAR.MULT*COV
      M.sigma <<- COVM(COV)
      update.max.var(M.sigma)
    }
  }
  else if(PROFILE) # unit max-variance filters
  {
    UNIT <- 1
    K.sigma <- scale.covm(M.sigma,1/max.var)

    # arg COV contains ML max.var estimate
    update.vars <- function(COV)
    {
      COV <- mean(diag(cbind(COV))) # diag is annoying
      max.var <<- VAR.MULT * c(COV) # relative in ratio to max var !!!
      # update sigma # M.sigma was in ratio to max.var
      M.sigma <<- scale.covm(K.sigma,max.var)
      # update error
      if(UERE && UERE<3) { CTMM$error <<- CTMM$error * sqrt(max.var) }
      # update.max.var(sigma)
    }
  }
  else # fixed-variance filters
  {
    UNIT <- FALSE
    update.vars <- function(COV) { }
  }

  #### RUN KALMAN FILTERS ###
  if(DIM==1/2) ### 2 separate 1D Kalman filters instead of 1 2D Kalman filter ###
  {
    SIGMA <- eigenvalues.covm(K.sigma) # (relative to max.var)

    # major axis likelihood
    CTMM$sigma <- SIGMA[1]
    KALMAN1 <- kalman(cbind(z[,1]),u,dt=dt,CTMM=CTMM,error=error[,1,1,drop=FALSE]) # errors are relative to max.var if PROFILE

    # minor axis likelihood
    CTMM$sigma <- SIGMA[2]
    KALMAN2 <- kalman(cbind(z[,2]),u,dt=dt,CTMM=CTMM,error=error[,2,2,drop=FALSE]) # errors are relative to max.var if PROFILE

    mu <- cbind(KALMAN1$mu,KALMAN2$mu)

    R.sigma <- c(KALMAN1$sigma + KALMAN2$sigma)/2 # residual variance (relative to M.sigma)
    if(PROFILE && !profile)
    { RATIO <- R.sigma/max.var } # ratio of residual variance to model variance
    else if(profile)
    { update.vars(R.sigma) }

    # combine uncorrelated estimates
    COV.mu <- array(c(KALMAN1$iW,diag(0,M),diag(0,M),KALMAN2$iW),c(M,M,2,2)) # -1/Hessian
    if(PROFILE) { COV.mu <- COV.mu * max.var }

    # put in first canonical form (2*M,2*M)
    COV.mu <- aperm(COV.mu,c(3,1,4,2))
    dim(COV.mu) <- c(2*M,2*M)

    logdetCOV <- KALMAN1$logdet + KALMAN2$logdet
    logdetcov <- -log(det(KALMAN1$W)) - log(det(KALMAN2$W))
  }
  else ### 1x 1D or 2D Kalman filter ###
  {
    # prepare variance/covariance of 1D/2D Kalman filters
    if(PROFILE==2 || DIM==1) # could be rotated & squeezed
    { CTMM$sigma <- K.sigma@par['major'] }
    else if(DIM==2) # circle, !isotropic, UERE=1,2
    { CTMM$sigma <- K.sigma } # else sigma is full covariance matrix

    if(DIM==1) { error <- error[,1,1,drop=FALSE] } # isotropic && UERE redundant error information

    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)

    mu <- KALMAN$mu

    R.sigma <- KALMAN$sigma # residual covariance (relative to K.sigma)

    if(profile)
    { update.vars(R.sigma) } # variance/covariance update based on residual covariance
    else if(PROFILE==2)
    { RATIO <- mean( diag( COVM(R.sigma) %*% solve.covm(M.sigma) ) ) }
    else if(PROFILE) # filter already partially standardized via eccentricity
    { RATIO <- mean(diag(cbind(R.sigma))) / var.covm(M.sigma,ave=TRUE) }

    logdetCOV <- (AXES/DIM)*KALMAN$logdet # log autocovariance / n

    ### mu covariance terms
    if(DIM==2 || circle) # cases where we have a u(t) for each dimension
    {
      logdetcov <- -log(det(KALMAN$W)) # log cov[beta] absolute

      COV.mu <- KALMAN$iW # (2*M,2*M) from riffle
      # variance was not included
      if(PROFILE) { COV.mu <- COV.mu * max.var } # iW=COV/area

      logdetcov <- log(det(COV.mu)) # (2*M,2*M)
    }
    else # lower or 1D filter return - cases where we have one u(t) for all dimensions
    {
      logdetcov <- -AXES*log(det(KALMAN$W)) # log cov[beta] absolute

      if(PROFILE) # unit covariance filter???
      { COV.mu <- KALMAN$iW %o% M.sigma } # (M,M,2,2)
      else if(!PROFILE)
      { COV.mu <- KALMAN$iW %o% diag(AXES) }
      # put either in first canonical form (AXES*M,AXES*M)
      COV.mu <- aperm(COV.mu,c(3,1,4,2))
      dim(COV.mu) <- c(AXES*M,AXES*M)
    }
  }   ### END KALMAN FILTER RUNS ###
  COV.mu <- nant(COV.mu,0)

  # discard infinite prior uncertainty in stationary mean for BM/IOU
  if(!CTMM$range) { logdetcov <- ifelse(M==1,0, -log(det(COV.mu[-(1:AXES),-(1:AXES)])) ) }

  # missing variances/covariances from profiling
  if(UNIT)
  {
    if(UNIT==1) { log.det.sigma <- AXES*log(max.var) }
    else if(UNIT==2) { log.det.sigma <- log(det.covm(M.sigma)) }

    logdetCOV <- logdetCOV + log.det.sigma # per n
    logdetcov <- logdetcov + M*log.det.sigma # absolute
  }

  if(SQUEEZE) # de-squeeze
  {
    M.sigma <- squeeze.covm(M.sigma,smgm=1/smgm)

    if(DIM==1) # DIM==2 squeezed u(t), so beta and COV[beta] have normal scale already
    {
      R <- c(smgm,1/smgm)
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

    M.sigma <- rotate.covm(M.sigma,theta)

    dim(COV.mu) <- c(2,M*2*M)
    COV.mu <- R %*% COV.mu
    dim(COV.mu) <- c(2,M,2,M)
    COV.mu <- aperm(COV.mu,c(3,4,1,2))
    dim(COV.mu) <- c(2,M*2*M)
    COV.mu <- R %*% COV.mu
    dim(COV.mu) <- c(2*M,2*M)
  }

  # restructure indices from x*n,y*m to x,m,n,y
  dim(COV.mu) <- c(AXES,M,AXES,M)
  COV.mu <- aperm(COV.mu,c(1,2,4,3))
  # should I drop the mean indices in COV.mu and DOF.mu if possible ?
  if(M==1) { dim(COV.mu) <- c(AXES,AXES) }

  # re-write all of this to calculate constant, divide constant by n || (n-1), and then subtract off from sum term by term?
  # likelihood constant/n: 2pi from det second term from variance-profiled quadratic term (which we will subtract if variance is not profiled)
  LL.CONST <- -AXES/2*log(2*pi)/VAR.MULT - AXES/2
  # not fixing this second term for REML yet... not wrong, but suboptimal maybe

  ### loglike: ( quadratic term of loglikelihood first )
  if(PROFILE && profile) { RATIO <- 1/VAR.MULT } # everything could be profiled - or transformed to variances whose mean could be profiled
  else if(!PROFILE) { RATIO <- mean(diag(cbind(R.sigma))) } # couldn't profile anything and didn't # residuals fully standardized by model
  # PROFILE && !profile was filter specific !
  # finish off the loglikelihood calculation # RATIO is variance ratio from above
  loglike <- -AXES/2*(RATIO-1) - logdetCOV/2 #
  loglike <- N * (loglike + (LL.CONST-zero/N)) # I expect the last part (all constants) to mostly cancel out
  logdetCOV <- N * logdetCOV # what is this for?

  # mean structure terms
  if(REML) { loglike <- loglike + CTMM$REML.loglike + logdetcov/2 }
  loglike <- nant(loglike,-Inf) # Inf - Inf

  if(verbose)
  {
    # assign variables
    if(profile && PROFILE) { CTMM$sigma <- M.sigma }
    else { CTMM$sigma <- sigma }

    CTMM <- ctmm.repair(CTMM,K=K)

    # mu <- drift@shift(mu,mu.center) # translate back to origin from center
    CTMM$mu <- mu
    CTMM$COV.mu <- COV.mu

    CTMM$loglike <- loglike
    attr(CTMM,"info") <- attr(data,"info")

    CTMM <- ctmm.ctmm(CTMM)
    return(CTMM)
  }
  else
  { return(loglike) }
}


# smallest resolutions in data (for soft bounding parameter guesstimates)
telemetry.mins <- function(data,axes=c('x','y'))
{
  dt <- stats::median(diff(data$t)) # median time difference

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
  method <- match.arg(method,c("ML","pREML","pHREML","HREML","REML"))
  axes <- CTMM$axes
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

  if(method=="REML") { REML <- TRUE }
  else { REML <- FALSE }

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
    # initial guess for optimization
    pars <<- get.parameters(CTMM,NAMES,linear.cov=linear.cov)
  }
  setup.parameters(CTMM,profile=TRUE)
  if("error" %nin% NAMES) { CTMM$error <- as.logical(CTMM$error) } # fix numeric error when it should be logical

  # degrees of freedom, including the mean, variance/covariance, tau, and error model
  k.mean <- ncol(CTMM$mean.vec)

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
      CTMM$sigma <- VAR.MULT * CTMM$sigma
      CTMM$sigma@par["area"] <- VAR.MULT * CTMM$sigma@par["area"]
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
    RESULT <- optimizer(par=pars,fn=fn,method=op.method,lower=lower,upper=upper,period=period,control=control)
    pars <- clean.parameters(RESULT$par)
    # copy over hessian from fit to COV.init ?

    # write best estimates over initial guess
    store.pars <- function(pars,profile=TRUE,finish=TRUE)
    {
      names(pars) <- NAMES
      pars <- clean.parameters(pars,linear.cov=linear.cov)

      CTMM <<- set.parameters(CTMM,pars,linear.cov=linear.cov)

      # this is a wasted evaluation !!! store verbose glob in environment?
      if(finish) { CTMM <<- ctmm.loglike(data,CTMM,REML=REML,verbose=TRUE,profile=profile) }
    }
    store.pars(pars,finish=TRUE)

    profile <- FALSE # no longer solving covariance analytically
    setup.parameters(CTMM,profile=FALSE)
    ### COV CALCULATION #############
    if(COV || method %in% c("pREML","pHREML"))
    {
      # if pREML, calculate Hessian in safe parameterization and then transform afterwards
      if(trace) { message("Calculating Hessian.") }
      DIFF <- genD(par=pars,fn=fn,zero=-CTMM$loglike,lower=lower,upper=upper,parscale=parscale,Richardson=2,mc.cores=1)
      hess <- DIFF$hessian
      grad <- DIFF$gradient

      # more robust covariance calculation than straight inverse
      CTMM$COV <- cov.loglike(hess,grad)
      dimnames(CTMM$COV) <- list(NAMES,NAMES)
    }

    # store MLE for faster re-optimization #
    MLE <- unscale.ctmm(CTMM)

    ### pREML correction ########## only do pREML if sufficiently away from boundaries
    if(method %in% c("pREML","pHREML") && mat.min(hess) > .Machine$double.eps*length(NAMES))
    {
      # parameter correction
      REML <- TRUE
      #ML.grad <- grad # save old ML gradient
      if(trace) { message("Calculating REML gradient.") }
      DIFF <- genD(par=pars,fn=fn,lower=lower,upper=upper,parscale=parscale,Richardson=2,order=1,mc.cores=1)

      J <- diag(length(pars))
      dimnames(J) <- list(NAMES,NAMES)
      # invoke Jacobian (major,minor/major,angle) -> linear parameterization !!!
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

      # increment transformed parameters
      # pars <- pars + d.pars
      # safety catch for bad models near boundaries
      pars <- line.boxer(d.pars,pars,lower=lower,upper=upper,period=period)
      names(pars) <- NAMES

      # store parameter correction
      profile <- FALSE
      store.pars(pars,profile=FALSE,finish=TRUE)

      linear.cov <- FALSE
      setup.parameters(CTMM,profile=FALSE)
    }
    else if(method %in% c("pREML",'pHREML'))
    {
      warning("pREML failure: indefinite ML Hessian.")
      if(method=='pREML') { method <- 'ML' }
      else if(method=='pHREML') { method <- 'HREML' }
    }
    ### end pREML correction ###

    ### profile linear REML parameters ###
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
        RESULT <- optimizer(par=pars,fn=fn,method=op.method,lower=lower,upper=upper,period=period,control=control)
        pars <- clean.parameters(RESULT$par)
      }
      # includes free profile
      store.pars(pars,profile=TRUE,finish=TRUE)
    }

    ### FINAL COVARIANCE ESTIMATE ###
    if(COV && method %in% c('pREML','pHREML','HREML')) ### CALCULATE COVARIANCE MATRIX ###
    {
      profile <- FALSE
      setup.parameters(CTMM,profile=FALSE)

      if(trace) { message("Calculating REML Hessian.") }
      # calcualte REML Hessian at pREML parameters
      DIFF <- genD(par=pars,fn=fn,lower=lower,upper=upper,parscale=parscale,Richardson=2,mc.cores=1)
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
