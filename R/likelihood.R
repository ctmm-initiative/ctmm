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
  ZERO <- which(dt<=.Machine$double.eps)
  if(length(ZERO) && length(CTMM$tau) && CTMM$tau[1])
  {
    if(CTMM$error==FALSE) { warning("Duplicate timestamps require an error model.") ; return(-Inf) }
    # check for HDOP==0 just in case
    ZERO <- error[ZERO,,,drop=FALSE]
    ZERO <- apply(ZERO,1,det)
    ZERO <- min(ZERO) * CTMM$error^4
    if(ZERO<=.Machine$double.eps^2) { warning("Duplicate timestamps require an error model.") ; return(-Inf) }
  }

  # check for bad variances
  TEST <- eigenvalues.covm(sigma)
  if(min(TEST)/max(TEST,CTMM$error^2,1)<=.Machine$double.eps)
  {
    ZERO <- apply(error,1,det)
    ZERO <- min(ZERO) * CTMM$error^4
    if(ZERO<=.Machine$double.eps^2) { return(-Inf) }
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

  # largest variance (to profile) --- used to be max, now really mean
  max.var <- 0
  update.max.var <- function(sigma)
  {
    max.var <<- mean(eigenvalues.covm(sigma)) # really profiling the variance with mean?
    if(UERE<3) { max.var <<- max.var + CTMM$error^2/AXES } # comparable error variance (@DOP==1)
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
    { CTMM$sigma <- var.covm(K.sigma,ave=TRUE) }
    else if(DIM==2) # circle, !isotropic, UERE=1,2
    { CTMM$sigma <- K.sigma } # else sigma is full covariance matrix

    if(DIM==1) { error <- error[,1,1,drop=FALSE] } # isotropic && UERE redundant error information

    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error,DIM=DIM)

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
