####################################
# log likelihood function
####################################
ctmm.loglike <- function(data,CTMM=ctmm(),REML=FALSE,profile=TRUE,zero=0,verbose=FALSE)
{
  # fail state - bad parameters or bad data
  if(verbose)
  {
    FAIL <- CTMM
    FAIL$loglike <- -Inf
  }
  else
  { FAIL <- -Inf }

  n <- length(data$t)
  AXES <- length(CTMM$axes)

  # save original tau length
  K <- length(CTMM$tau)
  # prepare model for numerics
  CTMM <- ctmm.prepare(data,CTMM)

  # drift <- get(CTMM$mean)

  range <- CTMM$range
  isotropic <- CTMM$isotropic

  ERROR <- CTMM$error

  COVM <- function(x)  # clean and enforce model structure
  {
    # 0/0 == Identity
    if(length(x)==1)
    { x <- nant(x,1) }
    else # length-4
    {
      SUB <- c(1,4)
      x[SUB] <- nant(x[SUB],1)
      SUB <- c(2,3)
      x[SUB] <- nant(x[SUB],0)
    }
    covm(x,isotropic=CTMM$isotropic,axes=CTMM$axes)
  }
  # sigma is current estimate and M.sigma is what we compute the Kalman filter with
  if(!is.null(CTMM$sigma))
  { sigma <- CTMM$sigma <- M.sigma <- COVM(CTMM$sigma) }
  else # better be profiling
  { sigma <- M.sigma <- COVM(diag(AXES)) }
  theta <- sigma@par["angle"] # NA in 1D

  STUFF <- squeezable.covm(CTMM)
  smgm <- STUFF$fact  # ratio of major axis to geometric mean axis
  ECC.EXT <- !STUFF$able # extreme eccentricity --- cannot squeeze data to match variances

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
  if(circle && ECC.EXT) { return(FAIL) } # can't squeeze !!! need 2D Langevin code

  n <- length(data$t)

  t <- data$t
  dt <- c(Inf,diff(t)) # time lags

  if(range) # timescale constant for profiling
  { DT <- 1 }
  else # IOU & BM
  {
    DT <- dt[-1] # not Inf, please
    DT <- DT[DT>0]

    if(length(DT))
    { DT <- stats::median(DT) }
    else
    { DT <- 1 }
  }

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
  class <- CTMM$class.mat
  ELLIPSE <- attr(error,"ellipse") # do we need error ellipses?
  # are we fitting the error?, then the above is not yet normalized
  TYPE <- DOP.match(CTMM$axes)
  UERE.RMS <- attr(data,"UERE")$UERE[,TYPE]
  UERE.DOF <- attr(data,"UERE")$DOF[,TYPE]
  names(UERE.DOF) <- names(UERE.RMS) <- rownames(attr(data,"UERE")$DOF) # R drops dimnames

  # only a subset of levels are in the data
  if(!is.null(names(CTMM$error)) && "class" %in% names(data))
  {
    LEVELS <- levels(data$class)
    UERE.RMS <- UERE.RMS[LEVELS]
    UERE.DOF <- UERE.DOF[LEVELS]
  }

  ## I don't recall what this was for, you can't profile from zero variance
  # if(is.null(CTMM$errors)) { CTMM$errors <- any(CTMM$error>0) }
  # UERE.FIT <- (CTMM$error | CTMM$errors) & !is.na(UERE.DOF) & UERE.DOF<Inf # will we be fitting error parameters?
  # UERE.FIX <- (CTMM$error | CTMM$errors) & (is.na(UERE.DOF) | UERE.DOF==Inf) # are there fixed error parameters
  UERE.FIT <- (CTMM$error) & !is.na(UERE.DOF) & UERE.DOF<Inf # will we be fitting error parameters?
  UERE.FIX <- (CTMM$error) & (is.na(UERE.DOF) | UERE.DOF==Inf) # are there fixed error parameters

  ### what kind of profiling is possible
  if((!any(CTMM$error>0) && !(circle && !isotropic)) || (!any(UERE.FIX) && isotropic)) # can profile full covariance matrix all at once
  { PROFILE <- 2 }
  else if(!any(UERE.FIX)) # can profile (max) variance with fixed-eccentricity 2D or 2x1D filters
  { PROFILE <- TRUE }
  else # no profiling possible - need exact-covariance filter (e.g. ELLIPSE FIXED)
  { PROFILE <- FALSE }

  ### 2D or 1D Kalman filters necessary?
  if(!circle && !isotropic && any(CTMM$error>0) && !ELLIPSE && !ECC.EXT) # can run 2x1D Kalman filters
  { DIM <- 1/2 }
  else if(ELLIPSE || (circle && !isotropic && any(CTMM$error>0)) || (ECC.EXT && AXES>1)) # full 2D filter necessary
  { DIM <- 2 }
  else # can run 1x1D Kalman filter
  { DIM <- 1 }

  # orient the data along the major & minor axes of sigma to run 2x1D filters - or to do basic circulation transformation after squeezing
  ROTATE <- !isotropic && (circle || (any(CTMM$error>0) && !ELLIPSE) && !ECC.EXT)
  if(ROTATE)
  {
    R <- rotate(-theta)
    z <- z %*% t(R)
    M.sigma <- rotate.covm(M.sigma,-theta)

    if(ELLIPSE) { error <- rotate.mat(error,-theta) } # rotate error ellipses
  }

  # squeeze from ellipse to circle for circulation model transformation
  # don't squeeze if some variances are zero
  SQUEEZE <- circle && !isotropic && !ECC.EXT
  if(SQUEEZE)
  {
    z <- squeeze(z,smgm)
    M.sigma <- squeeze.covm(M.sigma,circle=TRUE)

    # squeeze error circles into ellipses
    if(any(CTMM$error>0)) { error <- squeeze.mat(error,smgm) }

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
    if(ELLIPSE || (any(CTMM$error>0) && SQUEEZE)) { error <- rotates.mat(error,R) } # rotate error ellipses
    rm(R)
  }

  # Normalization of filters
  K.sigma <- M.sigma # default Kalman filter
  PRO.VAR <- 1 # value of variance multiplier to profile (default none)
  if(PROFILE==2 && !any(CTMM$error>0)) # unit-covariance filter # no error
  {
    UNIT <- 2

    # unit-COV KF sigma
    K.sigma <- COVM( diag(AXES) )

    # arg COV is ML COV matrix
    update.vars <- function(COV)
    {
      COV <- VAR.MULT*COV
      M.sigma <<- COVM(COV)
    }
  }
  else if(PROFILE) # unit-max-variance filters
  {
    UNIT <- 1

    # largest variance (to profile) --- used to be max, now really mean
    PRO.VAR <- profiled.var(CTMM,M.sigma,CTMM$error,DT=DT,AVE=TRUE)

    # unit-max-variance KF sigma # 1/DT for IOU/BM
    K.sigma <- scale.covm(M.sigma,1/PRO.VAR)

    if(any(UERE.FIT)) # unit-max-variance KF error
    { CTMM$error[UERE.FIT] <- CTMM$error[UERE.FIT] / sqrt(PRO.VAR) }

    # arg COV contains ML PRO.VAR estimate
    update.vars <- function(COV) # update PRO.VAR
    {
      COV <- mean(diag(cbind(COV))) # diag is annoying
      # relative in ratio to PRO.VAR
      if(any(UERE.FIT))
      { PRO.VAR <<- ( N*COV + sum( (UERE.DOF*(UERE.RMS/CTMM$error)^2)[UERE.FIT] ) )/(DOF+sum(UERE.DOF[UERE.FIT])) }
      else
      { PRO.VAR <<- VAR.MULT * c(COV) }


      # update sigma # M.sigma was in ratio to PRO.VAR
      M.sigma <<- scale.covm(K.sigma,PRO.VAR)
    }
  }
  else # fixed-variance filters
  {
    UNIT <- FALSE
    update.vars <- function(COV) { }
  }

  ### calibrate unknown errors given PROFILE state
  if(any(UERE.FIT)) # calibrate errors
  {
    class <- c( class %*% CTMM$error^2 )
    error[] <- class * error
  }
  rm(class)

  # check for bad time intervals # after evaluating location classes
  ZERO <- which(dt<=.Machine$double.eps)
  if(length(ZERO) && length(CTMM$tau) && CTMM$tau[1])
  {
    if(all(CTMM$error==FALSE)) { warning("Duplicate timestamps require an error model.") ; return(FAIL) }
    # check for HDOP==0 just in case
    ZERO <- error[ZERO,,,drop=FALSE]
    ZERO <- apply(ZERO,1,det) # AXES factors in product
    ZERO <- min(ZERO)
    if(ZERO<=.Machine$double.eps^AXES) { warning("Duplicate timestamps require an error model.") ; return(FAIL) }
  }

  # check for bad variances
  TEST <- eigenvalues.covm(sigma)
  if(min(TEST)/max(TEST,CTMM$error^2,1)<=.Machine$double.eps)
  {
    ZERO <- apply(error,1,det) # AXES factors in product
    ZERO <- min(ZERO)
    if(ZERO<=.Machine$double.eps^AXES) { return(FAIL) }
  }

  #### RUN KALMAN FILTERS ###
  if(DIM==1/2) ### 2 separate 1D Kalman filters instead of 1 2D Kalman filter ###
  {
    SIGMA <- eigenvalues.covm(K.sigma) # (relative to PRO.VAR)

    # major axis likelihood
    CTMM$sigma <- SIGMA[1]
    KALMAN1 <- kalman(cbind(z[,1]),u,dt=dt,CTMM=CTMM,error=error[,1,1,drop=FALSE]) # errors are relative to PRO.VAR if PROFILE

    # minor axis likelihood
    CTMM$sigma <- SIGMA[2]
    KALMAN2 <- kalman(cbind(z[,2]),u,dt=dt,CTMM=CTMM,error=error[,2,2,drop=FALSE]) # errors are relative to PRO.VAR if PROFILE

    mu <- cbind(KALMAN1$mu,KALMAN2$mu)

    R.sigma <- c(KALMAN1$sigma + KALMAN2$sigma)/2 # residual variance (relative to M.sigma)
    if(profile) { update.vars(R.sigma) } # update VARs before COV.mu!

    # combine uncorrelated estimates
    COV.mu <- array(c(KALMAN1$iW,diag(0,M),diag(0,M),KALMAN2$iW),c(M,M,2,2)) # -1/Hessian
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

    if(profile) { update.vars(R.sigma) }  # variance/covariance update based on residual covariance # update VARs before COV.mu!

    ### mu covariance terms
    if(DIM==2 || circle) # cases where we have a u(t) for each dimension
    {
      COV.mu <- KALMAN$iW # (2*M,2*M) from riffle

      logdetcov <- -log(det(KALMAN$W)) # log cov[beta] absolute # PRO.VAR handled below
    }
    else # lower or 1D filter return - cases where we have one u(t) for all dimensions
    {
      if(PROFILE) # unit covariance filter???
      { COV.mu <- KALMAN$iW %o% M.sigma } # (M,M,2,2)
      else if(!PROFILE)
      { COV.mu <- KALMAN$iW %o% diag(AXES) }

      # put either in first canonical form (AXES*M,AXES*M)
      COV.mu <- aperm(COV.mu,c(3,1,4,2))
      dim(COV.mu) <- c(AXES*M,AXES*M)

      logdetcov <- -AXES*log(det(KALMAN$W)) # log cov[beta] absolute
    }

    logdetCOV <- (AXES/DIM)*KALMAN$logdet # log autocovariance / n
  } ### END KALMAN FILTER RUNS ###

  if(PROFILE==1) { COV.mu <- COV.mu * PRO.VAR }
  COV.mu <- nant(COV.mu,0)

  # missing variances/covariances from profiling
  if(UNIT)
  {
    if(UNIT==1) { log.det.sigma <- AXES*log(PRO.VAR) } # unit-max-variance adjustment
    else if(UNIT==2) { log.det.sigma <- log(det.covm(M.sigma)) } # unit-COV adjustment

    logdetCOV <- logdetCOV + log.det.sigma # per n || n-1
    logdetcov <- logdetcov + M*log.det.sigma # absolute # !range handled below
  }

  # discard infinite prior uncertainty in stationary mean for BM/IOU
  if(!CTMM$range) { logdetcov <- ifelse(M==1,0, -log(det(COV.mu[-(1:AXES),-(1:AXES)])) ) }

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
  LL.CONST <- -AXES/2*log(2*pi) - AXES/2/VAR.MULT # why was VAR.MULT previously in the first term here?

  ## quadratic terms of log-likelihood
  if(profile && UNIT==2)
  { RATIO <- 1/VAR.MULT } # simple formula
  else if(UNIT==2) # M.sigma was only updated if profile
  { RATIO <- mean( diag( cbind(COVM(R.sigma) %*% solve.covm(M.sigma)) ) ) }
  else if(UNIT) # filter already partially standardized via eccentricity
  { RATIO <- mean(diag(cbind(R.sigma))) / PRO.VAR }
  else # couldn't profile anything and didn't # residuals standardized by model in KF
  { RATIO <- mean(diag(cbind(R.sigma))) }

  # else if(PROFILE && !profile) { filter specific above } !
  # finish off the loglikelihood calculation # RATIO is variance ratio from above
  loglike <- -AXES/2*(RATIO-1/VAR.MULT) - logdetCOV/2 # per n || n-1
  loglike <- N * (loglike + (LL.CONST-zero/N)) # I expect the last part (all constants) to mostly cancel out

  # mean structure terms
  if(REML) { loglike <- loglike + CTMM$REML.loglike + logdetcov/2 }

  if(any(UERE.FIT))
  {
    # update error from PROFILing
    if(PROFILE && !profile) # original error parameters
    { CTMM$error <- ERROR }
    else if(PROFILE && profile) # profiled error parameters
    { CTMM$error[UERE.FIT] <- CTMM$error[UERE.FIT] * sqrt(PRO.VAR) }

    # include calibration likelihood/prior
    SUB <- UERE.DOF>0 & UERE.FIT
    x <- (UERE.RMS/CTMM$error)[SUB]^2
    loglike <- loglike - AXES/2*sum(UERE.DOF[SUB]*(-log(x) + x - 1))
    # consider original calibration as log-likelihood as zero point
  }

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
