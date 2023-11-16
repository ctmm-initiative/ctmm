get.link <- function(CTMM)
{
  link <- CTMM$link

  if(is.null(link)) { link <- "identity" }

  link <- list(name=link,fn=get(link))

  if(link$name=="identity")
  { link$grad <- function(x){rep(1,length(x))} }
  else if(link$name=="log")
  { link$grad <- function(x){1/x} }
  else # numDeriv
  {} # TODO

  return(link)
}

ctmm.circulate <- function(CTMM,t,dt=diff(t))
{
  dt <- c(0,dt)
  dynamics <- CTMM$dynamics

  if(is.null(dynamics) || dynamics==FALSE || dynamics=="stationary")
  { circle <- CTMM$circle * (t-t[1]) }
  else
  {
    CP <- CTMM[[dynamics]] # change points
    circle <- rep(0,length(t))
    j <- 1
    s <- as.character(CP$state[j])
    for(i in 1:length(t))
    {
      if(i>1) { circle[i] <- circle[i-1] }

      while(t[i]>CP$stop[j]) # do we cross a change point?
      {
        DT <- CP$stop[j] - max(t[i-1],CP$start[j]) # time until change
        dt[i] <- dt[i] - DT
        circle[i] <- circle[i] + CTMM[[s]]$circle * DT
        j <- j + 1
        s <- as.character(CP$state[j])
      }

      circle[i] <- circle[i] + CTMM[[s]]$circle * dt[i]
    }
  }

  return(circle)
}

ctmm.apply <- function(CTMM,fn=identity,states=get.states(CTMM))
{
  CTMM <- fn(CTMM)
  for(s in states) { CTMM[[s]] <- fn(CTMM[[s]]) }
  return(CTMM)
}

sigma.apply <- function(CTMM,fn=identity,states=get.states(CTMM))
{
  CTMM$sigma <- fn(CTMM$sigma)
  for(s in states) { CTMM[[s]]$sigma <- fn(CTMM[[s]]$sigma) }
  return(CTMM)
}


####################################
# log likelihood function
####################################
ctmm.loglike <- function(data,CTMM=ctmm(),REML=FALSE,profile=TRUE,zero=0,verbose=FALSE,compute=TRUE,...)
{
  # fail state - bad parameters or bad data
  if(verbose)
  {
    FAIL <- CTMM
    FAIL$loglike <- -Inf
  }
  else
  { FAIL <- -Inf }

  # employ link function on time
  if(length(CTMM$timelink.par)) { data$t <- linktime(data,CTMM) }

  n <- length(data$t)
  AXES <- length(CTMM$axes)

  # save original tau length
  K <- length(CTMM$tau)
  # prepare model for numerics
  CTMM <- ctmm.prepare(data,CTMM,verbose=verbose)

  range <- CTMM$range
  isotropic <- CTMM$isotropic
  ERROR <- CTMM$error
  commute <- is.null(CTMM$commute) || ctmm.commute(CTMM)
  dynamics <- CTMM$dynamics
  states <- get.states(CTMM)

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

  # if(!drift.is.finite(CTMM,data)) { return(FAIL) }

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

  # get the error information
  error <- CTMM$error.mat # note for fitted errors, this is error matrix @ UERE=1 (CTMM$error)
  class <- CTMM$class.mat
  ELLIPSE <- attr(error,"ellipse") # do we need error ellipses?
  # are we fitting the error?, then the above is not yet normalized
  TYPE <- DOP.match(CTMM$axes)
  if(TYPE!="unknown")
  {
    UERE.RMS <- attr(data,"UERE")$UERE[,TYPE]
    UERE.DOF <- attr(data,"UERE")$DOF[,TYPE]
    names(UERE.DOF) <- names(UERE.RMS) <- rownames(attr(data,"UERE")$DOF) # R drops dimnames
  }
  else
  {
    UERE.RMS <- UERE.DOF <- 0
    names(UERE.RMS) <- names(UERE.DOF) <- "all"
  }

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
  if(commute && (!any(CTMM$error>0) && !(circle && !isotropic)) || (!any(UERE.FIX) && isotropic)) # can profile full covariance matrix all at once
  { PROFILE <- 2 }
  else if(!any(UERE.FIX)) # can profile (max) variance with fixed-eccentricity 2D or 2x1D filters
  { PROFILE <- TRUE }
  else # no profiling possible - need exact-covariance filter (e.g. ELLIPSE FIXED)
  { PROFILE <- FALSE }

  ### 2D or 1D Kalman filters necessary?
  if(commute && !circle && !isotropic && any(CTMM$error>0) && !ELLIPSE && !ECC.EXT) # can run 2x1D Kalman filters
  { DIM <- 1.5 }
  else if(!commute || ELLIPSE || (circle && !isotropic && any(CTMM$error>0)) || (ECC.EXT && AXES>1)) # full 2D filter necessary
  { DIM <- 2 }
  else # can run 1x1D Kalman filter
  { DIM <- 1 }

  # orient the data along the major & minor axes of sigma to run 2x1D filters - or to do basic circulation transformation after squeezing
  ROTATE <- commute && !isotropic && (circle || (any(CTMM$error>0) && !ELLIPSE) && !ECC.EXT)
  if(ROTATE)
  {
    R <- rotate(-theta)
    z <- z %*% t(R)
    fn <- function(sigma) { rotate.covm(sigma,-theta) }
    CTMM <- sigma.apply(CTMM,fn,states)
    M.sigma <- rotate.covm(M.sigma,-theta)

    if(ELLIPSE) { error <- rotate.mat(error,-theta) } # rotate error ellipses
  }

  # squeeze from ellipse to circle for circulation model transformation
  # don't squeeze if some variances are zero
  SQUEEZE <- commute && circle && !isotropic && !ECC.EXT
  if(SQUEEZE)
  {
    z <- squeeze(z,smgm)
    fn <- function(sigma) { squeeze.covm(sigma,circle=TRUE) }
    CTMM <- sigma.apply(CTMM,fn,states)
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
    # R <- circle*(t-t[1])
    R <- ctmm.circulate(CTMM,t,dt) # circle * (t-t[1])
    R <- rotates(-R) # rotation matrices
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
    scale.vars <- function(COV)
    {
      COV <- VAR.MULT*COV
      if(profile)
      {
        fn <- function(sigma) { COVM(COV*sigma) }
        CTMM <<- sigma.apply(CTMM,fn,states)
      }
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
    scale.vars <- function(COV) # update PRO.VAR
    {
      COV <- mean(diag(cbind(COV))) # diag is annoying
      # relative in ratio to PRO.VAR
      if(any(UERE.FIT))
      { PRO.VAR <<- ( N*COV + sum( (UERE.DOF*(UERE.RMS/CTMM$error)^2)[UERE.FIT] ) )/(DOF+sum(UERE.DOF[UERE.FIT])) }
      else
      { PRO.VAR <<- VAR.MULT * c(COV) }

      # update sigma # M.sigma was in ratio to PRO.VAR
      if(profile)
      {
        fn <- function(sigma) { scale.covm(sigma,PRO.VAR) }
        CTMM <<- sigma.apply(CTMM,fn,states)
      }
      M.sigma <<- scale.covm(K.sigma,PRO.VAR)
    }
  }
  else # fixed-variance filters
  {
    UNIT <- FALSE
    scale.vars <- function(COV) { }
  }
  KMR <- mean(diag(K.sigma/M.sigma))

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

  if(!compute) # return information on Kalman filter and parameters
  {
    R <- list()
    R$DIM <- DIM
    R$pars <- id.parameters(CTMM,UERE.FIT=UERE.FIT)$NAMES # numerically fitted parameters
    return(R)
  }

  #### RUN KALMAN FILTERS ###
  if(DIM==1.5) ### 2 separate 1D Kalman filters instead of 1 2D Kalman filter ###
  {
    SIGMA <- eigenvalues.covm(K.sigma) # (relative to PRO.VAR)

    # major axis likelihood
    fn <- function(sigma){ KMR*sigma[1] }
    CTMM <- sigma.apply(CTMM,fn,states)
    CTMM$sigma <- SIGMA[1] # should now be redundant
    KALMAN1 <- kalman(cbind(z[,1]),u,t=t,dt=dt,CTMM=CTMM,error=error[,1,1,drop=FALSE]) # errors are relative to PRO.VAR if PROFILE

    # minor axis likelihood
    fn <- function(sigma){ KMR*sigma[4] }
    CTMM <- sigma.apply(CTMM,fn,states)
    CTMM$sigma <- SIGMA[2] # should now be redundant
    KALMAN2 <- kalman(cbind(z[,2]),u,t=t,dt=dt,CTMM=CTMM,error=error[,2,2,drop=FALSE]) # errors are relative to PRO.VAR if PROFILE

    mu <- cbind(KALMAN1$mu,KALMAN2$mu)

    R.sigma <- c(KALMAN1$sigma + KALMAN2$sigma)/2 # residual variance (relative to M.sigma)
    if(profile) { scale.vars(R.sigma) } # update VARs before COV.mu!

    # combine uncorrelated estimates
    COV.mu <- array(c(KALMAN1$iW,diag(0,M),diag(0,M),KALMAN2$iW),c(M,M,2,2)) # -1/Hessian
    # put in first canonical form (2*M,2*M)
    COV.mu <- aperm(COV.mu,c(3,1,4,2))
    dim(COV.mu) <- c(2*M,2*M)

    logdetCOV <- KALMAN1$logdet + KALMAN2$logdet
    logdetcov <- -PDlogdet(KALMAN1$W) - PDlogdet(KALMAN2$W)
  }
  else ### 1x 1D or 2D Kalman filter ###
  {
    # prepare variance/covariance of 1D/2D Kalman filters
    if(PROFILE==2 || DIM==1) # could be rotated & squeezed
    {
      fn <- function(sigma) { var.covm(KMR*sigma,ave=TRUE) }
      CTMM <- sigma.apply(CTMM,fn,states)
      CTMM$sigma <- var.covm(K.sigma,ave=TRUE) # should be redundant
    }
    else if(DIM==2) # circle, !isotropic, UERE=1,2
    {
      fn <- function(sigma) { scale.covm(sigma,KMR) }
      CTMM <- sigma.apply(CTMM,fn,states)
      CTMM$sigma <- K.sigma # should be redundant
    } # else sigma is full covariance matrix

    if(DIM==1) { error <- error[,1,1,drop=FALSE] } # isotropic && UERE redundant error information

    KALMAN <- kalman(z,u,t=t,dt=dt,CTMM=CTMM,error=error,DIM=DIM)

    mu <- KALMAN$mu

    R.sigma <- KALMAN$sigma # residual covariance (relative to K.sigma)

    if(profile) { scale.vars(R.sigma) }  # variance/covariance update based on residual covariance # update VARs before COV.mu!

    ### mu covariance terms
    if(DIM==2 || circle) # cases where we have a u(t) for each dimension
    {
      COV.mu <- KALMAN$iW # (2*M,2*M) from riffle

      logdetcov <- -PDlogdet(KALMAN$W) # log cov[beta] absolute # PRO.VAR handled below
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

      logdetcov <- -AXES*PDlogdet(KALMAN$W) # log cov[beta] absolute
    }

    logdetCOV <- (AXES/DIM)*KALMAN$logdet # log autocovariance / n
  } ### END KALMAN FILTER RUNS ###

  if(PROFILE==1) { COV.mu <- COV.mu * PRO.VAR }
  COV.mu <- nant(COV.mu,0)

  # missing variances/covariances from profiling
  if(UNIT)
  {
    if(UNIT==1) { log.det.sigma <- AXES*log(PRO.VAR) } # unit-max-variance adjustment
    else if(UNIT==2) { log.det.sigma <- PDlogdet(M.sigma) } # unit-COV adjustment

    logdetCOV <- logdetCOV + log.det.sigma # per n || n-1
    logdetcov <- logdetcov + M*log.det.sigma # absolute # !range handled below
  }

  # discard infinite prior uncertainty in stationary mean for BM/IOU
  if(!length(states) && !CTMM$range)
  { logdetcov <- ifelse(M==1,0, -PDlogdet(COV.mu[-(1:AXES),-(1:AXES)]) ) }
  else if(length(states))
  {
    offset <- 0
    IND <- rep(TRUE,nrow(COV.mu))
    # fill matrix
    offset <- 0
    for(s in states)
    {
      if(!CTMM[[s]]$range) { IND[offset+1:AXES] <- FALSE }
      offset <- offset + nrow(CTMM[[s]]$mu)
    }
    if(any(IND))
    { logdetcov <- -PDlogdet(COV.mu[IND,IND]) }
    else
    { logdetcov <- 0 }
  }

  if(SQUEEZE) # de-squeeze
  {
    fn <- function(sigma) { squeeze.covm(sigma,smgm=1/smgm) }
    CTMM <- sigma.apply(CTMM,fn,states)
    M.sigma <- squeeze.covm(M.sigma,smgm=1/smgm) # should be redundant

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

    fn <- function(sigma) { rotate.covm(sigma,theta) }
    CTMM <- sigma.apply(CTMM,fn,states)
    M.sigma <- rotate.covm(M.sigma,theta) # should be redundant

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
  { RATIO <- mean( diag( cbind(COVM(R.sigma) %*% solve_covm(M.sigma)) ) ) }
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
    if(profile && PROFILE) { CTMM$sigma <- M.sigma } # states done in scale.vars()
    else { CTMM$sigma <- sigma }

    if(TYPE!="unknown") { CTMM$UERE <- attr(data,"UERE")$UERE[,TYPE] }

    CTMM <- ctmm.repair(CTMM,K=K)

    CTMM$mu <- mu
    CTMM$COV.mu <- COV.mu

    # store state means
    offset <- 0
    for(s in states)
    {
      IND <- offset + 1:nrow(CTMM[[s]]$mu)
      CTMM[[s]]$mu <- CTMM$mu[IND,]
      CTMM[[s]]$COV.mu <- CTMM$COV.mu[,IND,,IND]
      offset <- offset + nrow(CTMM[[s]]$mu)
    }

    CTMM$loglike <- loglike
    attr(CTMM,"info") <- attr(data,"info")

    CTMM <- ctmm.ctmm(CTMM)
    return(CTMM)
  }
  else
  { return(loglike) }
}
