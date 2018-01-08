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
  CTMM <- ctmm.prepare(data,CTMM,REML=REML)

  drift <- get(CTMM$mean)

  range <- CTMM$range
  isotropic <- CTMM$isotropic

  sigma <- CTMM$sigma
  if(!is.null(sigma))
  {
    area <- sigma@par[1]
    ecc <- sigma@par[2]
    theta <- sigma@par[3]
  }
  else
  {
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
  mu.center <- colMeans(z)
  z <- t(z) - mu.center
  z <- t(z)
  # add mu.center back to the mean value after kalman filter / mean profiling

  # variance debias factor
  if(REML)
  { VAR.MULT <- n/(n-M) }
  else # ML constant
  { VAR.MULT <- 1 }

  # make the errors relative to profile the variance
  if(profile) { CTMM$error <- CTMM$error / sqrt(area) }
  # get the error information
  error <- get.error(data,CTMM)
  # are we fitting the error
  UERE <- attr(error,"flag")

  # do we need to orient the data along the major an minor axes of sigma
  ROTATE <- !isotropic && (CTMM$error || circle)
  if(ROTATE) { z <- z %*% t(rotate(-theta)) }

  if(circle) # ONE KALMAN FILTER WITH COMPLEX SIGNAL (will always be 2D)
  {
    # proportional standardization from ellipse to circle
    if(!isotropic && ecc)
    {
      z[,1] <- z[,1] * exp(-ecc/4)
      z[,2] <- z[,2] * exp(+ecc/4)
    }
    z <- cbind(z[,1] + 1i*z[,2])

    # corotating frame
    R <- exp(-1i*circle*(t-t[1]))
    z <- R * z
    u <- R * u

    if(UERE<3 && profile)
    {
      CTMM$sigma <- 1
      KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error) # error is relative

      ML.area <- as.numeric(Re(KALMAN$sigma))/2
      # profile variance if unspecified
      # debias if REML
      area <- VAR.MULT*ML.area

      sigma <- covm(c(area,ecc,theta),isotropic=isotropic,axes=CTMM$axes)

      # convert from error ratio back to absolute error
      CTMM$error <- sqrt(area) * CTMM$error

      COV.mu <- KALMAN$iW * area
      DOF.mu <- KALMAN$W

      # terms needed for AICc
      logdetCOV <- 2*KALMAN$logdet + 2*log(area) # per n

      # this is loglikelihood/n quadratic term
      loglike <- -(1/VAR.MULT-1)
    }
    else
    {
      CTMM$sigma <- area
      KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)

      # what does this term here represent? relative variance or something
      R.sigma <- Re(KALMAN$sigma)/2

      COV.mu <- KALMAN$iW
      DOF.mu <- area * KALMAN$W

      # terms needed for AICc
      logdetCOV <- 2*KALMAN$logdet # per n

      # this is loglikelihood/n quadratic term
      loglike <- -(R.sigma-1)
    }

    # real array formatting
    mu <- KALMAN$mu
    mu <- cbind(Re(mu),Im(mu))
    # complex correlations are x-y correlations
    COV.mu <- array(c(Re(COV.mu),-Im(COV.mu),Im(COV.mu),Re(COV.mu)),c(M,M,2,2))
    DOF.mu <- array(c(Re(DOF.mu),-Im(DOF.mu),Im(DOF.mu),Re(DOF.mu)),c(M,M,2,2))

    # de-standardization
    R <- exp(c(+1,-1)*ecc/4)
    mu <- t(R * t(mu))
    COV.mu <- array(COV.mu,c(M^2*2,2))
    COV.mu <- t(R * t(COV.mu))
    COV.mu <- array(COV.mu,c(M,M,2,2))
    COV.mu <- aperm(COV.mu,c(1,2,4,3))
    COV.mu <- array(COV.mu,c(M^2*2,2))
    COV.mu <- t(R * t(COV.mu))
    COV.mu <- array(COV.mu,c(M,M,2,2))

    # now we can calculate the determinant
    logdetcov <- aperm(COV.mu,c(1,3,2,4))
    logdetcov <- array(logdetcov,c(2*M,2*M))
    logdetcov <- log(det(logdetcov))
  } # NOT CIRCLE
  else if(!CTMM$error) # ONE KALMAN FILTER WITH UNKNOWN COVARIANCE AND NO ERROR
  {
    CTMM$sigma <- 1
    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)

    ML.sigma <- KALMAN$sigma
    # profile covariance if sigma unspecified
    if(profile)
    {
      # debias if REML
      sigma <- VAR.MULT*ML.sigma

      sigma <- covm(sigma,isotropic=isotropic,axes=CTMM$axes)
    }

    mu <- KALMAN$mu
    COV.mu <- KALMAN$iW %o% sigma
    DOF.mu <- KALMAN$W %o% diag(AXES)

    # terms needed for AICc
    logdetCOV <- AXES*KALMAN$logdet + log(det(sigma)) # per n
    logdetcov <- log(det(sigma)^M/det(KALMAN$W)^AXES)

    # this is loglikelihood/n quadratic term
    if(!profile)
    { loglike <- -sum(diag(ML.sigma %*% PDsolve(sigma))-1)/2 }
    else
    { loglike <- -AXES*(1/VAR.MULT-1)/2 }
  } # ERROR
  else if(isotropic && UERE<3 && profile) # ONE KALMAN FILTER WITH UNKNOWN VARIANCE AND KNOWN ERROR RATIO
  {
    CTMM$sigma <- 1
    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error) # error is relative

    ML.area <- mean(diag(KALMAN$sigma))
    # profile variance if unspecified
    # debias if REML
    area <- VAR.MULT*ML.area

    # convert from error ratio back to absolute error
    CTMM$error <- sqrt(area) * CTMM$error

    sigma <- covm(area,isotropic=isotropic,axes=CTMM$axes)

    mu <- KALMAN$mu
    COV.mu <- KALMAN$iW %o% sigma
    DOF.mu <- KALMAN$W %o% diag(AXES)

    logdetCOV <- AXES*KALMAN$logdet + log(det(sigma))
    logdetcov <- log(det(sigma)^M/det(KALMAN$W)^AXES)

    # this is loglikelihood/n quadratic term
    loglike <- -AXES*(1/VAR.MULT-1)/2
  }
  else if(isotropic) # ONE KALMAN FILTER WITH KNOWN VARIANCE AND KNOWN ERROR
  {
    CTMM$sigma <- area
    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)

    mu <- KALMAN$mu
    COV.mu <- KALMAN$iW %o% diag(AXES)
    DOF.mu <- (area * KALMAN$W) %o% diag(AXES)

    # includes error and standardization
    R.sigma <- KALMAN$sigma
    R.sigma <- mean(diag(R.sigma))

    logdetCOV <- AXES*KALMAN$logdet
    logdetcov <- -AXES*log(det(KALMAN$W))

    # this is loglikelihood/n quadratic term
    loglike <- -AXES*(R.sigma-1)/2
  }
  else if(UERE<3 && profile) # TWO KALMAN FILTERS WITH UNKNOWN VARIANCE AND KNOWN ERROR RATIO
  {
    # eigen variances proportionally
    SIGMA <- exp(c(+1,-1)*ecc/2)

    # major axis likelihood
    CTMM$sigma <- SIGMA[1]
    KALMAN1 <- kalman(cbind(z[,1]),u,dt=dt,CTMM=CTMM,error=error) # error is relative

    # minor axis likelihood
    CTMM$sigma <- SIGMA[2]
    KALMAN2 <- kalman(cbind(z[,2]),u,dt=dt,CTMM=CTMM,error=error) # error is relative

    ML.area <- c(KALMAN1$sigma + KALMAN2$sigma)/2
    # profile variance if unspecified
    # debias if REML
    area <- VAR.MULT*ML.area

    # convert from error ratio back to absolute error
    CTMM$error <- sqrt(area) * CTMM$error

    sigma <- covm(c(area,ecc,theta),isotropic=isotropic,axes=CTMM$axes)

    logdet <- KALMAN1$logdet + KALMAN2$logdet

    mu <- cbind(KALMAN1$mu,KALMAN2$mu)
    COV.mu <- area * array(c(KALMAN1$iW,diag(0,M),diag(0,M),KALMAN2$iW),c(M,M,2,2)) # -1/Hessian
    DOF.mu <- array(c(SIGMA[1]*KALMAN1$W,diag(0,M),diag(0,M),SIGMA[2]*KALMAN2$W),c(M,M,2,2))

    logdetCOV <- logdet + 2*log(area)
    logdetcov <- -log(det(KALMAN1$W/area)*det(KALMAN2$W/area))

    # this is loglikelihood/n quadratic term
    loglike <- -(1/VAR.MULT-1)
  }
  else # TWO KALMAN FILTERS WITH KNOWN VARIANCE AND FIXED ERROR
  {
    # eigen variances
    SIGMA <- area * exp(c(+1,-1)*ecc/2)

    # major axis likelihood
    CTMM$sigma <- SIGMA[1]
    KALMAN1 <- kalman(cbind(z[,1]),u,dt=dt,CTMM=CTMM,error=error)

    # minor axis likelihood
    CTMM$sigma <- SIGMA[2]
    KALMAN2 <- kalman(cbind(z[,2]),u,dt=dt,CTMM=CTMM,error=error)

    mu <- cbind(KALMAN1$mu,KALMAN2$mu)
    COV.mu <- array(c(KALMAN1$iW,diag(0,M),diag(0,M),KALMAN2$iW),c(M,M,2,2)) # -1/Hessian
    DOF.mu <- array(c(SIGMA[1]*KALMAN1$W,diag(0,M),diag(0,M),SIGMA[2]*KALMAN2$W),c(M,M,2,2))

    R.sigma <- c(KALMAN1$sigma + KALMAN2$sigma)/2 # standardized residual variances
    logdetCOV <- KALMAN1$logdet + KALMAN2$logdet
    logdetcov <- -log(det(KALMAN1$W)*det(KALMAN2$W))

    # this is loglikelihood/n quadratic term
    loglike <- -(R.sigma-1)
  }

  # restructure indices from m,n,x,y to x,m,n,y
  COV.mu <- aperm( COV.mu , c(3,1,2,4))
  DOF.mu <- aperm( DOF.mu , c(3,1,2,4))

  # transform results back
  if(ROTATE)
  {
    R <- rotate(+theta)

    mu <- mu %*% t(R)

    COV.mu <- array(COV.mu,c(2,M^2*2))
    DOF.mu <- array(DOF.mu,c(2,M^2*2))

    COV.mu <- R %*% COV.mu
    DOF.mu <- R %*% DOF.mu

    COV.mu <- array(COV.mu,c(2*M^2,2))
    DOF.mu <- array(DOF.mu,c(2*M^2,2))

    COV.mu <- COV.mu %*% t(R)
    DOF.mu <- DOF.mu %*% t(R)

    COV.mu <- array(COV.mu,c(2,M,M,2))
    DOF.mu <- array(DOF.mu,c(2,M,M,2))
  }

  # should I drop the indices in COV.mu and DOF.mu if possible ?
  COV.mu <- drop(COV.mu)
  DOF.mu <- drop(DOF.mu)

  # isotropic reduction if possible
  if(length(dim(DOF.mu))==2 && AXES==2)
  { if(DOF.mu[1,1]==DOF.mu[2,2] && DOF.mu[1,2]==0) { DOF.mu <- mean(diag(DOF.mu)) } }

  # re-write all of this to calculate constant, divide constant by n, and then subtract off from sum term by term?

  # likelihood constant/n: 2pi from det second term from variance-profiled quadratic term (which we will subtract if variance is not profiled)
  if(REML)
  { LL.CONST <- -(n-M)/n *AXES/2*log(2*pi) - AXES/2 } # not fixing this second term for REML yet... not wrong, but suboptimal maybe }
  else # ML constant
  { LL.CONST <- -AXES/2*log(2*pi*exp(1)) }

  # finish off the loglikelihood calculation
  loglike <- loglike - logdetCOV/2
  loglike <- n * (loglike + (LL.CONST-zero/n)) # I epxect the last part (all constants) to mostly cancel out
  logdetCOV <- n * logdetCOV

  # mean structure terms
  if(REML) { loglike <- loglike + CTMM$REML.loglike + logdetcov }
  CTMM$REML.loglike <- NULL

  if(verbose)
  {
    # assign variables
    CTMM$sigma <- sigma
    CTMM <- ctmm.repair(CTMM,K=K)

    # if(range)
    {
      mu <- drift@shift(mu,mu.center) # translate back to origin from center
      CTMM$mu <- mu
      CTMM$COV.mu <- COV.mu
      CTMM$DOF.mu <- DOF.mu
    }

    CTMM$loglike <- loglike
    attr(CTMM,"info") <- attr(data,"info")

    CTMM <- ctmm.ctmm(CTMM)
    return(CTMM)
  }
  else  { return(loglike) }
}


###########################################################
# FIT MODEL WITH LIKELIHOOD FUNCTION (convenience wrapper to optim)
ctmm.fit <- function(data,CTMM=ctmm(),method="ML",COV=TRUE,control=list(),trace=FALSE)
{
  # used for minimum scale of parameter inspection
  n <- length(data$t)
  dt <- stats::median(diff(data$t))
  df <- 2*pi/(last(data$t)-data$t[1])

  method <- match.arg(method,c("ML","pREML","pHREML","HREML","REML"))

  default <- list(method="Nelder-Mead",precision=1/2,maxit=.Machine$integer.max)
  control <- replace(default,names(control),control)
  precision <- control$precision
  optimizer <- control$method
  control$method <- NULL

  if(method=="REML") { REML <- TRUE }
  else { REML <- FALSE }

  # fit from MLE and not perturbed MLE
  if(!is.null(CTMM$MLE)) { CTMM <- CTMM$MLE }

  # clean/validate
  CTMM <- ctmm.ctmm(CTMM)
  drift <- get(CTMM$mean)

  CTMM$mu <- NULL # can always profile mu analytically
  UERE <- get.error(data,CTMM,flag=TRUE) # do we fit the error?
  axes <- CTMM$axes
  range <- CTMM$range
  if(is.null(CTMM$sigma)) { CTMM$sigma <- covm(stats::cov(get.telemetry(data,axes))) }

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

  # evaluate mean function for this data set if no vector is provided
  if(is.null(CTMM$mean.vec))
  {
    U <- drift(data$t,CTMM)
    CTMM$mean.vec <- U

    UU <- t(U) %*% U
    CTMM$REML.loglike <- (length(axes)/2)*(log(det(UU)) + ncol(U)*log(2*pi))
  }

  # id and characterize parameters for profiling
  pars <- NAMES <- parscale <- lower <- upper <- period <- NULL
  setup.parameters <- function(CTMM,profile=TRUE,linear=FALSE)
  {
    STUFF <- id.parameters(CTMM,profile=profile,linear=linear,UERE=UERE,dt=dt,df=df)
    NAMES <<- STUFF$NAMES
    parscale <<- STUFF$parscale
    lower <<- STUFF$lower
    upper <<- STUFF$upper
    period <<- STUFF$period
    # initial guess for optimization
    pars <<- get.parameters(CTMM,NAMES)
  }
  setup.parameters(CTMM)

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

  # CAVEATS (NEED TO UPDATE CODE)
  # if(error && circle && !isotropic) { warning("error==TRUE & circle==TRUE & isotropic==FALSE distorts error circle for speed.") }
  # if(circle && !range) { stop("circle==TRUE & range==FALSE are incompatible.") }
  # if(length(axes)==1 & !isotropic) { stop("Only multidimensional models can be anisotropic.")  }
  # if(length(axes)==1 & circle) { stop("Only multidimensional models can circle.") }
  # if(length(axes)>2) { stop("Can only handle 1-2 dimensions at a time.") }

  # NOW OPTIMIZE
  profile <- TRUE
  if(length(NAMES)==0) # EXACT
  {
    if(method %in% c("pHREML","HREML")) { REML <- TRUE } # IID limit pHREML -> REML

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

    # fast calculation of sigma covariance
    COVSTUFF <- COV.covm(CTMM$sigma,n=n,k=k.mean)
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
      DIFF <- genD(par=pars,fn=fn,zero=-CTMM$loglike,lower=lower,upper=upper,covariance=covariance(),parscale=parscale,Richardson=2,mc.cores=1)
      hess <- DIFF$hessian
      grad <- DIFF$gradient

      # more robust covariance calculation than straight inverse
      CTMM$COV <- cov.loglike(hess,grad)
      dimnames(CTMM$COV) <- list(NAMES,NAMES)
    }

    # pREML correction ############################### only do pREML if sufficiently away from boundaries
    TEST <- method %in% c("pREML","pHREML")
    if(TEST) { TEST <- all( eigen(hess,only.values=TRUE)$values > .Machine$double.eps*length(NAMES) ) }
    if(TEST)
    {
      # store MLE for faster model selection (ML is what is optimized, not pREML or REML)
      CTMM$MLE <- CTMM
      CTMM$MLE$method <- "ML"

      # parameter correction
      REML <- TRUE
      #ML.grad <- grad # save old ML gradient
      if(trace) { message("Calculating REML gradient.") }
      DIFF <- genD(par=pars,fn=fn,zero=-CTMM$loglike,lower=lower,upper=upper,covariance=CTMM$COV,parscale=parscale,Richardson=2,order=1,mc.cores=1,jacobian=FALSE)

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
      DIFF <- genD(par=pars,fn=fn,zero=-CTMM$loglike,lower=lower,upper=upper,covariance=covariance(),parscale=parscale,Richardson=2,mc.cores=1)
      # Using MLE gradient, which should be zero off boundary
      CTMM$COV <- cov.loglike(DIFF$hessian,grad)
    }
    else if(TEST) ### don't confuse the ML COV with pREML COV
    { CTMM$COV <- NULL }

    if(COV) { dimnames(CTMM$COV) <- list(NAMES,NAMES) }
    if(COV && !is.null(CTMM$MLE)) { dimnames(CTMM$MLE$COV) <- list(NAMES,NAMES) }
  }

  # model likelihood
  if(method!='ML') { CTMM$loglike <- ctmm.loglike(data,CTMM=CTMM,REML=FALSE,profile=FALSE) }
  CTMM$method <- method

  # covariance parameters only
  setup.parameters(CTMM,profile=FALSE,linear=FALSE)
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

  return(CTMM)
}


####### calculate variance and variance-covaraince from area and eccentricity information
area2var <- function(CTMM,MEAN=TRUE)
{
  VAR <- mean(diag(CTMM$sigma))
  COV <- CTMM$COV

  if(!CTMM$isotropic)
  {
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
    colnames(grad) <- rownames(COV)
    rownames(grad) <- c("variance",rownames(COV)[-(1:3)])

    COV <- grad %*% COV %*% t(grad)
  }

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

  # estimate circulation period
  if(CTMM$circle==1)
  {
    n <- length(data$t)

    # residuals
    z <- get.telemetry(data,CTMM$axes)
    z <- z - (u %*% CTMM$mu)

    # velocities !!! update this to minimally filtered estimate
    v <- cbind(diff(z[,1]),diff(z[,2])) / diff(data$t)
    # midpoint locations during velocity v
    z <- cbind(z[-1,1]+z[-n,1],z[-1,2]+z[-n,2])/2

    # average angular momentum
    L <- c(z[,1]%*%v[,2] - z[,2]%*%v[,1]) / (n-1)

    circle <- L / mean(diag(CTMM$sigma))
    # circle <- 2*pi/circle

    CTMM$circle <- circle
  }

  variogram.fit(variogram,CTMM=CTMM,name=name,interactive=interactive)
}
