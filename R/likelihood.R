###############################
# Propagator/Green's function and Two-time correlation from Langevin equation for Kalman filter and simulations
langevin <- function(dt,CTMM)
{
  tau <- CTMM$tau
  CPF <- CTMM$CPF
  sigma <- methods::getDataPart(CTMM$sigma)
  K <- continuity(CTMM)

  if(K==0)
  {
    Green <- array(0,c(1,1))
    Sigma <- array(1,c(1,1))
  }
  else if(K==1)
  {
    c0 <- exp(-dt/tau)
    Green <- array(c0,c(1,1))
    Sigma <- array(1-c0^2,c(1,1))
  }
  else if(K==2)
  {
    f <- 1/tau

    DAMP <- FALSE
    if(!CPF) # overdamped to critically damped
    {
      Omega2 <- prod(f)

      nu <- diff(f)/2
      # should we use the exponential or hyperbolic representation?
      if(dt<Inf & (nu*dt)>0.8813736) # exponential
      { DAMP <- TRUE }
      else # hyperbolic
      { f <- mean(f) }
    }
    else if(CPF) # underdamped
    {
      nu <- 2*pi*f[1]
      f <- f[2]
      Omega2 <- f^2 + nu^2
    }

    TT <- 2*mean(f)/Omega2

    if(dt==Inf) # make this numerically relative in future
    {
      Green <- rbind( c(0,0) , c(0,0) )
      Sigma <- rbind( c(1,0) , c(0,Omega2) )
    }
    else
    {
      if(DAMP) # very overdamped
      {
        Exp <- exp(-dt/tau)/diff(tau)
        c0 <- diff(Exp*tau)
        c1 <- -diff(Exp)
        c2 <- diff(Exp/tau)
      }
      else # mildly overdamped to underdamped
      {
        Exp <- exp(-f*dt)

        if(!CPF) # overdamped
        {
          SincE <- sinch(nu*dt)*Exp
          CosE <- cosh(nu*dt)*Exp
        }
        else # underdamped
        {
          SincE <- sinc(nu*dt)*Exp
          CosE <- cos(nu*dt)*Exp
        }

        c0 <- CosE + (f*dt)*SincE
        c1 <- -(Omega2*dt)*SincE
        c2 <- -Omega2*(CosE - (f*dt)*SincE)
      }

      Green <- rbind( c(c0,-c1/Omega2) , c(c1,-c2/Omega2) )
      Sigma <- -TT*c1^2  #off-diagonal term
      Sigma <- rbind( c(1,0)-c(c0^2+c1^2/Omega2,Sigma) , c(0,Omega2)-c(Sigma,c1^2+c2^2/Omega2) )
    }
  }

  return(list(Green=Green, Sigma=sigma*Sigma))
}

#############################################################
# Internal Kalman filter/smoother for multiple derivatives, dimensions, trends
# Kalman filter/smoother for matrix P operator and multiple mean functions
# this is a lower level function
# more for generalizability/readability than speed at this point
kalman <- function(z,u,dt,CTMM,error=rep(0,nrow(z)),smooth=FALSE,sample=FALSE,weight=FALSE,residual=FALSE)
{
  n <- nrow(z)
  DATA <- 1:ncol(z)

  # glob data and mean functions together for Kalman filter
  z <- cbind(z,u)
  VEC <- ncol(z)

  # indices of mean functions
  MEAN <- (last(DATA)+1):VEC

  tau <- CTMM$tau
  K <- length(tau)  # dimension of hidden state per spatial dimension

  # observed dimensions
  OBS <- 1

  Id <- diag(OBS)
  IdH <- diag(K)

  # observable state projection operator (will eventially upgrade to use velocity telemetry)
  P <- array(0,c(K,OBS))
  P[1:OBS,] <- 1

  # forecast estimates for zero-mean, unit-variance filters
  zFor <- array(0,c(n,K,VEC))
  sFor <- array(0,c(n,K,K))

  # forcast residuals
  zRes <- array(0,c(n,OBS,VEC))
  sRes <- array(0,c(n,OBS,OBS))

  # concurrent/filtered estimates
  zCon <- array(0,c(n,K,VEC))
  sCon <- array(0,c(n,K,K))

  # propagation information
  Green <- array(0,c(n,K,K))
  Sigma <- array(0,c(n,K,K))

  # Propagators from Langevin equation
  for(i in 1:n)
  {
    # does the time lag change values? Then update the propagators.
    if(i==1 || dt[i] != dt[i-1])
    { Langevin <- langevin(dt=dt[i],CTMM=CTMM) }

    Green[i,,] <- Langevin$Green
    Sigma[i,,] <- Langevin$Sigma
  }

  # first zForcast is properly zeroed
  sFor[1,,] <- Sigma[1,,]

  for(i in 1:n)
  {
    # residual covariance
    sForP <- sFor[i,,] %*% P # why do I need this?
    ERROR <- error[i]*Id
    sRes[i,,] <- ((t(P) %*% sForP) + ERROR)

    # forcast residuals
    zRes[i,,] <- z[i,] - (t(P) %*% zFor[i,,])

    if(all(abs(sRes[i,,])<Inf)){ Gain <- sForP %*% solve(sRes[i,,]) }
    else { Gain <- sForP %*% (0*Id) } # solve() doesn't like this case...

    # concurrent estimates
    zCon[i,,] <- zFor[i,,] + (Gain %*% zRes[i,,])
    # manifestly symmetric form
    # sCon[i,,] <- (sFor[i,,] - (Gain %*% t(sForP)))
    # manifestly positive-definite form (Joseph)
    JOSEPH <- IdH - (Gain %*% t(P))
    sCon[i,,] <- (JOSEPH %*% sFor[i,,] %*% t(JOSEPH))
    if(error[i]<Inf & error[i]>0) { sCon[i,,] <- sCon[i,,] + (Gain %*% ERROR %*% t(Gain)) } # otherwise Gain==0
    # this is supposed to be more stable numerically

    # update forcast estimates for next iteration
    if(i<n)
    {
      #update forcast estimates now
      zFor[i+1,,] <- Green[i+1,,] %*% zCon[i,,]
      sFor[i+1,,] <- ((Green[i+1,,] %*% sCon[i,,] %*% t(Green[i+1,,])) + Sigma[i+1,,])
    }
  }

  # return residuals of Kalman filter only
  if(residual) { return(zRes[,1,]/sqrt(sRes[,1,1])) }

  if(!weight)
  {
    # general quadratic form, tracing over times
    # hard-code weights for location observation only
    M <- Adj(zRes[,1,]) %*% (zRes[,1,]/sRes[,1,1])

    # estimate mean parameter
    W <- as.matrix(M[MEAN,MEAN])
    mu <- solve(W) %*% M[MEAN,DATA]
  }

  # returned profiled mean
  if(!smooth && !sample && !weight)
  {
    # hard coded for position observations
    sigma <- (M[DATA,DATA] - (Adj(mu) %*% W %*% mu))/n

    # log det autocorrelation matrix == trace log autocorrelation matrix
    logdet <- mean(log(sRes)) # this is 1/n times the full term

    return(list(mu=mu,W=W,sigma=sigma,logdet=logdet))
  }

  #########################
  # INVERSE FILTER (DOESN'T SEEM TO WORK)
  #########################
  if(weight)
  {
    # would be standardized residuals in diagonalization
    zRes[,1,] <- zRes[,1,]/sRes[,1,1]

    # forecast estimates for zero-mean, unit-variance filters
    zFor <- array(0,c(n,K,VEC))

    # run Kalman filter loop in reverse but residual -> data instead of data -> residual
    for(i in 1:n)
    {
      # forcast residuals -> fake data
      zRes[i,,] + (t(P) %*% zFor[i,,]) -> z[i,]

      sForP <- sFor[i,,] %*% P # why do I need this?
      if(all(abs(sRes[i,,])<Inf)){ Gain <- sForP %*% solve(sRes[i,,]) }
      else { Gain <- sForP %*% (0*Id) } # solve() doesn't like this case...

      # concurrent estimates
      zCon[i,,] <- zFor[i,,] + (Gain %*% zRes[i,,])

      # update forcast estimates for next iteration
      if(i<n) { zFor[i+1,,] <- Green[i+1,,] %*% zCon[i,,] }
    }

    return(z)
  }

  # delete residuals
  rm(zRes,sRes)

  #####################
  # KALMAN SMOOTHER
  #####################
  # Finish detrending the effect of a stationary mean
  MU <- zFor[,,MEAN]
  dim(MU) <- c(n*K,length(MEAN))
  MU <- MU %*% mu
  dim(MU) <- c(n,K,length(DATA))
  zFor[,,DATA] <- zFor[,,DATA,drop=FALSE] - MU
  # there has to be a better way to do this?
  MU <- zCon[,,MEAN]
  dim(MU) <- c(n*K,length(MEAN))
  MU <- MU %*% mu
  dim(MU) <- c(n,K,length(DATA))
  zCon[,,DATA] <- zCon[,,DATA,drop=FALSE] - MU
  # why does R drop dimensions so randomly?

  # delete u(t)
  zCon <- zCon[,,DATA,drop=FALSE]
  zFor <- zFor[,,DATA,drop=FALSE]
  # drop=FALSE must be here for BM/OU and I don't fully understand why

  #################
  # RANDOM SAMPLER: end point
  #################
  if(sample)
  {
    zCon[n,,] <- sapply(DATA,function(d){MASS::mvrnorm(mu=zCon[n,,d],Sigma=sCon[n,,])})
    sCon[n,,] <- array(0,c(K,K))
  }

  # upgrade concurrent estimates to Kriged estimates
  for(i in (n-1):1)
  {
    # DOUBLE CHECK THIS
    L <- sCon[i,,] %*% t(Green[i+1,,]) %*% solve(sFor[i+1,,])

    # overwrite concurrent estimate with smoothed estimate
    # RTS smoother
    zCon[i,,] <- zCon[i,,] + L %*% (zCon[i+1,,]-zFor[i+1,,])
    # manifestly symmetric backsolver
    # sCon[i,,] <- sCon[i,,] + L %*% (sCon[i+1,,]-sFor[i+1,,]) %*% t(L)
    # manifestly PSD backsolver
    JOSEPH <- IdH - (L %*% Green[i+1,,])
    sCon[i,,] <- (JOSEPH %*% sCon[i,,] %*% t(JOSEPH)) + (L %*% (sCon[i+1,,]+Sigma[i+1,,]) %*% t(L))

    #################
    # RANDOM SAMPLER
    #################
    if(sample)
    {
      zCon[i,,] <- sapply(DATA,function(d){MASS::mvrnorm(mu=zCon[i,,d],Sigma=sCon[i,,])})
      sCon[i,,] <- array(0,c(K,K))
    }
  }

  # restore stationary mean to locations only
  zCon[,1,] <- zCon[,1,] + (cbind(u) %*% mu)

  zname <- c("position")
  if(K>1) { zname <- c(zname,"velocity") }
  dimnames(zCon) <- list(NULL,zname,c("x","y")[DATA])
  dimnames(sCon) <- list(NULL,zname,zname)

  # return smoothed states
  # this object is temporary
  state <- list(CTMM=CTMM,Z=zCon,S=sCon)
  class(state) <- "state"

  return(state)
}

####################################
# log likelihood function
####################################
ctmm.loglike <- function(data,CTMM=ctmm(),REML=FALSE,zero=0,verbose=FALSE)
{
  n <- length(data$t)
  AXES <- length(CTMM$axes)

  # prepare model for numerics
  CTMM <- ctmm.prepare(data,CTMM,REML=REML)

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
    area <- NA
    ecc <- 0
    theta <- 0
  }

  circle <- CTMM$circle
  if(circle) { circle <- 2*pi/circle }
  if(abs(circle) == Inf) { circle <- FALSE }

  n <- length(data$t)

  t <- data$t
  # time lags
  dt <- c(Inf,diff(t))

  # data z and mean vector u
  z <- get.telemetry(data,CTMM$axes)
  u <- CTMM$mean.vec
  M <- ncol(u) # number of linear parameters per spatial dimension
  if(REML)
  {
    # variance debias factor
    VAR.MULT <- n/(n-M)
    # likelihood constant: 2pi from det second term from variance-profiled quadratic term (which we will subtract if variance is not profiled)
    LL.CONST <- - (n-M)*AXES/2*log(2*pi) - n*AXES/2 # not fixing this second term for REML yet... not wrong, but suboptimal maybe
  }
  else # ML constant
  {
    VAR.MULT <- 1
    LL.CONST <- - n*AXES/2*log(2*pi*exp(1))
  }


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

    if(!CTMM$error || (isotropic && UERE<3))
    {
      CTMM$sigma <- 1
      KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)

      ML.area <- as.numeric(Re(KALMAN$sigma))/2
      # profile variance if unspecified
      if(is.na(area))
      {
        # debias if REML
        area <- VAR.MULT*ML.area

        sigma <- covm(c(area,ecc,theta),isotropic=isotropic,axes=CTMM$axes)

        # convert from error ratio to absolute error
        CTMM$error <- sqrt(area) * CTMM$error
      }

      COV.mu <- solve(KALMAN$W/area)
      DOF.mu <- KALMAN$W

      # this is loglikelihood/n
      loglike <- -KALMAN$logdet -log(area) - (ML.area/area-1)

      if(REML) { loglike <- loglike - log(det(KALMAN$W/area))/n }
    }
    else
    {
      CTMM$sigma <- area
      KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)

      # what does this term here represent? relative variance or something
      R.sigma <- KALMAN$sigma

      COV.mu <- solve(KALMAN$W)
      DOF.mu <- area * KALMAN$W

      # this is loglikelihood/n
      loglike <- -KALMAN$logdet - (1/2)*(R.sigma-2)

      if(REML) { loglike <- loglike - log(det(KALMAN$W))/n }
    }

    loglike <- Re(loglike)
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
  }
  else if(!CTMM$error) # ONE KALMAN FILTER WITH UNKNOWN COVARIANCE AND NO ERROR
  {
    CTMM$sigma <- 1
    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)

    ML.sigma <- KALMAN$sigma
    # profile covariance if sigma unspecified
    if(is.null(sigma))
    {
      # debias if REML
      sigma <- VAR.MULT*ML.sigma

      sigma <- covm(sigma,isotropic=isotropic,axes=CTMM$axes)
    }

    mu <- KALMAN$mu
    COV.mu <- solve(KALMAN$W) %o% sigma
    DOF.mu <- KALMAN$W %o% diag(AXES)

    # this is loglikelihood/n
    loglike <- -(AXES/2)*KALMAN$logdet -(1/2)*log(det(sigma)) - (1/2)*sum(diag(ML.sigma %*% solve(sigma))-1)

    # tensor product rule for determinants
    if(REML) { loglike <- loglike - (1/2)*log(det(KALMAN$W)^AXES/det(sigma)^M)/n }
  }
  else if(isotropic && UERE<3) # ONE KALMAN FILTER WITH UNKNOWN VARIANCE AND KNOWN ERROR RATIO
  {
    CTMM$sigma <- 1
    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)

    ML.area <- mean(diag(KALMAN$sigma))
    # profile variance if unspecified
    if(is.na(area))
    {
      # debias if REML
      area <- VAR.MULT*ML.area

      sigma <- covm(area,isotropic=isotropic,axes=CTMM$axes)

      # convert from error ratio to absolute error
      CTMM$error <- sqrt(area) * CTMM$error
    }

    mu <- KALMAN$mu
    COV.mu <- solve(KALMAN$W) %o% sigma
    DOF.mu <- KALMAN$W %o% diag(AXES)

    # this is loglikelihood/n
    loglike <- -(AXES/2)*KALMAN$logdet -(1/2)*log(det(sigma)) - (1/2)*AXES*(ML.area/area-1)

    # tensor product rule for determinants
    # not sure about this equation !!!
    if(REML) { loglike <- loglike - (1/2)*log(det(KALMAN$W)^AXES/det(sigma)^M)/n }
  }
  else if(isotropic) # ONE KALMAN FILTER WITH KNOWN VARIANCE AND KNOWN ERROR
  {
    CTMM$sigma <- area
    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)

    mu <- KALMAN$mu
    COV.mu <- solve(KALMAN$W) %o% diag(AXES)
    DOF.mu <- (area * KALMAN$W) %o% diag(AXES)

    # includes error and standardization
    R.sigma <- KALMAN$sigma
    R.sigma <- mean(diag(R.sigma))

    # this is loglikelihood/n
    loglike <- -(AXES/2)*KALMAN$logdet -(AXES/2)*(R.sigma-1)

    if(REML) { loglike <- loglike - (AXES/2)*log(det(KALMAN$W))/n }
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

    logdet <- KALMAN1$logdet + KALMAN2$logdet
    R.sigma <- KALMAN1$sigma + KALMAN2$sigma # standardized residual variances

    mu <- cbind(KALMAN1$mu,KALMAN2$mu)
    COV.mu <- array(c(solve(KALMAN1$W),diag(0,M),diag(0,M),solve(KALMAN2$W)),c(M,M,2,2)) # -1/Hessian
    DOF.mu <- array(c(SIGMA[1]*KALMAN1$W,diag(0,M),diag(0,M),SIGMA[2]*KALMAN2$W),c(M,M,2,2))

    # this is loglikelihood/n
    loglike <- -(1/2)*logdet - (1/2)*(R.sigma-2)

    if(REML) { loglike <- loglike - (1/2)*log(det(KALMAN1$W)*det(KALMAN2$W))/n }
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
  if(length(dim(DOF.mu))==2 && AXES==2) { if(DOF.mu[1,1]==DOF.mu[2,2] && DOF.mu[1,2]==0) { DOF.mu <- mean(diag(DOF.mu)) } }

  # mean structure terms
  if(REML) { loglike <- loglike + CTMM$REML.loglike }
  CTMM$REML.loglike <- NULL

  # finish off the loglikelihood calculation
  loglike <- n*(loglike-(zero-LL.CONST)/n)

  if(verbose)
  {
    # assign variables
    CTMM$sigma <- sigma
    CTMM <- ctmm.repair(CTMM)

    # if(range)
    {
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
ctmm.fit <- function(data,CTMM=ctmm(),method="ML",control=list(),trace=FALSE)
{
  method <- match.arg(method,c("ML","pREML","REML"))

  default <- list(method="Nelder-Mead",precision=1/2,maxit=.Machine$integer.max)
  control <- replace(default,names(control),control)
  precision <- control$precision
  optimizer <- control$method
  control$method <- NULL

  # clean/validate
  CTMM <- ctmm.ctmm(CTMM)

  if(method=="REML") { REML <- TRUE }
  else { REML <- FALSE }

  # basic info
  n <- length(data$t)
  tau <- CTMM$tau
  CPF <- CTMM$CPF
  circle <- CTMM$circle
  sigma <- CTMM$sigma
  if(!is.null(sigma)) { sigma <- sigma@par }
  CTMM$mu <- NULL # can always profile mu analytically
  isotropic <- CTMM$isotropic
  error <- CTMM$error
  UERE <- get.error(data,CTMM,flag=TRUE) # do we fit the error?
  range <- CTMM$range
  axes <- CTMM$axes

  # save for fitting
  COV.init <- CTMM$COV
  # erase previous fitting info if present
  CTMM$COV <- NULL
  CTMM$COV.mu <- NULL
  CTMM$DOF.mu <- NULL

  # evaluate mean function for this data set if no vector is provided
  if(is.null(CTMM$mean.vec))
  {
    drift <- get(CTMM$mean)
    U <- drift(data$t,CTMM)
    CTMM$mean.vec <- U

    UU <- t(U) %*% U
    CTMM$REML.loglike <- (length(axes)/2)*(log(det(UU)) + ncol(U)*log(2*pi))
  }

  # CAVEATS
  if(error && circle && !isotropic) { warning("error==TRUE & circle==TRUE & isotropic==FALSE distorts error circle for speed.") }
  if(circle && !range) { stop("circle==TRUE & range==FALSE are incompatible.") }
  if(length(axes)==1 & !isotropic) { stop("Only multidimensional models can be anisotropic.")  }
  if(length(axes)==1 & circle) { stop("Only multidimensional models can circle.") }
  if(length(axes)>2) { stop("Can only handle 1-2 dimensions at a time.") }

  # parameter indices for non-profiled parameters that we have to numerically optimize
  # PARAMETERS THAT WE MAY OR MAY NOT DIFFERENTIATE WRT
  if(error)
  {
    if(UERE<3 && isotropic) # can profile variance free when fitting error
    {
      SIGMA <- (if(isotropic){NULL}else{1:2})
      SIGMAV <- (if(isotropic){NULL}else{2:3})
    }
    else # must fit all 1-3 covariance parameters
    {
      SIGMA <- 1:(if(isotropic){1}else{3})
      SIGMAV <- 1:(if(isotropic){1}else{3})
    }
  }
  else if(circle) # must fit cov shape, but variance is free
  {
    SIGMA <- (if(isotropic){NULL}else{1:2})
    SIGMAV <- (if(isotropic){NULL}else{2:3})
  }
  else
  { SIGMA <- NULL }

  # PARAMETERS THAT WE WILL DIFFERENTIATE WRT
  NAMES <- NULL
  TAU <- NULL
  CIRCLE <- NULL
  ERROR <- NULL
  calPARS <- function()
  {
    if(length(SIGMA)==1) { NAMES <<- "area" }
    else if(length(SIGMA)==3) { NAMES <<- c("area","eccentricity","angle") }

    K <- length(tau)
    if(K+range==1) # IID, BM
    { TAU <<- NULL }
    else
    {
      NAMES <<- c(NAMES,paste("tau",names(tau)[(1+(1-range)):K]))
      TAU <<- length(SIGMA) + 1:(K-(1-range))
    }

    if(circle)
    {
      NAMES <<- c(NAMES,"circle")
      CIRCLE <<- length(SIGMA) + length(TAU) + 1
    }

    # If error is an error estimate rather than TRUE, and if there is no error annotated, then we will fit error
    if(error & UERE<3)
    {
      NAMES <<- c(NAMES,"error")
      ERROR <<- length(SIGMA) + length(TAU) + length(CIRCLE) + 1
    }
  }
  calPARS()

  # numerically fit parameters
  PARS <- length(SIGMA) + length(TAU) + length(CIRCLE) + length(ERROR)
  # degrees of freedom, including the mean, variance/covariance, tau, and error model
  k.mean <- ncol(CTMM$mean.vec)
  # covariance parameters only
  nu <- (if(isotropic){1}else{3}) + length(TAU) + length(CIRCLE) + length(ERROR)
  # all parameters
  k <- nu + length(axes)*(if(range){k.mean}else{k.mean-1})

  # OPTIMIZATION GUESS (pars)
  # also construct reasonable parscale
  pars <- NULL
  parscale <- NULL
  lower <- NULL
  upper <- NULL
  period <- NULL
  calpars <- function()
  {
    pars <<- NULL
    parscale <<- NULL
    period <<- NULL

    if(length(SIGMA))
    {
      # need some initial guess...
      if(is.null(sigma)) { sigma <<- drift@init(data,CTMM)$sigma@par }
      pars <<- sigma[SIGMAV]
      parscale <<- c(parscale,c(sigma[1],log(2),pi/2)[SIGMAV])
      LO <- c(0,0,-Inf)
      lower <<- LO[SIGMAV]
      UP <- c(Inf,Inf,Inf)
      upper <<- UP[SIGMAV]
      period <<- c(period,c(F,F,pi)[SIGMAV])

      # can we profile the variance, if so delete the guess
      if(!is.element(1,SIGMAV)) { sigma[1] <<- NA }
    }

    if(length(TAU))
    {
      PS <- pars.tauv(tau)
      pars <<- c(pars,PS) # pull out relevant tau elements
      parscale <<- c(parscale,PS)
      LO <- 0*PS
      lower <<- c(lower,LO)
      UP <- Inf*(PS+1)
      upper <<- c(upper,UP)
      period <<- c(period,F)
    }

    # use 1/T as circulation parameter
    if(length(CIRCLE))
    {
      pars <<- c(pars,1/circle)
      parscale <<- c(parscale,1/abs(circle))
      lower <<- c(lower,-Inf)
      upper <<- c(upper,Inf)
      period <<- c(period,F)
    }

    if(length(ERROR))
    {
      # fitting the error/variance ratio and not absolute error
      # in this case, one overall variance of motion + error is fit analytically
      if(isotropic)
      {
        pars <<- c(pars,error/sqrt(sigma[1]))
        parscale <<- c(parscale,error/sqrt(sigma[1]))
      }
      else # regular error fitting
      {
        pars <<- c(pars,error)
        parscale <<- c(parscale,error)
      }
      lower <<- c(lower,0)
      upper <<- c(upper,Inf)
      period <<- c(period,F)
    }
    # names(pars) <<- NAMES
  }
  calpars()

  cleanpars <- function(p)
  {

    if(length(SIGMA))
    {
      # enforce positivity
      if(isotropic)
      { p[SIGMA] <- abs(p[SIGMA]) }
      else
      {
        if(p[SIGMA][2]<0) { p[SIGMA][2] <- p[SIGMA][2] - pi/2 } # swap major and minor axis
        p[SIGMA][-3] <- abs(p[SIGMA][-3]) # positive area & eccentricity
        p[SIGMA][3] <- (((p[SIGMA][3]/pi+1/2) %% 1) - 1/2)*pi # wrap angle back
      }
    }

    if(length(TAU))
    {
      # enforce positivity
      p[TAU] <- abs(p[TAU])

      if(!CPF) { p[TAU] <- sort(p[TAU],decreasing=TRUE) }
    }

    if(length(ERROR))
    {
      # enforce positivity
      p[ERROR] <- abs(p[ERROR])
    }

    return(p)
  }

  # OPTIMIZATION FUNCTION (fn)
  # optional argument lengths: TAU, TAU+1, TAU+SIGMA
  fn <- function(p,zero=0)
  {
    p <- cleanpars(p)

    # Fix sigma if provided, up to degree provided
    if(length(SIGMA)==0)
    { sigma <- NULL }
    else
    {
      # write over inherited sigma with any specified parameters
      sigma[SIGMAV] <- p[SIGMA]

      # for some reason optim jumps to crazy big parameters
      if(!isotropic && !is.na(sigma[1]) && sigma[1]*exp(abs(sigma[2])/2)==Inf) { return(Inf) }
    }

    # fix tau from par
    if(length(TAU)==0)
    {
      if(range) { tau <- NULL }
      else { tau <- Inf }
    }
    else
    {
      tau <- p[TAU]

      if(!range) { tau <- c(Inf,tau) }
    }

    # fix circulation from par
    if(length(CIRCLE))
    { circle <- 1/p[CIRCLE] }

    # fix error from par
    if(length(ERROR))
    { error <- p[ERROR] }

    # overwrite modified parameters inside this function's environment only
    CTMM$tau <- tau
    CTMM$circle <- circle
    CTMM$sigma <- covm(sigma,isotropic=isotropic,axes=axes)
    CTMM$error <- error

    # negative log likelihood
    return(-ctmm.loglike(data,CTMM,REML=REML,zero=-zero))
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
  if(PARS==0) # EXACT
  {
    # Bi-variate Gaussian with zero error
    CTMM$sigma <- NULL
    CTMM <- ctmm.loglike(data,CTMM=CTMM,REML=REML,verbose=TRUE)

    # pREML perturbation
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
    pars <- cleanpars(RESULT$par)
    # copy over hessian from fit to COV.init ?

    # write best estimates over initial guess
    store.pars <- function(pars)
    {
      pars <- cleanpars(pars)

      tau <<- pars[TAU]
      if(!range){ tau <<- c(Inf,tau) }
      names(tau) <<- names(CTMM$tau)

      # save circulation if numerically optimized
      if(circle)
      {
        if(pars[CIRCLE])
        { circle <<- 1/pars[CIRCLE] }
        else
        { circle <<- FALSE }
        # In case ML circulation is zero, deactivate it in the model
      }

      # save sigma if numerically optimized
      if(length(SIGMA))
      { sigma[SIGMAV] <<- pars[SIGMA] }
      else
      { sigma <<- NULL }

      # save error magnitude if modeled
      if(length(ERROR)) { error <<- pars[ERROR] }

      CTMM$tau <<- tau
      CTMM$circle <<- circle
      CTMM$sigma <<- covm(sigma,isotropic=isotropic,axes=axes)
      CTMM$error <<- error

      # verbose ML information
      # this is a wasted evaluation !!! store verbose glob in environment?
      CTMM <<- ctmm.loglike(data,CTMM,REML=REML,verbose=TRUE)
      # zero-taus can be deleted above
      tau <<- CTMM$tau
    }
    store.pars(pars)

    # need to differentiate, no longer solve sigma analytically
    sigma <- CTMM$sigma@par
    SIGMA <- 1:(if(isotropic){1}else{3})
    SIGMAV <- 1:(if(isotropic){1}else{3})
    calPARS()
    # if error was just the relative error, then convert back
    if(length(ERROR) & isotropic) { error <- error * sqrt(sigma[1]) }

    calpars()
    if(trace) { message("Calculating covariance.") }
    DIFF <- genD(par=pars,fn=fn,zero=-CTMM$loglike,lower=lower,upper=upper,covariance=covariance(),parscale=parscale,Richardson=2,mc.cores=1)
    hess <- DIFF$hessian
    grad <- DIFF$gradient
    # robust covariance calculation
    CTMM$COV <- cov.loglike(hess,grad)
    dimnames(CTMM$COV) <- list(NAMES,NAMES)

    CTMM$method <- method
    # perturbative correction
    if(method=="pREML")
    {
      # store MLE for faster model selection (ML is what is optimized, not pREML or REML)
      CTMM$MLE <- CTMM
      CTMM$MLE$method <- "ML"
      loglike <- CTMM$loglike

      COV <- CTMM$COV
      pNAMES <- names(pars)

      # parameter correction
      REML <- TRUE
      #ML.grad <- grad # save old ML gradient
      if(trace) { message("Calculating pREML correction.") }
      DIFF <- genD(par=pars,fn=fn,zero=-CTMM$loglike,lower=lower,upper=upper,covariance=covariance(),parscale=parscale,Richardson=2,mc.cores=1)
      d.pars <- -c(COV %*% DIFF$gradient) # COV is -1/Hessian, grad is of -loglike

      # gradient transformations
      # Jacobian
      J <- diag(length(pars))
      dimnames(J) <- list(pNAMES,pNAMES)

      # final perturbation
      d.pars <- c(J %*% d.pars)

      # increment transformed parameters
      pars <- pars + d.pars

      # transform back
      names(pars) <- pNAMES

      # store parameter correction
      store.pars(pars)
      CTMM$loglike <- loglike #MLE

      # covariance correction
      # REML hessian & ML gradient (at boundaries)
      CTMM$COV <- cov.loglike(DIFF$hessian,grad)
      #CTMM$COV <- 2*COV - (COV %*% hess %*% COV)
    }

    # convert from circulation frequency to circulation period
    if(circle)
    {
      g <- -circle^2
      CTMM$COV[CIRCLE,] <- g*CTMM$COV[CIRCLE,]
      CTMM$COV[,CIRCLE] <- g*CTMM$COV[,CIRCLE]

      if(method=="pREML")
      {
        g <- -CTMM$MLE$circle^2
        CTMM$MLE$COV[CIRCLE,] <- g*CTMM$MLE$COV[CIRCLE,]
        CTMM$MLE$COV[,CIRCLE] <- g*CTMM$MLE$COV[,CIRCLE]
      }
      # before adding more junk here, just change the ctmm object to store frequency and not period for circulation
    }

    dimnames(CTMM$COV) <- list(NAMES,NAMES)
    if(method=="pREML") { dimnames(CTMM$MLE$COV) <- list(NAMES,NAMES) }
  }

  CTMM$AIC <- (2*k-2*CTMM$loglike)
  CTMM$AICc <- CTMM$AIC + 2*k*(k+nu)/(length(axes)*n-k-nu)

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

  # estimate circulation period
  if(CTMM$circle==1)
  {
    n <- length(data$t)

    # residuals
    z <- get.telemetry(data,CTMM$axes)
    z <- t(t(z)-mu)

    # velocities
    v <- cbind(diff(z[,1]),diff(z[,2])) / diff(data$t)
    # midpoint locations during velocity v
    z <- cbind(z[-1,1]+z[-n,1],z[-1,2]+z[-n,2])/2

    # average angular momentum
    L <- (z[,1]%*%v[,2] - z[,2]%*%v[,1]) / (n-1)

    circle <- L / mean(diag(CTMM$sigma))
    circle <- 2*pi/circle

    CTMM$circle <- circle
  }

  variogram.fit(variogram,CTMM=CTMM,name=name,interactive=interactive)
}
