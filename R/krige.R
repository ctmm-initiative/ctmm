# occurrence <- function(data,CTMM,H=0,res.time=10,res.space=10,grid=NULL,cor.min=0.5,dt.max=NULL) UseMethod("overlap") #S3 generic

################################
# Return hidden state estimates or simulations
################################
smoother <- function(DATA,CTMM,...)
{
  CTMM <- ctmm.prepare(DATA,CTMM)
  AXES <- length(CTMM$axes)

  t <- DATA$t

  dt <- c(Inf,diff(t))
  n <- length(t)

  u <- array(1,n)

  isotropic <- CTMM$isotropic
  sigma <- CTMM$sigma
  area <- sigma@par[1]
  if(AXES > 1)
  {
    ecc <- sigma@par[2]
    theta <- sigma@par[3]
  }
  else
  {
    ecc <- 0
    theta <- 0
  }

  K <- length(CTMM$tau)

  circle <- CTMM$circle
  if(circle) { circle <- 2*pi/circle }
  if(abs(circle) == Inf) { circle <- FALSE }

  R <- get.telemetry(DATA,CTMM$axes)
  # stationary mean function
  u <- array(1,n)

  V <- array(0,dim=c(n,AXES))
  COV <- array(0,dim=c(n,AXES,AXES))
  VCOV <- array(0,dim=c(n,AXES,AXES)) # velocity covariance

  # do we need to orient the data along the major an minor axes of sigma
  ROTATE <- !isotropic && (CTMM$error || circle)
  if(ROTATE)
  {
    M <- rotate(-theta)
    R <- R %*% t(M)
    mu <- M %*% as.numeric(CTMM$mu)
  }

  if(circle)
  {
    # proportional standardization from ellipse to circle
    if(ecc)
    {
      R[,1] <- R[,1] * exp(-ecc/4)
      R[,2] <- R[,2] * exp(+ecc/4)
    }
    R <- cbind(R[,1] + 1i*R[,2])

    # corotating frame
    M <- exp(-1i*circle*(t-t[1]))
    R <- M * R
    u <- M * u

    CTMM$sigma <- area
    KALMAN <- kalman(R,u,dt=dt,CTMM=CTMM,error=DATA$error,...)

    R <- KALMAN$Z[,"position",1]
    COV[,1,1] <- KALMAN$S[,"position","position"]
    COV[,2,2] <- KALMAN$S[,"position","position"]
    if(K>1)
    {
      V <- KALMAN$Z[,"velocity",1]
      VCOV[,1,1] <- KALMAN$S[,"velocity","velocity"]
      VCOV[,2,2] <- KALMAN$S[,"velocity","velocity"]
    }

    # spin back + velocity of frame
    M <- Conj(M)
    R <- M * R
    mu <- mu[1] + 1i*mu[2]
    R <- cbind(Re(R),Im(R))

    if(K>1)
    {
      V <- (M * V) + 1i*circle*(R-mu)
      V <- cbind(Re(V),Im(V))
    }

    # unstandardize
    if(ecc)
    {
      R[,1] <- R[,1] * exp(+ecc/4)
      R[,2] <- R[,2] * exp(-ecc/4)

      COV[,1,1] <- COV[,1,1] * exp(+ecc/2)
      COV[,2,2] <- COV[,2,2] * exp(-ecc/2)

      if(K>1)
      {
        V[,1] <- V[,1] * exp(+ecc/4)
        V[,2] <- V[,2] * exp(-ecc/4)

        VCOV[,1,1] <- VCOV[,1,1] * exp(+ecc/2)
        VCOV[,2,2] <- VCOV[,2,2] * exp(-ecc/2)
      }
    }
  }
  else if(CTMM$isotropic || !CTMM$error)
  {
    CTMM$sigma <- area
    KALMAN <- kalman(R,u,dt=dt,CTMM=CTMM,error=DATA$error,...)

    R <- KALMAN$Z[,"position",]
    R <- cbind(R) # R drops dimension-1
    COV <- KALMAN$S[,"position","position"] %o% covm(c(1,ecc,theta),isotropic,axes=CTMM$axes)
    if(K>1)
    {
      V <- KALMAN$Z[,"velocity",]
      V <- cbind(V) # R drops dimension-1
      VCOV <- KALMAN$S[,"velocity","velocity"] %o% covm(c(1,ecc,theta),isotropic,axes=CTMM$axes)
    }
  }
  else
  {
    #diagonalize data and then run two 1D Kalman filters with separate means
    # major axis likelihood
    CTMM$sigma <- area * exp(+ecc/2)
    KALMAN <- kalman(cbind(R[,1]),u,dt=dt,CTMM=CTMM,error=DATA$error,...)

    R[,1] <- KALMAN$Z[,"position",1]
    COV[,1,1] <- KALMAN$S[,"position","position"]
    if(K>1)
    {
      V[,1] <- KALMAN$Z[,"velocity",1]
      VCOV[,1,1] <- KALMAN$S[,"velocity","velocity"]
    }

    # minor axis likelihood
    CTMM$sigma <- area * exp(-ecc/2)
    KALMAN <- kalman(cbind(R[,2]),u,dt=dt,CTMM=CTMM,error=DATA$error,...)

    R[,2] <- KALMAN$Z[,"position",1]
    COV[,2,2] <- KALMAN$S[,"position","position"]
    if(K>1)
    {
      V[,2] <- KALMAN$Z[,"velocity",1]
      VCOV[,2,2] <- KALMAN$S[,"velocity","velocity"]
    }
  }

  if(ROTATE)
  {
    # transform results back
    M <- rotate(+theta)

    COV <- aperm(COV,perm=c(2,1,3))
    COV <- M %.% COV %.% t(M)
    COV <- aperm(COV,perm=c(2,1,3))

    R <- t(M %*% t(R))

    if(K>1)
    {
      V <- t(M %*% t(V))

      VCOV <- aperm(VCOV,perm=c(2,1,3))
      VCOV <- M %.% VCOV %.% t(M)
      VCOV <- aperm(VCOV,perm=c(2,1,3))
    }
  }

  colnames(R) <- CTMM$axes
  RETURN <- list(t=t,R=R,COV=COV)
  if(K>1) { RETURN$V <- V ; RETURN$VCOV <- VCOV }

  return(RETURN)
}


########################################
# fill in data gaps with missing observations of infinite error
########################################
fill.data <- function(data,CTMM=ctmm(tau=Inf),verbose=FALSE,t=NULL,dt=NULL,res=1,cor.min=0,dt.max=NULL)
{
  # prepare data error
  data$error <- get.error(data,CTMM)

  # FIX THE TIME GRID TO AVOID TINY DT
  if(is.null(t))
  {
    t <- data$t
    DT <- diff(t)

    # target resolution
    if(is.null(dt)){ dt <- stats::median(DT)/res }

    # maximum gap to bridge
    if(is.null(dt.max)) { dt.max <- -log(cor.min)*CTMM$tau[1] }

    # this regularization is not perfectly regular, but holds up to sampling drift in caribou data
    t.grid <- c()  # full (even-ish) grid
    dt.grid <- c() # local (numeric) sampling resolution
    t.new <- c()   # new times in this even grid
    for(i in 1:length(DT))
    {
      n.sub <- round(DT[i]/dt)+1
      n.sub <- max(n.sub,2) # fix for crazy small time-steps
      t.sub <- seq(from=t[i],to=t[i+1],length.out=n.sub)
      dt.sub <- DT[i]/(n.sub-1)

      # skip low correlation times in gaps
      INCLUDE <- (t.sub-t[i])<=dt.max/2 | (t[i+1]-t.sub)<=dt.max/2
      t.sub <- t.sub[INCLUDE]
      n.sub <- length(t.sub)

      dt.sub <- rep(dt.sub,n.sub)
      t.grid <- c(t.grid,t.sub)
      dt.grid <- c(dt.grid,dt.sub)

      t.new <- c(t.new,t.sub[c(-1,-n.sub)])
    }

    # half weight repeated endpoints in grid
    w.grid <- dt.grid
    REPEAT <- which(diff(t.grid)==0)
    w.grid[REPEAT] <- w.grid[REPEAT]/2
    w.grid[REPEAT+1] <- w.grid[REPEAT+1]/2
  }
  else # use a pre-specified time grid
  {
    t.new <- t[!(t %in% data$t)]
  }

  # empty observation row for these times
  blank <- data[1,]
  blank$error <- Inf
  blank <- blank[rep(1,length(t.new)),]
  blank$t <- t.new

  # attach empty measurements to data
  data <- rbind(data,blank)

  # sort times
  data <- data[sort.list(data$t,na.last=NA,method="quick"),]
  # this is now our fake data set to feed into the kalman smoother

  if(verbose)
  { data <- list(data=data,t.grid=t.grid,dt.grid=dt.grid,w.grid=w.grid) }

  return(data)
}


#################################
# Kriged Kernel Density Estimate
# H is your additional smoothing bandwidth matrix (zero by default)
# resolution is the number of kriged locations per median step
# cor.min is roughly the correlation required between locations to bridge them
# dt.max is (alternatively) the maximum gap allowed between locations to bridge them
#################################
occurrence <- function(data,CTMM,H=0,res.time=10,res.space=10,grid=NULL,cor.min=0.5,dt.max=NULL)
{
  dt <- stats::median(diff(data$t))

  if(length(H)==1) { H <- diag(H,2) }

  info <- attr(data,"info")
  SIGMA <- CTMM$sigma # diffusion matrix for later
  CTMM <- ctmm.prepare(data,CTMM)
  error <- get.error(data,CTMM)
  MIN.ERR <- min(error)

  # format data to be relatively evenly spaced with missing observations
  data <- fill.data(data,CTMM,verbose=TRUE,res=res.time,cor.min=cor.min,dt.max=dt.max)
  t.grid <- data$t.grid
  dt.grid <- data$dt.grid
  w.grid <- data$w.grid
  data <- data$data
  # run through smoother to get
  state <- smoother(data,CTMM,smooth=TRUE)

  # evenly sampled subset: data points (bridge ends) may be counted twice and weighted half
  GRID <- c( which(data$t %in% t.grid) , which(data$t %in% t.grid[which(diff(t.grid)==0)]) )
  GRID <- sort.int(GRID,method="quick")
  # GRID <- data$t %in% t.grid

  # t <- state$t[GRID]
  R <- state$R[GRID,]
  COV <- state$COV[GRID,,]
  n <- length(R[,1])

  # continuous velocities will give us more information to use
  if(length(CTMM$tau)>1)
  { V <- state$V[GRID,] }
  else # null velocity data otherwise
  { V <- array(0,c(n,2)) }

  # fake data
  data <- data.frame(x=R[,1],y=R[,2])

  # some covariances are slightly negative due to roundoff error
  # so here I clamp the negative singular values to zero
  COV <- lapply(1:n,function(i){ PDclamp(COV[i,,]) })
  COV <- unlist(COV)
  dim(COV) <- c(2,2,n)
  COV <- aperm(COV,perm=c(3,1,2))

  # uncertainties/bandwidths for this data
  h <- H # extra smoothing
  H <- array(0,c(n,2,2))
  for(i in 1:n)
  {
    # total covariance is bandwidth + krige uncertainty + kinetic numerical correction
    H[i,,] <- h + COV[i,,] + dt.grid[i]^2/12*(V[i,] %o% V[i,])
    # maybe publish the last part as a technical note
  }
  # there is probably a faster way to do that

  # estimate size of data blob
  dr <- diag(SIGMA)
  if(CTMM$range)
  {
    if(length(CTMM$tau)==1) #OU
    { dr <- dt/4 * dr/CTMM$tau[1] }
    else if(length(CTMM$tau)==2) #OUF
    { dr <- dt^2/24 * dr/prod(CTMM$tau)}
  }
  else # these sigmas are var/tau[1] diffusion limits
  {
    if(length(CTMM$tau)==1) #BM
    { dr <- dt/4 * dr }
    else if(length(CTMM$tau)==2) #IOU
    { dr <- dt^2/24 * dr/CTMM$tau[2] }
  }
  if(CTMM$error){ dr <- dr + MIN.ERR }
  dr <- sqrt(dr)

  # using the same data format as AKDE, but with only the ML estimate (alpha=1)
  KDE <- kde(data,H=H,W=w.grid,dr=dr/res.space)
  KDE$H <- diag(0,2)
  KDE <- new.UD(KDE,info=info)
  return(KDE)
}


##############################################
# SIMULATE DATA over time array t
simulate.ctmm <- function(object,nsim=1,seed=NULL,data=NULL,t=NULL,dt=NULL,res=1,...)
{
  if(!is.null(seed)){ set.seed(seed) }

  info <- attr(object,"info")
  axes <- object$axes
  AXES <- length(axes)

  CLASS <- class(data)
  CONDITIONAL <- FALSE
  if((CLASS=="telemetry" || CLASS=='data.frame'))
  {
    if(all(object$axes %in% names(data))) # condition off of data
    { CONDITIONAL <- TRUE }
    else # use data time & error for unconditional simulation
    { t <- data$t }
  }

  if(CONDITIONAL)
  {
    object <- ctmm.prepare(data,object)
    data <- fill.data(data,CTMM=object,t=t,dt=dt,res=res)
    object$error <- TRUE # avoids unit variance algorithm, data already contains fixed errors from fill.data
    data <- smoother(data,object,sample=TRUE)
    data <- cbind(t=data$t,data$R)
    data <- data.frame(data)

    # # the user probably only wants times t if t is specified
    # other stuff seems to be coded for uniform sampling...
    # if(!is.null(t))
    # {
    #   DEBUG <<- list(data=data,t=t)
    #   WHICH <- logical(length(data$t))
    #   i <- 1 ; j <- 1
    #   while(j < length(data$t))
    #   {
    #     if(t[i]==data$t[j])
    #     {
    #       WHICH[j] <- TRUE
    #       i <- i + 1
    #     }
    #     j <- j + 1
    #   }
    #   data <- data[WHICH,]
    # }

  }
  else # Gaussian simulation not conditioned off of any data
  {
    if(is.null(data)) { error <- FALSE } else { error <- get.error(data,object) } # get error if provided
    object <- ctmm.prepare(list(t=t),object)

    tau <- object$tau
    if(is.null(tau)) { tau = 0 }
    K <- length(tau)

    mu <- object$mu
    if(is.null(mu)) { mu <- array(0,c(1,AXES)) }
    sigma <- object$sigma
    if(is.null(sigma)) { sigma <- diag(1,AXES) }

    Lambda <- sqrtm(sigma)

    K <- length(tau)

    n <- length(t)
    dt <- c(Inf,diff(t)) # time lags

    # where we will store the data
    z <- array(0,c(n,AXES))

    # initial hidden state, for standardized process
    H <- array(0,c(K,AXES))

    object$sigma <- 1
    for(i in 1:n)
    {
      # tabulate propagators if necessary
      if((i==1)||(dt[i]!=dt[i-1]))
      {
        Langevin <- langevin(dt=dt[i],CTMM=object)
        Green <- Langevin$Green
        Sigma <- Langevin$Sigma
      }

      # generate standardized process
      # R drops dimensions ... so this code has to be a little weird
      H <- lapply(1:AXES,function(X){ MASS::mvrnorm(n=1,mu=as.vector(Green %*% H[,X]),Sigma=Sigma) })
      H <- do.call(cbind,H)
      # pull out location from hidden state
      z[i,] <- H[1,]
    }

    # rotate process if necessary
    circle <- object$circle
    if(circle)
    {
      z <- z[,1] + 1i*z[,2]

      circle <- 2*pi/circle
      R <- exp(1i*circle*(t-t[1]))
      z <- R * z

      z <- cbind(Re(z),Im(z))
    }

    z <- z %*% Lambda

    # calculate mean function
    mu <- object$mean.vec %*% mu
    z <- z + mu
    colnames(z) <- axes

    # throw in error
    if(any(as.logical(error)))
    {
      for(i in 1:length(axes))
      { z[,i] <- z[,i] + stats::rnorm(n,sd=sqrt(error)) }
    }

    data <- data.frame(z)
    data$t <- t
  }

  data <- new.telemetry(data,info=info)
  return(data)
}
#methods::setMethod("simulate",signature(object="ctmm"), function(object,...) simulate.ctmm(object,...))

simulate.telemetry <- function(object,nsim=1,seed=NULL,CTMM=NULL,t=NULL,dt=NULL,res=1,...)
{ simulate.ctmm(CTMM,nsim=nsim,seed=seed,data=object,t=t,dt=dt,res=res,...) }

##########################
# predict locations at certaint times !!! make times unique
##########################
predict.ctmm <- function(object,data=NULL,t=NULL,dt=NULL,res=1,...)
{
  info <- attr(object,"info")
  axes <- object$axes

  # Gaussian simulation not conditioned off of any data
  if(is.null(data))
  {
    object <- ctmm.prepare(list(t=t),object)

    mu <- object$mu
    if(is.null(mu)) { mu <- array(0,c(1,length(axes))) }

    # calculate mean function
    r <- object$mean.vec %*% mu
    colnames(r) <- axes

    data <- data.frame(r)
    data$t <- t
  }
  else # condition off of the data
  {
    object <- ctmm.prepare(data,object)
    K <- length(object$tau)
    data <- fill.data(data,CTMM=object,t=t,dt=dt,res=res)
    object$error <- TRUE # avoids unit variance algorithm
    data <- smoother(data,object,smooth=TRUE)

    NAMES <- colnames(data$R)
    COV <- data$COV

    VNAMES <- paste0("v.",NAMES)
    VNAMES -> colnames(data$V)
    VCOV <- data$VCOV

    data <- cbind(t=data$t,data$R,data$V)
    data <- data.frame(data)

    # flatten covariance matrix and include in data.frame
    for(i in 1:length(NAMES))
    {
      for(j in i:length(NAMES))
      {
        NAME <- paste("cov.",NAMES[i],".",NAMES[j],sep="")
        data[,NAME] <- COV[,i,j]
      }
    }

    if(K>1)
    {
      # flatten covariance matrix and include in data.frame
      for(i in 1:length(VNAMES))
      {
        for(j in i:length(VNAMES))
        {
          NAME <- paste("cov.",VNAMES[i],".",VNAMES[j],sep="")
          data[,NAME] <- VCOV[,i,j]
        }
      }
    }

    if(!is.null(t))
    {
      # pair down predictions only to those initially requested
      IN <- data$t %in% t
      data <- data[IN,]
      # remove duplicate predictions that arise with duplicate data (which is ok with error)
      IN <- which(diff(data$t)==0)
      if(length(IN)) { data <- data[-IN,] }
    }
  }

  data <- new.telemetry(data,info=info)
  return(data)
}

predict.telemetry <- function(object,CTMM=NULL,t=NULL,dt=NULL,res=1,...)
{ predict.ctmm(CTMM,data=object,t=t,dt=dt,res=res,...) }
