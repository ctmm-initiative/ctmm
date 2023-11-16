# occurrence <- function(data,CTMM,H=0,res.time=10,res.space=10,grid=NULL,cor.min=0.5,dt.max=NULL) UseMethod("overlap") #S3 generic

################################
# Return hidden state estimates or simulations
################################
smoother <- function(data,CTMM,precompute=FALSE,sample=FALSE,residual=FALSE,...)
{
  if(is.null(CTMM$error.mat)) { CTMM <- ctmm.prepare(data,CTMM) }
  if(is.null(data$record)) { data$record <- TRUE } # real recorded data or blank/empty timestamps from fill-data
  AXES <- length(CTMM$axes)

  t <- data$t

  dt <- c(Inf,diff(t))
  n <- length(t)

  isotropic <- CTMM$isotropic
  sigma <- CTMM$sigma
  COVM <- function(...) { covm(...,isotropic=CTMM$isotropic,axes=CTMM$axes) } # enforce model structure
  theta <- sigma@par['angle'] # NA in 1D

  STUFF <- squeezable.covm(CTMM)
  smgm <- STUFF$fact  # ratio of major axis to geometric mean axis
  ECC.EXT <- !STUFF$able # extreme eccentricity --- cannot squeeze data to match variances

  K <- max(length(CTMM$tau),1)

  circle <- CTMM$circle

  ####################################
  # PRECOMPUTE AVOIDS WASTEFUL ROTATIONS & TRANSFORMATIONS
  STUFF <- c("z","ROTATE","SQUEEZE","error","DIM","R")
  if(precompute>=0) # calculate new
  {
    # get the error information
    error <- CTMM$error.mat # note for fitted errors, this is error matrix @ UERE=1 (CTMM$error)
    class <- CTMM$class.mat
    ELLIPSE <- attr(error,"ellipse") # do we need error ellipses?
    TYPE <- DOP.match(CTMM$axes)
    UERE.DOF <- attr(data,"UERE")$DOF[,TYPE]
    names(UERE.DOF) <- rownames(attr(data,"UERE")$DOF)
    UERE.FIT <- CTMM$error>0 & !is.na(UERE.DOF) & UERE.DOF<Inf # will we be fitting any error parameters?

    # don't try to fit error class parameters absent from data
    if(any(CTMM$error>0) && "class" %in% names(data))
    {
      LEVELS <- levels(data$class)
      UERE.DOF <- UERE.DOF[LEVELS]
      UERE.FIT <- UERE.FIT[LEVELS]
      # CTMM$error <- CTMM$error[LEVELS]
    }

    # are we fitting the error, then the above is not yet normalized.
    if(any(UERE.FIT)) # calibrate errors
    {
      class <- c( class %*% CTMM$error^2 )
      error[] <- class * error
    }
    rm(class)

    if(ELLIPSE || (!isotropic && circle && any(CTMM$error>0)) || (ECC.EXT && AXES>1)) { DIM <- 2 } # requires 2D smoother
    else if(!isotropic & any(CTMM$error>0)) { DIM <- 1/2 } # requires 2x1D smoothers
    else { DIM <- 1 } # can use 1x1D smoother

    z <- get.telemetry(data,CTMM$axes)
    # u <- CTMM$mean.vec
    # mu <- CTMM$mu

    # orient the data along the major and minor axes of sigma
    ROTATE <- !isotropic && !ECC.EXT
    if(ROTATE)
    {
      z <- rotate.vec(z,-theta)
      sigma <- rotate.covm(sigma,-theta)

      if(ELLIPSE) { error <- rotate.mat(error,-theta) } # rotate error ellipses
    }

    # squeeze from ellipse to circle
    SQUEEZE <- !isotropic && (DIM<2 || circle) && !ECC.EXT
    if(SQUEEZE)
    {
      z <- squeeze(z,smgm)
      sigma <- squeeze.covm(sigma,circle=TRUE)

      if(any(CTMM$error>0)) { error <- squeeze.mat(error,smgm) } # squeeze error circles into ellipses
    }

    if(circle) ## COROTATING FRAME FOR circle=TRUE ##
    {
      R <- rotates(-circle*(t-t[1])) # rotation matrices
      z <- rotates.vec(z,R)
      if(ELLIPSE || (any(CTMM$error>0) && SQUEEZE)) { error <- rotates.mat(error,R) }
      # prepare R for inverse transformation
      R <- aperm(R,c(1,3,2))
    }
    else
    { R <- NULL }

    # fix variances of empty timestamps - set from fill.data
    if(!residual)
    {
      empty <- which(!data$record)
      if(length(empty)) { error[empty,,] <- aperm( array(diag(Inf,dim(error)[2]),c(dim(error)[2:3],length(empty))) ,c(3,1,2)) } # R is weird
      rm(empty)
    }

    # in case of SQUEEZE & rotate
    CTMM$sigma <- sigma
  }
  else # pull from old
  { for(thing in STUFF) { assign(thing,get(thing,pos=Kalman.env)) } }

  # store for later
  if(precompute>0) { for(thing in STUFF) { assign(thing,get(thing),pos=Kalman.env) } }
  # END PRECOMPUTE
  ################################

  if(!residual)
  {
    COV <- array(0,dim=c(n,AXES,AXES)) # position covariance
    if(K>1)
    {
      v <- array(0,dim=c(n,AXES))
      vCOV <- array(0,dim=c(n,AXES,AXES)) # velocity covariance
    }
  }

  SIGMA <- eigenvalues.covm(sigma)

  # rotated data with circular errors - 2x1D kalman smoother
  if(DIM==1/2) # diagonalize data and then run two 1D Kalman filters with separate means
  {
    # major axis likelihood
    CTMM$sigma <- SIGMA[1]
    KALMAN1 <- kalman(z[,1,drop=FALSE],u=NULL,t=t,dt=dt,CTMM=CTMM,error=error[,1,1,drop=FALSE],precompute=precompute,sample=sample,residual=residual,...)
    # minor axis likelihood
    CTMM$sigma <- SIGMA[2]
    KALMAN2 <- kalman(z[,2,drop=FALSE],u=NULL,t=t,dt=dt,CTMM=CTMM,error=error[,2,2,drop=FALSE],precompute=precompute,sample=sample,residual=residual,...)

    if(residual) { return(cbind(KALMAN1,KALMAN2)) }

    z[,1] <- KALMAN1$Z[,1,]
    z[,2] <- KALMAN2$Z[,1,]
    if(!sample)
    {
      COV[,1,1] <- KALMAN1$S[,1,1]
      COV[,2,2] <- KALMAN2$S[,1,1]
    }
    if(K>1)
    {
      v[,1] <- KALMAN1$Z[,2,]
      v[,2] <- KALMAN2$Z[,2,]
      if(!sample)
      {
        vCOV[,1,1] <- KALMAN1$S[,2,2]
        vCOV[,2,2] <- KALMAN2$S[,2,2]
      }
    }
  }
  else # use 1 Kalman filter - may be 1D or 2D
  {
    if(DIM==1)
    {
      CTMM$sigma <- SIGMA[1] # isotropic variance
      error <- error[,1,1,drop=FALSE] # isotropic && UERE redundant error information
    }

    KALMAN <- kalman(z,u=NULL,t=t,dt=dt,CTMM=CTMM,error=error,DIM=DIM,precompute=precompute,sample=sample,residual=residual,...)
    # point estimates will be correct but eccentricity is missing from variances

    if(residual) { return(KALMAN) }

    # position and velocity entries
    POS <- VEL <- array(FALSE,c(K,DIM))
    POS[1,] <- TRUE ; POS <- c(POS)
    if(K>1) { VEL[2,] <- TRUE ; VEL <- c(VEL) }

    z <- cbind(KALMAN$Z[,POS,])
    if(!sample) { COV <- KALMAN$S[,POS,POS,drop=FALSE] }
    if(K>1)
    {
      v <- cbind(KALMAN$Z[,VEL,])
      if(!sample) { vCOV <- KALMAN$S[,VEL,VEL,drop=FALSE] }
    }

    if(DIM<AXES && !sample) # promote from VAR to COV (2,2)
    {
      # fix for circular smoother # keeps track of SQUEEZE and ROTATE
      MAT <- diag(AXES)
      COV <- drop(COV) %o% MAT
      if(K>1) { vCOV <- drop(vCOV) %o% MAT }
    }
  } # end single filter

  if(circle) # circulate
  {
    z <- rotates.vec(z,R)
    if(!sample) { COV <- rotates.mat(COV,R) }
    if(K>1)
    {
      # includes non-inertial frame component: Omega x r
      v <- rotates.vec(v,R) + circle*cbind(-z[,2],z[,1])
      if(!sample) { vCOV <- rotates.mat(vCOV,R) }
    }
  }

  if(SQUEEZE) # unsqueeze the distribution
  {
    z <- squeeze(z,1/smgm)
    if(!sample) { COV <- squeeze.mat(COV,1/smgm) }
    if(K>1)
    {
      v <- squeeze(v,1/smgm)
      if(!sample) { vCOV <- squeeze.mat(vCOV,1/smgm) }
    }
  }

  if(ROTATE) # transform results back
  {
    z <- rotate.vec(z,+theta)
    if(!sample) { COV <- rotate.mat(COV,theta) }
    if(K>1)
    {
      v <- rotate.vec(v,theta)
      if(!sample) { vCOV <- rotate.mat(vCOV,theta) }
    }
  }

  colnames(z) <- CTMM$axes
  dimnames(COV) <- list(NULL,colnames(z),colnames(z))
  RETURN <- list(t=t,R=z,COV=COV)
  if(K>1)
  {
    colnames(v) <- paste0('v',CTMM$axes)
    dimnames(vCOV) <- list(NULL,colnames(v),colnames(v))
    RETURN$V <- v
    RETURN$VCOV <- vCOV
  }

  return(RETURN)
}


########################################
# fill in data gaps with missing observations of infinite error
########################################
fill.data <- function(data,CTMM=ctmm(tau=Inf),verbose=FALSE,t=NULL,dt=NULL,res=1,cor.min=0,dt.max=NULL,DT=diff(t),buffer=FALSE)
{
  DT.MAX <- dt.max # store for later

  # is this recorded data or empty gap
  data$record <- TRUE
  # repeated timestamps to skip in occurrence
  data$skip <- FALSE
  data$skip[diff(data$t)==0] <- TRUE

  if(is.null(t) && (!length(CTMM$tau) || CTMM$tau[1]==0)) { t <- data$t } # don't add further times

  # FIX THE TIME GRID TO AVOID TINY DT
  if(is.null(t))
  {
    t <- data$t
    if(is.null(DT)) { DT <- diff(t) }

    # target resolution
    if(is.null(dt)){ dt <- stats::median(DT)/res }

    # can cor.min argument be applied
    if(length(CTMM$tau)<2 || CTMM$tau[2]<=0) { cor.min <- FALSE }

    if(cor.min)
    {
      # convert from correlation to time
      cor.min <- -log(cor.min)*CTMM$tau[2] # need for buffer=TRUE
      # maximum gap to bridge
      dt.max <- max(DT.MAX,cor.min)
    }
    if(is.null(dt.max)) { dt.max <- Inf } # default don't skip gaps
    dt.max2 <- dt.max/2

    # this regularization is not perfectly regular, but holds up to sampling drift in caribou data
    t.grid <- c()  # full (locally) even grid
    dt.grid <- c() # local (numeric) sampling resolution
    t.new <- c()   # new times in this even grid
    for(i in which(DT>0)) # don't repeat timestamps, even if data does, unless dt changes
    {
      if(DT[i] <= dt.max)
      {
        n.sub <- round(DT[i]/dt)+1
        n.sub <- max(n.sub,2) # fix for crazy small time-steps
        t.sub <- seq(from=t[i],to=t[i+1],length.out=n.sub)
        dt.sub <- DT[i]/(n.sub-1)
        dt.sub <- rep(dt.sub,n.sub)
      }
      else # skip bulk of gap
      {
        t.sub <- seq(from=0,to=dt.max2,by=dt)
        t.sub <- c( t[i]+t.sub , t[i+1]-rev(t.sub) )
        n.sub <- length(t.sub)
        dt.sub <- rep(dt,n.sub)
      }

      t.grid <- c(t.grid,t.sub)
      dt.grid <- c(dt.grid,dt.sub)

      t.new <- c(t.new,t.sub[c(-1,-n.sub)])
    }

    # buffer observation period
    dt.max <- min(cor.min,DT.MAX) * buffer
    if(dt.max)
    {
      if(cor.min==Inf) { stop("buffer=TRUE incompatible with cor.min=0.") }
      dt.max2 <- dt.max/2

      dt.buffer <- seq(0,dt.max2,by=dt)
      buffer <- rev( t.grid[1] - dt.buffer )
      t.grid <- c(buffer,t.grid)
      t.new <- c(buffer[-length(buffer)],t.new)

      buffer <- last(t.grid) + dt.buffer
      t.grid <- c(t.grid,buffer)
      t.new <- c(t.new,buffer[-1])

      dt.buffer <- array(dt,length(dt.buffer))
      dt.grid <- c(dt.buffer,dt.grid,dt.buffer)
    }

    # don't need to repeat if dt doesn't change
    SAME <- diff(t.grid)==0 & diff(dt.grid)==0
    SAME <- c(SAME,FALSE) # keep last time
    t.grid <- t.grid[!SAME] # drop first same
    dt.grid <- dt.grid[!SAME] # drop first same

    # half weight repeated times
    w.grid <- dt.grid
    REPEAT <- which(diff(t.grid)==0)
    w.grid[REPEAT] <- w.grid[REPEAT]/2
    w.grid[REPEAT+1] <- w.grid[REPEAT+1]/2
  } # end if is.null(t)
  else # use a pre-specified time grid
  {
    t.new <- t[!(t %in% data$t)]
    t.grid <- t
    dt.grid <- diff(t)
    dt.grid <- pmin(c(Inf,dt.grid),c(dt.grid,Inf))
    w.grid <- rep(1,length(t))
  }

  # empty observation row for these times
  blank <- data[1,]
  blank$record <- FALSE # these are not TRUE records
  blank$skip <- FALSE # don't skip even if timestamp repeats
  blank <- blank[rep(1,length(t.new)),]
  blank$t <- t.new

  # attach empty measurements to data
  data <- rbind(data,blank)

  # sort times
  data <- data[sort.list(data$t,na.last=NA),]
  # this is now our fake data set to feed into the kalman smoother

  if(verbose) { data <- list(data=data,t.grid=t.grid,dt.grid=dt.grid,w.grid=w.grid) }

  return(data)
}


##############################################
# SIMULATE DATA over time array t
simulate.ctmm <- function(object,nsim=1,seed=NULL,data=NULL,VMM=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,precompute=FALSE,...)
{
  T.SPECIFIED <- !is.null(t)

  info <- attr(object,"info")
  if(!is.null(data)) { info$identity <- glue( attr(data,'info')$identity , info$identity ) }

  # have to do this becaues simulate is an S3 with 3 fixed arguments

  if(class(object)[1] %in% c("data.frame","telemetry"))
  {
    TEMP <- data
    data <- object
    object <- TEMP
  }

  if(class(nsim)[1] %in% c("data.frame","telemetry"))
  {
    TEMP <- data
    data <- nsim
    nsim <- TEMP
  }

  if(class(nsim)[1] == "ctmm")
  {
    TEMP <- object
    object <- nsim
    nsim <- TEMP
  }

  if(class(seed)[1] %in% c("data.frame","telemetry"))
  {
    TEMP <- data
    data <- seed
    seed <- TEMP
  }

  if(class(seed)[1] == "ctmm")
  {
    TEMP <- object
    object <- seed
    seed <- TEMP
  }

  if(is.null(nsim)) { nsim <- 1 }

  if(nsim>1)
  {
    S <- lapply(1:nsim,function(i){simulate.ctmm(object,seed,data,VMM=VMM,t=t,dt=dt,res=res,complete=complete,precompute=precompute,nsim=1,...)})
    return(S)
  }

  if(is.null(object) && !is.null(VMM)) # 1D
  {
    object <- VMM
    VMM <- NULL
  }
  else if(!is.null(object) && !is.null(VMM)) # 3D
  {
    # combine results (lazy code, calculates time grid twice)
    OUT <- simulate.ctmm(object,nsim=nsim,seed=seed,data=data,t=t,dt=dt,res=res,...)
    ZOUT <- simulate.ctmm(VMM,nsim=nsim,seed=seed,data=data,t=t,dt=dt,res=res,...)
    data <- cbind(OUT,ZOUT[,-1]) # drop redundant time column
    rm(OUT,ZOUT)

    data <- new.telemetry(data,info=info)
    if(complete) { data <- pseudonymize(data,tz=info$timezone,proj=info$projection,origin=EPOCH) }

    attr(data,"UERE") <- uere.null(data)
    attr(data,"UERE")$UERE[] <- 0
    attr(data,"UERE")$DOF[] <- Inf
    attr(data,"UERE")$N[] <- Inf

    return(data)
  }
  # 1-2D below

  axes <- object$axes
  AXES <- length(axes)

  # no movement model
  if(is.null(object)) { object <- ctmm(sigma=0,mu=rep(0,AXES),error=TRUE) }
  ZERO <- all(diag(object$sigma)==0) # no movement model logical

  if(!is.null(seed)){ set.seed(seed) }

  CLASS <- class(data)[1]
  CONDITIONAL <- FALSE
  if(CLASS=="telemetry" || CLASS=='data.frame')
  {
    if(!ZERO && all(object$axes %in% names(data))) # condition off of data
    { CONDITIONAL <- TRUE }
    else # use data time & error for unconditional simulation
    { t <- data$t }
  }

  if(CONDITIONAL)
  {
    STUFF <- c('object','data','drift','velocity')
    if(precompute>=0) # prepare model and data frame
    {
      object <- ctmm.prepare(data,object,precompute=FALSE) # u calculated here with unfilled t
      data <- fill.data(data,CTMM=object,t=t,dt=dt,res=res,...)
      # object$error <- TRUE # avoids unit variance algorithm - data contains fixed errors from fill.data

      # calculate trend
      velocity <- drift.velocity(object,data$t) %*% object$mu
      drift <- drift.mean(object,data$t) %*% object$mu

      # detrend for simulation - retrend later
      z <- get.telemetry(data,axes=object$axes)
      data[,object$axes] <- z - drift
    }
    else # recycle model and data frame
    { for(thing in STUFF){ assign(thing,get(thing,pos=Kalman.env)) } }

    # store prepared model and data frame
    if(precompute>0) { for(thing in STUFF){ assign(thing,get(thing),pos=Kalman.env) } }

    data <- smoother(data,object,sample=TRUE,precompute=precompute)

    # retrend data
    data$R <- data$R + drift
    rm(drift)
    # trend velocity
    if("V" %in% names(data)) { data$V <- data$V + velocity }

    data <- cbind(t=data$t,data$R,data$V)
    data <- data.frame(data)

    # the user probably only wants times t if t is specified
    if(T.SPECIFIED)
    {
      WHICH <- data$t %in% t # I'm assuming R is coded to do a respectable sort match
      data <- data[WHICH,]
    }

  } # conditional simulation
  else # Gaussian simulation not conditioned off of any data
  {
    STUFF <- c('Green','Sigma','error','ELLIPSE','object','mu','Lambda','n','K','z','v','circle','R')
    if(precompute>=0)
    {
      if(is.null(data))
      { error <- ELLIPSE <- FALSE }
      else # get error if provided
      {
        error <- get.error(data,object,calibrate=TRUE)
        ELLIPSE <- attr(error,'ellipse')
      }
      object <- ctmm.prepare(data.frame(t=t),object) # mean.vec calculated here

      n <- length(t)

      if(ZERO) # no model - simulation of errors only
      {
        z <- get.telemetry(data)
        v <- NULL
        Green <- Sigma <- mu <- Lambda <- K <- circle <- R <- NULL # make sure get() doesn't fail
      }
      else # have model
      {
        tau <- object$tau
        if(length(tau)==0) { tau = 0 }
        K <- length(tau)

        mu <- object$mu
        if(is.null(mu)) { mu <- array(0,c(1,AXES)) }
        sigma <- object$sigma
        if(is.null(sigma)) { sigma <- covm(1,axes=axes) }

        Lambda <- sqrtm.covm(sigma)

        K <- length(tau)

        dt <- c(Inf,diff(t)) # time lags

        # where we will store the data
        z <- array(0,c(n,AXES))
        if(K>1) { v <- array(0,c(n,AXES)) } else { v <- NULL }

        Green <- array(0,c(n,K,K))
        Sigma <- array(0,c(n,K,K))

        object$sigma <- 1
        object <- get.taus(object) # pre-compute stuff for Langevin equation solutions
        for(i in 1:n)
        {
          # tabulate propagators if necessary
          if((i==1)||(dt[i]!=dt[i-1])) { Langevin <- langevin(dt=dt[i],CTMM=object) }
          Green[i,,] <- Langevin$Green
          Sigma[i,,] <- Langevin$Sigma
        }
        if(!object$range) { Sigma[1,,] <- 0 } # start at first point instead of random point on Earth

        # Sigma is now standardization matrix
        Sigma <- vapply(1:n,function(i){PDfunc(Sigma[i,,],func=function(x){sqrt(abs(x))},pseudo=TRUE)},Sigma[1,,]) # (K,K,n)
        dim(Sigma) <- c(K,K,n)
        Sigma <- aperm(Sigma,c(3,1,2)) # (n,K,K)

        # circulation stuff
        circle <- object$circle
        R <- exp(1i*circle*(t-t[1]))
      } # end if(!is.null(object))

      # pre-compute error matrices
      if(any(object$error>0) && !ELLIPSE) # circular errors
      { error <- sqrt(error) }
      else if(ELLIPSE) # eliptical errors
      {
        error <- vapply(1:n, function(i){sqrtm(error[i,,])}, diag(2)) # (2,2,n)
        error <- aperm(error,c(3,1,2)) # (n,2,2)
      }
    } # END precompute
    else # precomputed objects from previous run
    { for(thing in STUFF) { assign(thing,get(thing,pos=Kalman.env)) } }

    # store precomputed objects for later
    if(precompute>0) { for(thing in STUFF) { assign(thing,get(thing),pos=Kalman.env) } }

    if(!ZERO) # have model
    {
      # initial hidden state, for standardized process
      H <- array(0,c(K,AXES))
      for(i in 1:n)
      {
        # generate standardized process - R arrays are something awful... awful
        H[] <- cbind(Green[i,,]) %*% H[,,drop=FALSE] + cbind(Sigma[i,,]) %*% array(stats::rnorm(K*AXES),c(K,AXES))
        # pull out location from hidden state
        z[i,] <- H[1,]
        if(K>1) { v[i,] <- H[2,] }
      }

      # rotate process if necessary
      if(circle)
      {
        z <- z[,1] + 1i*z[,2]
        z <- R * z

        if(K>1)
        {
          v <- v[,1] + 1i*v[,2]
          v <- (R*v) + 1i*circle*z # mean-zero z
          v <- cbind(Re(v),Im(v))
        }

        z <- cbind(Re(z),Im(z))
      }

      # calculate mean function
      z <- (z %*% Lambda) + (object$mean.vec %*% mu)
      colnames(z) <- axes
      if(K>1)
      {
        v <- (v %*% Lambda) + (drift.velocity(object,t) %*% mu)
        colnames(v) <- paste0("v",axes)
      }
    } # end if(!is.null(object))
    else # no process variance
    { z <- (object$mean.vec %*% mu) }

    # throw in error
    if(any(object$error>0))
    {
      if(!ELLIPSE) # circular errors
      { error <- error * array(stats::rnorm(n*length(axes)),c(n,length(axes))) }
      else # elliptical errors # can we do this with one 2n column product?
      {
        error <- vapply(1:n, function(i){error[i,,] %*% stats::rnorm(2)}, c(0,0) )
        error <- t(error)
      }
      z[] <- z + error # error became (dim,n)
      # velocity error?
    }

    # restore error columns if we simulated error
    if(is.null(data) || !any(object$error>0))
    {
      data <- cbind(t=t,z,v)
      data <- data.frame(data)
    }
    else
    {
      data[,axes] <- z
      # not storing velocity without error yet!
    }
  } # Gaussian simulation

  data <- new.telemetry(data,info=info)
  if(complete)
  {
    if(all(axes=='z')) { stop("(x,y) locations must also be simulated for complete=TRUE.") }
    data <- pseudonymize(data,tz=info$timezone,proj=info$projection,origin=EPOCH)
  }

  attr(data,"UERE") <- uere.null(data)
  attr(data,"UERE")$UERE[] <- 0
  attr(data,"UERE")$DOF[] <- Inf
  attr(data,"UERE")$N[] <- Inf

  return(data)
}
#methods::setMethod("simulate",signature(object="ctmm"), function(object,...) simulate.ctmm(object,...))


simulate.telemetry <- function(object,nsim=1,seed=NULL,CTMM=NULL,VMM=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,precompute=FALSE,...)
{ simulate.ctmm(CTMM,nsim=nsim,seed=seed,data=object,VMM=VMM,t=t,dt=dt,res=res,complete=complete,precompute=precompute,...) }


##########################
# predict locations at certain times !!! make times unique
##########################
predict.ctmm <- function(object,data=NULL,VMM=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,...)
{
  info <- attr(object,"info")
  if(!is.null(data)) { info$identity <- glue( attr(data,'info')$identity , info$identity ) }

  if(is.null(object) && !is.null(VMM)) # 1D
  {
    object <- VMM
    VMM <- NULL
  }
  else if(!is.null(object) && !is.null(VMM)) # 3D
  {
    # combine results (lazy code, calculates time grid twice)
    OUT <- predict.ctmm(object,data=data,t=t,dt=dt,res=res,...)
    ZOUT <- predict.ctmm(VMM,data=data,t=t,dt=dt,res=res,...)
    data <- cbind(OUT,ZOUT[,-1]) # drop redundant time column
    rm(OUT,ZOUT)

    data <- new.telemetry(data,info=info)
    if(complete) { data <- pseudonymize(data,tz=info$timezone,proj=info$projection,origin=EPOCH) }

    attr(data,"UERE") <- uere.null(data)
    attr(data,"UERE")$UERE[] <- 1
    attr(data,"UERE")$DOF[] <- Inf
    attr(data,"UERE")$N[] <- Inf

    return(data)
  }
  # 1-2D below

  axes <- object$axes

  # Gaussian simulation not conditioned off of any data
  if(is.null(data))
  {
    object <- ctmm.prepare(data.frame(t=t),object)

    mu <- object$mu
    if(is.null(mu)) { mu <- array(0,c(1,length(axes))) }

    # calculate mean function
    r <- object$mean.vec %*% mu
    colnames(r) <- axes

    v <- drift.velocity(object,t) %*% mu
    colnames(v) <- paste0("v",axes)

    data <- data.frame(r,v)
    data$t <- t

    # missing COVs !!!
    DOP <- DOP.match(axes)
    sigma <- methods::getDataPart(object$sigma)
    if(length(axes)==1 || object$isotropic)
    {
      if(length(sigma)>1) { sigma <- mean(diag(sigma)) }
      data[[DOP.LIST[[DOP]]$VAR]] <- sigma
      if(length(object$tau)>1)
      {
        if(length(axes)==1)
        { data[,paste0("VAR.v",axes)] <- sigma/prod(object$tau) }
        else
        { data[,DOP.LIST$speed$VAR] <- sigma/prod(object$tau) }
      }
    }
    else
    {
      sigma <- c(sigma)[-3]
      for(i in 1:3)
      {
        data[[DOP.LIST[[DOP]]$COV[i]]] <- sigma[i]
        if(length(object$tau)>1) { data[,DOP.LIST$speed$COV[i]] <- sigma[i]/prod(object$tau) }
      }
    }
  }
  else # condition off of the data
  {
    # object <- ctmm.prepare(data,object) # mean.vec here is calculated with pre-filled t
    K <- length(object$tau)
    data <- fill.data(data,CTMM=object,t=t,dt=dt,res=res)
    # object$error <- TRUE # avoids unit variance algorithm

    # calculate trend
    velocity <- drift.velocity(object,data$t) %*% object$mu
    drift <- drift.mean(object,data$t) %*% object$mu

    # detrend for simulation - retrend later
    z <- get.telemetry(data,axes=axes)
    data[,axes] <- z - drift

    # smooth mean-zero data
    data <- smoother(data,object,smooth=TRUE)

    # detrend for simulation - retrend later
    data$R <- data$R + drift
    # trend velocity
    if("V" %in% names(data)) { data$V <- data$V + velocity }

    NAMES <- colnames(data$R)
    COV <- data$COV

    if(K>1)
    {
      VNAMES <- colnames(data$V)
      VCOV <- data$VCOV
    }

    data <- cbind(t=data$t,data$R,data$V)
    data <- data.frame(data)

    # flatten covariance matrix and include in data.frame
    for(i in 1:length(NAMES))
    {
      for(j in i:length(NAMES))
      {
        NAME <- paste("COV.",NAMES[i],".",NAMES[j],sep="") # consistent with imported ARGOS error ellipse notation
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
          NAME <- paste("COV.",VNAMES[i],".",VNAMES[j],sep="") # consistent with above
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
  if(complete)
  {
    if(all(axes=='z')) { stop("(x,y) locations must also be predicted for complete=TRUE.") }
    data <- pseudonymize(data,tz=info$timezone,proj=info$projection,origin=EPOCH)
  }

  attr(data,"UERE") <- uere.null(data)
  attr(data,"UERE")$UERE[] <- 1
  attr(data,"UERE")$DOF[] <- Inf
  attr(data,"UERE")$N[] <- Inf

  return(data)
}

predict.telemetry <- function(object,CTMM=NULL,VMM=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,...)
{ predict.ctmm(CTMM,data=object,VMM=VMM,t=t,dt=dt,res=res,complete=complete,...) }
