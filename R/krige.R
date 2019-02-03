# occurrence <- function(data,CTMM,H=0,res.time=10,res.space=10,grid=NULL,cor.min=0.5,dt.max=NULL) UseMethod("overlap") #S3 generic

################################
# Return hidden state estimates or simulations
################################
smoother <- function(DATA,CTMM,precompute=FALSE,sample=FALSE,residual=FALSE,...)
{
  if(is.null(CTMM$error.mat)) { CTMM <- ctmm.prepare(DATA,CTMM) }
  if(is.null(DATA$record)) { DATA$record <- TRUE } # real recorded data or blank/empty timestamps from fill-data
  AXES <- length(CTMM$axes)

  t <- DATA$t

  dt <- c(Inf,diff(t))
  n <- length(t)

  isotropic <- CTMM$isotropic
  sigma <- CTMM$sigma
  area <- sigma@par[1]
  COVM <- function(...) { covm(...,isotropic=CTMM$isotropic,axes=CTMM$axes) } # enforce model structure
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

  K <- max(length(CTMM$tau),1)

  circle <- CTMM$circle

  ####################################
  # PRECOMPUTE AVOIDS WASTEFUL ROTATIONS & TRANSFORMATIONS
  STUFF <- c("z","ROTATE","SQUEEZE","error","UERE","DIM","R")
  if(precompute>=0) # calculate new
  {
    # get the error information
    error <- CTMM$error.mat # note for fitted errors, this is error matrix @ UERE=1 (CTMM$error)
    UERE <- attr(error,"flag")
    # are we fitting the error, then the above is not yet normalized
    if(UERE && UERE<3) { error <- error * CTMM$error } # set UERE value if necessary

    if(UERE>=4 || (!isotropic && circle && UERE)) { DIM <- 2 } # requires 2D smoother
    else if(!isotropic & UERE) { DIM <- 1/2 } # requires 2x1D smoothers
    else { DIM <- 1 } # can use 1x1D smoother

    z <- get.telemetry(DATA,CTMM$axes)
    # u <- CTMM$mean.vec
    # mu <- CTMM$mu

    # do we need to orient the data along the major an minor axes of sigma for 1D smoothers (or elliptical circulation)
    ROTATE <- !isotropic && (circle || (UERE && UERE<4))
    if(ROTATE)
    {
      z <- rotate.vec(z,-theta)

      sigma <- attr(sigma,"par")
      sigma["angle"] <- 0
      sigma <- COVM(sigma)

      if(UERE>=4) { error <- rotate.mat(error,-theta) } # rotate error ellipses
    }

    # squeeze from ellipse to circle for circulation model transformation
    SQUEEZE <- (!isotropic && circle)
    if(SQUEEZE)
    {
      z <- squeeze(z,ecc)

      sigma <- attr(sigma,"par")
      sigma["eccentricity"] <- 0
      sigma <- COVM(sigma)

      if(UERE) { error <- squeeze.mat(error,ecc) } # squeeze error circles into ellipses
    }

    if(circle) ## COROTATING FRAME FOR circle=TRUE ##
    {
      R <- rotates(-circle*(t-t[1])) # rotation matrices
      z <- rotates.vec(z,R)
      if(UERE>=4 || (UERE && SQUEEZE)) { error <- rotates.mat(error,R) }
      # prepare R for inverse transformation
      R <- aperm(R,c(1,3,2))
    }
    else
    { R <- NULL }

    # fix variances of empty timestamps - set from fill.data
    if(!residual)
    {
      empty <- which(!DATA$record)
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

  # rotated data with circular errors - 2x1D kalman smoother
  if(DIM==1/2) # diagonalize data and then run two 1D Kalman filters with separate means
  {
    # major axis likelihood
    CTMM$sigma <- area * exp(+ecc/2) # major variance
    KALMAN1 <- kalman(z[,1,drop=FALSE],u=NULL,dt=dt,CTMM=CTMM,error=error,precompute=precompute,sample=sample,residual=residual,...)
    # minor axis likelihood
    CTMM$sigma <- area * exp(-ecc/2) # minor variance
    KALMAN2 <- kalman(z[,2,drop=FALSE],u=NULL,dt=dt,CTMM=CTMM,error=error,precompute=precompute,sample=sample,residual=residual,...)

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
    if(DIM==1) { CTMM$sigma <- area } # isotropic variance
    KALMAN <- kalman(z,u=NULL,dt=dt,CTMM=CTMM,error=error,precompute=precompute,sample=sample,residual=residual,...)
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
      # fix for zero-error eccentric smoother
      MAT <- attr(sigma,'par') # keeps track of SQUEEZE and ROTATE
      MAT[1] <- 1
      MAT <- covm(MAT)
      MAT <- methods::getDataPart(MAT)
      COV <- drop(COV) %o% MAT
      if(K>1) { vCOV <- drop(vCOV) %o% MAT }
    }
  }

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
    z <- squeeze(z,-ecc)
    if(!sample) { COV <- squeeze.mat(COV,-ecc) }
    if(K>1)
    {
      v <- squeeze(v,-ecc)
      if(!sample) { vCOV <- squeeze.mat(vCOV,-ecc) }
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
fill.data <- function(data,CTMM=ctmm(tau=Inf),verbose=FALSE,t=NULL,dt=NULL,res=1,cor.min=0,dt.max=NULL,DT=diff(t))
{
  # is this recorded data or empty gap
  data$record <- TRUE

  if(!length(CTMM$tau) || CTMM$tau[1]==0) { t <- data$t } # don't add further times

  # FIX THE TIME GRID TO AVOID TINY DT
  if(is.null(t))
  {
    t <- data$t
    if(is.null(DT)) { DT <- diff(t) }

    # target resolution
    if(is.null(dt)){ dt <- stats::median(DT)/res }

    # maximum gap to bridge
    if(is.null(dt.max)) { dt.max <- -log(cor.min)*CTMM$tau[1] }
    dt.max2 <- dt.max/2

    # this regularization is not perfectly regular, but holds up to sampling drift in caribou data
    t.grid <- c()  # full (locally) even grid
    dt.grid <- c() # local (numeric) sampling resolution
    t.new <- c()   # new times in this even grid
    for(i in 1:length(DT))
    {
      if(DT[i]<=dt.max)
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

    # half weight repeated endpoints in grid
    w.grid <- dt.grid
    REPEAT <- which(diff(t.grid)==0)
    w.grid[REPEAT] <- w.grid[REPEAT]/2
    w.grid[REPEAT+1] <- w.grid[REPEAT+1]/2
  }
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
  blank <- blank[rep(1,length(t.new)),]
  blank$t <- t.new

  # attach empty measurements to data
  data <- rbind(data,blank)

  # sort times
  data <- data[sort.list(data$t,na.last=NA,method="quick"),]
  # this is now our fake data set to feed into the kalman smoother

  if(verbose) { data <- list(data=data,t.grid=t.grid,dt.grid=dt.grid,w.grid=w.grid) }

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
  CTMM0 <- CTMM
  dt <- stats::median(diff(data$t))

  if(length(H)==1) { H <- diag(H,2) }

  info <- attr(data,"info")
  SIGMA <- CTMM$sigma # diffusion matrix for later
  CTMM <- ctmm.prepare(data,CTMM,precompute=FALSE) # not the final t for calculating u
  error <- get.error(data,CTMM,circle=TRUE)
  MIN.ERR <- min(error) # Fix something here?

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
  COV <- vapply(1:n,function(i){ PDclamp(COV[i,,]) },COV[1,,]) # (d,d,n)
  COV <- aperm(COV,perm=c(3,1,2)) # (n,d,d)

  # uncertainties/bandwidths for this data
  H <- vapply(1:n,function(i){ H + COV[i,,] + dt.grid[i]^2/12*(V[i,] %o% V[i,]) },H) # (2,2,n)
  H <- aperm(H,c(3,1,2)) # (n,2,2)

  # estimate size of data blob
  dr <- diag(SIGMA)
  if(CTMM$range && length(CTMM$tau)) { dr <- dr/CTMM$tau[1] }
  # prefactors from mid-bridge variance
  if(length(CTMM$tau)==1) #BM/OU
  { dr <- dt/4 * dr }
  else if(length(CTMM$tau)==2) #IOU/OUF
  { dr <- dt^2/24 * dr/CTMM$tau[2] }

  if(CTMM$error){ dr <- dr + MIN.ERR }
  dr <- sqrt(dr)

  # using the same data format as AKDE, but with only the ML estimate (alpha=1)
  KDE <- kde(data,H=H,W=w.grid,dr=dr/res.space)
  KDE$H <- diag(0,2)
  KDE <- new.UD(KDE,info=info,type='occurrence',CTMM=CTMM0)
  return(KDE)
}


##############################################
# SIMULATE DATA over time array t
simulate.ctmm <- function(object,nsim=1,seed=NULL,data=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,precompute=FALSE,...)
{
  if(!is.null(seed)){ set.seed(seed) }

  info <- attr(object,"info")
  if(!is.null(data)) { info$identity <- glue( attr(data,'info')$identity , info$identity ) }
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
    STUFF <- c('object','data','drift','velocity')
    if(precompute>=0) # prepare model and data frame
    {
      object <- ctmm.prepare(data,object,precompute=FALSE) # u calculated here with unfilled t
      data <- fill.data(data,CTMM=object,t=t,dt=dt,res=res,...)
      # object$error <- TRUE # avoids unit variance algorithm - data contains fixed errors from fill.data

      # calculate trend
      drift <- get(object$mean)
      velocity <- drift@velocity(data$t,object) %*% object$mu
      drift <- drift(data$t,object) %*% object$mu

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

  } # conditional simulation
  else # Gaussian simulation not conditioned off of any data
  {
    STUFF <- c('Green','Sigma','error','object','mu','Lambda','n','K','z','v','circle','R','UERE')
    if(precompute>=0)
    {
      if(is.null(data)) { error <- UERE <- FALSE }
      else { error <- get.error(data,object) ; UERE <- attr(error,'flag') } # get error if provided
      object <- ctmm.prepare(data.frame(t=t),object) # mean.vec calculated here

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
        # Sigma is now standardization matrix
        Sigma[i,,] <- PDfunc(Sigma[i,,],func=function(x){sqrt(abs(x))},pseudo=TRUE)
      }

      # circulation stuff
      circle <- object$circle
      R <- exp(1i*circle*(t-t[1]))

      # pre-compute error matrices
      if(UERE && UERE<=3) # circular errors
      { error <- sqrt(error) }
      else if(UERE) # eliptical errors
      {
        error <- vapply(1:n, function(i){sqrtm(error[i,,])}, diag(2)) # (2,2,n)
        error <- aperm(error,c(3,1,2)) # (n,2,2)
      }
    } # END precompute
    else # precomputed objects from previous run
    { for(thing in STUFF) { assign(thing,get(thing,pos=Kalman.env)) } }

    # store precomputed objects for later
    if(precompute>0) { for(thing in STUFF) { assign(thing,get(thing),pos=Kalman.env) } }

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
      z <- cbind(Re(z),Im(z))

      if(K>1)
      {
        v <- v[,] + 1i*v[,2]
        v <- (R*v) + 1i*circle*z # mean zero here
        v <- cbind(Re(v),Im(v))
      }
    }

    # calculate mean function
    z <- (z %*% Lambda) + (object$mean.vec %*% mu)
    colnames(z) <- axes
    if(K>1)
    {
      v <- (v %*% Lambda) + (get(object$mean)@velocity(t,object) %*% mu)
      colnames(v) <- paste0("v",axes)
    }

    # throw in error
    if(UERE)
    {
      if(UERE<=3) # circular errors
      { error <- error * array(stats::rnorm(n*length(axes)),c(n,length(axes))) }
      else # eliptical errors # can we do this with one 2n column product?
      {
        error <- vapply(1:n, function(i){error[i,,] %*% stats::rnorm(2)}, c(0,0) )
        error <- t(error)
      }
      z[] <- z + error # error became (dim,n)
      # velocity error?
    }

    # restore error columns if we simulated error
    if(is.null(data) || !UERE)
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
  if(complete) { data <- pseudonymize(data,tz=info$timezone,proj=info$projection,origin=EPOCH)  }
  return(data)
}
#methods::setMethod("simulate",signature(object="ctmm"), function(object,...) simulate.ctmm(object,...))

simulate.telemetry <- function(object,nsim=1,seed=NULL,CTMM=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,precompute=FALSE,...)
{ simulate.ctmm(CTMM,nsim=nsim,seed=seed,data=object,t=t,dt=dt,res=res,complete=complete,precompute=precompute,...) }


##########################
# predict locations at certaint times !!! make times unique
##########################
predict.ctmm <- function(object,data=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,...)
{
  info <- attr(object,"info")
  if(!is.null(data)) { info$identity <- glue( attr(data,'info')$identity , info$identity ) }
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

    v <- get(object$mean)@velocity(t,object) %*% mu
    colnames(v) <- paste0("v",axes)

    data <- data.frame(r,v)
    data$t <- t

    # missing COVs !!!
    DOP <- DOP.match(axes)
    sigma <- methods::getDataPart(object$sigma)
    if(length(axes)==1 || object$isotropic)
    {
      sigma <- mean(diag(sigma,length(axes)))
      data[[DOP.LIST[[DOP]]$VAR]] <- sigma
      if(length(object$tau)>1) { data[[paste0("VAR.v",axes)]] <- sigma/prod(object$tau) }
    }
    else
    {
      sigma <- c(sigma)[-3]
      for(i in 1:3)
      {
        data[[DOP.LIST[[DOP]]$COV[i]]] <- sigma[i]
        if(length(object$tau)>1) { data[[DOP.LIST$speed$COV[i]]] <- sigma[i]/prod(object$tau) }
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
    drift <- get(object$mean)
    velocity <- drift@velocity(data$t,object) %*% object$mu
    drift <- drift(data$t,object) %*% object$mu

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
  if(complete) { data <- pseudonymize(data,tz=info$timezone,proj=info$projection,origin=EPOCH)  }
  return(data)
}

predict.telemetry <- function(object,CTMM=NULL,t=NULL,dt=NULL,res=1,complete=FALSE,...)
{ predict.ctmm(CTMM,data=object,t=t,dt=dt,res=res,complete=complete,...) }
