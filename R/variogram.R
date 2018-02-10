# extend subset method
subset.variogram <- function(x,...)
{
  info <- attr(x,"info")
  x <- data.frame(x)
  x <- subset.data.frame(x,...)
  x < - droplevels(x) # why is this here?
  new.variogram(x,info=info)
}

`[.variogram` <- function(x,...)
{
  info <- attr(x,"info")
  x <- data.frame(x)
  x <- "[.data.frame"(x,...)
  # if(class(x)=="data.frame") { x <- new.variogram(x,info=info) }
  x <- new.variogram(x,info=info)
  return(x)
}

#########################
# variogram funcion wrapper
variogram <- function(data,dt=NULL,fast=TRUE,res=1,CI="Markov",axes=c("x","y"))
{
  res <- round(res)

  if(length(dt)>1)
  {
    dt <- sort(dt,method="quick")
    res <- rep(res,length(dt))
    res[-1] <- pmax(res[-1],dt[-1]/dt[-length(dt)])

    # calculate a variogram at each dt
    SVF <- lapply(1:length(dt), function(i) { variogram.dt(data,dt=dt[i],fast=fast,res=res[i],CI=CI,axes=axes) } )

    # subset each variogram to relevant range of lags
    dt <- c(-dt[1],dt,Inf)
    for(i in 1:(length(dt)-2))
    {
      SUB <- (dt[i] < SVF[[i]]$lag) & (SVF[[i]]$lag < dt[i+2])
      SVF[[i]] <- SVF[[i]][SUB,]
    }

    # concatentate
    SVF <- do.call(rbind,SVF)
  }
  else
  {
    SVF <- variogram.dt(data,dt=dt,res=res,fast=fast,CI=CI,axes=axes)
  }

  # delete missing lags
  SVF <- SVF[SVF$DOF>0,]

  # average error variance at UERE=1, so it can be rescaled by model fit UERE if necessary
  if("COV.angle" %in% names(data)) # elliptical error
  { error <- rowSums(get.telemetry(data,axes=c('COV.major','COV.minor')))/2 }
  else # circular error - possibly not calibrated
  { error <- get.error(data,list(error=TRUE,axes=axes)) }
  error <- mean(error) # mean variance of error

  info=attr(data,"info")
  info$axes <- axes
  info$error <- error

  SVF <- new.variogram(SVF,info=info)
  return(SVF)
}

################################
# wrapper for fast and slow variogram codes, for a specified dt
variogram.dt <- function(data,dt=NULL,fast=NULL,res=1,CI="Markov",axes=c("x","y"))
{
  # intelligently select algorithm
  if(is.null(fast))
  {
    if(length(data$t)<1000) { fast <- FALSE }
    else { fast <- TRUE }
  }

  if(fast)
  { SVF <- variogram.fast(data=data,dt=dt,res=res,CI=CI,axes=axes) }
  else
  { SVF <- variogram.slow(data=data,dt=dt,CI=CI,axes=axes) }

  # skip missing data
  SVF <- SVF[which(SVF$DOF>0),]
  SVF <- stats::na.omit(SVF)
  return(SVF)
}

############################
# best initial time for a uniform grid
grid.init <- function(t,dt=stats::median(diff(t)),W=NULL)
{
  if(is.null(W)) { W <- array(1,length(t)) }

  # simple analytic periodic cost function
  # COST = sum_t w(t) sin^2(pi(t-t0)/dt)
  # maximized anaytically
  theta <- (2*pi/dt)*t
  SIN <- c(W %*% sin(theta))
  COS <- c(W %*% cos(theta))
  t0 <- -dt/(2*pi)*atan(SIN/COS)

  return(t0)
}

############################
# preliminary objects to smear data across a uniform grid
pregridder <- function(t,dt=NULL,W=NULL)
{
  n <- length(t)

  # default time step
  if(is.null(dt)) { dt <- stats::median(diff(t)) }

  # choose best grid alignment
  t <- t - grid.init(t,dt=dt,W=W)

  # fractional grid index -- starts at >=1
  index <- t/dt
  index <- index - round(index[1])

  while(index[1]<1) { index <- index + 1 }
  while(index[1]>=2) { index <- index - 1 }

  # data is smeared onto p at floor(index) and (1-p) at ceiling(index)
  FLOOR <- floor(index)
  p <- 1-(index-FLOOR)

  # uniform lag grid
  n <- ceiling(last(index)) + 1
  lag <- seq(0,n-1)*dt

  return(list(FLOOR=FLOOR,p=p,lag=lag))
}

############################
# smear data across a uniform grid
gridder <- function(t,z,dt=NULL,W=NULL,lag=NULL,p=NULL,FLOOR=NULL,finish=TRUE,res=1)
{
  n <- length(t)
  COL <- ncol(z)

  # time lags
  DT <- diff(t)

  # default time step
  if(is.null(dt)) { dt <- stats::median(DT) }

  # setup grid transformation
  if(is.null(lag))
  {
    # gap weights to prevent oversampling with coarse dt
    if(is.null(W))
    {
      W <- clamp(DT/dt) # inner weights
      W <- c(1,W) + c(W,1) # left & right weights
      W <- W/2 # average left & right
    }
    else if(!W)
    { W <- rep(1,length(t)) }

    FLOOR <- pregridder(t,dt=dt/res,W=W)
    p <- FLOOR$p
    lag <- FLOOR$lag
    FLOOR <- FLOOR$FLOOR
  }
  else if(!W)
  { W <- rep(1,length(t)) }

  n <- length(lag)
  # continuously distribute times over uniform grid
  W.grid <- numeric(n) # DONT USE ARRAY HERE :(
  Z.grid <- array(0,c(n,COL))
  # lower time spread
  P <- p*W
  W.grid[FLOOR] <- W.grid[FLOOR] + P
  Z.grid[FLOOR,] <- Z.grid[FLOOR,] + P*z
  # upper time spread
  FLOOR <- FLOOR+1
  P <- (1-p)*W
  W.grid[FLOOR] <- W.grid[FLOOR] + P
  Z.grid[FLOOR,] <- Z.grid[FLOOR,] + P*z

  if(finish)
  {
    # normalize distributed information
    POS <- (W.grid>0)
    Z.grid[POS,] <- Z.grid[POS,]/W.grid[POS]

    W <- sum(W) # now total DOF
  }

  return(list(w=W.grid,z=Z.grid,lag=lag,dt=dt))
}

############################
# FFT VARIOGRAM
# SLP sum of lagged product
variogram.fast <- function(data,dt=NULL,res=1,CI="Markov",axes=c("x","y"),ACF=FALSE)
{
  t <- data$t
  z <- get.telemetry(data,axes)
  COL <- ncol(z)

  # smear the data over an evenly spaced time grid (possibly at finer resolution)
  GRID <- gridder(t,z,dt=dt,res=res)
  dt <- GRID$dt
  lag <- GRID$lag
  df <- GRID$w
  # continuous weights eff up the FFT numerics so discretize weights
  w <- sign(GRID$w) # indicator function
  z <- GRID$z
  zz <- z^2

  n <- length(lag)
  #N <- 2*n
  N <- composite(2*n) # keep FFT fast

  df <- FFT(pad(df,N))
  w <- Conj(FFT(pad(w,N)))
  z <- FFT(rpad(z,N))
  zz <- FFT(rpad(zz,N))

  # indicator pair number (integer - stable FFT*)
  DOF <- COL*round(Re(IFFT(abs(w)^2)[1:n]))
  # SVF un-normalized
  if(!ACF) { SVF <- Re(IFFT(Re(w*rowSums(zz))-rowSums(abs(z)^2))[1:n]) }
  else { SVF <- Re(IFFT(rowSums(abs(z)^2))[1:n]) }

  # normalize SVF
  SVF <- SVF/DOF

  # prevent NaNs
  if(any(DOF==0)) { SVF[DOF==0] <- 0 }

  # weight pair number (more accurate DOF)
  DOF <- COL*(Re(IFFT(abs(df)^2)[1:n]))

  # aggregate to each dt in array if discretized at higher resolution
  if(res>1)
  {
    # resolution to aggregate away
    n <- dt/diff(lag[1:2])
    m <- ceiling(length(lag)/n)
    SVF <- pad(SVF,n*(m+1))
    DOF <- pad(DOF,n*(m+1))

    # window weights
    W <- floor(n)
    if(is.even(W)) { W <- W - 1 }
    # full weights
    W <- rep(1,floor(W))
    # partial weights on ends
    if(length(W)<n)
    {
      r <- n-length(W)
      r <- r/2 # half for each end
      W <- c(r,W,r) # half on each end
    }
    SUB <- (length(W)-1)/2
    SUB <- (-SUB):(+SUB)

    lag <- (0:m)*dt # aggregated lags > 0
    SVF <- c(SVF[1], sapply(1:m,function(j){ SUB <- j*n + SUB ; sum(W*DOF[SUB]*SVF[SUB]) }) )
    DOF <- c(DOF[1], sapply(1:m,function(j){ SUB <- j*n + SUB ; sum(W*DOF[SUB]) }) )

    SVF[-1] <- SVF[-1]/DOF[-1]
    if(any(DOF==0)) { SVF[DOF==0] <- 0 }
  }

  # only count non-overlapping lags... not perfect
  if(CI=="Markov")
  {
    # number of lags in the data
    dof <- COL*(last(t)-t[1])/lag
    dof[1] <- COL*length(t)

    # be conservative
    DOF <- pmin(DOF,dof)
  }
  else if(CI=="IID") # fix initial and total DOF
  {
    DOF[1] <- COL*length(t)
    DOF[-1] <- DOF[-1]*(1/sum(DOF[-1])*COL*(length(t)^2-length(t))/2)
  }

  SVF <- data.frame(SVF=SVF,DOF=DOF,lag=lag)

  # contribution to SVF from telemetry error when UERE=1
  #error <- get.error(data,ctmm(axes=axes,error=1))
  #error <- mean(error)
  #result$error <- error
  #result$error[1] <- 0

  return(SVF)
}

##################################
# LAG-WEIGHTED VARIOGRAM
variogram.slow <- function(data,dt=NULL,CI="Markov",axes=c("x","y"),ACF=FALSE)
{
  t <- data$t
  #error <- get.error(data,ctmm(axes=axes,error=1)) # telemetry error when UERE=1
  z <- get.telemetry(data,axes)
  COL <- ncol(z)

  n <- length(t)

  # time lags
  DT <- diff(t)

  # default time step
  if(is.null(dt)) { dt <- stats::median(DT) }

  # interval-weights for each time
  W.T <- clamp(DT/dt) # inner weights
  W.T <- c(1,W.T) + c(W.T,1) # left + right weights
  W.T <- W.T / 2 # average left & right

  # matrices (vectorizing all n^2 operations)
  # matrix of lags
  LAG <- abs(outer(t,t,'-'))
  # matrix of semi-variances (not normalized)
  if(!ACF) { VAR <- lapply(1:COL, function(i) { outer(z[,i],z[,i],'-')^2 }) }
  else { VAR <- lapply(1:COL, function(i) { outer(z[,i],z[,i],'*') }) }
  VAR <- Reduce("+",VAR)
  # matrix of weights
  W.T <- W.T %o% W.T

  # fractional index
  LAG <- 1 + LAG/dt
  # integer indices
  I1 <- floor(LAG)
  I2 <- I1 + 1 # not ceiling if LAG = round(LAG)

  # fraction to deposit at INT
  W1 <- 1 - (LAG-I1)
  # fraction to deposit at INT+1
  W2 <- 1 - W1
  rm(LAG)

  # final weights
  W1 <- W.T * W1
  W2 <- W.T * W2
  rm(W.T)

  # where we will store stuff
  lag <- seq(0,ceiling((t[n]-t[1])/dt)+1)*dt
  SVF <- numeric(length(lag))
  #ERR <- numeric(length(lag))
  DOF <- numeric(length(lag))
  DOF2 <- numeric(length(lag))

  # accumulate
  accumulate <- function(K,W,svf)
  {
    SVF[K] <<- SVF[K] + W*svf
    #ERR[K] <<- ERR[K] + W*err
    DOF[K] <<- DOF[K] + W
    DOF2[K] <<- DOF2[K] + W^2
  }

  pb <- utils::txtProgressBar(style=3)
  for(i in 1:n)
  {
    for(j in i:n)
    {
      accumulate(I1[i,j],W1[i,j],VAR[i,j])
      accumulate(I2[i,j],W2[i,j],VAR[i,j])
    }
    utils::setTxtProgressBar(pb,(i*(2*n-i))/(n^2))
  }
  rm(I1,I2,W1,W2,VAR)

  # delete missing lags
  SVF <- data.frame(SVF=SVF,DOF=DOF,DOF2=DOF2,lag=lag)
  SVF <- SVF[DOF>0,]
  rm(DOF,DOF2,lag)

  # normalize SVF before we correct DOF
  SVF$SVF <- SVF$SVF/((2*COL)*SVF$DOF)
  #error <- error/DOF

  # effective DOF from weights
  SVF$DOF <- pmin(SVF$DOF,SVF$DOF^2/SVF$DOF2)
  SVF$DOF2 <- NULL

  # only count non-overlapping lags... still not perfect
  if(CI=="Markov")
  {
    dof <- sapply(SVF$lag,function(lag){ sum(DT[DT<=lag])/lag })
    dof[1] <- length(t)

    SVF$DOF <- pmin(SVF$DOF,dof)
  }
  else if(CI=="IID") # fix initial and total DOF
  {
    SVF$DOF[1] <- length(t)
    SVF$DOF[-1] <- SVF$DOF[-1]/sum(SVF$DOF[-1])*(length(t)^2-length(t))/2
  }

  # finish off DOF, one for x and y
  SVF$DOF <- COL * SVF$DOF

  close(pb)

  return(SVF)
}

########################
# AVERAGE VARIOGRAMS
mean.variogram <- function(x,...)
{
  x <- c(x,list(...))
  IDS <- length(x)

  # assemble observed lag range
  lag.min <- numeric(IDS) # this will be dt
  lag.max <- numeric(IDS)
  for(id in 1:IDS)
  {
    n <- length(x[[id]]$lag)
    lag.max[id] <- x[[id]]$lag[n]
    lag.min[id] <- x[[id]]$lag[2]
  }
  lag.max <- max(lag.max)
  lag.min <- min(lag.min)
  lag <- seq(0,ceiling(lag.max/lag.min))*lag.min

  # where we will store everything
  n <- length(lag)
  SVF <- numeric(n)
  DOF <- numeric(n)

  # accumulate semivariance
  for(id in 1:IDS)
  {
    for(i in 1:length(x[[id]]$lag))
    {
      # lag index
      j <- 1 + round(x[[id]]$lag[i]/lag.min)
      # number weighted accumulation
      DOF[j] <- DOF[j] + x[[id]]$DOF[i]
      SVF[j] <- SVF[j] + x[[id]]$DOF[i]*x[[id]]$SVF[i]
    }
  }

  # delete missing lags
  variogram <- data.frame(SVF=SVF,DOF=DOF,lag=lag)
  variogram <- subset(variogram,DOF>0)

  # normalize SVF
  variogram$SVF <- variogram$SVF / variogram$DOF

  # drop unused levels
  variogram <- droplevels(variogram)

  # average average errors
  DOF <- sapply(1:length(x),function(i){ x[[i]]$DOF[1] })
  error <- sapply(1:length(x),function(i){ attr(x[[i]],"info")$error })
  error <- sum(DOF * error)/sum(DOF)

  info <- mean.info(x)
  info$error <- error

  variogram <- new.variogram(variogram,info=info)
  return(variogram)
}
#methods::setMethod("mean",signature(x="variogram"), function(x,...) mean.variogram(x,...))


#################
# consolodate info attributes from multiple datasets
mean.info <- function(x)
{
  if(class(x) != "list") { return( attr(x,"info")$identity ) }

  # mean identity
  identity <- sapply(x , function(v) { attr(v,"info")$identity } )
  identity <- unique(identity) # why did I do this?
  identity <- paste(identity,collapse=" ")

  info=attr(x[[1]],"info")
  info$identity <- identity

  return(info)
}


