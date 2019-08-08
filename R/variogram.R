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
variogram <- function(data,dt=NULL,fast=TRUE,res=1,CI="Markov",axes=c("x","y"),precision=1/8)
{
  res <- round(res)

  if(length(dt)>1)
  {
    dt <- sort(dt,method="quick")
    res <- rep(res,length(dt))
    res[-1] <- pmax(res[-1],dt[-1]/dt[-length(dt)])

    # calculate a variogram at each dt
    SVF <- lapply(1:length(dt), function(i) { variogram.dt(data,dt=dt[i],fast=fast,res=res[i],CI=CI,axes=axes,precision=precision) } )

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
  { SVF <- variogram.dt(data,dt=dt,res=res,fast=fast,CI=CI,axes=axes,precision=precision) }

  # delete missing lags
  SVF <- SVF[SVF$DOF>0,]

  info=attr(data,"info")
  info$axes <- axes
  # info$error <- error

  SVF <- new.variogram(SVF,info=info)
  return(SVF)
}


################################
# wrapper for fast and slow variogram codes, for a specified dt
variogram.dt <- function(data,dt=NULL,fast=NULL,res=1,CI="Markov",axes=c("x","y"),precision=1/2)
{
  # intelligently select algorithm
  if(is.null(fast))
  {
    if(length(data$t)<1000) { fast <- FALSE }
    else { fast <- TRUE }
  }

  # telemetry error when UERE=1
  error <- get.error(data,ctmm(axes=axes,error=TRUE),circle=TRUE)

  if(fast)
  { SVF <- variogram.fast(data=data,error=error,dt=dt,res=res,CI=CI,axes=axes) }
  else
  { SVF <- variogram.slow(data=data,error=error,dt=dt,CI=CI,axes=axes,precision=precision,calibrated=is.calibrated(data)) }

  # skip missing data
  SVF <- SVF[which(SVF$DOF>0),]
  SVF <- stats::na.omit(SVF)

  SVF$SVF[1] <- 0 # what about dupes?
  SVF$MSE[1] <- 0 # what about dupes?

  # clarify that this is not calibrated error
  if(attr(error,'flag')<=2) { SVF <- rename(SVF,"MSE","MSDOP") }

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
  if(!is.na(dt) && dt==0) { dt <- stats::median(DT[DT>0]) }
  if(is.na(dt)) { dt <- 1 } # doesn't really matter

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
variogram.fast <- function(data,error=NULL,dt=NULL,res=1,CI="Markov",axes=c("x","y"),ACF=FALSE)
{
  t <- data$t
  z <- get.telemetry(data,axes)
  COL <- ncol(z)

  # smear the data over an evenly spaced time grid (possibly at finer resolution)
  if(!ACF) { z <- cbind(error,z) }
  GRID <- gridder(t,z,dt=dt,res=res)
  dt <- GRID$dt
  lag <- GRID$lag
  df <- GRID$w
  # continuous weights eff up the FFT numerics so discretize weights
  w <- sign(GRID$w) # indicator function
  if(!ACF) { error <- GRID$z[,1] ; GRID$z <- GRID$z[,-1] }
  z <- GRID$z
  zz <- z^2

  n <- length(lag)
  #N <- 2*n
  N <- composite(2*n) # keep FFT fast

  df <- FFT(pad(df,N))
  w <- Conj(FFT(pad(w,N)))
  z <- FFT(rpad(z,N))
  if(!ACF) { error <- FFT(pad(error,N)) }
  zz <- FFT(rpad(zz,N))

  # indicator pair number (integer - stable FFT*)
  DOF <- round(Re(IFFT(abs(w)^2)[1:n]))
  # SVF un-normalized
  if(!ACF)
  {
    SVF <- Re(IFFT(Re(w*rowSums(zz))-rowSums(abs(z)^2))[1:n])
    error <- Re(IFFT(Re(w*error))[1:n])
  }
  else
  { SVF <- Re(IFFT(rowSums(abs(z)^2))[1:n]) }

  # normalize SVF
  if(!ACF) { error <- error/DOF } # already averaged over dimensions
  DOF <- COL*DOF
  SVF <- SVF/DOF

  # prevent NaNs
  ZERO <- DOF==0
  if(any(ZERO))
  {
    SVF[ZERO] <- 0
    if(!ACF) { error[ZERO] <- 0 }
  }
  rm(ZERO)

  # weight pair number (more accurate DOF)
  DOF <- COL*(Re(IFFT(abs(df)^2)[1:n]))

  rm(df,w,z,zz)

  # aggregate to each dt in array if discretized at higher resolution
  if(res>1)
  {
    # resolution to aggregate away
    n <- dt/diff(lag[1:2])
    m <- ceiling(length(lag)/n)
    SVF <- pad(SVF,n*(m+1))
    if(!ACF) { error <- pad(error,n*(m+1)) }
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
    if(!ACF) { error <- c(error[1], sapply(1:m,function(j){ SUB <- j*n + SUB ; sum(W*DOF[SUB]*error[SUB]) }) ) }
    DOF <- c(DOF[1], sapply(1:m,function(j){ SUB <- j*n + SUB ; sum(W*DOF[SUB]) }) )

    SVF[-1] <- SVF[-1]/DOF[-1]
    if(!ACF)
    { error[-1] <- error[-1]/DOF[-1] }

    ZERO <- DOF<=0
    if(any(ZERO))
    {
      SVF[ZERO] <- 0
      if(!ACF) { error[ZERO] <- 0 }
    }
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
  if(!ACF) { SVF$MSE <- error }

  return(SVF)
}

##################################
# LAG-WEIGHTED VARIOGRAM
variogram.slow <- function(data,error=NULL,dt=NULL,CI="Markov",axes=c("x","y"),ACF=FALSE,precision=1/2,calibrated=FALSE)
{
  t <- data$t
  z <- get.telemetry(data,axes)
  COL <- ncol(z)
  n <- length(t)

  # initial variance guess (for error accounting)
  MVAR <- mean( apply(z,2,stats::var) )

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
  if(!ACF)
  {
    VAR <- lapply(1:COL, function(i) { outer(z[,i],z[,i],'-')^2 })
    EVAR <- outer(error,error,'+')/2
  }
  else
  { VAR <- lapply(1:COL, function(i) { outer(z[,i],z[,i],'*') }) }
  VAR <- Reduce("+",VAR) / ifelse(!ACF,2*COL,COL) # semi-variance or variance
  # matrix of weights
  W.T <- W.T %o% W.T
  rm(z,error)

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
  SVF <- rep(MVAR,length(lag)) # IID initial guess for error accounting
  SVF[1] <- 0 # IID initial guess for error accounting
  if(!ACF & !calibrated) { ERR <- numeric(length(lag)) }
  SVF -> SVF.OLD

  # accumulate semi-variance
  accumulate <- function(K,W,svf,err=0,werr=1)
  {
    if(!ACF)
    {
      if(!calibrated) { ERR[K] <<- ERR[K] + W*err }
      else if(K>1) { werr <- 1/(1+err/SVF.OLD[K]) }
    }

    Werr <- W*werr
    Werr2 <- Werr*werr # W*werr^2

    SVF[K] <<- SVF[K] + Werr2*svf
    DOF[K] <<- DOF[K] + Werr # normalization constant & linear/trend weight
    DOF2[K] <<- DOF2[K] + Werr^2 # sum of squares of linear/trend weight
  }

  # iterative accounting of error
  ERROR <- Inf
  ERROR.OLD <- Inf
  TARGET <- .Machine$double.eps^precision
  pb <- utils::txtProgressBar(style=3) # time loops
  while(ERROR>TARGET && ERROR<=ERROR.OLD)
  {
    PG <- (TARGET/ERROR)^(1/4) # base progress
    SVF <- numeric(length(lag))
    DOF <- numeric(length(lag))
    DOF2 <- numeric(length(lag))

    for(i in 1:n)
    {
      for(j in i:n)
      {
        accumulate(I1[i,j],W1[i,j],VAR[i,j],EVAR[i,j])
        accumulate(I2[i,j],W2[i,j],VAR[i,j],EVAR[i,j])
      }
      utils::setTxtProgressBar(pb,PG+(1-PG)*(i*(2*n-i))/(n^2))
    }

    # normalize SVF before we compare to old and/or correct DOF
    SVF <- SVF/DOF
    if(!ACF && !calibrated) #error <- error/DOF
    { ERR <- ERR/DOF } # error already averaged over dimension
    else if(!ACF && calibrated)
    {
      ERROR <- abs(SVF-SVF.OLD)/pmax(SVF,min(1,MVAR))
      #plot(ERROR)
      ERROR <- ERROR[1:ceiling(length(ERROR)*precision)]
      ERROR <- max(ERROR,na.rm=TRUE)
      #title(ERROR)
      SVF -> SVF.OLD
      ERROR -> ERROR.OLD
    }
    else
    { ERROR <- 0 }

    if(ACF || !calibrated) { break }
  }
  rm(I1,I2,W1,W2,VAR,SVF.OLD)
  if(!ACF) { rm(EVAR) }

  # only count non-overlapping lags... still not perfect
  if(CI=="Markov")
  {
    DT <- sort(DT) # diff
    CDT <- cumsum(DT) # total lag for diff

    # for every lag, what is the max DT<=lag
    # DOF = CDT[max DT]/lag
    j <- length(DT)
    dof <- numeric(length(lag))
    for(i in length(lag):1)
    {
      # go down in DT until lag can fit <=DT
      while(DT[j]>lag[i] && j>1) { j <- j-1 }
      # quit if lag can no longer fit <=DT
      if(j==1 && DT[j]>lag[i]) { break }
      # store result
      dof[i] <- CDT[j]/lag[i]
    }

    dof[1] <- length(t)
    #DOF <- pmin(DOF,dof)
  }

  # delete missing lags
  SVF <- data.frame(SVF=SVF,DOF=DOF,DOF2=DOF2,lag=lag)
  if(!ACF)
  {
    if(calibrated) { SVF$MSE <- 0 }
    else { SVF$MSE <- ERR ; rm(ERR) }
  }
  if(CI=="Markov") { SVF$dof <- dof ; rm(dof) }
  SVF <- SVF[DOF>0,]
  rm(DOF,DOF2,lag)

  # effective DOF from weights, non-overlapping lags, take the most conservative of the 3 estimates
  SVF$DOF <- pmin(SVF$DOF,SVF$DOF^2/SVF$DOF2) ; SVF$DOF2 <- NULL
  if(CI=="Markov") { SVF$DOF <- pmin(SVF$DOF,SVF$dof) ; SVF$dof <- NULL }

  if(CI=="IID") # fix initial and total DOF
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
  MSE <- numeric(n)

  # error type to pull out -- assumes all are the same
  if("MSE" %in% names(x[[1]])) { ERR <- "MSE" } else { ERR <- "MSDOP" }

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
      if(ERR %in% names(x[[id]])) { MSE[j] <- MSE[j] + x[[id]]$DOF[i]*x[[id]][[ERR]][i] } # not ideal fix
    }
  }

  # delete missing lags
  variogram <- data.frame(lag=lag,SVF=SVF,DOF=DOF)
  if(ERR %in% names(x[[1]])) { variogram$MSE <- MSE }
  variogram <- variogram[DOF>0,]

  # normalize SVF
  variogram$SVF <- variogram$SVF / variogram$DOF
  if(ERR %in% names(x[[1]])) { variogram$MSE <- variogram$MSE / variogram$DOF }

  # drop unused levels
  variogram <- droplevels(variogram)

  # correct name if not calibrated
  if(ERR %in% names(x[[1]])) { rename(variogram,"MSE",ERR) }

  variogram <- new.variogram(variogram,info=mean.info(x))
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


