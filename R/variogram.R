# variogram class
new.variogram <- methods::setClass("variogram",representation("data.frame",info="list"))

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

  # average error when UERE=1
  error <- get.error(data,ctmm(error=TRUE,axes=axes))
  error <- mean(error)

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
  # mean identity
  identity <- sapply(x , function(v) { attr(v,"info")$identity } )
  identity <- unique(identity) # why did I do this?
  identity <- paste(identity,collapse=" ")

  info=attr(x[[1]],"info")
  info$identity <- identity

  return(info)
}


#########
# update to moment/cumulant with non-stationary mean
svf.func <- function(CTMM,moment=FALSE)
{
  # pull out relevant model parameters
  tau <- CTMM$tau

  # trace variance
  sigma <- mean(diag(CTMM$sigma)) # now AM.sigma

  CPF <- CTMM$CPF
  circle <- CTMM$circle

  range <- CTMM$range
  tau <- tau[tau>0]
  tau <- tau[tau<Inf]
  K <- length(tau)

  # no error considered if missing
  COV <- CTMM$COV
  if(is.null(COV))
  { COV <- diag(0,1+K+(if(circle){1}else{0})+(if(CTMM$error){1}else{0})) }
  else
  { COV <- area2var(CTMM,MEAN=TRUE) }

  # FIRST CONSTRUCT STANDARD ACF AND ITS PARAMTER GRADIENTS
  if(CPF) # Central place foraging
  {
    nu <- 2*pi/tau[1]
    f <- 1/tau[2]
    acf <- function(t){ (cos(nu*t)+f/nu*sin(nu*t))*exp(-f*t) }
    acf.grad <- function(t) { -c(2*pi,1)/tau^2 * c( (-(t+f/nu^2)*sin(nu*t)+f/nu*t*cos(nu*t))*exp(-f*t) , -t*acf(t) + 1/nu*sin(nu*t)*exp(-f*t) ) }
  }
  else if(K==0 && range) # Bivariate Gaussian
  {
    acf <- function(t){ if(t==0) {1} else {0} }
    acf.grad <- function(t){ NULL }
  }
  else if(K==0) # Brownian motion
  {
    acf <- function(t){ 1-t }
    acf.grad <- function(t){ NULL }
  }
  else if(K==1 && range) # OU motion
  {
    acf <- function(t){ exp(-t/tau) }
    acf.grad <- function(t){ t/tau^2*acf(t) }
  }
  else if(K==1) # IOU motion
  {
    acf <- function(t) { 1-(t-tau*(1-exp(-t/tau))) }
    acf.grad <- function(t){ 1-(1+t/tau)*exp(-t/tau) }
  }
  else if(K==2) # OUF motion
  {
    acf <- function(t){ diff(tau*exp(-t/tau))/diff(tau) }
    acf.grad <- function(t) { c(1,-1)*((1+t/tau)*exp(-t/tau)-acf(t))/diff(tau) }
  }

  # finish off svf function including circulation if present
  if(!circle)
  {
    ACF <- function(t) { acf(t) }
    svf <- function(t) { sigma*(1-acf(t)) }
    grad <- function(t) { c(svf(t)/sigma, -sigma*acf.grad(t)) }
  }
  else
  {
    f <- 2*pi/circle
    ACF <- function(t) { cos(f*t)*acf(t) }
    svf <- function(t) { sigma*(1-cos(f*t)*acf(t)) }
    grad <- function(t) { c(svf(t)/sigma, -sigma*cos(f*t)*acf.grad(t), -(f/circle)*sigma*t*sin(f*t)*acf(t)) }
  }

  # add error term
  if(CTMM$error && ("error" %in% dimnames(COV)[[1]]))
  {
    err.svf <- function(t) { (if(t==0) {0} else {CTMM$error^2/2}) }
    GRAD <- function(t) { c(grad(t) , (if(t==0) {0} else {CTMM$error}) ) }
  }
  else
  {
    err.svf <- function(t) { 0 }
    GRAD <- grad
  }

  if(moment)
  { drift <- get(CTMM$mean) }
  else
  { drift <- stationary }
  MEAN <- drift@svf(CTMM)

  SVF <- function(t) { svf(t) + err.svf(t) + MEAN$EST(t) }

  # variance of SVF
  VAR <- function(t)
  {
    g <- GRAD(t)
    return( c(g %*% COV %*% g) + MEAN$VAR(t) )
  }

  # chi-square effective degrees of freedom
  DOF <- function(t) { return( 2*SVF(t)^2/VAR(t) ) }

  return(list(svf=SVF,VAR=VAR,DOF=DOF,ACF=ACF))
}


##########
plot.svf <- function(lag,CTMM,error=0,alpha=0.05,col="red",type="l",...)
{
  # adjust model error to incorporate data HDOP average
  if(CTMM$error)
  {
    if(error==0) { error <- 1 } # default HDOP=1
    error <- sqrt(error)
    CTMM$error <- CTMM$error * error

    if(!is.null(CTMM$COV) && "error" %in% dimnames(CTMM$COV)[1])
    {
      CTMM$COV["error",] <- CTMM$COV["error",] * error
      CTMM$COV[,"error"] <- CTMM$COV[,"error"] * error
    }
  }

  SVF <- svf.func(CTMM,moment=TRUE)
  svf <- SVF$svf
  DOF <- SVF$DOF

  # point estimate plot
  SVF <- Vectorize(function(t) { svf(t) })
  graphics::curve(SVF,from=0,to=lag,n=1000,add=TRUE,col=col,type=type,...)

  # confidence intervals if COV provided
  if(any(diag(CTMM$COV)>0))
  {
    Lags <- seq(0,lag,lag/1000)

    for(j in 1:length(alpha))
    {
      svf.lower <- Vectorize(function(t){ svf(t) * CI.lower(DOF(t),alpha[j]) })
      svf.upper <- Vectorize(function(t){ svf(t) * CI.upper(DOF(t),alpha[j]) })

      graphics::polygon(c(Lags,rev(Lags)),c(svf.lower(Lags),rev(svf.upper(Lags))),col=scales::alpha(col,0.1/length(alpha)),border=NA,...)
    }
  }

}

###########################################################
# PLOT VARIOGRAM
###########################################################
plot.variogram <- function(x, CTMM=NULL, level=0.95, fraction=0.5, col="black", col.CTMM="red", xlim=NULL, ylim=NULL, ...)
{
  alpha <- 1-level

  # number of variograms
  if(class(x)=="variogram" || class(x)=="data.frame") { x <- list(x) }
  n <- length(x)

  # default single comparison model
  # if(is.null(CTMM) && n==1 && !is.null(attr(x[[1]],"info")$CTMM)) { CTMM <- attr(x[[1]],"info")$CTMM }
  ACF <- !is.null(attr(x[[1]],"info")$ACF)

  # maximum lag in data
  max.lag <- sapply(x, function(v){ last(v$lag) } )
  max.lag <- max(max.lag)
  # subset fraction of data
  max.lag <- fraction*max.lag

  # subset all data to fraction of total period
  x <- lapply(x, function(v) { subset.data.frame(v, lag <= max.lag) })

  if(!ACF) # SVF plot
  {
    min.SVF <- 0
    # maximum CI on SVF
    max.SVF <- max(sapply(x, function(v){ max(v$SVF * CI.upper(v$DOF,min(alpha))) } ))
    # limit plot range to twice max SVF point estimate (otherwise hard to see)
    max.cap <- 2*max(sapply(x, function(v){ max(v$SVF) } ))
    max.SVF <- min(max.SVF,max.cap)

    # choose SVF units
    SVF.scale <- unit(max.SVF,"area")
    SVF.name <- SVF.scale$name
    SVF.scale <- SVF.scale$scale

    SVF.name <- c(SVF.name,unit(max.SVF,"area",concise=TRUE)$name)
    SVF.name[3] <- SVF.name[2]

    ylab <- "Semi-variance"
    ylab <- c(ylab,ylab,"SVF")

    # range of possible ylabs with decreasing size
    ylab <- paste(ylab, " (", SVF.name, ")", sep="")

  }
  else # ACF plot
  {
    min.SVF <- sapply(x,function(v) { min(v$SVF) })
    min.SVF <- min(0,1.1*min.SVF)

    max.SVF <- sapply(x,function(v) { max(v$SVF[v$lag>0]) })
    max.SVF <- max(0,1.1*max.SVF)

    SVF.scale <- 1
    ylab <- "Autocorrelation"
    ylab <- c(ylab,ylab,"ACF")
  }

  # choose lag units
  lag.scale <- unit(max.lag,"time",2)
  lag.name <- lag.scale$name
  lag.scale <- lag.scale$scale

  lag.name <- c(lag.name,unit(max.lag,"time",thresh=2,concise=TRUE)$name)
  lag.name[3] <- lag.name[2]

  xlab <- "Time-lag"
  xlab <- c(xlab,xlab,"Lag")

  xlab <- paste(xlab, " (", lag.name, ")", sep="")

  # choose appropriately sized axis labels for base plot
  lab <- rbind(xlab,ylab)

  #work out string width max
  max.cex.w <- lab # copy dimensions and preserve below
  max.cex.w[] <- par('pin')/strwidth(lab,'inches')
  #work out string height max
  max.cex.h <- lab
  max.cex.h[] <- (par('mai')[1:2]/par('mar')[1:2])/strheight(lab,'inches')

  # min of x & y
  max.cex.w <- pmin(max.cex.w[1,],max.cex.w[2,])
  max.cex.h <- pmin(max.cex.h[1,],max.cex.h[2,])
  # min of width and height
  max.cex <- pmin(max.cex.w,max.cex.h)

  lab <- 1
  if(max.cex[lab]<1) { lab <- lab + 1 }
  if(max.cex[lab]<1) { lab <- lab + 1 }

  # unit convert scales if supplied
  if(!is.null(xlim)) { xlim <- xlim/lag.scale }
  if(!is.null(ylim)) { ylim <- ylim/SVF.scale }

  # fix base plot layer
  plot(c(0,max.lag/lag.scale),c(min.SVF/SVF.scale,max.SVF/SVF.scale), xlim=xlim, ylim=ylim, xlab=xlab[lab], ylab=ylab[lab], col=grDevices::rgb(1,1,1,0), ...)

  # color array for plots
  col <- array(col,n)

  for(i in 1:n)
  {
    lag <- x[[i]]$lag/lag.scale
    SVF <- x[[i]]$SVF/SVF.scale
    DOF <- x[[i]]$DOF

    # make sure plot looks nice and appropriate for data resolution
    type <- "l"
    if(length(lag) < 100) { type <- "p" }

    graphics::points(lag, SVF, type=type, col=col[[i]])

    for(j in 1:length(alpha))
    {
      # chi-square CIs for semi-variance
      if(!ACF)
      {
        SVF.lower <- SVF * CI.lower(DOF,alpha[j])
        SVF.upper <- SVF * CI.upper(DOF,alpha[j])
      }
      else # Fisher CIs for autocorrelation
      {
        # subset relevant data for Fisher transformation
        STUFF <- data.frame(lag=lag,SVF=SVF,DOF=DOF)
        STUFF <- STUFF[DOF>3,]
        lag <- STUFF$lag
        SVF <- STUFF$SVF
        DOF <- STUFF$DOF
        # Fisher transformation
        FISH <- atanh(SVF)
        SD <- 1/sqrt(DOF-3)
        SVF.lower <- tanh(stats::qnorm(alpha[j]/2,mean=FISH,sd=SD,lower.tail=TRUE))
        SVF.upper <- tanh(stats::qnorm(alpha[j]/2,mean=FISH,sd=SD,lower.tail=FALSE))

        graphics::abline(h=c(-1,0,1)/sqrt(DOF[1])*stats::qnorm(1-alpha[j]/2),col="red",lty=c(2,1,2))
      }

      graphics::polygon(c(lag,rev(lag)),c(SVF.lower,rev(SVF.upper)),col=scales::alpha(col[[i]],alpha=0.1),border=NA)
    }
  }

  # NOW PLOT THE MODELS
  if(!is.null(CTMM))
  {
    if(class(CTMM)=="ctmm") { CTMM <- list(CTMM) }
    n <- length(CTMM)

    # color array for plots
    col <- array(col.CTMM,n)
    type <- "l"

    for(i in 1:n)
    {
      # units conversion
      CTMM[[i]] <- unit.ctmm(CTMM[[i]],sqrt(SVF.scale),lag.scale)

      # include errors in svf plot
      error <- FALSE
      if(CTMM[[i]]$error)
      {
        j <- FALSE
        # choose relevant data
        if(length(x)==1) { j <- 1 }
        else if(length(x)==n) { j <- i }
        else warning("Ambiguity in who's error to plot.")

        if(j) { error <- attr(x[[j]],"info")$error }
      }

      plot.svf(max.lag/lag.scale,CTMM[[i]],error=error,alpha=alpha,type=type,col=col[[i]])
    }
  }

}
# PLOT.VARIOGRAM METHODS
#methods::setMethod("plot",signature(x="variogram",y="missing"), function(x,y,...) plot.variogram(x,...))
#methods::setMethod("plot",signature(x="variogram",y="variogram"), function(x,y,...) plot.variogram(list(x,y),...))
#methods::setMethod("plot",signature(x="variogram",y="ctmm"), function(x,y,...) plot.variogram(x,model=y,...))
#methods::setMethod("plot",signature(x="variogram"), function(x,...) plot.variogram(x,...))


#######################################
# plot a variogram with zoom slider
#######################################
zoom.variogram <- function(x, fraction=0.5, ...)
{
  # R CHECK CRAN BUG WORKAROUND
  z <- NULL

  # number of variograms
  n <- 1
  if(class(x)=="list") { n <- length(x) }
  else {x <- list(x) } # promote to list of one

  # maximum lag in data
  max.lag <- sapply(x, function(v){ last(v$lag) } )
  max.lag <- max(max.lag)

  min.lag <- sapply(x, function(v){ v$lag[2] } )
  min.lag <- min(min.lag)

  b <- 4
  min.step <- min(fraction,10*min.lag/max.lag)
  manipulate::manipulate( { plot.variogram(x, fraction=b^(z-1), ...) }, z=manipulate::slider(1+log(min.step,b),1,initial=1+log(fraction,b),label="zoom") )
}
methods::setMethod("zoom",signature(x="variogram"), function(x,fraction=0.5,...) zoom.variogram(x,fraction=fraction,...))


####################################
# guess variogram model parameters #
####################################
variogram.guess <- function(variogram,CTMM=ctmm())
{
  # guess at some missing parameters
  n <- length(variogram$lag)

  # variance estimate
  sigma <- mean(variogram$SVF[2:n])

  # peak curvature estimate
  # should occur at short lags
  v2 <- 2*max((variogram$SVF/variogram$lag^2)[2:n])

  # free frequency
  Omega2 <- v2/sigma
  Omega <- sqrt(Omega2)

  # peak diffusion rate estimate
  # should occur at intermediate lags
  # diffusion parameters
  D <- (variogram$SVF/variogram$lag)[2:n]
  # index of max diffusion
  tauD <- which.max(D)
  # max diffusion
  D <- D[tauD]
  # time lag of max diffusion
  tauD <- variogram$lag[tauD]

  # average f-rate
  f <- -log(D/(sigma*Omega))/tauD

  CPF <- CTMM$CPF
  if(CPF) # frequency, rate esitmate
  {
    omega2 <- Omega2 - f^2
    if(f>0 && omega2>0)
    { tau <- c(2*pi/sqrt(omega2),1/f) }
    else # bad backup estimate
    {
      tau <- sqrt(2)/Omega
      tau <- c(2*pi*tau , tau)
    }
  }
  else # position, velocity timescale estimate
  { tau <- c(sigma/D,D/v2)}

  if(!CTMM$range) { sigma <- D ; tau[1] <- Inf }

  if(length(CTMM$tau)==0) { CTMM$tau <- tau }
  #else if(length(CTMM$tau)==1) { CTMM$tau[2] <- tau[2] }

  # preserve orientation and eccentricity if available/necessary
  if(is.null(CTMM$sigma) || length(CTMM$axes)==1)
  { CTMM$sigma <- sigma }
  else
  {
    CTMM$sigma <- CTMM$sigma@par
    CTMM$sigma[1] <- sigma / cosh(CTMM$sigma[2]/2)
  }

  # don't overwrite or lose ctmm parameters not considered here
  model <- as.list(CTMM) # getDataPart deletes names()
  model$info <- attr(variogram,"info")
  model <- do.call("ctmm",model)
  return(model)
}


######################################################################
# visual fit of the variogram
######################################################################
variogram.fit <- function(variogram,CTMM=ctmm(),name="GUESS",fraction=0.5,interactive=TRUE,...)
{
  if(interactive && !manipulate::isAvailable()) { interactive <- FALSE }
  envir <- .GlobalEnv

  # R CHECK CRAN BUG WORKAROUNDS
  z <- NULL
  tau1 <- 1
  tau2 <- 0
  error <- CTMM$error
  store <- NULL ; rm(store)

  m <- 2 # slider length relative to point guestimate
  n <- length(variogram$lag)

  # fill in missing parameters non-destructively
  CTMM <- variogram.guess(variogram,CTMM)
  if(!interactive) { return(CTMM) }

  # parameters for logarithmic slider
  b <- 4
  min.step <- 10*variogram$lag[2]/variogram$lag[n]
  #min.step <- max(min.step,fraction)

  # manipulation controls
  manlist <- list(z = manipulate::slider(1+log(min.step,b),1,initial=1+log(fraction,b),label="zoom"))

  K <- length(CTMM$tau)
  if(K==1) { CTMM$tau[2] <- 0 }

  range <- CTMM$range
  sigma <- mean(diag(CTMM$sigma))
  if(range)
  {
    sigma.unit <- unit(sigma,"area",concise=TRUE)
    sigma <- sigma / sigma.unit$scale
    label <- paste("sigma variance (",sigma.unit$name,")",sep="")
    manlist <- c(manlist, list(sigma = manipulate::slider(0,m*sigma,initial=sigma,label=label)))
  }
  else
  {
    sigma.unit <- unit(sigma,"diffusion",concise=TRUE)
    sigma <- sigma / sigma.unit$scale
    label <- paste("sigma diffusion (",sigma.unit$name,")",sep="")
    manlist <- c(manlist, list(sigma = manipulate::slider(0,m*sigma,initial=sigma,label=label)))
  }

  CPF <- CTMM$CPF
  tau <- CTMM$tau
  tau1.unit <- unit(tau[1],"time",2,concise=TRUE)
  tau2.unit <- unit(tau[2],"time",2,concise=TRUE)
  tau[1] <- tau[1] / tau1.unit$scale
  tau[2] <- tau[2] / tau2.unit$scale
  if(CPF)
  {
    label <- paste("tau period (",tau1.unit$name,")",sep="")
    manlist <- c(manlist, list(tau1 = manipulate::slider(0,m*tau[1],initial=tau[1],label=label)))

    label <- paste("tau decay (",tau2.unit$name,")",sep="")
    manlist <- c(manlist, list(tau2 = manipulate::slider(0,m*tau[2],initial=tau[2],label=label)))

    tau2 <- NULL # not sure why necessary
  }
  else
  {
    label <- paste("tau position (",tau1.unit$name,")",sep="")
    manlist <- c(manlist, list(tau1 = manipulate::slider(0,m*tau[1],initial=tau[1],label=label)))

    label <- paste("tau velocity (",tau2.unit$name,")",sep="")
    manlist <- c(manlist, list(tau2 = manipulate::slider(0,m*tau[2],initial=tau[2],label=label)))
  }

  # circulation
  circle <- CTMM$circle
  if(circle)
  {
    circle.unit <- unit(circle,"time",concise=TRUE)
    circle <- circle / circle.unit$scale
    label <- paste("circulation (",circle.unit$name,")",sep="")
    c1 <- min(0,m*circle)
    c2 <- max(0,m*circle)
    manlist <- c(manlist, list(circle = manipulate::slider(c1,c2,initial=circle,label=label)))
  }

  # error
  e2 <- max(100,2*error)
  manlist <- c(manlist, list(error = manipulate::slider(0,e2,initial=as.numeric(error),step=0.1,label="error (m)")))

  # storage button
  manlist <- c(manlist, list(store = manipulate::button(paste("Save to",name))))

  if(!range)
  {
    manlist$tau1 <- NULL
    tau1 <- Inf
  }

  if(K==1)
  {
    manlist$tau2 <- NULL
    tau2 <- 0
  }

  # non-destructive parameter overwrite
  manipulate::manipulate(
    {
      # store trace, but preserve angle & eccentricity
      if(length(CTMM$axes)==2)
      {
        CTMM$sigma <- CTMM$sigma@par
        CTMM$sigma[1] <- sigma * sigma.unit$scale / cosh(CTMM$sigma[2]/2)
      }
      else
      { CTMM$sigma <- sigma }

      CTMM$tau <- c(tau1*tau1.unit$scale, tau2*tau2.unit$scale)
      if(circle) { CTMM$circle <- circle * circle.unit$scale }
      CTMM$error <- error

      CTMM <- as.list(CTMM)
      CTMM$info <- attr(variogram,"info")
      CTMM <- do.call("ctmm",CTMM)
      fraction <- b^(z-1)
      if(store) { assign(name,CTMM,envir=envir) }
      plot.variogram(variogram,CTMM=CTMM,fraction=fraction,...)
    },
    manlist
  )
}
