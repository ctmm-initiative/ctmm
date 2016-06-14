# variogram class
new.variogram <- methods::setClass("variogram",representation("data.frame",info="list"))

# extend subset method
subset.variogram <- function(x,...)
{
  info <- attr(x,"info")
  error <- attr(x,"error")
  x <- subset.data.frame(x,...)
  x < - droplevels(x)
  new.variogram(x,info=info)
}

`[.variogram` <- function(x,...)
{
  info <- attr(x,"info")
  error <- attr(x,"error")
  x <- utils::getS3method("[","data.frame")(x,...)
  if(class(x)=="data.frame") { x <- new.variogram(x,info=info) }
  return(x)
}

#########################
# variogram funcion wrapper
variogram <- function(data,dt=NULL,fast=TRUE,CI="Markov",axes=c("x","y"))
{
  if(length(dt)<2)
  { var <- variogram.dt(data,dt=dt,fast=fast,CI=CI,axes=axes) }
  else
  {
    # calculate a variograms at each dt
    vars <- lapply(dt, function(DT) { variogram.dt(data,dt=DT,fast=fast,CI=CI,axes=axes) } )
    
    # subset each variogram to relevant range of lags
    dt <- c(dt,Inf)
    lag <- vars[[1]]$lag
    vars[[1]] <- vars[[1]][lag<=dt[2],]
    for(i in 1:(length(dt)-1))
    {
      lag <- vars[[i]]$lag
      vars[[i]] <- vars[[i]][(dt[i]<lag)&(lag<=dt[i+1]),]
    }
    
    # coalate
    var <- vars[[1]]
    for(i in 2:(length(dt)-1)) { var <- rbind(var,vars[[i]]) }
  }
  
  # average error when UERE=1
  error <- get.error(data,ctmm(error=TRUE,axes=axes))
  error <- mean(error)
  
  info=attr(data,"info")
  info$axes <- axes
  info$error <- error
  
  var <- new.variogram(var,info=info)
  return(var)
} 
  
################################
# wrapper for fast and slow variogram codes, for a specified dt
variogram.dt <- function(data,dt=NULL,fast=NULL,CI="Markov",axes=c("x","y"))
{
  # intelligently select algorithm
  if(is.null(fast))
  {
    if(length(data$t)<1000) { fast <- FALSE }
    else { fast <- TRUE }
  }
  
  if(fast)
  { SVF <- variogram.fast(data=data,dt=dt,CI=CI,axes=axes) }
  else
  { SVF <- variogram.slow(data=data,dt=dt,CI=CI,axes=axes) }
  
  # skip missing data
  SVF <- SVF[where(SVF$DOF>0),]
  SVF <- stats::na.omit(SVF)
  return(SVF)
}

############################
# best initial time for a uniform grid
grid.init <- function(t,dt=stats::median(diff(t)),W=array(1,length(t)))
{
  # simple analytic periodic cost function
  # COST = sum_t w(t) sin^2(pi(t-t0)/dt)
  # maximized anaytically
  theta <- (2*pi/dt)*t
  SIN <- W %*% sin(theta)
  COS <- W %*% cos(theta)
  t0 <- -dt/(2*pi)*atan(SIN/COS)
  
  return(t0)   
}

############################
# smear data across a uniform grid
gridder <- function(t,z,dt=NULL)
{
  n <- length(t)
  COL <- ncol(z)
  
  # time lags
  DT <- diff(t)
  
  # default time step
  if(is.null(dt)) { dt <- stats::median(DT) }
  
  # gap weights to prevent oversampling with coarse dt
  W <- clamp(c(DT[1],DT)/dt) # left weights
  W <- W + clamp(c(DT,DT[n-1])/dt) # + right weights
  W <- W/2 # average left and right
  
  # choose best grid alignment
  t <- t - grid.init(t,dt,W)
  
  # fractional grid index -- starts at >=1
  index <- t/dt
  while(index[1]<1) { index <- index + 1 }
  while(index[1]>=2) { index <- index - 1 }
  
  # uniform lag grid
  n <- ceiling(max(index))
  lag <- seq(0,n-1)*dt
  
  # continuously distribute times over uniform grid
  W.grid <- numeric(n) # DONT USE ARRAY HERE :(
  Z.grid <- array(0,c(n,COL))
  for(i in 1:length(t))
  {
    j <- index[i]
    
    if(floor(j)==ceiling(j))
    { # trivial case
      J <- round(j)
      w <- W[i] # total weight
      W.grid[J] <- W.grid[J] + w
      Z.grid[J,] <- Z.grid[J,] + w*z[i,]
    }
    else
    { # distribution information between adjacent grids
      # left grid portion
      J <- floor(j)
      w <- W[i]*(1-(j-J))
      W.grid[J] <- W.grid[J] + w
      Z.grid[J,] <- Z.grid[J,] + w*z[i,]

      # right grid portion
      J <- ceiling(j)
      w <- W[i]*(1-(J-j))
      W.grid[J] <- W.grid[J] + w
      Z.grid[J,] <- Z.grid[J,] + w*z[i,]
    }
  }
  
  # normalize distributed information
  POS <- (W.grid>0)
  Z.grid[POS,] <- Z.grid[POS,]/W.grid[POS]
  
  
  # continuous weights eff up the FFT numerics so discretize weights
  W <- sum(W) # now total DOF
  W.grid <- sign(W.grid) # discrete weights
  
  return(list(w=W.grid,z=Z.grid,lag=lag,dt=dt))
}

############################
# FFT VARIOGRAM
# SLP sum of lagged product
variogram.fast <- function(data,dt=NULL,fast=fast,CI="Markov",axes=c("x","y"),SLP=FALSE,ACF=FALSE)
{
  t <- data$t
  z <- get.telemetry(data,axes)
  COL <- ncol(z)

  # smear the data over an evenly spaced time grid
  GRID <- gridder(t,z,dt)
  W.grid <- GRID$w
  Z.grid <- GRID$z
  lag <- GRID$lag
  
  n <- length(lag)
  #N <- 2*n
  N <- composite(2*n)

  W.grid <- Conj(FFT(pad(W.grid,N)))
  ZZ.grid <- FFT(rpad(Z.grid^2,N))
  Z.grid <- FFT(rpad(Z.grid,N))

  # pair number. one for x and y data
  DOF <- COL*round(Re(IFFT(abs(W.grid)^2)[1:n]))
  # SVF un-normalized
  SVF <- Re(IFFT(Re(W.grid*rowSums(ZZ.grid))-rowSums(abs(Z.grid)^2))[1:n])
  if(SLP) { slp <- Re(IFFT(rowSums(abs(Z.grid)^2))[1:n]) }
  
  # delete missing lags
  SVF <- data.frame(SVF=SVF,DOF=DOF,lag=lag)
  if(SLP) { SVF$SLP <- slp }
  SVF <- subset(SVF,DOF>0)
  if(SLP) { slp <- SVF$SLP }
  lag <- SVF$lag
  DOF <- SVF$DOF
  SVF <- SVF$SVF
  
  # normalize SVF
  SVF <- SVF/DOF
  
  # only count non-overlapping lags... not perfect
  if(CI=="Markov")
  {
    dof <- 2*(last(t)-t[1])/lag
    dof[1] <- 2*length(t)
  
    for(i in 1:length(lag))
    {
      if(dof[i]<DOF[i]) {DOF[i] <- dof[i] }
    }
  }
  else if(CI=="IID") # fix initial and total DOF
  {
    DOF[1] <- 2*length(t)
    DOF[-1] <- DOF[-1]/sum(DOF[-1])*(length(t)^2-length(t))
  }
  
  result <- data.frame(SVF=SVF,DOF=DOF,lag=lag)
  if(SLP) { result$SLP <- slp }
  
<<<<<<< HEAD
  # contribution to SVF from telemetry error when UERE=1
  #error <- get.error(data,ctmm(axes=axes,error=1))
  #error <- mean(error)
  #result$error <- error
  #result$error[1] <- 0
  
=======
>>>>>>> refs/remotes/origin/Chris
  return(result)
}

##################################
# LAG-WEIGHTED VARIOGRAM
variogram.slow <- function(data,dt=NULL,CI="Markov",axes=c("x","y"))
{
  t <- data$t
<<<<<<< HEAD
  error <- get.error(data,ctmm(axes=axes,error=1)) # telemetry error when UERE=1
=======
  #error <- get.error(data,ctmm(axes=axes,error=1)) # telemetry error when UERE=1
>>>>>>> refs/remotes/origin/Chris
  z <- get.telemetry(data,axes)
  COL <- ncol(z)
  
  n <- length(t)
  
  # time lags
  DT <- diff(t)
  DT.L <- c(DT[1],DT)
  DT.R <- c(DT,DT[n-1])
  
  # default time step
  if(is.null(dt)) { dt <- stats::median(DT) }

  # where we will store stuff
  lag <- seq(0,ceiling((t[n]-t[1])/dt))*dt
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
      tau <- t[j] - t[i]
      svf <- mean((z[j,]-z[i,])^2)/2
      #err <- (error[i]+error[j])/2 # telemetry error when UERE=1

      # gap weight
      if(tau==0) { w <- 1 }
      else { w <- (clamp(DT.L[j]/tau)+clamp(DT.R[j]/tau))*(clamp(DT.L[i]/tau)+clamp(DT.R[i]/tau)) }
      
      # fractional index
      k <- tau/dt + 1

      if(floor(k)==ceiling(k))
      { # even sampling
        # lag index
        K <- round(k)
        # total weight
        W <- w
        # accumulate        
        accumulate(K,W,svf)
      }
      else
      { # account for drift by distributing semi-variance
        
        # left index
        K <- floor(k)
        # left weight
        W <- w*(1-(k-K))
        # accumulate left portion
        accumulate(K,W,svf)
        
        # right index
        K <- ceiling(k)
        # right weight
        W <- w*(1-(K-k))
        # accumulate right portion
        accumulate(K,W,svf)
      }
    }
    utils::setTxtProgressBar(pb,(i*(2*n-i))/(n^2))
  }
  
  # delete missing lags
  SVF <- data.frame(SVF=SVF,DOF=DOF,DOF2=DOF2,lag=lag)
  SVF <- subset(SVF,DOF>0)
  lag <- SVF$lag
  DOF <- SVF$DOF
  DOF2 <- SVF$DOF2
  #error <- SVF$error
  SVF <- SVF$SVF
  
  # normalize SVF
  SVF <- SVF/DOF
  #error <- error/DOF
  # effective DOF from weights, one for x and y
  DOF <- COL*DOF^2/DOF2
  
  # only count non-overlapping lags... still not perfect
  if(CI=="Markov")
  {
    dof <- COL*length(t)
    if(dof<DOF[1]) { DOF[1] <- dof  }
    
    for(i in 2:length(lag))
    { # large gaps are missing data
      dof <- COL*sum(DT[DT<=lag[i]])/lag[i]
      if(dof<DOF[i]) { DOF[i] <- dof }
      
      utils::setTxtProgressBar(pb,i/length(lag))
    }
  }
  else if(CI=="IID") # fix initial and total DOF
  {
    DOF[1] <- COL*length(t)
    DOF[-1] <- DOF[-1]/sum(DOF[-1])*(length(t)^2-length(t))*COL/2
  }
  
  close(pb)
  
  result <- data.frame(SVF=SVF,DOF=DOF,lag=lag)
  return(result)
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
<<<<<<< HEAD
  error <- sapply(1:length(x),function(i){ attr(x[[i]],"error") })
=======
  error <- sapply(1:length(x),function(i){ attr(x[[i]],"info")$error })
>>>>>>> refs/remotes/origin/Chris
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
  { COV <- diag(0,1+K+(if(circle){1}else{0})) }
  else
  { COV <- area2var(CTMM) }
    
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
  
  if(moment)
  { drift <- get(CTMM$mean) }
  else
  { drift <- stationary }
  MEAN <- drift@svf(CTMM)
  
  SVF <- function(t) { svf(t) + MEAN$EST(t) }
  
  # variance of SVF
  VAR <- function(t)
  {
    g <- grad(t)
    return( (g %*% COV %*% g) + MEAN$VAR(t) )
  }
  
  # chi-square effective degrees of freedom
  DOF <- function(t) { return( 2*SVF(t)^2/VAR(t) ) }
  
  return(list(svf=SVF,VAR=VAR,DOF=DOF,ACF=ACF))
}


##########
plot.svf <- function(lag,CTMM,error=0,alpha=0.05,col="red",type="l",...)
{
  SVF <- svf.func(CTMM,moment=TRUE)
  svf <- SVF$svf
  DOF <- SVF$DOF

  # telemetry error svf
  error <- CTMM$error^2 * error
  errf <- function(t){ if(t==0) {0} else {error} }
    
  # point estimate plot
  SVF <- Vectorize(function(t) { svf(t)+errf(t) })
  graphics::curve(SVF,from=0,to=lag,n=1000,add=TRUE,col=col,type=type,...)
  
  # confidence intervals if COV provided
  if(any(diag(CTMM$COV)>0))
  {
    Lags <- seq(0,lag,lag/1000)
    
    for(j in 1:length(alpha))
    {
      svf.lower <- Vectorize(function(t){ svf(t) * CI.lower(DOF(t),alpha[j]) + errf(t) })
      svf.upper <- Vectorize(function(t){ svf(t) * CI.upper(DOF(t),alpha[j]) + errf(t) })
      
      graphics::polygon(c(Lags,rev(Lags)),c(svf.lower(Lags),rev(svf.upper(Lags))),col=scales::alpha(col,0.1/length(alpha)),border=NA,...)
    }
  }
  
}

###########################################################
# PLOT VARIOGRAM
###########################################################
plot.variogram <- function(x, CTMM=NULL, level=0.95, fraction=0.5, col="black", col.CTMM="red", ...)
{  
  alpha <- 1-level
  
  # number of variograms
  if(class(x)=="variogram" || class(x)=="data.frame") { x <- list(x) }
  n <- length(x)

  # default single comparison model
  if(is.null(CTMM) && n==1 && !is.null(attr(x[[1]],"info")$CTMM)) { CTMM <- attr(x[[1]],"info")$CTMM }
  
  # maximum lag in data
  max.lag <- sapply(x, function(v){ last(v$lag) } )
  max.lag <- max(max.lag)
  # subset fraction of data
  max.lag <- fraction*max.lag
  
  # subset all data to fraction of total period
  x <- lapply(x, function(v) { subset.data.frame(v, lag <= max.lag) })

  # maximum CI on SVF
  max.SVF <- max(sapply(x, function(v){ max(v$SVF * CI.upper(v$DOF,min(alpha))) } ))
  # limit plot range to twice max SVF point estimate (otherwise hard to see)
  max.cap <- 2*max(sapply(x, function(v){ max(v$SVF) } ))
  if(max.SVF>max.cap) { max.SVF <- max.cap }
  
  # choose SVF units
  SVF.scale <- unit(max.SVF,"area")
  SVF.name <- SVF.scale$name
  SVF.scale <- SVF.scale$scale
  
  # choose lag units
  lag.scale <- unit(max.lag,"time",2)
  lag.name <- lag.scale$name
  lag.scale <- lag.scale$scale
  
  xlab <- paste("Time-lag ", "(", lag.name, ")", sep="")
  ylab <- paste("Semi-variance ", "(", SVF.name, ")", sep="")
  
  # fix base plot layer
  plot(c(0,max.lag/lag.scale),c(0,max.SVF/SVF.scale), xlab=xlab, ylab=ylab, col=grDevices::rgb(1,1,1,0), ...)
  
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
      SVF.lower <- SVF * CI.lower(DOF,alpha[j])
      SVF.upper <- SVF * CI.upper(DOF,alpha[j])
      
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
  if(is.null(CTMM$tau) || is.null(CTMM$sigma))
  {
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
    
    # preserve orientation and eccentricity if available/necessary
    if(is.null(CTMM$sigma) || length(CTMM$axes)==1)
    { CTMM$sigma <- sigma }
    else
    {
      CTMM$sigma <- CTMM$sigma@par
      CTMM$sigma[1] <- sigma / cosh(CTMM$sigma[2]/2)
    }
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