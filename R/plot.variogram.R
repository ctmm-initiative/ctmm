#########
# error argument is <VAR> (fixed) or <VAR>/UERE^2 (fitted) passed from SVF object
# CTMM$error is logical (fixed) or UERE (fitted)
svf.func <- function(CTMM,moment=FALSE,error=0)
{
  # # adjust model error to incorporate data HDOP average
  # if(CTMM$error)
  # {
  #   error <- sqrt(length(CTMM$axes)*error) # <HDOP>
  #   # UERE adjustment
  #   CTMM$error <- CTMM$error * error # UERE * <HDOP> stored in CTMM$error
  #   if(!is.null(CTMM$COV) && "error" %in% dimnames(CTMM$COV)[1])
  #   {
  #     CTMM$COV["error",] <- CTMM$COV["error",] * error
  #     CTMM$COV[,"error"] <- CTMM$COV[,"error"] * error
  #   }
  # }

  # pull out relevant model parameters
  tau <- CTMM$tau

  # trace variance
  sigma <- mean(diag(CTMM$sigma)) # now AM.sigma

  circle <- CTMM$circle

  # no error considered if missing
  COV <- CTMM$COV
  if(!is.null(COV)) { COV <- area2var(CTMM,MEAN=TRUE) }

  range <- CTMM$range
  tau <- tau[tau<Inf]
  if(any(tau==0))
  {
    DEL <- paste("tau",names(tau[tau==0]))
    if(!is.null(COV)) { COV <- rm.name(COV,DEL) }
    tau <- tau[tau>0]
  }
  K <- length(tau)

  # FIRST CONSTRUCT STANDARD ACF AND ITS PARAMTER GRADIENTS
  if(K==0 && range) # Bivariate Gaussian
  {
    NAMES <- NULL
    acf <- function(t){ if(t==0) {1} else {0} }
    acf.grad <- function(t){ NULL }
  }
  else if(K==0) # Brownian motion
  {
    NAMES <- NULL
    acf <- function(t){ 1-t }
    acf.grad <- function(t){ NULL }
  }
  else if(K==1 && range) # OU motion
  {
    NAMES <- "tau position"
    acf <- function(t){ exp(-t/tau) }
    acf.grad <- function(t){ t/tau^2*acf(t) }
  }
  else if(K==1) # IOU motion
  {
    NAMES <- "tau velocity"
    acf <- function(t) { 1-(t-tau*(1-exp(-t/tau))) }
    acf.grad <- function(t){ 1-(1+t/tau)*exp(-t/tau) }
  }
  else if(K==2) # OUF motion
  {
    NAMES <- c("tau position","tau velocity")
    acf <- function(t){ diff(tau*exp(-t/tau))/diff(tau) }
    acf.grad <- function(t) { c(1,-1)*((1+t/tau)*exp(-t/tau)-acf(t))/diff(tau) }
  }

  # !!! this is all hard-coded dumb to the ordering of COV !!!

  # finish off svf function including circulation if present
  if(!circle)
  {
    ACF <- function(t) { acf(t) }
    svf <- function(t) { sigma*(1-acf(t)) }
    grad <- function(t) { c(svf(t)/sigma, -sigma*acf.grad(t)) }
  }
  else
  {
    NAMES <- c(NAMES,"circle")
    ACF <- function(t) { cos(circle*t)*acf(t) }
    svf <- function(t) { sigma*(1-cos(circle*t)*acf(t)) }
    grad <- function(t) { c(svf(t)/sigma, -sigma*cos(circle*t)*acf.grad(t), +sigma*t*sin(circle*t)*acf(t)) }
  }
  NAMES <- c("variance",NAMES)

  # add error term
  if(CTMM$error)
  {
    err.svf <- function(t) { (if(t==0) {0} else {CTMM$error^2 * error}) }
    if("error" %in% dimnames(CTMM$COV)[[1]]) # fit or fixed error?
    {
      NAMES <- c(NAMES,"error")
      GRAD <- function(t) { c(grad(t) , (if(t==0) {0} else {2 * CTMM$error * error}) ) }
    }
    else
    { GRAD <- grad }
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

  # no error provided
  if(is.null(COV)) { COV <- diag(0,nrow=length(GRAD(0))) }

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
# error argument is average variance of data's error passed from SVF object
plot.svf <- function(lag,CTMM,error=0,alpha=0.05,col="red",type="l",...)
{
  SVF <- svf.func(CTMM,moment=TRUE,error=error)
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
  x <- listify(x)
  n <- length(x)

  # default single comparison model
  # if(is.null(CTMM) && n==1 && !is.null(attr(x[[1]],"info")$CTMM)) { CTMM <- attr(x[[1]],"info")$CTMM }
  ACF <- !is.null(attr(x[[1]],"info")$ACF)

  # subset the data if xlim or fraction provided
  if(!is.null(xlim))
  {
    fraction <- 1 # xlim overrides fraction
    x <- lapply(x,function(y){ y[xlim[1]<=y$lag & y$lag<=xlim[2],] })
  }
  else
  {
    # maximum lag in data
    max.lag <- sapply(x, function(v){ last(v$lag) } )
    max.lag <- max(max.lag,xlim)
    # subset fraction of data
    max.lag <- fraction*max.lag

    # subset all data to fraction of total period
    if(fraction<1) { x <- lapply(x, function(y) { y[y$lag<=max.lag,] }) }

    xlim <- c(0,max.lag)
  }

  # calculate ylimits from all variograms
  if(is.null(ylim)) { ylim <- extent(x,level=max(level))$y }

  if(!ACF) # SVF plot
  {
    # choose SVF units
    SVF.scale <- unit(ylim,"area")
    SVF.name <- SVF.scale$name
    SVF.scale <- SVF.scale$scale

    SVF.name <- c(SVF.name,unit(ylim,"area",concise=TRUE)$name)
    SVF.name[3] <- SVF.name[2]

    ylab <- "Semi-variance"
    ylab <- c(ylab,ylab,"SVF")

    # range of possible ylabs with decreasing size
    ylab <- paste(ylab, " (", SVF.name, ")", sep="")

  }
  else # ACF plot
  {
    SVF.scale <- 1
    ylab <- "Autocorrelation"
    ylab <- c(ylab,ylab,"ACF")
  }

  # choose lag units
  lag.scale <- unit(xlim,"time",2)
  lag.name <- lag.scale$name
  lag.scale <- lag.scale$scale

  lag.name <- c(lag.name,unit(xlim,"time",thresh=2,concise=TRUE)$name)
  lag.name[3] <- lag.name[2]

  xlab <- "Time-lag"
  xlab <- c(xlab,xlab,"Lag")

  xlab <- paste(xlab, " (", lag.name, ")", sep="")

  # choose appropriately sized axis labels for base plot
  lab <- rbind(xlab,ylab)

  # string width max
  max.cex.w <- lab # copy dimensions and preserve below
  max.cex.w[] <- graphics::par('pin')/graphics::strwidth(lab,'inches')
  # string height max
  max.cex.h <- lab
  max.cex.h[] <- (graphics::par('mai')[1:2]/graphics::par('mar')[1:2])/graphics::strheight(lab,'inches')

  # min of x & y
  max.cex.w <- pmin(max.cex.w[1,],max.cex.w[2,])
  max.cex.h <- pmin(max.cex.h[1,],max.cex.h[2,])
  # min of width and height
  max.cex <- pmin(max.cex.w,max.cex.h)

  lab <- 1
  if(max.cex[lab]<1) { lab <- lab + 1 }
  if(max.cex[lab]<1) { lab <- lab + 1 }

  # unit convert scales if supplied
  xlim <- xlim/lag.scale
  ylim <- ylim/SVF.scale

  # fix base plot layer
  plot(xlim,ylim, xlim=xlim, ylim=ylim, xlab=xlab[lab], ylab=ylab[lab], col=grDevices::rgb(1,1,1,0), ...)

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
    CTMM <- listify(CTMM)
    n <- length(CTMM)

    # color array for plots
    col <- array(col.CTMM,n)
    type <- "l"

    for(i in 1:n)
    {
      # units conversion
      CTMM[[i]] <- unit.ctmm(CTMM[[i]],length=sqrt(SVF.scale),time=lag.scale)

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

      plot.svf(xlim[2],CTMM[[i]],error=error,alpha=alpha,type=type,col=col[[i]])
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
  x <- listify(x)
  n <- length(x)

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


