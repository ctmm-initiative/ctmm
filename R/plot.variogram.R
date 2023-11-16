#########
# CTMM$error is logical (fixed) or UERE (fitted)
# error (modulo UERE) is now an argument of the output functions in addition to lag
svf.func <- function(CTMM,moment=FALSE)
{
  CTMM <- get.taus(CTMM) # pre-calculate stuff
  # pull out relevant model parameters
  tau <- CTMM$tau

  # trace variance
  sigma <- var.covm(CTMM$sigma,ave=TRUE)

  circle <- CTMM$circle

  # no error considered if missing
  COV <- CTMM[["COV"]]
  if(!is.null(COV)) { COV <- axes2var(CTMM,MEAN=TRUE) }

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
  NAMES <- CTMM$tau.names
  if(K==0 && range) # Bivariate Gaussian
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
    if(tau[1]>tau[2]) # overdamped
    {
      acf <- function(t){ diff(tau*exp(-t/tau))/diff(tau) }
      acf.grad <- function(t) { c(1,-1)*((1+t/tau)*exp(-t/tau)-acf(t))/diff(tau) }
    }
    else if(!CTMM$omega) # critically damped
    {
      tau <- tau[1]
      acf <- function(t){ (1+t/tau)*exp(-t/tau) }
      acf.grad <- function(t) { (t^2/tau^3)*exp(-t/tau) }
    }
    else if(CTMM$omega) # underdamped
    {
      f <- CTMM$f.nu[1]
      nu <- CTMM$f.nu[2]
      acf <- function(t) { (cos(nu*t) + (f/nu)*sin(nu*t))*exp(-f*t) }
      acf.grad <- function(t) { c( exp(-f*t) * c( -t*cos(nu*t) + (1-f*t)/nu*sin(nu*t) , -(t+f/nu^2)*sin(nu*t) + (f/nu)*t*cos(nu*t) ) %*% CTMM$J.nu.tau ) }
    }
  }

  # finish off svf function including circulation if present
  if(!circle)
  {
    ACF <- function(t) { acf(t) }
    svf <- function(t) { sigma*(1-acf(t)) }
    grad <- function(t,...) { c(svf(t)/sigma, -sigma*acf.grad(t)) }
  }
  else
  {
    NAMES <- c(NAMES,"circle")
    ACF <- function(t) { cos(circle*t)*acf(t) }
    svf <- function(t) { sigma*(1-cos(circle*t)*acf(t)) }
    grad <- function(t,...) { c(svf(t)/sigma, -sigma*cos(circle*t)*acf.grad(t), +sigma*t*sin(circle*t)*acf(t)) }
  }
  NAMES <- c("variance",NAMES)

  # # add error terms?
  # if(any(CTMM$error))
  # {
  #   err.svf <- function(t,error=0) { ifelse(t>0,c(error %*% CTMM$error^2),0) }
  #   PARS <- dimnames(CTMM$COV)[[1]]
  #   PARS <- PARS[grepl("error",PARS)]
  #   if(length(PARS)) # fit or fixed error?
  #   {
  #     NAMES <- c(NAMES,PARS)
  #     PARS <- substr(PARS,nchar("error ?"),nchar(PARS))
  #     PARS <- which(names(CTMM$error) %in% PARS)
  #     GRAD <- function(t,error=0) { c(grad(t) , ifelse(rep(t>0,length(PARS)),2*c(error[PARS] * CTMM$error[PARS]),0) ) }
  #   }
  #   else
  #   { GRAD <- grad }
  # }
  # else
  # {
  #   err.svf <- function(t,error=0) { 0 }
  #   GRAD <- grad
  # }

  if(moment)
  { MEAN <- drift.svf(CTMM) }
  else
  { MEAN <- stationary.svf(CTMM) }

#  SVF <- function(t,error=0) { svf(t) + err.svf(t,error=error) + MEAN$EST(t) }
  SVF <- function(t) { svf(t) + MEAN$EST(t) }

  # no error provided
  if(is.null(COV)) { COV <- diag(0,nrow=length(grad(0))) }

  # empty covariance matrix
  BLANK <- array(0,length(NAMES)*c(1,1))
  dimnames(BLANK) <- list(NAMES,NAMES)

  # information that we have from CTMM
  NAMES <- NAMES[NAMES %in% dimnames(COV)[[1]]]
  # store that information appropriately
  if(length(NAMES)) { BLANK[NAMES,NAMES] <- COV[NAMES,NAMES] }
  # copy over
  BLANK -> COV

  # variance of SVF
  VAR <- function(t)
  {
    # g <- GRAD(t,error=error)
    g <- grad(t)
    return( c(g %*% COV %*% g) + MEAN$VAR(t) )
  }

  # chi-square effective degrees of freedom
  # DOF <- function(t,error=0) { return( 2*SVF(t,error=error)^2/VAR(t,error=error) ) }
  DOF <- function(t) { return( 2*SVF(t)^2/VAR(t) ) }

  return(list(svf=SVF,VAR=VAR,DOF=DOF,ACF=ACF))
}


##########
plot_svf <- function(lag,CTMM,alpha=0.05,col="red",type="l",...)
{
  # changed from max lag to all lags
  # changed from error=0 or number/logical to error=NULL or array

  # number of pixels across diagonal of display
  PX <- ceiling(sqrt(sum(grDevices::dev.size("px")^2)))

  # are we including errors?
  # ERROR <- !is.null(error) && CTMM$error
  # e0 <- 0
  # can we plot a smooth curve?
  # if(!ERROR)
  {
    lag <- seq(0,last(lag),length.out=PX)
    # error <- 0 -> e0
  }
  # else if(all(diff(error[-1])==0))
  # { error <- error[2] -> e0 } # can still plot curve because all errors it the same

  SVF <- svf.func(CTMM,moment=TRUE)
  svf <- SVF$svf
  DOF <- SVF$DOF

  # point estimate plot
  # SVF <- Vectorize(function(t,error=e0) { svf(t,error=error) })
  SVF <- Vectorize(function(t) { svf(t) })

  lag[1] <- lag[2]/1000 # almost go to origin, but not quite to avoid nugget

  # if(length(error)==1) # can plot curve
  { graphics::curve(SVF,from=0,to=last(lag),n=PX,add=TRUE,col=col,...) }
  # else
  # { graphics::points(lag,SVF(lag,error),col=col,type=type,...) }

  # confidence intervals if COV provided
  if(any(diag(CTMM$COV)>0))
  {
    # SVF <- Vectorize(function(t,error=e0){ svf(t,error=error) })(lag,error)
    SVF <- Vectorize(function(t){ svf(t) })(lag)

    for(j in 1:length(alpha))
    {
      # dof <- Vectorize(function(t,error=e0) { DOF(t,error=error) })(lag,error)
      dof <- Vectorize(function(t) { DOF(t) })(lag)
      svf.lower <- Vectorize(function(df){ CI.lower(df,alpha[j]) })(dof)
      svf.upper <- Vectorize(function(df){ CI.upper(df,alpha[j]) })(dof)

      graphics::polygon(c(lag,rev(lag)),c(SVF*svf.lower,rev(SVF*svf.upper)),col=malpha(col,0.1/length(alpha)),border=NA,...)
    }
  }

}

###########################################################
# PLOT VARIOGRAM
###########################################################
plot.variogram <- function(x, CTMM=NULL, level=0.95, units=TRUE, fraction=0.5, col="black", col.CTMM="red", xlim=NULL, ylim=NULL, ext=NULL, ...)
{
  alpha <- 1-level

  if(!is.null(ext))
  {
    xlim <- ext$x
    ylim <- ext$y
  }

  # empirical variograms
  x <- listify(x)
  n <- length(x)
  # theoretical models
  CTMM <- listify(CTMM)
  m <- length(CTMM)

  # default single comparison model
  # if(is.null(CTMM) && n==1 && !is.null(attr(x[[1]],"info")$CTMM)) { CTMM <- attr(x[[1]],"info")$CTMM }
  ACF <- attr(x[[1]],"info")$ACF
  ACF <- !is.null(ACF) && ACF

  RESIDUAL <- attr(x[[1]],"info")$residual
  RESIDUAL <- !is.null(RESIDUAL) && RESIDUAL

  ULAG <- attr(x[[1]],"info")$lags # NULL or 'time' by default
  if(is.null(ULAG)) { ULAG <- "time" }

  axes <- attr(x[[1]],"info")$axes
  TYPE <- DOP.match(axes)
  UNITS <- DOP.LIST[[TYPE]]$units # NA if unknown

  # don't plot if DOF<1
  x <- lapply(x,function(y){ y[y$DOF>=1,] })

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

  # clamp ACF estimates before calculating extent
  if(ACF) { x <- lapply(x,function(y) { y$SVF <- clamp(y$SVF,-1,1) ; y }) }

  # calculate ylimits from all variograms !!!
  if(is.null(ylim)) { ylim <- extent(x,level=max(level))$y }

  if(ACF) # ACF plot
  {
    ylab <- "Autocorrelation"
    ylab <- c(ylab,ylab,"ACF")
  }
  else # SVF plot
  {
    ylab <- "Semi-variance"
    ylab <- c(ylab,ylab,"SVF")
  }

  if(RESIDUAL || is.na(UNITS)) # unitless
  { SVF.scale <- 1 }
  else if(UNITS=="distance") # units of original data (axes)
  {
    # choose SVF units
    SVF.scale <- unit(ylim,"area",SI=!units)
    SVF.name <- SVF.scale$name
    SVF.scale <- SVF.scale$scale

    SVF.name <- c(SVF.name,unit(ylim,"area",concise=TRUE,SI=!units)$name)
    SVF.name[3] <- SVF.name[2]

    # range of possible ylabs with decreasing size
    ylab <- paste(ylab, " (", SVF.name, ")", sep="")
  }
  else if(UNITS=="frequency")
  {
    SVF.scale <- unit(sqrt(ylim),UNITS,SI=!units)
    SVF.name <- paste0(SVF.scale$name," squared") # e.g., per year squared
    SVF.scale <- SVF.scale$scale^2

    CONCISE <- unit(ylim,UNITS,concise=TRUE,SI=!units)$name
    substr(CONCISE,nchar(CONCISE),nchar(CONCISE)) <- "\u00B2"

    SVF.name <- c(SVF.name,CONCISE)
    SVF.name[3] <- SVF.name[2]

    # range of possible ylabs with decreasing size
    ylab <- paste(ylab, " (", SVF.name, ")", sep="")
  }
  else # some generic unit
  {
    SVF.scale <- unit(sqrt(ylim),UNITS,SI=!units)
    SVF.name <- paste0("square ",SVF.scale$name) # e.g., square kilograms
    SVF.scale <- SVF.scale$scale^2

    SVF.name <- c(SVF.name,paste0(unit(ylim,UNITS,concise=TRUE,SI=!units)$name,"\u00B2"))
    SVF.name[3] <- SVF.name[2]

    # range of possible ylabs with decreasing size
    ylab <- paste(ylab, " (", SVF.name, ")", sep="")
  }

  # choose lag units
  if(ULAG=="time")
  {
    lag.scale <- unit(xlim,"time",2,SI=!units)
    lag.name <- lag.scale$name
    lag.scale <- lag.scale$scale

    lag.name <- c(lag.name,unit(xlim,"time",thresh=2,concise=TRUE,SI=!units)$name)
    lag.name[3] <- lag.name[2]

    xlab <- "Time-lag"
    xlab <- c(xlab,xlab,"Lag")

    xlab <- paste(xlab, " (", lag.name, ")", sep="")
  }
  else # unitless lags
  {
    lag.scale <- 1
    xlab <- "Lags"
  }

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
  # color array for plots
  col.CTMM <- array(col.CTMM,m)

  # units conversions
  x <- lapply(x,function(X){unit.variogram(X,time=lag.scale,area=SVF.scale)})
  CTMM <- lapply(CTMM,function(M){unit.ctmm(M,length=sqrt(SVF.scale),time=lag.scale)})

  for(i in 1:n)
  {
    SVF <- x[[i]]$SVF
    # if(!ACF)
    # {
    #   if("MSE" %in% names(x[[i]])) # calibrated errors
    #   { MSE <- x[[i]]$MSE }
    #   else # uncalibrated errors - needs fitted error in model
    #   { MSE <- x[[i]]$MSDOP }
    # }
    lag <- x[[i]]$lag
    DOF <- x[[i]]$DOF

    # make sure plot looks nice and appropriate for data resolution
    if(length(lag) < 100) { type <- "p" } else { type <- "l" }

    graphics::points(lag, SVF, type=type, col=col[[i]],...)

    for(j in 1:length(alpha))
    {
      # chi-square CIs for semi-variance
      if(!ACF)
      {
        LOW <- CI.lower(DOF,alpha[j])
        HIGH <- CI.upper(DOF,alpha[j])

        SVF.lower <- SVF * LOW
        SVF.upper <- SVF * HIGH

        if(RESIDUAL) { graphics::abline(h=1,col="red",...) }
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
        SD <- 1/sqrt(max(DOF-3,0))
        SVF.lower <- tanh(stats::qnorm(alpha[j]/2,mean=FISH,sd=SD,lower.tail=TRUE))
        SVF.upper <- tanh(stats::qnorm(alpha[j]/2,mean=FISH,sd=SD,lower.tail=FALSE))

        if(RESIDUAL) { graphics::abline(h=c(-1,0,1)/sqrt(DOF[1])*stats::qnorm(1-alpha[j]/2),col="red",lty=c(2,1,2),...) }
      }

      graphics::polygon(c(lag,rev(lag)),c(SVF.lower,rev(SVF.upper)),col=malpha(col[[i]],alpha=0.1),border=NA,...)
    }

    # PLOT CORRESPONDING MODEL
    # if(i<=m) { plot_svf(lag,CTMM[[i]],error=MSE,alpha=alpha,type=type,col=col.CTMM[[i]]) }
    if(i<=m) { plot_svf(lag,CTMM[[i]],alpha=alpha,type=type,col=col.CTMM[[i]],...) }
  }
  # PLOT LEFTOVER MODELS USING THE LAST DATA
  # if(n<m) { for(i in n:m) { plot_svf(lag,CTMM[[i]],error=MSE,alpha=alpha,type=type,col=col.CTMM[[i]]) } }
  if(n<m) { for(i in n:m) { plot_svf(lag,CTMM[[i]],alpha=alpha,type=type,col=col.CTMM[[i]],...) } }

  # no projection for variograms
  assign("projection",NULL,pos=plot.env)
  # dimensional type
  assign("x.dim","time",pos=plot.env)
  assign("y.dim","area",pos=plot.env)
  # unit conversion
  assign("x.scale",lag.scale,pos=plot.env)
  assign("y.scale",SVF.scale,pos=plot.env)
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

  # pass second argument non-fraction
  if(class(fraction)[1]!="numeric")
  {
    arg2 <- fraction
    fraction <- 0.5
  }
  else
  {
    arg2 <- NULL
  }

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

  s1 <- 1 + log(min.step,b)
  initial <- 1 + log(fraction,b)

  if(is.null(arg2))
  { manipulate::manipulate( { plot.variogram(x,fraction=b^(z-1),...) }, z=manipulate::slider(s1,1,initial=initial,label="zoom") ) }
  else
  { manipulate::manipulate( { plot.variogram(x,arg2,fraction=b^(z-1),...) }, z=manipulate::slider(s1,1,initial=initial,label="zoom") ) }
}
methods::setMethod("zoom",signature(x="variogram"), function(x,fraction=0.5,...) zoom.variogram(x,fraction=fraction,...))


