intensity <- function(data,UD,RSF,R=list(),variable=NULL,empirical=TRUE,theoretical=TRUE,level=0.95,ticks=TRUE,smooth=TRUE,interpolate=TRUE,xlim=NULL,ylim=NULL,...)
{
  # how to sample rasters
  interpolate <- rep(interpolate,length(R))
  interpolate <- ifelse(interpolate,"bilinear","simple")
  names(interpolate) <- names(R)

  axes <- RSF$axes
  GEO <- c('longitude','latitude')
  CTMM <- UD@CTMM
  if(smooth && any(CTMM$error>0)) { data[,c(axes,GEO)] <- predict(data,CTMM=CTMM,t=data$t,complete=TRUE)[,c(axes,GEO)] }

  RVARS <- names(R)
  formula <- RSF$formula

  if(is.null(variable))
  {
    if(length(RVARS)==1)
    { variable <- RVARS[1] }
    else
    {
      for(r in RVARS)
      {
        intensity(data,UD=UD,RSF=RSF,R=R,variable=r,level=level,smooth=smooth,interpolate=interpolate,...)
        return()
      }
    }
  }

  ### USED #############################
  axes <- variable
  w <- UD$weights
  w <- w / sum(w)
  error <- 0.001
  res <- 10
  n <- UD$DOF.H[1]
  w2d <- 1/n
  w2o <- (n-1)/n
  MISE <- function(h) { w2d/sqrt(2*h^2) + w2o/sqrt(2+2*h^2) - 2/sqrt(2+h^2) + 1/sqrt(2) }
  # silverman's rule of thumb
  DIM <- 1
  h <- (4/(DIM+2)/n)^(1/(DIM+4))
  # find 1D bandwidth for same n
  control <- list(precision=.Machine$double.eps^(1/4))
  h <- optimizer(h,MISE,lower=0,control=control)$par
  # find bias for same h
  bias <- 1 - 1/n + h^2
  # integral K^2
  RK2 <- 1/sqrt(4*pi)

  # annotate data with variable of interest
  VARIABLE <- R[[variable]]
  PROJ <- raster::projection(VARIABLE)
  xy <- get.telemetry(data,GEO)
  xy <- project(xy,to=PROJ)
  data[,variable] <- raster::extract(VARIABLE,xy,method=interpolate[variable])

  # RFIT <- ctmm.guess(data,ctmm(axes=axes),interactive=FALSE)
  # RFIT <- ctmm.select(data,RFIT)
  RFIT <- ctmm.fit(data,ctmm(axes=axes))
  H <- h^2 * methods::getDataPart(RFIT$sigma)
  EXT <- extent(RFIT,level=1-error)[,axes,drop=FALSE] # Gaussian extent (includes uncertainty)
  dr <- c(sqrt(H)/res)
  grid <- format_grid(NULL,axes=axes) # grid not defined
  grid <- kde.grid(data,H=H,axes=axes,alpha=error,res=res,dr=dr,grid=grid,EXT.min=EXT)
  KDE <- kde(data,H=H,axes=axes,CTMM=RFIT,bias=bias,W=w,alpha=error,dr=dr,grid=grid,...)
  rm(grid)

  log.PDF <- log(KDE$PDF)
  VAR.log.PDF <- c(RK2/n/sqrt(H)) / KDE$PDF
  alpha <- 1-level
  z <- qnorm(1-alpha/2)
  SE.log.PDF <- z * sqrt(VAR.log.PDF)

  ### AVAILABLE ########################
  # zero out variable of interest for partial effect plot
  raster::values(R[[variable]]) <- 0 # interaction terms remain - variable could be an interaction
  if(projection(data)==raster::projection(VARIABLE))
  { grid <- VARIABLE }
  else
  { grid <- NULL }
  AGDE <- agde(data,RSF,R=R,grid=grid)
  rm(R,grid)

  # evaluate raster on AGDE grid
  VARIABLE <- R.grid(AGDE$r,projection(AGDE),VARIABLE)

  # extract variable availability
  R <- KDE$r[[axes]] # resource axis
  dR <- KDE$dr
  R1 <- R[1]
  P <- array(0,length(R)) # probability mass axis

  # flatten arrays
  AGDE <- c(AGDE$PDF)
  VARIABLE <- c(VARIABLE)

  # get indices
  VARIABLE <- (VARIABLE-R1)/dR
  FLOOR <- floor(VARIABLE)
  w <- 1-(VARIABLE-FLOOR)
  SUB <- which(FLOOR>=1 & FLOOR<=length(R))
  for(i in SUB)
  {
    I <- FLOOR[i]
    P[I] <- P[I] + w[i]*AGDE[i]
  }
  rm(FLOOR)

  CEIL <- ceiling(VARIABLE)
  w <- CEIL-VARIABLE
  SUB <- which(CEIL>=1 & CEIL<=length(R))
  for(i in SUB)
  {
    I <- CEIL[i]
    P[I] <- P[I] + w[i]*AGDE[i]
  }
  rm(CEIL,VARIABLE)

  # normalize the same as used
  P <- (sum(KDE$PDF)/sum(P))*P
  # this also fixes the density units (after log difference)
  log.P <- log(P)

  # Propagate Available uncertainty ?

  # combine everything
  SUB <- KDE$PDF>.Machine$double.eps & P>.Machine$double.eps
  R <- R[SUB]
  log.PDF <- log.PDF[SUB]
  SE.log.PDF <- SE.log.PDF[SUB]
  log.P <- log.P[SUB]
  log.UA <- log.PDF - log.P

  if(is.null(xlim)) { xlim <- range(R) }
  ylab <- '\u2206log(\u03BB)'
  if(is.null(ylim))
  {
    ylim <- c(min(log.UA),max(log.UA))
    dy <- ylim[2]-ylim[1]
    ylim <- c(max(ylim[1]-1.5*dy,min(log.UA-SE.log.PDF)) , min(ylim[2]+1.5*dy,max(log.UA+SE.log.PDF)))
  }
  if(empirical)
  {
    graphics::plot(R,log.UA,col='black',type='l',lwd=2,xlim=xlim,ylim=ylim,xlab=variable,ylab=ylab,...)
    graphics::polygon(c(R,rev(R)),c(log.UA-SE.log.PDF,rev(log.UA+SE.log.PDF)),border=NA,col=malpha('black',0.25))
  }

  if(theoretical)
  {
    if(empirical)
    { plot <- graphics::points }
    else
    { plot <- graphics::plot }

    # choose best intercept for data
    I <- which.min(SE.log.PDF)

    EST <- RSF$beta[variable]*(R-R[I])+log.UA[I]
    plot(R,EST,col='red',type='l',lwd=2,xlab=variable,ylab=ylab,xlim=xlim,ylim=ylim,...)
    SE <- sqrt(RSF$COV[variable,variable])*abs(R-R[I])
    graphics::polygon(c(R,rev(R)),c(EST-SE,rev(EST+SE)),border=NA,col=malpha('red',0.25))
  }

  # ticks at top corresponding to sampled resource values?
  if(ticks) { graphics::axis(side=3,at=data[[variable]],labels=FALSE,col=malpha('blue',0.5)) }
}
