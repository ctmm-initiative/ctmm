# how far you can go before R's bessel function breaks down
BESSEL_LIMIT <- 2^16

# estimate and assign speeds to times
outlie <- function(data,UERE=10,standardize=FALSE,plot=TRUE,...)
{
  # turn of extreme value stuff (too sensitive)
  mspeed=NULL;
  n=NULL

  if(class(data)[1]=="list")
  {
    if(is.null(n))
    {
      n <- sapply(data,nrow)
      n <- sum(n) - length(data)
    }
    return( lapply(1:length(data), function(i){OUT <- outlie(data[[i]],UERE=UERE,mspeed=mspeed,n=n,standardize=standardize,plot=plot); graphics::title(names(data)[i]); OUT} ) )
  }

  error <- get.error(data,ctmm(error=UERE,axes=c("x","y")),circle=TRUE) # VAR
  if(!is.null(mspeed)) { COV <- get.error(data,ctmm(error=UERE,axes=c("x","y")),DIM=2) } # COV
  else { COV <- NULL }

  Vs <- assign_speeds(data,UERE=error,COV=COV,mspeed=mspeed,n=n)
  v <- Vs$v.t
  VAR.v <- Vs$VAR.t
  p.EV <- Vs$p.t

  mu <- median.telemetry(data)

  d <- get.telemetry(data,axes=c("x","y"))
  mu <- get.telemetry(mu,axes=c("x","y"))
  mu <- c(mu)

  # detrend median
  d <- t(d) - mu
  # distances from median
  d <- colSums(d^2)
  d <- sqrt(d)
  D <- distanceMLE(d,error)
  VAR.d <- error/(2-(d^2-D^2)/error)
  d <- D

  if(plot)
  {
    # bounding box
    new.plot(data,...)
    # convert units on telemetry object to match base plot
    data <- unit.telemetry(data,length=get("x.scale",pos=plot.env))

    lwd <- Vs$v.dt
    lwd <- if(diff(range(lwd))) { lwd/max(lwd) } else { 0 }
    col <- grDevices::rgb(0,0,lwd,lwd)
    lwd <- 2*lwd
    n <- length(data$t)
    graphics::segments(x0=data$x[-n],y0=data$y[-n],x1=data$x[-1],y1=data$y[-1],col=col,lwd=lwd,asp=1,...)

    cex <- if(diff(range(d))) { d/max(d) } else { 0 }
    col <- grDevices::rgb(cex,0,0,cex)
    graphics::points(data$x,data$y,col=col,cex=cex,pch=20,...)
  }

  if(standardize)
  {
    MAD.v <- stats::mad(v)
    v <- v/MAD.v
    VAR.v <- VAR.v/MAD.v^2

    MAD.d <- stats::mad(d)
    d <- d/MAD.d
    VAR.d <- VAR.d/MAD.d^2
  }

  VAR.v <- nant(VAR.v,0)
  VAR.d <- nant(VAR.d,0)

  R <- data.frame(t=data$t,speed=v,distance=d,VAR.speed=VAR.v,VAR.distance=VAR.d)
  if(length(p.EV)) { R$p.EV <- p.EV }
  R <- new.outlie(R)
  return(R)
}


#########################
# plot outlier information with error bars
plot.outlie <- function(x,level=0.95,units=TRUE,axes=c('d','v'),...)
{
  n <- nrow(x)
  t <- x$t
  v <- x$speed
  d <- x$distance

  # CIs
  v <- sapply(1:n,function(i){tnorm.hdr(v[i],x$VAR.speed[i],level=level)})
  d <- sapply(1:n,function(i){tnorm.hdr(d[i],x$VAR.distance[i],level=level)})

  # unit conversions ... (simpler)
  UNITS <- unit(t,dimension='time',concise=FALSE,SI=!units)
  t <- t/UNITS$scale
  t.lab <- paste0("Time (",UNITS$name,")")

  UNITS <- unit(v,dimension='length',concise=TRUE,SI=!units)
  v <- v/UNITS$scale
  v.lab <- paste0("Minimum speed (",UNITS$name,"/s)")

  UNITS <- unit(d,dimension='length',concise=TRUE,SI=!units)
  d <- d/UNITS$scale
  d.lab <- paste0("Core deviation (",UNITS$name,")")

  mid <- function(x)
  {
    if(length(dim(x))) { x <- x[2,] }
    return(x)
  }

  x <- get(axes[1])
  y <- get(axes[2])

  xlab <- get(paste0(axes[1],".lab"))
  ylab <- get(paste0(axes[2],".lab"))

  # base plot
  plot(mid(x),mid(y),xlim=range(x),ylim=range(y),pch=19,xlab=xlab,ylab=ylab,...)

  # ERROR BAR PLOT
  # hack: we draw arrows but with very special "arrowheads"
  # (c) Laryx Decidua
  # very annoying warning when zero length --- cannot be suppressed with suppressWarnings()
  if(length(dim(y)))
  {
    SUB <- y[3,]-y[1,] > .Machine$double.eps # still does not avoid annoying warning
    suppressWarnings( graphics::arrows(mid(x)[SUB],y[1,SUB],mid(x)[SUB],y[3,SUB],length=0.05,angle=90,code=3,...) )
  }
  if(length(dim(x)))
  {
    SUB <- x[3,]-x[1,] > .Machine$double.eps # still does not avoid annoying warning
    suppressWarnings( graphics::arrows(x[1,SUB],mid(y)[SUB],x[3,SUB],mid(y)[SUB],length=0.05,angle=90,code=3,...) )
  }
  # will need to switch from arrows to segment to just avoid annoying warning...
}


# speeds assigned by blame
# dt[1] is recording interval
# dt[2] is minimum time between fixes, which can be smaller than dt[1]
# UERE is misnamed and is actually the per-time error covariance
assign_speeds <- function(data,dt=NULL,UERE=0,COV=NULL,mspeed=NULL,n=NULL,method="max")
{
  method <- match.arg(method,c("max","min"))
  p <- p.t <- NULL

  DT <- diff(data$t)
  # speed threshold is specified
  if(!is.null(mspeed))
  {
    if(is.null(n)) { n <- nrow(data) - 1 }

    if(class(mspeed)=="numeric")
    { SVF <- function(tau) { (1/pi*mspeed^2*tau^2) * diag(2) } }
    else if(class(mspeed)=="ctmm")
    {
      sigma <- mspeed$sigma # location covariance
      mspeed$sigma <- covm(1) # unit variance
      SVF <- function(tau) { svf.func(mspeed)$svf(tau) * sigma } # unit covariance SVF * covariance
    }
    else
    { stop("Argument 'mspeed' class not recognized.") }
  }
  else
  { SVF <- NULL }

  if(is.null(dt)) { dt <- time_res(DT) }

  # inner speed estimates
  v.dt <- speedMLE(data,dt=dt,UERE=UERE,DT=DT,SVF=SVF,COV=COV,n=n)
  VAR.dt <- v.dt$VAR
  p.dt <- v.dt$p
  v.dt <- v.dt$X
  if(length(v.dt)==1)
  {
    v <- c(v.dt,v.dt)
    if(length(p.dt)==1) { p <- c(p.dt,p.dt) }
    VAR <- c(VAR.dt,VAR.dt)
    return(list(v.t=v,VAR.t=VAR,v.dt=v.dt,VAR.dt=VAR.dt,p.t=p))
  }

  if(length(UERE)==1)
  {
    # end point contingency estimates
    v1 <- speedMLE(data[c(1,3),],dt=dt,UERE=UERE,SVF=SVF,COV=COV[c(1,3),,],n=n)
    N <- length(data$t)
    v2 <- speedMLE(data[c(N-2,N),],dt=dt,UERE=UERE,SVF=SVF,COV=COV[c(N-2,N),,],n=n)
  }
  else # pull out correct errors for calculation if fed all errors
  {
    v1 <- speedMLE(data[c(1,3),],dt=dt,UERE=UERE[c(1,3)],SVF=SVF,COV=COV[c(1,3),,],n=n)
    N <- length(data$t)
    v2 <- speedMLE(data[c(N-2,N),],dt=dt,UERE=UERE[c(N-2,N)],SVF=SVF,COV=COV[c(N-2,N),,],n=n)
  }
  VAR1 <- v1$VAR; p1 <- v1$p; v1 <- v1$X
  VAR2 <- v2$VAR; p2 <- v2$p; v2 <- v2$X

  # TODO PICK UP HERE

  # left and right estimates - n estimates, 1 lag apart
  v1 <- c(v1,v.dt)
  v2 <- c(v.dt,v2)

  VAR1 <- c(VAR1,VAR.dt)
  VAR2 <- c(VAR.dt,VAR2)

  p1 <- c(p1,p.dt)
  p2 <- c(p.dt,p2)

  if(method=="max")
  {
    # n-1 estimates, 2 lags apart
    v1 <- v1[-N]
    v2 <- v2[-1]

    # which side of lag index looks smaller
    LESS <- v1 < v2

    vs <- sort(v.dt,index.return=TRUE,decreasing=TRUE)
    is <- vs$ix
    vs <- vs$x

    v <- numeric(length(LESS))
    VAR <- numeric(length(LESS))
    # assign blame for smallest speeds, from greatest to least, last/least taking precedence in R
    LESS <- LESS[is]
    v[is+!LESS] <- vs # v.dt[is]
    VAR[is+!LESS] <- VAR.dt[is]

    # assign blame for largest speeds, from least to greatest, last/greatest taking precedence in R
    LESS <- rev(LESS)
    is <- rev(is)
    vs <- rev(vs)
    v[is+LESS] <- vs # v.dt[is]
    VAR[is+LESS] <- VAR.dt[is]

    rm(is,vs,LESS)

    if(length(p.dt))
    {
      # which side of lag index looks more common
      MORE <- p1 > p2

      p1 <- p1[-N]
      p2 <- p2[-1]

      ps <- sort(p.dt,index.return=TRUE,decreasing=FALSE)
      ip <- ps$ix
      ps <- ps$x

      p <- numeric(length(MORE))
      # assign blame for largest p's, from least to most, last/most taking precedence in R
      MORE <- MORE[ip]
      p[ip+!MORE] <- ps # p.dt[ip]

      # assign blame for smallest p's, from least to most, last/least taking precedence in R
      MORE <- rev(MORE)
      ip <- rev(ip)
      ps <- rev(ps)
      p[ip+MORE] <- ps # p.dt[ip]
    }
  }
  else if(method=="min")
  {
    # v <- pmin(v1,v2)

    v <- cbind(v1,v2)
    is <- apply(v,1,which.min)
    v <- vapply(1:nrow(v),function(i){v[i,is[i]]},v[,1])
    VAR <- vapply(1:nrow(VAR),function(i){VAR[i,is[i]]},v)

    if(length(p.dt))
    {
      # TODO
      # TODO
      # TODO

      p.dt <- NULL
    }
  }

  return(list(v.t=v,VAR.t=VAR,v.dt=v.dt,VAR.dt=VAR.dt,p.t=p))
}


# estimate the speed between the two rows of data with error UERE & temporal resolution dt
# dt[1] is recording interval
# dt[2] is minimum time between fixes, which can be smaller than dt[1]
speedMLE <- function(data,dt=NULL,UERE=0,CTMM=ctmm(error=UERE,axes=c("x","y"),circle=TRUE),DT=diff(data$t),SVF=NULL,COV=NULL,n=NULL)
{
  ######################
  # 1/diff(t) estimate

  ZERO <- DT==0
  if(any(ZERO)) { DT[ZERO] <- dt[2] }
  f <- 1/DT
  # if(dt) # truncate
  # {
  #   # result storage
  #   f <- numeric(length(DT))
  #
  #   # assuming two times are uniformly distributed with roundoff error
  #   SUB <- DT>dt
  #   DT.SUB <- DT[SUB]
  #   if(any(SUB)) { f[SUB] <- ( (DT.SUB+dt)*log(DT.SUB+dt) + (DT.SUB-dt)*log(DT.SUB-dt) - 2*DT.SUB*log(DT.SUB) )/dt^2 }
  #
  #   # limit of sampling interval == roundoff
  #   SUB <- DT==dt
  #   if(any(SUB)) { f[SUB] <- log(4)/dt }
  #
  #   # fall back - totally made up
  #   SUB <- DT<dt
  #   if(any(SUB)) { f[SUB] <- 2*log(4)/dt }
  # }
  # else
  # { f <- 1/DT }

  ######################
  # distance estimate

  # measured distances
  dx <- diff(data$x)
  dy <- diff(data$y)
  dr <- sqrt(dx^2+dy^2)
  dx <- cbind(dx,dy)
  rm(dy)

  # step covariances
  if(!is.null(COV))
  {
    for(i in 1:(nrow(data)-1)) { COV[i,,] <- COV[i,,] + COV[i+1,,] + 2*SVF(DT[i]) }

    p <- dr
    for(i in 1:nrow(dx)) { p[i] <- dx[i,] %*% PDsolve(COV[i,,]) %*% dx[i,] }
    p <- 1 - (-expm1(-p/2))^n
  }
  else
  { p <- NULL }
  rm(dx)

  # point estimate of distance with error>0
  if(length(UERE)>1 || UERE)
  {
    if(length(UERE)==1)
    {
      # 1-time errors
      error <- get.error(data,ctmm(error=UERE,axes=c("x","y")),circle=TRUE)
    }
    else # full error array was passed
    {
      UERE -> error
      rm(UERE)
    }

    # 2-time errors
    error <- error[-1] + error[-length(error)]

    DR <- distanceMLE(dr,error)
    VAR <- error/(2-(dr^2-DR^2)/error)
  }
  else
  {
    DR <- dr
    VAR <- numeric(length(dr))
  }

  RETURN <- data.frame(X=DR*f,VAR=VAR*f^2)
  if(length(p)) { RETURN$p <- p }
  return(RETURN)
}


####################
distanceMLE <- function(dr,error)
{
  SUB <- dr>0 & error>0

  # coefficient in transcendental Bessel equation
  # x I0(x) == y I1(x)
  y <- dr[SUB]^2/error[SUB]
  x <- BesselSolver(y)
  # x = dr*dR/error

  # fixed for Inf error
  if(any(SUB)) { dr[SUB] <- ifelse(error[SUB]<Inf,error[SUB]/dr[SUB]*x,0) }

  return(dr)
}


# x I0(x) == y I1(x)
BesselSolver <- function(y)
{
  # solution storage
  x <- numeric(length(y))

  # critical point, below which all point estimates are zero
  # SUB1 <- (y<=2)
  # if(any(SUB1)) { dr[SUB1] <- 0 }
  # x[SUB1] <- 0 # this is now default

  # perturbation of sqrt(EQ) of both sides (from x=0) and solving
  SUB2 <- (2<y) & (y<3.6)
  if(any(SUB2))
  {
    x.SUB <- sqrt(2*y[SUB2])
    x[SUB2] <- 4 * sqrt( (x.SUB-2)/(4-x.SUB) )
  }

  # expansion EQ/exp (from x=y) and solving
  SUB3 <- (3.6<=y) & (y<=BESSEL_LIMIT)
  if(any(SUB3))
  {
    y.SUB <- y[SUB3]

    BI0 <- besselI(y.SUB,0,expon.scaled=TRUE)
    BI1 <- besselI(y.SUB,1,expon.scaled=TRUE)

    x[SUB3] <- 2*y.SUB * ( y.SUB*BI0 - (y.SUB+1)*BI1 ) / ( (2*y.SUB-1)*BI0 - (2*y.SUB+1)*BI1 )
  }

  # expansion from x=y, solving, and then rational expansion from y=Inf
  SUB4 <- BESSEL_LIMIT<y
  if(any(SUB4))
  {
    y.SUB <- y[SUB4]

    x[SUB4] <- 2*y.SUB*(y.SUB-1)/(2*y.SUB-1)
  }

  # iterative to solution by expanding EQ/exp from x=x.old and solving
  SUB <- SUB2 & SUB3
  if(any(SUB))
  {
    x.SUB <- x[SUB]
    y.SUB <- y[SUB]

    BI0 <- besselI(x.SUB,0,expon.scaled=TRUE)
    BI1 <- besselI(x.SUB,1,expon.scaled=TRUE)

    x[SUB] <- x.SUB * ( x.SUB*(x.SUB+y.SUB)*BI0 - (x.SUB^2+(x.SUB+2)*y.SUB)*BI1 ) / ( x.SUB*(x.SUB+y.SUB-1)*BI0 - (x.SUB^2+(x.SUB+1)*y.SUB)*BI1 )
  }
  # one iteration is good enough to get everything below 1% relative error

  # this is now the point estimate, including telemetry error
  # SUB <- !SUB1

  return(x)
}


# estimate temporal resolution from data
# dt[1] = recording interval (truncation)
# dt[2] = fix time
time_res <- function(DT)
{
  ### recording interval ###
  dt <- DT[DT>0]
  # fallback values
  if(length(dt)==0) { return(c(1,1)) }
  if(length(dt)==1) { return(c(dt,min(1,dt/2))) }

  # assume decimal truncation
  M <- -log10(min(dt))
  M <- ceiling(M)
  M <- max(0,M)
  M <- 10^M

  dt <- round(M*dt) # shift decimal place to integers
  dt <- gcd.vec(dt)/M # shift back

  ### fix time ###
  # longest repeating 1
  DT <- DT==0
  if(any(DT))
  {
    for(i in 2:length(DT)) { if(DT[i]) { DT[i] <- DT[i] + DT[i-1] } }
    DT <- max(DT)
    DT <- dt/(1+DT)
  }
  else
  { DT <- 0 }

  return(c(dt,DT))
}


# greatest common divisor of an array
gcd.vec <- function(vec)
{
  vec <- sort(vec,method='quick')

  GCD <- gcd(vec[1],vec[2])
  for(i in vec[-(1:2)]) { GCD <- gcd(i,GCD) }
  # i > GCD because we sorted vec first

  return(GCD)
}

# greatest common divisor of 2 numbers
# fastest if x > y, I think...
gcd <- function(x,y)
{
  r <- x%%y;
  return(ifelse(r, gcd(y, r), y))
}
