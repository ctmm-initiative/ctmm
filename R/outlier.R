# how far you can go before R's bessel function breaks down
BESSEL_LIMIT <- 2^16

# estimate and assign speeds to times
outlie <- function(data,UERE=10,standardize=FALSE,plot=TRUE,by='d',...)
{
  if(class(data)[1]=="list")
  { return( lapply(1:length(data), function(i){OUT <- outlie(data[[i]],UERE=UERE,standardize=standardize,plot=plot); graphics::title(names(data)[i]); OUT} ) ) }

  error <- get.error(data,ctmm(error=UERE,axes=c("x","y")),circle=TRUE) # VAR

  # time resolution
  DT <- diff(data$t)
  time.res <- time_res(DT)
  ZERO <- DT==0
  if(any(ZERO)) { DT[ZERO] <- time.res[2] }

  a <- assign_accelerations(data,UERE=error,DT=DT)
  VAR.a <- a$VAR
  a <- a$a

  Vs <- assign_speeds(data,UERE=error,DT=DT)
  v <- Vs$v.t
  VAR.v <- Vs$VAR.t

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

    by <- get(by)
    cex <- if(diff(range(by))) { by/max(by) } else { 0 }
    col <- grDevices::rgb(cex,0,0,cex)
    graphics::points(data$x,data$y,col=col,cex=cex,pch=20,...)
  }

  if(standardize)
  {
    MAD.d <- stats::mad(d)
    d <- d/MAD.d
    VAR.d <- VAR.d/MAD.d^2

    MAD.v <- stats::mad(v)
    v <- v/MAD.v
    VAR.v <- VAR.v/MAD.v^2

    MAD.a <- stats::mad(a)
    a <- a/MAD.a
    VAR.a <- VAR.a/MAD.a^2
  }

  VAR.v <- nant(VAR.v,0)
  VAR.d <- nant(VAR.d,0)

  R <- data.frame(t=data$t,distance=d,speed=v,acceleration=a,VAR.speed=VAR.v,VAR.distance=VAR.d,VAR.acceleration=VAR.a)
  R <- new.outlie(R)
  return(R)
}


#########################
# plot outlier information with error bars
plot.outlie <- function(x,level=0.95,units=TRUE,axes=c('d','v'),...)
{
  n <- nrow(x)
  t <- x$t
  d <- x$distance
  v <- x$speed
  a <- x$acceleration

  # CIs
  d <- sapply(1:n,function(i){tnorm.hdr(d[i],x$VAR.distance[i],level=level)})
  v <- sapply(1:n,function(i){tnorm.hdr(v[i],x$VAR.speed[i],level=level)})
  a <- sapply(1:n,function(i){tnorm.hdr(a[i],x$VAR.acceleration[i],level=level)})

  # unit conversions ... (simpler)
  UNITS <- unit(t,dimension='time',concise=FALSE,SI=!units)
  t <- t/UNITS$scale
  t.lab <- paste0("Time (",UNITS$name,")")

  UNITS <- unit(d,dimension='length',concise=TRUE,SI=!units)
  d <- d/UNITS$scale
  d.lab <- paste0("Core deviation (",UNITS$name,")")

  UNITS <- unit(v,dimension='length',concise=TRUE,SI=!units)
  v <- v/UNITS$scale
  v.lab <- paste0("Minimum speed (",UNITS$name,"/s)")

  UNITS <- unit(a,dimension='length',concise=TRUE,SI=!units)
  a <- a/UNITS$scale
  a.lab <- paste0("Minimum acceleration (",UNITS$name,"/s\u00B2)")

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
assign_speeds <- function(data,DT=NULL,UERE=0,method="max")
{
  method <- match.arg(method,c("max","min"))

  if(is.null(DT))
  {
    DT <- diff(data$t)
    dt <- time_res(DT)
  }

  # inner speed estimates
  v.dt <- speedMLE(data,UERE=UERE,DT=DT)
  VAR.dt <- v.dt$VAR
  v.dt <- v.dt$X
  if(length(v.dt)==1)
  {
    v <- c(v.dt,v.dt)
    VAR <- c(VAR.dt,VAR.dt)
    return(list(v.t=v,VAR.t=VAR,v.dt=v.dt,VAR.dt=VAR.dt))
  }

  N <- length(data$t)
  if(length(UERE)==1)
  {
    # end point contingency estimates
    v1 <- speedMLE(data[c(1,3),],DT=sum(DT[1:2]),UERE=UERE)
    v2 <- speedMLE(data[c(N-2,N),],DT=sum(DT[(N-1):N]),UERE=UERE)
  }
  else # pull out correct errors for calculation if fed all errors
  {
    v1 <- speedMLE(data[c(1,3),],DT=sum(DT[1:2]),UERE=UERE[c(1,3)])
    v2 <- speedMLE(data[c(N-2,N),],DT=sum(DT[(N-2):(N-1)]),UERE=UERE[c(N-2,N)])
  }
  VAR1 <- v1$VAR; v1 <- v1$X
  VAR2 <- v2$VAR; v2 <- v2$X

  # left and right estimates - n estimates, 1 lag apart
  v1 <- c(v1,v.dt)
  v2 <- c(v.dt,v2)

  VAR1 <- c(VAR1,VAR.dt)
  VAR2 <- c(VAR.dt,VAR2)

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
  }
  else if(method=="min")
  {
    # v <- pmin(v1,v2)

    v <- cbind(v1,v2)
    is <- apply(v,1,which.min)
    v <- vapply(1:nrow(v),function(i){v[i,is[i]]},v[,1])
    VAR <- vapply(1:nrow(VAR),function(i){VAR[i,is[i]]},v)
  }

  return(list(v.t=v,VAR.t=VAR,v.dt=v.dt,VAR.dt=VAR.dt))
}


# accelerations assigned by midpoint time
assign_accelerations <- function(data,DT=NULL,UERE=0)
{
  if(is.null(DT))
  {
    DT <- diff(data$t)
    dt <- time_res(DT)
  }

  if(nrow(data)==2)
  {
    STUFF <- assign_speeds(data,DT=DT,UERE=UERE)
    a <- STUFF$v.t / DT
    VAR <- STUFF$VAR.t / DT^2
    return(list(a=a,VAR=VAR))
  }

  a <- accelerationMLE(data,DT=DT,UERE=UERE)
  VAR <- a$VAR
  a <- a$X

  # end point accelerations (from rest)
  N <- length(data$t)
  if(length(UERE)==1)
  {
    a1 <- assign_speeds(data[1:2,],DT=DT[1],UERE=UERE)
    a2 <- assign_speeds(data[(N-1):N,],DT=DT[N-1],UERE=UERE)
  }
  else
  {
    a1 <- assign_speeds(data[1:2,],DT=DT[1],UERE=UERE[1:2])
    a2 <- assign_speeds(data[(N-1):N,],DT=DT[N-1],UERE=UERE[(N-1):N])
  }
  VAR1 <- a1$VAR.dt / DT[1]^2;   a1 <- a1$v.dt / DT[1]
  VAR2 <- a2$VAR.dt / DT[N-1]^2; a2 <- a2$v.dt / DT[N-1]

  a <- c(a1,a,a2)
  VAR <- c(VAR1,VAR,VAR2)

  return(list(a=a,VAR=VAR))
}


# estimate the speed between the two rows of data with error UERE & temporal resolution dt
# dt[1] is recording interval
# dt[2] is minimum time between fixes, which can be smaller than dt[1]
speedMLE <- function(data,DT=NULL,UERE=0,CTMM=ctmm(error=UERE,axes=c("x","y"),circle=TRUE))
{
  ######################
  # 1/diff(t) estimate

  if(is.null(DT))
  {
    DT <- diff(data$t)
    dt <- time_res(DT)
  }
  f <- 1/DT

  ######################
  # distance estimate

  # measured distances
  dx <- diff(data$x)
  dy <- diff(data$y)
  dr <- sqrt(dx^2+dy^2)
  rm(dx,dy)

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
  return(RETURN)
}

accelerationMLE <- function(data,DT=NULL,UERE=0,CTMM=ctmm(error=UERE,axes=c("x","y"),circle=TRUE))
{
  ######################
  # 1/diff(t) estimate

  if(is.null(DT))
  {
    DT <- diff(data$t)
    dt <- time_res(DT)
  }

  # calculate velocities
  # DT is velocity time step
  f <- 1/DT
  vx <- diff(data$x) * f
  vy <- diff(data$y) * f
  w <- cbind(f[-length(f)],f[-length(f)]+f[-1],f[-1]) # velocity location weights

  # calculate accelerations
  DT <- (DT[-1] + DT[-length(DT)])/2 # acceleration time step
  f <- 1/DT
  ax <- diff(vx) * f; rm(vx)
  ay <- diff(vy) * f; rm(vy)

  a <- sqrt(ax^2+ay^2)
  rm(ax,ay)

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

    # 3-time errors
    ERROR <- array(0,length(error)-2)
    w <- w^2
    for(i in 1:(length(error)-2)) { ERROR[i] <- w[i,] %*% error[i + 0:2] }
    rm(error)
    ERROR <- ERROR * f^2

    A <- distanceMLE(a,ERROR)
    VAR <- ERROR/(2-(a^2-A^2)/ERROR)
  }
  else
  {
    A <- a
    VAR <- numeric(length(a))
  }

  RETURN <- data.frame(X=A,VAR=VAR)
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
