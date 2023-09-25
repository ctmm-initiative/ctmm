# how far you can go before R's bessel function breaks down
BESSEL_LIMIT <- 2^16

# estimate and assign speeds to times
outlie <- function(data,plot=TRUE,by='d',...)
{
  if(class(data)[1]=="list")
  { return( lapply(1:length(data), function(i){OUT <- outlie(data[[i]],plot=plot,by=by); if(plot){ graphics::title(names(data)[i]) }; OUT} ) ) }

  UERE <- uere(data)
  if(DOP.LIST$horizontal$VAR %nin% names(data)) { uere(data) <- UERE } # adds VAR guesstimate columns if missing (DOF==0)
  error <- UERE$UERE[,'horizontal']
  names(error) <- rownames(UERE$UERE) # R drops dimnames
  error <- ctmm(error=error,axes=c('x','y'))
  error <- get.error(data,error,calibrate=TRUE)

  # time resolution
  DT <- diff(data$t)
  time.res <- time_res(DT)
  ZERO <- DT==0
  if(any(ZERO)) { DT[ZERO] <- time.res[2] }

  Vs <- assign_speeds(data,UERE=error,DT=DT,axes=c('x','y'))
  v <- Vs$v.t
  VAR.v <- Vs$VAR.t

  mu <- median.telemetry(data)
  d <- get.telemetry(data,axes=c("x","y"))
  mu <- get.telemetry(mu,axes=c("x","y"))
  mu <- c(mu)

  # detrend median
  d <- t(d) - mu
  # distances from median
  if(length(dim(error))==3)
  { d <- t(d) }
  else
  {
    d <- colSums(d^2)
    d <- sqrt(d)
  }
  D <- distanceMLE(d,error,return.VAR=TRUE)
  d <- D[,1]
  VAR.d <- D[,2]
  rm(D)

  if('z' %in% names(data))
  {
    error <- UERE$UERE[,'vertical']
    names(error) <- rownames(UERE$UERE) # R drops dimnames
    error <- ctmm(error=error,axes=c('z'))
    error <- get.error(data,error,calibrate=TRUE) # VAR

    Vz <- assign_speeds(data,UERE=error,DT=DT,axes='z')
    vz <- Vz$v.t
    VAR.vz <- Vz$VAR.t

    dz <- get.telemetry(data,axes=c("z"))
    dz <- dz - stats::median(data$z)
    dz <- abs(dz)
    DZ <- distanceMLE(dz,error,axes='z',return.VAR=TRUE)
    dz <- DZ[,1]
    VAR.dz <- DZ[,2]
    rm(DZ)
  }

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

  VAR.v <- nant(VAR.v,0)
  VAR.d <- nant(VAR.d,0)
  if('z' %in% names(data))
  {
    VAR.vz <- nant(VAR.vz,0)
    VAR.dz <- nant(VAR.dz,0)
  }

  R <- data.frame(t=data$t,distance=d,VAR.distance=VAR.d,speed=v,VAR.speed=VAR.v)
  if('z' %in% names(data))
  {
    R$vertical.distance <- dz
    R$VAR.vertical.distance <- VAR.dz
    R$vertical.speed <- vz
    R$VAR.vertical.speed <- VAR.vz
  }
  R <- new.outlie(R)
  return(R)
}


#########################
# plot outlier information with error bars
plot.outlie <- function(x,level=0.95,units=TRUE,axes=c('d','v'),xlim=NULL,ylim=NULL,...)
{
  x <- listify(x)
  t <- NULL
  d <- v <- NULL
  dz <- vz <- NULL

  VERTICAL <- any(c('z','vz') %in% axes)

  for(i in 1:length(x))
  {
    n <- nrow(x[[i]])

    if(VERTICAL && "vertical.distance" %in% names(x[[i]]))
    {
      dzi <- sapply(1:n,function(j){tnorm.hdr(x[[i]]$vertical.distance[j],x[[i]]$VAR.vertical.distance[j],level=level)})
      dz <- cbind(dz,dzi)

      vzi <- sapply(1:n,function(j){tnorm.hdr(x[[i]]$vertical.speed[j],x[[i]]$VAR.vertical.speed[j],level=level)})
      vz <- cbind(vz,vzi)
    }
    else if(VERTICAL)
    { next }

    t <- c(t, x[[i]]$t )

    di <- sapply(1:n,function(j){tnorm.hdr(x[[i]]$distance[j],x[[i]]$VAR.distance[j],level=level)})
    d <- cbind(d,di)

    vi <- sapply(1:n,function(j){tnorm.hdr(x[[i]]$speed[j],x[[i]]$VAR.speed[j],level=level)})
    v <- cbind(v,vi)
  }
  n <- length(t)

  # unit conversions ... (simpler)
  UNITS <- unit(t,dimension='time',concise=FALSE,SI=!units)
  t.scale <- UNITS$scale
  t <- t/UNITS$scale
  t.lab <- paste0("Time (",UNITS$name,")")

  UNITS <- unit(d,dimension='length',concise=TRUE,SI=!units)
  d.scale <- UNITS$scale
  d <- d/UNITS$scale
  d.lab <- paste0("Median deviation (",UNITS$name,")")

  UNITS <- unit(v,dimension='length',concise=TRUE,SI=!units)
  v.scale <- UNITS$scale
  v <- v/UNITS$scale
  v.lab <- paste0("Minimum speed (",UNITS$name,"/s)")

  if(VERTICAL)
  {
    UNITS <- unit(dz,dimension='length',concise=TRUE,SI=!units)
    dz.scale <- UNITS$scale
    dz <- dz/UNITS$scale
    dz.lab <- paste0("Median vertical deviation (",UNITS$name,")")

    UNITS <- unit(vz,dimension='length',concise=TRUE,SI=!units)
    vz.scale <- UNITS$scale
    vz <- vz/UNITS$scale
    vz.lab <- paste0("Minimum vertical speed (",UNITS$name,"/s)")
  }

  mid <- function(x)
  {
    if(length(dim(x))) { x <- x[2,] }
    return(x)
  }

  x <- get(axes[1])
  y <- get(axes[2])

  xlab <- get(paste0(axes[1],".lab"))
  ylab <- get(paste0(axes[2],".lab"))

  if(is.null(xlim))
  {
    xlim <- x
    xlim[3,] <- pmin(x[3,],x[2,]*2)
    xlim <- range(xlim)
  }
  else
  {
    SCALE <- get(paste0(axes[1],".scale"))
    xlim <- xlim/SCALE
  }

  if(is.null(ylim))
  {
    ylim <- y
    ylim[3,] <- pmin(y[3,],y[2,]*2)
    ylim <- range(ylim)
  }
  else
  {
    SCALE <- get(paste0(axes[2],".scale"))
    ylim <- ylim/SCALE
  }

  # base plot
  plot(mid(x),mid(y),xlim=xlim,ylim=ylim,pch=19,xlab=xlab,ylab=ylab,...)

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
assign_speeds <- function(data,DT=NULL,UERE=0,method="max",axes=c('x','y'))
{
  method <- match.arg(method,c("max","min"))

  if(is.null(DT))
  {
    DT <- diff(data$t)
    dt <- time_res(DT)
  }

  # inner speed estimates
  v.dt <- speedMLE(data,UERE=UERE,DT=DT,axes=axes)
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
  else if(length(dim(UERE))==3) # pull out correct errors for calculation if fed all errors
  {
    v1 <- speedMLE(data[c(1,3),],DT=sum(DT[1:2]),UERE=UERE[c(1,3),,])
    v2 <- speedMLE(data[c(N-2,N),],DT=sum(DT[(N-2):(N-1)]),UERE=UERE[c(N-2,N),,])
  }
  else # pull out correct errors for calculation if fed all errors (circular)
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


# estimate the speed between the two rows of data with error UERE & temporal resolution dt
# dt[1] is recording interval
# dt[2] is minimum time between fixes, which can be smaller than dt[1]
speedMLE <- function(data,DT=NULL,UERE=0,axes=c('x','y'),CTMM=ctmm(error=UERE,axes=axes))
{
  AXES <- length(axes)

  ######################
  # 1/diff(t) estimate

  if(is.null(DT))
  {
    DT <- diff(data$t)
    dt <- time_res(DT)
  }
  f <- 1/DT

  ###################
  if(length(UERE)==1) # 1-time errors
  { error <- get.error(data,ctmm(error=UERE,axes=axes)) }
  else # full error array was passed
  { error <- UERE }

  ######################
  # distance estimate

  # measured distances
  if(AXES==2)
  {
    dx <- diff(data$x)
    dy <- diff(data$y)

    if(length(dim(error))==3) # error ellipse
    { dr <- cbind(dx,dy) } # need full vector
    else # error circle or interval
    { dr <- sqrt(dx^2+dy^2) }
    rm(dx,dy)
  }
  else if(AXES==1)
  {
    dr <- diff(data$z)
    dr <- abs(dr)
  }

  # point estimate of distance with error>0
  if(length(UERE)>1 || UERE)
  {
    # 2-time errors
    if(length(dim(error))==3) # error ellipse
    { error <- error[-1,,,drop=FALSE] + error[-nrow(error),,,drop=FALSE] }
    else # error circle or interval
    { error <- error[-1] + error[-length(error)] }

    DR <- distanceMLE(dr,error,axes=axes,return.VAR=TRUE)
    VAR <- DR[,2]
    DR <- DR[,1]
  }
  else
  {
    DR <- dr
    VAR <- numeric(length(dr))
  }

  RETURN <- data.frame(X=DR*f,VAR=VAR*f^2)
  return(RETURN)
}


####################
distanceMLE <- function(dr,error,axes=c('x','y'),return.VAR=FALSE)
{
  if(length(dim(error))==3) # error ellipse
  {
    dr <- sapply(1:nrow(dr),function(i){ abs_bivar(dr[i,],error[i,,],return.VAR=TRUE) })
    dr <- t(dr) #
  }
  else # error circle or interval
  {
    AXES <- length(axes)
    DR <- dr
    SUB <- dr>0 & error>0

    if(any(SUB))
    {
      if(AXES==2) # error circle
      {
        # coefficient in transcendental Bessel equation
        # x I0(x) == y I1(x)
        y <- DR[SUB]^2/error[SUB]
        x <- BesselSolver(y)
        # x = DR*dR/error

        # fixed for Inf error
        DR[SUB] <- ifelse(error[SUB]<Inf,error[SUB]/DR[SUB]*x,0)
      }
      else if(AXES==1) # error interval
      {
        error[SUB] <- sqrt(error[SUB]) # now standard deviation

        # x = y tanh(x y)
        y <- DR[SUB]/error[SUB]
        x <- TanhSolver(y)

        DR[SUB] <- ifelse(error[SUB]<Inf,error[SUB]*x,0)
      } # end error interval
    } # end SUB

    VAR <- error/(AXES-(dr^2-DR^2)/error)
    dr <- cbind(DR,VAR)
  } # end error circle or error interval

  if(!return.VAR) { dr <- dr[,1] }
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


# x == y tanh(x y)
TanhSolver <- function(y)
{
  x <- y

  SUB <- y<1
  if(any(SUB)) { x[SUB] <- 0 }

  # Taylor expansion from x=0
  SUB <- 1<y & y<=1.1106677241426588
  if(any(SUB))
  {
    x[SUB] <- sqrt((5*y[SUB]-sqrt(120-95*y[SUB]^2))/y[SUB]^3)/2
  }

  # Taylor expansion from x=y
  SUB <- 1.1106677241426588<=y & y<Inf
  if(any(SUB))
  {
    TEMP <- 2*y[SUB]^2
    x[SUB] <- nant((sinh(TEMP)-TEMP)/(cosh(TEMP)-TEMP+1),1) * y[SUB]
  }

  # One Newton-Raphson iteration and we are within 0.5% error max
  SUB <- 1<y & y<Inf
  if(any(SUB))
  {
    TANH <- tanh(x[SUB]*y[SUB])
    GRAD <- (1-TANH^2) * y[SUB]
    # x + dx == y (TANH + GRAD*dx) + O(dx^2)
    # (1 - y*GRAD) dx == y*TANH - x + O(dx^2)
    dx <- (y[SUB]*TANH-x[SUB])/(1-y[SUB]*GRAD)
    x[SUB] <- x[SUB] + dx
  }

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
