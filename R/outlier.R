# how far you can go before R's bessel function breaks down
BESSEL_LIMIT <- 2^16

# estimate and assign speeds to times
outlie <- function(data,UERE=10,standardize=FALSE,plot=TRUE,...)
{
  if(class(data)=="list") { return( lapply(data, function(d){outlie(d,UERE=UERE,standardize=standardize,plot=plot)} ) ) }

  error <- get.error(data,ctmm(error=UERE,axes=c("x","y")),circle=TRUE)

  Vs <- assign_speeds(data,UERE=error)
  v <- Vs$v.t

  mu <- median.telemetry(data)

  d <- get.telemetry(data,axes=c("x","y"))
  mu <- get.telemetry(mu,axes=c("x","y"))
  mu <- c(mu)

  # detrend median
  d <- t(d) - mu
  # distances from median
  d <- colSums(d^2)
  d <- sqrt(d)
  d <- distanceMLE(d,error)

  if(plot)
  {
    # bounding box
    new.plot(data,...)
    # convert units on telemetry object to match base plot
    data <- unit.telemetry(data,length=get("x.scale",pos=plot.env))

    lwd <- Vs$v.dt
    lwd <- (lwd/max(lwd))
    col <- grDevices::rgb(0,0,lwd,lwd)
    lwd <- 2*lwd
    n <- length(data$t)
    graphics::segments(x0=data$x[-n],y0=data$y[-n],x1=data$x[-1],y1=data$y[-1],col=col,lwd=lwd,asp=1,...)

    cex <- (d/max(d))
    col <- grDevices::rgb(cex,0,0,cex)
    graphics::points(data$x,data$y,col=col,cex=cex,pch=20,...)
  }

  if(standardize)
  {
    v <- v/stats::mad(v)
    d <- d/stats::mad(d)
  }

  return(data.frame(speed=v,distance=d))
}


# speeds assigned by blame
# dt[1] is recording interval
# dt[2] is minimum time between fixes, which can be smaller than dt[1]
assign_speeds <- function(data,dt=NULL,UERE=0,method="max")
{
  method <- match.arg(method,c("max","min"))

  DT <- diff(data$t)
  if(is.null(dt)) { dt <- time_res(DT) }

  # inner speed estimates
  v.dt <- speedMLE(data,dt=dt,UERE=UERE,DT=DT)
  if(length(v.dt)==1)
  {
    v <- c(v.dt,v.dt)
    return(list(v.t=v,v.dt=v.dt))
  }

  if(length(UERE)==1)
  {
    # end point contingency estimates
    v1 <- speedMLE(data[c(1,3),],dt=dt,UERE=UERE)
    n <- length(data$t)
    v2 <- speedMLE(data[c(n-2,n),],dt=dt,UERE=UERE)
  }
  else # pull out correct errors for calculation if fed all errors
  {
    v1 <- speedMLE(data[c(1,3),],dt=dt,UERE=UERE[c(1,3)])
    n <- length(data$t)
    v2 <- speedMLE(data[c(n-2,n),],dt=dt,UERE=UERE[c(n-2,n)])
  }

  # left and right estimates - n estimates, 1 lag apart
  v1 <- c(v1,v.dt)
  v2 <- c(v.dt,v2)

  if(method=="max")
  {
    # n-1 estimates, 2 lags apart
    v1 <- v1[-n]
    v2 <- v2[-1]

    # which side of lag index looks smaller
    LESS <- v1 < v2

    vs <- sort(v.dt,method='quick',index.return=TRUE,decreasing=TRUE)
    is <- vs$ix
    vs <- vs$x

    v <- numeric(length(LESS))
    # assign blame for smallest speeds, from greatest to least, last/least taking precidence in R
    LESS <- LESS[is]
    v[is+!LESS] <- vs

    # assign blame for largest speeds, from least to greatest, last/greatest taking precidence in R
    LESS <- rev(LESS)
    is <- rev(is)
    vs <- rev(vs)
    v[is+LESS] <- vs
  }
  else if(method=="min")
  {
    v <- pmin(v1,v2)
  }

  return(list(v.t=v,v.dt=v.dt))
}


# estimate the speed between the two rows of data with error UERE & temporal resolution dt
# dt[1] is recording interval
# dt[2] is minimum time between fixes, which can be smaller than dt[1]
speedMLE <- function(data,dt=NULL,UERE=0,CTMM=ctmm(error=UERE,axes=c("x","y"),circle=TRUE),DT=diff(data$t))
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
  dr <- sqrt(diff(data$x)^2+diff(data$y)^2)

  # point esitmate of distance with error>0
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

    dr <- distanceMLE(dr,error)
  }

  return(dr*f)
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

  if(any(SUB)) { dr[SUB] <- error[SUB]/dr[SUB] * x }

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
