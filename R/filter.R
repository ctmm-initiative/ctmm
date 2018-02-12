median.telemetry <- function(x,na.rm=FALSE,...)
{
  x <- listify(x)

  id <- paste0("median of ",mean.info(x)$identity)
  proj <- projection(x[[1]])
  tz <- attr(x[[1]],'info')$timezone
  x <- do.call(rbind,x)

  long <- stats::median(x$longitude)
  lat <- stats::median(x$latitude)
  t <- stats::median(x$t)
  timestamp <- as.character(as.POSIXct(t,tz=tz,origin="1970/01/01"))

  x <- data.frame(timestamp=timestamp,t=t,longitude=long,latitude=lat)
  x <- new.telemetry(x,info=list(identity=id,projection=proj,timezone=tz))
  projection(x) <- proj

  return(x)
}


# UNFINISHED
distances <- function(x,y=median(x),error=NULL)
{
  x <- listify(x)
  y <- listify(y)
  if(length(x)>length(y))
  y <- rep(y,length(x)/length(y))

  d <- list()
  for(i in 1:length(x))
  {
    n <- length(x[[i]])

  }


}


#
speeds <- function(data,error=NULL,dt=NULL)
{
  if(is.null(error)) { error <- 10 } # GPS default 10 meter error
  if(is.null(dt)) { dt <- time_res(data) }
  assign_speeds(data,UERE=error,dt=dt)$v.t
}


# speeds assigned by blame
assign_speeds <- function(data,dt=time_res(data),UERE=0,method=c("max","min"))
{
  method <- match.arg(method,c("max","min"))

  CTMM <- ctmm(error=UERE,axes=c("x","y"))

  # inner speed estimates
  v.dt <- speedMLE(data,dt=dt,UERE=UERE,CTMM=CTMM)
  if(length(v.dt)==1)
  {
    v <- c(v.dt,v.dt)
    return(list(v.t=v,v.dt=v.dt))
  }

  # end point contingency estimates
  v1 <- speedMLE(data[c(1,3),],dt=dt,UERE=UERE,CTMM=CTMM)
  n <- length(data$t)
  v2 <- speedMLE(data[c(n-2,n),],dt=dt,UERE=UERE,CTMM=CTMM)

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
speedMLE <- function(data,dt=time_res(data),UERE=0,CTMM=ctmm(error=UERE,axes=c("x","y")))
{
  ######################
  # 1/diff(t) estimate

  DT <- diff(data$t)
  if(dt) # truncate
  {
    # result storage
    f <- numeric(length(DT))

    # assuming two times are uniformly distributed with roundoff error
    SUB <- DT>dt
    DT.SUB <- DT[SUB]
    if(any(SUB)) { f[SUB] <- ( (DT.SUB+dt)*log(DT.SUB+dt) + (DT.SUB-dt)*log(DT.SUB-dt) - 2*DT.SUB*log(DT.SUB) )/dt^2 }

    # limit of sampling interval == roundoff
    SUB <- DT==dt
    if(any(SUB)) { f[SUB] <- log(4)/dt }

    # fall back - totally made up
    SUB <- DT<dt
    if(any(SUB)) { f[SUB] <- 2*log(4)/dt }
  }
  else
  { f <- 1/DT }

  ######################
  # distance estimate

  # measured distances
  dr <- sqrt(diff(data$x)^2+diff(data$y)^2)

  # point esitmate of distance with error>0
  if(UERE>0)
  {
    # 1-time errors
    error <- get.error(data,CTMM)
    # 2-time errors
    error <- error[-1] + error[-length(error)]

    # coefficient in transcendental Bessel equation
    # x I0(x) == y I1(x)
    y <- dr^2/error

    x <- BesselSolver(y)

    SUB <- dr>0
    if(any(SUB)) { dr[SUB] <- error[SUB]/dr[SUB] * x[SUB] }
  }

  return(dr*f)
}


# x I0(x) == y I1(x)
BesselSolver <- function(y)
{
  # solution storage
  x <- numeric(length(y))
  # how far you can go before R's bessel function breaks down
  BESSEL_LIMIT <- 2^16

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
time_res <- function(data)
{
  dt <- diff(data$t)
  dt <- dt[dt>0]
  if(length(dt)==0) { return(0) }
  if(length(dt)==1) { return(dt) }

  # assume decimal truncation
  M <- -log10(min(dt))
  M <- ceiling(M)
  M <- max(0,M)
  M <- 10^M

  dt <- round(M*dt) # shift decimal place to integers
  dt <- gcd.vec(dt)/M # shift back

  return(dt)
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
