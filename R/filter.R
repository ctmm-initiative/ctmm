# net speeds for the line segments between points
mid_speeds <- function(data,dt=time_res(data),UERE=0)
{
  v <- sapply(1:(length(data$t)-1) , function(i) { speedMLE(data[i+0:1,],dt=dt,UERE=UERE) } )  
  
  return(v)
}

# minimum nearby speed
min_speeds <- function(data,dt=time_res(data),UERE=0)
{
  # inner speed estimates
  v <- mid_speeds(data,dt=dt,UERE=UERE)
  
  # end point contingency
  v1 <- speedMLE(data[c(1,3),],dt=dt,UERE=UERE)
  n <- length(data$t)
  v2 <- speedMLE(data[c(n-2,n),],dt=dt,UERE=UERE)
  
  # before and after estimates
  v1 <- c(v1,v)
  v2 <- c(v,v2)
  
  # minimum adjacent estimate
  v <- pmin(v1,v2)
  
  return(v)
}

# estimate the speed between the two rows of data with error UERE & temporal resolution dt
# data is assumed to have 2 rows !!! will try to vectorize later
# TOL is relative error tolerance for Bessel equation
speedMLE <- function(data,dt=0,UERE=0,CTMM=ctmm(error=UERE,axes=c("x","y")),TOL=0.01)
{
  # how far you can go before R's bessel function breaks down
  BESSEL_LIMIT <- 2^16
  
  # measured distance
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
    
    if(y<=2) { return(0) } # critical point, below which all point estimates are zero
    
    # first guess
    if(y<3.6) # perturbation of sqrt(EQ) of both sides (from x=0) and solving
    {
      x <- sqrt(2*y)
      x <- 4 * sqrt( (x-2)/(4-x) )
      
      BI0 <- BI1 <- 1 # initialize to pass test
    }
    else if(y<BESSEL_LIMIT) # expansion EQ/exp (from x=y) and solving
    {
      BI0 <- besselI(y,0,expon.scaled=TRUE)
      BI1 <- besselI(y,1,expon.scaled=TRUE)
      
      x <- 2*y * ( y*BI0 - (y+1)*BI1 ) / ( (2*y-1)*BI0 - (2*y+1)*BI1 )
    }
    else # expansion from x=y, solving, and then rational expansion from y=Inf
    {
      x <- 2*y*(y-1)/(2*y-1)
    }
   
    # iterative to solution by expanding EQ/exp from x=x.old and solving
    ERROR <- abs(x-y)/x
    while(y<BESSEL_LIMIT && ERROR>TOL)
    {
      x.old <- x
      
      BI0 <- besselI(x,0,expon.scaled=TRUE)
      BI1 <- besselI(x,1,expon.scaled=TRUE)
      
      x <- x * ( x*(x+y)*BI0 - (x^2+(x+2)*y)*BI1 ) / ( x*(x+y-1)*BI0 - (x^2+(x+1)*y)*BI1 )
      
      ERROR <- abs(x-x.old)/x
    }
    
    # this is now the point estimate, including telemetry error
    dr <- error/dr * x
  }
  
  DT <- diff(data$t)
  # point estimate of frequency without/with truncation error accounted for
  if(dt==0) # no truncation error
  { f <- 1/DT }
  else # truncation error
  { f <- log((DT+dt/2)/(DT-dt/2))/dt }
  
  return(dr*f)
}


# estimate temporal resolution from data
time_res <- function(data)
{
  dt <- diff(data$t)
  dt <- dt[dt>0]
  
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
