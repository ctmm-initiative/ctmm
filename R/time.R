# accumulated time information assumes <12hr sampling
get.sundial <- function(object,CTMM=NULL,twilight="civil",dt.max=6 %#% 'hr')
{
  twilight <- match.arg(twilight,c('nautical','civil','none'))
  if(twilight=='civil')
  { keep <- c("dawn","solarNoon","dusk","nadir") }
  else if(twilight=='nautical')
  { keep <- c("nauticalDawn","solarNoon","nauticalDusk","nadir") }
  else if(twilight=='none')
  { keep <- c("sunrise","solarNoon","sunset","nadir") }

  INTERPOLATE <- max(diff(object$t)) > dt.max
  if(INTERPOLATE)
  {
    # this will provide linear interpolation with predict
    if(is.null(CTMM))
    {
      CTMM <- ctmm(range=FALSE,tau=Inf,isotropic=TRUE,sigma=1)
      CTMM <- ctmm.loglike(object,CTMM,verbose=TRUE)
    }

    # keep track of new times introduced
    KEEP <- object$t

    # interpolate data based on CTMM if !NULL
    object <- predict(object,CTMM,dt=dt.max,complete=TRUE)

    KEEP <- (object$t %in% KEEP)
  }

  TODAY <- data.frame(date=object$timestamp,lat=object$latitude,lon=object$longitude)
  # calculate local timezones
  tz <- round(object$lon/15)
  tz <- ifelse(tz>=0,paste0("Etc/GMT+",tz),paste0("Etc/GMT",tz))
  # specify date in local timezone
  n <- length(tz)
  TODAY$date <- sapply(1:n,function(i){as.Date(TODAY$date[i],tz=tz[i])})
  # R data.frames are supposed to keep class information like Date...?
  TODAY$date <- as.Date(TODAY$date,origin=EPOCH)

  days <- 2 # buffer of days to run suncalc, because its behavior is unpredictable
  TIMES <- lapply((-days):days,function(d){
    TODAY$date <- TODAY$dat + d;
    TODAY <- suncalc::getSunlightTimes(data=TODAY,keep=keep,tz="UTC")[,keep]
    })
  TIMES <- do.call(cbind,TIMES)

  # so annoying that as.numeric does not work on data.frames
  TIMES <- sapply(1:ncol(TIMES),function(col){as.numeric(TIMES[[col]])})

  # corresponding circular variables
  ANGLE <- rep(c(0,pi/2,pi,3*pi/2),2*days+1)

  # convert data time into circular variable
  angle <- rep(0,n)
  day <- rep(1,n)
  dawn <- dusk <- array(0,c(n,2)) # previous and next
  for(i in 1:n)
  {
    # the problem with the suncalc info is that it comes out unsorted
    IND <- sort(TIMES[i,],index.return=TRUE)$ix
    times <- TIMES[i,IND]
    theta <- ANGLE[IND]

    # encapsulating 2 knots
    IND <- which.min(abs(times-object$t[i]))
    if(object$t[i] < times[IND])
    { IND <- c(IND-1,IND) }
    else
    { IND <- c(IND,IND+1) }

    t <- times[IND]
    a <- theta[IND]

    # sort circular variable for interpolation
    while(a[1]>a[2]) { a[2] <- a[2] + 2*pi }

    if(t[2]==object$t[i]) # exact match
    { angle[i] <- a[2] }
    else # interpolate
    { angle[i] <- a[1] + (object$t[i]-t[1])/diff(t)*diff(a) }
    angle[i] <- angle[i] %% (2*pi)

    # previous dawn
    d <- which(theta==0 & times<=object$t[i])
    dawn[i,1] <- times[last(d)]
    # next dawn
    d <- which(theta==0 & times>object$t[i])
    dawn[i,2] <- times[first(d)]

    # previous dusk
    d <- which(theta==pi & times<=object$t[i])
    dusk[i,1] <- times[last(d)]
    # next dusk
    d <- which(theta==pi & times>object$t[i])
    dusk[i,2] <- times[first(d)]
  }

  # accumulated light & darkness from the beginning
  light <- dark <- rep(0,n)
  for(i in 2:n)
  {
    # accumulate from past accumulated light and darkness
    light[i] <- light[i-1]
    dark[i] <- dark[i-1]

    dt <- object$t[i]-object$t[i-1]
    # da <- (angle[i]-angle[i-1]) %% (2*pi)

    if(angle[i]<pi) # -to-light
    {
      if(angle[i-1]<=pi) # light-to-light
      {
        light[i] <- light[i] + dt
      }
      else # dark-to-light transition
      {
        w <- c( dawn[i-1,2]-object$t[i-1] , object$t[i]-dawn[i,1] )
        w <- w/sum(w)
        # linearly interpolated local dawn time
        t <- w[2]*dawn[i-1,2] + w[1]*dawn[i,1] # reciprocal weighting
        # partition accumulated times
        dark[i] <- dark[i] + (t-object$t[i-1])
        light[i] <- light[i] + (object$t[i]-t)
      }
    } # end -to-light
    else # -to-dark
    {
      if(angle[i-1]>pi) # dark-to-dark
      {
        dark[i] <- dark[i] + dt
      }
      else # light-to-dark transition
      {
        w <- c( dusk[i-1,2]-object$t[i-1] , object$t[i]-dusk[i,1] )
        w <- w/sum(w)
        # linearly interpolated local dusk time
        t <- w[2]*dusk[i-1,2] + w[1]*dusk[i,1] # reciprocal weighting
        # partition accumulated times
        dark[i] <- dark[i] + (t-object$t[i-1])
        light[i] <- light[i] + (object$t[i]-t)
      }
    } # end -to-dark
  } # end accumulation loop

  R <- data.frame(light.time=light,dark.time=dark)
  R$light <- angle>0 & angle<=pi # already modulo 2*pi

  # information for spline timelinks
  R$sundial <- angle
  R$suntime <- ifelse(0<=angle & angle<pi,dusk[,2]-dawn[,1],dawn[,2]-dusk[,1])

  # remove new times
  if(INTERPOLATE) { R <- R[KEEP,] }

  return(R)
}


# get time-linked times
linktime <- function(data,CTMM)
{
  timelink <- CTMM$timelink
  p <- length(CTMM$timelink.par)
  t <- data$t

  if(is.null(timelink) || timelink=="identity" || p==0)
  {  }
  else if(timelink=="switch")
  {
    par <- CTMM$timelink.par
    t <- data$light.time*(1+par) + data$dark.time*(1-par)
  }
  else
  {
    R <- timelink.fn(CTMM)

    # handle initial offset
    if(data$sundial[1]<=pi) # counting time since previous sunrise
    { data$light.time <- data$light.time + data$sundial[1]/pi*data$suntime[1] }
    else # counting time since previous sunset
    { data$dark.time <- data$dark.time + (data$sundial[1]-pi)/pi*data$suntime[1] }

    n <- length(t)
    t <- 0
    for(i in 1:n)
    {
      if(data$sundial[i]<=pi)
      {
        # whole periods accumulated
        t <- R$acc.dark * data$dark.time[i]
        #
        #
      }
      else
      {
        # whole periods accumulated
        t <- R$acc.light * data$light.time[i]
        #
        #
      }
    }

    # t <- t + int(data$sundial)/data$sundial.rate

  }
  # TODO diurnal & nocturnal spline models

  return(t)
}


#
timelink.parinfo <- function(CTMM)
{
  timelink <- CTMM$timelink
  p <- length(CTMM$timelink.par)

  if(is.null(timelink) || timelink=="identity" || p==0)
  { return(list()) }

  fn <- get(paste0(timelink,".timelink.parinfo"))
  fn(CTMM)
}


#
timelink.clean <- function(par,timelink="identity")
{
  p <- length(par)

  if(is.null(timelink) || timelink=="identity" || p==0)
  { return(par) }

  fn <- get(paste0(timelink,".timelink.clean"))
  fn(par)
}


# how fast is time moving (factor for speed) for data
timelink.rate <- function(data,CTMM)
{
  timelink <- CTMM$timelink
  p <- length(CTMM$timelink.par)

  if(is.null(timelink) || timelink=="identity" || p==0)
  { return(rep(1,nrow(data))) }

  if(timelink %in% "switch")
  {
    fn <- get(paste0(timelink,".timelink.rate"))
    r <- fn(data,CTMM)
  }
  else
  {
    R <- timelink.fn(CTMM)
    r <- R$fn(data$sundial)
  }
  return(r)
}


# function for sundial plot and CIs for time argument
timelink.fn <- function(CTMM)
{
  timelink <- CTMM$timelink
  p <- length(CTMM$timelink.par)

  if(is.null(timelink) || timelink=="identity" || p==0)
  {
    R <- list()
    R$fn <- function(angle) { rep(1,length(angle)) }
    R$grad <- function(angle) { rep(0,length(angle)) }
  }
  else
  {
    fn <- get(paste0(timelink,".timelink.fn"))
    R <- fn(CTMM)
  }

  return(R)
}


# increase timelink complexity
timelink.complexify <- function(CTMM)
{
  timelink <- CTMM$timelink

  if(is.null(timelink) || timelink=="identity")
  { return(CTMM) }

  fn <- get(paste0(timelink,".timelink.complexify"))
  fn(CTMM)
}


# decrease timelink complexity
timelink.simplify <- function(CTMM)
{
  timelink <- CTMM$timelink
  p <- length(CTMM$timelink.par)

  if(is.null(timelink) || timelink=="identity" || p==0)
  { return(CTMM) }

  fn <- get(paste0(timelink,".timelink.simplify"))
  fn(CTMM)
}


# return timelink name
timelink.name <- function(CTMM)
{
  timelink <- CTMM$timelink
  p <- length(CTMM$timelink.par)

  if(is.null(timelink) || timelink=="identity" || p==0)
  { return(NULL) }

  fn <- get(paste0(timelink,".timelink.name"))
  fn(CTMM)
}


# summarize parameters
timelink.summary <- function(CTMM,level=0.95)
{
  timelink <- CTMM$timelink
  p <- length(CTMM$timelink.par)

  if(is.null(timelink) || timelink=="identity" || p==0)
  { return(NULL) }

  fn <- get(paste0(timelink,".timelink.summary"))
  fn(CTMM,level=level)
}


###############################
# day & night activity rate time-link model
# day   rate = 1+par
# night rate = 1-par

switch.timelink.parinfo <- function(CTMM)
{
  R <- list()
  R$lower <- -1
  R$upper <- 1
  R$parscale <- 1
  return(R)
}

switch.timelink.clean <- function(par)
{ clamp(par,-1,1) }

# dt_internal/dt_clock
switch.timelink.rate <- function(data,CTMM)
{
  par <- CTMM$timelink.par
  ifelse(data$light,1+par,1-par)
}

# rate = dt_internal/dt_clock
# drate/dpar
switch.timelink.fn <- function(par)
{
  R <- list()

  R$fn <- function(angle) { ifelse(0<=angle & angle<pi,1+par,1-par) }
  R$grad <- function(angle) { ifelse(0<=angle & angle<pi,+1,-1) }

  return(R)
}

switch.timelink.complexify <- function(object)
{
  # if(length(object$timelink.par)==0)
  # { object$timelink.par <- c(object$timelink.par,0) }
  object$timelink.par <- 0
  return(object)
}

switch.timelink.simplify <- function(object)
{
  # n <- length(object$timelink.par)
  # object$timelink.par <- object$timelink.par[-n]
  object$timelink.par <- NULL
  return(object)
}

switch.timelink.name <- function(object)
{
  return("diel-switch")
}

switch.timelink.summary <- function(object,level=0.95)
{
  par <- object$timelink.par

  CI <- c(0,0.5,1)

  if('timelink-1' %in% dimnames(object$COV)[[1]])
  {
    VAR <- object$COV['timelink-1','timelink-1']
    VAR <- VAR/4 # beta ~ 1/2 +/- par/2
  }
  else
  { VAR <- Inf }

  if(par==0)
  {
    NAME <- "% cathemeral"
    CI <- beta.ci(0.5,VAR,level=level)
    CI <- CI * 2 # one sided CI
    CI[3] <- 1
  }
  else if(par>0)
  {
    NAME <- "% diurnal"
    beta <- (1+par)/2
    CI <- beta.ci(beta,VAR,level=level)
  }
  else if(par<0)
  {
    NAME <- "% nocturnal"
    beta <- (1-par)/2
    CI <- beta.ci(beta,VAR,level=level)
  }

  CI <- 100*CI # fraction -> %

  CI <- rbind(CI)
  rownames(CI) <- NAME
  colnames(CI) <- NAMES.CI

  return(CI)
}

##################
# asymmetric periodic cubic splines

spline.timelink.fn <- function(CTMM,even=FALSE,half=FALSE,fast=FALSE)
{
  y <- CTMM$timelink.par
  p0 <- length(y)
  # missing y[0] is noon rate - inferred from mean 100% activity
  y <- c(p0+1-sum(y),y)
  p <- length(y)
  h <- 2*pi/p # angle between knots

  yn <- c(y[-1],y[1]) # next index
  yp <- c(y[p],y[-p]) # previous index

  if(fast) # fast solver (in progress)
  {
    # tri-band circulant matrix to solve for quadratic terms
    M <- array(0,c(p,p))
    for(i in 1:p)
    {
      M[i,i] <- 4
      j <- 1+i%%p
      M[i,j] <- M[i,j] + 1
      j <- 1+(i-2)%%p
      M[i,j] <- M[i,j] + 1
    }
    b <- 3/h^2*( yn - 2*y + yp )
    M <- solve(M) # swap this out for tri-band circulant solver
    C <- M %*% b

    # GRADIENT UNFINISHED

    rm(M,b)
    Cn <- c(C[-1],C[1])

    # linear term
    B <- ( yn - y )/h - h/3*( Cn + 2*C )

    # cubic term
    D <- ( Cn - C )/(3*h)

    # all coefficients
    Q <- cbind(y,B,C,D)
    rm(y,B,C,D)
  }
  else # slow solver
  {
    # Jacobian for gradient
    J <- diag(p)[,-1]
    dim(J) <- c(p,p0)
    J[1,] <- -1

    Jn <- rbind(J[-1,],J[1,])

    M1 <- M2 <- M3 <- array(0,c(p,p,3))

    for(i in 1:p)
    {
      # continuity
      M1[i,i,] <- c(h,h^2,h^3)

      # continuous derivative
      M2[i,i,] <- +c(1,2*h,3*h^2)
      M2[i,1+i%%p,] <- -c(1,0,0)

      # continuous curvature
      M3[i,i,] <- +c(0,2,6*h)
      M3[i,1+i%%p,] <- -c(0,2,0)
    }
    dim(M1) <- dim(M2) <- dim(M3) <- c(p,p*3)
    b <- yn - y
    Jb <- Jn - J

    M <- rbind(M1,M2,M3)
    rm(M1,M2,M3)

    M <- solve(M)
    M <- M[,1:length(b)] # non-zero b

    Q <- M %*% b
    grad <- M %*% Jb

    dim(Q) <- c(p,3)
    Q <- cbind(y,Q) # [p,4]

    # THIS IS ALL WRONG
    dim(grad) <- c(p,3,p0)
    grad <- aperm(grad,c(3,1,2)) # [p0,p,3]
    dim(grad) <- c(p0*p,3)
    grad <- cbind(c(t(J)),grad) # [p0*p,4]
    dim(grad) <- c(p0,p,4)

    rm(M,b)
  }

  R <- list()

  R$fn <- Vectorize( function(angle)
  {
    # angle oriented from noon==0
    angle <- (angle-pi/2) %% (2*pi)
    # angle from knot
    da <- angle %% h
    # knot index
    i <- 1 + round((angle-da)/h)
    # rate
    c(Q[i,] %*% c(1,da,da^2,da^3))
  } )

  R$grad <- Vectorize( function(angle)
  {
    # angle oriented from noon==0
    angle <- (angle-pi/2) %% (2*pi)
    # angle from knot
    da <- angle %% h
    # knot index
    i <- 1 + round((angle-da)/h)
    # gradient of rate w.r.t timelink.par
    c( grad[,i,] %*% c(1,da,da^2,da^3) ) # [p0]
  } )

  # segment-wise integrals
  INT <- Q %*% c(h,h^2/2,h^3/3,h^4/4)
  INT <- cumsum(INT)
  INT <- c(0,INT[-p])

  int <- function(angle)
  {
    # angle oriented from noon==0
    angle <- (angle-pi/2) %% (2*pi)
    # angle from knot
    da <- angle %% h
    # knot index
    i <- 1 + round((angle-da)/h)
    # integral of rate - time t
    INT[i] + c(Q[i,] %*% c(da,da^2/2,da^3/3,da^4/4))
  }

  # accumulated activity over a whole light|dark period
  R$acc.light <- ( int(pi) + int(5/2*pi)-int(2*pi) ) / pi
  R$acc.dark <- ( int(2*pi)-int(pi) ) / pi

  # accumulated activity over a fraction of a light|dark period
  R$acc <- Vectorize( function(angle)
  {
    if(angle<=pi) #day
    {
      if(angle<=pi/2)
      { M <- int(2*pi+angle)-int(2*pi) }
      else
      { M <- int(5/2*pi)-int(2*pi) + int(angle) }
      M <- M/angle
    }
    else # night
    {
      M <- int(angle)-int(pi)
      M <- M/(angle-pi)
    }
    return(M)
  } )

  # will need derivatives for speed code in future periodic spline mean function
  R$deriv <- Vectorize( function(angle)
  {
    # angle oriented from noon==0
    angle <- (angle-pi/2) %% (2*pi)
    # angle from knot
    da <- angle %% h
    # knot index
    i <- 1 + round((angle-da)/h)
    # gradient of rate w.r.t timelink.par
    c( grad[,i,] %*% c(0,1,2*da,3*da^2) )
  })

  # MSD
  #
  #

  return(R)
}

spline.timelink.parinfo <- function(CTMM)
{
  p0 <- length(CTMM$timelink.par)

  R <- list()
  R$lower <- 0
  R$upper <- p0+1
  R$parscale <- 1
  return(R)
}

# par reset function for fitting
spline.timelink.clean <- function(par)
{
  p0 <- length(par)

  par <- clamp(par,0,Inf)

  # p[0] = p0+1 - sum(par)
  if(sum(par)>p0+1)  { par <- par * (p0+1)/sum(par) }
  # now p[0] = 0

  return(par)
}

spline.timelink.complexify <- function(CTMM)
{
  par <- CTMM$timelink.par
  p <- length(par) + 1 # new length, still doesn't include noon knot
  # new knot locations
  theta <- seq(1/2*pi,5/2*pi,length.out=1+p+1)[-(p+2)]

  # evaluate old spline function at new knots
  fn <- spline.timelink.fn(CTMM)$fn
  par <- fn(theta)
  par <- par/mean(par) # might not be mean-1
  par <- par[-1] # drop noon knot

  CTMM$timelink.par <- par
  return(CTMM)
}

spline.timelink.simplify <- function(CTMM)
{
  par <- CTMM$timelink.par
  p <- length(par)

  if(p==1)
  { CTMM$timelink.par <- NULL }
  else
  {
    p <- p - 1 # new length, does not include noon
    # new knot locations
    theta <- seq(1/2*pi,5/2*pi,length.out=1+p+1)[-(p+2)]

    # evaluate old spline function at new knots
    fn <- spline.timelink.fn(CTMM)$fn
    par <- fn(theta)
    par <- par/mean(par) # might not be mean-1
    par <- par[-1] # drop noon knot

    CTMM$timelink.par <- par
  }

  return(CTMM)
}

spline.timelink.name <- function(object)
{
  n <- length(object$timelink.par)
  if(n)
  { NAME <- paste("spline-timelink",n) }
  else
  { NAME <- NULL }
  return(NAME)
}

spline.timelink.summary <- function(object,level=0.95)
{ return(NULL) }


#########################################
# generic periodic function that is only a function of the (local) time of day
fourier.timelink.fn <- function(CTMM)
{
  par <- CTMM$timelink.par
  p <- length(par)/2
  R <- list()

  # rate
  R$fn <- function(sundial)
  {
    r <- 1
    for(i in 1%:%p)
    {
      if(is.odd(i)) # symmetric              # anti-symmetric
      { r <- r + par[2*i-1]*sin(i*sundial) + par[2*i]*cos(i*sundial) }
      else          # symmetric              # anti-symmetric
      { r <- r + par[2*i-1]*cos(i*sundial) + par[2*i]*sin(i*sundial) }
    }
    return(r)
  }

  # gradient of rate
  R$grad <- function(sundial)
  {
    r <- NULL
    for(i in 1%:%p)
    {
      if(is.odd(i))
      { r <- c(r,sin(i*sundial),cos(i*sundial)) }
      else
      { r <- c(r,cos(i*sundial),sin(i*sundial)) }
    }
    return(r)
  }

  # angular integral
  R$int <- function(sundial)
  {
    r <- 0
    for(i in 1%:%p)
    {
      if(is.odd(i))
      { r <- r - par[2*i-1]/i*cos(i*sundial) + par[2*i]/i*sin(i*sundial) }
      else
      { r <- r + par[2*i-1]/i*sin(i*sundial) - par[2*i]/i*cos(i*sundial) }
    }
    return(r)
  }

  return(R)
}

fourier.timelink.complexify <- function(object)
{
  object$timelink.par <- c(object$timelink.par,0,0)
  return(object)
}

fourier.timelink.simplify <- function(object)
{
  n <- length(object$timelink.par)
  object$timelink.par <- object$timelink.par[1:(n-2)]
  return(object)
}

fourier.timelink.name <- function(object)
{
  n <- length(object$timelink.par)/2
  NAME <- paste("fourier-timelink",n)
  return(NAME)
}

fourier.timelink.summary <- function(object,level=0.95)
{ return(NULL) }


#############################
# generic periodic function that is only a function of light level (symmetric/cosine w.r.t. noon/nadir)
cosine.timelink.fn <- function(CTMM)
{
  par <- CTMM$timelink.par
  p <- length(par)
  R <- list()

  # rate
  R$fn <- function(sundial)
  {
    r <- 1
    for(i in 1%:%p)
    {
      if(is.odd(i))
      { r <- r + par[i]*sin(i*sundial) }
      else
      { r <- r + par[i]*cos(i*sundial) }
    }
    return(r)
  }

  # gradient of rate
  R$grad <- function(sundial)
  {
    r <- NULL
    for(i in 1%:%p)
    {
      if(is.odd(i))
      { r <- c(r,sin(i*sundial)) }
      else
      { r <- c(r,cos(i*sundial)) }
    }
    return(r)
  }

  # angular integral
  R$int <- function(sundial)
  {
    r <- 0
    for(i in 1%:%p)
    {
      if(is.odd(i))
      { r <- r - par[i]/i*cos(i*sundial) }
      else
      { r <- r + par[i]/i*sin(i*sundial) }
    }
    return(r)
  }

  return(R)
}

cosine.timelink.complexify <- function(object)
{
  object$timelink.par <- c(object$timelink.par,0)
  return(object)
}

cosine.timelink.simplify <- switch.timelink.simplify

cosine.timelink.name <- function(object)
{
  n <- length(object$timelink.par)
  NAME <- paste("cosine-timelink",n)
  return(NAME)
}

cosine.timelink.summary <- function(object,level=0.95)
{ return(NULL) }


#################
circadian <- function(CTMM,level=0.95,n=100,...)
{
  FACT <- 1
  # HFACT <- 1 + (FACT-1)/2

  COV <- CTMM$COV
  PAR <- which(grepl("timelink",rownames(COV)))
  COV <- COV[PAR,PAR]

  n <- ceiling(n/2)
  theta <- c(seq(0,pi,length.out=n),seq(pi,2*pi,length.out=n))
  DAY <- 1:n
  NIT <- n + 1:n
  n <- 2*n

  R <- timelink.fn(CTMM)
  GRAD <- sapply(theta,R$grad)
  dim(GRAD) <- c(length(PAR),n)
  R <- R$fn(theta)
  VAR <- sapply(1:n,function(i){GRAD[,i] %*% COV %*% GRAD[,i]})

  CI <- sapply(1:n,function(i){lognorm.ci(R[i],VAR[i],level=level)})

  MAX <- max(CI[3,])
  xlim <- c(-MAX,MAX) * FACT
  ylim <- c(-MAX,MAX) * FACT

  plot(0,0,xlim=xlim,ylim=ylim,xlab=NA,ylab=NA,asp=1,col=c(0,0,0,0),...)
  # lim <- par("usr")
  # xlim <- lim[1:2]
  # ylim <- lim[3:4]
  # X <- c(xlim,rev(xlim))
  # Y <- c(0,0,ylim[2],ylim[2])
  # polygon(X,Y,border=NA,col=rgb(1,1,0,0.1)) # day shade
  # Y <- c(0,0,ylim[1],ylim[1])
  # polygon(X,Y,border=NA,col=rgb(0,0,1,0.1)) # night shade

  COS <- cos(theta)
  SIN <- sin(theta)

  # point estimate - 100% circle - colored
  X <- c( (CI[2,]*COS)[DAY] , rev(COS[DAY]) )
  Y <- c( (CI[2,]*SIN)[DAY] , rev(SIN[DAY]) )
  graphics::polygon(X,Y,border=NA,col=grDevices::rgb(1,1,0,1/2),...)
  X <- c( (CI[2,]*COS)[NIT] , rev(COS[NIT]) )
  Y <- c( (CI[2,]*SIN)[NIT] , rev(SIN[NIT]) )
  graphics::polygon(X,Y,border=NA,col=grDevices::rgb(0,0,1,1/2),...)

  # CIs
  X <- c(CI[1,]*COS,rev(CI[3,]*COS))
  Y <- c(CI[1,]*SIN,rev(CI[3,]*SIN))
  graphics::polygon(X,Y,border=NA,col=grDevices::rgb(0,0,0,1/4),...)

  # point estimates
  X <- CI[2,]*COS
  Y <- CI[2,]*SIN
  graphics::lines(X[DAY],Y[DAY],...)
  graphics::lines(X[NIT],Y[NIT],...)

  # normal circle
  #shape::plotcircle(r=1,lcol=rgb(0,0,0,1/4),lwd=1)
  #shape::plotcircle(r=1/2,lcol=rgb(0,0,0,1/4),lwd=1)
  #shape::plotcircle(r=3/2,lcol=rgb(0,0,0,1/4),lwd=1)
  graphics::abline(h=0,col=grDevices::rgb(0,0,0,1/4))
  graphics::abline(v=0,col=grDevices::rgb(0,0,0,1/4))

  # sun "\u2600"
  # text(0,FACT*MAX,labels="\u2600",cex=2)

  # star "\u2605"
  # text(0,-FACT*MAX,labels="\u2605")
}
