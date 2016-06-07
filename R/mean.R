# drift/mean function models
new.drift <- methods::setClass("drift", representation("function",init="function",scale="function",svf="function",refine="function",name="function",speed="function",summary="function"))

#################################
# stationary mean/drift function
#################################
stationary.drift <- function(t,CTMM) { cbind( array(1,length(t)) ) }

# guess some parameters and check the model parameter sanity
stationary.init <- function(data,CTMM)
{
  z <- get.telemetry(data,CTMM$axes)
  
  # weights from errors
  error <- get.error(data,CTMM)
  if(CTMM$error) { w <- 1/error }
  else { w <- rep(1,length(data$t)) }
  # normalize weights
  w <- w/sum(w)
  
  CTMM$mu <- c(w %*% z)

  z <- t(t(z)-CTMM$mu)
  
  CTMM$sigma <- t(z) %*% z
  
  # remove error from variability
  if(CTMM$error) { CTMM$sigma <- CTMM$sigma - (sum(error) * diag(length(CTMM$axes))) }
  
  n <- length(data$t)
  CTMM$sigma <- CTMM$sigma / (n-1)  
  
  CTMM$sigma <- covm(CTMM$sigma,isotropic=CTMM$isotropic,axes=CTMM$axes)
  
  return(CTMM)
}

# how do we rescale time
stationary.scale <- function(CTMM,time) { CTMM }

# svf of mean function
stationary.svf <- function(CTMM)
{
  # default
  EST <- function(t) { 0 }
  VAR <- function(t) { 0 }
  
  return(list(EST=EST,VAR=VAR))
}

# increase number of linear model parameters
stationary.refine <- function(CTMM) { list() }

# name of mean function
stationary.name <- function(CTMM) { NULL }

# mean square speed function: point estimate and variance
stationary.speed <- function(CTMM) { list(EST=0,VAR=0) }

# place to put optional summary information
stationary.summary <- function(CTMM,level,level.UD) { NULL }

# combine this all together for convenience
stationary <- new.drift(stationary.drift,init=stationary.init,scale=stationary.scale,svf=stationary.svf,refine=stationary.refine,name=stationary.name,speed=stationary.speed,summary=stationary.summary)

############################
# Periodic drift function
############################
# periodic mean/drift function
periodic.drift <- function(t,CTMM)
{
  harmonic <- CTMM$harmonic
  period <- CTMM$period
  
  # constant term
  U <- cbind( stationary.drift(t,CTMM) )
  
  omega <- periodic.omega(CTMM)
  if(length(omega))
  {
    omega <- t %o% omega
    U <- cbind( U , cos(omega) , sin(omega) )
  }
  
  return(U)
}

# guess parameters
periodic.init <- function(data,CTMM)
{
  # default period of 1 day
  if(is.null(CTMM$period)) { CTMM$period <- 24*60^2 } 

  # default harmonics is none
  if(is.null(CTMM$harmonic)) { CTMM$harmonic <- numeric(length(CTMM$period)) }

  if(is.null(CTMM$mu)) { CTMM <- stationary.init(data,CTMM) }

  return(CTMM)
}

# how do we rescale time
periodic.scale <- function(CTMM,time)
{
  CTMM$period <- CTMM$period / time
  
  return(CTMM)
}

# SVF of mean function
periodic.svf <- function(CTMM)
{
  if(sum(CTMM$harmonic)==0) { return( stationary.svf(CTMM) ) }
  
  STUFF <- periodic.stuff(CTMM)
  omega <- STUFF$omega
  A <- STUFF$A
  COV <- STUFF$COV
  
  EST <- function(t) { sum( 1/4 * A^2 * (1-cos(omega*t) )) }
  
  VAR <- function(t)
  {
    grad <- 1/2 * A * (1-cos(omega*t))
    return(grad %*% COV %*% grad)
  }
  
  return(list(EST=EST,VAR=VAR))
}


# name of mean function
periodic.name <- function(CTMM)
{
  NAME <- CTMM$harmonic
  NAME <- paste(NAME,collapse=" ")
  NAME <- paste("harmonic",NAME)
  return(NAME)
}

# calculate deterministic mean square speed and its variance
periodic.speed <- function(CTMM)
{
  if(sum(CTMM$harmonic)==0) { return(stationary.speed(CTMM)) }
  
  STUFF <- periodic.stuff(CTMM)
  omega <- STUFF$omega
  A <- STUFF$A
  COV <- STUFF$COV
  
  if(is.null(omega)) { return( stationary.speed(CTMM) ) }

  EST <- sum(omega^2*A^2)/2
  grad <- omega^2*A
  VAR <- grad %*% COV %*% grad

  return(list(EST=EST,VAR=VAR))
}

# calculate rotational indices
periodic.summary <- function(CTMM,level,level.UD)
{
  alpha <- 1-level
  SUM <- stationary.summary(CTMM,level,level.UD)
  
  if(sum(CTMM$harmonic)==0) { return(SUM) }
  
  STUFF <- periodic.stuff(CTMM)
  omega <- STUFF$omega
  A <- STUFF$A
  COV <- STUFF$COV

  # ROTATIONAL VARIANCE INDEX
  ROT <- sum(A^2)/2 # sine and cosine average 1/2
  area <- CTMM$sigma@par["area"]
  ecc <- CTMM$sigma@par["eccentricity"]
  RAN <- 2 * area * cosh(ecc/2)
  MLE <- ROT/(ROT+RAN)
  
  GRAD <- A
  COV.ROT <- GRAD %*% COV %*% GRAD

  GRAD <- RAN/area
  PARS <- "area"
  if(!CTMM$isotropic)
  {
    GRAD <- c(GRAD,RAN*tanh(ecc/2)/2)
    PARS <- c(PARS,"eccentricity")
  }
  COV.RAN <- CTMM$COV[PARS,PARS]
  COV.RAN <- GRAD %*% COV.RAN %*% GRAD

  GRAD <- c(RAN,-ROT)/(ROT+RAN)^2
  VAR <- sum(GRAD^2 * c(COV.ROT,COV.RAN))

  CI <- ci.tau(MLE,VAR,alpha=alpha,min=0,max=1)
  CI <- 100*sqrt(rbind(CI))
  rownames(CI) <- "rotation/deviation %"

  SUM <- rbind(SUM,CI)
  
  # ROTATIONAL SPEED INDEX
  # CIs are too big to be useful
  if(length(CTMM$tau)>1)
  {
    ROT <- sum((omega*A)^2)/2 # sine and cosine average 1/2
    Omega <- if(CTMM$circle) { 2*pi/CTMM$circle } else { 0 }
    RAN <- 2*area*cosh(ecc/2) * ( 1/prod(CTMM$tau) + Omega^2 )
    MLE <- ROT/(ROT+RAN)

    GRAD <- omega^2*A
    COV.ROT <- GRAD %*% COV %*% GRAD

    GRAD <- c(RAN/area,-2*area*cosh(ecc/2)/prod(CTMM$tau)/CTMM$tau)
    PARS <- c("area","tau position","tau velocity")
    if(!CTMM$isotropic)
    {
      GRAD <- c(GRAD,RAN*tanh(ecc/2)/2)
      PARS <- c(PARS,"eccentricity")
    }
    if(CTMM$circle)
    {
      GRAD <- c(GRAD,-4*area*cosh(ecc/2)*Omega^2/CTMM$circle)
      PARS <- c(PARS,"circle")
    }
    COV.RAN <- CTMM$COV[PARS,PARS]
    COV.RAN <- GRAD %*% COV.RAN %*% GRAD

    GRAD <- c(RAN,-ROT)/(ROT+RAN)^2
    VAR <- sum(GRAD^2 * c(COV.ROT,COV.RAN))

    CI <- ci.tau(MLE,VAR,alpha=alpha,min=0,max=1)
    CI <- 100*sqrt(rbind(CI))
    rownames(CI) <- "rotation/speed %"

    SUM <- rbind(SUM,CI)
  }
          
  return(SUM)
}

# increase number of harmonics in model
periodic.refine <- function(CTMM)
{
  period <- CTMM$period
  harmonic <- CTMM$harmonic
  
  GUESS <- list()
  
  # dt <- stats::median(diff(data$t))
  # P <- attr(CTMM$mean,"par")$P
  # # max harmonics for Nyquist frequency
  # kmax <- P/(2*dt)
  
  for(i in 1:length(period))
  {
    GUESS[[length(GUESS)+1]] <- CTMM
    GUESS[[length(GUESS)]]$harmonic[i] <- harmonic[i] + 1
  }
  
  return(GUESS)
}

# combine this all together for convenience
periodic <- new.drift(periodic.drift,init=periodic.init,scale=periodic.scale,svf=periodic.svf,refine=periodic.refine,name=periodic.name,speed=periodic.speed,summary=periodic.summary)

# list of non-zero frequencies of model
periodic.omega <- function(CTMM)
{
  harmonic <- CTMM$harmonic
  period <- CTMM$period
  
  omega <- NULL
  for(i in 1:length(harmonic))
  {
    if(harmonic[i] > 0)
    {
      w <- (2*pi/period[i]) * (1:harmonic[i])
      omega <- c( omega , w )
    }
  }
  
  return(omega)
}

# useful info from periodic mean function parameters
periodic.stuff <- function(CTMM)
{
  harmonic <- CTMM$harmonic
  period <- CTMM$period
  omega <- periodic.omega(CTMM)
  
  # amplitudes and covariance
  A <- CTMM$mu
  COV <- CTMM$COV.mu

  # total number of harmonics
  k <- length(omega)

  # default uncertainty of none
  if(is.null(COV)) { COV <- 0 }
  # default structure
  COV <- array(COV,c(2,2*k+1,2*k+1,2))
  
  # delete stationary coefficients
  A <- A[-1,,drop=FALSE]
  COV <- COV[,-1,,,drop=FALSE]
  COV <- COV[,,-1,,drop=FALSE]
  
  # structure omega like A
  omega <- c(omega,omega)
  omega <- cbind(omega,omega)
  # this is now (2k)*2
  
  # flatten block-vectors and block-matrices
  A <- array(A,(2*k)*2)
  omega <- array(omega,(2*k)*2)
  COV <- aperm(COV,c(2,1,3,4))
  COV <- array(COV,c((2*k)*2,(2*k)*2))
  
  return(list(A=A,COV=COV,omega=omega))
}

#################################
# continuous uniform spline mean functions
#################################
uspline.drift <- function(t,CTMM)
{
  degree <- CTMM$degree
  knot <- CTMM$knot

  STUFF <- uspline.stuff(CTMM)
  tknot <- STUFF$tknot
  dt <- STUFF$dt
  
  # midpoint / no real grid / degenerate case
  if(knot==1)
  {
    U <- stationary.drift(t,CTMM)
    if(degree==1) { return( U ) }
    
    U <- cbind( U , (t-tknot) )
    return(U)
  }
  
  U <- array(0,c(length(t),knot,degree))
  
  # current starting knot
  k <- 1
  for(i in 1:length(t))
  {
    s <- t[i]
    # iterate knot number k until we are between the appropriate knots
    while(s>tknot[k+1] && k<knot) { k <- k+1 }
    
    # relative distance betwee knots
    s <- (s - tknot[k])/dt
    
    if(degree==1) # linear interpolant
    {
      U[i,k,1] <- 1-s
      U[i,k+1,1] <- s
    }
    else # cubic Hermite interpolant
    {
      U[i,k,1] <- (2*s^3 - 3*s^2 + 1)
      U[i,k,2] <- (s^3 - 2*s^2 + s)*dt
      U[i,k+1,1] <- (-2*s^3 + 3*s^2)
      U[i,k+1,2] <- (s^3 - s^2)*dt
    }
  }

  U <- array(U,c(length(t),knot*degree))
  
  return(U)  
}

# initialize default parameters
uspline.init <- function(data,CTMM)
{
  # degree of continuity: 1,2 - location,velocity
  if(is.null(CTMM$degree)) { CTMM$degree <- 1 }
  # default number of knots 
  if(is.null(CTMM$knot)) { CTMM$knot <- 1 }
  # default domain of spline grid
  if(is.null(CTMM$domain)) { CTMM$domain <- c(data$t[1],last(data$t)) }
  
  if(is.null(CTMM$mu)) { CTMM <- stationary.init(CTMM) }
  
  return(CTMM)  
}

# scale times
uspline.scale <- function(CTMM,time)
{
  knot <- CTMM$knot
  degree <- CTMM$degree
  
  if(CTMM$degree>1)
  {
    scale <- array(1,c(knot,degree))
    scale[,2] <- time
    scale <- array(scale,knot*degree)
    
    CTMM$mu[,1] <- CTMM$mu[,1] * scale
    CTMM$mu[,2] <- CTMM$mu[,2] * scale
  }
  
  return(CTMM)
}

# svf

# refine
uspline.refine <- function(CTMM)
{
  CTMM$knot <- CTMM$knot + 1
  return(list(CTMM))
}


# name of mean function
uspline.name <- function(CTMM)
{
  NAME <- paste("degree",CTMM$degree,"knot",CTMM$knot)
  return(NAME)
}

# ms speed of spline function
uspline.speed <- function(CTMM)
{
  degree <- CTMM$degree
  knot <- CTMM$knot
  
  STUFF <- uspline.stuff()
  dt <- STUFF$dt
  
  if(degree==1 && knot==1)
  {
    EST <- 0
    VAR <- 0
  }
  else if(degree==1)
  {
    EST <- (diff(CTMM$mu[,1])^2 + diff(CTMM$mu[,2])^2)/dt^2/(knot-1)
    VAR <- NA
  }
  else
  {
    EST <- NA
    VAR <- NA
  }
  
  return(list(EST=EST,VAR=VAR))
}

# summary

# convenience function
uspline <- new.drift(uspline.drift,init=uspline.init,scale=uspline.scale,refine=uspline.refine,name=uspline.name)

# calculate spline grid
uspline.stuff <- function(CTMM)
{
  knot <- CTMM$knot
  domain <- CTMM$domain
  
  if(knot>1)
  {
    # location of knots
    tknot <- seq(domain[1],domain[2],length.out=knot)
    # grid spacing
    dt <- diff(domain)/(knot-1)
  }
  else
  {
    tknot <- mean(domain)
    dt <- diff(domain)
  }
    
  return(list(tknot=tknot,dt=dt))
}