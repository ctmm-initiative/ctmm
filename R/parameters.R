# clean up parameter arrays
clean.parameters <- function(par)
{
  NAMES <- names(par)

  # enforce positivity
  if("area" %in% NAMES) { par["area"] <- abs(par["area"]) }

  if("eccentricity" %in% NAMES)
  {
    # swap major and minor axis
    if(par["eccentricity"]<0)
    {
      par["eccentricity"] <- abs(par["eccentricity"])
      par["angle"] <- par["angle"] - pi/2
    }

    # wrap angle back
    par["angle"] <- (((par["angle"]/pi+1/2) %% 1) - 1/2)*pi
  }

  # enforce positivity
  TAUS <- startsWith(NAMES,"tau")
  if(length(TAUS)) { par[TAUS] <- abs(par[TAUS]) }

  # enforce positivity
  if("error" %in% NAMES) { par["error"] <- abs(par["error"]) }

  return(par)
}


# identify autocovariance parameters in ctmm object
# if profile=TRUE, some parameters can be solved exactly and so aren't identified
# if linear=TRUE, only return linear non-problematic parameters
id.parameters <- function(CTMM,profile=TRUE,linear=FALSE,UERE=FALSE,dt=0,df=0,dz=10)
{
  # identify and name autocovariance parameters
  NAMES <- NULL

  sigma <- attr(CTMM$sigma,"par")
  if(!profile || (CTMM$error && UERE>=3)) # must fit all 1-3 covariance parameters
  { if(CTMM$isotropic) { NAMES <- c(NAMES,"area") } else { NAMES <- c(NAMES,names(sigma)) } }
  else if(CTMM$error || CTMM$circle) # must fit shape, but scale/area/var can be profiled for free
  { if(!CTMM$isotropic) { NAMES <- c(NAMES,names(sigma[2:3])) } }

  if(!linear) # nonlinear autocorrelation parameters
  {
    TAU <- CTMM$tau
    if(!CTMM$range) { TAU <- TAU[-1] }
    if(length(TAU)) { NAMES <- c(NAMES,paste("tau",names(TAU))) }

    if(CTMM$circle) { NAMES <- c(NAMES,"circle") }
  }

  if(CTMM$error && UERE<3) { NAMES <- c(NAMES,"error") }

  # setup parameter information
  parscale <- NULL
  lower <- NULL
  upper <- NULL
  period <- NULL

  if("area" %in% NAMES)
  {
    parscale <- c(parscale,sigma[1])
    lower <- c(lower,0)
    upper <- c(upper,Inf)
    period <- c(period,FALSE)
  }

  # eccentricity and angle
  if("eccentricity" %in% NAMES)
  {
    parscale <- c(parscale,log(2),pi/2)
    lower <- c(lower,-Inf,-Inf)
    upper <- c(upper,Inf,Inf)
    period <- c(period,FALSE,pi)
  }

  if(!linear) # nonlinear autocorrelation parameters
  {
    if(length(TAU))
    {
      parscale <- c(parscale,pmax(TAU,dt))
      lower <- c(lower,rep(0,length(TAU)))
      upper <- c(upper,rep(Inf,length(TAU)))
      period <- c(period,rep(FALSE,length(TAU)))
    }

    if(CTMM$circle)
    {
      parscale <- c(parscale,max(abs(CTMM$circle),df))
      lower <- c(lower,-Inf)
      upper <- c(upper,Inf)
      period <- c(period,FALSE)
    }
  }

  if(CTMM$error && UERE<3)
  {
    parscale <- c(parscale,max(dz,CTMM$error)) # minimum of dz error parscale
    lower <- c(lower,0)
    upper <- c(upper,Inf)
    period <- c(period,FALSE)
  }

  if(length(parscale))
  {
    names(parscale) <- NAMES
    names(lower) <- NAMES
    names(upper) <- NAMES
    names(period) <- NAMES
  }

  return(list(NAMES=NAMES,parscale=parscale,lower=lower,upper=upper,period=period))
}


# pull out array of named autocovariance parameters
get.parameters <- function(CTMM,NAMES)
{
  par <- numeric(length(NAMES))
  names(par) <- NAMES

  # can't loop this easily because R collapses NULL values
  getter <- function(NAME,VALUE=CTMM[[NAME]]) { if(NAME %in% NAMES) { par[NAME] <<- VALUE } }

  sigma <- attr(CTMM$sigma,"par")
  getter("area",sigma[1])
  getter("eccentricity",sigma[2])
  getter("angle",sigma[3])

  tau <- CTMM$tau
  getter("tau position",tau[1])
  getter("tau velocity",tau[2])

  getter("circle")
  getter("error")

  return(par)
}


# write autocovariance parameter array into a CTMM object
set.parameters <- function(CTMM,par)
{
  NAMES <- names(par)

  sigma <- CTMM$sigma@par
  NAME <- "area"; if(NAME %in% NAMES) { sigma[NAME] <- par[NAME] }
  NAME <- "eccentricity"; if(NAME %in% NAMES) { sigma[NAME] <- par[NAME] }
  NAME <- "angle"; if(NAME %in% NAMES) { sigma[NAME] <- par[NAME] }
  CTMM$sigma <- covm(sigma)

  NAME <- "tau position"; if(NAME %in% NAMES) { CTMM$tau[1] <- par[NAME] ; names(CTMM$tau)[1] <- 'position' }
  NAME <- "tau velocity"; if(NAME %in% NAMES) { CTMM$tau[2] <- par[NAME] ; names(CTMM$tau)[2] <- 'velocity' }

  NAME <- "circle"; if(NAME %in% NAMES) { CTMM[[NAME]] <- par[NAME] }
  NAME <- "error"; if(NAME %in% NAMES) { CTMM[[NAME]] <- par[NAME] }

  return(CTMM)
}
