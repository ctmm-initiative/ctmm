POSITIVE.PARAMETERS <- c("major","minor","tau position","tau velocity","tau","omega","error")

# clean up parameter arrays
clean.parameters <- function(par,linear.cov=FALSE)
{
  NAMES <- names(par)

  for(P in POSITIVE.PARAMETERS) { if(P %in% NAMES) { par[P] <- clamp(par[P],0,Inf) } }

  # covariance parameters
  if("minor" %in% NAMES)
  {
    if(!linear.cov)
    {
      # swap major and minor axis
      if(par["minor"]>1)
      {
        par["angle"] <- par["angle"] - pi/2

        # not profiled
        if('major' %in% NAMES) { par['major'] <- par['major'] * par['minor'] }
        # else 'major' will be profiled

        par['minor'] <- 1/par['minor']
      }

      # wrap angle back
      par["angle"] <- (((par["angle"]/pi+1/2) %% 1) - 1/2)*pi
    }
    else
    {
      # sigma[x,x] sigma[y,y]
      MAX <- sqrt(prod(par[c('major','minor')]))
      # sigma[x,y]
      par['angle'] <- clamp(par['angle'],-MAX,MAX)
    }
  }

  P <- c("tau position","tau velocity")
  if(all(P %in% NAMES)) { par[P] <- sort(par[P],decreasing=TRUE) }

  return(par)
}


# identify autocovariance parameters in ctmm object
# if profile==TRUE, some parameters can be solved exactly and so aren't identified
# if linear==TRUE, only return linear non-problematic parameters
# if linear.cov==TRUE, use linear covariance paramters for pREML
# STRUCT determines the model structure incase some estimated parameters in CTMM are zero
id.parameters <- function(CTMM,profile=TRUE,linear=FALSE,linear.cov=FALSE,UERE=FALSE,dt=0,df=0,dz=0,STRUCT=CTMM)
{
  # identify and name autocovariance parameters
  NAMES <- NULL

  if(STRUCT$sigma@par['major'] || length(STRUCT$tau))
  {
    sigma <- attr(CTMM$sigma,"par")
    if(!profile || (STRUCT$error && UERE>=3)) # must fit all 1-3 covariance parameters
    { if(STRUCT$isotropic) { NAMES <- c(NAMES,"major") } else { NAMES <- c(NAMES,names(sigma)) } }
    else if(STRUCT$error || STRUCT$circle) # must fit shape, but scale/area/var can be profiled for free
    { if(!STRUCT$isotropic) { NAMES <- c(NAMES,names(sigma[2:3])) } } # the latter are ratios to the former
  }

  if(!linear) # nonlinear autocorrelation parameters
  {
    TAU <- CTMM$tau
    if(!STRUCT$range) { TAU <- TAU[-1] }
    if(length(TAU))
    {
      if(length(TAU)==2 && TAU[1]==TAU[2]) { TAU <- TAU[1] ; NAMES <- c(NAMES,"tau") } # identical timescales
      else { NAMES <- c(NAMES,paste("tau",names(TAU))) } # distinct timescales
    }

    if(STRUCT$omega) { NAMES <- c(NAMES,"omega") }

    if(STRUCT$circle) { NAMES <- c(NAMES,"circle") }
  }

  if((STRUCT$error) && UERE<3) { NAMES <- c(NAMES,"error") }

  # setup parameter information
  parscale <- NULL
  lower <- NULL
  upper <- NULL
  period <- NULL

  SIGMIN <- dz^2*ifelse(CTMM$range,1,df)

  if("major" %in% NAMES)
  {
    parscale <- c(parscale,max(sigma['major'],SIGMIN))
    lower <- c(lower,0)
    upper <- c(upper,Inf)
    period <- c(period,FALSE)
  }

  # minor and angle
  if("minor" %in% NAMES)
  {
    if(!linear.cov)
    {
      parscale <- c(parscale,1,pi/2)
      period <- c(period,FALSE,pi)
    }
    else
    {
      parscale <- c(parscale,pmax(abs(sigma[c('minor','angle')]),SIGMIN))
      period <- c(period,FALSE,FALSE)
    }
    lower <- c(lower,0,-Inf)
    upper <- c(upper,Inf,Inf) # could make 1 for minor/major
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

    if(STRUCT$omega)
    {
      parscale <- c(parscale,max(CTMM$omega,df))
      lower <- c(lower,0)
      upper <- c(upper,Inf)
      period <- c(period,FALSE)
    }

    if(STRUCT$circle)
    {
      parscale <- c(parscale,max(abs(CTMM$circle),df))
      lower <- c(lower,-Inf)
      upper <- c(upper,Inf)
      period <- c(period,FALSE)
    }
  }

  if(STRUCT$error && UERE<3)
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
# linear.cov uses linear representation of covariance matrix
get.parameters <- function(CTMM,NAMES,linear.cov=FALSE)
{
  par <- numeric(length(NAMES))
  names(par) <- NAMES

  # can't loop this easily because R collapses NULL values
  getter <- function(NAME,VALUE=CTMM[[NAME]]) { if(NAME %in% NAMES) { par[NAME] <<- VALUE } }

  sigma <- CTMM$sigma
  if(!linear.cov)
  {
    sigma <- attr(sigma,"par")
    getter("major",sigma[1])
    getter("minor",sigma[2])
    getter("angle",sigma[3])
  }
  else
  {
    getter("major",sigma[1]) # sigma[x,x]
    getter("minor",sigma[4]) # sigma[y,y]
    getter("angle",sigma[2]) # sigma[x,y]
  }

  tau <- CTMM$tau
  getter("tau position",if(length(tau)>0) { tau[1] } else { 0 })
  getter("tau velocity",if(length(tau)>1) { tau[2] } else { 0 })
  getter("tau",if(length(tau)>0) { min(tau) } else { 0 }) # identical timescales
  getter("omega")

  getter("circle")
  getter("error")

  return(par)
}


# write autocovariance parameter array into a CTMM object
set.parameters <- function(CTMM,par,linear.cov=FALSE)
{
  NAMES <- names(par)

  sigma <- CTMM$sigma
  if(!linear.cov)
  {
    sigma <- sigma@par
    NAME <- "major"; if(NAME %in% NAMES) { sigma[NAME] <- par[NAME] }
    NAME <- "minor"; if(NAME %in% NAMES) { sigma[NAME] <- par[NAME] }
    NAME <- "angle"; if(NAME %in% NAMES) { sigma[NAME] <- par[NAME] }
  }
  else
  {
    sigma <- methods::getDataPart(sigma)
    NAME <- "major"; if(NAME %in% NAMES) { diag(sigma) <- par[c(NAME,NAME)] }
    NAME <- "minor"; if(NAME %in% NAMES) { sigma[4] <- par[NAME] }
    NAME <- "angle"; if(NAME %in% NAMES) { sigma[2:3] <- par[c(NAME,NAME)] }
  }
  CTMM$sigma <- covm(sigma)

  NAME <- "tau position"; if(NAME %in% NAMES) { CTMM$tau[1] <- par[NAME] ; names(CTMM$tau)[1] <- 'position' }
  NAME <- "tau velocity"; if(NAME %in% NAMES) { CTMM$tau[2] <- par[NAME] ; names(CTMM$tau)[2] <- 'velocity' }
  NAME <- "tau"; if(NAME %in% NAMES) { CTMM$tau <- c(1,1)*par[NAME] ; names(CTMM$tau) <- c('position','velocity') } # identical timescales
  NAME <- "omega"; if(NAME %in% NAMES) { CTMM$omega <- par[NAME] }

  NAME <- "circle"; if(NAME %in% NAMES) { CTMM[[NAME]] <- par[NAME] }
  NAME <- "error"; if(NAME %in% NAMES) { CTMM[[NAME]] <- par[NAME] }

  return(CTMM)
}
