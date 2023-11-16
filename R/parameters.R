POSITIVE.PARAMETERS <- c("major","minor","tau position","tau velocity","tau","omega","error")


# identify autocovariance parameters in ctmm object
# if profile==TRUE, some parameters can be solved exactly and so aren't identified -- in particular, only relative variances matter
# if linear==TRUE, only return linear non-problematic parameters
# if linear.cov==TRUE, use linear covariance parameters for pREML
# STRUCT determines the model structure incase some estimated parameters in CTMM are zero
id.parameters <- function(CTMM,profile=TRUE,linear=FALSE,linear.cov=FALSE,UERE.FIT=NA,dt=0,df=0,dz=0,STRUCT=CTMM,fit=FALSE)
{
  AXES <- length(CTMM$axes)
  # identify and name autocovariance parameters
  NAMES <- NULL

  states <- get.states(CTMM)
  if(length(states))
  {
    R <- list()
    for(s in states)
    { R[[s]] <- id.parameters(CTMM[[s]],profile=profile,linear=linear,linear.cov=linear.cov,UERE.FIT=UERE.FIT,dt=dt,dz=dz,STRUCT=STRUCT) }

    NAMES <- list()
    parscale <- lower <- upper <- period <- NULL
    for(s in states)
    {
      NAMES[[s]] <- R[[s]]$NAMES
      parscale <- c(parscale,R[[s]]$parscale)
      lower <- c(lower,R[[s]]$lower)
      upper <- c(upper,R[[s]]$upper)
      period <- c(period,R[[s]]$period)
    }

    return(list(NAMES=NAMES,parscale=parscale,lower=lower,upper=upper,period=period))
  }

  if(STRUCT$sigma@par['major'] || length(STRUCT$tau))
  {
    sigma <- attr(CTMM$sigma,"par")
    if(!profile || any(STRUCT$error & !UERE.FIT)) # must fit all 1-3 covariance parameters - not profiling or fixed errors
    { if(STRUCT$isotropic) { NAMES <- c(NAMES,"major") } else { NAMES <- c(NAMES,names(sigma)) } }
    else if(any(UERE.FIT) || STRUCT$circle) # must fit shape, but overall scale/var can be profiled for free - circulation or relative errors
    { if(!STRUCT$isotropic) { NAMES <- c(NAMES,names(sigma[2:3])) } }
  }

  if(any(UERE.FIT)) { NAMES <- c(NAMES, paste("error",names(UERE.FIT)[UERE.FIT>0]) ) }

  if(!linear) # nonlinear mean and autocovariance parameters
  {
    DRIFT <- drift.pars(CTMM,all=TRUE,fit=fit)
    if(length(DRIFT$pars)) { NAMES <- c(NAMES,DRIFT$NAMES) }

    TIMELINK <- length(CTMM$timelink.par)
    if(TIMELINK) { NAMES <- c(NAMES,paste0("timelink-",1:TIMELINK)) }

    TAU <- STRUCT$tau # for structure here
    if(!STRUCT$range) { TAU <- TAU[-1] }
    if(length(TAU))
    {
      if(length(TAU)==2 && TAU[1]==TAU[2]) { TAU <- TAU[1] ; NAMES <- c(NAMES,"tau") } # identical timescales
      else { NAMES <- c(NAMES,paste("tau",names(TAU))) } # distinct timescales
    }
    TAU <- CTMM$tau # for parscale later
    if(!STRUCT$range) { TAU <- TAU[-1] }
    if("tau" %in% NAMES) { TAU <- mean(TAU) } # just the one

    if(STRUCT$omega) { NAMES <- c(NAMES,"omega") }

    if(STRUCT$circle) { NAMES <- c(NAMES,"circle") }
  }

  ## ORDER MATTERS BETWEEN THESE TWO CODE BLOCKS!

  # setup parameter information
  parscale <- NULL
  lower <- NULL
  upper <- NULL
  period <- NULL
  # reflect <- NULL

  SIGMIN <- dz^2*ifelse(CTMM$range,1,df)

  MAX <- eigenvalues.covm(CTMM$sigma)[1]
  MAX <- max(MAX,SIGMIN)

  if(!linear.cov) # nonlinear sigma: major, major/minor, angle
  {
    if("major" %in% NAMES)
    {
      parscale <- c(parscale,MAX)
      lower <- c(lower,0)
      upper <- c(upper,Inf)
      period <- c(period,FALSE)
    }

    # minor and angle
    if("minor" %in% NAMES)
    {
      parscale <- c(parscale,MAX,pi/2)
      lower <- c(lower,0,-Inf)
      upper <- c(upper,Inf,Inf)
      period <- c(period,FALSE,pi)
    }

    # RMS UERE - error parameters
    if(any(UERE.FIT))
    {
      parscale <- c(parscale, pmax(CTMM$error[UERE.FIT],dz) ) # minimum of dz error parscale
      lower <- c(lower, rep(0,sum(UERE.FIT)) )
      upper <- c(upper, rep(Inf,sum(UERE.FIT)) )
      period <- c(period, rep(FALSE,sum(UERE.FIT)) )
    }
  }
  else # linear sigma: xx, yy, xy
  {
    if("major" %in% NAMES) # xx || yy
    {
      parscale <- c(parscale,MAX)
      lower <- c(lower,0)
      upper <- c(upper,Inf)
      period <- c(period,FALSE)
    }

    if("minor" %in% NAMES) # yy & xy
    {
      parscale <- c(parscale,MAX,MAX)
      lower <- c(lower,0,-Inf)
      upper <- c(upper,Inf,Inf)
      period <- c(period,FALSE,FALSE)
    }

    # MS UERE - error parameters
    if(any(UERE.FIT))
    {
      parscale <- c(parscale, pmax(CTMM$error[UERE.FIT],dz)^2 ) # minimum of dz error parscale
      lower <- c(lower, rep(0,sum(UERE.FIT)) )
      upper <- c(upper, rep(Inf,sum(UERE.FIT)) )
      period <- c(period, rep(FALSE,sum(UERE.FIT)) )
    }
  } # end linear sigma

  if(!linear) # nonlinear autocorrelation parameters
  {
    if(length(DRIFT$pars))
    {
      parscale <- c(parscale,DRIFT$parscale)
      lower <- c(lower,DRIFT$lower)
      upper <- c(upper,DRIFT$upper)
      period <- c(period,rep(FALSE,length(DRIFT$pars)))
    }

    if(TIMELINK)
    {
      TLI <- timelink.parinfo(CTMM)

      parscale <- c(parscale,rep(TLI$parscale,TIMELINK))
      lower <- c(lower,rep(TLI$lower,TIMELINK))
      upper <- c(upper,rep(TLI$upper,TIMELINK))
      period <- c(period,rep(FALSE,TIMELINK))
    }

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

  if(length(parscale))
  {
    names(parscale) <- NAMES
    names(lower) <- NAMES
    names(upper) <- NAMES
    names(period) <- NAMES
  }

  # only relative variance matters - easier optimization (more constrained, well behaved)
  if(profile)
  {
    PARS <- c("major","minor")
    PARS <- PARS[PARS %in% NAMES]

    if(length(PARS))
    {
      parscale[PARS] <- 1
      upper[PARS] <- 1
    }

    # error fraction could be much smaller
    if(any(UERE.FIT))
    {
      PARS <- paste("error",names(UERE.FIT)[UERE.FIT])

      parscale[PARS] <- min(1, pmax(CTMM$error[UERE.FIT],dz) / sqrt(MAX) )
      upper[PARS] <- 1
    }
  } # end profile

  return(list(NAMES=NAMES,parscale=parscale,lower=lower,upper=upper,period=period))
}


# pull out array of named autocovariance parameters
# linear.cov uses linear representation of covariance matrix
get.parameters <- function(CTMM,NAMES,profile=FALSE,linear.cov=FALSE)
{
  par <- numeric(length(NAMES))
  names(par) <- NAMES
  states <- get.states(CTMM)
  if(length(states))
  {
    par <- list()
    for(s in states) { par[[s]] <- get.parameters(CTMM[[s]],NAMES[[s]],profile=profile,linear.cov=linear.cov) }
    names(par) <- states
    return(par)
  }

  # can't loop this easily because R collapses NULL values
  getter <- function(NAME,VALUE=CTMM[[NAME]])
  { if(NAME %in% NAMES) { par[NAME] <<- VALUE } }

  beta <- CTMM$beta
  for(B in names(beta)) { getter(B,beta[B]) }

  sigma <- CTMM$sigma
  if(!linear.cov)
  {
    sigma <- attr(sigma,"par")
    getter("major",sigma[1])
    getter("minor",sigma[2])
    getter("angle",sigma[3])

    PARS <- NAMES[ grepl("error",NAMES) ]
    for(P in PARS) { par[P] <- CTMM$error[ substr(P,nchar("error ?"),nchar(P)) ] }
  }
  else
  {
    getter("major",sigma[1]) # sigma[x,x]
    getter("minor",sigma[4]) # sigma[y,y]
    getter("angle",sigma[2]) # sigma[x,y]

    PARS <- NAMES[ grepl("error",NAMES) ]
    for(P in PARS) { par[P] <- CTMM$error[ substr(P,nchar("error ?"),nchar(P)) ]^2 }
  }

  DRIFT <- drift.pars(CTMM)
  if(length(DRIFT))
  {
    DPARS <- names(DRIFT)
    DPARS <- DPARS[ DPARS %in% NAMES ]
    par[DPARS] <- DRIFT[DPARS]
  }

  timelink <- CTMM$timelink.par
  TIMELINK <- which(grepl("timelink",NAMES))
  if(length(timelink) && length(TIMELINK)) { par[TIMELINK] <- timelink }

  tau <- CTMM$tau
  getter("tau position",if(length(tau)>0) { tau[1] } else { 0 })
  getter("tau velocity",if(length(tau)>1) { tau[2] } else { 0 })
  getter("tau",if(length(tau)>0) { min(tau) } else { 0 }) # identical timescales
  getter("omega")

  getter("circle")

  # only relative variances matter - return variance fraction [0,1]
  if(profile)
  {
    UERE.RMS <- paste("error",names(CTMM$error)) %in% PARS
    UERE.RMS <- CTMM$error[UERE.RMS]

    MAX <- profiled.var(CTMM,UERE.RMS=UERE.RMS,AVE=FALSE)

    # this needs to be consistent with profiled.var() and set.parameters()
    if(length(PARS))
    { par[PARS] <- par[PARS] / sqrt(MAX) }

    SPARS <- c("major","minor")
    SPARS <- SPARS[SPARS %in% NAMES]

    # this needs to be consistent with profiled.var() and set.parameters()
    if(length(SPARS)) { par[SPARS] <- par[SPARS] / MAX }
  } # end profile

  return(par)
}


# total variance to profile
profiled.var <- function(CTMM,sigma=CTMM$sigma,UERE.RMS=CTMM$error,DT=1,AVE=FALSE)
{
  PRO.VAR <- sum(eigenvalues.covm(sigma))*DT # really profiling the variance with mean?

  if(any(UERE.RMS>0)) { PRO.VAR <- PRO.VAR + sum(UERE.RMS^2) } # comparable error variance (@DOP==1)

  if(AVE) # variance per dimension
  {
    AXES <- length(CTMM$axes)
    PRO.VAR <- PRO.VAR / AXES
  }

  return(PRO.VAR)
}


# clean up parameter arrays
clean.parameters <- function(par,profile=FALSE,linear.cov=FALSE,timelink="identity")
{
  # multiple states
  if(class(par)[1]=='list')
  {
    par <- lapply(par,function(p){clean.parameters(p,profile=profile,linear.cov=linear.cov,timelink=timelink)})
    return(par)
  }

  NAMES <- names(par)

  for(P in POSITIVE.PARAMETERS) { if(P %in% NAMES) { par[P] <- clamp(par[P],0,Inf) } }

  # arbitrary error parameters > 0
  EP <- NAMES[ grepl("error",NAMES) ]
  par[EP] <- clamp(par[EP],0,Inf)

  # make sure relative variance is constrained properly
  if(profile)
  {
    SP <- c('major','minor') # major can't really be here
    SP <- SP[SP %in% NAMES]

    # 0 <= VBOUND <= 1
    VBOUND <- sum(par[SP]) + sum(par[EP]^2)
    if(VBOUND>1)
    {
      if(length(EP)) { par[EP] <- par[EP] / sqrt(VBOUND) }

      if(length(SP))
      {
        # par[SP] <- par[SP] / VBOUND # new major
        par['minor'] <- 0 # old major
        par["angle"] <- par["angle"] - pi/2
      }
    } # else 'major' > 0
    else if(length(SP) && VBOUND<par['minor'])
    {
      par['minor'] <- VBOUND # old major is true minor
      par["angle"] <- par["angle"] - pi/2 # rotate old major to true major
    }
  }
  else if(!linear.cov && all(c("major","minor") %in% NAMES) && par["minor"]>par['major'])
  {
    # swap major and minor axis --- if not profiling major axis
    VARS <- par[c('minor','major')]
    names(VARS) <- NULL
    par[c('major','minor')] <- VARS

    par["angle"] <- par["angle"] - pi/2
  }
  else if(linear.cov && all(c("major","minor") %in% NAMES))
  {
    # sigma[x,x] sigma[y,y] --- positive semi-definite
    MAX <- sqrt(prod(par[c('major','minor')]))
    # sigma[x,y]
    par['angle'] <- clamp(par['angle'],-MAX,MAX)
  }

  # shift back to (-pi/2,pi/2)
  if("angle" %in% NAMES && !linear.cov)
  { par["angle"] <- (((par["angle"]/pi+1/2) %% 1) - 1/2)*pi }

  # swap tau velocity and tau position if they become reversed (to keep parscale sane)
  P <- c("tau position","tau velocity")
  if(all(P %in% NAMES)) { par[P] <- sort(par[P],decreasing=TRUE) }

  # clean timelink parameters
  if(timelink!="identity")
  {
    PAR <- grepl("timelink",NAMES)
    if(any(PAR)) { par[PAR] <- timelink.clean(par[PAR],timelink=timelink) }
  }

  return(par)
}


# write autocovariance parameter array into a CTMM object
set.parameters <- function(CTMM,par,profile=FALSE,linear.cov=FALSE)
{
  NAMES <- names(par)
  AXES <- length(CTMM$axes)

  beta <- CTMM$beta
  BETA <- names(beta)
  beta[] <- NA
  for(B in BETA) { beta[B] <- par[B] }
  CTMM$beta <- beta

  sigma <- CTMM$sigma
  if(!linear.cov)
  {
    sigma <- attr(sigma,"par")
    NAME <- "major"; if(NAME %in% NAMES) { sigma[NAME] <- par[NAME]; if(AXES>1) { sigma['minor'] <- par[NAME] } } # in case of isotropic
    NAME <- "minor"; if(NAME %in% NAMES) { sigma[NAME] <- par[NAME] }
    NAME <- "angle"; if(NAME %in% NAMES) { sigma[NAME] <- par[NAME] }

    EP <- NAMES[ grepl("error",NAMES) ]
    for(P in EP) { CTMM$error[ substr(P,nchar("error ?"),nchar(P)) ] <- par[P] }

    if(profile)
    {
      SP <- c("major","minor") # major shouldn't be here
      SP <- SP[SP %in% NAMES]

      MAX <- sum(par[SP]) + sum(par[EP]^2)
      sigma['major'] <- 1 - MAX # remaining fraction of total variance
    }
  }
  else # major,minor,angle is storing xx,yy,xy in linear representation # no profiling ever here
  {
    sigma <- methods::getDataPart(sigma)
    NAME <- "major"; if(NAME %in% NAMES) { if(AXES>1) { diag(sigma) <- par[c(NAME,NAME)] } else { sigma <- par[NAME] } }
    NAME <- "minor"; if(NAME %in% NAMES) { sigma[4] <- par[NAME] }
    NAME <- "angle"; if(NAME %in% NAMES) { sigma[2:3] <- par[c(NAME,NAME)] }

    PARS <- NAMES[ grepl("error",NAMES) ]
    for(P in PARS) { CTMM$error[ substr(P,nchar("error ?"),nchar(P)) ] <- sqrt(par[P]) }
  }
  CTMM$sigma <- covm(sigma,axes=CTMM$axes)

  DRIFT <- drift.pars(CTMM)
  if(length(DRIFT))
  {
    DPARS <- names(DRIFT)
    DPARS <- DPARS[ DPARS %in% NAMES ]
    DRIFT[DPARS] <- par[DPARS]
    CTMM <- drift.assign(CTMM,DRIFT)
  }

  timelink <- par[grepl("timelink",NAMES)]
  if(length(timelink)) { CTMM$timelink.par <- timelink }

  NAME <- "tau position"; if(NAME %in% NAMES) { CTMM$tau[1] <- par[NAME] ; names(CTMM$tau)[1] <- 'position' }
  NAME <- "tau velocity"; if(NAME %in% NAMES) { CTMM$tau[2] <- par[NAME] ; names(CTMM$tau)[2] <- 'velocity' }
  NAME <- "tau"; if(NAME %in% NAMES) { CTMM$tau <- c(1,1)*par[NAME] ; names(CTMM$tau) <- c('position','velocity') } # identical timescales
  NAME <- "omega"; if(NAME %in% NAMES) { CTMM$omega <- par[NAME] }

  NAME <- "circle"; if(NAME %in% NAMES) { CTMM[[NAME]] <- par[NAME] }

  return(CTMM)
}


# copy autocovariance parameters 'par' from 'value' to 'x', which might have more parameters
copy.parameters <- function(x,value,par=value$features,destructive=TRUE)
{
  Pv <- par
  Px <- x$features

  if('minor' %in% Pv) { x$isotropic <- value$isotropic }
  if('minor' %in% Pv || 'minor' %nin% Px) # copy over full covariance matrix
  { x$sigma <- value$sigma }
  else # copy over variance only
  {
    sigma <- scale.covm(x$sigma) # var-1
    sigma <- scale.covm( sigma, var.covm(value$sigma,ave=TRUE) )
    sigma@isotropic=FALSE
    x$sigma <- sigma
  }

  timelink <- which(grepl('timelink',Pv))
  if(length(timelink))
  {
    x$timelink <- value$timelink
    x$timelink.par[1:length(value$timelink.par)] <- value$timelink.par
  }

  if("tau position" %in% Pv) { x$tau[1] <- value$tau[1] }
  if("tau velocity" %in% Pv) { x$tau[2] <- value$tau[2] }
  if("tau" %in% Pv) { x$tau <- value$tau }
  if("omega" %in% Pv) { x$omega <- value$omega }
  if("circle" %in% Pv) { x$circle <- value$circle }

  PARS <- Pv[ grepl("error",Pv) ]
  for(P in PARS)
  {
    P <- substr(P,nchar("error ?"),nchar(P))
    x$error[P] <- value$error[P]
  }

  DRIFT <- drift.pars(value)
  if(length(DRIFT))
  {
    DRIFTO <- drift.pars(x)
    DPARS <- names(DRIFT)
    DPARS <- DPARS[ DPARS %in% names(DRIFTO) ]
    DRIFTO[DPARS] <- DRIFT[DPARS]
    x <- drift.assign(x,DRIFTO)
  }

  # don't think I need this
  if("MLE" %in% names(x))
  {
    if("MLE" %in% names(value)) { value <- value$MLE }
    x$MLE <- copy.parameters(x$MLE,value,par=par,destructive=destructive)
  }

  return(x)
}


get.states <- function(CTMM)
{
  dynamics <- CTMM$dynamics
  if(is.null(dynamics) || dynamics=="stationary")
  { states <- NULL }
  else if(dynamics=="change.point")
  { states <- levels(CTMM[[dynamics]]$state) }

  return(states)
}
