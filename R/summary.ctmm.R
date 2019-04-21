###################################################
# Calculate good CIs for other functions
confint.ctmm <- function(model,alpha=0.05,UNICODE=FALSE)
{
  # make sure parameters are named correctly
  model <- ctmm.ctmm(model)

  model <- get.taus(model,zeroes=TRUE)
  tau <- model$tau
  tau <- tau[tau<Inf]
  K <- length(tau)

  COV <- model$COV

  par <- NULL
  NAME <- NULL

  # timescale uncertainty: can hit 0 and Inf
  if(K>0)
  {
    NAME <- model$tau.names # canonical parameter names
    if(K==2 && model$omega) # oscillatory model
    {
      TAU <- model$TAU # human-readible parameters
      COV.TAU <- model$J.TAU.tau %*% COV[NAME,NAME] %*% t(model$J.TAU.tau)
      for(k in 1:K) { par <- rbind(par,ci.tau(TAU[k],COV.TAU[k,k],alpha=alpha)) }

      if(UNICODE) { NAME <- c("\u03C4[decay]","\u03C4[period]") }
      else { NAME <- c("tau decay","tau period") }
    }
    else if(K>1 && tau[1]==tau[2]) # identical timescales
    {
      par <- rbind(par,ci.tau(tau[1],COV[NAME,NAME],alpha=alpha))
      if(UNICODE) { NAME <- "\u03C4" }
    }
    else # OU, OUF, IOU
    {
      for(k in 1:K) { par <- rbind(par,ci.tau(tau[k],COV[NAME[k],NAME[k]],alpha=alpha)) }
      if(UNICODE) { NAME <- paste0("\u03C4[",names(tau),"]") }
    }
  }
  # else IID, BM

  # circulation period
  circle <- model$circle
  if(circle)
  {
    NAME <- c(NAME,"circle")
    VAR <- COV["circle","circle"]
    CI <- stats::qnorm(c(alpha/2,0.5,1-alpha/2),mean=circle,sd=sqrt(VAR))

    # does frequency CI cross zero?
    if(CI[1]*CI[3]<0)
    {
      # quick fix
      CI <- abs(CI)
      CI[3] <- max(CI)
      CI[1] <- 0
    }

    # frequency to period
    CI <- sort(2*pi/abs(CI))

    par <- rbind(par,CI)
  }

  # standard area uncertainty: chi-square
  NAME <- c("area",NAME)
  AREA <- area.covm(model$sigma)
  DOF <- 2*DOF.area(model)
  par <- rbind(chisq.ci(AREA,DOF=DOF,alpha=alpha),par)

  rownames(par) <- NAME

  return(par)
}


# timescale uncertainty: can hit 0 and Inf
ci.tau <- function(tau,COV,alpha=0.05,min=0,max=Inf)
{
  z <- stats::qnorm(1-alpha/2)
  if(is.nan(COV)) { COV <- Inf }

  # tau normal for lower CI
  CI <- tau + c(1,-1)*z*sqrt(COV)

  # lower CI of f==1/tau normal for upper CI of tau
  CI <- c(CI, (1/tau + c(1,-1)*z*sqrt(COV/tau^4))^-1)

  # take most conservative estimates
  CI <- sort(range(CI,na.rm=TRUE))

  # enforce boundary constraints
  CI <- c(max(CI[1],min),min(CI[2],max))

  CI <- c(CI[1],tau,CI[2])
  return(CI)
}


###
summary.ctmm <- function(object,level=0.95,level.UD=0.95,units=TRUE,IC="AICc",MSPE="position",...)
{
  CLASS <- class(object)
  if(CLASS=="ctmm")
  { return(summary.ctmm.single(object,level=level,level.UD=level.UD,units=units)) }
  else if(CLASS=="list")
  { return(summary.ctmm.list(object,level=level,level.UD=level.UD,IC=IC,MSPE=MSPE,units=units)) }
}


# stochastic square speed contribution
rand.speed <- function(CTMM)
{
  MSD <- var.covm(CTMM$sigma,ave=FALSE)
  Omega <- CTMM$circle
  Omega2 <- CTMM$Omega2 + Omega^2
  MLE <- MSD * Omega2

  COV.RAN <- axes2var(CTMM,MEAN=FALSE)

  GRAD <- Omega2
  PARS <- "variance"

  PARS <- c(PARS,CTMM$tau.names)
  GRAD <- c(GRAD,MSD*CTMM$J.Omega2)

  if(CTMM$circle)
  {
    PARS <- c(PARS,"circle")
    GRAD <- c(GRAD,2*MSD*Omega)
  }

  COV.RAN <- COV.RAN[PARS,PARS]
  COV.RAN <- c(GRAD %*% COV.RAN %*% GRAD)

  return(list(MLE=MLE,COV=COV.RAN))
}


######################################################
# Summarize results
summary.ctmm.single <- function(object, level=0.95, level.UD=0.95, units=TRUE, ...)
{
  alpha <- 1-level
  alpha.UD <- 1-level.UD

  object <- get.taus(object) # populate with parameter stuff
  drift <- get(object$mean)

  circle <- object$circle
  tau <- object$tau
  range <- object$range
  tau <- tau[tau<Inf]

  if("COV" %in% names(object))
  { COV <- object$COV }
  else # fill in with infinite covariance
  {
    P <- id.parameters(object,profile=FALSE)$NAMES
    COV <- diag(Inf,nrow=length(P))
    dimnames(COV) <- list(P,P)
    object$COV <- COV
  }
  P <- nrow(COV)

  if(is.null(object$COV.mu))
  {
    object$COV.mu <- diag(Inf,nrow=length(object$mu))
    object$DOF.mu <- diag(0,nrow=length(object$mu))
  }

  par <- confint.ctmm(object,alpha=alpha,UNICODE=TRUE)
  # where to store unit information
  K <- nrow(par)
  name <- character(K)
  scale <- numeric(K)

  # standard area to home-range area
  par[1,] <- -2*log(alpha.UD)*pi*par[1,]

  # pretty area units   # do we convert units
  unit.list <- unit.par(par[1,],"area",SI=!units)
  name[1] <- unit.list$name
  scale[1] <- unit.list$scale

  # pretty time units
  P <- nrow(par)
  if(P>1)
  {
    for(i in 2:P)
    {
      # set scale by upper CI if point estimate is zero
      unit.list <- unit.par(par[i,],"time",SI=!units)
      name[i] <- unit.list$name
      scale[i] <- unit.list$scale
    }
  }

  # can we estimate speed?
  if(continuity(object)>1)
  {
    # stochastic contribution
    STUFF <- rand.speed(object)
    ms <- STUFF$MLE
    var.ms <- STUFF$COV

    # include mean
    MSPEED <- drift@speed(object)
    ms <- ms + MSPEED$EST
    var.ms <- var.ms + MSPEED$VAR

    # chi^2 degrees of freedom
    DOF.speed <- 2*ms^2/var.ms

    # root mean square velocity
    # pretty units
    rms <- sqrt(chisq.ci(ms,COV=var.ms,alpha=alpha))
    unit.list <- unit.par(rms,"speed",SI=!units)
    name <- c(name,unit.list$name)
    scale <- c(scale,unit.list$scale)

    par <- rbind(par,rms)
    rownames(par)[nrow(par)] <- "speed"
  }
  else
  { DOF.speed <- 0 }

  # did we estimate errors?
  if("error" %in% rownames(object$COV))
  {
    error <- object$error
    VAR <- object$COV["error","error"]
    # convert to chi^2
    VAR <- (2*error)^2 * VAR
    error <- error^2
    # CIs
    error <- chisq.ci(error,COV=VAR,alpha=alpha)
    # back to meters/distance
    error <- sqrt(error)
    unit.list <- unit.par(error,"length",SI=!units)
    name <- c(name,unit.list$name)
    scale <- c(scale,unit.list$scale)

    par <- rbind(par,error)
    rownames(par)[nrow(par)] <- "error"
  }

  # Fix unit choice
  par <- par/scale

  # affix units
  rownames(par) <- paste(rownames(par)," (",name,")",sep="")

  if(!range) { par <- par[-1,] } # delete off "area" (really diffusion)

  # anything else interesting from the mean function
  par <- rbind(drift@summary(object,level,level.UD),par)

  colnames(par) <- c("low","ML","high")

  SUM <- list(name=name.ctmm(object))

  # affix DOF info
  # only valid for processes with a stationary mean
  SUM$DOF <- c( DOF.mean(object) , DOF.area(object) , DOF.speed/2 )
  names(SUM$DOF) <- c("mean","area","speed")

  SUM$CI <- par

  return(SUM)
}
#methods::setMethod("summary",signature(object="ctmm"), function(object,...) summary.ctmm(object,...))

#DOF of area
DOF.area <- function(CTMM)
{
  sigma <- CTMM$sigma@par
  AREA <- area.covm(CTMM$sigma)

  if(CTMM$isotropic)
  { VAR <- CTMM$COV['major','major'] }
  else
  {
    P <- c('major','minor')
    GRAD <- c( sqrt(sigma['minor']) , sigma['major']/sqrt(sigma['minor'])/2 )
    VAR <- c( GRAD %*% CTMM$COV[P,P] %*% GRAD )
  }

  if(CTMM$range) { DOF <- AREA^2/abs(VAR) } else { DOF <- 0 }

  if(is.nan(DOF)) { DOF <- 0 }

  return(DOF)
}

#########
DOF.mean <- function(CTMM)
{
  if(!CTMM$range || "COV.mu" %nin% names(CTMM)) { return(0) }

  sigma <- CTMM$sigma
  COV <- CTMM$COV.mu

  if(length(dim(COV))==4) { COV <- COV[,1,1,] } # take only stationary mean COV

  sigma <- sqrtm.covm(sigma)
  # symmetric under trace and det
  DOF <- sigma %*% PDsolve(COV) %*% sigma
  DOF <- mean(diag(DOF))

  return(DOF)
}

#######
DOF.speed <- function(CTMM) { return( summary(CTMM)$DOF['speed'] ) }
