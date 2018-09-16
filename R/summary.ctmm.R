###################################################
# Calculate good CIs for other functions
confint.ctmm <- function(model,alpha=0.05)
{
  tau <- model$tau
  tau <- tau[tau<Inf]
  K <- length(tau)

  COV <- model$COV

  par <- NULL
  NAME <- NULL

  # timescale uncertainty: can hit 0 and Inf
  if(K>0)
  {
    NAME <- paste("tau",names(tau))
    for(k in 1:K) { par <- rbind(par,ci.tau(tau[k],COV[NAME[k],NAME[k]],alpha=alpha)) }
  }

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
  GM.sigma <- model$sigma@par["area"]
  COV <- COV["area","area"]
  par <- rbind(chisq.ci(GM.sigma,COV=COV,alpha=alpha),par)

  rownames(par) <- NAME

  return(par)
}

# timescale uncertainty: can hit 0 and Inf
ci.tau <- function(tau,COV,alpha=0.05,min=0,max=Inf)
{
  z <- stats::qnorm(1-alpha/2)

  # tau normal for lower CI
  CI <- tau + c(1,-1)*z*sqrt(COV)

  # lower CI of f==1/tau normal for upper CI of tau
  CI <- c(CI, (1/tau + c(1,-1)*z*sqrt(COV/tau^4))^-1)

  # take most conservative estimates
  CI <- range(CI,na.rm=TRUE)

  # enforce boundary constraints
  CI <- c(max(CI[1],min),min(CI[2],max))

  CI <- c(CI[1],tau,CI[2])
  return(CI)
}

###
summary.ctmm <- function(object,level=0.95,level.UD=0.95,units=TRUE,IC="AICc",MSPE=TRUE,...)
{
  CLASS <- class(object)
  if(CLASS=="ctmm")
  { return(summary.ctmm.single(object,level=level,level.UD=level.UD,units=units)) }
  else if(CLASS=="list")
  { return(summary.ctmm.list(object,level=level,level.UD=level.UD,IC=IC,MSPE=MSPE,units=units)) }
}

######################################################
# Summarize results
summary.ctmm.single <- function(object, level=0.95, level.UD=0.95, units=TRUE, ...)
{
  alpha <- 1-level
  alpha.UD <- 1-level.UD

  drift <- get(object$mean)

  circle <- object$circle
  tau <- object$tau
  range <- object$range
  tau <- tau[tau<Inf]
  K <- length(tau)

  AM.sigma <- mean(diag(object$sigma))
  GM.sigma <- object$sigma@par["area"]
  ecc <- object$sigma@par["eccentricity"]

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

  # where to store unit information
  name <- character(K+1)
  scale <- numeric(K+1)

  par <- confint.ctmm(object,alpha=alpha)

  # standard area to home-range area
  par[1,] <- -2*log(alpha.UD)*pi*par[1,]

  # pretty area units   # do we convert units
  unit.list <- unit.par(par,"area",SI=!units)
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
    COV <- area2var(object,MEAN=TRUE)
    COV <- rm.name(COV,"error")

    # RMS velocity
    Omega2 <- 1/prod(tau)
    grad <- -Omega2/tau

    # contribution from circulation
    omega2 <- 0
    if(circle)
    {
      omega2 <- circle^2
      grad <- c(grad,2*circle)
    }

    # contribution from sigma
    # GM.sigma <- cosh(ecc)*AM.sigma
    ms <- (2*AM.sigma)*(Omega2+omega2)
    grad <- 2*c(Omega2+omega2, AM.sigma*grad)
    var.ms <- c((grad) %*% COV %*% (grad))

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

  SUM <- list()

  # affix DOF info
  # only valid for processes with a stationary mean
  if(object$range)
  {
    SUM$DOF <- c( DOF.mean(object) , DOF.area(object) , DOF.speed/2 )
    names(SUM$DOF) <- c("mean","area","speed")
  }

  SUM$CI <- par

  return(SUM)
}
#methods::setMethod("summary",signature(object="ctmm"), function(object,...) summary.ctmm(object,...))

#DOF of area
DOF.area <- function(CTMM) { CTMM$sigma@par["area"]^2/abs(CTMM$COV["area","area"]) }

#########
DOF.mean <- function(CTMM)
{
  DOF <- CTMM$DOF.mu
  DIM <- dim(DOF)

  if(!is.null(DIM))
  {
    # stationary contribution
    if(length(DIM)==4) { DOF <- DOF[,1,1,] }
    # average of x-y DOFs
    DOF <- mean(diag(DOF))
  }

  return(DOF)
}
