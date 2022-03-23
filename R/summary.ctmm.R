###################################################
# Calculate good CIs for other functions
confint.ctmm <- function(model,alpha=0.05,UNICODE=FALSE)
{
  # make sure parameters are named correctly
  model <- ctmm.ctmm(model)
  isotropic <- model$isotropic
  sigma <- model$sigma

  model <- get.taus(model,zeroes=TRUE)
  tau <- model$tau
  tau <- tau[tau<Inf]
  K <- length(tau)

  COV <- model[["COV"]]
  POV <- model[["POV"]]
  COV.POV <- model[["COV.POV"]]
  VARS <- dimnames(COV)[[1]]

  par <- NULL
  NAMES <- NULL

  # standard area uncertainty: chi-square
  NAMES <- c(NAMES,"area")
  AREA <- area.covm(model$sigma)
  DOF <- 2*DOF.area(model)
  VAR <- 2*AREA^2/DOF
  AREA <- chisq.ci(AREA,DOF=DOF,alpha=alpha)
  par <- rbind(par,AREA)

  if(!is.null(POV))
  {
    J <- J.zero(POV)

    if(isotropic[1])
    {
      P <- 'major'
      AREA <- sigma@par[P]
      # CoV^2
      PAR <- POV[P,P] / AREA^2
      # dArea/dpar
      J[P,P] <- 1
      J <- quad2lin(J,diag=TRUE)
    }
    else
    {
      P <- c('major','minor')
      AREA <- sqrt(prod(sigma@par[P]))
      # dArea/dpar
      J0 <- 1/2 * AREA / sigma@par[P]
      # CoV^2
      PAR <- (J0 %*% POV[P,P] %*% t(J0))/AREA^2
      # Jacobian for VAR[Area]
      J[P,P] <- J0
      J <- quad2lin(J,diag=TRUE)
    }
    # VAR[PAR|POV] + VAR[PAR|Area]
    VAR <- tr(J %*% COV.POV %*% t(J))/AREA^4 + (-2*PAR/AREA)^2*VAR
    PAR <- chisq.ci(PAR,VAR=VAR,alpha=alpha)
    PAR <- sqrt(PAR)

    NAMES <- c(NAMES,"CoV[area]")
    par <- rbind(par,PAR)
  }

  # timescale uncertainty: can hit 0 and Inf
  if(K>0)
  {
    NAME <- LABEL <- model$tau.names # canonical parameter names
    TAU <- model$tau # human-readable parameters
    J0 <- diag(1,length(NAME))

    if(K==2 && "omega" %in% VARS) # OUO
    {
      if(UNICODE) { LABEL <- c("\u03C4[decay]","\u03C4[period]") }
      else { LABEL <- c("tau decay","tau period") }

      J0 <- model$J.TAU.tau
    } # end OUO
    else if(K>1 && "tau" %in% VARS) # OUf
    {
      if(UNICODE) { LABEL <- "\u03C4" }
    } # end OUf
    else # OU, OUF, IOU
    {
      if(UNICODE) { LABEL <- paste0("\u03C4[",names(tau),"]") }
    } # end OU, OUF, IOU

    VAR.TAU <- J0 %*% COV[NAME,NAME] %*% t(J0)
    VAR.TAU <- diag(VAR.TAU)

    if(!is.null(POV))
    {
      PAR <- diag( J0 %*% POV[NAME,NAME] %*% t(J0) )/TAU^2
      J <- J.zero(POV)
      J[NAME,NAME] <- J0
      J <- quad2lin(J,diag=TRUE)
      VAR <- diag( J %*% COV.POV %*% t(J) )
      names(VAR) <- model$features
      VAR <- VAR[NAME]/TAU^4 + (-2*PAR/TAU)^2*VAR.TAU
    }

    for(k in 1:length(NAME))
    {
      par <- rbind(par,ci.tau(TAU[k],VAR.TAU[k],alpha=alpha))
      NAMES <- c(NAMES,LABEL[k])

      if(!is.null(POV))
      {
        par <- rbind(par,sqrt(chisq.ci(PAR[k],VAR[k],alpha=alpha)))
        NAMES <- c(NAMES,paste0("CoV[",LABEL[k],"]"))
      }
    } # end for k
  }
  # else IID, BM

  # circulation period
  circle <- model$circle
  if("circle" %in% VARS)
  {
    VAR <- COV["circle","circle"]
    if(VAR<Inf)
    { CI <- stats::qnorm(c(alpha/2,0.5,1-alpha/2),mean=circle,sd=sqrt(VAR)) }
    else
    { CI <- c(-Inf,circle,Inf) }

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
    NAMES <- c(NAMES,"circle")

    if(!is.null(POV))
    {
      f <- 2*pi/circle
      J0 <- -2*pi/circle^2
      PAR <- (J0 * POV['circle','circle'] * J0)/f^2

      J <- J.zero(POV)
      J["circle","circle"] <- J0
      J <- quad2lin(J,diag=TRUE)
      VAR <- tr(J %*% COV.POV %*% t(J))/f^4 + (-2*PAR/circle)^2*COV["circle","circle"]

      CI <- sqrt( chisq.ci(PAR,VAR=VAR,alpha=alpha) )
      par <- rbind(par,CI)
      NAMES <- c(NAMES,"CoV[circle]")
    }
  }

  rownames(par) <- NAMES
  return(par)
}

#######
J.zero <- function(POV)
{
  J <- array(0,dim(POV))
  dimnames(J) <- dimnames(POV)
  return(J)
}

# timescale uncertainty: can hit 0 and Inf
ci.tau <- function(tau,COV,alpha=0.05,min=0,max=Inf)
{
  z <- stats::qnorm(1-alpha/2)
  if(is.nan(COV)) { COV <- Inf }

  # tau normal for lower CI
  # CI <- tau + c(1,-1)*z*sqrt(COV)
  # lower CI of f==1/tau normal for upper CI of tau
  # CI <- c(CI, (1/tau + c(1,-1)*z*sqrt(COV/tau^4))^-1)

  # tau chi^2 for lower CI
  CI <- chisq.ci(tau,VAR=COV,alpha=alpha)
  # tau chi^2 for upper CI
  CI <- c(CI,1/chisq.ci(1/tau,VAR=COV/tau^4,alpha=alpha))

  # take most conservative estimates
  CI <- sort(range(CI,na.rm=TRUE))

  # enforce boundary constraints
  CI <- c(max(CI[1],min),min(CI[2],max))

  CI <- c(CI[1],tau,CI[2])
  return(CI)
}


###
summary.ctmm <- function(object,level=0.95,level.UD=0.95,units=TRUE,IC=NULL,MSPE=NULL,...)
{
  CLASS <- class(object)[1]
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
  GRAD <- rep(0,nrow(COV.RAN))
  names(GRAD) <- rownames(COV.RAN)
  # somewhat replicated wrt old GRAD code, but need this for mean.ctmm summary
  J <- rep(0,length(CTMM$features))
  names(J) <- CTMM$features

  GRAD["variance"] <- Omega2
  if(CTMM$isotropic[1])
  { J["major"] <- MLE/CTMM$sigma@par['major'] }
  else
  { J[c("major","minor")] <- MLE/MSD }

  GRAD[CTMM$tau.names] <- J[CTMM$tau.names] <- MSD*CTMM$J.Omega2

  if(CTMM$circle)
  { GRAD["circle"] <- J["circle"] <- 2*MSD*Omega }

  COV.RAN <- c(GRAD %*% COV.RAN %*% GRAD)

  return(list(MLE=MLE,COV=COV.RAN,J=J))
}


######################################################
# Summarize results
summary.ctmm.single <- function(object, level=0.95, level.UD=0.95, units=TRUE, ...)
{
  alpha <- 1-level
  alpha.UD <- 1-level.UD
  axes <- object$axes

  object <- get.taus(object) # populate with parameter stuff
  drift <- get(object$mean)

  circle <- object$circle
  tau <- object$tau
  range <- object$range
  tau <- tau[tau<Inf]

  if("COV" %in% names(object))
  {
    COV <- object$COV
    POV <- object$POV
    COV.POV <- object$COV.POV
  }
  else # fill in with infinite covariance
  {
    UERE.FIT <- object$error>0
    P <- id.parameters(object,profile=FALSE,UERE.FIT=UERE.FIT)$NAMES
    COV <- diag(Inf,nrow=length(P))
    dimnames(COV) <- list(P,P)
    object$COV <- COV
    POV <- NULL
  }
  P <- nrow(COV)

  if(is.null(object$COV.mu))
  {
    object$COV.mu <- diag(Inf,nrow=length(axes))
    object$DOF.mu <- diag(0,nrow=length(axes))
    dimnames(object$COV.mu) <- dimnames(object$DOF.mu) <- list(axes,axes)
  }

  ### AREA
  par <- confint.ctmm(object,alpha=alpha,UNICODE=TRUE)
  # where to store unit information
  K <- nrow(par)
  name <- rep("",K)
  scale <- rep(1,K)
  names(name) <- names(scale) <- rownames(par)

  # standard area to home-range area
  par["area",] <- -2*log(alpha.UD)*pi*par[1,]

  # pretty area units   # do we convert units
  unit.list <- unit.par(par["area",],"area",SI=!units)
  name["area"] <- unit.list$name
  scale["area"] <- unit.list$scale

  # pretty time units
  I <- !grepl("CoV",rownames(par),fixed=TRUE)
  I <- which(I)[-1] # all but area
  for(i in I)
  {
    # set scale by upper CI if point estimate is zero
    unit.list <- unit.par(par[i,],"time",SI=!units)
    name[i] <- unit.list$name
    scale[i] <- unit.list$scale
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
    rms <- sqrt(chisq.ci(ms,VAR=var.ms,alpha=alpha))
    rms <- rms / chi.bias(DOF.speed) # if chi^2 estimate was unbiased, then chi estimate is now unbiased too
    unit.list <- unit.par(rms,"speed",SI=!units)
    name <- c(name,unit.list$name)
    scale <- c(scale,unit.list$scale)

    par <- rbind(par,rms)
    rownames(par)[nrow(par)] <- "speed"

    if(!is.null(POV))
    {
      J0 <- STUFF$J
      PAR <- c(J0 %*% POV %*% J0)/ms^2

      J0 <- rbind(J0,array(0,length(J0)-1:0))
      J <- quad2lin(J0,diag=TRUE)
      VAR <- tr(J %*% COV.POV %*% t(J))/ms^4 + (-2*PAR/ms)^2*var.ms

      # CoV^2 of MS[speed]
      PAR <- chisq.ci(PAR,VAR=VAR,alpha=alpha)
      # CoV^2 of RMS[speed]
      PAR <- PAR/4
      # CoV of RMS[speed]
      PAR <- sqrt(PAR)

      par <- rbind(par,PAR)
      rownames(par)[nrow(par)] <- "CoV[speed]"
      name <- c(name,"")
      scale <- c(scale,1)
    }
  }
  else
  { DOF.speed <- 0 }

  # can we estimate diffusion?
  if(length(tau) && tau[1]>0)
  {
    STUFF <- diffusion(object,finish=FALSE)
    DOF.diffusion <- STUFF$DOF
    D <- chisq.ci(STUFF$D,VAR=STUFF$VAR,level=level)
    unit.list <- unit.par(D,"diffusion",SI=!units)
    name <- c(name,unit.list$name)
    scale <- c(scale,unit.list$scale)

    par <- rbind(par,D)
    rownames(par)[nrow(par)] <- "diffusion"

    if(!is.null(POV))
    {
      D <- D[2]
      J0 <- STUFF$J
      PAR <- c(J0 %*% POV %*% J0)/D^2

      J0 <- rbind(J0,array(0,length(J0)-1:0))
      J <- quad2lin(J0,diag=TRUE)
      VAR <- tr(J %*% COV.POV %*% t(J))/D^4 + (-2*PAR/D)^2*STUFF$VAR

      # CoV^2
      PAR <- chisq.ci(PAR,VAR=VAR,alpha=alpha)
      PAR <- sqrt(PAR)
      par <- rbind(par,PAR)
      rownames(par)[nrow(par)] <- "CoV[diffusion]"
      name <- c(name,"")
      scale <- c(scale,1)
    }
  }
  else
  { DOF.diffusion <- 0 }

  # did we estimate errors?
  PARS <- rownames(object$COV)
  PARS <- PARS[grepl("error",PARS)]
  for(P in PARS)
  {
    error <- object$error[substr(P,nchar("error #"),nchar(P))]
    VAR <- object$COV[P,P]
    # convert to chi^2
    VAR <- (2*error)^2 * VAR
    error <- error^2
    # CIs
    error <- chisq.ci(error,VAR=VAR,alpha=alpha)
    # back to meters/distance
    error <- sqrt(error)
    unit.list <- unit.par(error,"length",SI=!units)
    name <- c(name,unit.list$name)
    scale <- c(scale,unit.list$scale)

    par <- rbind(par,error)
    rownames(par)[nrow(par)] <- P
  }

  # did we estimate RSF betas?
  PARS <- names(object$beta)
  for(i in 1%:%length(PARS)) # stack in front
  {
    P <- PARS[i]
    RSF <- object$beta[P]
    VAR <- diag(object$COV)[P]
    VAR <- nant(VAR,Inf)
    RSF <- norm.ci(RSF,VAR,level=level)

    NAME <- paste0("1/",P)
    SCALE <- 1

    if(!is.null(POV))
    {
      Q <- paste0(P,"-",P)
      PAR <- chisq.ci(POV[P,P],VAR=COV.POV[Q,Q],alpha=alpha)
      PAR <- sqrt(PAR)
      RSF <- rbind(RSF,PAR)

      NAME <- c(NAME,NAME)
      SCALE <- c(1,1)
      P <- c(P,paste0("SD[",P,"]"))
    }

    # PU <- substr(P,nchar("RSF *"),nchar(P))
    # name <- c(paste0("1/",PU),name)
    name <- c(NAME,name)
    scale <- c(SCALE,scale)

    par <- rbind(RSF,par)
    rownames(par)[1:length(P)] <- P
  }

  # Fix unit choice
  par <- par/scale

  # affix units if given units
  DO <- name!=""
  rownames(par)[DO] <- paste(rownames(par)[DO]," (",name[DO],")",sep="")

  if(!range) { par <- par[-1,,drop=FALSE] } # delete off "area" (really diffusion)... stop R drop :(

  # anything interesting from the timelink function
  par <- rbind(timelink.summary(object,level=level),par)

  # anything else interesting from the mean function
  par <- rbind(drift@summary(object,level,level.UD),par)

  colnames(par) <- NAMES.CI

  SUM <- list(name=name.ctmm(object))

  # affix DOF info
  # only valid for processes with a stationary mean
  SUM$DOF <- c( nant(DOF.mean(object),0) , DOF.area(object) , DOF.diffusion/2 , DOF.speed/2 )
  names(SUM$DOF) <- c("mean","area","diffusion","speed")

  SUM$CI <- par

  return(SUM)
}
#methods::setMethod("summary",signature(object="ctmm"), function(object,...) summary.ctmm(object,...))

#DOF of area
DOF.area <- function(CTMM)
{
  sigma <- CTMM$sigma@par
  AREA <- area.covm(CTMM$sigma)

  if("major" %nin% dimnames(CTMM$COV)[[1]])
  { VAR <- Inf }
  else if(CTMM$isotropic[1])
  { VAR <- CTMM$COV['major','major'] }
  else
  {
    P <- c('major','minor')
    GRAD <- AREA / sigma[P] / 2
    VAR <- c( GRAD %*% CTMM$COV[P,P] %*% GRAD )
  }

  if(CTMM$range) { DOF <- AREA^2/abs(VAR) } else { DOF <- 0 }

  if(is.nan(DOF)) { DOF <- 0 }

  return(DOF)
}


# DOF of trace variance
DOF.var <- function(CTMM)
{
  MSD <- var.covm(CTMM$sigma,ave=FALSE)
  COV.RAN <- axes2var(CTMM,MEAN=FALSE)

  2*MSD^2/COV.RAN
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
