speed.telemetry <- function(object,CTMM,t=NULL,level=0.95,robust=FALSE,units=TRUE,prior=TRUE,fast=TRUE,cor.min=0.5,dt.max=NULL,error=0.01,cores=1,trace=TRUE,...)
{ speed.ctmm(CTMM,data=object,t=t,level=level,robust=robust,units=units,prior=prior,fast=fast,cor.min=cor.min,dt.max=dt.max,error=error,cores=cores,trace=trace,...) }

speed.ctmm <- function(object,data=NULL,t=NULL,level=0.95,robust=FALSE,units=TRUE,prior=TRUE,fast=TRUE,cor.min=0.5,dt.max=NULL,error=0.01,cores=1,trace=TRUE,...)
{
  # bad return value
  INF <- c(0,Inf,Inf)
  names(INF) <- NAMES.CI
  INF <- rbind(INF)
  rownames(INF) <- "speed (meters/second)"
  INF <- list(DOF=0,CI=INF)
  names(INF$DOF) <- "speed"
  class(INF) <- "speed"

  if(length(object$tau)<2 || object$tau[2]<=.Machine$double.eps)
  {
    warning("Movement model is fractal.")
    return(INF)
  }

  if(prior && fast)
  {
    TEST <- try(any(eigen(object$COV)$values<=.Machine$double.eps),silent=TRUE)
    TEST <- class(TEST)[1]!="logical" || TEST
    if(TEST)
    {
      warning("Indefinite covariance matrix estimate. Consider fast=FALSE.")
      return(INF)
    }
  }

  if(is.null(data) && prior && !fast)
  { stop("Parametric boostrap cannot be performed without data.") }

  # analytically solvable cases
  if(!robust && is.null(data) && (!prior || fast))
  {
    if(object$isotropic[1] && object$mean!="stationary") # chi_2 : circular velocity distribution
    {
      CI <- summary(object,level=level,units=FALSE)
      DOF <- CI$DOF['speed']
      CI <- CI$CI['speed (meters/second)',] * sqrt(pi/2/2)
    }
    else # elliptical velocity distribution
    {
      # propagate error uncertainty
      UERE.FIT <- object$error>0 & ( paste("error",names(object$error)) %in% dimnames(object)[[1]] )
      STUFF <- id.parameters(object,profile=FALSE,linear=FALSE,UERE.FIT=UERE.FIT)
      NAMES <- STUFF$NAMES
      parscale <- STUFF$parscale
      lower <- STUFF$lower
      upper <- STUFF$upper

      fn <- function(p)
      {
        object <- set.parameters(object,p)
        speed_deterministic(object)
      }

      PAR <- get.parameters(object,NAMES)
      MEAN <- fn(PAR)

      if(prior && fast)
      {
        GRAD <- genD(par=PAR,fn=fn,lower=lower,upper=upper,parscale=parscale,Richardson=2,mc.cores=1)$grad
        # GRAD <- numDeriv::grad(fn,PAR) # don't step outside of (0,1)
        if("COV" %in% names(object))
        { VAR <- c(GRAD %*% object$COV[NAMES,NAMES] %*% GRAD) } # emulate can lose features
        else
        { VAR <- Inf }

        if(object$mean!="stationary")
        {
          fn <- function(p)
          {
            object$mu[] <- p
            speed_deterministic(object)
          }

          GRAD <- genD(par=c(object$mu),fn=fn,Richardson=2,mc.cores=1)$grad
          VAR2 <- aperm(object$COV.mu,c(2,1,3,4))
          dim(VAR2) <- c(1,1)*length(GRAD)
          VAR2 <- c(GRAD %*% VAR2 %*% GRAD)
          VAR <- VAR + VAR2
        }

        # propagate errors (chi)
        M2 <- VAR + MEAN^2
        DOF <- chi.dof(MEAN,M2) # chi DOF
        CI <- sqrt(chisq.ci(M2,DOF=DOF,alpha=1-level)) # chi=sqrt(chi^2) CI
        CI <- CI / chi.bias(DOF)
        CI[2] <- MEAN # shouldn't be necessary?
      }
      else # ! prior
      {
        CI <- c(1,1,1)*MEAN
        names(CI) <- NAMES.CI
        DOF <- 0
      }
    }
  }
  else # simulation based evaluation
  {
    # pilot estimate of DOF
    DOF <- 2*speed(object)$DOF

    # time range to consider
    if(is.null(t)) { t <- data$t }
    t <- t[c(1,length(t))]

    cores <- resolveCores(cores,fast=FALSE)

    # random speed calculation
    DT <- diff(data$t)
    spd.fn <- function(i=0) { speed_rand(object,data=data,prior=prior,fast=fast,cor.min=cor.min,dt.max=dt.max,error=error,precompute=precompute,DT=DT,TP=t,DOF=DOF,...) }

    # setup precompute stuff
    if(prior==FALSE)
    {
      precompute <- TRUE # precompute matrices
      SPEEDS <- spd.fn()
      precompute <- -1 # from now on use precomputed matrices
    }
    else
    {
      precompute <- FALSE
      SPEEDS <- numeric(0)
    }

    # keep replicating until error target
    if(trace) { pb <- utils::txtProgressBar(style=3) }
    ERROR <- Inf
    N <- length(SPEEDS)
    if(!robust)
    {
      S1 <- sum(SPEEDS)
      S2 <- sum(SPEEDS^2)
    }
    while(ERROR>=error || length(SPEEDS)<=20)
    {
      ADD <- unlist(plapply(1:cores,spd.fn,cores=cores,fast=FALSE))
      SPEEDS <- c(SPEEDS,ADD)
      # rolling mean & variance
      N <- N + cores
      if(!robust) # rolling sums to keep complexity O(N)
      {
        if(any(ADD==Inf))
        {
          warning("Sampling distribution does not always resolve velocity. Try robust=TRUE.")
          return(INF)
        }
        S1 <- S1 + sum(ADD)
        S2 <- S2 + sum(ADD^2)
      }
      else # use insert sort to keep SPEEDS sorted
      { SPEEDS <- sort(SPEEDS,method="radix") }

      # standard_error(mean) / mean
      if(N>1)
      {
        # calculate averages and their standard errors
        if(!robust)
        {
          AVE <- S1/N
          VAR <- abs(S2 - N*AVE^2)/(N-1)
          ERROR <- nant( sqrt(VAR/N) / AVE ,Inf)
        }
        else
        {
          # median
          AVE <- vint(SPEEDS,(N+1)/2)
          # standard error on the median
          Q1 <- vint(SPEEDS,(N+1-sqrt(N))/2)
          Q2 <- vint(SPEEDS,(N+1+sqrt(N))/2)
          ERROR <- nant( max(AVE-Q1,Q2-AVE) / AVE ,Inf)

          if(N>1/error^2)
          {
            warning("Expectation values did not converge after ",N," iterations.")
            break
          }
        }

        # update progress bar
        if(trace) { utils::setTxtProgressBar(pb,clamp(min(length(SPEEDS)/20,(error/ERROR)^2))) }
      } # end N>1 ERROR calc
    } # end while ERROR

    # return raw data (undocumented)
    if(is.na(level) || is.null(level)) { return(SPEEDS) }

    # calculate averages and variances
    if(!robust)
    {
      # more correct chi^1 statistics
      M1 <- mean(SPEEDS)
      M2 <- mean(SPEEDS^2)
      DOF <- chi.dof(M1,M2)
      CI <- chisq.ci(M2,DOF=DOF,alpha=1-level) # correct mean and variance
      CI <- sqrt(CI)
      CI[2] <- M1
    }
    else
    {
      alpha <- (1-level)/2
      CI <- stats::quantile(SPEEDS,c(alpha,0.5,1-alpha))
      names(CI) <- NAMES.CI

      M1 <- AVE
      M2 <- ((Q2-Q1)/2)^2 + M1^2
      DOF <- chi.dof(M1,M2)
    }

    if(trace) { close(pb) }
  } # end simulations

  UNITS <- unit(CI,"speed",SI=!units)
  CI <- rbind(CI)/UNITS$scale
  rownames(CI) <- paste0("speed (",UNITS$name,")")
  #attr(CI,"DOF") <- DOF
  if(CI[1]==Inf) { CI[1] <- 0 } # sampled all Inf
  CI <- list(DOF=DOF/2,CI=CI)
  names(CI$DOF) <- "speed"
  class(CI) <- "speed"
  return(CI)
}


# calculate speed of one random trajectory
speed_rand <- function(CTMM,data=NULL,prior=TRUE,fast=TRUE,cor.min=0.5,dt.max=NULL,error=0.01,precompute=FALSE,DT=diff(data$t),TP=range(data$t),DOF=Inf,...)
{
  # capture model uncertainty
  if(prior) { CTMM <- emulate(CTMM,data=data,fast=fast,...) }
  # fail state for fractal process
  if(length(CTMM$tau)<2 || CTMM$tau[2]<=.Machine$double.eps) { return(Inf) }
  if(CTMM$tau[2]==Inf) { return(0) }

  dt <- CTMM$tau[2]*(error/10)^(1/3) # this gives O(error/10) instantaneous error

  if(is.null(data))
  {
    # analytic result possible here
    if(CTMM$mean=="stationary") { return(speed_deterministic(CTMM)/chi.bias(DOF)) }
    # else do a sufficient length simulation
    t <- seq(0,CTMM$tau[2]/error^2,dt) # this should give O(error) estimation error
    data <- simulate(CTMM,t=t,precompute=precompute)
  }
  else
  {
    # check for rare cases in sampling disribution where motion is almost but not fractal relative to sampling schedule
    # first check if non-stationary mean is irrelevant
    MSV <- sum(diag(CTMM$sigma))/ ifelse(CTMM$range,prod(CTMM$tau),CTMM$tau[2])
    if(nant(drift.speed(CTMM)$EST/MSV,0)<error^2)
    {
      # now check if there are many steps per sampled interval
      FRAC <- CTMM$tau[2]/DT
      if(all(FRAC<error))
      {
        # finally check if the simulated distance would be much greater than sampled net displacement
        FAKE <- CTMM
        FAKE$mean <- "stationary"
        SPD <- speed(FAKE,prior=FALSE)[2]
        FRAC <- sqrt(diff(data$x)^2+diff(data$y)^2)/DT / SPD
        if(all(FRAC<error)) { return(SPD) }
      }
    }
    if(MSV==0) # need better case handling for tiny sigma estimate?
    { return(0) }

    if(is.null(dt.max)) { dt.max <- -log(cor.min)*CTMM$tau[2] }
    data <- simulate(CTMM,data=data,dt=dt,precompute=precompute,dt.max=dt.max,DT=DT)
    t <- data$t
  }
  if(is.null(dt.max)) { dt.max <- Inf }

  # instantaneous speeds
  data <- sqrt(data$vx^2+data$vy^2)

  # weights for Riemmann sum
  w <- diff(t)
  w <- w * (w<=dt.max) # skip gaps
  w <- c(0,w) + c(w,0)
  w <- w * (t>=TP[1] & t<=TP[2]) # only consider time-interval of interest

  # weighted average speed
  v <- sum(w*data)/sum(w)

  return(v)
}


# only for stationary processes
speed_deterministic <- function(CTMM,sigma=CTMM$sigma)
{
  if(CTMM$range)
  { sigma <- sigma/prod(CTMM$tau) } # OUF
  else
  { sigma <- sigma/CTMM$tau[2] } # IOU

  if(CTMM$mean!="stationary") # this is approximate
  {
    E <- drift.energy(CTMM)$VV
    E <- t(CTMM$mu) %*% E %*% CTMM$mu
    sigma <- sigma + E
  }

  sigma <- eigen(sigma)$values

  if(CTMM$isotropic || sigma[1]==sigma[2])
  { v <- sqrt(sigma[1] * pi/2) }
  else
  { v <- sqrt(2/pi) * sqrt(sigma[1]) * pracma::ellipke(1-clamp(sigma[2]/sigma[1]))$e }

  return(v)
}
