speed.telemetry <- function(object,CTMM,level=0.95,prior=TRUE,fast=TRUE,cor.min=0.5,dt.max=NULL,error=0.01,cores=1,...)
{ speed.ctmm(CTMM,data=object,level=level,prior=prior,error=error,cores=cores,...) }

speed.ctmm <- function(object,data=NULL,level=0.95,prior=TRUE,fast=TRUE,cor.min=0.5,dt.max=NULL,error=0.01,cores=1,...)
{
  if(length(object$tau)<2 || object$tau[2]<=.Machine$double.eps)
  { stop("Movement model is fractal. Speed cannot be estimated.") }

  if(prior && fast && any(eigen(object$COV,only.values=TRUE)$values<=.Machine$double.eps))
  { stop("Indefinite covariance matrix in sampling distribution.") }

  # analytically solvable cases
  if(is.null(data) && object$mean=="stationary" && (!prior || fast))
  {
    if(object$isotropic) # chi_2 : circular velocity distribution
    { CI <- summary(object,level=level,units=FALSE)$CI['speed (meters/second)',] * sqrt(pi/2/2) }
    else # elliptical velocity distribution
    {
      NAMES <- id.parameters(object,profile=FALSE)$NAMES

      fn <- function(p)
      {
        CTMM <- set.parameters(object,p)
        speed.deterministic(CTMM)
      }

      PAR <- get.parameters(object,NAMES)
      MEAN <- fn(PAR)

      if(prior && fast)
      {
        GRAD <- numDeriv::grad(fn,PAR) # don't step outside of (0,1)
        VAR <- c(GRAD %*% object$COV %*% GRAD)

        # propagate errors chi -> chi^2
        CI <- chisq.ci(MEAN^2,(2*MEAN)^2*VAR,alpha=1-level)
        # transform back
        CI <- sqrt(CI)
      }
      else # ! prior
      {
        CI <- c(1,1,1)*MEAN
        names(CI) <- c("low","ML","high")
      }
    }
  }
  else # simulation based evaluation
  {
    cores <- resolveCores(cores,fast=FALSE)

    # random speed calculation
    spd.fn <- function(i=0) { speed.rand(object,data=data,prior=prior,fast=fast,cor.min=cor.min,dt.max=dt.max,error=error,precompute=precompute,...) }

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
    pb <- utils::txtProgressBar(style=3)
    ERROR <- Inf
    N <- length(SPEEDS)
    S1 <- sum(SPEEDS)
    S2 <- sum(SPEEDS^2)
    while(ERROR>=error || length(SPEEDS)<=20)
    {
      ADD <- unlist(plapply(1:cores,spd.fn,cores=cores,fast=FALSE))
      SPEEDS <- c(SPEEDS,ADD)
      # rolling mean & variance
      N <- N + cores
      S1 <- S1 + sum(ADD)
      S2 <- S2 + sum(ADD^2)

      # standard_error(mean) / mean
      if(length(SPEEDS)>1)
      {
        MEAN <- S1/N
        VAR <- abs(S2 - N*MEAN^2)/(N-1)
        ERROR <- sqrt(VAR/N) / MEAN

        # fallback on quantiles if some Inf speeds
        if((is.na(level) || is.null(level)))
        {
          SPEEDS <- sort(SPEEDS,method='quick')
          # median
          M <- SPEEDS[round(N/2)]
          # standard error on the median
          Q1 <- SPEEDS[round((N-sqrt(N))/2)]
          Q2 <- SPEEDS[round((1+N+sqrt(N))/2)]
          ERROR <- max(M-Q1,Q2-M) / M
        }

        # update progress bar
        utils::setTxtProgressBar(pb,clamp(min(length(SPEEDS)/20,(error/ERROR)^2)))
      }
    }

    # return raw data (undocumented)
    if(is.na(level) || is.null(level)) { return(SPEEDS) }

    MEAN <- mean(SPEEDS)
    # calculate variance under log transform, where data is more normal
    VAR <- stats::var(log(SPEEDS))
    # transform back
    VAR <- MEAN^2 * VAR
    # chi^2 square speed
    CI <- chisq.ci(MEAN^2,(2*MEAN)^2*VAR,alpha=1-level)
    CI <- sqrt(CI)

    close(pb)
  }

  UNITS <- unit(CI,"speed")
  CI <- rbind(CI)/UNITS$scale
  rownames(CI) <- paste0("speed (",UNITS$name,")")

  return(CI)
}

# calculate speed of one random trajectory
speed.rand <- function(CTMM,data=NULL,prior=TRUE,fast=TRUE,cor.min=0.5,dt.max=NULL,error=0.01,precompute=FALSE,...)
{
  # capture model uncertainty
  if(prior) { CTMM <- emulate(CTMM,data=data,fast=fast,...) }
  # fail state for fractal process
  if(length(CTMM$tau)==1 || CTMM$tau[2]<=.Machine$double.eps) { return(Inf) }

  dt <- CTMM$tau[2]*sqrt(error) # how thin should this be, exactly?

  if(is.null(data))
  {
    # analytic result possible here
    if(data$mean=="stationary") { return(speed.deterministic(CTMM)) }
    # else do a sufficient length simulation
    t <- seq(0,CTMM$tau[2]/error^2,dt)
    data <- simulate(CTMM,t=t,precompute=precompute)
  }
  else
  {
    # check for rare cases in sampling disribution where motion is almost but not fractal relative to sampling schedule
    if(prior && !fast)
    {
      # first check if non-stationary mean is irrelevant
      drift <- get(CTMM$mean)
      if(drift@speed(CTMM)$EST/(sum(diag(CTMM$sigma))/prod(CTMM$tau))<error^2)
      {
        # now check if there are many steps per sampled interval
        DT <- diff(data$t)
        FRAC <- CTMM$tau[2]/DT
        if(all(FRAC<error))
        {
          # finally check if the simulated distance would be much greater than sampled net displacement
          FAKE <- CTMM
          FAKE$mean <- "stationary"
          SPD <- speed(FAKE)[2]
          FRAC <- sqrt(diff(data$x)^2+diff(data$y)^2)/DT / SPD
          if(all(FRAC<error)) { return(SPD) }
        }
      }
    }

    if(cor.min) { dt.max <- -log(cor.min)*CTMM$tau[2] }
    data <- simulate(CTMM,data=data,dt=dt,precompute=precompute,dt.max=dt.max)
    t <- data$t
  }
  if(is.null(dt.max)) { dt.max <- Inf }

  # instantaneous speeds
  data <- sqrt(data$v.x^2+data$v.y^2)

  # weights for Riemmann sum
  w <- diff(t)
  w <- w * (w<=dt.max) # skip gaps
  w <- c(0,w) + c(w,0)

  # weighted average speed
  v <- sum(w*data)/sum(w)

  return(v)
}

# only for stationary processes
speed.deterministic <- function(CTMM)
{
  if(CTMM$isotropic)
  { v <- sqrt(sum(diag(CTMM$sigma))/prod(CTMM$tau) * pi/2/2) }
  else
  {
    sigma <- attr(CTMM$sigma,"par")
    # eigen values of velocity variance
    sigma <- sigma['area'] / prod(CTMM$tau) * exp(c(1,-1)/2*sigma['eccentricity'])
    v <- sqrt(2/pi) * sqrt(sigma[1]) * pracma::ellipke(-diff(sigma)/sigma[1])$e
  }
  return(v)
}
