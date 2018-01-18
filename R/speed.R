speed.telemetry <- function(object,CTMM,level=0.95,prior=TRUE,fast=TRUE,error=0.01,mc.cores=NULL,...)
{ speed.ctmm(CTMM,data=object,level=level,prior=prior,error=error,mc.cores=mc.cores,...) }

speed.ctmm <- function(object,data=NULL,level=0.95,prior=TRUE,fast=TRUE,error=0.01,mc.cores=NULL,...)
{
  if(length(object$tau)<2 || object$tau[2]<=.Machine$double.eps)
  { stop("Movement model is fractal. Speed cannot be estimated.") }

  if(prior && fast && any(eigen(object$COV,only.values=TRUE)$values<=.Machine$double.eps))
  { stop("Indefinite covariance matrix in sampling distribution.") }

  # analytically solvable cases
  if(is.null(data) && object$mean=="stationary")
  {
    if(object$isotropic) # chi_2 : circular velocity distribution
    { CI <- summary(object,level=level,units=FALSE)$CI['speed (meters/second)',] * sqrt(pi/2/2) }
    else # elliptical velocity distribution
    {
      NAMES <- id.parameters(object,profile=FALSE)$NAMES

      fn <- function(p)
      {
        CTMM <- set.parameters(object,p)
        sigma <- attr(CTMM$sigma,"par")
        # eigen values of velocity variance
        sigma <- sigma['area'] / p["tau position"] / p["tau velocity"] * exp(c(1,-1)/2*sigma['eccentricity'])
        sqrt(2/pi) * sqrt(sigma[1]) * pracma::ellipke(-diff(sigma)/sigma[1])$e
      }

      PAR <- get.parameters(object,NAMES)
      GRAD <- numDeriv::grad(fn,PAR) # don't step outside of (0,1)
      VAR <- c(GRAD %*% object$COV %*% GRAD)
      MEAN <- fn(PAR)

      # propagate errors chi -> chi^2
      CI <- chisq.ci(MEAN^2,(2*MEAN)^2*VAR,alpha=1-level)
      # transform back
      CI <- sqrt(CI)
    }
  }
  else # simulation based evaluation
  {
    if(is.null(mc.cores)) { mc.cores <- detectCores() }

    # random speed calculation
    spd.fn <- function(i) { speed.rand(object,data=data,prior=prior,fast=fast,error=error,precompute=precompute,...) }

    # setup precompute stuff
    if(prior==FALSE)
    {
      precompute <- TRUE # precompute matrices
      SPEEDS <- spd.fn(0)
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
    MAX <- 0
    N <- length(SPEEDS)
    S1 <- sum(SPEEDS)
    S2 <- sum(SPEEDS^2)
    while(ERROR>=error || length(SPEEDS)<=10)
    {
      ADD <- unlist(mclapply(1:mc.cores,spd.fn,mc.cores=mc.cores))
      SPEEDS <- c(SPEEDS,ADD)
      # rolling mean & variance
      N <- N + mc.cores
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
        MAX <- c(MAX,ERROR,error)
        MAX <- max(MAX[MAX<Inf])
        if(error<MAX) { PROG <- clamp(log(MAX/ERROR)/log(MAX/error)) } else { PROG <- 1 }
        utils::setTxtProgressBar(pb,PROG)
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
speed.rand <- function(CTMM,data=NULL,prior=TRUE,fast=TRUE,error=0.01,precompute=FALSE,...)
{
  # capture model uncertainty
  if(prior) { CTMM <- emulate(CTMM,data=data,fast=fast,...) }
  # fail state for fractal process
  if(length(CTMM$tau)==1 || CTMM$tau[2]<=.Machine$double.eps) { return(Inf) }

  # SIMPSON'S RULE MIGHT NOT BE CORRECT FOR SQRT ERRORS
  # time range
  # t <- data$t
  # # how long should this be to minimize within-track uncertainty?
  # if(is.null(data)) { t <- c(0,CTMM$tau[2]/error^2) }
  #
  # dt <- CTMM$tau[2]/5 # how thin should this be?
  # RANGE <- last(t)-t[1]
  # n <- RANGE/dt
  # n <- ceiling(n)
  # if(is.even(n)) { n <- n + 1 }
  # dt <- RANGE/n
  #
  # t <- seq(t[1],last(t),length.out=n)
  #
  # data <- simulate(CTMM,data=data,t=t)
  #
  # # pull out just the times we are interested in (not the data times)
  # data <- data[data$t %in% t,]
  #
  # v <- sqrt( data$v.x^2 + data$v.y^2 )
  # rm(data)
  #
  # # simpson weights
  # w <- 2 + 2*is.even(1:n) # 2,4,2,...,2,4,2
  # w[1] <- w[n] <- 1
  #
  # v <- (dt/3) * sum(w*v) / RANGE

  dt <- CTMM$tau[2]/10 # how thin should this be?

  if(is.null(data))
  {
    t <- seq(0,CTMM$tau[2]/error^2,dt)
    data <- simulate(CTMM,t=t,precompute=precompute)
  }
  else
  {
    data <- simulate(CTMM,data=data,dt=dt,precompute=precompute)
    t <- data$t
  }

  # instantaneous speeds
  data <- sqrt(data$v.x^2+data$v.y^2)

  # weights for Riemmann sum
  w <- diff(t)
  w <- c(0,w) + c(w,0)

  # weighted average speed
  v <- sum(w*data)/sum(w)

  return(v)
}
