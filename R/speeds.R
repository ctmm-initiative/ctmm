# predicted speeds over time
# what is the easiest way to include parameter uncertainty?
# add independent variance from RMS speed?

####
speeds.telemetry <- function(object,CTMM,t=NULL,cycle=Inf,level=0.95,robust=FALSE,prior=FALSE,fast=TRUE,error=0.01,cores=1,...)
{
  data <- object
  # check conflicting conditions

  if(is.null(t)) { t <- data$t }

  if(!prior && fast) { SPEEDS <- speeds.fast(data,CTMM=CTMM,t=t,cycle=cycle,level=level,robust=robust,...) }
  else { SPEEDS <- speeds.slow(data,CTMM=CTMM,t=t,cycle=cycle,level=level,robust=robust,prior=prior,fast=fast,error=error,cores=cores,...) }

  SPEEDS <- as.data.frame(SPEEDS)
  SPEEDS$t <- t
  # include timestamps if possible
  timezone <- attr(data,'info')$timezone
  if(!is.null(timezone))
  { SPEEDS$timestamp <- as.POSIXct(SPEEDS$t,tz=timezone,origin=EPOCH) }

  return(SPEEDS)
}

speeds.ctmm <- function(object,data=NULL,t=NULL,cycle=Inf,level=0.95,robust=FALSE,prior=FALSE,fast=TRUE,error=0.01,cores=1,...)
{ speeds.telemetry(data,CTMM=object,t=t,cycle=cycle,level=level,robust=robust,prior=prior,fast=fast,error=error,cores=cores,...) }


# emulate then simulate
speeds.slow <- function(data,CTMM=NULL,t=NULL,level=0.95,robust=FALSE,prior=FALSE,fast=TRUE,error=0.01,cores=1,...)
{
  n <- length(t)

  # function to evaluate random speeds
  SUB <- NULL
  spds.fn <- function(i=0)
  {
    # capture model uncertainty
    if(prior) { CTMM <- emulate(CTMM,data=data,fast=fast,...) }
    # fail state for fractal process
    if(length(CTMM$tau)<2 || CTMM$tau[2]<=.Machine$double.eps) { return(rep(Inf,n)) }
    if(CTMM$tau[2]==Inf) { return(rep(0,n)) }

    data <- simulate(CTMM,data=data,t=t,precompute=precompute)

    if(nrow(data)>n)
    {
      if(is.null(SUB)) { SUB <<- (data$t %in% t) }
      data <- data[SUB,]
    }

    # instantaneous speeds
    data <- sqrt(data$vx^2+data$vy^2)

    return(data)
  }

  # bad return value
  INF <- array(c(0,Inf,Inf),c(3,n))
  INF <- t(INF)
  colnames(INF) <- c("low","ML","high")

  if(!robust)
  {
    S1 <- array(0,n)
    S2 <- array(0,n)
  }

  if(!prior) # precompute kalman filter
  {
    precompute <- TRUE
    SPEEDS <- spds.fn()
    dim(SPEEDS) <- c(n,1)
    N <- 1
    precompute <- -1

    S1 <- rowSums(SPEEDS)
    S2 <- rowSums(SPEEDS^2)

    if(any(SPEEDS==Inf))
    {
      warning("Sampling distribution does not always resolve velocity. Try robust=TRUE.")
      return(INF)
    }
  }
  else
  {
    SPEEDS <- NULL
    N <- 0
    precompute <- FALSE
  }

  # loop over emulations
  ERROR <- Inf
  pb <- utils::txtProgressBar(style=3)
  while(ERROR>=error || N<=20)
  {
    ADD <- plapply(1:cores,spds.fn,cores=cores,fast=FALSE)
    ADD <- simplify2array(ADD) # (t,cores)
    dim(ADD) <- c(n,cores)
    SPEEDS <- cbind(SPEEDS,ADD) # (t,N)
    N <- N + cores

    if(!robust)
    {
      S1 <- S1 + rowSums(ADD)
      S2 <- S2 + rowSums(ADD^2)

      if(any(ADD==Inf))
      {
        warning("Sampling distribution does not always resolve velocity. Try robust=TRUE.")
        return(INF)
      }
    }
    else # use insert sort to keep SPEEDS sorted
    {
      SPEEDS <- vapply(1:nrow(SPEEDS),function(i){sort(SPEEDS[i,],method="radix")},SPEEDS[1,]) # (N,t)
      dim(SPEEDS) <- c(N,n) # R drops indices for YOUR CONVENIENCE :)
      SPEEDS <- t(SPEEDS) # (t,N)
    }

    if(N>1)
    {
      if(!robust)
      {
        AVE <- S1/N
        VAR <- abs(S2 - N*AVE^2)/(N-1)
        ERROR <- sqrt(VAR/N) / AVE
        ERROR <- max(ERROR)
      }
      else
      {
        # median
        AVE <- mint(SPEEDS,(N+1)/2)
        # standard error on the median
        Q1 <- mint(SPEEDS,(N+1-sqrt(N))/2)
        Q2 <- mint(SPEEDS,(N+1+sqrt(N))/2)
        ERROR <- max(1-Q1/AVE,Q2/AVE-1)
        # correct for Inf AVE
        if(is.nan(ERROR))
        {
          ERROR <- Inf
          warning("Speeds accumulating on boundary and expectation value may not converge.")
        }
      }

      # update progress bar
      utils::setTxtProgressBar(pb,clamp(min(length(SPEEDS)/20,(error/ERROR)^2)))
    }
  } # end while

  # return raw data (undocumented)
  if(is.na(level) || is.null(level)) { return(SPEEDS) }

  # calculate averages and variances
  if(!robust)
  {
    M1 <- rowMeans(SPEEDS) # (n)
    M2 <- rowMeans(SPEEDS^2)

    # chi^1 DOF consistent with 1-2 moments
    DOF <- vapply(1:n,function(i){chi.dof(M1[i],M2[i])},numeric(1))

    CI <- vapply(1:n,function(i){ chisq.ci(M2[i],DOF=DOF[i],level=level) },numeric(3)) # (3,t)
    CI <- sqrt(CI)
    CI[2,] <- M1
  }
  else ### start here !!!
  {
    alpha <- (1-level)/2
    CI <- vapply(1:nrow(SPEEDS),function(i){stats::quantile(SPEEDS[i,],c(alpha,0.5,1-alpha))},numeric(3)) # (3,t)
    rownames(CI) <- c("low","ML","high")
  }
  CI <- t(CI)

  close(pb)
  return(CI)
}


####
speeds.fast <- function(data,CTMM=NULL,t=NULL,level=0.95,robust=FALSE,...)
{
  n <- length(t)

  DOF <- summary(CTMM)$DOF['speed']
  if(!DOF)
  {
    M1 <- rep(Inf,n) # upgrade to outlie estimates?
    M2 <- rep(Inf,n)
    DOF <- numeric(n)
  }
  else
  {
    if(!is.null(CTMM)) { data <- predict(CTMM,data=data,t=t,...) }

    axes <- c("vx","vy")
    v <- get.telemetry(data,axes=axes) # (n,2)
    VAR <- get.error(data,list(axes=axes,error=TRUE),DIM=2) # (n,2,2)

    # exact second moment - not used currently
    M2 <- vapply(1:n,function(i){ sum( v[i,]^2 + diag(VAR[i,,]) ) },numeric(1))

    # good approximation to mean speed (delta method fails when velocity estimate is small)
    M1 <- vapply(1:n,function(i){abs.bivar(v[i,],VAR[i,,])},numeric(1)) # n

    # variance of square speed # exact?
    # VAR <- 4 * vapply(1:n,function(i){v[i,] %*% VAR[i,,] %*% v[i,]},numeric(1)) # (n)
    # variance of speed - delta method (M2 might be inconsistent with M1)
    # VAR <- VAR/AVE^2

    # chi^1 DOF consistent with M1 & VAR about M1 (M2 might be inconsistent with M1)
    DOF <- vapply(1:n,function(i){chi.dof(M1[i],M2[i])},numeric(1))
  }


  # if no level, return point estimate and DOF
  if(is.null(level))
  { v <- cbind(speed=M1,DOF=DOF,VAR=pmax(0,M2-M1^2)) } # output v and chi DOF
  else
  {
    v <- vapply(1:n,function(i){ chisq.ci(M2[i],DOF=DOF[i],level=level,robust=robust) },numeric(3)) # (3,n)
    v <- sqrt(t(v)) # (n,3)
    v[,2] <- M1
  }

  return(v)
}

# approximate <|r|> of bi-variate Gaussian distribution
# bessel function stuff (exact for equal variance)
# elliptical functions stuff (exact for zero mean)
# anistotropic + nonzero combination is approximate
# returns first two moments
abs.bivar <- function(mu,Sigma)
{
  sigma0 <- mean(diag(Sigma))
  stdev0 <- sqrt(sigma0)
  mu2 <- sum(mu^2)
  mu <- sqrt(mu2)

  sigma <- eigen(Sigma)$values

  Barg <- mu2/(4*sigma0)
  if(Barg >= BESSEL_LIMIT) { return(mu) }

  B0 <- besselI(Barg,0,expon.scaled=TRUE)
  B1 <- besselI(Barg,1,expon.scaled=TRUE)

  # contains deterministic limit
  sqrtpi2 <- sqrt(pi/2)
  Bv <-  sqrtpi2 * sqrt(Barg) * ( B0 + B1 ) * mu
  # contains stochastic limit
  Bs <- B0 / sqrtpi2 * sqrt(sigma[1]) * pracma::ellipke(1-sigma[2]/sigma[1])$e

  M1 <- Bv + Bs

  return(M1)
}
