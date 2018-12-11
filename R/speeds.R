# predicted speeds over time
# what is the easiest way to include parameter uncertainty?
# add independent variance from RMS speed?

####
speeds.telemetry <- function(object,CTMM,t=NULL,level=0.95,robust=FALSE,prior=FALSE,fast=TRUE,error=0.01,cores=1,...)
{
  data <- object
  # check conflicting conditions

  # return Infs if OU/BM

  if(is.null(t)) { t <- data$t }

  if(!prior) { SPEEDS <- speeds.fast(data,CTMM=CTMM,t=t,level=level,robust=robust,...) }
  else { SPEEDS <- speeds.slow(data,CTMM=CTMM,t=t,level=level,robust=robust,fast=fast,error=error,cores=cores,...) }

  SPEEDS <- as.data.frame(SPEEDS)
  SPEEDS$t <- t
  # include timestamps if possible
  timezone <- attr(data,'info')$timezone
  if(!is.null(timezone))
  { SPEEDS$timestamp <- as.POSIXct(SPEEDS$t,tz=timezone,origin=EPOCH) }

  return(SPEEDS)
}

speeds.ctmm <- function(object,data,t=NULL,level=0.95,robust=FALSE,prior=FALSE,fast=TRUE,error=0.01,cores=1,...)
{ speeds.telemetry(data,CTMM=object,t=t,level=level,robust=robust,prior=prior,fast=fast,error=error,cores=cores,...) }


# emulate then simulate
speeds.slow <- function(data,CTMM=NULL,t=NULL,level=0.95,robust=FALSE,fast=TRUE,error=0.01,cores=1,...)
{
  n <- length(t)

  # function to evaluate random speeds
  SUB <- NULL
  spds.fn <- function(i)
  {
    # capture model uncertainty
    CTMM <- emulate(CTMM,data=data,fast=fast,...)
    # fail state for fractal process
    if(length(CTMM$tau)<2 || CTMM$tau[2]<=.Machine$double.eps) { return(rep(Inf,n)) }
    if(CTMM$tau[2]==Inf) { return(rep(0,n)) }

    data <- simulate(data,CTMM=CTMM,t=t)

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

  # loop over emulations
  N <- 0
  SPEEDS <- NULL
  if(!robust)
  {
    S1 <- array(0,n)
    S2 <- array(0,n)
  }
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
    AVE <- rowMeans(SPEEDS) # (n)
    # calculate MVU variance under log transformation
    VAR <- apply(log(SPEEDS),1,stats::var)
    # transform back
    VAR <- VAR * AVE^2
    # chi^2 speed^2
    DOF <- 2*AVE^2/VAR

    CI <- sapply(1:n,function(i){ chisq.ci(AVE[i]^2,DOF=DOF[i],level=level) }) # (t,3)
    CI <- sqrt(CI)
  }
  else ### start here !!!
  {
    alpha <- (1-level)/2
    CI <- sapply(1:nrow(SPEEDS),function(i){stats::quantile(SPEEDS[i,],c(alpha,0.5,1-alpha))}) # (3,t)
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
    AVE <- rep(Inf,n) # upgrade to outlie estimates?
    DOF <- numeric(n)
  }
  else
  {
    if(!is.null(CTMM)) { data <- predict(data,CTMM=CTMM,t=t,...) }

    axes <- c("vx","vy")
    v <- get.telemetry(data,axes=axes) # (n,2)
    VAR <- get.error(data,list(axes=axes,error=TRUE),DIM=2) # (n,2,2)

    # delta method doesn't work very well for small speed estimates
    # # delta method - variance of square speed
    # VAR <- 4 * vapply(1:n,function(i){v[i,] %*% VAR[i,,] %*% v[i,]},numeric(1)) # (n)
    # # speeds
    # v2 <- rowSums(v^2)
    # DOF <- 2*v2^2/VAR

    AVE <- numeric(n)
    DOF <- numeric(n)
    for(i in 1:n)
    {
      EIGEN <- eigen(VAR[i,,])
      CHI2.VAR <- sum(EIGEN$values) # variance of velocity
      AVE[i] <- sum( (t(EIGEN$vectors) %*% v[i,])^2 ) + CHI2.VAR # mean square velocity
      CHI2.VAR <- 2*sum(EIGEN$values^2) # variance of square speed in 2D
      DOF[i] <- 2*AVE[i]^2/CHI2.VAR # chi^2 DOF
    }
  }


  # if no level, return point estimate and DOF
  if(is.null(level))
  {
    # output v and chi DOF
    v <- sqrt(AVE)
    v <- cbind(v,DOF)
    colnames(v) <- c("speed","DOF")
  }
  else
  {
    v <- vapply(1:n,function(i){ chisq.ci(AVE[i],DOF=DOF[i],level=level,robust=robust) },numeric(3)) # (3,n)
    v <- sqrt(t(v))
  }

  return(v)
}
