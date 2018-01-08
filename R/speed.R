# new generic
speed <- function(object,...) UseMethod("speed")

speed.telemetry <- function(object,CTMM,level=0.95,prior=TRUE,fast=TRUE,error=0.01,mc.cores=NULL,...)
{ speed.ctmm(CTMM,data=object,level=level,prior=prior,error=error,mc.cores=mc.cores,...) }

speed.ctmm <- function(object,data=NULL,level=0.95,prior=TRUE,fast=TRUE,error=0.01,mc.cores=NULL,...)
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
  ERROR.MEAN <- Inf
  ERROR.SE <- Inf
  while(ERROR.MEAN>=error || ERROR.SE>=error || length(SPEEDS)<=10)
  {
    SPEEDS <- c(SPEEDS,unlist(mclapply(1:mc.cores,spd.fn,mc.cores=mc.cores)))

    MEAN <- mean(SPEEDS)

    if(length(SPEEDS)>1) { VAR <- stats::var(SPEEDS) }
    else { VAR <- Inf }

    # standard error(mean) / mean
    ERROR.MEAN <- sqrt(VAR/length(SPEEDS)) / MEAN

    # standard error(standard deviation) / mean
    ERROR.SE <- 1/2 * sqrt(2/(length(SPEEDS)-1)*VAR) / MEAN
  }

  # calculate 95% CIs (log-normal)
  SPEEDS <- log(SPEEDS)
  MEAN <- mean(SPEEDS)
  VAR <- stats::var(SPEEDS)
  CI <- norm.ci(MEAN,VAR,alpha=1-level)
  CI <- exp(CI)

  UNITS <- unit(CI,"speed")
  CI <- rbind(CI)/UNITS$scale
  rownames(CI) <- paste0("speed (",UNITS$name,")")

  return(CI)
}

# calculate speed of one random trajectory
speed.rand <- function(CTMM,data=NULL,prior=TRUE,fast=TRUE,error=0.01,precompute=FALSE,...)
{
  # capture model uncertainty
  if(prior) {  CTMM <- emulate(CTMM,data=data,fast=fast,...) }

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
