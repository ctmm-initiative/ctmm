# R refuses to call log.list for list objects outside of ctmm
Log <- function(x,variable="area",debias=TRUE,...)
{
  x <- listify(x)
  NAMES <- names(x)
  x <- log.list(x,variable=variable,debias=debias,...)

  x <- data.frame(log=x$log,VAR.log=x$VAR.log)
  rownames(x) <- NAMES

  INF <- x$VAR.log==Inf
  x$log[INF] <- 0

  return(x)
}


log_ctmms <- function(x,variable="area",debias=TRUE,level.UD=0.95,...)
{
  x <- listify(x)
  x <- import.variable(x,variable=variable,level.UD=level.UD)
  # list(ID=ID,AREA=AREA,DOF=DOF,variable=variable)
  y <- list()
  y$log <- log(x$AREA)

  # 2D correction made in import.variables
  if(!debias)
  { y$VAR.log <- 2/x$DOF }
  else
  {
    y$VAR.log <- trigamma(x$DOF/2)

    BIAS <- log_chi2_bias(x$DOF)
    y$log <- y$log - BIAS
  }

  return(y)
}


log_area <- function(x,variable="area",debias=TRUE,...)
{ log_ctmms(x,variable=variable,debias=debias,...) }


log_UD <- function(x,variable="area",debias=TRUE,level.UD=0.95,...)
{
  x <- listify(x)
  x <- lapply(x,function(y){summary(y,level.UD=level.UD,units=FALSE)})
  x <- log_area(x,debias=debias,...)
  return(x)
}


# speed
log_speed <- function(x,variable="speed",debias=TRUE,...)
{
  x <- listify(x)
  x <- import.variable(x,variable="speed",chi=TRUE)
  # list(ID=ID,AREA=AREA,DOF=DOF,variable=variable)
  y <- list()
  y$log <- log(x$AREA)

  # 2D correction made in import.variables
  if(!debias)
  { y$VAR.log <- 2/x$DOF/4 }
  else
  {
    y$VAR.log <- trigamma(x$DOF/2)
    y$VAR.log <- y$VAR.log/4

    BIAS <- digamma(x$DOF/2) - log(x$DOF/2)
    BIAS <- BIAS/2
    BIAS <- nant(BIAS,0)
    y$log <- y$log - BIAS
  }

  return(y)
}


Exp <- function(est,VAR.est=0,VAR=0,VAR.VAR=0,variable="area",debias=TRUE,level=0.95,units=TRUE,...)
{
  # convert from log-chi to log-chi^2
  R <- exp_log(est=est,VAR.est=VAR.est,VAR=VAR,VAR.VAR=VAR.VAR,...)
  est <- R$mu
  VAR <- R$VAR

  if(variable=="speed")
  { DOF <- chi.dof(est,est^2+VAR) }
  else # chi^2
  { DOF <- 2*est^2/VAR }

  if(debias)
  {
    BIAS <- digamma(DOF/2) - log(DOF/2)
    BIAS <- nant(BIAS,0)
    if(variable=="speed") { BIAS <- BIAS/2 }
    est <- est + BIAS
  }

  CI <- chisq.ci(est,VAR=VAR,level=level)

  # apply units and name
  CI <- summary_UD_format(CI,DOF/2,units=units)

  return(CI)
}


exp_log <- function(est,VAR.est=0,VAR=0,VAR.VAR=0,...)
{
  mu <- exp(est + VAR/2)
  # grad <- c(1,1/2) %o% mu
  VAR.mu <- (mu)^2*VAR.est + (mu/2)^2*VAR.VAR

  R <- list(mu=mu,VAR=VAR.mu)
  return(R)
}
