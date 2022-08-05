# R refuses to call log.list for list objects outside of ctmm
Log <- function(x,debias=TRUE,...)
{
  x <- listify(x)
  NAMES <- names(x)
  x <- log.list(x,debias=debias,...)
  names(x$log) <- names(x$VAR.log) <- NAMES
  return(x)
}


Exp <- function(est,VAR.est=0,VAR=0,VAR.VAR=0,variable="area",debias=TRUE,level=0.95,units=TRUE,...)
{
  fn <- paste0("exp.log.",variable)
  fn <- get(fn)

  fn(est,VAR.est=VAR.est,VAR=VAR,VAR.VAR=VAR.VAR,debias=debias,level=level,units=units,...)
}

exp.log <- function(est,VAR.est=0,VAR=0,VAR.VAR=0,...)
{
  mu <- exp(est + VAR/2)
  grad <- c(1,1/2) * mu
  VAR.mu <- sum( grad^2 * c(VAR.est,VAR.VAR) )

  R <- list(mu=mu,VAR=VAR.mu)
  return(R)
}


log.UD <- function(x,debias=TRUE,level.UD=0.95,...)
{
  x <- listify(x)
  x <- lapply(x,function(y){summary(y,level.UD=level.UD,units=FALSE)})
  x <- log.area(x,debias=debias,...)
  return(x)
}


log.area <- function(x,debias=TRUE,...)
{
  x <- listify(x)
  x <- import.variable(x,variable="area")
  # list(ID=ID,AREA=AREA,DOF=DOF,variable=variable)
  y <- list()
  y$log <- log(x$AREA)

  # 2D correction made in import.variables
  if(!debias)
  { y$VAR.log <- 2/x$DOF }
  else
  {
    y$VAR.log <- trigamma(x$DOF/2)

    BIAS <- digamma(x$DOF/2) - log(x$DOF/2)
    BIAS <- nant(BIAS,0)
    y$log <- y$log - BIAS
  }

  return(y)
}


exp.log.area <- function(est,VAR.est=0,VAR=0,VAR.VAR=0,debias=TRUE,level=0.95,units=TRUE,...)
{
  R <- exp.log(est=est,VAR.est=VAR.est,VAR=VAR,VAR.VAR=VAR.VAR,...)
  est <- R$est
  VAR <- R$VAR

  DOF <- 2*est^2/VAR

  if(debias)
  {
    BIAS <- digamma(DOF/2) - log(DOF/2)
    est <- est + BIAS
  }

  CI <- chisq.ci(est,VAR=VAR,level=level)

  # apply units and name
  CI <- summary.UD.format(CI,DOF/2,units=units)

  return(CI)
}


# speed
log.speed <- function(x,debias=TRUE,...)
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

exp.log.speed <- function(est,VAR.est=0,VAR=0,VAR.VAR=0,debias=TRUE,level=0.95,units=TRUE,...)
{
  # convert from log-chi to log-chi^2
  R <- exp.log(est=est,VAR.est=VAR.est,VAR=VAR,VAR.VAR=VAR.VAR,...)
  est <- R$est
  VAR <- R$VAR

  DOF <- chi.dof(est,est^2+VAR)

  if(debias)
  {
    BIAS <- digamma(DOF/2) - log(DOF/2)
    BIAS <- BIAS/2
    est <- est + BIAS
  }

  CI <- chisq.ci(est,VAR=VAR,level=level)

  # apply units and name
  CI <- summary.UD.format(CI,DOF/2,units=units)

  return(CI)
}

# models
