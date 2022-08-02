# R refuses to call log.list for list objects outside of ctmm
Log <- function(x,debias=TRUE,...)
{
  x <- listify(x)
  NAMES <- names(x)
  x <- log.list(x,debias=debias,...)
  names(x$log) <- names(x$VAR.log) <- NAMES
  return(x)
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

exp.area <- function()
{
  # DOF <- 2*itrigamma(VAR)
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
    y$VAR.log <- trigamma(x$DOF/2)/4

    BIAS <- ( digamma(x$DOF/2) - log(x$DOF/2) )/2
    BIAS <- nant(BIAS,0)
    y$log <- y$log - BIAS
  }

  return(y)
}

# models
