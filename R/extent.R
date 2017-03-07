# 2D range of ctmm object
# I can't get the raster generic to work?? It won't accept my classes?
# extent <- function(x,...) UseMethod("extent")

# what is the range of a list of objects of any type
extent.list <- function(x,...)
{
  # list of element ranges
  RANGE <- lapply(x,function(d){extent(d,...)})
  # concatenate ranges
  RANGE <- Reduce(rbind,RANGE)
  # final range x & y
  RANGE <- extent.telemetry(RANGE,level=1)
  # level has already been applied
  return(RANGE)
}
# raster doesn't really need this?
setMethod('extent', signature(x='list'), extent.list)


# range of telemetry data
extent.telemetry <- function(x,level=1,...)
{
  alpha <- (1-level)/2
  X <- stats::quantile(x$x,probs=c(alpha,1-alpha))
  Y <- stats::quantile(x$y,probs=c(alpha,1-alpha))
  RANGE <- data.frame(x=X,y=Y)
  row.names(RANGE) <- c("min","max")
  return(RANGE)
}
setMethod('extent', signature(x='telemetry'), extent.telemetry)


# range of Gaussian contours
extent.ctmm <- function(x,level=0.95,level.UD=0.95,...)
{
  if(is.null(x$mu)) { stop("This model has no mean location. Try ctmm.guess.") }
  
  alpha <- 1 - level
  alpha.UD <- 1 - level.UD
  # standard Gaussian quantile stuff
  z <- sqrt(-2*log(alpha.UD))
  
  # proportionality constants for outer CIs
  sigma <- x$sigma
  
  # capture outer contour if present
  const <- 1
  if(is.na(level)) { alpha <- 1 } # will use ML contour
  if(!is.null(x$COV))
  {
    K <- length(x$tau)
    const <- confint.ctmm(x,alpha)[1,3]/sqrt(det(sigma))
  }
  
  buff <- z*sqrt(const*diag(sigma))
  
  X <- x$mu[1] + c(-1,1)*buff[1]
  Y <- x$mu[2] + c(-1,1)*buff[2]
  
  RANGE <- data.frame(x=X,y=Y)
  row.names(RANGE) <- c("min","max")
  return(RANGE)
}
setMethod('extent', signature(x='ctmm'), extent.ctmm)


# range of UD contours
extent.UD <- function(x,level=0.95,level.UD=0.95,...)
{
  if(level.UD==1 || (level==1 && !is.null(x$DOF.area)))
  {
    RANGE <- data.frame(x=c(-1,1),y=c(-1,1))*Inf
    row.names(RANGE) <- c("min","max")
    return(RANGE)
  }
  
  # capture ML contour
  if(is.null(x$DOF.area))
  { P <- level.UD }
  else # capture outer contour
  { P <- CI.UD(x,level.UD,level,P=TRUE)[3] }
  
  # do we extend this far?
  MAT <- (x$CDF <= P)
  X <- apply(MAT,1,any)
  Y <- apply(MAT,2,any)
  rm(MAT)
  
  # now indices
  X <- which(X)
  Y <- which(Y)
  
  # range indices
  X <- c(X[1],last(X))
  Y <- c(Y[1],last(Y))
  
  # location ranges
  X <- x$r$x[X]
  Y <- x$r$y[Y]
  
  RANGE <- data.frame(x=X,y=Y)
  row.names(RANGE) <- c("min","max")
  return(RANGE)
}
setMethod('extent', signature(x='UD'), extent.UD)
