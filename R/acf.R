############
residuals.ctmm <- function(object,data,...)
{
  object <- ctmm.prepare(data,object)
  error <- get.error(data,object)
  t <- data$t
  dt <- c(Inf,diff(t))
  
  # effective DOF of variance estimate
  n <- 2*object$sigma@par["area"]^2/object$COV["area","area"]
  
  z <- cbind(data$x,data$y)
  # subtract off mean function
  drift <- get(object$mean)
  mu <- drift(data$t,object) %*% object$mu
  z <- z - mu
  
  sigma <- object$sigma
  area <- sigma@par["area"]
  theta <- sigma@par["angle"]
  ecc <- sigma@par["eccentricity"]
  
  circle <- object$circle
  if(circle) { circle <- 2*pi/circle }
  
  
  # rotate to eigen-basis
  z <- z %*% t(rotate(-theta))
  
  if(circle)
  {
    # proportional standardization from ellipse to circle
    if(ecc)
    {
      z[,1] <- z[,1] * exp(-ecc/4)
      z[,2] <- z[,2] * exp(+ecc/4)
    }
    z <- cbind(z[,1] + 1i*z[,2])
    
    # corotating frame
    R <- exp(-1i*circle*(t-t[1]))
    z <- R * z
    
    # preserves error volume, but distorts error ellipse
    object$sigma <- area
    z <- kalman(z,u=NULL,dt=dt,object,error=error,residual=TRUE)
    
    # de-complexify result
    z <- cbind(Re(z),Im(z))
  }
  else
  {
    # send each axis through with correct variance and error
    object$sigma <- area*exp(+ecc/2)
    z[,1] <- kalman(z[,1,drop=FALSE],u=NULL,dt=dt,object,error=error,residual=TRUE)
    
    object$sigma <- area*exp(-ecc/2)
    z[,2] <- kalman(z[,2,drop=FALSE],u=NULL,dt=dt,object,error=error,residual=TRUE)
  }
  
  attr(data,"info")$UERE <- NULL
  attr(data,"info")$residual <- TRUE
  COV <- cbind(2/n)
  dimnames(COV) <- list("area","area")
  # data@info$CTMM <- ctmm(mu=rbind(c(0,0)),sigma=1,isotropic=T,COV=COV)
  data <- data[,c("t","x","y")]
  data$x <- z[,1]
  data$y <- z[,2]
  
  return(data) # check class of this object: data.frame or telemetry
}
residuals.telemetry <- function(object,CTMM,...) { residuals.ctmm(CTMM,object,...) }

# acf function
correlogram <- function(data,dt=NULL,axes=c("x","y"))
{
  ACF <- variogram.fast(data,dt=NULL,CI="IID",axes=c("x","y"),ACF=TRUE)
  ACF <- new.variogram(ACF,info=attr(data,"info"))
  attr(ACF,"info")$ACF <- TRUE
  # attr(ACF,"info")$CTMM <- NULL
  return(ACF)
}