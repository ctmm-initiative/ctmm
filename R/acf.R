############
residuals.ctmm <- function(object,data,...)
{
  object <- ctmm.prepare(data,object)
  axes <- object$axes

  # detrend mean
  drift <- get(object$mean)
  drift <- drift(data$t,object) %*% object$mu
  data[,axes] <- get.telemetry(data,axes=axes) - drift
  rm(drift)

  # calculate residuals
  data[,axes] <- smoother(data,object,residual=TRUE)
  data <- data[,c("t",axes)]

  attr(data,"info")$UERE <- NULL
  attr(data,"info")$residual <- TRUE

  # effective DOF of variance estimate
  # n <- 2*object$sigma@par["area"]^2/object$COV["area","area"]
  # COV <- cbind(2/n)
  # dimnames(COV) <- list("area","area")
  # data@info$CTMM <- ctmm(mu=rbind(c(0,0)),sigma=1,isotropic=T,COV=COV)

  return(data) # check class of this object: data.frame or telemetry
}
residuals.telemetry <- function(object,CTMM,...) { residuals.ctmm(CTMM,object,...) }

# acf function
correlogram <- function(data,dt=NULL,fast=TRUE,res=1,axes=c("x","y"))
{
  if(fast) { ACF <- variogram.fast(data,dt=dt,res=res,CI="IID",axes=axes,ACF=TRUE) }
  else { ACF <- variogram.slow(data,dt=dt,CI="IID",axes=axes,ACF=TRUE) }

  # normalize ACF to 1 at lag 0
  ACF$SVF <- ACF$SVF / ACF$SVF[1]

  ACF <- new.variogram(ACF,info=attr(data,"info"))
  attr(ACF,"info")$ACF <- TRUE
  # attr(ACF,"info")$CTMM <- NULL
  return(ACF)
}
