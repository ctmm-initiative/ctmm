# predicted speeds over time
# what is the easiest way to include parameter uncertainty?
# add independent variance from RMS speed?
speeds <- function(data,CTMM=NULL,level=0.95,HDR=FALSE,fast=TRUE,...)
{
  n <- nrow(data)

  DOF <- summary(CTMM)$DOF['speed']
  if(!DOF)
  {
    v2 <- rep(Inf,n) # upgrade to outlie estimates?
    DOF <- numeric(n)
  }
  else
  {
    if(!is.null(CTMM)) { data <- predict(data,CTMM=CTMM,t=data$t,...) }

    axes <- c("vx","vy")
    v <- get.telemetry(data,axes=axes) # (n,2)
    VAR <- get.error(data,list(axes=axes,error=TRUE),DIM=2) # (n,2,2)

    # delta method - variance of square speed
    VAR <- 4 * vapply(1:n,function(i){v[i,] %*% VAR[i,,] %*% v[i,]},numeric(1)) # (n)

    # speeds
    v2 <- rowSums(v^2)

    # chi-squared DOF -- VAR sum rule
    DOF <- 1/( 1/DOF + VAR/(2*v2^2) )
  }


  # if no level, return point estimate and DOF
  if(is.null(level))
  {
    # output v and chi DOF
    v <- sqrt(v2)
    v <- cbind(v,DOF)
    colnames(v) <- c("speed","DOF")
  }
  else
  {
    v2 <- vapply(1:n,function(i){ chisq.ci(v2[i],DOF=DOF[i],level=level,HDR=HDR) },numeric(3)) # (3,n)
    v <- sqrt(t(v2))
  }

  return(v)
}
