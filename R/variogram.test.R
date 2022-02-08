svf.test <- function(x,y,test="F",level=0.95)
{
  if(class(x)=="ctmm")
  {
    z <- x
    x <- y # variogram
    y <- z # ctmm
  }

  if(class(x)=="variogram" && class(y)=="ctmm")
  {
    # # expected location error contribution
    # if("MSE" %in% names(x)) { y$MSE <- x$MSE }
    # else { y$MSE <- y$error^2 * x$MSDOP }

    svf <- svf.func(y,moment=TRUE)
    y$SVF <- Vectorize(svf$svf)(x$lag)
    y$DOF <- Vectorize(svf$DOF)(x$lag)
  }

  # if(class(x)=="variogram" && class(y)=="variogram")
  # {
  #   # TODO consider only common lags
  # }

  Fstat <- x$SVF/y$SVF
  Fstat <- pmin( stats::pf(Fstat,x$DOF,y$DOF,lower.tail=TRUE) , stats::pf(Fstat,x$DOF,y$DOF,lower.tail=FALSE) )
  R <- which.min(Fstat)
  Fstat <- Fstat[R]
  if(x$SVF[R]>y$SVF[R])
  {
    S1 <- x$SVF[R]
    V1 <- 2*S1^2/x$DOF[R]
    S2 <- 1/y$SVF[R]
    V2 <- 2*S2^2/y$DOF[R]
    R <- F.CI(S1,V1,S2,V2,level=level)
  }
  else
  {
    S1 <- y$SVF[R]
    V1 <- 2*S1^2/y$DOF[R]
    S2 <- 1/x$SVF[R]
    V2 <- 2*S2^2/x$DOF[R]
    R <- F.CI(S1,V1,S2,V2,level=level)
  }
  names(R) <- NAMES.CI

  # else if(test=="overlap")
  # {
  #   x$VAR <- 2*x$SVF^2/x$DOF
  #   y$VAR <- 2*y$SVF^2/y$DOF
  #   VAR.F <- 1/y$SVF^2*x$VAR + (-x$SVF/y$SVF^2)^2*y$VAR
  #
  #   # Gaussian overlap argument
  #   VAR <- (1/4-1/Fstat^2/4)^2*VAR.F
  #   R <- (1/2+Fstat/4+1/Fstat/4)
  #   # distance
  #   VAR <- (1/R/4)^2*VAR
  #   R <- log(R)/4
  #
  #   # maximum distance
  #   MAX <- which.max(R)
  #   R <- R[MAX]
  #   VAR <- VAR[MAX]
  #
  #   R <- chisq.ci(R,VAR,level=level)
  #   R <- rev(exp(-R)) # overlap (minimum)
  #   names(R) <- NAMES.CI
  # }

  return(R)
}
