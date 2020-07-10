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
    # expected location error contribution
    if("MSE" %in% names(x)) { y$MSE <- x$MSE }
    else { y$MSE <- y$error^2 * x$MSDOP }

    svf <- svf.func(y,moment=TRUE)
    y$SVF <- Vectorize(svf$svf)(x$lag,error=y$MSE)
    y$DOF <- Vectorize(svf$DOF)(x$lag,error=y$MSE)
  }

  if(class(x)=="variogram" && class(y)=="variogram")
  {
    # TODO consider only common lags
  }

  Fstat <- x$SVF/y$SVF
  if(test=="F")
  {
    TAIL <- Fstat < 1
    Fstat <- stats::pf(Fstat,x$DOF,y$DOF,lower.tail=TAIL)
    R <- min(Fstat[-1])
  }
  else if(test=="overlap")
  {
    x$VAR <- 2*x$SVF^2/x$DOF
    y$VAR <- 2*y$SVF^2/y$DOF
    VAR.F <- 1/y$SVF^2*x$VAR + (-x$SVF/y$SVF^2)^2*y$VAR

    # Gaussian overlap argument
    VAR <- (1/4-1/Fstat^2/4)^2*VAR.F
    R <- (1/2+Fstat/4+1/Fstat/4)
    # distance
    VAR <- (1/R/4)^2*VAR
    R <- log(R)/4

    # maximum distance
    MAX <- which.max(R)
    R <- R[MAX]
    VAR <- VAR[MAX]

    R <- chisq.ci(R,VAR,level=level)
    R <- rev(exp(-R)) # overlap (minimum)
    names(R) <- NAMES.CI
  }

  return(R)
}
