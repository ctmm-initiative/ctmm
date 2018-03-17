svf.test <- function(variogram,CTMM)
{
  svf <- svf.func(CTMM,moment=TRUE,error=as.logical(CTMM$error))
  SVF <- Vectorize(svf$svf)(variogram$lag[-1])
  DOF <- Vectorize(svf$DOF)(variogram$lag[-1])

  Fstat <- SVF/variogram$SVF[-1]
  Fstat <- stats::pf(Fstat,DOF,variogram$DOF[-1],lower.tail=FALSE)

  RETURN <- min(Fstat)

  return(RETURN)
}

# signature (CTMM,variogram)

# signature (variogram,variogram)
