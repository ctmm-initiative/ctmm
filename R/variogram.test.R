svf.test <- function(variogram,CTMM)
{
  if("MSE" %in% names(variogram)) { MSE <- variogram$MSE }
  else { MSE <- CTMM$error * variogram$MSDOP }

  svf <- svf.func(CTMM,moment=TRUE)
  SVF <- Vectorize(svf$svf)(variogram$lag[-1],error=MSE[-1])
  DOF <- Vectorize(svf$DOF)(variogram$lag[-1],error=MSE[-1])

  Fstat <- SVF/variogram$SVF[-1]
  Fstat <- pmax(Fstat,1/Fstat)
  Fstat <- stats::pf(Fstat,DOF,variogram$DOF[-1],lower.tail=FALSE)

  RETURN <- min(Fstat)

  return(RETURN)
}
