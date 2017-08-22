# simulate an estimated model
emulate <- function(CTMM,data=NULL)
{
  if(is.null(data)) { return(emulate.ctmm(CTMM)) }
  else { return(emulate.telemetry(data,CTMM)) }
}

# simulate an estimated model from the parameter estimate covariances
emulate.ctmm <- function(CTMM)
{
  # needs to be updated for periodic stuff
  CTMM$mu <- MASS::mvrnorm(mu=CTMM$mu,Sigma=CTMM$COV.mu)

  COV <- CTMM$COV
  NAMES <- dimnames(COV)[[1]]
  par <- get.parameters(CTMM,NAMES)

  # positive variables (potential)
  PAR <- c("area","tau position","tau velocity","error")
  # included variables that are positive
  PAR <- PAR[PAR %in% NAMES]
  # log transform positive parameters
  for(P in PAR)
  {
    COV[P,] <- COV[P,]/par[P]
    COV[,P] <- COV[,P]/par[P]

    par[P] <- log(par[P])
  }

  par <- MASS::mvrnorm(mu=par,Sigma=COV)

  # transform log back to positive parameters
  par[PAR] <- exp(par[PAR])

  CTMM <- set.parameters(CTMM,par)

  return(CTMM)
}

# simulate an estimated model by simulating data and then fitting a model to that data
emulate.telemetry <- function(data,CTMM)
{
  data <- simulate(CTMM,t=data$t) #!!! NEED TO HAVE AN HDOP/HERE option here???
}
