#################################
# Nelder-Mead parallel optimizer for multiple CPU cores
mcoptim.NM <- function(par,fn,control=list())
{
  if(is.null(control$mc.cores)) { mc.cores <- parallel::detectCores() }
  if(is.null(control$alpha)) { alpha <- 1 }
  if(is.null(control$alpha)) { beta <- 2 }
  if(is.null(control$alpha)) { gamma <- 1/2 }
  if(is.null(control$alpha)) { sigma <- 1/2 }
  
  
  
  DIM <- length(par)
  lower <- rep(lower,DIM)
  upper <- rep(upper,DIM)
  
  # default COV initialization
  COV <- control$COV
  if(is.null(COV))
  {
    COV <- pmax(par-lower,upper-par)
    COV[COV==Inf] <- par
    COV[COV==0] <- 1
    COV <- diag(COV^2)
  }
  
  # while !two.in.a.row.tol & maxit
  
  
  # elipsoid transformation
  #STD <- sqrtm(COV)
  #R <- STD %*% R + par
  
  # mclapply
  
  # Bayesian update initially
  # quadratic regression eventually
  
  
  detach(control)
  return()
}