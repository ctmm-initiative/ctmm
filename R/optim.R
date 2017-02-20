#################################
# quasi-Newton parallel optimizer for multiple CPU cores
mcoptim <- function(par,fn,gr=NULL,lower=-Inf,upper=Inf,control=list(COV=NULL,mc.cores=parallel::detectCores()))
{
  # CRAN DOESN'T LIKE ATTACH
  #attach(control)
  
  DIM <- length(par)
  lower <- rep(lower,DIM)
  upper <- rep(upper,DIM)
  
  # default COV initialization
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