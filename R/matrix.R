##### det shouldn't fail because R dropped indices
det.numeric <- function(x,...) { x }
determinant.numeric <- function(x,logarithm=TRUE,...)
{
  SIGN <- sign(x)
  if(logarithm)
  { x <- log(abs(x)) }

  RESULT <- list(modulus=x,sign=SIGN)
  attr(RESULT$modulus,"logarithm") <- logarithm

  class(RESULT) <- "det"

  return(det)
}

# dyadic product default
outer <- function(X,Y=X,FUN="*",...) { base::outer(X,Y,FUN=FUN,...) }

# adjoint of matrix
Adj <- function(M) { t(Conj(M)) }

# Hermitian part of matrix
He <- function(M) { (M + Adj(M))/2 }


# map function for PSD matrices
PDfunc <-function(M,func=function(m){1/m},force=FALSE,pseudo=FALSE)
{
  DIM <- dim(M)
  if(is.null(DIM)) { DIM <- 1 }
  TOL <- DIM[1]*.Machine$double.eps

  M <- eigen(M)
  V <- M$vectors
  M <- M$values

  if(any(M<=0) && !force && !pseudo) { stop("Matrix not positive definite.") }

  FORCE <- (M < TOL)
  PSEUDO <- (abs(M) < TOL)

  if(any(FORCE) && force) { M[FORCE] <- TOL }
  M <- func(M)
  if(any(PSEUDO) && pseudo) { M[PSEUDO] <- 0 }

  # add up from smallest contribution to largest contribution
  INDEX <- sort(M,method="quick",index.return=TRUE)$ix
  M <- lapply(INDEX,function(i) M[i]*(V[,i] %o% Conj(V[,i])) )
  M <- Reduce('+',M)

  return(M)
}


# Positive definite solver
PDsolve <- function(M,force=FALSE,pseudo=FALSE)
{
  DIM <- dim(M)
  if(is.null(DIM)) { DIM <- 1 }

  if(DIM[1]==1) { return(1/M) }
  if(DIM[1]==2)
  {
    DET <- M[1,1]*M[2,2]-M[1,2]*M[2,1]
    SWP <- M[1,1] ; M[1,1] <- M[2,2] ; M[2,2] <- SWP
    M[1,2] <- -M[1,2]
    M[2,1] <- -M[2,1]
    return( M/DET )
  }

  TOL <- DIM[1]*.Machine$double.eps

  # symmetrize
  M <- He(M)

  # rescale
  W <- abs(diag(M))
  W <- sqrt(W)
  W <- W %o% W

  # now a correlation matrix that is easier to invert
  M <- M/W

  # try ordinary inverse
  M.try <- try(qr.solve(M,tol=TOL))
  # fall back on decomposition
  if(class(M.try) == "matrix")
  { M <- M.try }
  else
  { M <- PDfunc(M,func=function(m){1/m},force=force,pseudo=pseudo) }

  # back to covariance matrix
  M <- M/W

  # symmetrize
  M <- He(M)

  return(M)
}


# sqrtm fails on 1x1 matrix
# it also gives annoying notes if I don't cast it right
# also, my matrices are PSD
sqrtm <- function(M,force=FALSE,pseudo=FALSE)
{
  DIM <- dim(M)
  if(is.null(DIM)) { DIM <- 1 }
  TOL <- DIM[1]*.Machine$double.eps

  if(all(DIM==1))
  {
    if(M>=0)
    { M <- sqrt(M) }
    else
    {
      if(force || pseudo) { M <- 0 } # round off error
      else { stop("Matrix is not positive definite.") } # complex sqrt
    }
  }
  else
  {
    if(all(diag(M)>=-TOL))
    {
      R <- Matrix::Matrix(M,sparse=FALSE,doDiag=FALSE) # complains still.... ???
      R <- expm::sqrtm(M) # reduces back to class "matrix" ?
    }
    else
    { R <- diag(-1,nrow=DIM[1]) }

    if(all(Re(diag(R))>=-TOL && abs(Im(diag(R)))<=TOL))
    {
      M <- Re(R)
      TEST <- (diag(M)<=0)
      if(any(TEST)) { diag(M)[TEST] <- 0 }
    }
    else
    {
      if(force || pseudo)
      { M <- PDfunc(M,func=function(m){sqrt(m)},force=force,pseudo=pseudo) }
      else
      { stop("Matrix is not positive definite.") }
    }
  }

  return(M)
}


# condition number
conditionNumber <- function(M)
{
  M <- eigen(M)$values
  M <- range(M)
  return(M[2]/M[1])
}


# Positive definite part of matrix
PDclamp <- function(M)
{
  # symmetrize
  M <- He(M)

  M <- PDfunc(M,function(m){clamp(m,0,Inf)},force=TRUE)

  # symmetrize
  M <- He(M)

  return(M)
}
