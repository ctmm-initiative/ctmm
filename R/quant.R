quant <- function(x,p=0.5,low=-Inf,high=Inf)
{
  n <- length(x)
  x <- sort(x)
  q <- (1:n)/(n+1)
  # monotonically increasing ecdf
  fn <- stats::splinefun(x,q,method="hyman")
  # encapsulate p (indices)
  I <- p * (n+1)
  I <- c(floor(I),ceiling(I))
  # enclosing knots
  if(I[1]>0) { I[1] <- x[I[1]] } else { I[1] <- low }
  if(I[2]<=n) { I[2] <- x[I[2]] } else {I[2] <- high }
  # arbitrary point to extract cubic coefficients
  if(I[1]==-Inf) { X <- x[1] - 1 } # left of first knot
  else if(I[2]==Inf) { X <- x[n] + 1 } # right of last knot
  else { X <- mean(I) } # mid-point
  # extract coefficients at closest spline
  coef <- c( fn(X)-p , fn(X,deriv=1) , fn(X,deriv=2)/2 , fn(X,deriv=3)/6 )
  # solve cubic equation (fractional index)
  dx <- polyroot(coef)
  X <- X+dx # possible solution
  X <- X[abs(Im(X))<=.Machine$double.eps*n]
  X <- Re(X)
  X <- X[X>=I[1] & X<=I[2]]
  return(X)
}
