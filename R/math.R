# inverse psigamma function
ipsigamma <- function(x,deriv=0,precision=1/2)
{
  if(deriv==0)
  { y <- exp(x) }
  else # y == (-1)^(deriv+1) / x^deriv
  { y <- (-1)^(deriv+1)/x^(1/deriv) }

  tol <- .Machine$double.eps^precision
  ERROR <- Inf
  while(ERROR>tol)
  {
    # x == psigamma(y)
    # x == x0 + psigamma'(y0)*(y-y0) + ...
    # dx == psigamma'(y0)*dy
    # dy == dx/psigamma'(y0)
    dx <- x - psigamma(y,deriv=deriv)
    dy <- dx/psigamma(y,deriv=deriv+1)
    dy <- nant(dy,1)
    y <- y + dy
    ERROR <- max(abs(dy/y))
  }

  return(y)
}

# inverse trigamma
itrigamma <- function(x,precision=1/2)
{ ipsigamma(x,deriv=1,precision=precision) }
