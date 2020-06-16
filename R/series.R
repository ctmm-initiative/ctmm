series <- function(x,coef)
{
  n <- length(coef)
  x <- c(1,x^(1:(n-1)))
  x <- coef*x
  x <- sum(x)
  return(x)
}


# log(1+x)/x --- to machine precision (including small x)
log1pxdx <- Vectorize( function(x)
{
  if(x>0.05)
  { x <- log(1+x)/x }
  else
  {
    # [1/11] Pade approximate @ x=0 to avoid numerical underflow in above representation
    coef <- c(1, 1/2, -(1/12), 1/24, -(19/720), 3/160, -(863/60480), 275/24192, -(33953/3628800), 8183/1036800, -(3250433/479001600), 4671/788480)
    x <- 1/series(x,coef)
  }
  return(x)
} )


# log(beta(a,b)) + a*log(b+a*x) --- to machine precision (including large b)
lbetaplog <- Vectorize( function(a,b,x)
{
  if(b<1/sqrt(.Machine$double.eps))
  { x <- lbeta(a,b) + a*log(b+a*x) }
  else
  {
    # Laurent series WRT 1/b @ 1/b=0 to avoid underflow
    coef <- c(lgamma(a),a/2 - a^2/2 + a^2*x,a/12 - a^2/4 + a^3/6 - (a^3*x^2)/2)
    x <- series(1/b,coef)
  }

  return(x)
})
