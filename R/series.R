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


# (sqrt(x+1)-1)/x
sqrtxp1 <- Vectorize( function(x)
{
  if(x>0.5) { return((sqrt(x+1)-1)/x) }

  x <- sqrt(x)
  cn <- c(1/2, 0, 19/8, 0, 153/32, 0, 85/16, 0, 455/128, 0, 3003/2048, 0, 3003/8192, 0, 429/8192, 0, 495/131072, 0, 55/524288, 0, 1/2097152)
  cd <- c(1, 0, 5, 0, 171/16, 0, 51/4, 0, 595/64, 0, 273/64, 0, 5005/4096, 0, 429/2048, 0, 1287/65536, 0, 55/65536, 0, 11/1048576)
  series(x,cn)/series(x,cd)
})


# log K[n/2+r/2](t*sqrt(1+s/t))/K[r/2](t)
lKK <- function(r,n,t,s)
{
  # log( besselK(sqrt(1+s/t)*t,n/2+r/2)/besselK(t,r/2) )
  # 0/0 limits to 1
  -s*sqrtxp1(s/t) + nant( BesselK(sqrt(1+s/t)*t,n/2+r/2,expon.scaled=TRUE,log=TRUE) - BesselK(t,r/2,expon.scaled=TRUE,log=TRUE), 0)
}


BesselK <- function(x,nu,expon.scaled=FALSE,log=FALSE)
{
  nu <- abs(nu)

  # try base besselK
  y <- log(besselK(x,nu,expon.scaled=expon.scaled))
  MAX <- log(.Machine$double.xmax/2)
  # overflow cases
  SUB <- y>=MAX & x>0 & nu>1 & nu<Inf
  if(any(SUB)) # nu too big for base Bessel
  { y[SUB] <- Bessel::besselK.nuAsym(x[SUB],nu[SUB],k.max=5,expon.scaled=expon.scaled,log=TRUE) }

  if(!log) { y <- exp(y) }

  return(y)
}
