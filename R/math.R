EulerGamma <- 0.57721566490153286061
Zeta3 <- 1.2020569031595942854

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


# shifted Legendre polynomials
legendre <- function(n,x)
{
  if(n==0)
  { 1 }
  else if(n==1)
  { 2*x - 1 }
  else if(n==2)
  { 6*x^2 - 6*x + 1 }
  else if(n==3)
  { 20*x^3 - 30*x^2 + 12*x - 1 }
  else if(n==4)
  { 70*x^4 - 140*x^3 + 90*x^2 - 20*x + 1 }
  else if(n==5)
  { 252*x^5 - 630*x^4 + 560*x^3 - 210*x^2 + 30*x - 1 }
}

# sinc functions
sinc <- function(x,SIN=sin(x))
{
  return(ifelse(x==0, 1, SIN/x))
}

sinch <- function(x,SINH=sinh(x))
{
  return(ifelse(x==0, 1, SINH/x))
}


# multivariate polygamma function
mpsigamma <- function(x,deriv=0,dim=1)
{
  PSI <- 1 - 1:dim
  PSI <- x + PSI/2
  if(deriv>=0) { PSI <- sapply(PSI,function(p) psigamma(p,deriv=deriv)) }
  else if(deriv==-1) { PSI <- sapply(PSI,function(p) lgamma(p)) }
  else { stop("Derivative ",deriv+1," of log(Gamma(x)) not supported.") }
  PSI <- sum(PSI)
  return(PSI)
}


logit <- function(p)
{
  p <- clamp(p,0,1)
  log(p/(1-p))
}

ilogit <- function(z) { 1/(1+exp(-z)) }

binom <- function(x,y) { 1/(x+1)/beta(y+1,x-y+1) }

lbinom <- function(x,y) { -log(x+1) - lbeta(y+1,x-y+1) }
