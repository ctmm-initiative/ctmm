# solve for sample size of VAR[log(det(Wishart))] of dim c(d,d)
n.varlogdet <- function(varlogdet,d,TOL=.Machine$double.eps^(3/4))
{
  # zeroth-order term of expansion
  arg <- (1:(d-1))/2
  V0 <- sum(trigamma(arg)) + pi^2/6
  if(varlogdet>V0) # pertubative solution from n = d-1 + epsilon
  { n <- d-1 + 2/sqrt(varlogdet-V0) }
  else if(varlogdet<2*d/(d-1)) # asymptotic solution from 1/n = epsilon
  { n <- 2*d/varlogdet }
  else
  { n <- d } # blind guess

  dn <- Inf
  while(abs(dn)>TOL)
  {
    arg <- (n+1-(1:d))/2
    # variance of the logarithm of the determinant of a Wishart matrix
    V0 <- sum(trigamma(arg))
    # derivative of that
    V1 <- sum(psigamma(arg,2))/2
    # Newton-Raphson iterate
    dn <- (varlogdet-V0)/V1
    # Newton-Raphson iteration under log transform: n = d-1 + exp(z)
    n <- (d-1) + (n-d+1)*exp(dn/(n-d+1))
  }

  return(n)
}
