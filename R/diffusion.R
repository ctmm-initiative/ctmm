# calculate the max_lag diffusion rate for a movement model
diffusion <- function(CTMM,level=0.95,finish=TRUE)
{
  CTMM <- get.taus(CTMM) # pre-calculate stuff

  range <- CTMM$range
  tau <- CTMM$tau
  omega <- CTMM$omega
  f <- CTMM$f.nu[1]
  nu <- CTMM$f.nu[2]
  Omega2 <- CTMM$Omega2
  circle <- abs(CTMM$circle) # only magnitude matters for this

  sigma <- var.covm(CTMM$sigma,ave=TRUE) # variance
  COV <- CTMM$COV
  if(!is.null(COV)) { COV <- axes2var(CTMM,MEAN=TRUE) }

  NAMES <- CTMM$tau.names
  if(!length(tau)) # IID
  { return( c(0,Inf,Inf) ) }
  else if(!range) # Brownian motion # IOU # max diffusion rate at infinite lag
  {
    D <- 1
    D.grad <- NULL
    NAMES <- NULL
  }
  else if(length(tau)==1) # OU
  {
    if(circle*tau <= 1) # max diffusion rate at zero lag
    {
      D <- 1/tau
      D.grad <- -1/tau^2
    }
    else # circulation enhanced max diffusion rate
    {
      D1 <- circle*tau
      D <- atan((D1^2-1)/(2*D1))/circle
      D.grad <- 2/(1+D1^2)/circle * c(circle,tau) - c(0,D/circle)
      NAMES <- c(NAMES,"circle")
    }
  }
  else if(length(tau)==2 && !omega && tau[1]!=tau[2] && !circle) # OUF
  {
    D1 <- tau[2]/tau[1]
    D <- D1^(D1/(1-D1)) / tau[1]
    D.grad <- D * (1/(1-D1) + log(D1)/(1-D1)^2 ) * c(-D1/tau[1],1/tau[1]) - c(D/tau[1],0)
  }
  else if(length(tau)==2 && omega && tau[1]==tau[2] && !circle) # OUO
  {
    D1 <- tau[1]*omega
    D <- sqrt(1+D1)*exp(-atan(D1)/D1) / tau[1]
    D.grad <- D * ( (D1^1-1)/(D1^1+1)/D1 + atan(D1)/D1^2 ) * c(omega,tau[1]) - c(D/tau[1],0)
  }
  else if(length(tau)==2 && !omega && tau[1]==tau[2] && !circle) # OUf
  {
    D <- exp(-1)/tau[1]
    D.grad <- -D/tau[1]
  }
  else if(circle) # OUF, OUO, OUf with circulation can't be solved analytically
  {
    NAMES <- c(NAMES,"circle")

    if(!omega && tau[1]!=tau[2]) # OUF
    {
      A0.fn <- function(t) { diff(exp(-t/tau)/tau)/diff(tau) }
      D0.fn <- function(t) { diff(exp(-t/tau))/diff(tau) }
      J <- diag(2)
      LAG0 <- log(tau[1]/tau[2])/(1/tau[2]-1/tau[1])
    }
    else if(!omega && tau[1]==tau[2]) # OUf
    {
      A0.fn <- function(t) { (1+f*t) * exp(-f*t) }
      D0.fn <- function(t) { Omega2/nu * sin(nu*t) * exp(-f*t) }
      J <- CTMM$J.f.tau
      LAG0 <- tau[1]
    }
    else if(omega) # OUO
    {
      A0.fn <- function(t) { (cos(nu*t)+f/nu*sin(nu*t)) * exp(-f*t) }
      D0.fn <- function(t) { Omega2/nu * sin(nu*t) * exp(-f*t) }
      J <- CTMM$J.nu.tau
      LAG0 <- atan(nu*t)/nu
    }

    # total diffusion rate function of lag
    D.fn <- function(t) { D0.fn(t)*cos(circle*t) + A0.fn(t)*circle*sin(circle*t) }

    upper <- max(LAG0,pi/circle)

    MAX <- stats::optimize(D.fn,lower=0,upper=upper,maximum=TRUE,tol=.Machine$double.eps^0.5)
    LAG <- MAX$objective
    D <- MAX$maximum

    # zero circulation result
    R0 <- CTMM
    R0$circle <- FALSE
    R0 <- diffusion(R0,finish=FALSE)

    D.grad <- c(R0$D.grad,(D-R0$D)/circle) # crude calculation of gradient if circle>>0
    # very annoying to calculate this better (optimize inside gradient)
  }

  # return information for CIs but not completed CIs
  if(!finish) { return(list(D=D,D.grad=D.grad,NAMES=NAMES)) }

  NAMES <- c(NAMES,"variance")
  D.grad <- c( sigma * D.grad, D )
  D <- sigma * D

  names(D.grad) <- NAMES
  if(!is.null(COV)) { VAR <- c(D.grad %*% COV[NAMES,NAMES] %*% D.grad) }
  else { VAR <- Inf }

  DOF <- 2*D^2/VAR / length(CTMM$axes)
  return(list(D=D,VAR=VAR,DOF=DOF))

  CI <- chisq.ci(D,VAR,level=level)
  return(CI)
}
