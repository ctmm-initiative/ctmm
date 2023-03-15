# NECESSARY NOT-SO-SPECIAL FUNCTIONS

# the z->1 limit is numerically unstable
Diff.OUF.fn <- function(z,deriv=0)
{
  w <- 1-z
  if(w > 0.004)
  {
    if(deriv==0)
    { z <- z^(z/w) }
    else if(deriv==1)
    { z <- Diff.OUF.fn(z) * ( w + log(z) )/w^2 }
  }
  else
  {
    if(deriv==0)
    { z <- exp(-1) * series(w,c(1,1/2,7/24,3/16,743/5760,215/2304)) }
    else if(deriv==1)
    { z <- -Diff.OUF.fn(z) * series(w,1/2:10) }
  }
  return(z)
}

# NaN at 0
Diff.OUO.fn <- function(x,deriv=0)
{
  if(deriv==0)
  { x <- sqrt(1+x^2) * exp(-ifelse(x==0,1,atan(x)/x)) }
  else if(deriv==1)
  {
    if(x>0.1)
    { x <- Diff.OUO.fn(x)/x^2 * ( x - 2*x/(1+x^2) + atan(x) ) }
    else
    {
      n <- 10
      coef <- seq(5,by=4,length.out=n)/seq(3,by=2,length.out=n) * (-1)^(1+1:n)
      x <- Diff.OUO.fn(x) * x*series(x^2,coef)
    }
  }
  return(x)
}


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
  circle <- FALSE # turning circulation off, because it reduces diffusion rate

  sigma <- var.covm(CTMM$sigma,ave=FALSE) # variance
  COV <- CTMM$COV
  if(!is.null(COV))
  {
    COV <- axes2var(CTMM,MEAN=FALSE)
    D.grad <- rep(0,nrow(COV))
    names(D.grad) <- rownames(COV)
  }
  else
  { D.grad <- NULL }
  # this is redundant with D.grad, but needed for mean.ctmm summary
  J <- rep(0,length(CTMM$features))
  names(J) <- CTMM$features

  if(!length(tau) || all(tau==0)) # IID
  {
    if(!finish) { return(list(D=Inf,grad=0,VAR=Inf,DOF=0)) }
    return( c(0,Inf,Inf) )
  }
  else if(!range) # Brownian motion # IOU # max diffusion rate at infinite lag
  { D <- 1 }
  else if(length(tau)==1 || tau[2]==0) # OU
  {
    tau <- tau[1]
    NAMES <- CTMM$tau.names[1]

    if(circle*tau <= 1) # max diffusion rate at zero lag
    {
      D <- 1/tau
      D.grad[NAMES] <- J[NAMES] <- -1/tau^2
    }
    else # circulation enhanced max diffusion rate
    {
      NAMES <- c(NAMES,"circle")
      z <- circle*tau
      z.grad <- c(circle,tau)
      # D <- atan((z^2-1)/(2*z))/circle
      D <- (2*atan(z)-pi/2)/circle
      D.grad[NAMES] <- J[NAMES] <- 2/(1+z^2)/circle * z.grad - c(0,D/circle)
    }
  }
  else if(length(tau)==2)
  {
    NAMES <- CTMM$tau.names

    if(!omega && tau[1]!=tau[2] && !circle) # OUF
    {
      z <- tau[2]/tau[1]
      z.grad <- c(-1,1)*z/tau
      D <- Diff.OUF.fn(z) / tau[1]
      D.grad[NAMES] <- J[NAMES] <- Diff.OUF.fn(z,deriv=1)/tau[1]*z.grad - c(D/tau[1],0)
    }
    else if(omega && tau[1]==tau[2] && !circle) # OUO
    {
      z <- tau[1]*omega
      z.grad <- c(omega,tau[1])
      D <- Diff.OUO.fn(z) / tau[1]
      D.grad[NAMES] <- J[NAMES] <- Diff.OUO.fn(z,deriv=1)/tau[1]*z.grad - c(D/tau[1],0)
    }
    else if(!omega && tau[1]==tau[2] && !circle) # OUf
    {
      NAMES <- NAMES[1]
      D <- exp(-1)/tau[1]
      D.grad[NAMES] <- J[NAMES] <- -D/tau[1]
    }
  }
  # else if(circle) # OUF, OUO, OUf with circulation can't be solved analytically
  # {
  #   if(!omega && tau[1]!=tau[2]) # OUF
  #   {
  #     NAMES <- paste("tau",NAMES)
  #     S0.fn <- function(t) { -diff(exp(-t/tau)*tau)/diff(tau) }
  #     D0.fn <- function(t) { diff(exp(-t/tau))/diff(tau) }
  #     D1.fn <- function(t) { -diff(exp(-t/tau)/tau)/diff(tau) }
  #     J <- diag(2)
  #     LAG0 <- log(tau[1]/tau[2])/(1/tau[2]-1/tau[1])
  #   }
  #   else if(!omega && tau[1]==tau[2]) # OUf
  #   {
  #     S0.fn <- function(t) { -(1+f*t) * exp(-f*t) }
  #     D0.fn <- function(t) { f^2*t * exp(-f*t) }
  #     D1.fn <- function(t) { f^2*(1-f*t) * exp(-f*t) }
  #     J <- CTMM$J.f.tau
  #     LAG0 <- tau[1]
  #   }
  #   else if(omega) # OUO
  #   {
  #     S0.fn <- function(t) { -(cos(nu*t)+f/nu*sin(nu*t)) * exp(-f*t) }
  #     D0.fn <- function(t) { Omega2/nu * sin(nu*t) * exp(-f*t) }
  #     D1.fn <- function(t) { Omega2 * ( cos(nu*t) - f/nu*sin(nu*t) ) * exp(-f*t) }
  #     J <- CTMM$J.nu.tau
  #     LAG0 <- atan(nu*t)/nu
  #   }
  #
  #   # negative diffusion rate function of lag
  #   nD.fn <- function(t) { -(D0.fn(t)*cos(circle*t) - S0.fn(t)*circle*sin(circle*t)) }
  #
  #   MAX <- optimizer(LAG0,nD.fn,lower=0)
  #   LAG <- MAX$par
  #   D <- -MAX$value
  #
  #   # UNFINISHED
  #   # UNFINISHED
  #   # UNFINISHED
  #
  #   # zero circulation result
  #   R0 <- CTMM
  #   R0$circle <- FALSE
  #   R0 <- diffusion(R0,finish=FALSE)
  #
  #   D.grad <- c(R0$D.grad,(D-R0$D)/circle) # crude calculation of gradient if circle>>0
  #   # very annoying to calculate this better (optimize inside gradient)
  # }

  D.grad <- sigma * D.grad
  J <- sigma * J

  D.grad["variance"] <- D
  if(CTMM$isotropic[1])
  { J["major"] <- 2*D }
  else
  { J["major"] <- J["minor"] <- D }
  D <- sigma * D

  if(!is.null(COV))
  {
    VAR <- c(D.grad %*% COV %*% D.grad)
    VAR <- nant(VAR,Inf)
  }
  else
  { VAR <- Inf }

  # this is for chi^2 CIs
  DOF <- 2*D^2/VAR # / length(CTMM$axes) # /2 is for counting

  # return information for CIs but not completed CIs
  if(!finish) { return(list(D=D,grad=D.grad,VAR=VAR,DOF=DOF,J=J)) }

  CI <- chisq.ci(D,VAR,level=level)
  return(CI)
}
