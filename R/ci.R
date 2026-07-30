ci.beta <- function(est,var,se=sqrt(var),level=0.95,...)
{
  if(missing(var)) { var = se^2 }
  CI <- t(Vectorize(beta.ci)(MLE=est,VAR=var,level=level,...))
  if(nrow(CI)==1) { CI <- CI[1,] }
  return(CI)
}

# beta distributed CI given mean and variance estimates
beta.ci <- function(MLE,VAR,level=0.95,alpha=1-level)
{
  MLE <- nant(MLE,0)
  VAR <- nant(VAR,Inf)
  if(VAR == 0) {
    # handle the degenerate case (the beta distribution becomes a delta)
    CI <- c(MLE, MLE, MLE)
  } else {
    n <- MLE*(1-MLE)/VAR - 1
    if(n<=0)
    { CI <- c(0,MLE,1) }
    else
    {
      a <- n * MLE
      b <- n * (1-MLE)
      CI <- c(0, MLE, 1)
      aux <- stats::qbeta(c(alpha/2, 1-alpha/2), a, b, ncp=0)
      CI[1] <- aux[1]
      CI[3] <- aux[2]
    }
  }
  names(CI) <- NAMES.CI
  return(CI)
}


# chisq fast confidence interval functions
CI.upper <- Vectorize(function(k,Alpha){stats::qchisq(Alpha/2,k,lower.tail=FALSE)/k})
CI.lower <- Vectorize(function(k,Alpha){stats::qchisq(Alpha/2,k,lower.tail=TRUE)/k})

ci.chisq <- function(est,var,se=sqrt(var),level=0.95,...)
{
  if(missing(var)) { var = se^2 }
  CI <- t(Vectorize(chisq.ci)(MLE=est,VAR=var,level=level,...))
  if(nrow(CI)==1) { CI <- CI[1,] }
  return(CI)
}

# calculate chi^2 confidence intervals from MLE and COV estimates
chisq.ci <- function(MLE,VAR=NULL,level=0.95,alpha=1-level,DOF=2*MLE^2/VAR,robust=FALSE,HDR=FALSE)
{
  #DEBUG <<- list(MLE=MLE,COV=COV,level=level,alpha=alpha,DOF=DOF,robust=robust,HDR=HDR)
  # try to do something reasonable on failure cases
  if(is.nan(DOF) || is.na(DOF)) { DOF <- 0 } # NaN comes from infinite variance divsion
  if(is.na(MLE)) { MLE <- Inf }

  if(is.na(level))
  {
    VAR <- 2*MLE^2/DOF
    CI <- c(DOF,MLE,VAR)
    names(CI) <- c("DOF","est","VAR")
    return(CI)
  }

  if(DOF==Inf)
  { CI <- c(1,1,1)*MLE }
  else if(DOF==0)
  { CI <- c(0,MLE,Inf) }
  else if(MLE==0)
  { CI <- c(0,0,0) }
  else if(MLE==Inf)
  { CI <- c(Inf,Inf,Inf) }
  else if(!is.null(VAR) && VAR<0) # try an exponential distribution?
  {
    warning("VAR[Area] = ",VAR," < 0")
    if(!HDR)
    {
      CI <- c(1,1,1)*MLE
      CI[1] <- stats::qexp(alpha/2,rate=1/min(sqrt(-VAR),MLE))
      CI[3] <- stats::qexp(1-alpha/2,rate=1/max(sqrt(-VAR),MLE))
    }
    else
    { CI <- c(0,0,MLE) * stats::qexp(1-alpha,rate=1/min(sqrt(-VAR),MLE)) }
  }
  else     # regular estimate
  {
    if(HDR) # highest density region
    { CI <- MLE * chisq.hdr(df=DOF,level=level,pow=HDR)/DOF }
    else if(robust) # quantile CIs # not sure how well this works for k<<1
    { CI <- MLE * c(CI.lower(DOF,alpha),stats::qchisq(0.5,DOF)/DOF,CI.upper(DOF,alpha)) }
    else # traditional confidence intervals
    { CI <- MLE * c(CI.lower(DOF,alpha),1,CI.upper(DOF,alpha)) }

    # Normal backup for upper.tail
    if(is.null(VAR)) { VAR <- 2*MLE^2/DOF }
    UPPER <- norm.ci(CI[2],VAR,alpha=alpha)[3]
    # qchisq upper.tail is too small when DOF<<1
    # qchisq bug that no regular use of chi-square/gamma would come across
    if(is.nan(CI[3]) || CI[3]<UPPER) { CI[3] <- UPPER }
  }

  names(CI) <- NAMES.CI
  return(CI)
}


ci.invgauss <- function(est,var,se=sqrt(var),level=0.95,...)
{
  if(missing(var)) { var = se^2 }
  CI <- t(Vectorize(IG.ci)(est,VAR=var,level=level,...))
  if(nrow(CI)==1) { CI <- CI[1,] }
  return(CI)
}


## F-distribution CIs with exact means and variances for the ratio, numerator, and denominator
# E1 == E[numerator]
# VAR1 == VAR[numerator]
# E2 == E[1/denominator]
# VAR2 == VAR[1/denominator]
F.CI <- function(E1,VAR1,E2,VAR2,level=0.95)
{
  EST <- nant(E1*E2,0)
  N1 <- 2*E1^2/VAR1 # chi^2 DOF
  N2 <- nant(2*E2^2/VAR2 + 4,0) # inverse-chi^2 DOF

  if(N1<=0 || N2<=0)
  { CI <- c(0,EST,Inf) }
  else
  {
    BIAS <- N2/(N2-2) # F-distribution mean bias factor
    if(BIAS<=0) { BIAS <- Inf }

    alpha <- (1-level)/2
    CI <- stats::qf(c(alpha,1-alpha),N1,N2)
    CI <- CI / BIAS # debiased ratio corresponding to EST=1
    CI <- c(CI[1],1,CI[2]) * EST
  }

  names(CI) <- NAMES.CI
  return(CI)
}

## log(F) CIs
# E1 == E[numerator]
# VAR1 == VAR[numerator]
# E2 == E[denominator]
# VAR2 == VAR[denominator]
Log.F.CI <- function(E1,VAR1,E2,VAR2,level=0.95)
{
  N1 <- 2*E1^2/VAR1 # chi^2 DOF
  N2 <- 2*E2^2/VAR2 # chi^2 DOF

  alpha <- (1-level)/2
  CI <- numeric(3)
  CI[c(1,3)] <- log( stats::qf(c(alpha,1-alpha),N1,N2) )
  CI[2] <- CI[2] + (log(E1) - log_chi2_bias(N1)) - (log(E2) - log_chi2_bias(N2))

  names(CI) <- NAMES.CI
  return(CI)
}


# inverse Gaussian CIs
IG.ci <- function(mu,VAR,k=VAR/mu^3,level=0.95,precision=1/2)
{
  if(is.na(level))
  {
    CI <- c(2*mu^2/VAR,mu,VAR)
    names(CI) <- c("DOF","est","VAR")
    return(CI)
  }

  if(k==Inf)
  { CI <- c(0,mu,Inf) }
  else if(k>0)
  {
    tol <- .Machine$double.eps^precision
    p <- (1-level)/2
    p <- c(p,1-p)

    CI <- rep(NA_real_,3)
    CI[2] <- mu # mean
    CI[c(1,3)] <- statmod::qinvgauss(p,mean=mu,shape=1/k,tol=min(tol,1E-14),maxit=.Machine$integer.max)
  }
  else
  { CI <- c(mu,mu,mu) }
  names(CI) <- NAMES.CI
  return(CI)
}


# CI that are chi^2 when VAR[M]>>VAR.pop
chisq.IG.ci <- function(M,VAR,w,level=0.95,precision=1/2)
{
  w <- c(1-w,w)

  # approximation that bridges the two with correct mean and variance
  CI.1 <- chisq.ci(M,VAR=VAR,level=level)
  CI.2 <- IG.ci(M,VAR=VAR,level=level,precision=precision)

  CI <- w[1]*CI.1 + w[2]*CI.2

  # TODO make the second term some kind of GIG?
  CI[3] <- nant(CI[3],Inf)

  return(CI)
}


ci.lognorm <- function(est,var,se=sqrt(var),level=0.95,...)
{
  if(missing(var)) { var = se^2 }
  CI <- t(Vectorize(lognorm.ci)(MLE=est,VAR=var,level=level,...))
  if(nrow(CI)==1) { CI <- CI[1,] }
  return(CI)
}

# calculate log-normal confidence intervals from MLE and COV estimates
lognorm.ci <- function(MLE,VAR,level=0.95,alpha=1-level)
{
  # MLE = exp(mu + sigma/2)
  # VAR = (exp(sigma)-1) * MLE^2

  # exp(sigma)-1 = VAR/MLE^2
  # exp(sigma) = 1+VAR/MLE^2
  sigma <- log(1 + VAR/MLE^2)

  # mu+sigma/2 = log(MLE)
  mu <- log(MLE) - sigma/2

  CI <- norm.ci(mu,sigma,alpha=alpha)

  # transform back
  CI <- exp(CI)
  CI[2] <- MLE

  return(CI)
}


ci.norm <- function(est,var,se=sqrt(var),level=0.95,...)
{
  if(missing(var)) { var = se^2 }
  CI <- t(Vectorize(norm.ci)(MLE=est,VAR=var,level=level,...))
  if(nrow(CI)==1) { CI <- CI[1,] }
  return(CI)
}

# normal confidence intervals
norm.ci <- function(MLE,VAR,level=0.95,alpha=1-level)
{
  # z-values for low, ML, high estimates
  z <- stats::qnorm(1-alpha/2)*c(-1,0,1)

  # normal ci
  CI <- MLE + z*sqrt(VAR)

  names(CI) <- NAMES.CI
  return(CI)
}
