# iterpolate vector by continuous index
# vec is a vector, ind is a continuous index
vint <- function(vec,ind,return.ind=FALSE)
{
  n <- length(vec)

  lo <- floor(ind)
  hi <- ceiling(ind)

  # extrapolate
  lo <- max(1,lo)
  hi <- max(2,hi)
  # extrapolate
  lo <- min(n-1,lo)
  hi <- min(n,hi)

  if(return.ind) { return(c(lo,hi)) }

  IND <- abs(vec[c(lo,hi)])
  if(any(IND==Inf))
  {
    IND <- which.max(IND)
    return(vec[c(lo,hi)[IND]])
  }

  # linear interpolation
  vec <- vec[lo] + (vec[hi]-vec[lo])*(ind-lo)

  return(vec)
}
# same thing as above but with a block-vector
mint <- function(mat,ind)
{
  IND <- vint(mat[1,],ind,return.ind=TRUE)
  mat <- mat[,IND[1]] + (mat[,IND[2]]-mat[,IND[1]])*(ind-IND[1])
  return(mat)
}


# confidence interval functions
CI.upper <- Vectorize(function(k,Alpha){stats::qchisq(Alpha/2,k,lower.tail=FALSE)/k})
CI.lower <- Vectorize(function(k,Alpha){stats::qchisq(Alpha/2,k,lower.tail=TRUE)/k})


# calculate chi^2 confidence intervals from MLE and COV estimates
chisq.ci <- function(MLE,COV=NULL,level=0.95,alpha=1-level,DOF=2*MLE^2/COV,robust=FALSE,HDR=FALSE)
{
  #DEBUG <<- list(MLE=MLE,COV=COV,level=level,alpha=alpha,DOF=DOF,robust=robust,HDR=HDR)
  # try to do something reasonable on failure cases
  if(is.nan(DOF) || is.na(DOF)) { DOF <- 0 } # NaN comes from infinite variance divsion
  if(is.na(MLE)) { MLE <- Inf }

  if(DOF==Inf)
  { CI <- c(1,1,1)*MLE }
  else if(DOF==0)
  { CI <- c(0,MLE,Inf) }
  else if(MLE==0)
  { CI <- c(0,0,0) }
  else if(MLE==Inf)
  { CI <- c(Inf,Inf,Inf) }
  else if(!is.null(COV) && COV<0) # try an exponential distribution?
  {
    warning("VAR[Area] = ",COV," < 0")
    if(!HDR)
    {
      CI <- c(1,1,1)*MLE
      CI[1] <- stats::qexp(alpha/2,rate=1/min(sqrt(-COV),MLE))
      CI[3] <- stats::qexp(1-alpha/2,rate=1/max(sqrt(-COV),MLE))
    }
    else
    { CI <- c(0,0,MLE) * stats::qexp(1-alpha,rate=1/min(sqrt(-COV),MLE)) }
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
    if(is.null(COV)) { COV <- 2*MLE^2/DOF }
    UPPER <- norm.ci(CI[2],COV,alpha=alpha)[3]
    # qchisq upper.tail is too small when DOF<<1
    # qchisq bug that no regular use of chi-square/gamma would come across
    if(is.nan(CI[3]) || CI[3]<UPPER) { CI[3] <- UPPER }
  }

  names(CI) <- NAMES.CI
  return(CI)
}


# highest density region (HDR) for chi/chi^2 distribution
# pow=1 gives chi HDR (in terms of chi^2 values)
# pow=2 gives chi^2 HDR
chisq.hdr <- function(df,level=0.95,alpha=1-level,pow=1)
{
  # mode == 0
  if(df <= pow)
  { CI <- c(0,0,stats::qchisq(level,df,lower.tail=TRUE)) } # this goes badly when df<<1
  else # mode == df - pow
  {
    # chi and chi^2 modes
    mode <- df-pow

    # given some chi^2 value under the mode, solve for the equiprobability-density chi^2 value over the mode
    X2 <- function(X1) { if(X1==0){return(Inf)} ; X1 <- -X1/mode ; return( -mode*gsl::lambert_Wm1(X1*exp(X1)) ) }
    # how far off are we from the desired coverage level
    dP <- function(X1) {stats::pchisq(X2(X1),df,lower.tail=TRUE)-stats::pchisq(X1,df,lower.tail=TRUE)-level}
    # solve for X1 to get the desired coverage level
    X1 <- stats::uniroot(dP,c(0,mode),f.lower=alpha,f.upper=-level,extendInt="downX",tol=.Machine$double.eps,maxiter=.Machine$integer.max)$root
    # uniroot is not reliable when root is very near boundary
    if(X1==0)
    {
      # cost function with logit link
      dP2 <- function(z) { dP(mode/(1+exp(z)))^2 }
      X1 <- stats::nlm(dP2,1,ndigit=16,gradtol=.Machine$double.eps,stepmax=100,steptol=.Machine$double.eps,iterlim=.Machine$integer.max)$minimum
      X1 <- mode/(1+exp(X1))
      # nlm is not reliable in some other cases...
    }
    # HDR CI
    CI <- c(X1,mode,X2(X1))
  }
  return(CI)
}


# inverse density function for chi-square distribution
# p == dchisq(x,df)
idchisq <- function(p,df)
{
  # constant
  p <- 2^(df/2) * gamma(df/2) * p
  # unit power
  p <- p^2
  p <- p^(1/(df-2))
  # unit scale
  p <- -p/(df-2)
  # solve
  # x <- c( lamW::lambertW0(p) , lamW::lambertWm1(p) ) ### CRAN check won't like without depends
  # transform back
  x <- -(df-2)*x
  x <- sort(x)
  return(x)
}


# normal confidence intervals
norm.ci <- function(MLE,COV,alpha=0.05)
{
  # z-values for low, ML, high estimates
  z <- stats::qnorm(1-alpha/2)*c(-1,0,1)

  # normal ci
  CI <- MLE + z*sqrt(COV)

  names(CI) <- NAMES.CI
  return(CI)
}

# calculate log-normal confidence intervals from MLE and COV estimates
lognorm.ci <- function(MLE,COV,alpha=0.05)
{
  # log transform of variance
  COV <- COV/MLE^2
  # log transform of point estimate
  MLE <- log(MLE)

  CI <- norm.ci(MLE,COV,alpha=alpha)

  # transform back
  CI <- exp(CI)

  return(CI)
}

# robust central tendency & dispersal estimates with high breakdown threshold
# nstart should probably be increased from default of 1
rcov <- function(x,...)
{
  DIM <- dim(x)
  if(DIM[1]>DIM[2]) { x <- t(x) }

  # first estimate
  AVE <- apply(x,1,stats::median)
  # MAD esitmate of standard deviation (robust)
  MAD <- abs(x-AVE)
  CONST <- 1/stats::qnorm(3/4)
  MAD <- CONST * apply(MAD,1,stats::median)

  # standardize before calculating geometric median
  if(length(AVE)>1)
  {
    STUFF <- Gmedian::GmedianCov(t((x-AVE)/MAD),init=rep(0,length(AVE)),scores=0,nstart=10,...)
    AVE <- c(STUFF$median * MAD + AVE)
    COV <- t(STUFF$covmedian * MAD) * MAD
  }
  else
  { COV <- MAD^2 * diag(length(AVE)) }

  # too many infinities for Gmedian to handle --- fall back to MAD
  NANS <- is.nan(diag(COV)) | is.na(diag(COV))
  if(any(NANS))
  {
    COV[NANS,] <- COV[,NANS] <- 0
    COV[NANS,NANS] <- MAD[NANS]^2
  }

  return(list(median=AVE,COV=COV))
}


# minimally trimmed mean
# could make O(n) without full sort
mtmean <- function(x,lower=-Inf,upper=Inf,func=mean)
{
  x <- sort(x)

  # lower trim necessary
  n <- sum(x<=lower)
  # fallbacks
  if(n==length(x)) { return(lower) }
  if(2*n>=length(x)) { return(x[n+1]) }

  # upper trim necessary
  m <- sum(x>=upper)
  # fallbacks
  if(m==length(x)) { return(upper) }
  if(2*m>=length(x)) { return(x[length(x)-m-1]) }

  n <- max(n,m)
  x <- x[(1+n):(length(x)-n)]

  x <- func(x)
  return(x)
}


# degrees-of-freedom of (proportionally) chi distribution with specified moments
chi.dof <- function(M1,M2,error=1/2)
{
  # solve for chi^2 DOF consistent with M1 & M2
  R <- M1^2/M2 # == 2*pi/DOF / Beta(DOF/2,1/2)^2 # 0 <= R <= 1
  if(R>=1) { return(Inf) } # purely deterministic
  if(R<=0) { return(0) }

  DOF <- M1^2/(M2-M1^2)/2 # initial guess - asymptotic limit
  error <- .Machine$double.eps^error
  ERROR <- Inf
  while(ERROR>=error)
  {
    # current value at guess
    R0 <- 2*pi/DOF/beta(DOF/2,1/2)^2
    # current value of gradient
    G0 <- ( digamma((DOF+1)/2) - digamma(DOF/2) - 1/DOF )*R0
    # correction
    delta <- (R-R0)/G0
    # make sure DOF remains positive
    if(DOF+delta<=0)
    { DOF <- DOF/2 }
    else
    { DOF <- DOF + delta }
    ERROR <- abs(delta)/DOF
  }

  return(DOF)
}


# highest density region of truncated normal distribution
# mu is mode
# lower <= mu <= upper
tnorm.hdr <- function(mu=0,VAR=1,lower=0,upper=Inf,level=0.95)
{
  sd <- sqrt(VAR)
  MASS <- stats::pnorm(upper,mean=mu,sd=sd) - stats::pnorm(lower,mean=mu,sd=sd)
  CDF <- function(a,b) { ( stats::pnorm(b,mean=mu,sd=sd) - stats::pnorm(a,mean=mu,sd=sd) )/MASS }

  CI <- c(lower,mu,upper)
  DIFF <- diff(CI)

  if(DIFF[1]<=DIFF[2] && CDF(lower,mu+DIFF[1])<=level) # we hit the lower boundary
  {
    level <- level - CDF(lower,mu) # upper half mass remaining
    level <- level * MASS # mass remaining (when not truncated)
    CI[3] <- stats::qnorm(level+0.5,mean=mu,sd=sd,lower.tail=TRUE)
  }
  else if(DIFF[2]<=DIFF[1] && CDF(mu-DIFF[2],upper)<=level) # we hit the upper boundary
  {
    level <- level - CDF(mu,upper) # lower half mass remaining
    level <- level * MASS # mass remining (when not truncated)
    CI[1] <- stats::qnorm(level+0.5,mean=mu,sd=sd,lower.tail=FALSE)
  }
  else # no boundary hit
  {
    level <- level/2 # probability mass on each side
    level <- level * MASS # mass on each side (when not truncated)
    CI[1] <- stats::qnorm(level+0.5,mean=mu,sd=sd,lower.tail=FALSE)
    CI[3] <- stats::qnorm(level+0.5,mean=mu,sd=sd,lower.tail=TRUE)
  }

  names(CI) <- NAMES.CI
  return(CI)
}
