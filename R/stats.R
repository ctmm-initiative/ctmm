# statistical mode
Mode <- function(x)
{
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
  mean(ux)
}


# confidence interval functions
CI.upper <- Vectorize(function(k,Alpha){stats::qchisq(Alpha/2,k,lower.tail=FALSE)/k})
CI.lower <- Vectorize(function(k,Alpha){stats::qchisq(Alpha/2,k,lower.tail=TRUE)/k})


# calculate chi^2 confidence intervals from MLE and COV estimates
chisq.ci <- function(MLE,COV=NULL,level=0.95,alpha=1-level,DOF=2*MLE^2/COV,HDR=FALSE)
{
  # try to do something reasonable on failure cases
  if(is.nan(DOF)) { DOF <- 0 } # this comes from infinite variance divsion
  if(DOF==0)
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
    if(!HDR) # traditional confidence intervals
    { CI <- MLE * c(CI.lower(DOF,alpha),1,CI.upper(DOF,alpha)) }
    else # highest density region
    { CI <- MLE * chisq.hdr(df=DOF,level=level,pow=HDR)/DOF }

    # Normal backup for upper.tail
    if(is.null(COV)) { COV <- 2*MLE^2/DOF }
    UPPER <- norm.ci(CI[2],COV,alpha=alpha)[3]
    # qchisq upper.tail is too small when DOF<<1
    # probably an R bug that no regular use of chi-square/gamma would come across
    if(CI[3]<UPPER) { CI[3] <- UPPER }
  }

  names(CI) <- c("low","ML","high")
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

  names(CI) <- c("low","ML","high")
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
