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
chisq.ci <- function(MLE,COV=NULL,alpha=0.05,DOF=2*MLE^2/COV,fast=TRUE)
{
  # try to do something reasonable on failure cases
  if(DOF==0)
  { CI <- c(0,MLE,Inf) }
  else if(MLE==0)
  { CI <- c(0,0,0) }
  else if(MLE==Inf)
  { CI <- c(Inf,Inf,Inf) }
  else if(!is.null(COV) && COV<0) # try an exponential distribution?
  {
    warning("VAR[Area] = ",COV," < 0")
    CI <- c(1,1,1)*MLE
    CI[1] <- stats::qexp(alpha/2,rate=1/min(sqrt(-COV),MLE))
    CI[3] <- stats::qexp(1-alpha/2,rate=1/max(sqrt(-COV),MLE))
  }
  else     # regular estimate
  {
    # traditional confidence intervals
    if(fast)
    { CI <- MLE * c(CI.lower(DOF,alpha),1,CI.upper(DOF,alpha)) }
    else # more probable confidence intervals
    {
      CI <- MLE * q2chisq(1-alpha,DOF)/DOF
      CI[3] <- CI[2]
      CI[2] <- MLE
    }

    # qchisq upper.tail is too small when DOF<<1
    # probably an R bug that no regular use of chi-square/gamma would come across
    if(is.null(COV)) { COV <- 2*MLE^2/DOF }
    # Normal backup for upper.tail
    UPPER <- norm.ci(MLE,COV,alpha=alpha)[3]
    if(CI[3]<UPPER) { CI[3] <- UPPER }
  }

  names(CI) <- c("low","ML","high")
  return(CI)
}


# proper 2-sided quantile function for chi-squared distribution
q2chisq <- function(p,df)
{
  # mode == 0
  if(df <= 2) { CI <- c(0,0,stats::qchisq(p,df,lower.tail=TRUE)) }
  else # mode == df
  {
    # conventional CIs
    x1 <- stats::qchisq((1-p)/2,df,lower.tail=TRUE)
    x2 <- stats::qchisq((1-p)/2,df,lower.taul=FALSE)
    # mismatched density
    p1 <- stats::dchisq(x1,df)
    p2 <- stats::dchisq(x2,df)
    # average density as first guess
    p0 <- sqrt(p1*p2)

    # MORE TO COME !!!

    # start loop

    # quantiles for this density
    x <- idchisq(p0,df)

    # total probability for these quantiles

    # newton iteration towards target probability

    # end loop

    CI <- c(x1,df-2,x2)
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

  return(list(median=AVE,COV=COV))
}
