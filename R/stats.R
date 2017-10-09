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
chisq.ci <- function(MLE,COV=NULL,alpha=0.05,DOF=2*MLE^2/COV)
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
    CI <- MLE * c(CI.lower(DOF,alpha),1,CI.upper(DOF,alpha))

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
