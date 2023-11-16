# generalized covariance from negative-log-likelihood derivatives
cov.loglike <- function(hess,grad=rep(0,sqrt(length(hess))),tol=.Machine$double.eps,WARN=TRUE)
{
  EXCLUDE <- c("ctmm.boot","cv.like","ctmm.select")

  # in case of bad derivatives, use worst-case numbers
  grad <- nant(grad,Inf)
  hess <- nant(hess,0)

  # if hessian is likely to be positive definite
  if(all(diag(hess)>0))
  {
    COV <- try(PDsolve(hess),silent=TRUE)
    if(class(COV)[1]=="matrix" && all(diag(COV)>0)) { return(COV) }
  }
  # one of the curvatures is negative or close to negative
  # return something sensible just in case we are on a boundary and this makes sense

  # normalize parameter scales by curvature or gradient (whatever is larger)
  V <- abs(diag(hess))
  V <- sqrt(V)
  V <- pmax(V,abs(grad))

  # don't divide by zero
  TEST <- V<=tol
  if(any(TEST)) { V[TEST] <- 1 }

  W <- V %o% V

  grad <- nant(grad/V,1)
  hess <- nant(hess/W,1)

  # clamp off-diagonal
  MAX <- sqrt(abs(diag(hess)))
  MAX <- MAX %o% MAX
  hess[] <- clamp(hess,-MAX,MAX)

  EIGEN <- eigen(hess)
  values <- EIGEN$values
  if(any(values<=0))
  {
    # shouldn't need to warn if using ctmm.select, but do if using ctmm.fit directly
    if(WARN)
    {
      N <- sys.nframe()
      if(N>=2)
      {
        for(i in 2:N)
        {
          CALL <- deparse(sys.call(-i))[1]
          CALL <- sapply(EXCLUDE,function(E){grepl(E,CALL,fixed=TRUE)})
          CALL <- any(CALL)
          if(CALL)
          {
            WARN <- FALSE
            break
          }
        } # 2:N
      } # N >= 2
    }
    # warn if weren't using ctmm.select
    if(WARN) { warning("MLE is near a boundary or optimizer failed.") }
  }
  values <- clamp(values,0,Inf)
  vectors <- EIGEN$vectors

  # transform gradient to hess' coordinate system
  grad <- t(vectors) %*% grad

  # generalized Wald-like formula with zero-curvature limit
  # VAR ~ square change required to decrease log-likelihood by 1/2
  for(i in 1:length(values))
  {
    DET <- values[i]+grad[i]^2

    if(values[i]==0.0) # Wald limit of below
    { values[i] <- 1/(2*grad[i])^2 }
    else if(DET>=0.0) # Wald
    { values[i] <- ((sqrt(DET)-grad[i])/values[i])^2 }
    else # minimum loglike? optim probably failed or hit a boundary
    {
      # (parameter distance to worst parameter * 1/2 / loglike difference to worst parameter)^2
      # values[i] <- 1/grad[i]^2
      # pretty close to the other formula, so just using that
      values[i] <- 1/(2*grad[i])^2
    }
  }
  values <- nant(values,0)

  COV <- array(0,dim(hess))
  values <- nant(values,Inf) # worst case NaN fix
  SUB <- values<Inf
  if(any(SUB)) # separate out the finite part
  { COV <- COV + Reduce("+",lapply((1:length(grad))[SUB],function(i){ values[i] * outer(vectors[,i]) })) }
  SUB <- !SUB
  if(any(SUB)) # toss out the off-diagonal NaNs
  { COV <- COV + Reduce("+",lapply((1:length(grad))[SUB],function(i){ D <- diag(outer(vectors[,i])) ; D[D>0] <- Inf ; diag(D,length(D)) })) }

  COV <- COV/W

  return(COV)
}


# confidence interval functions
CI.upper <- Vectorize(function(k,Alpha){stats::qchisq(Alpha/2,k,lower.tail=FALSE)/k})
CI.lower <- Vectorize(function(k,Alpha){stats::qchisq(Alpha/2,k,lower.tail=TRUE)/k})


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
norm.ci <- function(MLE,VAR,level=0.95,alpha=1-level)
{
  # z-values for low, ML, high estimates
  z <- stats::qnorm(1-alpha/2)*c(-1,0,1)

  # normal ci
  CI <- MLE + z*sqrt(VAR)

  names(CI) <- NAMES.CI
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


# beta distributed CI given mean and variance estimates
beta.ci <- function(MLE,VAR,level=0.95,alpha=1-level)
{
  MLE <- nant(MLE,0)
  VAR <- nant(VAR,Inf)
  n <- MLE*(1-MLE)/VAR - 1
  if(n<=0)
  { CI <- c(0,MLE,1) }
  else
  {
    a <- n * MLE
    b <- n * (1-MLE)
    CI <- stats::qbeta(c(alpha/2,0.5,1-alpha/2),a,b)
    CI[2] <- MLE # replace median with mean
  }
  names(CI) <- NAMES.CI
  return(CI)
}


# Binomial CDF defined for effective size (n) and outcome fraction (q)
pfbinom <- function(q,size,prob,lower.tail=TRUE,log.p=FALSE)
{
  x <- (1-prob)/prob * (q+1/size)/(1-q)
  df1 <- 2*size*(1-q)
  df2 <- 2*size*(q+1/size)
  stats::pf(x,df1,df2,lower.tail=lower.tail,log.p=log.p)
}


# Binomial quantile function defined for effective size (n) and outcome fraction (q)
qfbinom <- function(p,size,prob)
{
  parscale <- c(p,1-p,prob,1-prob)
  parscale <- min(parscale[parscale>0])

  # solve p==pfbinom(q)
  TEST <- pfbinom(1/2,size,prob) # median quantile
  ## prevent underflow... not sure if this is necessary
  if(p<=TEST) # below median -- lower-tail probability optimization
  {
    lower.tail <- TRUE
    par <- 1/4
    lower <- 0
    upper <- 1/2
  }
  else # above median -- upper-tail probability optimization
  {
    lower.tail <- FALSE
    p <- 1-p # upper-tail probability
    par <- 3/4
    lower <- 1/2
    upper <- 1
  }
  cost <- function(q) { (p-pfbinom(q,size,prob,lower.tail=lower.tail))^2 }
  FIT <- optimizer(par,cost,lower=lower,upper=upper,control=list(parscale=parscale))
  q <- FIT$par

  return(q)
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
chi.dof <- function(M1,M2,precision=1/2)
{
  # DEBUG <<- list(M1=M1,M2=M2,precision=precision)
  error <- .Machine$double.eps^precision

  # solve for chi^2 DOF consistent with M1 & M2
  R <- M1^2/M2 # == 2*pi/DOF / Beta(DOF/2,1/2)^2 # 0 <= R <= 1
  if(M1==0 && M2==0) { R <- 1 }
  if(1-R <= 0) { return(Inf) } # purely deterministic
  if(R <= 0) { return(0) }

  DOF <- M1^2/(M2-M1^2)/2 # initial guess - asymptotic limit
  ERROR <- Inf
  while(ERROR>=error)
  {
    DOF.old <- DOF
    ERROR.old <- ERROR

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

    if(ERROR>ERROR.old) { return(DOF.old) }
  }

  return(DOF)
}

# variance of chi variable, given dof
chi.var <- function(DOF,M1=1)
{
  R <- DOF

  fn <- function(DOF){ 2*pi/DOF/beta(DOF/2,1/2)^2 }

  # Laurent expansion (large DOF)
  MAX <- 26 # switch over point in numerical accuracy
  coef <- c( 1, -(1/2), 1/8, 1/16, -(5/128), -(23/256), 53/1024, 593/2048, -(5165/32768), -(110123/65536), 231743/262144 )
  pn <- Vectorize( function(DOF) { series(1/DOF,coef) } )
  SUB <- DOF>=MAX
  if(any(SUB)) { R[SUB] <- pn(DOF[SUB]) }

  # Taylor expansion (small DOF)
  MIN <- 0.000002
  coef <- c(0, 1/2, -log(2), pi^2/24 + log(2)^2) * pi # further terms require gsl::zeta()
  pn <- Vectorize( function(DOF) { series(DOF,coef) } )
  SUB <- DOF<=MIN
  if(any(SUB)) { R[SUB] <- pn(DOF[SUB]) }

  SUB <- DOF>MIN & DOF<MAX
  if(any(SUB)) { R[SUB] <- fn(DOF[SUB]) }

  VAR <- (1/R-1)*M1^2
  # VAR <- nant(VAR,1/DOF)

  return(VAR)
}

# relative bias in chi estimate, given unbiased chi^2 estimate
# 1 is no bias
chi.bias <- function(DOF)
{
  fn <- function(DOF) { sqrt(2/DOF)*exp(lgamma((DOF+1)/2)-lgamma(DOF/2)) }
  # Laurent expansion
  coef <- c(1, -(1/4), 1/32, 5/128, -(21/2048), -(399/8192), 869/65536, 39325/262144, -(334477/8388608), -(28717403/33554432), 59697183/268435456)
  pn <- Vectorize( function(DOF) { series(1/DOF,coef) } )
  # switch over
  MAX <- 45

  BIAS <- rep(1,length(DOF))
  SUB <- DOF<MAX
  if(any(SUB)) { BIAS[SUB] <- fn(DOF[SUB]) }
  SUB <- DOF>=MAX
  if(any(SUB)) { BIAS[SUB] <- pn(DOF[SUB]) }

  BIAS <- ifelse(DOF==0,0,BIAS)
  # BIAS <- ifelse(DOF==Inf,1,BIAS)
  return(BIAS)
}

# (scaled) chi^2 degrees-of-freedom from median and interquartile range
chisq.dof <- function(MED,IQR,alpha=0.25)
{
  if(IQR==0) { return(Inf) }
  if(IQR==Inf) { return(0) }
  if(MED==0) { return(0) }

  cost <- function(nu,zero=0)
  {
    if(nu==0) { return(Inf) }
    if(nu==Inf) { return(Inf) }

    Q1 <- stats::qchisq(1/4,df=nu)/nu # Q1/mean
    M <- stats::qchisq(1/2,df=nu)/nu # median/mean
    Q2 <- stats::qchisq(3/4,df=nu)/nu # Q2/mean
    R <- M/(Q2-Q1) # IQR/MED from nu
    r <- MED/IQR # IQR/MED from data

    return((R-r)^2)
  }

  # initial guess (asymptotic relation)
  nu <- (MED/IQR)^2 / 0.2747632133101263
  R <- optimizer(nu,cost,lower=0,upper=Inf,control=list(parscale=1))

  return(R$par)
}


# highest density region of truncated normal distribution
# mu is mode
# lower <= mu <= upper
tnorm.hdr <- function(mu=0,VAR=1,lower=0,upper=Inf,level=0.95)
{
  sd <- sqrt(VAR)
  MASS <- stats::pnorm(upper,mean=mu,sd=sd) - stats::pnorm(lower,mean=mu,sd=sd)
  CDF <- function(a,b) { nant( (stats::pnorm(b,mean=mu,sd=sd)-stats::pnorm(a,mean=mu,sd=sd))/MASS ,1) }

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

# Dirac-delta to inverse-Gaussian ratio in sampling distribution
# par = (mu,k=1/lambda)
# VAR = VAR[mu]
DD.IG.ratio <- function(par,VAR,n)
{
  mu <- par[1]
  # population moments
  theta <- 1/prod(par[1:2]); rho <- 1
  M.pop <- mu * (besselK(theta,rho/2-1,expon.scaled=TRUE)/besselK(theta,rho/2,expon.scaled=TRUE))
  VAR.pop <- mu^2 * (besselK(theta,rho/2-2,expon.scaled=TRUE)/besselK(theta,rho/2,expon.scaled=TRUE)) - M.pop^2
  # chi^2 VAR.pop==0
  # IG when VAR==VAR.pop/n
  DRATIO <- VAR.pop/VAR # (0,n) : (chi^2,IG)
  DRATIO <- clamp(DRATIO,0,n) / n # (0,1) : (chi^2,IG)
  DRATIO
}


# generalized inverse Gaussian CIs
# GIG.CI <- function(mu,k,rho,level=0.95,precision=1/2)
# {
#   tol <- .Machine$double.eps^precision
#   p <- (1-level)/2
#   p <- c(p,1-p)
#
#   chi <- theta*eta
#   psi <- theta/eta
#   lambda <- -rho/2
#
#   CI <- rep(NA_real_,3)
#   CI[2] <- eta*(besselK(theta,1-rho/2,expon.scaled=TRUE)/besselK(theta,-rho/2,expon.scaled=TRUE))
#   CI[c(1,3)] <- GeneralizedHyperbolic::qgig(p,chi=chi,psi=psi,lambda=lambda,method="integrate",ibfTol=min(tol,.Machine$double.eps^.85),uniTol=min(tol,1E-7))
#   names(CI) <- NAMES.CI
#   return(CI)
# }


# F-distribution CIs with exact means and variances for the ratio, numerator, and denominator
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


# reduced chi^2 log likelihood for s=1
loglike.chisq <- function(sigma,dof,constant=FALSE)
{
  df2 <- dof/2
  R <- - df2*(log(sigma)+1/sigma)
  if(constant) { R <- R + df2*log(df2) - lgamma(df2) }
  return(R)
}


#
qmvnorm <- function(p,dim=1,tol=1/2)
{
  tol <- .Machine$double.eps^tol
  alpha <- 1-p

  if(dim==1)
  { z <- stats::qnorm(1-alpha/2) }
  else if(dim==2)
  { z <- sqrt(-2*log(alpha)) }
  else
  {
    if(dim==3)
    {
      # distribution function
      p.fn <- function(z) { stats::pchisq(z^2,df=1) - sqrt(2/pi)*z*exp(-z^2/2) }
      # density function (derivative)
      dp.fn <- function(z) { sqrt(2/pi)*z^2*exp(-z^2/2) }
      # initial guess # dimensional extrapolation
      z <- 2*qmvnorm(p,2) - qmvnorm(p,1)
    }

    ERROR <- Inf
    while(ERROR>tol)
    {
      # p == p.fn(z) + dp.fn(z)*dz + O(dz^2)
      dz <- (p-p.fn(z))/dp.fn(z)
      z <- z + dz
      ERROR <- abs(dz)
    }
  }

  return(z)
}
