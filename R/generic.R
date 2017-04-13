# is a package installed?
is.installed <- function(pkg) is.element(pkg, utils::installed.packages()[,1]) 

# generic FFT functions
FFT <- function(X,inverse=FALSE)
{
  if(is.null(ncol(X)) || is.na(ncol(X)))
  { 
    if(!inverse) { X <- stats::fft(X) }
    else { X <- stats::fft(X,inverse=TRUE)/length(X) }
  }
  else
  {
    if(!inverse) { X <- mvfft(X) }
    else { X <- stats::mvfft(X,inverse=TRUE)/nrow(X) }
  }
  
  return(X)
}

# fastest FFT functions... don't use on integers
FFTW <- function(X,inverse=FALSE)
{
  if(is.null(ncol(X)) || is.na(ncol(X)))
  { 
    if(!inverse) { X <- fftw::FFT(X) }
    else { X <- fftw::IFFT(X) }
  }
  else
  {
    if(!inverse) { X <- sapply(1:ncol(X),function(j){ fftw::FFT(X[,j]) }) }
    else { X <- sapply(1:ncol(X),function(j){ fftw::IFFT(X[,j]) }) }
  }
  
  return(X)
}

# parallel functions
detectCores <- parallel::detectCores
mclapply <- parallel::mclapply

.onLoad <- function(...)
{
  # choose FFTW if installed
  if(is.installed("fftw")) { utils::assignInMyNamespace("FFT", FFTW) }

  # don't try to fork if winblows
  if(.Platform$OS.type=="windows")
  { 
    utils::assignInMyNamespace("detectCores", function(...) { 1 })
    utils::assignInMyNamespace("mclapply", function(X,FUN,mc.cores=1,...) { lapply(X,FUN,...) })
  }
}
.onAttach <- .onLoad

IFFT <- function(X,plan=NULL) { FFT(X,inverse=TRUE) }

composite <- function(n) { 2^ceiling(log(n,2)) }

# sinc functions
sinc <- Vectorize( function(x)
{
  if(x==0)
  { return(1) }
  else
  { return(sin(x)/x) }
} )

sinch <- Vectorize( function(x)
{
  if(x==0)
  { return(1) }
  else
  { return(sinh(x)/x) }
} )


# forwarding function for list of a particular datatype
zoom.list <- function(x,...)
{
  CLASS <- class(x[[1]])
  #utils::getS3method("zoom",CLASS)(x,...)
  methods::getMethod("zoom",signature=CLASS)(x,...)
}
methods::setMethod("zoom",signature(x="list"), function(x,...) zoom.list(x,...))


# forwarding function for list of a particular datatype
mean.list <- function(x,...)
{
  CLASS <- class(x[[1]])
  utils::getS3method("mean",CLASS)(x,...)
}
#methods::setMethod("mean",signature(x="list"), function(x,...) mean.list(x,...))


# forwarding function for list of a particular datatype
plot.list <- function(x,...)
{
  CLASS <- class(x[[1]])
  utils::getS3method("plot",CLASS)(x,...)
}
#methods::setMethod("plot",signature(x="list"), function(x,...) plot.list(x,...))

# forwarding function for list of a particular datatype
summary.list <- function(object,...)
{
  CLASS <- class(object[[1]])
  utils::getS3method("summary",CLASS)(object,...)
}


# replace NA elements
na.replace <- function(x,rep)
{
  REP <- is.na(x)
  x[REP] <- rep[REP]
  return(x)
}


# parity tests
is.even <- Vectorize(function(x) {x %% 2 == 0})

is.odd <- Vectorize(function(x) {x %% 2 != 0})


# 2D rotation matrix
rotate <- function(theta)
{
  COS <- cos(theta)
  SIN <- sin(theta)
  R <- rbind( c(COS,-SIN), c(SIN,COS) )
  return(R)
}


# statistical mode
Mode <- function(x)
{
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
  mean(ux)
}


# generalized covariance from -likelihood derivatives
cov.loglike <- function(hess,grad=rep(0,nrow(hess)))
{
  # if hessian is likely to be positive definite
  if(all(diag(hess)>0))
  {
    COV <- try(PDsolve(hess))
    if(class(COV)=="matrix") { return(COV) }
  }
  # one of the curvatures is negative or close to negative
  # return something sensible just in case we are on a boundary and this makes sense
  
  # normalize parameter scales by curvature or gradient (whatever is larger)
  V <- abs(diag(hess))
  V <- sqrt(V)
  V <- pmax(V,abs(grad))
  W <- V %o% V
  
  grad <- grad/V
  hess <- hess/W
  
  EIGEN <- eigen(hess)
  values <- EIGEN$values
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
    { values[i] <- min(((c(1,-1)*sqrt(DET)-grad[i])/values[i])^2) }
    else # minimum loglike? optim probably failed or hit a boundary
    {
      warning("Likelihood has not been maximized. MLE could be near a boundary.")
      # distance to worst parameter
      values[i] <- (grad[i]/values[i])^2
      # not sure what else to do...
    }
  }
  
  COV <- vectors %*% diag(values,length(values)) %*% t(vectors)
  COV <- COV/W
    
  return(COV)
}

# confidence interval functions
CI.upper <- Vectorize(function(k,Alpha){stats::qchisq(Alpha/2,k,lower.tail=FALSE)/k})
CI.lower <- Vectorize(function(k,Alpha){stats::qchisq(Alpha/2,k,lower.tail=TRUE)/k})


# calculate chi^2 confidence intervals from MLE and COV estimates
chisq.ci <- function(MLE,COV=NULL,alpha=0.05,DOF=2*MLE^2/COV)
{
  # try to do something reasonable on failure cases
  if(MLE==0)
  { CI <- c(0,0,0) }
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


# last element of array
last <- function(vec) { vec[length(vec)] }
first <- function(vec) { vec[1] }
# assign to last element... doesn't work
# "last<-" <- function(vec,ass)
# {
#   vec[length(vec)] <- ass
#   return(vec)
# }

# CLAMP A NUMBER
clamp <- Vectorize(function(num,min=0,max=1) { if(num<min) {min} else if(num<max) {num} else {max} })


# PAD VECTOR
pad <- function(vec,size,padding=0,side="right")
{
  # this is now the pad length instead of total length
  size <- size - length(vec)
  padding <- array(padding,size)
  
  if(side=="right"||side=="r")
  { vec <- c(vec,padding) }
  else if(side=="left"||side=="l")
  { vec <- c(padding,vec) }
  
  return(vec)
}

# row pad for data frames / matrices
rpad <- function(mat,size,padding=0,side="right")
{
  size <- size - nrow(mat)
  COL <- ncol(mat)
  padding <- array(padding,c(size,COL))
    
  if(side=="right"||side=="r")
  { mat <- rbind(mat,padding) }
  else if(side=="left" || side=="l")
  { mat <- rbind(padding,mat) }
  
  return(mat)
}

#remove rows and columns by name
rm.name <- function(object,name)
{
  object[!rownames(object) %in% name,!colnames(object) %in% name] 
}
