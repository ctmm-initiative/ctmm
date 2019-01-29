# base match.arg cannot match NA ???
match.arg <- function(arg,choices,...)
{
  if(is.na(arg) && any(is.na(choices))) { return(NA) }
  else { return(base::match.arg(arg,choices,...)) }
}

# is a package installed?
is.installed <- function(pkg) is.element(pkg, utils::installed.packages()[,1])

# does this thing exist and, if so, is it true
true <- function(x) { !is.null(x) & !is.na(x) & x }

# not in #
"%nin%" <- function(x, table) { match(x, table, nomatch = 0) == 0 }

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
    if(!inverse) { X <- stats::mvfft(X) }
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

.onLoad <- function(...)
{
  # choose FFTW if installed
  if(is.installed("fftw")) { utils::assignInMyNamespace("FFT", FFTW) }
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


# multivariate polygamma function
mpsigamma <- function(x,deriv=0,dim=1)
{
  PSI <- 1 - 1:dim
  PSI <- x + PSI/2
  PSI <- sapply(PSI,function(p) psigamma(p,deriv=deriv))
  PSI <- sum(PSI)
  return(PSI)
}


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
median.list <- function(x,na.rm=FALSE,...)
{
  CLASS <- class(x[[1]])
  utils::getS3method("median",CLASS)(x,na.rm=na.rm,...)
}


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

# forwarding function for list of a particular datatype
writeShapefile.list <- function(object,folder,file=NULL,...)
{
  CLASS <- class(object[[1]])
  utils::getS3method("writeShapefile",CLASS)(object,folder,file=file,...)
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


# generalized covariance from -likelihood derivatives
cov.loglike <- function(hess,grad=rep(0,nrow(hess)))
{
  # if hessian is likely to be positive definite
  if(all(diag(hess)>0))
  {
    COV <- try(PDsolve(hess))
    if(class(COV)=="matrix" && all(diag(COV)>0)) { return(COV) }
  }
  # one of the curvatures is negative or close to negative
  # return something sensible just in case we are on a boundary and this makes sense

  # normalize parameter scales by curvature or gradient (whatever is larger)
  V <- abs(diag(hess))
  V <- sqrt(V)
  V <- pmax(V,abs(grad))

  # don't divide by zero
  # TOL <- .Machine$double.eps * length(V)
  TEST <- V==0
  if(any(TEST)) { V[TEST] <- 1 }

  W <- V %o% V

  grad <- grad/V
  hess <- hess/W

  EIGEN <- eigen(hess)
  values <- EIGEN$values
  if(any(values<=0)) { warning("MLE is near a boundary or optim() failed.") }
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
      # warning("MLE is near a boundary.")
      # (parameter distance to worst parameter * 1/2 / loglike difference to worst parameter)^2
      # values[i] <- 1/grad[i]^2
      # pretty close to the other formula, so just using that
      values[i] <- 1/(2*grad[i])^2
    }
  }

  COV <- array(0,dim(hess))
  SUB <- values<Inf
  if(any(SUB)) # separate out the finite part
  { COV <- COV + Reduce("+",lapply((1:length(grad))[SUB],function(i){ values[i] * outer(vectors[,i]) })) }
  SUB <- !SUB
  if(any(SUB)) # toss out the off-diagonal NaNs
  { COV <- COV + Reduce("+",lapply((1:length(grad))[SUB],function(i){ D <- diag(outer(vectors[,i])) ; D[D>0] <- Inf ; diag(D) })) }

  COV <- COV/W

  return(COV)
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


# prepend to a vector
prepend <- function(x,values,before=1)
{ append(x,values,after=before-1) }


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
  mat <- cbind(mat)
  size <- size - nrow(mat)
  COL <- ncol(mat)
  padding <- array(padding,c(size,COL))
  colnames(padding) <- colnames(mat)

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


# put in a list if not in a list already
listify <- function(x)
{
  if(is.null(x)) { return(x) }
  if(class(x) != "list") { x <- list(x) }
  return(x)
}


# rename elements of an object
rename <- function(object,name1,name2)
{
  IND <- which(names(object)==name1)
  names(object)[IND] <- name2
  return(object)
}


# glue strings together if they different
glue <- function(...)
{
  x <- c(...) # removes NULLs
  x <- unique(x) # removes matchs
  x <- paste(x) # paste different strings together
  return(x)
}
