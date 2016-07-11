# existing S4 generic functions
methods::setGeneric("raster", getGeneric("raster", package="raster"))
methods::setGeneric("zoom", getGeneric("zoom", package="raster"))

# existing functions -> S4 generics
# this doesn't work
#methods::setGeneric("SpatialPoints",package="sp",signature=signature("coords",...))
#methods::setGeneric("SpatialPolygonsDataFrame",package="sp",signature="Sr")

# existing funtions -> S3 generics
# this works but is masked if you load sp
#SpatialPoints <- function(object,...) UseMethod("SpatialPoints")
#SpatialPoints.matrix <- function(object,...) sp::SpatialPoints(coords=object,...)
#SpatialPoints.data.frame <- function(object,...) sp::SpatialPoints(coords=object,...)

#SpatialPolygonsDataFrame <- function(object,...) UseMethod("SpatialPolygonsDataFrame")
#SpatialPolygonsDataFrame.SpatialPolygons <- function(object,...) sp::SpatialPolygonsDataFrame(Sr=object,...)

# new S3 generic functions
writeShapefile <- function(object,...) UseMethod("writeShapefile")


# is a package installed?
is.installed <- function(pkg) is.element(pkg, utils::installed.packages()[,1]) 

# generic FFT functions
FFT <- function(X,inverse=FALSE)
{
  if(is.null(ncol(X)))
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
  if(is.null(ncol(X)))
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

# choose FFTW if installed
.onLoad <- function(...)
{
  if(is.installed("fftw")) { utils::assignInMyNamespace("FFT", FFTW) }

  #if(!isNamespaceLoaded("raster")) { assign("zoom",function(x,...) UseMethod("zoom"),envir=envir) }
}
.onAttach <- .onLoad

IFFT <- function(X,plan=NULL) { FFT(X,inverse=TRUE) }

composite <- function(n) { 2^ceiling(log(n,2)) }

##### det shouldn't fail because R dropped indices
det.numeric <- function(x,...) { x }
determinant.numeric <- function(x,logarithm=TRUE,...)
{
  SIGN <- sign(x)
  if(logarithm)
  { x <- log(abs(x)) }
  
  RESULT <- list(modulus=x,sign=SIGN)
  attr(RESULT$modulus,"logarithm") <- logarithm

  class(RESULT) <- "det"
  
  return(det)
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

# indices where a condition is met
where <- function(x)
{
  (1:length(x))[x]
}


# statistical mode
Mode <- function(x)
{
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
  mean(ux)
}

# adjoint of matrix
Adj <- function(M) { t(Conj(M)) }

# Hermitian part of matrix
He <- function(M) { (M + Adj(M))/2 }

# Positive definite solver
PDsolve <- function(M)
{
  # symmetrize
  M <- He(M)
  
  # rescale
  W <- abs(diag(M))
  W <- sqrt(W)
  W <- W %o% W
  
  # now a correlation matrix that is easy to invert
  M <- M/W
  M <- qr.solve(M)
  M <- M/W

  # symmetrize
  M <- He(M)

  return(M)
}

# sqrtm fails on 1x1 matrix
# I cannot figure out how to make this "Note" go away!
# x <- Matrix::Matrix(x,sparse=FALSE,doDiag=FALSE)
sqrtm <- function(x)
{
  DIM <- dim(x)
  if(all(DIM==c(1,1)))
  { return ( sqrt(x) ) }
  else
  { return( expm::sqrtm(x) ) }
}


# generalized covariance from likelihood derivatives
cov.loglike <- function(hess,grad)
{
  # if hessian is likely to be positive definite
  if(all(diag(hess)>0))
  {
    COV <- try(PDsolve(hess))
    if(class(COV)=="matrix") { return(COV) }
  }
  # one of the curvatures is negative
  # return something sensible just in case we are on a boundary and this makes sense
  
  # normalize parameter scales by curvature or gradient (whatever is larger)
  V <- abs(diag(hess))
  V <- sqrt(V)
  V <- sapply(1:length(grad), function(i) { max(V[i],abs(grad[i])) } )
  W <- V %o% V
  
  grad <- grad/V
  hess <- hess/W
  
  EIGEN <- eigen(hess)
  values <- EIGEN$values
  vectors <- EIGEN$vectors

  # transform gradient to hess' coordinate system
  grad <- t(vectors) %*% grad

  # generalized Wald-like formula with zero-curvature limit
  for(i in 1:length(values))
  {
    DET <- values[i]+grad[i]^2
    
    if(values[i]==0.0) # Wald limit of below
    { values[i] <- 1/(2*grad[i])^2 }
    else if(DET>=0.0) # Wald
    { values[i] <- min(((c(1,-1)*sqrt(DET)-grad[i])/values[i])^2) }
    else # minimum loglike? optim probably failed
    {
      warning("Likelihood has not been maximized.")
      values[i] <- (grad[i]/values[i])^2
    }
  }
  
  COV <- vectors %*% diag(values,length(values)) %*% t(vectors)
  COV <- COV/W
    
  return(COV)
}

# confidence interval functions
CI.upper <- Vectorize(function(k,Alpha){qchisq(Alpha/2,k,lower.tail=FALSE)/k})
CI.lower <- Vectorize(function(k,Alpha){qchisq(Alpha/2,k,lower.tail=TRUE)/k})


# calculate chi^2 confidence intervals from MLE and COV estimates
chisq.ci <- function(MLE,COV=NULL,alpha=0.05,DOF=2*MLE^2/COV)
{
  CI <- MLE * c(CI.lower(DOF,alpha),1,CI.upper(DOF,alpha))
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

# CLAMP A NUMBER
clamp <- Vectorize(function(num,min=0,max=1) { if(num<min) {min} else if(num<max) {num} else {max} })


# Positive definite part of matrix
PDclamp <- function(M)
{ 
  
  # singular value decomposition method
  M <- svd(M)
  M$d <- clamp(M$d,max=Inf) # toss out small negative values
  M$u <- (M$u + M$v)/2 # symmetrize
  
  M <- lapply(1:length(M$d),function(i){M$d[i]*(M$u[i,]%o%Conj(M$u[i,]))})
  M <- Reduce("+",M)
  
  return(M)

  # simple method
  M[1,1] <- clamp(M[1,1],0,Inf)
  M[2,2] <- clamp(M[2,2],0,Inf)
  m <- sqrt(M[1,1]*M[2,2])
  M[1,2] <- clamp(M[1,2],-m,m)
  M[2,1] <- M[1,2]
  
  return(M)
  
}


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


# CHOOSE BEST UNITS FOR A LIST OF DATA
unit <- function(data,dimension,thresh=1,concise=FALSE)
{
  if(dimension=="length")
  {
    name.list <- c("meters","kilometers")
    abrv.list <- c("m","km")
    scale.list <- c(1,1000)
  }
  else if(dimension=="area")
  {
    name.list <- c("square meters","hectares","square kilometers")
    abrv.list <- c("m^2","hm^2","km^2")
    scale.list <- c(1,100^2,1000^2) 
  }
  else if(dimension=="time")
  {
    name.list <- c("seconds","minutes","hours","days","months","years")
    abrv.list <- c("sec","min","hr","day","mon","yr")
    scale.list <- c(1,60*c(1,60*c(1,24*c(1,29.53059,365.24))))
  }
  else if(dimension=="speed")
  {
    name.list <- c("meters/day","kilometers/day")
    abrv.list <- c("m/day","km/day")
    scale.list <- c(1,1000)/(60*60*24)
  }
  else if(dimension=="diffusion")
  {
    name.list <- c("square meters/day","hectares/day","square kilometers/day")
    abrv.list <- c("m^2/day","hm^2/day","km^2/day")
    scale.list <- c(1,100^2,1000^2)/(60*60*24)
  }
    
  max.data <- max(abs(data))
  
  if(concise) { name.list <- abrv.list }
  
  # choose most parsimonious units
  I <- max.data > thresh * scale.list
  if(any(I))
  {
    I <- (1:length(I))[I]
    I <- last(I)
  }
  else { I <- 1 }
  
  name <- name.list[I]
  scale <- scale.list[I]
  
  return(list(scale=scale,name=name))
}


# convert units
setUnits <- function(arg1,arg2)
{
  return(arg1 %#% arg2) 
}


# convert units
`%#%` <- function(arg1,arg2)
{
  # convert to si units
  if(is.numeric(arg1))
  {
    num <- arg1
    name <- arg2
    pow <- +1
  }
  else # convert from si units
  {
    num <- arg2
    name <- arg1
    pow <- -1
  }
  
  alias <- list()
  scale <- c()
  
  add <- function(a,s)
  {
    n <- length(alias)
    alias[[n+1]] <<- a 
    scale[n+1] <<- s
  }
  
  # TIME
  add(c("s","s.","sec","sec.","second","seconds"),1)
  add(c("min","min.","minute","minutes"),60)
  add(c("h","h.","hr","hr.","hour","hours"),60^2)
  add(c("day","days"),24*60^2)
  add(c("week","weeks"),7*24*60^2)
  add(c("month","months"),365.24/12*7*24*60^2)
  add(c("yr","yr.","year","years"),365.24*7*24*60^2)
  
  # Distance conversions
  add(c("m","m.","meter","meters"),1)
  add(c("km","km.","kilometer","kilometers"),1000)
  
  # Area conversions
  add(c("m^2","m.^2","meter^2","meters^2","square meter","square meters","meter squared","meters squared"),1)
  add(c("ha","hectare","hectares"),100^2)
  add(c("km^2","km.^2","kilometer^2","kilometers^2","square kilometer","square kilometers","kilometer squared","kilometers squared"),1000^2)
  
  for(i in 1:length(alias))
  {
    if(name %in% alias[[i]]) { return(num*scale[i]^pow) }
  }
}


## rescale the units of dimensionful parameters
unit.ctmm <- function(CTMM,length=1,time=1)
{
  if(length(CTMM$tau)>0){ CTMM$tau <- CTMM$tau/time }
  CTMM$circle <- CTMM$circle/time
  
  # all means scale with length the same way... but not time
  CTMM$mu <- CTMM$mu/length
  drift <- get(CTMM$mean)
  CTMM <- drift@scale(CTMM,time)
  
  CTMM$error <- CTMM$error/length
  CTMM$sigma <- CTMM$sigma/length^2
  CTMM$sigma@par["area"] <- CTMM$sigma@par["area"]/length^2
  
  # variance -> diffusion adjustment
  if(!CTMM$range)
  {
    CTMM$sigma <- CTMM$sigma*time
    CTMM$sigma@par["area"] <- CTMM$sigma@par["area"]*time
  }
  
  if(!is.null(CTMM$COV))
  {
    CTMM$COV.mu <- CTMM$COV.mu/length^2
    
    CTMM$COV["area",] <- CTMM$COV["area",]/length^2
    CTMM$COV[,"area"] <- CTMM$COV[,"area"]/length^2
    
    if(!CTMM$range)
    {
      CTMM$COV["area",] <- CTMM$COV["area",]*time
      CTMM$COV[,"area"] <- CTMM$COV[,"area"]*time
    }
    
    tau <- CTMM$tau
    tau <- tau[tau<Inf]
    if(length(tau))
    {
      tau <- names(tau)
      tau <- paste("tau",tau)

      CTMM$COV[tau,] <- CTMM$COV[tau,]/time
      CTMM$COV[,tau] <- CTMM$COV[,tau]/time
    }
    
    if(CTMM$circle)
    {
      CTMM$COV["circle",] <- CTMM$COV["circle",]/time
      CTMM$COV[,"circle"] <- CTMM$COV[,"circle"]/time
    }
  }

  return(CTMM)
}

######################
unit.UD <- function(UD,length=1)
{
  UD$r <- lapply(UD$r,function(x){ x/length })
  UD$PDF <- UD$PDF * length^2
  UD$dr <- UD$dr / length
  UD$H <- UD$H / length^2
  
  return(UD)
}