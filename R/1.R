# universal names for confidence intervals & point estimate
NAMES.CI <- c("low","est","high")

# this is stuff that needs to be run first (and in the right order) for S4 crap to work

methods::setOldClass("UERE")
new.UERE <- methods::setClass("UERE",contains="list",representation=methods::representation(info="list"),
                              prototype=methods::prototype(list(UERE=cbind(0),DOF=cbind(0),AICc=Inf,Zsq=Inf,VAR.Zsq=Inf,N=0),info=list()))
#DOF="matrix",AICc="numeric",Zsq="numeric",VAR.Zsq="numeric",N="numeric"

methods::setOldClass("telemetry")
new.telemetry <- methods::setClass("telemetry",contains="data.frame",representation=methods::representation(info="list",UERE="UERE"),
                                   prototype=methods::prototype(data.frame(),info=list(),UERE=new.UERE()) )

methods::setOldClass("ctmm")
new.ctmm <- methods::setClass("ctmm",contains="list",representation=methods::representation(info="list"),
                              prototype=methods::prototype(list(),info=list()))

methods::setOldClass("UD")
new.UD <- methods::setClass("UD",contains="list",representation=methods::representation(info="list",type="character",variable="character",CTMM="ctmm"),
                            prototype=methods::prototype(list(),info=list(),type=character(),variable=character(),CTMM=new.ctmm()) )

methods::setOldClass("variogram")
new.variogram <- methods::setClass("variogram",representation=methods::representation("data.frame",info="list",UERE="UERE"),
                                   prototype=methods::prototype(data.frame(),info=list(),UERE=new.UERE()) )

methods::setOldClass("outlie")
new.outlie <- methods::setClass("outlie",representation=methods::representation("data.frame"),prototype=methods::prototype(data.frame()))

# R drop is very annoying and yet this doesn't do anything despite the error on det(1)
#setMethod('determinant', signature(x='numeric'), identity)

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

# existing S4 generic functions
methods::setGeneric("projection", getGeneric("projection", package="raster"))
methods::setGeneric("projection<-", getGeneric("projection<-", package="raster"))
methods::setGeneric("raster", getGeneric("raster", package="raster"))
methods::setGeneric("zoom", getGeneric("zoom", package="raster"))
methods::setGeneric("writeVector", getGeneric("writeVector", package="terra"))


# new S3 generic functions
emulate <- function(object,...) UseMethod("emulate")
AICc <- function(object,...) UseMethod("AICc")
speed <- function(object,...) UseMethod("speed")
speeds <- function(object,...) UseMethod("speeds")
mag <- function(x,...) UseMethod("mag")
#modes <- function(object,...) UseMethod("modes") # CRAN is annoying: do not declare generics for non-exported functions
ridges <- function(object,...) UseMethod("ridges")

# internal S3 generic function
pars <- function(...) { UseMethod("pars") }

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

IFFT <- function(X,plan=NULL) { FFT(X,inverse=TRUE) }

# is a package installed?
is.installed <- function(pkg) is.element(pkg, utils::installed.packages()[,1])

.onLoad <- function(...)
{
  # new global options
  if(is.null(getOption("time.units"))) { options(time.units='mean') }
  utils::assignInMyNamespace("UNIT", generate.units())

  # choose FFTW if installed
  if(is.installed("fftw")) { utils::assignInMyNamespace("FFT", FFTW) }
}
.onAttach <- .onLoad
