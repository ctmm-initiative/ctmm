# this is stuff that needs to be run first (and in the right order) for S4 crap to work

methods::setOldClass("UERE")
new.UERE <- methods::setClass("UERE",contains="matrix",representation=methods::representation(DOF="matrix"),
                              prototype=methods::prototype(matrix(),DOF=matrix()))

methods::setOldClass("telemetry")
new.telemetry <- methods::setClass("telemetry",contains="data.frame",representation=methods::representation(info="list",UERE="UERE"),
                                   prototype=methods::prototype(data.frame(),info=list(),UERE=new.UERE()) )

methods::setOldClass("ctmm")
new.ctmm <- methods::setClass("ctmm",contains="list",representation=methods::representation(info="list"),
                              prototype=methods::prototype(list(),info=list()))

methods::setOldClass("UD")
new.UD <- methods::setClass("UD",contains="list",representation=methods::representation(info="list",type="character",CTMM="ctmm"),
                            prototype=methods::prototype(list(),info=list(),type=character(),CTMM=new.ctmm()))

methods::setOldClass("variogram")
new.variogram <- methods::setClass("variogram",representation=methods::representation("data.frame",info="list"),
                                   prototype=methods::prototype(data.frame(),info=list()))

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

# new S3 generic functions
writeShapefile <- function(object,folder,file=NULL,...) UseMethod("writeShapefile")
emulate <- function(object,...) UseMethod("emulate")
AICc <- function(object,...) UseMethod("AICc")
speed <- function(object,...) UseMethod("speed")
