# this is stuff that needs to be run first and in order for S4 crap to work
methods::setOldClass("telemetry")
new.telemetry <- methods::setClass("telemetry", representation(info="list"), contains="data.frame")

methods::setOldClass("UD")
new.UD <- methods::setClass("UD", representation(info="list"), contains="list")

methods::setOldClass("ctmm")
new.ctmm <- methods::setClass("ctmm", representation(info="list"), contains="list")

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
methods::setGeneric("raster", getGeneric("raster", package="raster"))
methods::setGeneric("zoom", getGeneric("zoom", package="raster"))

# new S3 generic functions
writeShapefile <- function(object,...) UseMethod("writeShapefile")
