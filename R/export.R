################################
# create a raster of the ML akde
raster.UD <- function(x,DF="CDF",...)
{
  UD <- x

  dx <- UD$dr[1]
  dy <- UD$dr[2]

  xmn <- UD$r$x[1]-dx/2
  xmx <- last(UD$r$x)+dx/2

  ymn <- UD$r$y[1]-dy/2
  ymx <- last(UD$r$y)+dy/2

  # probability mass for the cells
  if(DF=="PMF")
  {
    DF <- "PDF"
    UD[[DF]] <- UD[[DF]] * prod(UD$dr)
  }

  Raster <- raster::raster(t(UD[[DF]][,dim(UD[[DF]])[2]:1]),xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx,crs=attr(UD,"info")$projection)

  return(Raster)
}
methods::setMethod("raster",signature(x="UD"), function(x,DF="CDF",...) raster.UD(x,DF=DF,...))


##########################
writeRaster.UD <- function(x,filename,format,DF="CDF",...)
{
  x <- raster(x,DF="CDF") # ... here seems much less useful
  writeRaster(x,filename,format,...)
}
methods::setMethod("writeRaster",signature(x="UD",filename="character"), function(x,filename,format,DF="CDF",...) writeRaster.UD(x,filename,format,DF=DF,...))


################
# Is contour A inside contour B
inside <- function(A,B)
{
  result <- mode(sp::point.in.polygon(A$x,A$y,B$x,B$y))
  if(1<=result && result<=2) { return(1) } else { return(0) }
}


##############
SpatialPolygonsDataFrame.UD <- function(object,level.UD=0.95,level=0.95,...)
{
  UD <- object

  # populate arrays: level.UD versus level
  P <- NULL
  ID <- NULL
  for(i in 1:length(level.UD))
  {
    p <- CI.UD(UD,level.UD[i],level,P=TRUE)
    P <- cbind(P,p)
    ID <- cbind(ID,paste(UD@info$identity," ",round(100*level.UD[i]),"% ",names(p),sep=""))
  }

  # flatten arrays
  P <- c(P)
  ID <- c(ID)

  polygons <- list()
  for(i in 1:length(P))
  {
    CL <- grDevices::contourLines(UD$r,z=UD$CDF,levels=P[i])

    if(length(CL)==0) # nuge sp to make a contour
    { CL <- grDevices::contourLines(UD$r,z=UD$CDF,levels=P[i]*(1+.Machine$double.eps)) }

    # create contour heirarchy matrix (half of it)
    H <- array(0,c(1,1)*length(CL))
    for(row in 1:length(CL))
    {
      for(col in row:length(CL))
      {
        H[row,col] <- inside(CL[[row]],CL[[col]])
      }
    }

    # number of contours that this contour is inside
    I <- rowSums(H)

    # if I is odd, then you are a hole inside a positive area
    hole <- is.odd(I)

    # polygon
    polygons[[i]] <- list()
    for(j in 1:length(CL))
    {
      polygons[[i]][[j]] <- sp::Polygon( cbind( CL[[j]]$x , CL[[j]]$y ) , hole=hole[j] )
    }

    # polygonS
    polygons[[i]] <- sp::Polygons(polygons[[i]],ID=ID[i])
  }
  names(polygons) <- ID

  # spatial polygons
  polygons <- sp::SpatialPolygons(polygons, proj4string=sp::CRS(attr(UD,"info")$projection))

  # spatial polygons data frame
  data <- data.frame(name=rev(ID))
  rownames(data) <- rev(ID)
  polygons <- sp::SpatialPolygonsDataFrame(polygons,data)

  return(polygons)
}
#methods::setMethod("SpatialPolygonsDataFrame",signature(Sr="UD"), function(Sr,level.UD=0.95,level=0.95) SpatialPolygonsDataFrame.UD(Sr,level.UD=level.UD,level=level))


################
writeShapefile.UD <- function(object,folder,file=NULL,level.UD=0.95,level=0.95,...)
{
  UD <- object
  if(is.null(file)) { file <- attr(object,"info")$identity }

  SP <- SpatialPolygonsDataFrame.UD(UD,level.UD=level.UD,level=level)

  rgdal::writeOGR(SP, dsn=folder, layer=file, driver="ESRI Shapefile",...)
}


#########################
# convert to spatialpoints object
SpatialPoints.telemetry <- function(object,...)
{
  CLASS <- class(object)
  # promote to list if not already
  object <- listify(object)

  SP <- lapply( object, function(d) { sp::SpatialPoints( "[.data.frame"(d,c("x","y")), proj4string=sp::CRS(attr(d,"info")$projection) ) } )

  # rbind all together
  SP <- do.call(sp::rbind.SpatialPoints,SP)

  return(SP)
}
#methods::setMethod("SpatialPoints",signature(coords="telemetry"), function(coords) SpatialPoints.telemetry(coords))


##############
# convert to spatialpoints data.frame object
SpatialPointsDataFrame.telemetry <- function(object,...)
{
  # promote to list if not already
  object <- listify(object)

  # make identity array
  identity <- unlist(sapply(object,function(o){ rep(attr(o,"info")$identity,length(o$t)) }))

  SP <- SpatialPoints.telemetry(object,...)

  # make SPDF with identity information
  SP <- sp::SpatialPointsDataFrame(SP,data.frame(identity=identity),match.ID=FALSE)

  return(SP)
}

################
writeShapefile.telemetry <- function(object,folder,file=NULL, ...)
{
  if(is.null(file)) { file <- mean.info(object) }

  # make one long SPDF
  SP <- SpatialPointsDataFrame.telemetry(object)

  # make shape file from points
  rgdal::writeOGR(SP, dsn=folder, layer=file, driver="ESRI Shapefile",...)

  # no return value
}
