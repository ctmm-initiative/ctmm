################################
# create a raster of the ML akde
raster.UD <- function(x,DF="CDF",...)
{
  proj <- attr(x,"info")$projection
  UD <- x

  dx <- UD$dr[1]
  dy <- UD$dr[2]

  xmn <- UD$r$x[1]-dx/2
  xmx <- last(UD$r$x)+dx/2

  ymn <- UD$r$y[1]-dy/2
  ymx <- last(UD$r$y)+dy/2

  z <- UD$r$z

  # probability mass for the cells
  if(DF=="PMF")
  { UD <- UD[["PDF"]] * prod(UD$dr) }
  else
  { UD <- UD[[DF]] }

  if(length(dim(UD))==2)
  {
    UD <- t(UD[,dim(UD)[2]:1])
    R <- raster::raster(UD,xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx,crs=proj)
  }
  else
  {
    UD <- aperm(UD[,dim(UD)[2]:1,],c(2,1,3))
    R <- raster::brick(UD,xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx,crs=proj)
    R <- raster::setZ(R,z,name="height")
  }

  return(R)
}
methods::setMethod("raster",methods::signature(x="UD"), function(x,DF="CDF",...) raster.UD(x,DF=DF,...))


##########################
writeRaster.UD <- function(x,filename,format,DF="CDF",...)
{
  if(missing(filename)) { filename <- attr(x,"info")$identity }

  x <- raster(x,DF=DF)
  writeRaster(x,filename,format,...)
}
methods::setMethod("writeRaster",methods::signature(x="UD",filename="character"), function(x,filename,format,DF="CDF",...) writeRaster.UD(x,filename,format,DF=DF,...))


################
# Is contour A inside contour B
inside <- function(A,B)
{
  result <- mode(sp::point.in.polygon(A$x,A$y,B$x,B$y))
  if(1<=result && result<=2) { return(1) } else { return(0) }
}


##############
SpatialPolygonsDataFrame.UD <- function(object,level.UD=0.95,level=0.95,convex=FALSE,...)
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

  seg.id <- function(id) { ifelse(id<=length(CL),id,id-length(CL)) }
  seg.rv <- function(id) { id <- seg.id(id); CL[[id]]$x <<- rev(CL[[id]]$x) ; CL[[id]]$y <<- rev(CL[[id]]$y) }

  polygons <- list()
  for(i in 1:length(P))
  {
    if(!convex)
    {
      options(max.contour.segments=.Machine$integer.max)
      CL <- contourLines(UD,levels=P[i])

      CHANGE <- TRUE
      while(CHANGE && length(CL)>1) # segment join loop
      {
        ## join contour segments ## contourLines can fail to join segments into complete polygons
        DIST <- array(0,c(2*length(CL),2))
        DIST[1:length(CL),] <- t( sapply(CL,function(cl){c(cl$x[1],cl$y[1])}) )
        DIST[length(CL)+1:length(CL),] <- t( sapply(CL,function(cl){c(last(cl$x),last(cl$y))}) )
        DIST <- outer(DIST[,1],FUN='-')^2 + outer(DIST[,2],FUN='-')^2
        diag(DIST) <- Inf

        MINS <- sort(DIST,index.return=TRUE)$ix
        MINS <- arrayInd(MINS,dim(DIST))

        MATCH <- array(0,2*length(CL))
        for(j in nrow(MINS):1)
        {
          MATCH[MINS[j,1]] <- MINS[j,2]
          MATCH[MINS[j,2]] <- MINS[j,1]
        }

        CHANGE <- FALSE
        for(j in 1:nrow(MATCH))
        {
          k <- MATCH[j]
          if(seg.id(j)!=seg.id(k))
          {
            # connect to the tail of i
            if(j<=seg.id(j)) { seg.rv(j) }
            j <- seg.id(j)

            # connect to the beginning of j
            if(k>seg.id(k)) { seg.rv(k) }
            k <- seg.id(k)

            CL[[j]]$x <- c( CL[[j]]$x , CL[[k]]$x )
            CL[[j]]$y <- c( CL[[j]]$y , CL[[k]]$y )
            CL <- CL[-k]
            CHANGE <- TRUE
            break
          } # made change and re-calculated distances
        } # no changes made, consider next match
      } # no more disconnected segments

      # create contour heirarchy matrix (half of it)
      H <- array(0,c(1,1)*length(CL))
      for(row in 1:length(CL))
      {
        for(col in row:length(CL))
        { H[row,col] <- inside(CL[[row]],CL[[col]]) }
      }

      # number of contours that this contour is inside
      I <- rowSums(H)

      # if I is odd, then you are a hole inside a positive area
      hole <- is.odd(I)

      # polygon
      polygons[[i]] <- list()
      for(j in 1:length(CL))
      { polygons[[i]][[j]] <- sp::Polygon( cbind( CL[[j]]$x , CL[[j]]$y ) , hole=hole[j] ) }

      # polygonS
      polygons[[i]] <- sp::Polygons(polygons[[i]],ID=ID[i])
    }
    else # if(convex)
    {
      CL <- convex(UD,level=P[i],SP=FALSE,ID=ID[i],convex=convex)

      polygons[[i]] <- CL
    }
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
#methods::setMethod("SpatialPolygonsDataFrame",methods::signature(Sr="UD"), function(Sr,level.UD=0.95,level=0.95) SpatialPolygonsDataFrame.UD(Sr,level.UD=level.UD,level=level))


###########
writeVector.UD <- function(x,filename,filetype="ESRI Shapefile",level.UD=0.95,level=0.95,convex=FALSE,...)
{
  if(missing(filename)) { filename <- attr(x,"info")$identity }
  x <- SpatialPolygonsDataFrame.UD(x,convex=convex,level.UD=level.UD,level=level)
  x <- terra::vect(x)
  terra::writeVector(x,filename,filetype=filetype,...)
}
methods::setMethod("writeVector",methods::signature(x="UD",filename="character"), function(x,filename,filetype="ESRI Shapefile",level.UD=0.95,level=0.95,convex=FALSE,...) writeVector.UD(x,filename,filetype=filetype,convex=convex,level.UD=level.UD,level=level,...) )
methods::setMethod("writeVector",methods::signature(x="UD",filename="missing"), function(x,filename,filetype="ESRI Shapefile",level.UD=0.95,level=0.95,convex=FALSE,...) writeVector.UD(x,filename,filetype=filetype,convex=convex,level.UD=level.UD,level=level,...) )


################
writeVector.telemetry <- function(x,filename,filetype="ESRI Shapefile",error=TRUE,level.UD=0.95,...)
{
  if(missing(filename)) { filename <- mean_info(x)$identity }

  # make one long SPDF
  if(error) { x <- SpatialPolygonsDataFrame.telemetry(x,level.UD=level.UD) }
  else { x <- SpatialPointsDataFrame.telemetry(x) }

  # make shape file from points
  x <- terra::vect(x)
  terra::writeVector(x,filename,filetype=filetype,...)

  # no return value
}
methods::setMethod("writeVector",methods::signature(x="telemetry",filename="character"), function(x,filename,filetype="ESRI Shapefile",error=TRUE,level.UD=0.95,...) writeVector.telemetry(x,filename,filetype=filetype,error=error,level.UD=level.UD,...) )
methods::setMethod("writeVector",methods::signature(x="telemetry",filename="missing"), function(x,filename,filetype="ESRI Shapefile",error=TRUE,level.UD=0.95,...) writeVector.telemetry(x,filename,filetype=filetype,error=error,level.UD=level.UD,...) )


#########################
# convert to spatialpoints object
SpatialPoints.telemetry <- function(object,...)
{
  CLASS <- class(object)[1]
  # promote to list if not already
  object <- listify(object)

  SP <- lapply( object, function(d) { sp::SpatialPoints( "[.data.frame"(d,c("x","y")), proj4string=sp::CRS(attr(d,"info")$projection) ) } )

  # rbind all together
  SP <- do.call(sp::rbind.SpatialPoints,SP)

  return(SP)
}
#methods::setMethod("SpatialPoints",methods::signature(coords="telemetry"), function(coords) SpatialPoints.telemetry(coords))


##############
# convert to spatialpoints data.frame object
SpatialPointsDataFrame.telemetry <- function(object,...)
{
  # promote to list if not already
  object <- listify(object)

  SP <- SpatialPoints.telemetry(object,...)

  # make identity array
  identity <- unlist(sapply(object,function(o){ rep(attr(o,"info")$identity,length(o$t)) }))
  # sp does something weird to POSIXct columns -> character
  timestamp <- do.call(c,lapply(object,function(o){ paste(o$timestamp,attr(o,"info")$timezone) }))

  DF <- data.frame(identity=identity,timestamp=timestamp,row.names=1:length(identity))
  names(DF) <- c("identity","timestamp") # weird namespace collision with identity()

  # make SPDF with identity information
  SP <- sp::SpatialPointsDataFrame(SP,DF,match.ID=FALSE)

  return(SP)
}


########
SpatialPolygonsDataFrame.telemetry <- function(object,level.UD=0.95,...)
{
  object <- listify(object)
  identity <- unlist(sapply(object,function(o){ rep(attr(o,"info")$identity,length(o$t)) }))
  if(!length(identity)) { identity <- 1:length(object) }
  # sp does something weird to POSIXct columns -> character
  timestamp <- do.call(c,lapply(object,function(o){ paste(o$timestamp,attr(o,"info")$timezone) }))

  proj <- projection(object)
  if(length(proj)!=1) { stop("Inconsistent projections.") }

  polygons <- list()
  ID <- 1
  for(i in 1:length(object))
  {
    R <- get.telemetry(object[[i]])
    S <- get.error(object[[i]],ctmm(error=10),DIM=2) # UERE of 10 default

    polygons[[i]] <- list()
    for(j in 1:nrow(object[[i]]))
    {
      P <- ellipsograph(mu=R[j,],sigma=S[j,,],level=level.UD,PLOT=FALSE)
      P <- rbind(P,P[1,]) # first=last
      P <- sp::Polygon(P,hole=FALSE)
      P <- sp::Polygons(list(P),ID=ID)
      ID <- ID + 1 # annoying sp requirement
      polygons[[i]][[j]] <- P
    }
  }
  polygons <- do.call(c,polygons)
  polygons <- sp::SpatialPolygons(polygons, proj4string=sp::CRS(proj))
  # names(polygons) <- 1:length(polygons)

  n <- length(polygons)
  if(n>length(identity))
  {
    identity <- rep(identity,n)
    timestamp <- rep(timestamp,n)
  }
  # spatial polygons data frame
  DF <- data.frame(identity=identity,timestamp=timestamp,row.names=names(polygons))
  names(DF) <- c("identity","timestamp") # weird namespace collision with identity()
  polygons <- sp::SpatialPolygonsDataFrame(polygons,DF)

  return(polygons)
}


# Export: ctmm -> sp -> sf
as.sf <- function(x,error=FALSE,...)
{
  CLASS <- class(x)[1]

  if(CLASS=="UD")
  { x <- SpatialPolygonsDataFrame.UD(x,...) }
  else if(CLASS=="telemetry")
  {
    if(error)
    { x <- SpatialPolygonsDataFrame.telemetry(x,...) }
    else
    { x <- SpatialPointsDataFrame.telemetry(x,...) }
  }

  sf::st_as_sf(x)
}
