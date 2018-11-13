# global variable
DATA.EARTH <- list(R.EQ=6378137,R.PL=6356752.3142) # equatorial & polar radii


# range of telemetry data
projection.telemetry <- function(x,asText=TRUE)
{
  proj <- attr(x,"info")$projection
  if(!asText) { proj <- sp::CRS(proj,doCheckCRSArgs=FALSE) }
  return(proj)
}
setMethod('projection', signature(x='telemetry'), projection.telemetry)
setMethod('projection', signature(x='ctmm'), projection.telemetry)
setMethod('projection', signature(x='UD'), projection.telemetry)


projection.list <- function(x,asText=TRUE)
{
  PROJS <- sapply(x,projection)
  PROJ <- unique(PROJS)

  if(is.null(PROJS[[1]])) { return(NULL) }
  if(length(PROJ)==1) { return(PROJ) }
  return(PROJS)
}
setMethod('projection',signature(x='list'),projection.list)


# change the projection on a list of objects
"projection<-.list" <- function(x,value)
{
  lapply(x,function(y){methods::getMethod("projection<-",signature=class(y))(y,value)})
}
setMethod('projection<-', signature(x='list'), `projection<-.list`)


# change the projection of one telemetry object
"projection<-.telemetry" <- function(x,value)
{
  # delete projection information
  if(is.null(value))
  {
    DELETE <- c(DOP.LIST$horizontal$axes,DOP.LIST$horizontal$VAR,DOP.LIST$horizontal$COV,DOP.LIST$speed$axes,DOP.LIST$speed$VAR,DOP.LIST$speed$COV)
    for(NAME in DELETE) { x[[NAME]] <- NULL }

    attr(x,'info')$projection <- NULL
    return(x)
  }

  # convert to PROJ4 format if location
  value <- format.projection(value)

  NAMES <- names(x)
  n <- nrow(x)

  ### project locations ###
  R <- cbind(x$longitude,x$latitude)
  colnames(R) <- c("x","y")
  R <- rgdal::project(R,value)
  x[c('x','y')] <- R
  rm(R)

  if(any(c('speed','COV.angle') %in% NAMES))
  {
    # local vectors pointing north
    NORTH <- northing(x,value)

    ### project error ellipses ###
    if('COV.angle' %in%  NAMES)
    {
      # major axis eigen vector
      COV <- rotate.north(NORTH,x$COV.angle) # [n,2]
      # major axis eigen matrices
      COV <- sapply(1:n,function(i){outer(COV[i,])},simplify="array") # [2,2,n]
      ID <- diag(2)
      # covariance matrices
      COV <- sapply(1:n,function(i){x$COV.major[i]*COV[,,i] + x$COV.minor[i]*(ID-COV[,,i])},simplify="array") # [2,2,n]
      # flatten spatial indices
      dim(COV) <- c(2*2,n)
      # unique entries
      ID <- c(upper.tri(ID,diag=TRUE))
      x[DOP.LIST$horizontal$COV] <- t(COV[ID,]) # R is so confusing

      rm(COV)
    }

    ### project velocities ###
    if('speed' %in% NAMES)
    {
      # set magnitude and heading of velocity
      v <- rotate.north(x$speed*NORTH,x$heading)
      # velocity components
      x[DOP.LIST$speed$axes] <- v

      rm(v)
    }
  }

  attr(x,"info")$projection <- value

  return(x)
}
setMethod('projection<-', signature(x='telemetry'), `projection<-.telemetry`)


# unit vector pointing north at projected locations R
# x = partially projected telemetry data dim(n,2)
# return dim(n,2)
northing <- function(x,proj)
{
  # WGS-84 ellipsoid
  R.EQ <- DATA.EARTH$R.EQ
  R.PL <- DATA.EARTH$R.PL
  # approximate 1-meter-North latitude displacements
  d.lambda <- 1/sqrt((R.EQ*sin(x$latitude))^2+(R.PL*cos(x$latitude))^2)
  # could use grad() but would be slowwwww....
  u <- cbind(x$longitude,x$latitude + d.lambda*(360/2/pi)) # arg, degrees!
  colnames(u) <- c("x","y")
  u <- rgdal::project(u,proj)
  # difference vectors pointing North ~1 meters
  u <- u - get.telemetry(x) # [n,2]
  # difference vectors pointing North 1 meters exact
  u <- u / sqrt(rowSums(u^2)) # [n,2]
  return(u)
}

# rotate northing to heading
# u = dim(2,n)
# return = dim(n,2)
rotate.north <- function(u,heading)
{
  heading <- heading * (2*pi/360) # ack, degrees!

  u <- u[,1] + 1i*u[,2] # velocity vector
  R <- exp(-1i*heading) # rotation matrix
  u <- R*u
  u <- cbind(Re(u),Im(u))

  return(u)
}


# put projection into character format
format.projection <- function(proj,datum="WGS84")
{
  if(class(proj)=="CRS")
  { proj <- as.character(proj) }
  else if(class(proj)!="character")
  {
    # pull out geodesic coordinates and format into matrix
    proj <- as.matrix(rbind(proj)[,c("longitude","latitude")])

    if(nrow(proj)==1)
    { proj <- paste0("+proj=aeqd  +lon_0=",proj[1,1]," +lat_0=",proj[1,2]," +datum=",datum) }
    else if(nrow(proj)==2)
    { proj <- paste0("+proj=tpeqd +lon_1=",proj[1,1]," +lat_1=",proj[1,2]," +lon_2=",proj[2,1]," +lat_2=",proj[2,2]," +datum=",datum) }
    else if(nrow(proj)==3)
    { proj <- paste0("+proj=chamb +lon_1=",proj[1,1]," +lat_1=",proj[1,2]," +lon_2=",proj[2,1]," +lat_2=",proj[2,2]," +lon_3=",proj[3,1]," +lat_3=",proj[3,2]," +datum=",datum) }
    else
    { stop("PROJ4 does not support ",nrow(proj)," foci projections.") }
  }

  validate.projection(proj)
  return(proj)
}


# only allow compatible projections
validate.projection <- function(projection)
{
  if(class(projection)=="character") { projection <- sp::CRS(projection) } # this adds missing longlat specification
  if(class(projection)=="CRS") { projection <- as.character(projection) }

  if(grepl("longlat",projection,fixed=TRUE) || grepl("latlong",projection,fixed=TRUE))
  { stop("A projected coordinate system must be specified.") }

  if(grepl("units=",projection,fixed=TRUE) && !grepl("units=m",projection,fixed=TRUE))
  { stop("Units of distance other than meters not supported.") }
}
