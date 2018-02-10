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


"projection<-.telemetry" <- function(x,value)
{
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
  R.EQ <- 6378137
  R.PL <- 6356752.3142
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


########################################
# Suggest a good projection
########################################
suggest.projection <- function(data,PROJ="tpeqd",datum="WGS84")
{
  # projections that we will populate/automate
  PROJS <- c("aeqd","tpeqd")
  if(!(PROJ %in% PROJS)) { return(PROJ) }

  # assume Movebank data.frame
  lon <- data$longitude
  lat <- data$latitude

  # as a first approximation use one-point equidistant at average geolocation
  lon_0 <- stats::median(lon)
  lat_0 <- stats::median(lat)
  proj <- paste0("+proj=aeqd +lon_0=",lon_0," +lat_0=",lat_0," +datum=",datum)

  if(PROJ=="aeqd") { return(proj) }

  xy <- rgdal::project(cbind(lon,lat),proj)

  # calculate and detrend average
  mu <- c(stats::median(xy[,1]),stats::median(xy[,2]))
  xy <- xy - mu
  colnames(xy) <- c("x","y")

  # cross correlation
  cov <- mean(xy[,1]*xy[,2])
  # covariance matrix
  cov <- rbind( c( mean(xy[,1]^2) , cov ) , c( cov , mean(xy[,2]^2) ) )

  # figure out long axis (always first dim)
  R <- eigen(cov)$vectors
  # rotate data to long axis
  xy <- xy %*% R

  # bi-modal split of data
  xy1 <- xy[xy[,1]<0,]
  xy2 <- xy[xy[,1]>0,]

  # bi-modal modes
  mu1 <- c(stats::median(xy1[,1]),stats::median(xy1[,2]))
  mu2 <- c(stats::median(xy2[,1]),stats::median(xy2[,2]))

  # reverse rotation
  R <- solve(R)
  mu1 <- mu1 %*% R
  mu2 <- mu2 %*% R

  # re-trend mean
  mu1 <- mu1 + mu
  mu2 <- mu2 + mu

  # get long lat
  mu1 <- rgdal::project(mu1,proj,inv=TRUE)[1,]
  mu2 <- rgdal::project(mu2,proj,inv=TRUE)[1,]

  # did east and west get mixed up?
  if(mu1[1] > mu2[1])
  {
    mu <- mu1
    mu1 <- mu2
    mu2 <- mu
  }

  proj <- paste0("+proj=tpeqd +lon_1=",mu1[1]," +lat_1=",mu1[2]," +lon_2=",mu2[1]," +lat_2=",mu2[2]," +datum=",datum)

  return(proj)
  #################
  #STOP HERE

  # NON-FUNCTIONAL NORTH ROTATION CODE
  # This doesn't seem to do anything. I don't think the +axis and +towgs84 options are fully implemented in PROJ4.
  # keeping this code here for later.

  # project origin back
  mu <- rgdal::project(rbind(c(0,0)),proj,inv=TRUE)[1,]

  # add a touch of north
  mu <- mu + rbind(c(0,0.001))

  # project forward
  mu <- rgdal::project(mu,proj)[1,]

  # solve for rotation angle to get north vector
  theta <- atan2(mu[2],mu[1])

  # generate rotated projection...
  # abusing PROj4 small-angle rotation + rescaling = true rotation
  # this would ruin z data if we had any
  proj <- paste(proj," +towgs84=0,0,0,0,0,",tan(theta),",",cos(theta),sep="")
}


validate.projection <- function(projection)
{
  if(grepl("latlong",projection,fixed=TRUE) || grepl("latlong",projection,fixed=TRUE))
  { stop("A projected coordinate system must be specified.") }

  if(grepl("units=",projection,fixed=TRUE) && !grepl("units=m",projection,fixed=TRUE))
  { stop("Units of distance other than meters not supported.") }
}
