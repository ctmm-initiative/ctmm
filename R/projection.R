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
