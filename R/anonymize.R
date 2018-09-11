# unfinished: speed & heading, error ellipse
# unfinished: per individual || per list

# strip projection, origin, orientation, epoch
anonymize <- function(data)
{
  DROP <- class(data)=="telemetry"
  data <- listify(data)

  for(i in 1:length(data))
  {
    projection(data[[i]]) <- median(data[[i]],k=2)
    data[[i]]$t <- data[[i]]$t - data[[i]]$t[1]

    COLS <- c("timestamp","longitude","latitude")
    for(cl in COLS) { data[[i]][[cl]] <- NULL }

    SLOTS <- c("timezone","projection")
    for(sl in SLOTS) { attr(data[[i]],'info')[[sl]] <- NULL }
  }

  if(DROP) { data <- data[[1]] }
  return(data)
}


# give false origin, orientation, dispatch epoch
pseudonymize <- function(data,center=c(0,0),datum="WGS84",origin="0000-01-01 00:00.00 UTC",tz="GMT")
{
  DROP <- class(data)=="telemetry"
  data <- listify(data)

  proj <- paste0("+proj=aeqd +lon_0=",center[1]," +lat_0=",center[2]," +datum=",datum)

  for(i in 1:length(data))
  {
    xy <- get.telemetry(data[[i]])
    xy <- rgdal::project(xy,proj,inv=TRUE)
    data[[i]]$longitude <- xy[,1]
    data[[i]]$latitude <- xy[,2]
    attr(data[[i]],"info")$projection <- proj

    data[[i]]$timestamp <- as.POSIXct(data[[i]]$t,tz=tz,origin=origin)
    attr(data[[i]],"info")$timezone <- tz
  }

  if(DROP) { data <- data[[1]] }
  return(data)
}
