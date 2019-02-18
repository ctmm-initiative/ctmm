# append lunar phase to telemetry object
sunmoon.telemetry <- function(object,...)
{
  ADD <- suncalc::getMoonIllumination(object$timestamp,keep="fraction")$fraction
  object$moon <- ADD

  ADD <- suncalc::getSunlightPosition(date=object$timestamp,lat=object$latitude,lon=object$longitude,keep="altitude")
  object$sun <- ADD/pi + 1/2

  return(object)
}



# color by annual season, lunar phase, sunlight, or timespan
color.telemetry.list <- function(object,col="rainbow",...)
{
  n <- length(object)
  col <- col.fn(n,col=col)

  if(n<3) { return(col) }

  MED <- sapply(object,median.telemetry)
  dist <- function(m1,m2) { sqrt(sum((m1-m2)^2)) }
  D <- outer(MED,MED,dist) # distance matrix

  #

  #G <- igraph::graph_from_adjacency_matrix(D,mode="undirected",weigthed=TRUE)
  #W <- edge_attr(G,"weight")

  # can always remove a few edges and stay biconnected
  if(n>3)
  {
    #
  }

}


col.fn <- function(n,col="rainbow",alpha=1)
{
  if(col=="rainbow")
  { col <- grDevices::rainbow(n,alpha=alpha) }
  else if(col=="temperature")
  {
    col <- (0:(n-1))/max(1,n-1)
    col <- grDevices::rgb(col,0,1-col,alpha)
  }

  return(col)
}


color <- function(object,col="rainbow",alpha=1,cycle=FALSE,recursive=FALSE,...)
{
  col <- match.arg(col,c("rainbow","temperature"))

  CLASS <- class(object)
  if(CLASS=="list")
  {
    n <- length(object)

    if(recursive)
    {
      t1 <- min( sapply(object, function(o) { o$t[1] } ) )
      t2 <- max( sapply(object, function(o) { last(o$t) } ) )

      # normalize data opacity by sampling frequency
      dt <- sapply(object, function(o) { stats::median(diff(o$t)) } )
      alpha <- (dt/max(dt)) * alpha

      col <- lapply(1:n,function(i) { color.telemetry(object[[i]],col=col,alpha=alpha[i],cycle=cycle,t1=t1,t2=t2,...) })
    }
    else
    {
      if(col=="rainbow")
      { col <- grDevices::rainbow(n,alpha=alpha) }
      else if(col=="temperature")
      {
        col <- (0:(n-1))/max(1,n-1)
        col <- grDevices::rgb(col,0,1-col,alpha)
      }
    }
  }
  else if(CLASS=="telemetry")
  {
    col <- color.telemetry(object,col=col,alpha=alpha,cycle=cycle,...)
  }
  else
  { stop("Class ",CLASS," not yet supported.") }

  return(col)
}

color.telemetry <- function(object,col="rainbow",alpha=1,cycle=FALSE,t1=object$t[1],t2=last(object$t),...)
{
  t <- object$t
  n <- length(t)
  # normalize opacity by sampling frequency
  dt <- diff(t)
  dt <- dt/stats::median(dt)
  dt <- clamp(dt)
  dt <- (c(1,dt) + c(dt,1))/2
  alpha <- dt * alpha

  if(!cycle)
  {
    t <- (t-t1)/(t2-t1)

    if(col=="rainbow")
    { col <- grDevices::hsv(4/6*t,alpha=alpha) }
    else if(col=="temperature")
    { col <- grDevices::rgb(t,0,1-t,alpha) }
  }
  else
  {
    ts <- object$timestamp
    ts <- as.POSIXlt(ts)
    year <- sapply(1:n, function(i) { ts[i]$year })
    zone <- attr(object,"info")$timezone
    # bounding new years
    t1 <- as.numeric(as.POSIXct(sprintf("%d-01-01 00:00:00 %s",1900+year,zone)))
    t2 <- as.numeric(as.POSIXct(sprintf("%d-01-01 00:00:00 %s",1901+year,zone)))
    # time from/since new year in seconds
    t <- pmin(abs(t-t1),abs(t2-t))
    # normalize to 0-1 scale
    t <- t/max(t,365.24*24*60^2/2)

    if(col=="rainbow")
    { col <- grDevices::hsv(4/6*(1-t),alpha=alpha) }
    else if(col=="temperature")
    { col <- grDevices::rgb(t,0,1-t,alpha) }
  }

  return(col)
}
