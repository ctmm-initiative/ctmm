# circular quantile functions for longitudes
# assuming no global outliers in the data

median.longitude <- function(x,na.rm=FALSE,...)
{ quantile.longitude(x,probs=0.5,na.rm=na.rm,...) }

quantile.longitude <- function(x,probs=c(0,1),na.rm=FALSE,...)
{
  if(na.rm) { x <- x[!is.na(x)] }
  n <- length(x)
  angle <- sort(x)
  angle <- c(angle,angle[1]) # wrap around
  delta <- diff(angle) %% 360 # gap widths
  MAX <- which.max(delta) # gap complement of the data # assuming no global outliers
  x <- angle[1:MAX]
  if(MAX<n) { x <- c(angle[(MAX+1):n],x) } # x is now compacted into a single interval (not yet manifest)
  MIN <- x[1]
  x <- (x - MIN) %% 360 # [0,360) interval
  x <- stats::quantile(x,probs=probs,...)
  x <- x + MIN
  x <- ((x + 180) %% 360) - 180 # (-180,180]
  return(x)
}
