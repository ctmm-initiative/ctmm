# return the most frequent convex hull containing P probability mass
convex <- function(UD,level=0.95,convex=TRUE)
{
  PMF <- raster(UD,DF="PMF")
  M <- max(UD$CDF)

  # DIM <- dim(UD$CDF)
  # xy <- array(0,c(DIM,2))
  # xy[,,1] <- UD$r$x
  # xy <- aperm(xy,c(2,1,3))
  # xy[,,2] <- UD$r$y
  # xy <- aperm(xy,c(2,1,3))
  # dim(xy) <- c(prod(DIM),2)
  # colnames(xy) <- c('x','y')

  xy <- NULL
  cost <- function(p)
  {
    if(p==0)
    { R <- 0-level }
    else if(p==1)
    { R <- M-level }
    else
    {
      # SUB <- c(UD$CDF<=p)
      # xy <- xy[SUB,] # points within p

      CL <- grDevices::contourLines(UD$r,z=UD$CDF,levels=p)
      for(cl in CL) { xy <- rbind(xy, cbind(x=cl$x,y=cl$y) ) }

      SUB <- grDevices::chull(xy) # convex hull indices
      xy <- xy[SUB,] # convex hull points

      xy <- sp::Polygon(xy) # convex hull
      xy <- sp::Polygons(list(xy),ID="ID")
      xy <- sp::SpatialPolygons(list(xy))

      PMF <- PMF * raster::rasterize(xy,PMF,background=0)
      M <- raster::cellStats(PMF,stat='sum')
      R <- M-level
    }

    return(R^2)
  }

  # this is not a continuous function
  # RESULT <- optimizer(P,cost,lower=0,upper=1,parscale=0.1)
  # p <- RESULT$par

  if(convex)
  {
    p <- c(0,level)
    RESULT <- stats::optimize(cost,p,tol = .Machine$double.eps^0.5)
    p <- RESULT$minimum
  }
  else
  { p <- level }

  # SUB <- c(UD$CDF<=p)
  # xy <- xy[SUB,] # points within p

  CL <- grDevices::contourLines(UD$r,z=UD$CDF,levels=p)
  for(cl in CL) { xy <- rbind(xy, cbind(x=cl$x,y=cl$y) ) }

  if(convex)
  {
    SUB <- grDevices::chull(xy) # convex hull indices
    xy <- xy[SUB,] # convex hull points
  }

  xy <- sp::Polygon(xy) # convex hull
  xy <- sp::Polygons(list(xy),ID="ID")
  xy <- sp::SpatialPolygons(list(xy))

  return(xy)
}
