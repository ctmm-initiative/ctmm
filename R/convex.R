# safe contour line function
contourLines <- function(UD=NULL,r=UD$r,CDF=UD$CDF,levels=0.95)
{
  CL <- grDevices::contourLines(r,z=CDF,levels=levels)

  if(length(CL)==0) # nudge sp to make a contour
  { CL <- grDevices::contourLines(r,z=CDF,levels=levels*(1+.Machine$double.eps)) }

  if(length(CL)==0) # fallback on the mode
  {
    CL <- list()
    CL$level <- levels
    MODE <- which.min(CDF)
    MODE <- arrayInd(MODE,dim(CDF))
    MODE <- c( r$x[ MODE[1,1] ] , r$y[ MODE[1,2] ] )
    AREA <- summary(UD,level.UD=levels,units=FALSE)$CI[2]
    RAD <- sqrt(AREA/pi)
    RAD <- max(RAD,.Machine$double.eps)
    CL$x <- MODE[1] + c(1,0,-1,0)*RAD
    CL$y <- MODE[2] + c(0,1,0,-1)*RAD
    CL <- list(CL)
  }

  # # sp complains when less than 4 vertices
  # GOOD <- rep(TRUE,length(CL))
  # for(i in 1:length(CL))
  # { if(length(CL[[i]]$x)<4) { GOOD[i] <- FALSE } }
  # CL <- CL[GOOD]
  # # doing this, sp then gives an error about the area slot missing
  # # TODO check if sf also complains about this

  return(CL)
}

# return the most frequent convex hull containing P probability mass
convex <- function(UD,level=0.95,convex=TRUE,SP=TRUE,ID="ID")
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

      CL <- contourLines(UD,levels=p)
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

  CL <- contourLines(UD,levels=p)
  for(cl in CL) { xy <- rbind(xy, cbind(x=cl$x,y=cl$y) ) }

  if(convex)
  {
    SUB <- grDevices::chull(xy) # convex hull indices
    xy <- xy[SUB,] # convex hull points
  }

  xy <- sp::Polygon(xy) # convex hull
  xy <- sp::Polygons(list(xy),ID=ID)

  if(SP)
  { xy <- sp::SpatialPolygons(list(xy)) }

  return(xy)
}
