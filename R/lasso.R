# polygon lasso selector
lasso <- function(object,...) { selector(object,tool="lasso",...) }

# rectangular marquee selector
marquee <- function(object,...) { selector(object,tool="marquee",...) }

# shared code
selector <- function(object,tool,...)
{
  object <- listify(object)
  plot(object,col='blue',...)
  SCALE <- get('x.scale',pos=plot.env) # get units of base plot

  if(tool=='lasso')
  {
    if(.Platform$GUI=="RStudio") { message("Click to place vertices and press [ESC] when done.") }
    else { message("Left-click to place vertices and right-click for menu when done.") }
    SP <- raster::drawPoly(sp=FALSE,col='black')
  }
  else if(tool=='marquee')
  {
    message("Click opposite corners.")
    SP <- raster::drawExtent(col='black')
    SP <- rbind( c(SP@xmin,SP@ymin), c(SP@xmax,SP@ymin), c(SP@xmax,SP@ymax), c(SP@xmin,SP@ymax) )
  }
  SP <- SCALE * SP # convert to SI units
  SP <- lapply(object,function(o){ as.logical( sp::point.in.polygon(o$x,o$y,SP[,1],SP[,2]) ) })

  IN <- list()
  OUT <- list()
  for(i in 1:length(object))
  {
    IN[[i]] <- object[[i]][SP[[i]],]
    OUT[[i]] <- object[[i]][!SP[[i]],]
  }
  names(IN) <- paste("interior",names(object))
  names(OUT) <- paste("exterior",names(object))
  rm(object)

  ## replace with generic points() function ##
  # set units of OUT (for plot only)
  # highlight selected points
  #PIN <- lapply(IN,function(O){ unit.telemetry(O,length=SCALE) })
  plot(IN,col='red',add=TRUE,...) #!!!
  #rm(PIN)

  # return everything
  return( c(IN,OUT) )
}
