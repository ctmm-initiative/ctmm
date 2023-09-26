######
is.grid.complete <- function(grid)
{
  if(is.null(grid) || class(grid)[1]=="Extent") { return(FALSE) }
  if(class(grid)[1] %in% c("UD","RasterLayer")) { return(TRUE) }
  if(class(grid)[1]=="list") # list of options
  { return( "r" %in% names(grid) || all(c("extent","dr") %in% names(grid)) ) }
  else # assuming extent matrix or dr numeric
  { return(FALSE) }
}


is.consistent <- function(grids)
{
  TEST <- logical(length(grids)-1)
  for(i in 1:(length(TEST)))
  { TEST[i] <- is.consistent.pair(grids[[i]],grids[[i+1]]) }
  TEST <- all(TEST)
  return(TEST)
}


is.consistent.pair <- function(grid1,grid2)
{
  dr1 <- grid1$dr
  dr2 <- grid2$dr

  MAX <- pmax(dr1,dr2)
  MIN <- pmin(dr1,dr2)

  TEST <- log2(MAX/MIN) # should be an integer
  TEST <- abs( TEST - round(TEST) )
  TEST <- any( TEST > .Machine$double.eps )
  if(TEST) { return(FALSE) }

  r1 <- c(grid1$r$x[1],grid1$r$y[1])
  r2 <- c(grid2$r$x[1],grid2$r$y[1])

  TEST <- (r1-r2)/MIN
  TEST <- TEST - round(TEST)
  TEST <- any( TEST > .Machine$double.eps )
  return(!TEST)
}


grid.comp <- function(grids)
{
  DR <- sapply(grids,function(g){g$dr}) # [2,n]
  dr <- apply(DR,1,min) # smallest (dx,dy)
  DR <- t(t(DR)/dr)
}


###########
# returns grid information for union
# assumes origin aligned grids with same resolution
grid.union <- function(UD)
{
  n <- length(UD)

  r <- list()
  for(i in 1:n) { r[[i]] <- UD[[i]]$r }

  dr <- sapply(UD,function(ud){ud$dr}) # (axes,individuals)
  if(any(diff(dr[1,])!=0 | diff(dr[2,])!=0)) { stop("Inconsistent grid resolutions.") }
  dr <- UD[[1]]$dr

  # are all grids the same?
  TEST <- TRUE
  for(i in 1%:%(n-1)) { TEST <- TEST && identical(r[[i]],r[[i+1]]) }

  if(TEST)
  { R <- r[[1]] }
  else
  {
    x.min <- min( sapply(r,function(R){first(R$x)}) )
    x.max <- max( sapply(r,function(R){last(R$x)}) )

    y.min <- min( sapply(r,function(R){first(R$y)}) )
    y.max <- max( sapply(r,function(R){last(R$y)}) )

    R <- list()
    if(x.min<=x.max) { R$x <- seq(x.min,x.max,length.out=1+round((x.max-x.min)/dr['x'])) } # by fails with numerical error
    if(y.min<=y.max) { R$y <- seq(y.min,y.max,length.out=1+round((y.max-y.min)/dr['y'])) } # by fails with numerical error
  }

  return(list(r=R,dr=dr))
}


###########
# returns list of subset indices for intersection
# assumes origin aligned grids with same resolution
grid.intersection <- function(UD)
{
  tol <- .Machine$double.eps^(1/2)
  n <- length(UD)

  r <- list()
  for(i in 1:n) { r[[i]] <- UD[[i]]$r }

  dr <- sapply(UD,function(ud){ud$dr}) # (axes,individuals)
  if(any(diff(dr[1,])!=0 | diff(dr[2,])!=0)) { stop("Inconsistent grid resolutions.") }
  dr <- dr[,1]

  x.min <- max( sapply(r,function(R){first(R$x)}) )
  x.max <- min( sapply(r,function(R){last(R$x)}) )

  y.min <- max( sapply(r,function(R){first(R$y)}) )
  y.max <- min( sapply(r,function(R){last(R$y)}) )

  # R <- list()
  # if(x.min<=x.max) { R$x <- seq(x.min,x.max,length.out=1+round((x.max-x.min)/dr['x'])) } # by fails with numerical error
  # if(y.min<=y.max) { R$y <- seq(y.min,y.max,length.out=1+round((y.max-y.min)/dr['y'])) } # by fails with numerical error

  SUB <- list()
  for(i in 1:n)
  {
    SUB[[i]] <- list()
    SUB[[i]]$x <- (r[[i]]$x-x.min)/dr['x'] >= -tol & (x.max-r[[i]]$x)/dr['x'] >= -tol # careful with round-off error
    SUB[[i]]$y <- (r[[i]]$y-y.min)/dr['y'] >= -tol & (y.max-r[[i]]$y)/dr['y'] >= -tol # careful with round-off error
  }

  return(SUB)
}


########
# UDs have origin-aligned grids with identical resolutions
# put their PDFs on the same grid
same.grids <- function(UD)
{
  n <- length(UD)

  SUB <- grid.intersection(UD)

  for(i in 1:n)
  {
    UD[[i]]$r$x <- UD[[i]]$r$x[SUB[[i]]$x]
    UD[[i]]$r$y <- UD[[i]]$r$y[SUB[[i]]$y]
    UD[[i]]$PDF <- UD[[i]]$PDF[SUB[[i]]$x,SUB[[i]]$y]
  }

  return(UD)
}


##############
# give grid argument 1canonical formatting - returning list(r,dr,extent)
format_grid <- function(grid,axes=c('x','y'))
{
  if(is.null(grid)) { grid <- list(axes=axes) }

  # format lazy grid arguments
  if(!is.null(grid) && class(grid)[1]=="list" && all(axes %in% names(grid))) # assuming coordinate list
  { grid <- list(r=grid) }
  else if(!is.null(grid) && class(grid)[1] %nin% c("list","UD","RasterLayer"))
  {
    if(class(grid)[1]=="Extent")
    { grid <- list(extent=grid) }
    else # assuming extent or dr
    {
      if(length(dim(grid))==2) # assuming extent
      { grid <- list(extent=grid) }
      else if(length(grid)==1 || length(grid)==length(axes)) # assuming dr
      { grid <- list(dr=grid) }
      else
      { stop("Malformed grid argument.") }
    }
  }

  # recycle resolution if necessary
  if("dr" %in% names(grid)) { grid$dr <- array(grid$dr,length(axes)); names(grid$dr) <- axes }

  if(class(grid)[1]=="RasterLayer") ### grid defined by raster ###
  {
    # this includes outer margin of pixels (don't use!)
    EXT <- raster::extent(grid)
    EXT <- cbind( c(EXT@xmin,EXT@xmax) , c(EXT@ymin,EXT@ymax) )
    dEXT <- EXT[2,]-EXT[1,]

    res <- dim(grid)[2:1]
    dr <- raster::res(grid)

    # grid locations (pixel centers)
    r <- list(x=raster::xFromCol(grid),y=rev(raster::yFromRow(grid)))

    grid <- list(dr=dr,r=r)
  } # end raster
  else if(!is.null(grid$r)) ### grid fully pre-specified ###
  {
    r <- grid$r
    # UD object will also have dr specified
    if(is.null(grid$dr)) { grid$dr <- sapply(r,function(r){mean(diff(r))}) }
  }

  # default resolution for multiple individuals with different resolutions
  if("dr.fn" %nin% names(grid)) { grid$dr.fn <- min }

  return(grid)
}


############################################
# construct a grid for the density function
# non-destructive WRT argument 'grid'
# non-grid arguments are merely suggestions
kde.grid <- function(data,H,axes=c("x","y"),alpha=0.001,res=NULL,dr=NULL,EXT=NULL,EXT.min=NULL,grid=NULL)
{
  DIM <- length(axes)
  H <- prepare.H(H,n=length(data$t),axes=axes) # (times,dim,dim)

  # how far to extend range from data as to ensure alpha significance in total probability
  z <- qmvnorm(1-alpha,DIM)

  dH <- z * apply(H,1,function(h){sqrt(diag(h))}) # (dim,times)
  dim(dH) <- c(DIM,nrow(data))
  dH <- t(dH) # (times,dim)

  if(!is.null(grid$r)) ### grid fully pre-specified ###
  {
    R <- grid$r
    # UD object will also have dr
    if("dr"%in%names(grid))
    { dr <- grid$dr }
    else
    { dr <- sapply(R,function(r){mean(diff(r))}) }
  }
  else if("dr"%in%names(grid) && !is.null(grid$extent)) ### grid fully pre-specified... with possible conflicts ###
  {
    # raster extents include pixel margins
    MARGIN <- class(grid$extent)[1]=="Extent" || DIM==1

    dr <- grid$dr
    EXT <- as.matrix(grid$extent)

    # raster extents include pixel margins
    if(MARGIN)
    {
      EXT[1,] <- EXT[1,] + dr/2
      EXT[2,] <- EXT[2,] - dr/2
    }

    if(isTRUE(grid$align.to.origin))
    {
      EXT[1,] <- floor(EXT[1,]/dr)*dr
      EXT[2,] <- ceiling(EXT[2,]/dr)*dr
    }

    dEXT <- EXT[2,]-EXT[1,]
    res <- round(dEXT/dr) # extent could be slightly off --- assuming mostly correct

    # preserve dr
    R <- lapply(1:DIM,function(i){seq(EXT[1,i],EXT[2,i],length.out=1+res[i])})
  }
  else if(!is.null(grid$extent)) ### grid extent specified, but not resolution ###
  {
    # raster extents include pixel margins
    MARGIN <- class(grid$extent)[1]=="Extent" || DIM==1

    # align.to.origin doesn't make sense without dr specified
    EXT <- as.matrix(grid$extent)
    dEXT <- EXT[2,]-EXT[1,]

    if(is.null(dr)) { dr <- dEXT/(res+MARGIN) }
    res <- ceiling(dEXT/dr) - MARGIN
    res <- max(res,1)
    dr <- dEXT/(res+MARGIN) # correct dr if necessary

    # raster extents include pixel margins
    if(MARGIN)
    {
      EXT[1,] <- EXT[1,] + dr/2
      EXT[2,] <- EXT[2,] - dr/2
    }

    R <- lapply(1:DIM,function(i){seq(EXT[1,i],EXT[2,i],length.out=1+res[i])})
  } ### end grid extent specified ###
  else if("dr"%in%names(grid)) ### grid resolution specified, but not extent ###
  {
    dr <- grid$dr
    dr <- array(dr,DIM)

    R <- get.telemetry(data,axes) # (times,dim)
    # minimum extent
    if(is.null(EXT)) { EXT <- rbind( apply(R-dH,2,min) , apply(R+dH,2,max) ) } # (ext,dim) }

    dEXT <- EXT[2,]-EXT[1,]
    res <- ceiling(dEXT/dr)
    dEXT <- res*dr

    # grid center
    mu <- apply(EXT,2,mean)

    # preserve dr and mu
    EXT[1,] <- mu - dEXT/2
    EXT[2,] <- mu + dEXT/2

    if(isTRUE(grid$align.to.origin))
    {
      EXT[1,] <- floor(EXT[1,]/dr)*dr
      EXT[2,] <- ceiling(EXT[2,]/dr)*dr
      dEXT <- EXT[2,]-EXT[1,]
      res <- round(dEXT/dr)
    }

    # add one cell buffer on all sides for occurrence with zero-error IID models
    R <- lapply(1:DIM,function(i){seq(EXT[1,i]-dr[i],EXT[2,i]+dr[i],length.out=1+res[i]+2)})
  } ### end grid resolution specified ###
  else ### grid not specified at all ###
  {
    # align.to.origin doesn't make sense without dr
    R <- get.telemetry(data,axes) # (times,dim)

    if(is.null(EXT)) { EXT <- rbind( apply(R-dH,2,min) , apply(R+dH,2,max) ) } # (ext,dim) }
    colnames(EXT) <- axes
    if(!is.null(EXT.min)) { EXT <- as.matrix(extent(list(EXT,EXT.min),level=1))[,axes,drop=FALSE] }
    dEXT <- EXT[2,]-EXT[1,]

    # grid center
    mu <- apply(EXT,2,mean)

    # complete suggestions
    if(is.null(dr)) { dr <- dEXT/res }
    res <- ceiling(dEXT/dr)
    dr <- dEXT/res

    # grid locations
    # add one cell buffer on all sides for occurrence with zero-error IID models
    R <- lapply(1:DIM,function(i){ seq(EXT[1,i]-dr[i],EXT[2,i]+dr[i],length.out=1+res[i]+2) } ) # (grid,dim)
  } ### end not grid specification at all ###
  ### END SPECIFY GRID ###

  names(R) <- axes
  names(dr) <- axes
  #colnames(EXT) <- axes
  #rownames(EXT) <- c('min','max')
  #names(dEXT) <- axes
  #grid <- list(r=R,dH=dH,dr=dr,extent=EXT,dEXT=dEXT,axes=axes)
  grid <- list(r=R,dH=dH,dr=dr,axes=axes)
  return(grid)
}
