#########################################
# telemetry class definition
#########################################
new.telemetry <- methods::setClass("telemetry", representation(info="list"), contains="data.frame")

subset.telemetry <- function(x,...)
{
   info <- attr(x,"info")
   x <- utils::getS3method("subset","data.frame")(x,...)
   x <- new.telemetry(x,info=info)
   return(x)
}

`[.telemetry` <- function(x,...)
{
  info <- attr(x,"info")
  x <- utils::getS3method("[","data.frame")(x,...)
  if(class(x)=="data.frame") { x <- new.telemetry(x,info=info) }
  return(x)
}

get.telemetry <- function(data,axes=c("x","y"))
{
  z <- "[.data.frame"(data,axes)
  z <- as.matrix(z)
  colnames(z) <- axes
  return(z)
}

#######################
# Generic import function
as.telemetry <- function(CSV,timezone="GMT",projection=NULL,UERE=NULL,...) UseMethod("as.telemetry")

# Move object
as.telemetry.Move <- function(CSV,timezone="GMT",projection=NULL,UERE=NULL,...)
{
  # clean this up to just dump idData into columns
  DATA <- data.frame(timestamp=CSV$timestamp,
                    location.long=CSV@coords[,1],
                    location.lat=CSV@coords[,2]
                    )
  
  # possibly empty columns
  DATA$individual.local.identifier <- CSV@idData$individual.local.identifier
  DATA$tag.local.identifier <- CSV@idData$tag.local.identifier
  DATA$eobs.horizontal.accuracy.estimate <- CSV$eobs.horizontal.accuracy.estimate
  DATA$GPS.HDOP <- CSV$GPS.HDOP
  DATA$height.above.ellipsoid <- CSV$height.above.ellipsoid
  DATA$GPS.VDOP <- CSV$GPS.VDOP
  
  DATA <- as.telemetry.data.frame(DATA,timezone=timezone,projection=projection,UERE=UERE)
  return(DATA)
}
  
# this assumes a MoveBank data.frame
as.telemetry.data.frame <- function(CSV,timezone="GMT",projection=NULL,UERE=NULL,...)
{
  # choose id strings from what's present
  id <- as.factor(CSV$individual.local.identifier)
  if(length(id)==0) { id <- as.factor(CSV$tag.local.identifier) }
  
  DATA <- data.frame(id=id,
                    timestamp=as.POSIXct(CSV$timestamp,tz=timezone),
                    longitude=as.numeric(CSV$location.long),
                    latitude=as.numeric(CSV$location.lat)
                    )

  # Import and use HDOP if available
  HDOP <- CSV$GPS.HDOP
  if(!is.null(HDOP))
  { 
    if(is.null(UERE))
    {
      warning("HDOP values found but UERE not specified. See help(\"uere\").")
      DATA$HDOP <- HDOP
    }
    else
    { DATA$HERE <- (HDOP*UERE) }
  }
  
  # Import and use e-obs accuracy if available
  HSTD <- CSV$eobs.horizontal.accuracy.estimate
  if(!is.null(HSTD)) { DATA$HERE <- sqrt(2)*HSTD }
  # I emailed them, but they didn't know if there needed to be a sqrt(2) factor here
  # Do I assume this is a sigma_H ?
  # Do I assume this is an x-y standard deviation?
  # Scott's calibration data is more like the latter
    
  # Import third axis if available
  DATA$z <- CSV$height.above.ellipsoid
  VDOP <- CSV$GPS.VDOP
  if(!is.null(VDOP))
  { 
    if(is.null(UERE))
    {
      warning("VDOP values found but UERE not specified. See help(\"uere\").")
      DATA$VDOP <- VDOP
    }
    else
    { DATA$VERE <- (VDOP*UERE) }
  }
  # need to know where the ground is too
  
  DATA <- stats::na.omit(DATA)
  DATA$t <- as.numeric(DATA$timestamp)
  
  xy <- cbind(DATA$longitude,DATA$latitude)
  colnames(xy) <- c("x","y")
  
  if(is.null(projection)) { projection <- suggest.projection(DATA) }
  else { validate.projection(projection) }
  xy <- rgdal::project(xy,projection)
  
  DATA$x <- xy[,1]
  DATA$y <- xy[,2]
  
  # do this or possibly get empty animals from subset
  DATA <- droplevels(DATA)
  
  id <- levels(DATA$id)
  n <- length(id)
  
  telist <- list()  
  for(i in 1:n)
  {
    telist[[i]] <- DATA[DATA$id==id[i],]
    telist[[i]]$id <- NULL
    
    # clean through duplicates, etc..
    telist[[i]] <- telemetry.clean(telist[[i]],id=id[i])
        
    # combine data.frame with ancillary info
    info <- list(identity=id[i], timezone=timezone, projection=projection, UERE=UERE)
    telist[[i]] <- new.telemetry( telist[[i]] , info=info )
  }
  names(telist) <- id
  
  if (n>1) { return(telist) }
  else { return(telist[[1]]) }
}

# read in a MoveBank CSV file
as.telemetry.character <- function(CSV,timezone="GMT",projection=NULL,UERE=NULL,...)
{
  data <- utils::read.csv(CSV,...)
  data <- as.telemetry.data.frame(data,timezone=timezone,projection=projection,UERE=UERE)
  return(data)
}

#################
# clean up data
telemetry.clean <- function(data,id)
{  
  # sort in time
  ORDER <- sort.list(data$t,na.last=NA,method="quick")
  data <- data[ORDER,]
  if(any(ORDER != 1:length(ORDER))) { warning("Times out of order in ",id," sorted") }
  
  # remove duplicate observations
  ORDER <- length(data$t)
  data <- unique(data)
  if(ORDER != length(data$t)) { warning("Duplicate data in ",id," removed") }
  
  # exit with warning on duplicate times
  if(anyDuplicated(data$t)) { stop("Duplicate times in ",id) }
  
  # remove old level information
  data <- droplevels(data)
  
  dt <- diff(data$t)
  v <- sqrt(diff(data$x)^2+diff(data$y)^2)/dt
  v <- max(v)
  message("Maximum speed of ",v," m/s observed in ",id)
  dt <- min(dt)
  message("Minimum sampling interval of ",dt," s in ",id)
  
  return(data)
}

########################################
# Suggest a good projection
########################################
suggest.projection <- function(data,datum="WGS84")
{
  # assume Movebank data.frame
  lon <- data$longitude
  lat <- data$latitude
  
  # as a first approximation use one-point equidistant at average geolocation
  lon_0 <- stats::median(lon)
  lat_0 <- stats::median(lat)
  proj <- paste("+proj=aeqd +lon_0=",lon_0," +lat_0=",lat_0," +datum=",datum,sep="")
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
  
  proj <- paste("+proj=tpeqd +lon_1=",mu1[1]," +lat_1=",mu1[2]," +lon_2=",mu2[1]," +lat_2=",mu2[2]," +datum=",datum,sep="")

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

################################
# ZOOM INTO TELEMETRY DATA
zoom.telemetry <- function(x,fraction=1,...)
{
  manipulate::manipulate(
  { plot(x,fraction=fraction,...) },
  fraction = manipulate::slider(0, 2.0,initial=fraction,step=1/100)
  ) 
}
methods::setMethod("zoom",signature(x="telemetry"), function(x,fraction=1,...) zoom.telemetry(x,fraction=fraction,...))
methods::setMethod("zoom",signature(x="UD"), function(x,fraction=1,...) zoom.telemetry(x,fraction=fraction,...))

##############
new.plot <- function(data=NULL,CTMM=NULL,UD=NULL,level.UD=0.95,level=0.95,fraction=1,add=FALSE,xlim=NULL,ylim=NULL,...)
{
  alpha.UD <- 1-level.UD
  alpha <- 1-level
  
  dist <- list()
  dist$name <- "meters"
  dist$scale <- 1

  if(!add)
  {
    ext.x <- NULL
    ext.y <- NULL
    
    # bounding locations from data
    if(!is.null(data))
    {
      ext.x <- range(sapply(data, function(d){ range(d$x) } ))
      ext.y <- range(sapply(data, function(d){ range(d$y) } ))
    }
    
    # bounding locations from UDs
    if(!is.null(UD))
    {
      for(i in 1:length(UD))
      {
        EXT <- rowSums(UD[[i]]$PDF) > 0
        EXT <- UD[[i]]$r$x[EXT]
        EXT <- c(EXT[1],last(EXT))
        ext.x <- range(ext.x,EXT)
        
        EXT <- colSums(UD[[i]]$PDF) > 0
        EXT <- UD[[i]]$r$y[EXT]
        EXT <- c(EXT[1],last(EXT))
        ext.y <- range(ext.y,EXT)
      }
    }
    
    # bounding locations from Gaussian CTMM
    if(!is.null(CTMM))
    {
      if(is.na(alpha.UD)) { alpha.UD <- exp(-1) } # mean area
      z <- sqrt(-2*log(alpha.UD))
      
      for(i in 1:length(CTMM))
      {
        # proportionality constants for outer CIs
        sigma <- CTMM[[i]]$sigma
        
        # capture outer contour if present
        const <- 1
        if(!is.null(CTMM[[i]]$COV))
        {
          K <- length(CTMM[[i]]$tau)
          const <- confint.ctmm(CTMM[[i]],alpha)[1,3]/sqrt(det(sigma))
        }
        
        buff <- z*sqrt(const*diag(sigma))
        
        ext.x[1] <- min(ext.x[1], CTMM[[i]]$mu[1] - buff[1])
        ext.x[2] <- max(ext.x[2], CTMM[[i]]$mu[1] + buff[1])
        
        ext.y[1] <- min(ext.y[1], CTMM[[i]]$mu[2] - buff[2])
        ext.y[2] <- max(ext.y[2], CTMM[[i]]$mu[2] + buff[2])
      }
    }
    
    # bounding box
    mu <- c(mean(ext.x),mean(ext.y))
    buff <- c(diff(ext.x),diff(ext.y))/2
    
    # now zoom in/out to some fraction of the grid
    buff <- fraction*buff
    
    ext.x <- mu[1] + buff[1]*c(-1,1)
    ext.y <- mu[2] + buff[2]*c(-1,1)
    
    # try to obey xlim/ylim if provided
    if(!is.null(xlim) || !is.null(ylim))
    {
      max.diff <- max(diff(xlim),diff(ylim))*c(-1,1)/2
      
      if(is.null(ylim))
      { ylim <- mu[2] + max.diff }
      else if(is.null(xlim))
      { xlim <- mu[1] + max.diff }
      
      ext.x <- xlim
      ext.y <- ylim
    }
    
    # Get best unit scale
    dist <- unit(abs(c(ext.x,ext.y)),"length")

    xlab <- paste("x ", "(", dist$name, ")", sep="")
    ylab <- paste("y ", "(", dist$name, ")", sep="")
    
    mu <- mu/dist$scale
    
    ext.x <- ext.x/dist$scale
    ext.y <- ext.y/dist$scale
    
    # empty base layer plot
    plot(ext.x,ext.y, xlab=xlab, ylab=ylab, col=grDevices::rgb(1,1,1,0), asp=1, ...)
  }
  
  return(dist)
}

#######################################
# PLOT TELEMETRY DATA
#######################################
plot.telemetry <- function(x,CTMM=NULL,UD=NULL,level.UD=0.95,level=0.95,DF="CDF",col="red",col.level="black",col.DF="blue",col.grid="grey",fraction=1,add=FALSE,xlim=NULL,ylim=NULL,cex=1,lwd=1,...)
{
  alpha <- 1-level
  alpha.UD <- 1-level.UD
  if(is.na(alpha.UD)) { alpha.UD <- exp(-1) } # mean area
  
  # listify everything for generality
  if(class(x)=="telemetry" || class(x)=="data.frame") { x <- list(x)  }
  if(!is.null(CTMM)) { if(class(CTMM)=="ctmm") { CTMM <- list(CTMM) } }
  if(!is.null(UD)) { if(class(UD)=="UD") { UD <- list(UD) } }
  
  # median time step of data
  dt <- lapply(x,function(X){diff(X$t)})
  dt <- stats::median(unlist(dt))
  
  dist <- new.plot(data=x,CTMM=CTMM,UD=UD,level.UD=level.UD,level=level,fraction=fraction,add=add,xlim=xlim,ylim=ylim,...)
  
  # plot cm per unit of distance plotted (m or km)
  cmpkm <- 2.54*mean(graphics::par("fin")*diff(graphics::par("plt"))[-2]/diff(graphics::par("usr"))[-2])
  # plot px per unit of distance plotted (m or km)
  pxpkm <- mean(grDevices::dev.size("px")*diff(graphics::par("plt"))[-2]/diff(graphics::par("usr"))[-2])
  
  #########################
  # PLOT GAUSSIAN CONTOURS AND DENSITY
  if(!is.null(CTMM))
  {
    # contours colour
    col.level <- array(col.level,length(CTMM))
    col.DF <- array(col.DF,length(CTMM))
    
    for(i in 1:length(CTMM))
    {
      # scale units
      CTMM[[i]] <- unit.ctmm(CTMM[[i]],dist$scale)
      
      # plot denisty function lazily reusing KDE code
      pdf <- kde(data.frame(CTMM[[i]]$mu[1,,drop=FALSE]),H=methods::getDataPart(CTMM[[i]]$sigma),axes=c("x","y"),res=500)
      plot.df(pdf,DF=DF,col=col.DF[[i]],...)
      
      # plot ML estimate, regular style
      plot.ctmm(CTMM[[i]],alpha.UD,col=col.level[[i]],lwd=lwd,...)
      
      # plot CIs dashed if present
      if(!is.null(CTMM[[i]]$COV))
      {
        # proportionality constants for outer CIs
        const <- confint.ctmm(CTMM[[i]],alpha)["area",]
        const <- const[c(1,3)]/const[2]
        sigma <- CTMM[[i]]$sigma
        
        for(j in 1:2)
        {
          CTMM[[i]]$sigma <- const[j]*sigma
          plot.ctmm(CTMM[[i]],alpha.UD,col=scales::alpha(col.level[[i]],0.5),lwd=lwd/2,...)
        }
      }
    }
  }

  ##########################################
  # PLOT KDE CONTOURS... AND DENSITY
  if(!is.null(UD))
  {
    UD <- lapply(UD,function(ud){ unit.UD(ud,length=dist$scale) })
    plot.UD(UD,level.UD=level.UD,level=level,DF=DF,col.level=col.level,col.DF=col.DF,col.grid=col.grid,fraction=fraction,add=TRUE,xlim=xlim,ylim=ylim,cex=cex,lwd=lwd,...)
  }
  
  #########################
  # PLOT TELEMETRY DATA
  
  # color array for plots
  col <- array(col,length(x))
  
  # automagic the plot point size
  p <- sum(sapply(x, function(d) { length(d$t) } ))
  if(p>1000) { cex <- 1000/p * cex }
  
  # minimum error
  suppressWarnings(MIN <- min(sapply(x,function(X){min(X$HERE)})))
  # minimum area
  MIN <- pi*MIN^2/2/dist$scale^2
  
  # scale error to level.UD radius
  z <- sqrt(-2*log(alpha.UD))
  MIN <- z*MIN
  
  for(i in 1:length(x))
  { 
    r <- x[[i]][,c("x","y")]/dist$scale
    
    if(is.null(x[[i]]$HERE))
    {
      graphics::points(r, cex=cex, col=col[[i]],...)
    }
    else 
    {
      # circle radius
      x[[i]]$HERE <- z*x[[i]]$HERE/sqrt(2)/dist$scale
      # color density proportional to true density
      alpha <- clamp(max((1/pxpkm)^2,MIN)/(pi*x[[i]]$HERE^2))
      bg <- scales::alpha(col[[i]],alpha)
      fg <- scales::alpha(col[[i]],alpha/2)
      graphics::symbols(x=r$x,y=r$y,circles=x[[i]]$HERE,fg=fg,bg=bg,inches=FALSE,add=TRUE,...)
    }
    
    # also plot velocity vectors at dt scale
    #     if(all(c("vx","vy") %in% names(x[[i]])))
    #     {
    #       dr <- x[[i]][,c("vx","vy")]/dist$scale*dt
    #       
    #       arr.length <- dr^2
    #       arr.length <- sqrt(arr.length[,1]+arr.length[,2])
    #       arr.length <- 0.1*cmpkm*arr.length
    # 
    #       shape::Arrows(x0=r$x, y0=r$y, x1=(r$x+dr$vx), y1=(r$y+dr$vy), col=col[[i]], code=2, segment=T, arr.adj=1, arr.length=arr.length, arr.type="curved")
    #     }
  }
  
}
# SET METHODS FOR PLOT.TELEMETRY
#methods::setMethod("plot",signature(x="telemetry",y="missing"), function(x,y,...) plot.telemetry(x,...))
#methods::setMethod("plot",signature(x="telemetry",y="telemetry"), function(x,y,...) plot.telemetry(list(x,y),...))
#methods::setMethod("plot",signature(x="telemetry",y="ctmm"), function(x,y,...) plot.telemetry(x,model=y,...))
#methods::setMethod("plot",signature(x="telemetry",y="UD"), function(x,y,...) plot.telemetry(x,akde=y,...))
#methods::setMethod("plot",signature(x="telemetry"), function(x,...) plot.telemetry(x,...))


##############
plot.UD <- function(x,level.UD=0.95,level=0.95,DF="CDF",col.level="black",col.DF="blue",col.grid="grey",fraction=1,add=FALSE,xlim=NULL,ylim=NULL,cex=1,lwd=1,...)
{
  if(!is.null(x)) { if(class(x)=="UD") { x <- list(x) } }
  
  dist <- new.plot(UD=x,fraction=fraction,add=add,xlim=xlim,ylim=ylim,...)
  
  # contours colour
  col.level <- array(col.level,length(x))
  col.DF <- array(col.DF,length(x))
  col.grid <- array(col.grid,length(x))
  
  # UNIT CONVERSIONS
  for(i in 1:length(x))
  {
    # unit conversion
    x[[i]] <- unit.UD(x[[i]],length=dist$scale)
    
    # ML DENSITY PLOTS
    plot.df(x[[i]],DF=DF,col=col.DF[[i]],...)
  }
  
  # CONTOURS
  for(i in 1:length(x))
  {
    if(!is.na(col.level[[i]]))
    {
      # make sure that correct style is used for low,ML,high even in absence of lows and highs
      if(is.na(level.UD)) { level.UD <- CI.UD(x[[i]],level.UD,level,P=TRUE)[2] } # ugly/redundant code for now!!!
      plot.kde(x[[i]],level=level.UD,col=scales::alpha(col.level[[i]],1),lwd=lwd,...)
      
      if(!is.null(x[[i]]$DOF.H))
      {
        P <- CI.UD(x[[i]],level.UD,level,P=TRUE)
        plot.kde(x[[i]],level=P[-2],labels=round(100*P[2]),col=scales::alpha(col.level[[i]],0.5),lwd=lwd/2,...)
      }
    }
  }
  
  # plot grid
  for(i in 1:length(x))
  {
    if(!is.null(x[[i]]$DOF.H))
    {    
      H <- covm(x[[i]]$H)
      theta <- H@par["angle"]
      ecc <- H@par["eccentricity"]
      sigma <- H@par["area"]
      
      X <- x[[i]]$r$x
      Y <- x[[i]]$r$y
      
      COS <- cos(theta)
      SIN <- sin(theta)
      
      # grid spacing
      du <- sqrt(sigma*exp(+ecc/2))
      dv <- sqrt(sigma*exp(-ecc/2))
      
      # bandwidth axes
      u <- outer(+X*COS,+Y*SIN,"+")
      v <- outer(-X*SIN,+Y*COS,"+")
      
      # extent of data
      B <- (x[[i]]$PDF > 0)
      u <- u[B]
      v <- v[B]
      
      ex.u <- range(u)
      ex.v <- range(v)
      
      mu.u <- mean(ex.u)
      mu.v <- mean(ex.v)
      
      # grid numbers
      n.u <- diff(ex.u)/du
      n.v <- diff(ex.v)/dv
      
      n.u <- ceiling(n.u/2)
      n.v <- ceiling(n.v/2)
     
      # grid nodes
      u <- mu.u + du*(-n.u):n.u
      v <- mu.v + dv*(-n.v):n.v
      
      # transform back
      X <- outer(u*COS,-v*SIN,"+")
      Y <- outer(u*SIN,+v*COS,"+")
      
      for(j in 1:length(u)) { graphics::segments(x0=X[j,1],y0=Y[j,1],x1=last(X[j,]),y1=last(Y[j,]),col=col.grid[i],...) }
      for(j in 1:length(v)) { graphics::segments(x0=X[1,j],y0=Y[1,j],x1=last(X[,j]),y1=last(Y[,j]),col=col.grid[i],...) }
    }
  }
  
}

##################################
# plot PDF stored as KDE object
plot.df <- function(kde,DF="CDF",col="blue",...)
{
  col <- scales::alpha(col,(0:255)/255)

  if(DF=="PDF")
  {
    zlim <- c(0,max(kde$PDF))
  }
  else if(DF=="CDF")
  {
    zlim <- c(0,1)
    kde$CDF <- 1 - kde$CDF
  }

  graphics::image(kde$r,z=kde[[DF]],useRaster=TRUE,zlim=zlim,col=col,add=TRUE,...)
}


#############################
# Plot a KDE object's contours
plot.kde <- function(kde,level=0.95,labels=round(level*100),col="black",...)
{
  # record current option
  # MAX <- getOption("max.contour.segments")

  # do something that works
  options(max.contour.segments=.Machine$integer.max)
  graphics::contour(kde$r,z=kde$CDF,levels=level,labels=labels,labelcex=1,col=col,add=TRUE,...)
  
  # reinstate initial option (or default if was NULL--can't set back to NULL???)
  # if(is.null(MAX)) { MAX <- 25000 }
  # options(max.contour.segments=MAX)
}


##############################
# Plot Gaussian ctmm contours
plot.ctmm <- function(model,alpha=0.05,col="blue",...)
{
  mu <- model$mu
  sigma <- model$sigma
  
  Eigen <- eigen(sigma)
  std <- sqrt(Eigen$values)
  vec <- Eigen$vectors
  
  z <- sqrt(-2*log(alpha))
  
  num <- 100
  theta <- 2*pi*(0:num)/(num+1)
  Sin <- sin(theta)
  Cos <- cos(theta)
  
  x <- mu[1] + z*(Cos*std[1]*vec[1,1] + Sin*std[2]*vec[1,2])
  y <- mu[2] + z*(Cos*std[1]*vec[2,1] + Sin*std[2]*vec[2,2])
  
  graphics::xspline(x, y=y, shape=-1, open=FALSE, border=col, ...)
}


########################
# summarize telemetry data
summary.telemetry <- function(object,...)
{
  result <- attr(object,"info")
  
  dt <- stats::median(diff(object$t))
  units <- unit(dt,"time",thresh=1)
  result <- c(result,dt/units$scale)
  names(result)[length(result)] <- paste("sampling interval (",units$name,")",sep="")
  
  T <- last(object$t)-object$t[1]
  units <- unit(T,"time",thresh=1)
  result <- c(result,T/units$scale)
  names(result)[length(result)] <- paste("sampling period (",units$name,")",sep="")
  
  lon <- c(min(object$longitude),max(object$longitude))
  result <- c(result,list(lon=lon))
  names(result)[length(result)] <- paste("longitude range")
  
  lat <- c(min(object$latitude),max(object$latitude))
  result <- c(result,list(lat=lat))
  names(result)[length(result)] <- paste("latitude range")
  
  return(result)
}
#methods::setMethod("summary",signature(object="telemetry"), function(object,...) summary.telemetry(object,...))


#########################
# convert to spatialpoints object
SpatialPoints.telemetry <- function(object,...)
{
  data <- object
  
  if(class(data)=="telemetry" || class(data)=="data.frame")
  {
    return( sp::SpatialPoints( "[.data.frame"(data,c("x","y")), proj4string=sp::CRS(attr(data,"info")$projection) ) )
  }
  else if(class(data)=="list")
  {
    SPL <- lapply( data, function(d) { sp::SpatialPoints( "[.data.frame"(d,c("x","y")), proj4string=sp::CRS(attr(d,"info")$projection) ) } )
    return(SPL)
  }
}
#methods::setMethod("SpatialPoints",signature(coords="telemetry"), function(coords) SpatialPoints.telemetry(coords))


##############
# BUFFALO DATA
##############
# buffalo <- as.telemetry("../Data/buffalo/Kruger African Buffalo, GPS tracking, South Africa.csv")
## some bad data points
# buffalo[[6]] <- ctmm:::new.telemetry(buffalo[[6]][-5720,],info=attr(buffalo[[6]],"info"))
# buffalo[[5]] <- ctmm:::new.telemetry(buffalo[[5]][-869,],info=attr(buffalo[[5]],"info"))
# buffalo[[4]] <- ctmm:::new.telemetry(buffalo[[4]][-606,],info=attr(buffalo[[4]],"info"))
# save(buffalo,file="data/buffalo.rda",compress="xz")


##############
# GAZELLE DATA
##############
# gazelle <- read.csv("../Data/gazelles/data_semiVarianceToIdentifyingMovementModes.csv")
# gazelle$t <- gazelle$time ; gazelle$time <- NULL
# gazelle$id <- gazelle$gazelle ; gazelle$gazelle <- NULL
# ids <- levels(gazelle$id)
# g <- list()
# for(i in 1:length(ids)) { g[[i]] <- ctmm:::new.telemetry(droplevels(gazelle[gazelle$id==ids[i],][,1:3]),info=list(identity=ids[i])) }
# names(g) <- ids
# gazelle <- g ; rm(g)
# save(gazelle,file="data/gazelle.rda",compress="xz")
