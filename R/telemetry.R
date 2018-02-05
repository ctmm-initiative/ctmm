subset.telemetry <- function(x,...)
{
   info <- attr(x,"info")
   x <- data.frame(x)
   x <- subset.data.frame(x,...)
   x <- new.telemetry(x,info=info)
   return(x)
}

`[.telemetry` <- function(x,...)
{
  info <- attr(x,"info")
  x <- data.frame(x)
  x <- "[.data.frame"(x,...)
  # if(class(x)=="data.frame") { x <- new.telemetry(x,info=info) }
  x <- new.telemetry(x,info=info)
  return(x)
}


# rbind track segments
# rbind.telemetry <- function(...,deparse.level=1,make.row.names=TRUE,stringsAsFactors=default.stringsAsFactors())
# {
#   # sorting checks
#   x <- list(...)
#   t1 <- sapply(x,function(y){ y$t[1] })
#   t2 <- sapply(x,function(y){ last(y$t) })
#   IND1 <- sort(t1,method="quick",index.return=TRUE)$ix
#   IND2 <- sort(t2,method="quick",index.return=TRUE)$ix
#   if(any(IND1!=IND2)) { warning("Segments are overlapping and will be intermixed.") }
#   # all other sorting is segment-segment sorting, which is valid
#
#   # combine data
#   info <- mean.info(x)
#   x <- rbind.data.frame(...,stringsAsFactors=FALSE)
#   x <- new.telemetry(x,info=info)
#
#   # sort times
#   IND <- sort(x$t,method="quick",index.return=TRUE)$ix
#   x <- x[IND,]
#
#   return(x)
# }


# time-ordered sample of a track
# sample.telemetry <- function(x,size,replace=FALSE,prob=NULL)
# {
#   n <- length(x$t)
#
#   SAM <- sample(n,size,replace=replace,prob=prob)
#
#   SAM <- sort.int(SAM,method="quick")
#
#   x <- x[SAM,]
#
#   return(x)
# }


# pull out axis data
get.telemetry <- function(data,axes=c("x","y"))
{
  # z <- "[.data.frame"(data,axes)
  # z <- as.matrix(z)
  # colnames(z) <- axes
  # return(z)
  return(as.matrix(data.frame(data)[, axes], dimnames = axes))
}

#######################
# Generic import function
as.telemetry <- function(object,timeformat="",timezone="UTC",projection=NULL,UERE=NULL,drop=TRUE,...) UseMethod("as.telemetry")

# MoveStack object
as.telemetry.MoveStack <- function(object,timeformat="",timezone="UTC",projection=NULL,UERE=NULL,drop=TRUE,...)
{
  # need to first conglomerate to MoveBank format, then run as.telemetry
  object <- move::split(object)
  DATA <- lapply(object,function(mv){ Move2CSV(mv,timeformat=timeformat,timezone=timezone,projection=projection,UERE=UERE) })
  DATA <- do.call(rbind,DATA)
  DATA <- as.telemetry.data.frame(DATA,timeformat=timeformat,timezone=timezone,projection=projection,UERE=UERE,drop=drop)
  return(DATA)
}

# Move object
as.telemetry.Move <- function(object,timeformat="",timezone="UTC",projection=NULL,UERE=NULL,drop=TRUE,...)
{
  DATA <- Move2CSV(object,timeformat=timeformat,timezone=timezone,projection=projection,UERE=UERE)
  # can now treat this as a MoveBank object
  DATA <- as.telemetry.data.frame(DATA,timeformat=timeformat,timezone=timezone,projection=projection,UERE=UERE,drop=drop)
  return(DATA)
}

# convert Move object back to MoveBank CSV
Move2CSV <- function(object,timeformat="",timezone="UTC",projection=NULL,UERE=NULL,...)
{
  DATA <- data.frame(timestamp=move::timestamps(object))
  if(raster::isLonLat(object))
  { DATA[,c('location.long','location.lat')] <- sp::coordinates(object) }
  else
  { DATA[,c('location.long','location.lat')] <- sp::coordinates(sp::spTransform(object,sp::CRS("+proj=longlat +datum=WGS84"))) }

  # break Move object up into data.frame and idData
  idData <- move::idData(object)
  object <- as.data.frame(object)

  # add data.frame columns to data
  DATA <- cbind(DATA,object)

  # add idData to data
  # DATA[,names(idData)] <- t(array(idData,c(length(idData),nrow(DATA))))
  for(i in 1:length(idData)) { DATA[,names(idData)[i]] <- idData[i] }

  return(DATA)
}

# pull out a column with different possible names
pull.column <- function(object,NAMES,FUNC=as.numeric)
{
  # consider alternative spellings of NAMES, but preserve order (preference) of NAMES
  COPY <- NULL
  for(NAME in NAMES)
  { COPY <- c(COPY,unique(c(NAME,gsub("[.: ]","_",NAME)))) }
  NAMES <- tolower(COPY) # all lower case

  COLS <- names(object) # already lower case

  for(NAME in NAMES)
  {
    if(NAME %in% COLS)
    {
      # preference non-empty columns
      COL <- FUNC(object[,NAME])
      if(any(!is.na(COL))) {return(COL) }
    }
  }
  # nothing non-empty matched
  return(NULL)
}

# read in a MoveBank object file
as.telemetry.character <- function(object,timeformat="",timezone="UTC",projection=NULL,UERE=NULL,drop=TRUE,...)
{
  # read with 3 methods: fread, temp_unzip, read.csv, fall back to next if have error.
  # fread error message is lost, we can use print(e) for debugging.
  data <- tryCatch(data.table::fread(object,data.table=FALSE,check.names=TRUE,nrows=5),
                   error = function(e) "error")
  # if fread fails, then decompress zip to temp file, read data, remove temp file
  if (class(data) == "data.frame") { data <- data.table::fread(object,data.table=FALSE,check.names=TRUE,...) }
  else {
    data <- tryCatch(temp_unzip(object, data.table::fread, data.table=FALSE,check.names=TRUE,...),
                     error = function(e) "error")
    if (identical(data, "error")) {
      cat("fread failed, fall back to read.csv", "\n")
      data <- utils::read.csv(object,...)
    }
  }
  data <- as.telemetry.data.frame(data,timeformat=timeformat,timezone=timezone,projection=projection,UERE=UERE,drop=drop)
  return(data)
}

# this assumes a MoveBank data.frame
as.telemetry.data.frame <- function(object,timeformat="",timezone="UTC",projection=NULL,UERE=NULL,drop=TRUE,...)
{
  # make column names canonicalish
  names(object) <- tolower(names(object))

  # fastPOSIXct doesn't have a timeformat argument and as.POSIXct doesn't accept this argument if empty/NA/NULL ???
  if(class(object$timestamp)=="character" && timeformat=="") { DATA <- fasttime::fastPOSIXct(object$timestamp,tz=timezone) }
  else if(timeformat=="") { DATA <- as.POSIXct(object$timestamp,tz=timezone) }
  else { DATA <- as.POSIXct(object$timestamp,tz=timezone,format=timeformat) }
  DATA <- data.frame(timestamp=DATA)

  COL <- c("animal.ID","individual.local.identifier","local.identifier","Name","ID","tag.local.identifier","tag.ID","deployment.ID","trackId")
  COL <- pull.column(object,COL,as.factor)
  if(length(COL)==0)
  {
    warning("No MoveBank identification found. Assuming data corresponds to one indiviudal.")
    COL <- factor(rep('unknown',nrow(object)))
  }
  DATA$id <- COL

  COL <- c("location.long","Longitude","long","lon")
  COL <- pull.column(object,COL)
  DATA$longitude <- COL

  COL <- c("location.lat","Latitude","lat")
  COL <- pull.column(object,COL)
  DATA$latitude <- COL

  ###############################################
  # PROJECTION
  # delete missing rows for necessary information
  COLS <- c("timestamp","latitude","longitude")
  for(COL in COLS)
  {
    NAS <- which(is.na(DATA[,COL]))
    if(length(NAS))
    {
      DATA <- DATA[-NAS,]
      object <- object[-NAS,]
    }
  }

  DATA$t <- as.numeric(DATA$timestamp)

  xy <- cbind(DATA$longitude,DATA$latitude)
  colnames(xy) <- c("x","y")

  if(class(projection)=="CRS") { projection <- as.character(projection) }
  if(is.null(projection)) { projection <- suggest.projection(DATA) }
  else { validate.projection(projection) }
  xy <- rgdal::project(xy,projection)

  DATA$x <- xy[,1]
  DATA$y <- xy[,2]

  rm(xy)

  ###################################
  # HDOP
  # Import and use e-obs accuracy if available
  COL <- "eobs.horizontal.accuracy.estimate"
  COL <- pull.column(object,COL)
  if(length(COL)) { DATA$HERE <- sqrt(2)*COL }
  # I emailed them, but they didn't know if there needed to be a sqrt(2) factor here
  # Do I assume this is an HDOP sigma_H ?
  # Do I assume this is an x-y standard deviation?
  # Scott's calibration data is more like the latter

  # ARGOS error ellipse goes here
  # !!

  # Import and use HDOP if available
  COL <- c("GPS.HDOP","HDOP","DOP")
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    DATA$HDOP <- COL
    if(is.null(UERE)) { warning("HDOP values found, but UERE not specified and will have to be fit. See help(\"uere\").") }
  }

  # ARGOS error categories go here
  # !!

  # approximate DOP from # satellites if necessary
  if(!any(c("HDOP","HERE") %in% names(DATA)))
  {
    COL <- c("GPS.satellite.count","satellite.count","NumSats","Sats") # Counts?
    COL <- pull.column(object,COL)

    if(length(COL))
    {
      warning("HDOP values not found. Approximating with # satellites.")
      COL <- 10/(COL-2)
      DATA$HDOP <- COL

      if(is.null(UERE)) { warning("HDOP values approximated, but UERE not specified and will have to be fit. See help(\"uere\").") }
    }
  }

  ################################
  # HEIGHT
  # Import third axis if available
  COL <- c("height.above.ellipsoid","height.above.msl")
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    DATA$z <- COL

    COL <- c("GPS.VDOP","VDOP")
    COL <- pull.column(object,COL)
    if(length(COL))
    {
      DATA$VDOP <- COL
      if(is.null(UERE)) { warning("VDOP values found, but UERE not specified and will have to be fit. See help(\"uere\").") }
    }
    # need to know where the ground is too

    # if no VDOP, USE HDOP/HERE as approximate VDOP
    COL <- intersect(c("HERE","HDOP"),names(DATA))
    if(!any(c("VDOP","VERE") %in% names(DATA)) && length(COL))
    {
      COL <- COL[1]
      DATA$VDOP <- DATA[,COL]
      warning("VDOP not found. ",COL," used as an approximate VDOP, which will require a separate UERE. See help(\"uere\").")
    }
  }

  ########################################
  # VELOCITY
  # Import velocity information if present
  COL <- c("ground.speed","speed","GPS.speed")
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    v <- COL
    COL <- c("heading","GPS.heading")
    d.theta <- pull.column(object,COL)

    # WGS-84
    R.EQ <- 6378137
    R.PL <- 6356752.3142
    # approximate 1-meter-North latitude displacements
    d.lambda <- 1/sqrt((R.EQ*sin(DATA$latitude))^2+(R.PL*cos(DATA$latitude))^2)
    # could use grad() but would be slowwwww....
    xy <- cbind(DATA$longitude,DATA$latitude + d.lambda)
    colnames(xy) <- c("x","y")
    xy <- rgdal::project(xy,projection)
    # difference vector
    xy <- xy - DATA[,c('x','y')]
    # North heading
    theta <- atan2(xy[,1],xy[,2])
    # velocity heading
    theta <- theta + d.theta
    # velocity components
    DATA$v.x <- v * cos(theta)
    DATA$v.y <- v * sin(theta)

    rm(v,theta,d.theta,d.lambda,xy)

    ####################
    # SPEED ERE
    COL <- "eobs.speed.accuracy.estimate"
    COL <- pull.column(object,COL)
    if(length(COL)) { DATA$SERE <- sqrt(2)*COL } # assuming same form as EOBS horizontal accuracy estimate

    # if no VDOP, USE HDOP/HERE as approximate SDOP
    COL <- intersect(c("HERE","HDOP"),names(DATA))
    if(!any(c("SDOP","SERE") %in% names(DATA)) && length(COL))
    {
      COL <- COL[1]
      DATA$SDOP <- DATA[,COL]
      warning(COL," used as an approximate SDOP, which will require a separate UERE. See help(\"uere\").")
    }
  } # END SPEED IMPORT

  #######################################
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

    # delete empty columns in case collars/tags are different
    for(COL in names(telist[[i]]))
    { if(any(is.na(telist[[i]][[COL]]))) { telist[[i]][[COL]] <- NULL } }
  }
  names(telist) <- id

  # finally set the UERE if present
  if(!is.null(UERE)) { uere(telist) <- UERE }

  # return single or list
  if (n>1 || !drop) { return(telist) }
  else { return(telist[[1]]) }
}

#################
# clean up data
telemetry.clean <- function(data,id)
{
  # sort in time
  ORDER <- sort.list(data$t,na.last=NA,method="quick")
  data <- data[ORDER,]
  if(any(ORDER != 1:length(ORDER))) { warning("Times might be out of order or duplicated in ",id,". Make sure that timeformat and timezone are correctly specified.") }

  # remove duplicate observations
  ORDER <- length(data$t)
  data <- unique(data)
  if(ORDER != length(data$t)) { warning("Duplicate data in ",id," removed.") }

  # exit with warning on duplicate times
  if(anyDuplicated(data$t)) { warning("Duplicate times in ",id,". Data cannot be fit without an error model.") }

  # remove old level information
  data <- droplevels(data)

  dt <- diff(data$t)
  v <- sqrt(diff(data$x)^2+diff(data$y)^2)/dt
  v <- max(v)
  message("Maximum speed of ",format(v,digits=3)," m/s observed in ",id)
  dt <- min(dt)
  units <- unit(dt,dimension='time')
  message("Minimum sampling interval of ",format(dt/units$scale,digits=3)," ",units$name," in ",id)

  return(data)
}


################################
# ZOOM INTO TELEMETRY DATA
# this needs some work
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
  dist <- list()
  dist$name <- "meters"
  dist$scale <- 1

  if(!add)
  {
    ext <- NULL

    # bounding locations from data
    if(!is.null(data))
    { ext <- rbind(ext,extent(data)) }

    # bounding locations from UDs
    if(!is.null(UD))
    { ext <- rbind(ext,extent(UD,level=level,level.UD=level.UD)) }

    # bounding locations from Gaussian CTMM
    if(!is.null(CTMM) & !any(is.na(level.UD)))
    { ext <- rbind(ext,extent(CTMM,level=level,level.UD=level.UD)) }

    # combine ranges
    ext <- extent.telemetry(ext)

    # bounding box
    mu <- c(mean(ext$x),mean(ext$y))
    buff <- c(diff(ext$x),diff(ext$y))/2

    # now zoom in/out to some fraction of the grid
    buff <- fraction*buff

    ext$x <- mu[1] + buff[1]*c(-1,1)
    ext$y <- mu[2] + buff[2]*c(-1,1)

    # try to obey xlim/ylim if provided
    if(!is.null(xlim) || !is.null(ylim))
    {
      max.diff <- max(diff(xlim),diff(ylim))*c(-1,1)/2

      if(is.null(ylim))
      { ylim <- mu[2] + max.diff }
      else if(is.null(xlim))
      { xlim <- mu[1] + max.diff }

      ext$x <- xlim
      ext$y <- ylim
    }

    # Get best unit scale
    dist <- unit(unlist(ext),"length")

    xlab <- paste("x ", "(", dist$name, ")", sep="")
    ylab <- paste("y ", "(", dist$name, ")", sep="")

    # residuals have no units
    if(!is.null(data) && !is.null(attr(data[[1]],"info")$residual))
    {
      xlab <- "x"
      ylab <- "y"
    }

    ext <- ext/dist$scale

    # empty base layer plot
    plot(ext, xlab=xlab, ylab=ylab, col=grDevices::rgb(1,1,1,0), asp=1, ...)
  }

  return(dist)
}

#######################################
# PLOT TELEMETRY DATA
#######################################
plot.telemetry <- function(x,CTMM=NULL,UD=NULL,level.UD=0.95,level=0.95,DF="CDF",velocity=FALSE,col="red",col.level="black",col.DF="blue",col.grid="white",pch=1,type='p',labels=NULL,fraction=1,add=FALSE,xlim=NULL,ylim=NULL,cex=NULL,lwd=1,...)
{
  alpha <- 1-level
  alpha.UD <- 1-level.UD

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
      if(!is.null(CTMM[[i]]$COV) & !is.na(level))
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
    plot.UD(UD,level.UD=level.UD,level=level,DF=DF,col.level=col.level,col.DF=col.DF,col.grid=col.grid,labels=labels,fraction=fraction,add=TRUE,xlim=xlim,ylim=ylim,cex=cex,lwd=lwd,...)
  }

  #########################
  # PLOT TELEMETRY DATA

  # prepare point characteristics
  prepare.p <- function(pchar,all=FALSE)
  {
    if(!is.list(pchar))
    {
      if(length(x)>1)
      { pchar <- array(pchar,length(x)) }
      else if(!all)
      { pchar <- list(array(pchar,length(x[[1]]$t))) }
    }
    return(pchar)
  }

  col <- prepare.p(col)
  pch <- prepare.p(pch)
  type <- prepare.p(type,all=TRUE)

  # automagic the plot point size
  if(is.null(cex))
  {
    p <- sum(sapply(x, function(d) { length(d$t) } ))
    if(p>1000) { cex <- 1000/p } else { cex <- 1 }
  }
  cex <- prepare.p(cex)

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
      graphics::points(r, cex=cex[[i]], col=col[[i]], pch=pch[[i]], type=type[[i]],...)
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
    if(velocity && all(c("v.x","v.y") %in% names(x[[i]])))
    {
      dr <- x[[i]][,c("v.x","v.y")]/dist$scale*dt

      arr.length <- dr^2
      arr.length <- sqrt(arr.length[,1]+arr.length[,2])
      arr.length <- 0.1*cmpkm*arr.length

      shape::Arrows(x0=r$x, y0=r$y, x1=(r$x+dr$v.x), y1=(r$y+dr$v.y), col=col[[i]], code=2, segment=T, arr.adj=1, arr.length=arr.length, arr.type="curved")
    }
  }

}
# SET METHODS FOR PLOT.TELEMETRY
#methods::setMethod("plot",signature(x="telemetry",y="missing"), function(x,y,...) plot.telemetry(x,...))
#methods::setMethod("plot",signature(x="telemetry",y="telemetry"), function(x,y,...) plot.telemetry(list(x,y),...))
#methods::setMethod("plot",signature(x="telemetry",y="ctmm"), function(x,y,...) plot.telemetry(x,model=y,...))
#methods::setMethod("plot",signature(x="telemetry",y="UD"), function(x,y,...) plot.telemetry(x,akde=y,...))
#methods::setMethod("plot",signature(x="telemetry"), function(x,...) plot.telemetry(x,...))


##############
plot.UD <- function(x,level.UD=0.95,level=0.95,DF="CDF",col.level="black",col.DF="blue",col.grid="white",labels=NULL,fraction=1,add=FALSE,xlim=NULL,ylim=NULL,cex=NULL,lwd=1,...)
{
  if(!is.null(x)) { if(class(x)=="UD") { x <- list(x) } }

  dist <- new.plot(UD=x,fraction=fraction,add=add,xlim=xlim,ylim=ylim,...)

  # contours colour
  if(length(col.level)==length(level.UD) && length(col.level) != length(x))
  { col.level <- t(array(col.level,c(length(level.UD),length(x)))) }
  else
  { col.level <- array(col.level,c(length(x),length(level.UD))) }
  col.level <- array(col.level,c(length(x),length(level.UD),3))

  # contour labels
  if(is.null(labels)) { labels <- round(100*level.UD) }
  if((length(labels)==length(level.UD) || length(labels)==3*length(level.UD)) && length(labels) != length(x))
  {
    labels <- array(labels,c(length(level.UD),3,length(x)))
    labels <- aperm(labels,c(3,1,2))
  }
  else
  { labels <- array(labels,c(length(x),length(level.UD),3)) }

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

  # plot grid
  for(i in 1:length(x))
  {
    if(sum(diag(x[[i]]$H)>0))
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

    # not sure why this is necessary
    graphics::box(lwd=lwd,...)
  }

  # CONTOURS
  for(i in 1:length(x))
  {
    if(!any(is.na(col.level[i,,])) && !any(is.na(level.UD)))
    {
      # make sure that correct style is used for low,ML,high even in absence of lows and highs
      plot.kde(x[[i]],level=level.UD,labels=labels[i,,2],col=scales::alpha(col.level[i,,2],1),lwd=lwd,...)

      if(!is.na(level) && !is.null(x[[i]]$DOF.area))
      {
        P <- sapply(level.UD, function(l) { CI.UD(x[[i]],l,level,P=TRUE)[-2] } )
        plot.kde(x[[i]],level=P,labels=c(t(labels[i,,c(1,3)])),col=scales::alpha(c(t(col.level[i,,c(1,3)])),0.5),lwd=lwd/2,...)
      }
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
  drawlabels <- !(labels==FALSE | is.na(labels))
  graphics::contour(kde$r,z=kde$CDF,levels=level,labels=labels,labelcex=1,drawlabels=drawlabels,col=col,add=TRUE,...)

  # reinstate initial option (or default if was NULL--can't set back to NULL???)
  # if(is.null(MAX)) { MAX <- 25000 }
  # options(max.contour.segments=MAX)
}


##############################
# Plot Gaussian ctmm contours
plot.ctmm <- function(model,alpha=0.05,col="blue",...)
{
  mu <- model$mu # mean vector
  sigma <- model$sigma # covariance matrix

  Eigen <- eigen(sigma)
  std <- sqrt(Eigen$values)
  vec <- Eigen$vectors

  # confidence level = 1-alpha
  z <- sqrt(-2*log(alpha))

  num <- 100 # number of plotted points
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
  if(class(object)=="telemetry")
  {
    result <- attr(object,"info")

    dt <- stats::median(diff(object$t))
    units <- unit(dt,"time",thresh=1)
    result <- c(result,dt/units$scale)
    names(result)[length(result)] <- paste0("sampling interval (",units$name,")")

    P <- last(object$t)-object$t[1]
    units <- unit(P,"time",thresh=1)
    result <- c(result,P/units$scale)
    names(result)[length(result)] <- paste0("sampling period (",units$name,")")

    lon <- c(min(object$longitude),max(object$longitude))
    result <- c(result,list(lon=lon))
    names(result)[length(result)] <- "longitude range"

    lat <- c(min(object$latitude),max(object$latitude))
    result <- c(result,list(lat=lat))
    names(result)[length(result)] <- "latitude range"
  }
  else if(class(object)=="list")
  {
    # NAME <- sapply(object,function(o){ attr(o,"info")$identity })
    DT <- sapply(object,function(o){ stats::median(diff(o$t)) })
    P <- sapply(object,function(o){ last(o$t) - o$t[1] })
    lon <- sapply(object, function(o){ stats::median(o$longitude) })
    lat <- sapply(object, function(o){ stats::median(o$latitude) })

    result <- data.frame(interval=DT,period=P,longitude=lon,latitude=lat)

    # unit conversions
    COL <- 1
    units <- unit(stats::median(DT),"time",thresh=1,concise=TRUE)
    result[,COL] <- result[,COL]/units$scale
    colnames(result)[COL] <- paste0("interval (",units$name,")")

    COL <- 2
    units <- unit(stats::median(P),"time",thresh=1,concise=TRUE)
    result[,COL] <- result[,COL]/units$scale
    colnames(result)[COL] <- paste0("period (",units$name,")")
  }

  return(result)
}
#methods::setMethod("summary",signature(object="telemetry"), function(object,...) summary.telemetry(object,...))


##############
# BUFFALO DATA
##############
# buffalo <- as.telemetry("../DATA/Paul Cross/Kruger African Buffalo, GPS tracking, South Africa.csv")
## some bad data points
# buffalo[[6]] <- buffalo[[6]][-5720,]
# buffalo[[5]] <- buffalo[[5]][-869,]
# buffalo[[4]] <- buffalo[[4]][-606,]
# buffalo[[1]] <- buffalo[[1]][-1,]
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
