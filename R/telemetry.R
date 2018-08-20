subset.telemetry <- function(x,...)
{
   info <- attr(x,"info")
   UERE <- attr(x,"UERE")
   x <- data.frame(x)
   x <- subset.data.frame(x,...)
   x <- new.telemetry(x,info=info,UERE=UERE)
   return(x)
}

`[.telemetry` <- function(x,...)
{
  info <- attr(x,"info")
  UERE <- attr(x,"UERE")
  x <- data.frame(x)
  x <- "[.data.frame"(x,...)
  # if(class(x)=="data.frame") { x <- new.telemetry(x,info=info) }
  x <- new.telemetry(x,info=info,UERE=UERE)
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
  if(all(axes %in% names(data)))
  { data <- as.matrix(data.frame(data)[, axes], dimnames = axes) }
  else
  { data <- numeric(0) }

  return(data)
}

#######################
# Generic import function
as.telemetry <- function(object,timeformat="",timezone="UTC",projection=NULL,timeout=Inf,na.rm="row",drop=TRUE,...) UseMethod("as.telemetry")

# MoveStack object
as.telemetry.MoveStack <- function(object,timeformat="",timezone="UTC",projection=NULL,timeout=Inf,na.rm="row",drop=TRUE,...)
{
  # need to first conglomerate to MoveBank format, then run as.telemetry
  object <- move::split(object)
  DATA <- lapply(object,function(mv){ Move2CSV(mv,timeformat=timeformat,timezone=timezone,projection=projection) })
  DATA <- do.call(rbind,DATA)
  DATA <- as.telemetry.data.frame(DATA,timeformat=timeformat,timezone=timezone,projection=projection,timeout=timeout,na.rm=na.rm,drop=drop)
  return(DATA)
}

# Move object
as.telemetry.Move <- function(object,timeformat="",timezone="UTC",projection=NULL,timeout=Inf,na.rm="row",drop=TRUE,...)
{
  DATA <- Move2CSV(object,timeformat=timeformat,timezone=timezone,projection=projection)
  # can now treat this as a MoveBank object
  DATA <- as.telemetry.data.frame(DATA,timeformat=timeformat,timezone=timezone,projection=projection,timeout=timeout,na.rm=na.rm,drop=drop)
  return(DATA)
}

# convert Move object back to MoveBank CSV
Move2CSV <- function(object,timeformat="",timezone="UTC",projection=NULL,...)
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
# consider alternative spellings of NAMES, but preserve order (preference) of NAMES
pull.column <- function(object,NAMES,FUNC=as.numeric)
{
  canonical <- function(NAME) { unique(tolower(gsub("[.:_ ]","",NAME))) }

  NAMES <- canonical(NAMES)
  names(object) <- canonical(names(object)) -> COLS

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

# some columns are missing when fix is inferior class
# do we need a new location class for this
missing.class <- function(DATA,TYPE)
{
  # column to check for NAs
  if(TYPE=="speed") { COL <- "speed" }
  else if(TYPE=="vertical") { COL <- "z" }

  # some location classes are missing TYPE
  NAS <- is.na(DATA[[COL]])
  if(any(NAS))
  {
    # are these NAs specific to a location class
    if('class' %in% names(DATA))
    {
      IN <- unique(DATA$class[NAS])
      OUT <- unique(DATA$class[!NAS])
      OVER <- intersect(IN,OUT)
      # do we need an additional location class for missing TYPE?
      if(length(OVER))
      {
        OVER <- as.factor(NAS)
        levels(OVER) <- c(TYPE,"missing")
        DATA$class <- paste(as.character(DATA$class),as.character(OVER))
        DATA$class <- as.factor(DATA$class)
        rm(OVER)
      }
      # otherwise, current location classes are sufficient
    }
    else # we need a location class for missing TYPE
    {
      DATA$class <- as.factor(NAS)
      levels(DATA$class) <- c(TYPE,"missing")
    }

    # zero missing TYPE
    DATA[[COL]][NAS] <- 0
    if(TYPE=="speed") { DATA$heading[NAS] <- 0 }

    # do we have TYPE DOP
    DOP <- DOP.LIST[[TYPE]]$DOP
    VAR <- DOP.LIST[[TYPE]]$VAR

    if(!(DOP %in% names(DATA))) { DATA[[DOP]] <- 1 }
    # adjust DOP values for missing TYPE
    DATA[[DOP]][NAS] <- Inf
    # adjust calibrated errors --- shouldn't be necessary
    if(VAR %in% names(DATA)) { DATA[[VAR]][NAS] <- Inf }
  }

  return(DATA)
}


# read in a MoveBank object file
as.telemetry.character <- function(object,timeformat="",timezone="UTC",projection=NULL,timeout=Inf,na.rm="row",drop=TRUE,...)
{
  # read with 3 methods: fread, temp_unzip, read.csv, fall back to next if have error.
  # fread error message is lost, we can use print(e) for debugging.
  data <- tryCatch( suppressWarnings( data.table::fread(object,data.table=FALSE,check.names=TRUE,nrows=5) ) , error=function(e){"error"} )
  # if fread fails, then decompress zip to temp file, read data, remove temp file
  # previous data.table will generate error when reading zip, now it's warning and result is an empty data.frame.
  if(class(data) == "data.frame" && nrow(data) > 0) { data <- data.table::fread(object,data.table=FALSE,check.names=TRUE,...) }
  else {
    data <- tryCatch( temp_unzip(object, data.table::fread, data.table=FALSE,check.names=TRUE,...) , error=function(e){"error"} )
    if(identical(data,"error"))
    {
      cat("fread failed, fall back to read.csv","\n")
      data <- utils::read.csv(object,...)
    }
  }
  data <- as.telemetry.data.frame(data,timeformat=timeformat,timezone=timezone,projection=projection,timeout=timeout,na.rm=na.rm,drop=drop)
  return(data)
}


# this assumes a MoveBank data.frame
as.telemetry.data.frame <- function(object,timeformat="",timezone="UTC",projection=NULL,timeout=Inf,na.rm="row",drop=TRUE,...)
{
  na.rm <- match.arg(na.rm,c("row","col"))

  # make column names canonicalish
  names(object) <- tolower(names(object))

  # timestamp column
  COL <- c('timestamp','Acquisition.Start.Time','Acquisition.Time','time')
  COL <- pull.column(object,COL,FUNC=as.character)
  # fastPOSIXct doesn't have a timeformat argument and as.POSIXct doesn't accept this argument if empty/NA/NULL ???
  if(class(COL)=="character" && timeformat=="") { COL <- fasttime::fastPOSIXct(COL,tz=timezone) }
  else if(timeformat=="") { COL <- as.POSIXct(COL,tz=timezone) }
  else { COL <- as.POSIXct(COL,tz=timezone,format=timeformat) }
  DATA <- data.frame(timestamp=COL)

  COL <- c("animal.ID","individual.local.identifier","local.identifier","individual.ID","Name","ID","tag.local.identifier","tag.ID","deployment.ID","track.ID")
  COL <- pull.column(object,COL,as.factor)
  if(length(COL)==0)
  {
    warning("No MoveBank identification found. Assuming data corresponds to one indiviudal.")
    COL <- factor(rep('unknown',nrow(object)))
  }
  DATA$id <- COL

  COL <- c("location.long","Longitude","long","lon","GPS.Longitude")
  COL <- pull.column(object,COL)
  DATA$longitude <- COL

  COL <- c("location.lat","Latitude","lat","GPS.Latitude")
  COL <- pull.column(object,COL)
  DATA$latitude <- COL

  ###############################################
  # TIME & PROJECTION
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

  ##################################
  # ARGOS error ellipse/circle (newer ARGOS data >2011)
  COL <- "Argos.orientation"
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    DATA$COV.angle <- COL
    DATA$HDOP <- pull.column(object,"Argos.GDOP")
    # according to ARGOS, the following can be missing on <4 message data... but it seems present regardless
    DATA[[DOP.LIST$horizontal$VAR]] <- pull.column(object,"Argos.error.radius")^2/2
    DATA$COV.major <- pull.column(object,"Argos.semi.major")^2/2
    DATA$COV.minor <- pull.column(object,"Argos.semi.minor")^2/2
    # 1/2 from McClintock et al (2014) & in line with HDOP conventions
  }

  # ARGOS error categories (older ARGOS data <2011)
  # converted to error ellipses from ...
  COL <- c("Argos.location.class","Argos.lc")
  COL <- pull.column(object,COL,as.factor)
  if(!all(c("COV.angle","COV.major","COV.major") %in% names(DATA)) && length(COL))
  {
    # major axis always longitude
    DATA$COV.angle <- 90
    # error eigen variances
    ARGOS.minor <- c('3'=157,'2'=259,'1'=494, '0'=2271,'A'=762, 'B'=4596,'Z'=Inf)^2
    ARGOS.major <- c('3'=295,'2'=485,'1'=1021,'0'=3308,'A'=1244,'B'=7214,'Z'=Inf)^2
    # numbers from McClintock et al (2015)

    # error radii (average variance)
    ARGOS.radii <- (ARGOS.major+ARGOS.minor)/2

    message("ARGOS error ellipses not found. Using location class estimates from McClintock et al (2015).")
    DATA$class <- COL
    COL <- as.character(COL) # factors are weird
    DATA$COV.minor <- ARGOS.minor[COL]
    DATA$COV.major <- ARGOS.major[COL]
    DATA[[DOP.LIST$horizontal$VAR]] <- ARGOS.radii[COL]
    DATA$HDOP <- sqrt(2*ARGOS.radii)[COL]
  }

  ############
  # EOBS calibrated GPS errors
  COL <- c("eobs.horizontal.accuracy.estimate")
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    COL <- 1.1778678310260233 * COL # estimated from Scott's calibration data
    DATA$HDOP <- sqrt(2)*COL
    DATA[[DOP.LIST$horizontal$VAR]] <- COL^2
  }

  ###########
  # Telonics Gen4 GPS errors
  COL <- c("Horizontal.Error","GPS.Horizontal.Error","Telonics.Horizontal.Error")
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    # COL <- 0.14699275951173810 * COL # estimated from Patricia's calibration data

    DATA$HDOP <- sqrt(2)*COL
    # DATA[[DOP.LIST$horizontal$VAR]] <- COL^2

    # approximate UERE lower bound for cases where NA error and !NA HDOP
    COL <- c("GPS.HDOP","HDOP","DOP","Horizontal.Dilution","GPS.Horizontal.Dilution")
    COL <- pull.column(object,COL)
    if(length(COL))
    {
      NAS <- is.na(DATA$HDOP) & !is.na(COL) # errors that need to be calibrated
      if(any(NAS))
      {
        # DATA[[DOP.LIST$horizontal$VAR]][NAS] <- Inf # not yet calibrated
        DATA$HDOP[NAS] <- COL[NAS] # could use further calibration data to rescale
        # message("Telonics quick-fix HDOPs found; UERE will have to be fit. See help(\"uere\").")
      }
      rm(NAS)
    }
    else
    { stop("Telonics data detected, but missing HDOP column.") }
  }

  ###################################
  # HDOP
  if(!("HDOP" %in% names(DATA)))
  {
    COL <- c("GPS.HDOP","HDOP","DOP","GPS.Horizontal.Dilution")
    COL <- pull.column(object,COL)
    if(length(COL))
    {
      DATA$HDOP <- COL
      # message("HDOP values found; UERE will have to be fit. See help(\"uere\").")
    }
  }

  # approximate DOP from # satellites if necessary
  if(!("HDOP" %in% names(DATA)))
  {
    COL <- c("GPS.satellite.count","satellite.count","NumSats","Sats") # Counts? Messages?
    COL <- pull.column(object,COL)
    if(length(COL))
    {
      message("HDOP values not found. Approximating via # satellites.")
      COL <- 10/(COL-2)
      DATA$HDOP <- COL
      # message("HDOP values approximated; UERE will have to be fit. See help(\"uere\").")
    }
  }

  ###########################
  # generic location classes
  # includes Telonics Gen4 location classes (use with HDOP information)
  COL <- c("GPS.fix.type","fix.type","Fix.Attempt","GPS.Fix.Attempt","Telonics.Fix.Attempt")
  COL <- pull.column(object,COL,FUNC=as.factor)
  if(length(COL)) { DATA$class <- COL }

  #######################
  # timed-out fixes
  if(timeout<Inf)
  {
    COL <- c("GPS.time.to.fix","time.to.fix","fix.time","time.to.get.fix")
    COL <- pull.column(object,COL)
    if(length(COL))
    {
      if(class(timeout)=="function") { timeout <- timeout(COL) }
      COL <- (COL<timeout)
      COL <- as.factor(COL)
      levels(COL) <- c("timeout","in-time")

      if("class" %in% names(DATA)) # combine with existing class information
      {
        DATA$class <- paste(as.character(DATA$class),as.character(COL))
        DATA$class <- as.factor(DATA$class)
      }
      else # only class information so far
      { DATA$class <- COL }
    }
  }

  ################################
  # HEIGHT
  # Import third axis if available
  COL <- c("height.above.ellipsoid","height.above.msl","height.above.mean.sea.level","height.raw","height.(raw)","barometric.height","height","Argos.altitude","GPS.Altitude","altitude","barometric.depth","depth")
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    DATA$z <- COL

    COL <- c("GPS.VDOP","VDOP")
    COL <- pull.column(object,COL)
    if(length(COL))
    {
      DATA$VDOP <- COL
      # message("VDOP values found; UERE will have to be fit. See help(\"uere\").")
    }
    # need to know where the ground is too

    # if no VDOP, USE HDOP/HERE as approximate VDOP
    if(!("VDOP" %in% names(DATA)) && ("HDOP" %in% names(DATA)))
    {
      DATA$VDOP <- DATA$HDOP
      message("VDOP not found. HDOP used as an approximate VDOP.")
    }

    # do we need a location class for missing heights?
    DATA <- missing.class(DATA,"vertical")
  }

  ########################################
  # VELOCITY
  # Import velocity information if present
  COL <- c("ground.speed","speed","GPS.speed")
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    DATA$speed <- COL
    COL <- c("heading","GPS.heading","Course")
    COL <- pull.column(object,COL)
    if(length(COL)) { DATA$heading <- COL }
    else { DATA$speed <- NULL }

    ####################
    # SPEED ERE
    COL <- "eobs.speed.accuracy.estimate"
    COL <- pull.column(object,COL)
    if(length(COL)) # assuming same form as EOBS horizontal accuracy estimate
    {
      # UERE from Scott's calibration data
      COL <- 0.1089484 * COL
      DATA[[DOP.LIST$speed$DOP]] <- sqrt(2) * COL
      DATA[[DOP.LIST$speed$VAR]] <- COL^2
    }
    else if("HDOP" %in% names(DATA)) # USE HDOP as approximate SDOP
    {
      DATA$SDOP <- DATA$HDOP
      message("HDOP used as an approximate speed DOP.")
    }

    # do we need a location class for missing speeds?
    if("speed" %in% names(DATA))
    { DATA <- missing.class(DATA,"speed") }
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

    # delete empty columns in case collars/tags are different
    for(COL in names(telist[[i]]))
    {
      BAD <- is.na(telist[[i]][[COL]])
      if(all(BAD) || (na.rm=="col" && any(BAD))) { telist[[i]][[COL]] <- NULL } # delete the whole column
      else if(na.rm=="row" && any(BAD)) { telist[[i]] <- telist[[i]][!BAD,] } # only delete the rows
    }

    # delete incomplete vector data (caused by partial NAs)
    # i've seen speeds without headings...
    # COL <- c("speed","heading")
    # telist[[i]] <- rm.incomplete(telist[[i]],COL)
    # should rather return an error below?
    # COL <- c("COV.angle","COV.major","COV.minor")
    # telist[[i]] <- rm.incomplete(telist[[i]],COL)

    # drop missing levels and set UERE according to calibration
    telist[[i]] <- droplevels(telist[[i]])
    # NA UERE object given data.frame columns
    UERE <- uere.null(telist[[i]])
    # which UERE values were specified
    for(TYPE in colnames(UERE))
    {
      # is data calibrated
      VAR <- DOP.LIST[[TYPE]]$VAR
      if(VAR %in% names(telist[[i]]))
      {
        VAR <- (telist[[i]][[VAR]] < Inf) # is calibrated

        for(CLASS in levels(telist[[i]]$class))
        { if(all(VAR[telist[[i]]$class==CLASS])) { UERE[CLASS,TYPE] <- 1 } }

        if(!nlevels(telist[[i]]$class) && all(VAR)) { UERE["all",TYPE] <- 1 }
      }
    }

    # combine data.frame with ancillary info
    info <- list(identity=id[i], timezone=timezone, projection=projection)
    telist[[i]] <- new.telemetry( telist[[i]] , info=info, UERE=UERE )
  }
  rm(DATA)
  names(telist) <- id

  # determine projection without bias towards individuals
  if(is.null(projection)) { projection <- median.telemetry(telist,k=2) }
  # enforce projection
  telist <- "projection<-.list"(telist,projection)

  # return single or list
  if (n>1 || !drop) { return(telist) }
  else { return(telist[[1]]) }
}


# delete all columns of set if some columns of set are missing
rm.incomplete <- function(DF,COL)
{
  IN <- sum(COL %in% names(DF))
  if(IN==0 || IN==length(COL)) { return(DF) }
  # need to delete present COLs
  IN <- names(DF) %in% COL # COLs to delete
  DF <- DF[,!IN]
  return(DF)
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
  # v <- sqrt(diff(data$x)^2+diff(data$y)^2)/dt
  # v <- max(v)
  # message("Maximum speed of ",format(v,digits=3)," m/s observed in ",id)
  dt <- min(dt)
  units <- unit(dt,dimension='time')
  message("Minimum sampling interval of ",format(dt/units$scale,digits=3)," ",units$name," in ",id)

  return(data)
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
