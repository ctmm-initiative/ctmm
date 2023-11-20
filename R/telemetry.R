DATUM <- "+proj=longlat +datum=WGS84"

ATTRIBUTE <- list()
ATTRIBUTE$timestamp <- c('timestamp','timestamp.of.fix','Acquisition.Time',
                         'Date.Time','Date.Time.GMT','UTC.Date.Time',"DT.TM",'Ser.Local','GPS_YYYY.MM.DD_HH.MM.SS',
                         'Acquisition.Start.Time','start.timestamp',
                         'Time.GMT','GMT.Time','Local.Time','time',"\u6642\u523B",
                         'Date.Time','Date.GMT','Date','Date.Local',"\u65E5\u4ED8",
                         't','t_dat','use_date',"event.Date","observation.Date")
ATTRIBUTE$id <- c("animal.ID","individual.local.identifier","local.identifier","individual.ID","Name","ID","ID.Names","Animal","Full.ID",
                  "tag.local.identifier","tag.ID","band.number","band.num","device.info.serial","Device.ID","collar.id","Logger","Logger.ID",
                  "Deployment","deployment.ID","track.ID")
ATTRIBUTE$taxa <- c("verbatim.Scientific.Name")
ATTRIBUTE$long <- c("location.longitude","location.long","Longitude","longitude.WGS84","Longitude.deg","long","lon","lng","GPS.Longitude","\u7D4C\u5EA6","decimal.Longitude")
ATTRIBUTE$lat <- c("location.latitude","location.lat","Latitude","latitude.WGS84","Latitude.deg","latt","lat","GPS.Latitude","\u7DEF\u5EA6","decimal.Latitude")
ATTRIBUTE$zone <- c("GPS.UTM.zone","UTM.zone","zone")
ATTRIBUTE$east <- c("GPS.UTM.Easting","GPS.UTM.East","GPS.UTM.x","UTM.Easting","UTM.East","UTM.E","UTM.x","Easting","East","x")
ATTRIBUTE$north <- c("GPS.UTM.Northing","GPS.UTM.North","GPS.UTM.y","UTM.Northing","UTM.North","UTM.N","UTM.y","Northing","North","y")
ATTRIBUTE$error <- c("eobs.horizontal.accuracy.estimate","eobs.horizontal.accuracy.estimate.m","eobs.horizontal.accuracy",
                     "gps.horizontal.accuracy.estimate","gps.horizontal.accuracy.estimate.m","gps.horizontal.accuracy",
                     "horizontal.accuracy.estimate","horizontal.accuracy.estimate.m","horizontal.accuracy",
                     "error","error.m","3D.error.m","location.error","location.error.m","HEPE","EPE","EHPE",
                     "\u8AA4\u5DEE","\u8AA4\u5DEE.m","\u8AA4\u5DEE\uFF08m\uFF09","coordinate.Uncertainty.In.Meters")
ATTRIBUTE$Telonics <- c("Horizontal.Error","GPS.Horizontal.Error","Telonics.Horizontal.Error")
ATTRIBUTE$HDOP <- c("GPS.HDOP","HDOP","Horizontal.DOP","GPS.Horizontal.Dilution","Horizontal.Dilution","Hor.Dil","Hor.DOP","HPE","coordinate.Precision")
ATTRIBUTE$DOP <- c("GPS.DOP","DOP","GPS.Dilution","Dilution","Dil")
ATTRIBUTE$PDOP <- c("GPS.PDOP","PDOP","Position.DOP","GPS.Position.Dilution","Position.Dilution","Pos.Dil","Pos.DOP")
ATTRIBUTE$GDOP <- c("GPS.GDOP","GDOP","Geometric.DOP","GPS.Geometric.Dilution","Geometric.Dilution","Geo.Dil","Geo.DOP")
ATTRIBUTE$VDOP <- c("GPS.VDOP","VDOP","Vertical.DOP","GPS.Vertical.Dilusion","Vertical.Dilution","Ver.Dil","Ver.DOP","elevation.Accuracy","depth.Accuracy")
ATTRIBUTE$nsat <- c("GPS.satellite.count","satellite.count","Sat.Count","Number.of.Sats","Num.Sats","Nr.Sat","NSat","NSats","Sat.Num","satellites.used","Satellites","Sats","SVs.in.use") # Counts? Messages?
ATTRIBUTE$FIX <- c("GPS.fix.type","GPS.fix.type.raw","fix.type","type.of.fix","e.obs.type.of.fix","Fix.Attempt","GPS.Fix.Attempt","Telonics.Fix.Attempt","Fix.Status","sensor.type","Fix","eobs.type.of.fix","2D/3D","X3.equals.3dfix.2.equals.2dfix","Nav","Validated","VALID")
ATTRIBUTE$TTF <- c("GPS.time.to.fix","time.to.fix","time.to.GPS.fix","time.to.GPS.fix.s","GPS.TTF","TTF","GPS.fix.time","fix.time","time.to.get.fix","used.time.to.get.fix","e.obs.used.time.to.get.fix","Duration","GPS.navigation.time","navigation.time","Time.On","Searching.Time")
ATTRIBUTE$z <- c("height.above.ellipsoid","height.above.elipsoid","height.above.ellipsoid.m","height.above.elipsoid.m","height.above.msl","height.above.mean.sea.level","height.raw","height","height.m","barometric.height","altimeter","altimeter.m","Argos.altitude","GPS.Altitude","MSL_altitude_m","altitude","altitude.m","Alt","barometric.depth","depth","elevation","elevation.m","elev")
ATTRIBUTE$v <- c("ground.speed",'speed.over.ground','speed.over.ground.m.s',"speed","GPS.speed")
ATTRIBUTE$heading <- c("heading","heading.degree","heading.degrees","GPS.heading","Course","direction","direction.deg")
ATTRIBUTE$outliers <- c("manually.marked.outlier","algorithm.marked.outlier","import.marked.outlier","marked.outlier","outlier")
ATTRIBUTE$COV.angle <- c("Argos.orientation","Error.ellipse.orientation")
ATTRIBUTE$COV.major <- c("Argos.semi.major","Error.semi-major.axis")
ATTRIBUTE$COV.minor <- c("Argos.semi.minor","Error.semi-minor.axis")
ATTRIBUTE$COV.mean <- c("Argos.error.radius","Error.radius")
ATTRIBUTE$COV.xx <- c("VAR.x","COV.xx")
ATTRIBUTE$COV.yy <- c("VAR.y","COV.yy")
ATTRIBUTE$COV.xy <- c("COV.xy")


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
  # if(class(x)[1]=="data.frame") { x <- new.telemetry(x,info=info) }
  x <- new.telemetry(x,info=info,UERE=UERE)
  return(x)
}

head.telemetry <- function(x, n = 6L, ...) { utils::head(data.frame(x),n=n,...) }
tail.telemetry <- function(x, n = 6L, ...) { utils::tail(data.frame(x),n=n,...) }

# rbind track segments
tbind <- function(...)
{
  x <- list(...)
  if(length(x)==1)
  {
    x <- x[[1]]
    if(class(x)[1]=="telemetry") { return(x) }
  }

  PROJS <- length(projection(x))
  info <- mean_info(x)

  # unique names - just in case
  if(is.null(names(x)))
  { names(x) <- sapply(1:length(x),function(i){paste(attr(x[[i]],'info')$identity,i)}) }

  for(i in 1:length(x)) { if("class" %nin% names(x[[i]])) { x[[i]]$class <- as.factor('all') } }

  x <- ubind(x)
  UERE <- uere(x)

  n <- sapply(x,nrow)
  n <- c(0,cumsum(n))

  COLS <- lapply(x,names)
  COLS <- do.call(c,COLS)
  COLS <- unique(COLS)

  # somehow R does not have built-in functionality to rbind data.frames with extra/missing columns...
  y <- data.frame(stringsAsFactors=FALSE)
  for(col in COLS)
  {
    for(i in 1:length(x))
    {
      r1 <- 1 + n[i]
      r2 <- n[i+1]
      if(col %in% names(x[[i]]))
      {
        if(col=="class") # factor
        { y[r1:r2,col] <- as.character(x[[i]][[col]]) }
        else
        { y[r1:r2,col] <- x[[i]][[col]] }
      }
      else
      { y[r1:r2,col] <- NA }
    }
  }
  rm(x)

  # make sure class and UERE are in the same order
  if("class" %in% names(y))
  {
    CLASS <- unique(y$class)
    CLASS <- sort(CLASS) # sort levels
    y$class <- factor(y$class,levels=CLASS)

    UERE$UERE <- UERE$UERE[CLASS,,drop=FALSE]
    UERE$DOF <- UERE$DOF[CLASS,,drop=FALSE]
  }

  # handle missing error information
  for(i in 2:length(DOP.LIST))
  {
    # missing DOP <- 1
    if(DOP.LIST[[i]]$DOP %in% COLS)
    {
      NAS <- is.na(y[[DOP.LIST[[i]]$DOP]])
      if(any(NAS)) { y[[DOP.LIST[[i]]$DOP]][NAS] <- 1 }
    }

    # missing variance <- infinite
    if(DOP.LIST[[i]]$VAR %in% COLS)
    {
      NAS <- is.na(y[[DOP.LIST[[i]]$VAR]])
      if(any(NAS)) { y[[DOP.LIST[[i]]$VAR]][NAS] <- Inf }
    }

    # missing error ellipse <- error circle
    if(DOP.LIST[[i]]$COV[1] %in% COLS)
    {
      NAS <- is.na(y[[DOP.LIST[[i]]$COV[1]]])
      if(any(NAS))
      {
        y[[DOP.LIST[[i]]$COV[1]]][NAS] <- y[[DOP.LIST[[i]]$VAR]][NAS] # COV.x.x
        y[[DOP.LIST[[i]]$COV[2]]][NAS] <- 0 # COV.x.y
        y[[DOP.LIST[[i]]$COV[3]]][NAS] <- y[[DOP.LIST[[i]]$VAR]][NAS] # COV.x.y

        if(DOP.LIST[[i]]$COV.geo[1] %in% COLS)
        {
          y[[DOP.LIST[[i]]$COV.geo[1]]][NAS] <- y[[DOP.LIST[[i]]$VAR]][NAS] # COV.major
          y[[DOP.LIST[[i]]$COV.geo[2]]][NAS] <- y[[DOP.LIST[[i]]$VAR]][NAS] # COV.minor
          y[[DOP.LIST[[i]]$COV.geo[3]]][NAS] <- 0 # COV.angle
        }
      }
    }
  }

  # sort times
  IND <- sort(y$t,index.return=TRUE)$ix
  y <- y[IND,]

  y <- new.telemetry(y,info=info,UERE=UERE)

  # inconsistent projections
  if(PROJS>1)
  {
    warning("Re-projecting due to inconsistent projections.")
    projection(y) <- median(y,k=2)
  }

  return(y)
}

# bind UERE calibrations
ubind <- function(x)
{
  UERE <- uere(x) # initial UERE information
  if(class(UERE)[1]=="list") # non-unique calibration information
  {
    # which calibrations are equivalent
    # EQUAL <- outer(UERE,identical) # annoyingly, this does not work
    EQUAL <- array(TRUE,c(length(x),length(x)))
    for(i in 1:length(x))
    {
      for(j in 1%:%(i-1))
      { EQUAL[i,j] <- EQUAL[j,i] <- identical(UERE[[i]],UERE[[j]]) }
    }

    # are any in-equivalent datasets sharing the same class names?
    RENAME <- rep(FALSE,length(x))
    for(i in 1:length(x))
    {
      for(j in 1%:%(i-1))
      {
        if(!EQUAL[i,j] && any(classnames(UERE[[i]]) %in% classnames(UERE[[j]])))
        { RENAME[i] <- RENAME[j] <- TRUE }
      }
    }

    # some devices have the same class names, but different calibration information
    if(any(RENAME))
    {
      warning("Inconsistent location classes renamed.")
      for(i in 1:length(x))
      {
        LEVELS <- classnames(UERE[[i]]) %in% levels(x[[i]]$class)
        CLASS <- paste(names(x)[i],classnames(UERE[[i]]))
        classnames(UERE[[i]]) <- CLASS

        if('class' %nin% names(x[[i]])) { x[[i]]$class <- as.factor(CLASS) }
        else { levels(x[[i]]$class) <- CLASS[LEVELS] }
      }
    }

    # fix missing type information
    TYPES <- lapply(UERE,function(U){colnames(U$UERE)})
    TYPES <- do.call(c,TYPES)
    TYPES <- unique(TYPES)
    for(i in 1:length(UERE))
    {
      for(TYPE in TYPES)
      {
        if(TYPE %nin% colnames(UERE[[i]]$UERE))
        {
          NAS <- UERE[[i]]$UERE[,1,drop=FALSE]
          NAS[] <- Inf # no precisionf or missing data type
          colnames(NAS) <- TYPE
          UERE[[i]]$UERE <- cbind(UERE[[i]]$UERE,NAS)
          UERE[[i]]$DOF <- cbind(UERE[[i]]$DOF,NAS)

          NAS <- c(NAS)
          names(NAS) <- TYPE
          UERE[[i]]$AICc <- c(UERE[[i]]$AICc,NAS)
          UERE[[i]]$Zsq <- c(UERE[[i]]$Zsq,NAS)
          UERE[[i]]$VAR.Zsq <- c(UERE[[i]]$VAR.Zsq,NAS)
          NAS[] <- 0
          UERE[[i]]$N <- c(UERE[[i]]$N,NAS)
        }
      } # end TYPE loop

      # enforce consistent type sorting for rbind()
      UERE[[i]]$UERE <- UERE[[i]]$UERE[,TYPES,drop=FALSE] # R drops rownames here too?
      UERE[[i]]$DOF <- UERE[[i]]$DOF[,TYPES,drop=FALSE]
      UERE[[i]]$AICc <- UERE[[i]]$AICc[TYPES,drop=FALSE]
      UERE[[i]]$Zsq <- UERE[[i]]$Zsq[TYPES,drop=FALSE]
      UERE[[i]]$VAR.Zsq <- UERE[[i]]$VAR.Zsq[TYPES,drop=FALSE]
      UERE[[i]]$N <- UERE[[i]]$N[TYPES,drop=FALSE]
    } # end individual loop

    # combine UEREs
    if(class(UERE)[1]=="list")
    {
      MUERE <- UERE[[1]]
      for(i in 2:length(UERE))
      {
        if(!any(EQUAL[i,1%:%(i-1)])) # not equal to any previous UERE already included
        {
          MUERE$UERE <- rbind(MUERE$UERE,UERE[[i]]$UERE)
          MUERE$DOF <- rbind(MUERE$DOF,UERE[[i]]$DOF)
          MUERE$AICc <- MUERE$AICc + UERE[[i]]$AICc
          MUERE$Zsq <- nant(MUERE$N*MUERE$Zsq,Inf) + nant(UERE[[i]]$N*UERE[[i]]$Zsq,Inf)
          MUERE$VAR.Zsq <- nant(MUERE$N*MUERE$VAR.Zsq,Inf) + nant(UERE[[i]]$N*UERE[[i]]$VAR.Zsq,Inf)
          MUERE$N <- MUERE$N + UERE[[i]]$N
          MUERE$Zsq <- nant(MUERE$Zsq / MUERE$N,Inf)
          MUERE$VAR.Zsq <- nant(MUERE$VAR.Zsq / MUERE$N,Inf)
        }
      }
      UERE <- MUERE
      rm(MUERE)
    } # end UERE merger
  } # end UERE list

  uere(x) <- UERE # assign
  for (i in 1:length(x)) { x[[i]]@UERE <- UERE } # but preserve class order

  return(x)
}


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


########
set.telemetry <- function(data,value,axes=colnames(value))
{
  if(length(axes)==1) { data[[axes]] <- value }
  else if(length(axes)>1) { data[,axes] <- value }

  return(data)
}

######
set.name <- function(x,NAMES=NULL,drop=TRUE)
{
  x <- listify(x)
  for(i in 1:length(x))
  { names(x)[i] <- ifelse(is.null(NAMES),attr(x[[i]],"info")$identity,NAMES[i]) }
  if(length(x)==1 && drop) { x <- x[[1]] }
  return(x)
}

#######################
# Generic import function
as.telemetry <- function(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",dt.hot=NA,timeout=Inf,na.rm="row",mark.rm=FALSE,keep=FALSE,drop=TRUE,...) UseMethod("as.telemetry")


# MoveStack object
as.telemetry.MoveStack <- function(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",dt.hot=NA,timeout=Inf,na.rm="row",mark.rm=FALSE,keep=FALSE,drop=TRUE,...)
{
  # get individual names
  NAMES <- move::trackId(object)
  NAMES <- levels(NAMES)

  # preserve Move object projection if possible
  if(is.null(projection) && !raster::isLonLat(object)) { projection <- raster::projection(object) }

  # convert individually
  object <- move::split(object)
  object <- lapply(object,function(mv){ as.telemetry.Move(mv,timeformat=timeformat,timezone=timezone,projection=projection,datum=datum,dt.hot=dt.hot,timeout=timeout,na.rm=na.rm,mark.rm=mark.rm,keep=keep,...) })
  # name by MoveStack convention
  names(object) <- NAMES
  for(i in 1:length(object)) { attr(object,"info")$identity <- NAMES[i] }

  # unified projection if missing (somewhat redundant, but whatever)
  if(is.null(projection)) { projection(object) <- median(object,k=2) }

  if(drop && length(object)==1) { object <- object[[1]] }

  return(object)
}


# Move object
as.telemetry.Move <- function(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",dt.hot=NA,timeout=Inf,na.rm="row",mark.rm=FALSE,keep=FALSE,drop=TRUE,...)
{
  # preserve Move object projection if possible
  if(is.null(projection) && !raster::isLonLat(object)) { projection <- raster::projection(object) }

  DATA <- Move2CSV(object,timeformat=timeformat,timezone=timezone,projection=projection)
  # can now treat this as a MoveBank object
  DATA <- as.telemetry.data.frame(DATA,timeformat=timeformat,timezone=timezone,projection=projection,dt.hot=dt.hot,timeout=timeout,na.rm=na.rm,mark.rm=mark.rm,keep=keep,drop=drop)
  return(DATA)
}


# convert Move object back to MoveBank CSV
Move2CSV <- function(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",...)
{
  DATA <- data.frame(timestamp=move::timestamps(object))

  if(raster::isLonLat(object)) # projection will be automated
  {
    DATA[,c('location.long','location.lat')] <- sp::coordinates(object)
    if(is.null(projection)) { warning("Move objects in geographic coordinates are automatically projected.") }
  }
  else # projection will be preserved, but still need long-lat
  { DATA[,c('location.long','location.lat')] <- sp::coordinates(sp::spTransform(object,sp::CRS(DATUM))) }

  # break Move object up into data.frame and idData
  idData <- move::idData(object)
  NAMES <- names(idData)
  idData2 <- attr(object,"idData") # backup
  NAMES2 <- names(idData2)

  object <- as.data.frame(object)

  # add data.frame columns to data
  DATA <- cbind(DATA,object)

  # add idData to data
  if(length(idData))
  {
    for(i in 1:length(idData))
    {
      NAME <- NAMES[i]
      if(is.null(NAME)) { NAME <- NAMES2[i] } # idData() method can drop name from attribute
      if(is.null(NAME)) { NAME <- "Animal.ID" } # just in case
      DATA[,NAME] <- idData[i]
    }
  }

  return(DATA)
}


# sf object (C) jfsmenezes & CHF [UNFINISHED]
UNFINISHED.as.telemetry.sf = function(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",dt.hot=NA,timeout=Inf,na.rm="row",mark.rm=FALSE,keep=FALSE,drop=TRUE,...)
{
  if(sf::st_geometry_type(object,by_geometry=FALSE)!="POINT")
  { stop("only point features supported at the moment") }

  coords <- sf::st_coordinates(object)

  # PROJ4 is no longer supported
  from <- sf::st_crs(object)$proj4string

  # preserve object projection if not long-lat

  # re-project long-lat if necessary

  object <- sf::st_drop_geometry(object)
  as.telemetry(object,timeformat=timeformat,timezone=timezone,projection=projection,dt.hot=dt.hot,timeout=timeout,na.rm=na.rm,mark.rm=mark.rm,keep=keep,drop=drop)
}

# make names canonical
canonical.name <- function(NAME,UNIQUE=TRUE)
{
  NAME <- gsub("[\uef\ubb\ubf]","",NAME) # remove UTF-8 Byte Order Mark (BOM)
  NAME <- gsub("[.:_ -]","",NAME) # remove separators
  NAME <- gsub("\\(|\\)","",NAME) # remove parentheses
  NAME <- gsub("\\[|\\]", "", NAME) # remove square brackets
  NAME <- gsub("\\\uFF08|\\\uFF09","",NAME) # THIS DOES NOT WORK!
  NAME <- tolower(NAME)
  if(UNIQUE) { NAME <- unique(NAME) }
  return(NAME)
}


# pull out a column with different possible names
# consider alternative spellings of NAMES, but preserve order (preference) of NAMES
pull.column <- function(object,NAMES,FUNC=as.numeric,name.only=FALSE)
{
  NAMES <- canonical.name(NAMES)
  names(object) <- canonical.name(names(object),FALSE) -> COLS

  for(NAME in NAMES)
  {
    if(NAME %in% COLS)
    {
      # preference non-empty columns
      COL <- FUNC(object[,NAME])
      NAS <- is.na(COL)
      if(!all(NAS))
      {
        if(name.only) { return(NAME) }

        # otherwise NA is not a level... which is the whole point of this
        if(class(COL)[1]=="factor" && any(NAS))
        {
          COL <- addNA(COL)
          NAS <- is.na(levels(COL))
          levels(COL)[NAS] <- "NA" # level can't literally be NA, because NA is not a value
        }

        return(COL)
      }
    }
  }
  # nothing non-empty matched
  return(NULL)
}


# some columns are missing when fix is inferior class
# do we need a new location class for this
missing.class <- function(DATA,TYPE)
{
  LEVELS <- paste0("[",TYPE,"]")
  LEVELS[2] <- paste0("[NA-",TYPE,"]")

  # column to check for NAs
  if(TYPE=="speed") { COL <- "speed" }
  else if(TYPE=="vertical") { COL <- "z" }
  else { COL <- TYPE }

  # some location classes are missing TYPE
  NAS <- is.na(DATA[[COL]])
  if(any(NAS))
  {
    # are these NAs specific to a location class
    if('class' %in% names(DATA))
    {
      MISS <- as.factor(NAS)
      levels(MISS) <- LEVELS
      # do we need an additional location class for missing TYPE?
      DATA$class <- merge.class(DATA$class,MISS)
      # otherwise, current location classes are sufficient
    }
    else # we need a location class for missing TYPE
    {
      DATA$class <- as.factor(NAS)
      levels(DATA$class) <- LEVELS
    }

    # zero missing TYPE
    DATA[[COL]][NAS] <- 0
    if(TYPE=="speed") { DATA$heading[NAS] <- 0 }

    # do we have an associated TYPE DOP
    if(TYPE %in% names(DOP.LIST))
    {
      DOP <- DOP.LIST[[TYPE]]$DOP
      VAR <- DOP.LIST[[TYPE]]$VAR

      if(!(DOP %in% names(DATA))) { DATA[[DOP]] <- 1 }
      # adjust DOP values for missing TYPE
      DATA[[DOP]][NAS] <- Inf
      # adjust calibrated errors --- shouldn't be necessary
      if(VAR %in% names(DATA)) { DATA[[VAR]][NAS] <- Inf }
    }
    else if(TYPE %in% c("HDOP","VDOP"))
    { DATA[[COL]][NAS] <- 100 }
  } # end if(any(NAS))

  return(DATA)
}


# merge two location classes [NOT with minimal levels]
merge.class <- function(class1,class2)
{
  if(is.null(class2)) { return(class1) }
  if(is.null(class1)) { return(class2) }

  # LEVEL2 <- levels(class2)
  # N <- length(LEVEL2)
  # CLASS12 <- list()
  # for(i in 1:N) { CLASS12[[i]] <- unique( class1[ class2==LEVEL2[i] ] ) }
  #
  # # do we need additional location classes
  # OVER <- matrix(0,N,N)
  # for(i in 1:N) { for(j in 1:N) { OVER[i,j] <- sum( length( intersect(CLASS12[[i]],CLASS12[[j]]) ) ) } }
  # diag(OVER) <- 0 # don't count self similarity
  # OVER <- sum(OVER)/2

  # if(OVER)
  # {
  class1 <- as.character(class1)
  class2 <- as.character(class2)
  class <- paste(class1,class2)
  class <- as.factor(class)
  # }
  # else
  # { class <- class1 }

  return(class)
}


# read in a MoveBank object file
as.telemetry.character <- function(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",dt.hot=NA,timeout=Inf,na.rm="row",mark.rm=FALSE,keep=FALSE,drop=TRUE,...)
{
  # read with 3 methods: fread, temp_unzip, read.csv, fall back to next if have error.
  # fread error message is lost, we can use print(e) for debugging.
  data <- tryCatch( suppressWarnings( data.table::fread(object,data.table=FALSE,check.names=TRUE,nrows=5) ) , error=function(e){"error"} )
  # if fread fails, then decompress zip to temp file, read data, remove temp file
  # previous data.table will generate error when reading zip, now it's warning and result is an empty data.frame.
  if(class(data)[1] == "data.frame" && nrow(data) > 0) { data <- data.table::fread(object,data.table=FALSE,check.names=TRUE,...) }
  else {
    data <- tryCatch( temp_unzip(object, data.table::fread, data.table=FALSE,check.names=TRUE,...) , error=function(e){"error"} )
    if(identical(data,"error"))
    {
      cat("fread failed, fall back to read.csv","\n")
      data <- utils::read.csv(object,...)
    }
  }
  data <- as.telemetry.data.frame(data,timeformat=timeformat,timezone=timezone,projection=projection,datum=datum,dt.hot=dt.hot,timeout=timeout,na.rm=na.rm,mark.rm=mark.rm,keep=keep,drop=drop)
  return(data)
}


# fastPOSIXct doesn't have a timeformat argument
# fastPOSIXct requires GMT timezone
# fastPOSIXct only works on times after the 1970 epoch !!!
# as.POSIXct doesn't accept this argument if empty/NA/NULL ???
asPOSIXct <- function(x,timeformat="auto",timezone="UTC",...)
{
  if(timeformat=="auto")
  {
    x <- parsedate::parse_date(x,approx=FALSE,default_tz=timezone)
    return(x)
  }

  # try fastPOSIXct
  if(class(x)[1]=="character" && timeformat=="" && timezone %in% c("UTC","GMT"))
  {
    y <- fasttime::fastPOSIXct(x,tz=timezone)
    # did fastPOSIXct fail?
    if(!any(!is.na(x) & is.na(y))) { return(y) }
  }

  if(timeformat=="") { x <- as.POSIXct(x,tz=timezone) }
  else { x <- as.POSIXct(x,tz=timezone,format=timeformat) }

  return(x)
}


# this assumes a MoveBank data.frame
as.telemetry.data.frame <- function(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",dt.hot=NA,timeout=Inf,na.rm="row",mark.rm=FALSE,keep=FALSE,drop=TRUE,...)
{
  # ARGOS have error-ellipse estimates in long-lat
  # ATLAS have error-ellipse estimates in x-y
  ARGOS <- ATLAS <- FALSE

  if(grepl("+datum=",datum,fixed=TRUE))
  {
    datum <- strsplit(datum,'+datum=',fixed=TRUE)[[1]][2]
    datum <- strsplit(datum,' +',fixed=TRUE)[[1]][1]
  }

  # get rid of tibble classes # they don't extend data.frame objects quite right
  if(class(object)[1] == "grouped_df")
  { object <- dplyr::ungroup(object) }
  if(class(object)[1] %in% c("tbl_df","tbl"))
  { object <- as.data.frame(object) }

  # UNIT CONVERSIONS
  COL <- c("speed.km.h")
  COL <- pull.column(object,COL)
  if(length(COL)) { object$speed <- COL * (1000/60^2) }

  na.rm <- match.arg(na.rm,c("row","col"))

  # make column names canonicalish
  # names(object) <- tolower(names(object))

  # marked outliers
  COL <- pull.column(object,ATTRIBUTE$outliers,as.logical)
  if(length(COL))
  {
    if(mark.rm)
    {
      COL <- is.na(COL) | !COL # ignore outlier=NA
      object <- object[COL,]
      message(sum(!COL,na.rm=T)," marked outliers removed.")
    }
    else
    { message(sum(COL,na.rm=T)," outlier markings ignored.") }
  }

  # timestamp column
  options(digits.secs=6) # otherwise, R will truncate to seconds...
  COL <- ATTRIBUTE$timestamp
  COL <- pull.column(object,COL,FUNC=as.character,name.only=TRUE)
  if("POSIXct" %in% class(object[1,COL])) # numeric timestamp
  {
    COL <- pull.column(object,COL,FUNC=identity)
    # overwrite timezone argument with existing timezone
    timezone <- format(COL,format="%Z")
    timezone <- unique(timezone)
  }
  else # character timestamp
  {
    COL <- pull.column(object,COL,FUNC=as.character)
    COL <- asPOSIXct(COL,timeformat=timeformat,timezone=timezone)
  }
  DATA <- data.frame(timestamp=COL)

  COL <- ATTRIBUTE$id
  COL <- pull.column(object,COL,as.factor)
  OCCURRENCE <- FALSE # default assumption
  if(length(COL)==0)
  {
    COL <- ATTRIBUTE$taxa
    COL <- pull.column(object,COL,as.factor)
    if(length(COL))
    { OCCURRENCE <- TRUE }
    else
    {
      warning("No MoveBank or GBIF identification found. Assuming data corresponds to one indiviudal.")
      COL <- factor(rep('unknown',nrow(object)))
    }
  }
  DATA$id <- COL

  COL <- ATTRIBUTE$long
  COL <- pull.column(object,COL)
  DATA$longitude <- COL

  COL <- ATTRIBUTE$lat
  COL <- pull.column(object,COL)
  DATA$latitude <- COL

  # UTM locations if long-lat not present
  if(all(c('longitude','latitude') %nin% names(DATA)))
  {
    message("Geocentric coordinates not found. Looking for UTM coordinates.")

    COL <- ATTRIBUTE$zone
    COL <- pull.column(object,COL,FUNC=as.character)
    zone <- COL

    COL <- ATTRIBUTE$east
    COL <- pull.column(object,COL)
    XY <- COL

    COL <- ATTRIBUTE$north
    COL <- pull.column(object,COL)
    XY <- cbind(XY,COL)

    # XY information present but zone not present
    if(!is.null(XY) && ncol(XY)==2)
    {
      DATA$x <- XY[,1]
      DATA$y <- XY[,2]

      # missing zone
      if(is.null(zone))
      {
        message('UTM zone missing; assuming UTM zone="1 +north".')
        zone <- rep("1 +north",nrow(XY))
      }

      # missing hemisphere / latitude
      # if(any(grepl("^[0-9]*$",zone)))
      # { message('UTM zone missing lattitude bands; assuming UTM hemisphere="north". Alternatively, format zone column "# +south".') }

      zone <- paste0("+proj=utm +zone=",zone," +datum=",datum)
      zone <- factor(zone)
      for(z in levels(zone))
      {
        SUB <- zone==z
        XY[SUB,] <- project(XY[SUB,,drop=FALSE],from=z)
      }
      colnames(XY) <- c("x","y")

      # # convert to long-lat
      # colnames(XY) <- c("x","y")
      # XY <- sapply(1:nrow(XY),function(i){rgdal::project(XY[i,,drop=FALSE],zone[i],inv=TRUE)})
      # XY <- t(XY)

      # work with long-lat as if imported directly
      DATA$longitude <- XY[,1]
      DATA$latitude <- XY[,2]

      ## VHF error ellipses
      COL <- ATTRIBUTE$COV.xx
      COL <- pull.column(object,COL)
      XY <- COL

      COL <- ATTRIBUTE$COV.yy
      COL <- pull.column(object,COL)
      XY <- cbind(XY,COL)

      COL <- ATTRIBUTE$COV.xy
      COL <- pull.column(object,COL)
      XY <- cbind(XY,COL)

      if(!is.null(XY) && ncol(XY)==3)
      {
        DATA$COV.x.x <- XY[,1]
        DATA$COV.y.y <- XY[,2]
        DATA$COV.x.y <- XY[,3]
        DATA$VAR.xy <- (DATA$COV.x.x+DATA$COV.y.y)/2

        DATA$COV.major <- DATA$COV.minor <- DATA$COV.angle <- NA
        for(z in levels(zone))
        {
          SUB <- zone==z
          DATA[SUB,] <- cov.xy2geo(DATA[SUB,],z)
        }

        ATLAS <- TRUE
      }
    }
    else
    { stop("Could not identify location columns.") }

    rm(zone,XY)
  } # end UTM
  else if(datum!="WGS84") # convert from input datum to default DATUM
  {
    ROWS <- is.na(DATA$longitude) | is.na(DATA$latitude)
    ROWS <- !ROWS

    FROM <- paste0("+proj=longlat +datum=",datum)

    XY <- DATA[ROWS,c('longitude','latitude')]
    XY <- project(XY,from=FROM)
    DATA[ROWS,c('longitude','latitude')] <- XY

    rm(ROWS,XY)
  }

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

  #################################
  # ESTABLISH BASE LOCATION CLASSES BEFORE MISSINGNESS CLASSES
  ###########################
  # generic location classes
  # includes Telonics Gen4 location classes (use with HDOP information)
  # unlike other data, don't use first choice but loop over all columns
  # retain ARGOS location classes if mixed, and any previous class
  for(COL in canonical.name(ATTRIBUTE$FIX))
  {
    PULL <- which(COL==canonical.name(names(object)))
    if(any(PULL))
    {
      COL <- names(object)[PULL[1]]
      CLASS <- as.factor(object[[COL]])
      DATA$class <- merge.class(CLASS,DATA$class)
    }
  }
  #COL <- pull.column(object,ATTRIBUTE$FIX,FUNC=as.factor)
  #if(length(COL)) { DATA$class <- merge.class(COL,DATA$class) }

  ERROR <- rep(10,length(DATA$t)) # default 10 meters GPS error
  CAL <- rep(0,length(DATA$t)) # default calibration degrees-of-freedom (none) for above estimate

  # VHF error ellipses provided
  if(ATLAS && all(DOP.LIST$horizontal$COV %in% names(DATA)))
  {
    NAS <- is.na(DATA[,DOP.LIST$horizontal$COV])
    NAS <- apply(NAS,1,any)
    ERROR[!NAS] <- TRUE
    CAL[!NAS] <- Inf
  }

  ##################################
  # ARGOS error ellipse/circle (newer ARGOS data >2011)
  COL <- ATTRIBUTE$COV.angle
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    DATA$COV.angle <- COL
    DATA$HDOP <- pull.column(object,"Argos.GDOP")

    # according to ARGOS, the following can be missing on <4 message data... but it seems present regardless
    DATA$COV.major <- pull.column(object,ATTRIBUTE$COV.major)^2/2
    DATA$COV.minor <- pull.column(object,ATTRIBUTE$COV.minor)^2/2

    if(DOP.LIST$horizontal$VAR %in% names(object))
    { DATA[[DOP.LIST$horizontal$VAR]] <- pull.column(object,ATTRIBUTE$COV.mean)^2/2 }
    else
    { DATA[[DOP.LIST$horizontal$VAR]] <- (DATA$COV.minor + DATA$COV.major)/2 }

    ARGOS <- TRUE
    NAS <- is.na(COL) # missing error ellipse rows
    ERROR[!NAS] <- TRUE # Argos KF calibrated
    CAL[!NAS] <- Inf
  }
  else  # no error ellipse information
  { NAS <- rep(TRUE,nrow(object)) }

  # ARGOS error categories (older ARGOS data <2011)
  # converted to error ellipses from ...
  COL <- c("Argos.location.class","Argos.lc","lc")
  COL <- pull.column(object,COL,as.factor)
  # ARGOS location classes present, but no error ellipses
  if(length(COL) && any(NAS) && any(c(3:0,'A','B','Z')%in%levels(COL)))
  {
    # major axis always longitude
    DATA$COV.angle <- 90
    # error eigen-variances
    ARGOS.minor <- c('3'=157,'2'=259,'1'=494, '0'=2271,'A'=762, 'B'=4596,'Z'=Inf)^2
    ARGOS.major <- c('3'=295,'2'=485,'1'=1021,'0'=3308,'A'=1244,'B'=7214,'Z'=Inf)^2
    ARGOS.DOF <- c('3'=25,'2'=51,'1'=55, '0'=60,'A'=103, 'B'=132,'Z'=0)
    # numbers from Vincent et al (2002)

    # error radii (average variance)
    ARGOS.radii <- (ARGOS.major+ARGOS.minor)/2

    message(sum(NAS)," ARGOS error ellipses missing. Using location class estimates from Vincent et al (2002).")
    COL <- as.character(COL) # factors are weird

    DATA$COV.minor[NAS] <- ARGOS.minor[COL][NAS]
    DATA$COV.major[NAS] <- ARGOS.major[COL][NAS]
    DATA[[DOP.LIST$horizontal$VAR]][NAS] <- ARGOS.radii[COL][NAS]
    DATA$HDOP[NAS] <- sqrt(2*ARGOS.radii)[COL][NAS]

    # classify Kalman filtered locations separately
    if(any(!NAS)) { COL[!NAS] <- "KF" }

    DATA$class <- as.factor(COL)

    # remove Z class (no calibration data)
    SUB <- (NAS & COL!="Z") | !NAS
    object <- object[SUB,]
    DATA <- DATA[SUB,]

    ARGOS <- TRUE

    COL <- NAS & !is.na(DATA$COV.major) # !KF and LC
    ERROR[COL] <- TRUE # Argos LC calibrated
    CAL[COL] <- Inf
  }

  # get dop values, but don't overwrite ARGOS GDOP if also present in ARGOS/GPS data
  try.dop <- function(DOPS,MESS=NULL,FN=identity,UERE=10)
  {
    if(ARGOS || !("HDOP" %in% names(DATA)))
    {
      COL <- pull.column(object,DOPS)
      if(length(COL))
      {
        # HDOPS to assign
        if(ARGOS) { NAS <- is.na(DATA$HDOP) }
        else { NAS <- rep(TRUE,length(COL)) }

        # don't overwrite ARGOS GDOPs
        if(any(NAS))
        {
          if(length(MESS)) { message(MESS) }
          DATA$HDOP[NAS] <<- FN(COL[NAS])
          ERROR[NAS] <<- UERE # assign error scale
          CAL[NAS] <<- 0 # assign degrees of freedom for above estimate
        }

        # don't make ARGOS exception again (2 DOP types)
        ARGOS <<- FALSE
      }
    }
  }
  # try to get dop values from best to worst
  try.dop(ATTRIBUTE$error,UERE=1)
  try.dop(ATTRIBUTE$HDOP)
  try.dop(ATTRIBUTE$DOP,"HDOP values not found. Using ambiguous DOP.")
  try.dop(ATTRIBUTE$PDOP,"HDOP values not found. Using PDOP.")
  try.dop(ATTRIBUTE$GDOP,"HDOP values not found. Using GDOP.")
  try.dop(ATTRIBUTE$Telonics,"HDOP values not found. Using Telonics error estimates.",UERE=1)
  try.dop(ATTRIBUTE$nsat,"HDOP values not found. Approximating via # satellites.",FN=function(x){(12-2)/(x-2)})

  # GPS-ARGOS hybrid data clean-up
  COL <- "sensor.type"
  COL <- pull.column(object,COL,FUNC=as.factor)
  if(length(COL))
  {
    levels(COL) <- tolower(levels(COL))
    LEVELS <- levels(COL)
    if(('gps' %in% LEVELS) && any(grepl('argos',LEVELS)))
    {
      GPS <- (COL=='gps')
      # missing GPS HDOP fix
      if(all(is.na(DATA$HDOP[GPS]))) { DATA$HDOP[GPS] <- 1 }
      # GPS errors are relatively small
      DATA$VAR.xy[GPS] <- 0
      DATA$COV.angle[GPS] <- 0
      DATA$COV.major[GPS] <- 0
      DATA$COV.minor[GPS] <- 0
      rm(GPS)
    }
    rm(LEVELS)
  }

  ###########
  # Telonics Gen4 GPS errors
  COL <- pull.column(object,ATTRIBUTE$Telonics)
  TELONICS <- length(COL)

  # detect if Telonics by location classes
  if(!TELONICS && "class" %in% names(DATA))
  {
    if(all(levels(DATA$class) %in% c("Succeeded (3D)","Succeeded (2D)","Resolved QFP","Resolved QFP (Uncertain)","Unresolved QFP","Failed")))
    { TELONICS <- TRUE }
  }

  # consolidate Telonics location classes as per manual
  if(TELONICS && "class" %in% names(DATA))
  {
    CLASS <- levels(DATA$class)
    SUB <- CLASS %in% c("Succeeded (3D)","Succeeded (2D)")
    if(any(SUB)) { CLASS[SUB] <- "Succeeded" }
    SUB <- CLASS %in% c("Resolved QFP","Resolved QFP (Uncertain)","Unresolved QFP")
    if(any(SUB)) { CLASS[SUB] <- "QFP" }
    levels(DATA$class) <- CLASS

    # some QFP locations can be missing HDOP, etc.
    DATA <- missing.class(DATA,"HDOP")
  }

  # account for missing DOP values
  if("HDOP" %in% names(DATA)) { DATA <- missing.class(DATA,"HDOP") }

  #######################
  # timed-out fixes
  if(timeout<Inf)
  {
    COL <- ATTRIBUTE$TTF
    COL <- pull.column(object,COL)
    if(length(COL))
    {
      if(class(timeout)[1]=="function") { timeout <- timeout(COL) }
      COL <- (COL<timeout)
      # factor levels are based on what's present and not what's possible
      COL <- c(TRUE,FALSE,COL)
      COL <- as.factor(COL)
      levels(COL) <- c("timeout","in-time")
      COL <- COL[-(1:2)]

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
  COL <- ATTRIBUTE$z
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    DATA$z <- COL

    COL <- ATTRIBUTE$VDOP
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

    # account for missing DOP values
    if("VDOP" %in% names(DATA)) { DATA <- missing.class(DATA,"VDOP") }
  }

  ########################################
  # VELOCITY
  # Import velocity information if present
  COL <- ATTRIBUTE$v
  COL <- pull.column(object,COL)
  if(length(COL))
  {
    DATA$speed <- COL
    COL <- ATTRIBUTE$heading
    COL <- pull.column(object,COL)
    if(length(COL))
    {
      DATA$heading <- COL

      ####################
      # SPEED ERE
      COL <- "eobs.speed.accuracy.estimate"
      COL <- pull.column(object,COL)
      if(length(COL) && FALSE) # e-obs column is terrible for error estimation, location error estimates are better # but they do not cross validate
      {
        DATA[[DOP.LIST$speed$DOP]] <- COL
      }
      else if("HDOP" %in% names(DATA)) # USE HDOP as approximate SDOP
      {
        DATA$SDOP <- DATA$HDOP
        # message("HDOP used as an approximate speed DOP.")
      }

      # do we need a location class for missing speeds?
      if("speed" %in% names(DATA))
      { DATA <- missing.class(DATA,"speed") }
    }
    else { DATA$speed <- NULL }
  } # END SPEED IMPORT

  #######################################

  # keep everything from original data.frame
  if(class(keep)[1]=="logical" && keep)
  { DATA <- cbind(DATA,object) }
  else if(class(keep)[1]=="character") # keep specified  columns
  { DATA <- cbind(DATA,object[,keep,drop=FALSE]) }

  # do this or possibly get empty animals from subset
  DATA <- droplevels(DATA)

  id <- levels(DATA$id)
  n <- length(id)

  # sort and clean individuals
  telist <- list()
  for(i in 1:n)
  {
    telist[[i]] <- DATA[DATA$id==id[i],]
    telist[[i]]$id <- NULL

    # clean through duplicates, sort times, etc..
    telist[[i]] <- telemetry.clean(telist[[i]],id=id[i],OCCURRENCE=OCCURRENCE)

    # can only check for hot/cold class after sorting times
    if(!is.na(dt.hot))
    {
      HOT <- diff(telist[[i]]$t) <= dt.hot
      HOT <- c(FALSE,HOT) # first fix is assumed to be cold
      HOT <- factor(HOT,labels=c("cold","hot"))

      if("class" %in% names(telist[[i]])) # combine with existing class information
      {
        HOT <- as.character(HOT)
        telist[[i]]$class <- paste(as.character(telist[[i]]$class),HOT)
        telist[[i]]$class <- as.factor(telist[[i]]$class)
      }
      else # only class information so far
      { telist[[i]]$class <- HOT }
    }

    # delete empty columns in case collars/tags are different
    for(COL in names(telist[[i]]))
    {
      if(COL %nin% keep)
      {
        BAD <- is.na(telist[[i]][[COL]])
        if(all(BAD) || (na.rm=="col" && any(BAD))) { telist[[i]][[COL]] <- NULL } # delete the whole column
        else if(na.rm=="row" && any(BAD)) { telist[[i]] <- telist[[i]][!BAD,] } # only delete the rows
      }
    }

    # delete incomplete vector data (caused by partial NAs)
    # i've seen speeds without headings...
    # COL <- c("speed","heading")
    # telist[[i]] <- rm.incomplete(telist[[i]],COL)
    # should rather return an error below?
    # COL <- c("COV.angle","COV.major","COV.minor")
    # telist[[i]] <- rm.incomplete(telist[[i]],COL)
  } # for(i in 1:n)

  # format error structure
  for(i in 1:n)
  {
    # drop missing levels and set UERE according to calibration
    telist[[i]] <- droplevels(telist[[i]])

    # NA UERE object given data.frame columns
    UERE <- uere.null(telist[[i]])
    DOF <- UERE$DOF
    UERE <- UERE$UERE

    # # which UERE values were specified
    # for(TYPE in colnames(UERE))
    # {
    #   # is data calibrated
    #   VAR <- DOP.LIST[[TYPE]]$VAR
    #   if(VAR %in% names(telist[[i]]))
    #   {
    #     VAR <- (telist[[i]][[VAR]] < Inf) # is calibrated
    #
    #     for(CLASS in levels(telist[[i]]$class))
    #     { if(all(VAR[telist[[i]]$class==CLASS])) { UERE[CLASS,TYPE] <- 1 } }
    #
    #     if(!nlevels(telist[[i]]$class) && all(VAR)) { UERE["all",TYPE] <- 1 }
    #   }
    # }

    # assign UERE and DOF estimates
    CLASSES <- levels(telist[[i]]$class)
    if(length(CLASSES))
    {
      CLASSES <- sort(CLASSES) # sort levels alphabetically (canonical)
      telist[[i]]$class <- factor(telist[[i]]$class,levels=CLASSES)

      for(CLASS in CLASSES)
      {
        COL <- telist[[i]]$class==CLASS
        UERE[CLASS,] <- median(ERROR[COL]) # these numbers should all be the same
        DOF[CLASS,'horizontal'] <- median(CAL[COL]) # these numbers should all be the same
      }

      # consistent ordering
      UERE <- UERE[CLASSES,,drop=FALSE]
      DOF <- DOF[CLASSES,,drop=FALSE]
    }
    else
    {
      UERE["all",] <- median(ERROR) # these numbers should all be the same
      DOF["all",'horizontal'] <- median(CAL) # these numbers should all be the same
    }

    # combine data.frame with ancillary info
    info <- list(identity=id[i], timezone=timezone, projection=projection)
    INF <- pmax(UERE[1,],Inf)
    names(INF) <- colnames(UERE) # R drops dimnames...
    N <- DOF[1,]
    names(N) <- colnames(UERE) # R drops dimnames...
    UERE <- list(UERE=UERE,DOF=DOF,AICc=INF,Zsq=INF,VAR.Zsq=INF,N=N)
    UERE <- new.UERE(UERE)
    telist[[i]] <- new.telemetry( telist[[i]] , info=info, UERE=UERE )
  } # for(i in 1:n)
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
telemetry.clean <- function(data,id,OCCURRENCE=FALSE)
{
  # test for out of order (or out of reverse order)
  if(nrow(data)>1)
  {
    DIFF <- diff(data$t) # all positive if in forward order
    DIFF <- DIFF * DIFF[DIFF!=0][1] # reverse order becomes all positive too
    if(!OCCURRENCE) # tracking data shouldn't have repeating times
    {
      DIFF <- length(DIFF) && all(DIFF>0)
      if(!DIFF) { warning("Times might be out of order or duplicated in ",id,". Make sure that timeformat and timezone are correctly specified.") }
    }
    else # occurrence data can have repeating times, but they should still be in order
    {
      DIFF <- length(DIFF) && all(DIFF>=0)
      if(!DIFF) { warning("Times might be out of order in ",id,". Make sure that timeformat and timezone are correctly specified.") }
    }
  }

  # sort in time
  ORDER <- sort.list(data$t,na.last=NA,method="quick")
  data <- data[ORDER,]

  if(!OCCURRENCE)
  {
    # remove duplicate observations
    ORDER <- length(data$t)
    data <- unique(data)
    if(ORDER != length(data$t)) { warning("Duplicate data in ",id," removed.") }

    # exit with warning on duplicate times
    if(anyDuplicated(data$t)) { warning("Duplicate times in ",id,". Data cannot be fit without an error model.") }
  }

  # remove old level information
  data <- droplevels(data)

  if(!OCCURRENCE)
  {
    dt <- diff(data$t)
    # v <- sqrt(diff(data$x)^2+diff(data$y)^2)/dt
    # v <- max(v)
    # message("Maximum speed of ",format(v,digits=3)," m/s observed in ",id)
    dt <- min(dt,Inf)
    if(dt<Inf)
    {
      units <- unit(dt,dimension='time')
      message("Minimum sampling interval of ",format(dt/units$scale,digits=3)," ",units$name," in ",id)
    }
    else
    { message("Minimum sampling interval of NA in ",id) }
  }

  return(data)
}


########################
# summarize telemetry data
summary.telemetry <- function(object,...)
{
  if(class(object)[1]=="telemetry")
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
  else if(class(object)[1]=="list")
  {
    # NAME <- sapply(object,function(o){ attr(o,"info")$identity })
    DT <- sapply(object,function(o){ stats::median(diff(o$t)) })
    P <- sapply(object,function(o){ last(o$t) - o$t[1] })
    lon <- sapply(object, function(o){ stats::median(o$longitude) })
    lat <- sapply(object, function(o){ stats::median(o$latitude) })

    result <- data.frame(interval=DT,period=P,longitude=lon,latitude=lat)

    # unit conversions
    COL <- 1
    units <- unit(stats::median(DT,na.rm=TRUE),"time",thresh=1,concise=TRUE)
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


check.class <- function(data,target=c("telemetry","data.frame"))
{
  CLASS <- class(data)[1]
  if(CLASS %nin% target)
  {
    if(CLASS=="list")
    {
      CLASS.1 <- class(data[[1]])[1]
      if(CLASS.1  %in% target)
      {
        n <- length(data)
        stop("Argument is a list of ",n," objects, and not a single class ",target[1]," object")
      }
    }

    stop("Argument is of class ",CLASS,", and not class ",target[1])
  }
}

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
