# some default and non-default units
UNITS <- list()
UNITS$mean <- list()
UNITS$mean$DAY <- 86400.002 # mean solar day
UNITS$mean$YEAR <- 365.24217 * UNITS$mean$DAY # mean tropical year
UNITS$mean$MONTH <- 2.8 + 60*( 44 + 60*( 12 + 24*29 ) ) # mean synodic month
UNITS$calendar <- list()
UNITS$calendar$DAY <- 86400
UNITS$calendar$YEAR <- 365 * UNITS$calendar$DAY
UNITS$calendar$MONTH <- 30 * UNITS$calendar$DAY

# create units dictionary
generate.units <- function()
{
  alias <- list()
  scale <- c()

  add <- function(a,s)
  {
    n <- length(alias)
    alias[[n+1]] <<- canonical.name(a)
    scale[n+1] <<- s
  }

  OP <- getOption("time.units")

  # TIME
  add(c("\u03BCs","microsecond","microseconds"),1E-6)
  add(c("ms","milisecond","miliseconds"),1/1000)
  add(c("s","sec","sec.","second","seconds"),1)
  add(c("min","minute","minutes"),60)
  add(c("h","hr","hour","hours"),60^2)
  add(c("day","days"),UNITS[[OP]]$DAY) #day
  add(c("wk","week","weeks"),7*UNITS[[OP]]$DAY) # week
  add(c("mon","month","months"),UNITS[[OP]]$MONTH) # month
  add(c("yr","year","years"),UNITS[[OP]]$YEAR) # year
  add(c("ka","ky","millennium","millenniums","millennia","kiloannum","kiloannums","kiloyear","kiloyears"),1000*UNITS[[OP]]$YEAR)
  add(c("ma","my","megaannum","megaannums","megaanna","megayear","megayears","millionennium","millionenniums","millionennia"),1000^2*UNITS[[OP]]$YEAR)
  add(c("ae","ga","gy","gyr","aeon","aeons","eon","eons","gigayear","gigayears","gigaannum","gigaannums","giggaanna"),1000^3*UNITS[[OP]]$YEAR)

  # Distance conversions
  add(c("\u03BCm","micron","microns","micrometer","micrometers"),1E-6)
  add(c("mm","milimeter","milimeters"),1/1000)
  add(c("cm","centimeter","centimeters"),1/100)
  add(c("m","meter","meters"),1)
  add(c("km","kilometer","kilometers"),1000)
  add(c("in","inch","inches"),0.3048/12)
  add(c("ft","foot","feet"),0.3048)
  add(c("yd","yard","yards"),0.3048*3)
  add(c("mi","mile","miles"),0.3048*5280)

  # Area conversions
  add(c("\u03BCm\u00B2","micron\u00B2","microns\u00B2","micrometer\u00B2","micrometers\u00B2","\u03BCm^2","\u03BCm.^2","micron^2","microns^2","micrometer^2","micrometers^2","square micron","square microns","square micrometer","square micrometers","micron squared","microns squared","micrometer squared","micrometers squared"),1E-12)
  add(c("mm\u00B2","milimeter\u00B2","milimeters\u00B2","mm^2","mm.^2","milimeter^2","milimeters^2","square milimeter","square milimeters","milimeter squared","milimeters squared"),1/1000^2)
  add(c("cm\u00B2","centimeter\u00B2","centimeters\u00B2","cm^2","cm.^2","centimeter^2","centimeters^2","square centimeter","square centimeters","centimeter squared","centimeters squared"),1/100^2)
  add(c("m\u00B2","meter\u00B2","meters\u00B2","m^2","m.^2","meter^2","meters^2","square meter","square meters","meter squared","meters squared"),1)
  add(c("ha","hectare","hectares","hm\u00B2","hectometer\u00B2","hectometre\u00B2","hectometers\u00B2","hectometres\u00B2","hm^2","hectometer^2","hectometre^2","hectometers^2","hectometres^2","square hm","square hectometer","square hectometre","square hectometers","square hectometres"),100^2)
  add(c("km\u00B2","kilometer\u00B2","kilometers\u00B2","km^2","km.^2","kilometer^2","kilometers^2","square kilometer","square kilometers","kilometer squared","kilometers squared"),1000^2)
  add(c("in\u00B2","inch\u00B2","inches\u00B2","in^2","in.^2","inch^2","inches^2"),(0.3048/12)^2)
  add(c("ft\u00B2","foot\u00B2","feet\u00B2","ft^2","ft.^2","foot^2","feet^2","square foot","square feet","foot squared","feet squared"),0.3048^2)
  add(c("yd\u00B2","yard\u00B2","yards\u00B2","yd^2","yd.^2","yard^2","yards^2","square yard","square yards","yard squared","yards squared"),(0.3048*3)^2)
  add(c("mi\u00B2","mile\u00B2","miles\u00B2","mi^2","mi.^2","mile^2","miles^2","square mile","square miles","mile squared","miles squared"),(0.3048*5280)^2)

  # speed
  add(c("mps","m/s","m/sec","meter/sec","meter/second","meters/second"),1)
  add(c("kmph","kph","km/h","km/h","km/hr","kilometer/hour","kilometers/hour"),0.277777777777777777777)
  add(c("mph","mi/h","mi/hr","mile/h","mile/hr","mile/hour","miles/hour"),0.44704)
  add(c("fps","ft/s","ft/sec","feet/second"),0.3048)
  add(c('kt','kn','knot','knots'),1.852 * 0.277777777777777777777)
  add(c("km/s","kmps","km/sec"),1000)
  add(c("cm/s","cmps","cm/sec"),1/100)
  add(c("mm/s","mmps","mm/sec"),1/1000)
  add(c("\u03BCm/s","\u03BCmps","\u03BCm/sec"),1E-6)

  # frequency
  add(c("Hz","hertz"),1)
  add(c("kHz","kilohertz"),1000)
  add(c("MHz","megahertz"),1000^2)
  add(c("GHz","gigahertz"),1000^3)
  add(c("THz","terahertz"),1000^4)
  add(c("per min","1/min","min^-1","min\u207B\u00B9","per minute","1/minute","minute^-1","minute\u207B\u00B9"),1/60)
  add(c("per hr","1/hr","hr^-1","hr\u207B\u00B9","per hour","1/hour","hour^-1","hour\u207B\u00B9"),1/60^2)
  add(c("per day","1/day","day^-1","day\u207B\u00B9"),1/UNITS[[OP]]$DAY)
  add(c("per mon","1/mon","mon^-1","mon\u207B\u00B9","per month","1/month","month^-1","month\u207B\u00B9"),1/UNITS[[OP]]$MONTH)
  add(c("per yr","1/yr","yr^-1","yr\u207B\u00B9","per year","1/year","year^-1","year\u207B\u00B9"),1/UNITS[[OP]]$YEAR)
  add(c("per ky","1/ky","ky^-1","ky\u207B\u00B9","per ka","1/ka","ka^-1","ka\u207B\u00B9","per millennium","1/millennium","millennium^-1","millennium\u207B\u00B9","per kiloannum","1/kiloannum","kiloannum^-1","kiloannum\u207B\u00B9","per kiloyear","1/kiloyear","kiloyear^-1","kiloyear\u207B\u00B9"),1/UNITS[[OP]]$YEAR/1000)
  add(c("per my","1/my","my^-1","my\u207B\u00B9","per ma","1/ma","ma^-1","ma\u207B\u00B9","per megaannum","1/megaannum","megaannum^-1","megaannum\u207B\u00B9","per megayear","1/megayear","megayear^-1","megayear\u207B\u00B9","per millionennium","1/millionennium","millionennium^-1","millionennium\u207B\u00B9"),1/UNITS[[OP]]$YEAR/1000^2)
  add(c("per ae","1/ae","ae^-1","ae\u207B\u00B9","per ga","1/ga","ga^-1","ga\u207B\u00B9","per gy","1/gy","gy^-1","gy\u207B\u00B9","per gyr","1/gyr","gyr^-1","gyr\u207B\u00B9","per aeon","1/aeon","aeon^-1","aeon\u207B\u00B9","per eon","1/eon","eon^-1","eon\u207B\u00B9","per gigayear","1/gigayear","gigayear^-1","gigayear\u207B\u00B9","per gigaannum","1/gigaannum","gigaannum^-1","gigaannum\u207B\u00B9"),1/UNITS[[OP]]$YEAR/1000^3)

  # mass
  add(c("g","gm","gram","grams"),1/1000) # kg is SI
  add(c("kg","kilogram","kilograms"),1) # kg is SI
  add(c("Mg","t","tonne","tonnes","mt","m ton","m tons","metric ton","metric tons"),1000)
  add(c("mg","milligram","milligrams"),1/1000)
  add(c("\u03BCg","microgram","micrograms"),1/1000^2)
  add(c("ng","nanogram","nanograms"),1/1000^3)
  add(c("lb","lbs","pound","pounds"),0.45359237)
  add(c("oz","ounce","ounces"),0.45359237/16)
  add(c("slug","slugs"),14.59390)
  add(c("st","stone","stones"),0.45359237*14)
  add(c("ton","tons"),0.45359237*2000) # NA ton (not UK)

  # memory
  add(c("byte","bytes","B"),1)
  add(c("Kb","KiB"),1024)
  add(c("Mb","MiB"),1024^2)
  add(c("Gb","GiB"),1024^3)
  add(c("Tb","TiB"),1024^4)
  add(c("Pb","PiB"),1024^5)
  add(c("Eb","EiB"),1024^6)
  add(c("Zb","ZiB"),1024^7)
  add(c("Yb","YiB"),1024^8)

  return(list(alias=alias,scale=scale))
}
UNIT <- list() # generated onLoad
#UNIT <- generate.units()


dimfig <- function(data,dimension,thresh=1,...)
{
  UNITS <- unit(data,dimension=dimension,thresh=thresh,...)
  data <- data/UNITS$scale
  R <- list(data=data,units=c(UNITS$name,UNITS$abrv))
  return(R)
}

# CHOOSE BEST UNITS FOR A LIST OF DATA
# thresh is threshold for switching to coarser unit
# concise gives abbreviated names
unit <- function(data,dimension,thresh=1,concise=FALSE,SI=FALSE,...)
{
  if(SI) { data <- 1.01 ; thresh <- 1 } # will always choose base units
  OP <- getOption("time.units")

  if(dimension %in% c("length",'distance'))
  {
    name.list <- c("microns","milimeters","centimeters","meters","kilometers")
    abrv.list <- c("\u03BCm","mm","cm","m","km")
    scale.list <- c(1E-6,1/1000,1/100,1,1000)
  }
  else if(dimension=="area")
  {
    name.list <- c("square microns","square milimeters","square centimeters","square meters","hectares","square kilometers")
    abrv.list <- c("\u03BCm\u00B2","mm\u00B2","cm\u00B2","m\u00B2","hm\u00B2","km\u00B2")
    scale.list <- c(1E-12,1/1000^2,1/100^2,1,100^2,1000^2)
  }
  else if(dimension=="time")
  {
    name.list <- c("microseconds","miliseconds","seconds","minutes","hours","days","months","years","millennia","mega-anna","aeons")
    abrv.list <- c("\u03BCs","ms","sec","min","hr","day","mon","yr","ka","Ma","AE")
    scale.list <- c(1E-6,1/1000,1,60*c(1,60)) # through minutes
    scale.list[6] <- UNITS[[OP]]$DAY
    scale.list[7] <- UNITS[[OP]]$MONTH
    scale.list[8] <- UNITS[[OP]]$YEAR
    scale.list[9] <- 1000 * UNITS[[OP]]$YEAR
    scale.list[10] <- 1000^2 * UNITS[[OP]]$YEAR
    scale.list[11] <- 1000^3 * UNITS[[OP]]$YEAR
  }
  else if(dimension %in% c("speed",'velocity'))
  {
    name.list <- c("microns/day","milimeters/day","centimeters/day","meters/day","kilometers/day")
    abrv.list <- c("\u03BCm/day","mm/day","cm/day","m/day","km/day")
    scale.list <- c(1E-6,1/1000,1/100,1,1000)/UNITS[[OP]]$DAY

    # SI units fix
    if(SI)
    {
      name.list <- "meters/second"
      abrv.list <- "m/s"
      scale.list <- 1
    }
  }
  else if(dimension=="diffusion")
  {
    name.list <- c("square microns/day","square milimeters/day","square centimeters/day","square meters/day","hectares/day","square kilometers/day")
    abrv.list <- c("\u03BCm\u00B2/day","mm\u00B2/day","cm\u00B2/day","m\u00B2/day","hm\u00B2/day","km\u00B2/day")
    scale.list <- c(1E-12,1/1000^2,1/100^2,1,100^2,1000^2)/UNITS[[OP]]$DAY

    # SI units fix
    if(SI)
    {
      name.list <- "square meters/second"
      abrv.list <- "m\u00B2/s"
      scale.list <- 1
    }
  }
  else if(dimension=="frequency")
  {
    name.list <- c("per aeon","per mega-anna","per millennia","per year","per month","per day","per hour","per minute","hertz","kilohertz","megahertz","gigahertz","terahertz")
    abrv.list <- c("AE\u207B\u00B9","Ma\u207B\u00B9","ka\u207B\u00B9","yr\u207B\u00B9","mon\u207B\u00B9","day\u207B\u00B9","hr\u207B\u00B9","min\u207B\u00B9","Hz","kHz","MHz","GHz","THz")
    scale.list <- c(1/UNITS[[OP]]$YEAR/1000^3,1/UNITS[[OP]]$YEAR/1000^2,1/UNITS[[OP]]$YEAR/1000,1/UNITS[[OP]]$YEAR,1/UNITS[[OP]]$MONTH,1/UNITS[[OP]]$DAY,1/60^2,1/60,1,1000,1000^3,1000^6,1000^9,1000^12)
  }
  else if(dimension=="mass")
  {
    name.list <- c("nanograms","micrograms","milligrams","grams","kilograms","tonnes")
    abrv.list <- c("ng","\u03BCg","mg","gm","kg","Mg")
    scale.list <- c(1/1000^4,1/1000^3,1/1000^2,1/1000,1,1000)
  }
  else # units not recognized
  {
    R <- list(scale=1,name=NULL)
    return(R)
  }

  data <- data[!is.na(data)]
  data <- data[abs(data)<Inf]
  if(length(data)) { max.data <- max(abs(data)) } else { max.data <- 1 }

  if(concise) { name.list <- abrv.list }

  # choose most parsimonious units
  I <- (max.data >= thresh * scale.list)
  if(any(I))
  {
    I <- (1:length(I))[I]
    I <- last(I)
  }
  else { I <- 1 }

  name <- name.list[I]
  abrv <- abrv.list[I]
  scale <- scale.list[I]

  return(list(scale=scale,name=name,abrv=abrv))
}


### determine parsimonious units for parameter CIs
# preference point estimate, but fall back on high-CI if point estimate is zero
unit.par <- function(par,...)
{
  PAR <- par[2:3]
  PAR <- PAR[PAR>.Machine$double.eps]
  if(length(PAR)) { PAR <- min(PAR) } else { PAR <- 0 }

  return( unit(PAR,...) )
}


## rescale the units of telemetry object
unit.telemetry <- function(data,length=1,time=1,axes=c('x','y'))
{
  if(class(data)[1]=="phylometry")
  {
    lag <- attr(data,"lag")/time
    data[,axes] <- data[,axes]/length
    attr(data,"lag") <- lag
    # class(data) <- "phylometry"
    return(data)
  }
  # TELEMETRY CLASS BELOW

  convert <- function(NAMES,scale) { for(NAME in NAMES) { if(NAME %in% names(data)) { data[[NAME]] <<- data[[NAME]]/scale } } }

  convert('t',time)
  convert('light.time',time)
  convert('dark.time',time)
  convert('sundial.rate',1/time)

  if(any(axes %in% c('x','y','z')))
  {
    convert(DOP.LIST$horizontal$axes,length)
    convert(DOP.LIST$vertical$axes,length)
    convert(DOP.LIST$speed$axes,length/time)

    convert(DOP.LIST$horizontal$VAR,length^2)
    convert(DOP.LIST$horizontal$COV,length^2)

    convert(DOP.LIST$vertical$VAR,length^2)

    convert(DOP.LIST$speed$VAR,(length/time)^2)
    convert(DOP.LIST$speed$COV,(length/time)^2)
  }
  else
  { for(axis in axes) { convert(axis,length) } }

  # calibration constants
  attr(data,"UERE")$UERE <- attr(data,"UERE")$UERE/length # don't logicals get divided this way?

  return(data)
}


## rescale the units of dimensionful parameters
unit.ctmm <- function(CTMM,length=1,time=1)
{
  if(length(CTMM$tau)){ CTMM$tau <- CTMM$tau/time }
  CTMM$omega <- CTMM$omega * time
  CTMM$circle <- CTMM$circle * time

  # all means scale with length the same way... but not time
  if("mu" %in% names(CTMM))
  {
    CTMM$mu <- CTMM$mu/length
    CTMM <- drift.scale(CTMM,time)
  }

  if(!is.null(CTMM$timelink.cycle))
  { CTMM$timelink.cycle <- CTMM$timelink.cycle/time }

  # if(class(CTMM$error)[1]=='numeric')
  { CTMM$error <- CTMM$error/length } # don't divide logicals

  if("sigma" %in% names(CTMM))
  {
    CTMM$sigma <- scale.covm(CTMM$sigma,1/length^2)

    # variance -> diffusion adjustment
    if(!CTMM$range)
    { CTMM$sigma <- scale.covm(CTMM$sigma,time) }
  }

  if("COV.mu" %in% names(CTMM)) { CTMM$COV.mu <- CTMM$COV.mu/length^2 }

  if("COV" %in% names(CTMM))
  {
    NAMES <- dimnames(CTMM$COV)[[1]]

    PAR <- c('major','minor')
    PAR <- PAR[PAR %in% NAMES]
    if(length(PAR))
    {
      CTMM$COV[PAR,] <- CTMM$COV[PAR,]/length^2
      CTMM$COV[,PAR] <- CTMM$COV[,PAR]/length^2

      if(!CTMM$range)
      {
        CTMM$COV[PAR,] <- CTMM$COV[PAR,]*time
        CTMM$COV[,PAR] <- CTMM$COV[,PAR]*time
      }
    }

    tau <- CTMM$tau
    tau <- tau[tau<Inf]
    if(length(tau))
    {
      tau <- NAMES[grepl("tau",NAMES)]
      if(length(tau))
      {
        CTMM$COV[tau,] <- CTMM$COV[tau,]/time
        CTMM$COV[,tau] <- CTMM$COV[,tau]/time
      }

      if("omega" %in% NAMES)
      {
        CTMM$COV["omega",] <- CTMM$COV["omega",] * time
        CTMM$COV[,"omega"] <- CTMM$COV[,"omega"] * time
      }
    }

    ERROR <- NAMES[grepl("error",NAMES)] # error estimate covariance
    if(length(ERROR))
    {
      CTMM$COV[ERROR,] <- CTMM$COV[ERROR,]/length
      CTMM$COV[,ERROR] <- CTMM$COV[,ERROR]/length
    }

    if("circle" %in% NAMES)
    {
      CTMM$COV["circle",] <- CTMM$COV["circle",] * time
      CTMM$COV[,"circle"] <- CTMM$COV[,"circle"] * time
    }
  }

  return(CTMM)
}


######################
unit.UD <- function(UD,length=1)
{
  UD$r <- lapply(UD$r,function(x){ x/length })
  UD$PDF <- UD$PDF * length^2
  UD$dr <- UD$dr / length
  UD$H <- UD$H / length^2

  return(UD)
}


##################
unit.variogram <- function(SVF,time=1,area=1)
{
  SVF$lag <- SVF$lag / time
  SVF$SVF <- SVF$SVF / area
  if("MSE" %in% names(SVF)) { SVF$MSE <- SVF$MSE / area }

  return(SVF)
}


# convert units
`%#%` <- function(x,y)
{
  # convert to si units
  if(is.numeric(x))
  {
    num <- x
    name <- y
    pow <- +1
  }
  else # convert from si units
  {
    if(!is.numeric(y)) { return(x %#% 1 %#% y) }

    num <- y
    name <- x
    pow <- -1
  }

  name <- canonical.name(name)
  if(name=="") { return(num) }

  name <- strsplit(name,'*',fixed=TRUE)[[1]]
  if(length(name)>1)
  {
    if(pow==1)
    { for(i in 1:length(name)) { num <- num %#% name[i] } }
    else if(pow==-1)
    { for(i in 1:length(name)) { num <- name[i] %#% num } }
    return(num)
  }

  name <- strsplit(name,"/",fixed=TRUE)[[1]]
  if(length(name)>1)
  {
    if(pow==1)
    {
      num <- num %#% name[1]
      for(i in 2:length(name)) { num <- name[i] %#% num }
    }
    else if(pow==-1)
    {
      num <- name[1] %#% num
      for(i in 2:length(name)) { num <- num %#% name[i] }
    }
    return(num)
  }

  alias <- UNIT$alias
  scale <- UNIT$scale
  for(i in 1:length(alias))
  {
    if(name %in% alias[[i]]) { return(num*scale[i]^pow) }
  }
  stop(paste("Unit",name,"unknown."))
}

# interpret a string as number with units
ustring <- function(x)
{
  x <- canonical.name(x)

  y <- strsplit(x,"[a-z,A-Z]")[[1]][1]
  n <- nchar(y)+1
  y <- as.numeric(y)
  x <- substr(x,n,nchar(x))

  x <- y %#% x
  return(x)
}
