UNITS <- list()
UNITS$mean <- list()
UNITS$mean$DAY <- 86400.002 # mean solar day
UNITS$mean$YEAR <- 365.24217 * UNITS$mean$DAY # mean tropical year
UNITS$mean$MONTH <- 2.8 + 60*( 44 + 60*( 12 + 24*29 ) ) # mean synodic month
UNITS$calendar <- list()
UNITS$calendar$DAY <- 86400
UNITS$calendar$YEAR <- 365 * UNITS$calendar$DAY
UNITS$calendar$MONTH <- 30 * UNITS$calendar$DAY

# CHOOSE BEST UNITS FOR A LIST OF DATA
# thresh is threshold for switching to coarser unit
# concise gives abbreviated names
unit <- function(data,dimension,thresh=1,concise=FALSE,SI=FALSE)
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
    name.list <- c("microseconds","miliseconds","seconds","minutes","hours","days","months","years")
    abrv.list <- c("\u03BCs","ms","sec","min","hr","day","mon","yr")
    scale.list <- c(1E-6,1/1000,1,60*c(1,60)) # through minutes
    scale.list[6] <- UNITS[[OP]]$DAY
    scale.list[7] <- UNITS[[OP]]$MONTH
    scale.list[8] <- UNITS[[OP]]$YEAR
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

  data <- data[!is.na(data)]
  max.data <- max(abs(data))

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
  scale <- scale.list[I]

  return(list(scale=scale,name=name))
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
unit.telemetry <- function(data,length=1,time=1)
{
  convert <- function(NAMES,scale) { for(NAME in NAMES) { if(NAME %in% names(data)) { data[[NAME]] <<- data[[NAME]]/scale } } }

  convert(DOP.LIST$horizontal$axes,length)
  convert('t',time)
  convert(DOP.LIST$speed$axes,length/time)

  convert(DOP.LIST$horizontal$VAR,length^2)
  convert(DOP.LIST$horizontal$COV,length^2)

  convert(DOP.LIST$speed$VAR,(length/time)^2)
  convert(DOP.LIST$speed$COV,(length/time)^2)

  # HDOP is unitless

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
    drift <- get(CTMM$mean)
    CTMM <- drift@scale(CTMM,time)
  }

  if(class(CTMM$error)=='numeric') { CTMM$error <- CTMM$error/length } # don't divide logicals

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

    if("major" %in% NAMES)
    {
      CTMM$COV["major",] <- CTMM$COV["major",]/length^2
      CTMM$COV[,"major"] <- CTMM$COV[,"major"]/length^2

      if(!CTMM$range)
      {
        CTMM$COV["major",] <- CTMM$COV["major",]*time
        CTMM$COV[,"major"] <- CTMM$COV[,"major"]*time
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

    if("error" %in% NAMES)
    {
      CTMM$COV["error",] <- CTMM$COV["error",]/length
      CTMM$COV[,"error"] <- CTMM$COV[,"error"]/length
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
    num <- y
    name <- x
    pow <- -1
  }

  name <- tolower(name)

  # DIV <- grepl("/",name)
  # if(DIV)
  # {
  #   #
  # }

  OP <- getOption("time.units")

  alias <- list()
  scale <- c()

  add <- function(a,s)
  {
    n <- length(alias)
    alias[[n+1]] <<- a
    scale[n+1] <<- s
  }

  # TIME
  add(c("\u03BCs","\u03BCs.","microsecond","microseconds"),1E-6)
  add(c("ms","ms.","milisecond","miliseconds"),1/1000)
  add(c("s","s.","sec","sec.","second","seconds"),1)
  add(c("min","min.","minute","minutes"),60)
  add(c("h","h.","hr","hr.","hour","hours"),60^2)
  add(c("day","days"),UNITS[[OP]]$DAY) #day
  add(c("wk","wk.","week","weeks"),7*UNITS[[OP]]$DAY) # week
  add(c("mon","mon.","month","months"),UNITS[[OP]]$MONTH) # month
  add(c("yr","yr.","year","years"),UNITS[[OP]]$YEAR) # year

  # Distance conversions
  add(c("\u03BCm","\u03BCm.","micron","microns","micrometer","micrometers"),1E-6)
  add(c("mm","mm.","milimeter","milimeters"),1/1000)
  add(c("cm","cm.","centimeter","centimeters"),1/100)
  add(c("m","m.","meter","meters"),1)
  add(c("km","km.","kilometer","kilometers"),1000)
  add(c("in","in.","inch","inches"),0.3048/12)
  add(c("ft","ft.","foot","feet"),0.3048)
  add(c("yd","yd.","yard","yards"),0.3048*3)
  add(c("mi","mi.","mile","miles"),0.3048*5280)

  # Area conversions
  add(c("\u03BCm\u00B2","\u03BCm.\u00B2","micron\u00B2","microns\u00B2","micrometer\u00B2","micrometers\u00B2","\u03BCm^2","\u03BCm.^2","micron^2","microns^2","micrometer^2","micrometers^2","square micron","square microns","square micrometer","square micrometers","micron squared","microns squared","micrometer squared","micrometers squared"),1E-12)
  add(c("mm\u00B2","mm.\u00B2","milimeter\u00B2","milimeters\u00B2","mm^2","mm.^2","milimeter^2","milimeters^2","square milimeter","square milimeters","milimeter squared","milimeters squared"),1/1000^2)
  add(c("cm\u00B2","cm.\u00B2","centimeter\u00B2","centimeters\u00B2","cm^2","cm.^2","centimeter^2","centimeters^2","square centimeter","square centimeters","centimeter squared","centimeters squared"),1/100^2)
  add(c("m\u00B2","m.\u00B2","meter\u00B2","meters\u00B2","m^2","m.^2","meter^2","meters^2","square meter","square meters","meter squared","meters squared"),1)
  add(c("ha","hectare","hectares","hm\u00B2","hectometer\u00B2","hectometre\u00B2","hectometers\u00B2","hectometres\u00B2","hm^2","hectometer^2","hectometre^2","hectometers^2","hectometres^2","square hm","square hectometer","square hectometre","square hectometers","square hectometres"),100^2)
  add(c("km\u00B2","km.\u00B2","kilometer\u00B2","kilometers\u00B2","km^2","km.^2","kilometer^2","kilometers^2","square kilometer","square kilometers","kilometer squared","kilometers squared"),1000^2)
  add(c("in\u00B2","in.\u00B2","inch\u00B2","inches\u00B2","in^2","in.^2","inch^2","inches^2"),(0.3048/12)^2)
  add(c("ft\u00B2","ft.\u00B2","foot\u00B2","feet\u00B2","ft^2","ft.^2","foot^2","feet^2","square foot","square feet","foot squared","feet squared"),0.3048^2)
  add(c("yd\u00B2","yd.\u00B2","yard\u00B2","yards\u00B2","yd^2","yd.^2","yard^2","yards^2","square yard","square yards","yard squared","yards squared"),(0.3048*3)^2)
  add(c("mi\u00B2","mi.\u00B2","mile\u00B2","miles\u00B2","mi^2","mi.^2","mile^2","miles^2","square mile","square miles","mile squared","miles squared"),(0.3048*5280)^2)

  # speed
  add(c("mps","m/s","m/sec","meter/sec","meter/second","meters/second"),1)
  add(c("kmph","kph","km/h","km/h","km/hr","kilometer/hour","kilometers/hour"),0.277777777777777777777)
  add(c("mph","mi/h","mi/hr","mile/h","mile/hr","mile/hour","miles/hour"),0.44704)
  add(c("fps","ft/s","ft/sec","feet/second"),0.3048)
  add(c('kt','kn','knot','knots'),1.852 * 0.277777777777777777777)

  for(i in 1:length(alias))
  {
    if(name %in% alias[[i]]) { return(num*scale[i]^pow) }
  }
  stop(paste("Unit",name,"unknown."))
}
