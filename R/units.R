# CHOOSE BEST UNITS FOR A LIST OF DATA
# thresh is threshold for switching to coarser unit
# concise gives abbreviated names
unit <- function(data,dimension,thresh=1,concise=FALSE)
{
  if(dimension=="length")
  {
    name.list <- c("meters","kilometers")
    abrv.list <- c("m","km")
    scale.list <- c(1,1000)
  }
  else if(dimension=="area")
  {
    name.list <- c("square meters","hectares","square kilometers")
    abrv.list <- c("m\u00B2","hm\u00B2","km\u00B2")
    scale.list <- c(1,100^2,1000^2)
  }
  else if(dimension=="time")
  {
    name.list <- c("seconds","minutes","hours","days","months","years")
    abrv.list <- c("sec","min","hr","day","mon","yr")
    scale.list <- c(1,60*c(1,60*c(1,24*c(1,29.53059,365.24))))
  }
  else if(dimension=="speed")
  {
    name.list <- c("meters/day","kilometers/day")
    abrv.list <- c("m/day","km/day")
    scale.list <- c(1,1000)/(60*60*24)
  }
  else if(dimension=="diffusion")
  {
    name.list <- c("square meters/day","hectares/day","square kilometers/day")
    abrv.list <- c("m\u00B2/day","hm\u00B2/day","km\u00B2/day")
    scale.list <- c(1,100^2,1000^2)/(60*60*24)
  }

  max.data <- max(abs(data))

  if(concise) { name.list <- abrv.list }

  # choose most parsimonious units
  I <- max.data > thresh * scale.list
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


## rescale the units of telemetry object
unit.telemetry <- function(data,length=1,time=1)
{
  data$x <- data$x/length
  data$y <- data$y/length
  data$t <- data$t/length

  # error in meters
  if("HERE" %in% names(data)) { data$HERE <- data$HERE/length }
  # HDOP is unitless

  return(data)
}


## rescale the units of dimensionful parameters
unit.ctmm <- function(CTMM,length=1,time=1)
{
  if(length(CTMM$tau)){ CTMM$tau <- CTMM$tau/time }
  CTMM$circle <- CTMM$circle * time

  # all means scale with length the same way... but not time
  CTMM$mu <- CTMM$mu/length
  drift <- get(CTMM$mean)
  CTMM <- drift@scale(CTMM,time)

  CTMM$error <- CTMM$error/length
  CTMM$sigma <- CTMM$sigma/length^2
  CTMM$sigma@par["area"] <- CTMM$sigma@par["area"]/length^2

  # variance -> diffusion adjustment
  if(!CTMM$range)
  {
    CTMM$sigma <- CTMM$sigma*time
    CTMM$sigma@par["area"] <- CTMM$sigma@par["area"]*time
  }

  if(!is.null(CTMM$COV))
  {
    CTMM$COV.mu <- CTMM$COV.mu/length^2

    CTMM$COV["area",] <- CTMM$COV["area",]/length^2
    CTMM$COV[,"area"] <- CTMM$COV[,"area"]/length^2

    if(!CTMM$range)
    {
      CTMM$COV["area",] <- CTMM$COV["area",]*time
      CTMM$COV[,"area"] <- CTMM$COV[,"area"]*time
    }

    tau <- CTMM$tau
    tau <- tau[tau<Inf]
    if(length(tau))
    {
      tau <- names(tau)
      tau <- paste("tau",tau)

      CTMM$COV[tau,] <- CTMM$COV[tau,]/time
      CTMM$COV[,tau] <- CTMM$COV[,tau]/time
    }

    if(CTMM$circle)
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

  alias <- list()
  scale <- c()

  add <- function(a,s)
  {
    n <- length(alias)
    alias[[n+1]] <<- a
    scale[n+1] <<- s
  }

  # TIME
  add(c("s","s.","sec","sec.","second","seconds"),1)
  add(c("min","min.","minute","minutes"),60)
  add(c("h","h.","hr","hr.","hour","hours"),60^2)
  add(c("day","days"),24*60^2)
  add(c("wk","wk.","week","weeks"),7*24*60^2)
  add(c("mon","mon.","month","months"),((29*24+12)*60+44)*60+2.8)
  add(c("yr","yr.","year","years"),365.24*7*24*60^2)

  # Distance conversions
  add(c("m","m.","meter","meters"),1)
  add(c("km","km.","kilometer","kilometers"),1000)
  add(c("in","in.","inch","inches"),0.3048/12)
  add(c("ft","ft.","foot","feet"),0.3048)
  add(c("yd","yd.","yard","yards"),0.3048*3)
  add(c("mi","mi.","mile","miles"),0.3048*5280)

  # Area conversions
  add(c("m\u00B2","m.\u00B2","meter\u00B2","meters\u00B2","m^2","m.^2","meter^2","meters^2","square meter","square meters","meter squared","meters squared"),1)
  add(c("ha","hectare","hectares"),100^2)
  add(c("km\u00B2","km.\u00B2","kilometer\u00B2","kilometers\u00B2","km^2","km.^2","kilometer^2","kilometers^2","square kilometer","square kilometers","kilometer squared","kilometers squared"),1000^2)
  add(c("in\u00B2","in.\u00B2","inch\u00B2","inches\u00B2","in^2","in.^2","inch^2","inches^2"),(0.3048/12)^2)
  add(c("ft\u00B2","ft.\u00B2","foot\u00B2","feet\u00B2","ft^2","ft.^2","foot^2","feet^2","square foot","square feet","foot squared","feet squared"),0.3048^2)
  add(c("yd\u00B2","yd.\u00B2","yard\u00B2","yards\u00B2","yd^2","yd.^2","yard^2","yards^2","square yard","square yards","yard squared","yards squared"),(0.3048*3)^2)
  add(c("mi\u00B2","mi.\u00B2","mile\u00B2","miles\u00B2","mi^2","mi.^2","mile^2","miles^2","square mile","square miles","mile squared","miles squared"),(0.3048*5280)^2)

  for(i in 1:length(alias))
  {
    if(name %in% alias[[i]]) { return(num*scale[i]^pow) }
  }
  stop(paste("Unit",name,"unknown."))
}
