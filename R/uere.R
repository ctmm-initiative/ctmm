# global variables for dop/uere functions (special axes)
DOP.LIST <- list(horizontal=list(axes=c("x","y"),DOP="HDOP",ERE="HERE"),
                 vertical=list(axes="z",DOP="VDOP",ERE="VERE"),
                 speed=list(axes=c("v.x","v.y"),DOP="SDOP",ERE="SERE"))


# match DOP type for DOP.LIST by axes argument
DOP.match <- function(axes)
{
  NAMES <- names(DOP.LIST)
  for(i in 1:length(DOP.LIST)) { if(all(axes==DOP.LIST[[i]]$axes)) { return(NAMES[i]) } }
  # match was not found
  stop("axes=",paste(axes,collapse=",")," not of known DOP type.")
}


# try to assign one UERE to one data
try.assign.uere <- function(data,UERE,TYPE=names(UERE))
{
  if(is.null(TYPE)) # universal assignment from unnamed UERE
  {
    # try every possible UERE type assignment
    for(TYPE in names(DOP.LIST)) { data <- try.assign.uere(data,UERE,TYPE) }
    attr(data,"info")$UERE <- UERE # overwrite individual UEREs
  }
  else # non-destructive assignment from named UERE
  {
    NAMES <- names(attr(data,"info")$UERE)
    # global variable DOP.LIST
    LIST <- DOP.LIST[[TYPE]]
    ERE <- LIST$ERE
    DOP <- LIST$DOP

    if(is.null(UERE)) # null assignment
    {
      data[,ERE] <- NULL

      # code to delete a possible element seems complicated in R ???
      if(TYPE %in% NAMES) { attr(data,"info")$UERE <- attr(data,"info")$UERE[-which(NAMES==TYPE)] }
      # else { attr(data,"info")$UERE <- NULL }
    }
    else # numeric assignment
    {
      if(DOP %in% names(data)) { data[[ERE]] <- UERE * data[[DOP]] }

      # store specific UERE if it exists already with names
      if(is.null(attr(data,"info")$UERE) || is.null(NAMES)) { attr(data,"info")$UERE <- UERE }
      else { attr(data,"info")$UERE[TYPE] <- UERE }
    }
  } # END non-destructive assingment

  return(data)
}

# store the UERE
# axes determined by name(value) consistent with DOP.LIST global variable above
# NULL is for universal UERE, else "horizontal", "vertical", "speed"
"uere<-" <- function(data,value)
{
  UERE <- value

  # promote to list and revert back if DROP
  DROP <- FALSE
  if(class(data)=="telemetry" || class(data)=="data.frame")
  {
    data <- list(data)
    DROP <- TRUE
  }

  for(i in 1:length(data)) { for(j in 1:max(1,length(UERE))) { data[[i]] <- try.assign.uere(data[[i]],UERE[j]) } }

  if(DROP) { data <- data[[1]] }

  # demote to data.frame
  return(data)
}

# return the UERE from set data or calculate it from calibration data
uere <- function(data,axes=c("x","y"),diagnostic=FALSE)
{
  if(class(data)=="telemetry" || class(data)=="data.frame")
  {
    # return present UERE
    UERE <- attr(data,"info")$UERE
    if(!is.null(UERE)) { return(UERE) }

    # promote to list of calibration data and proceed
    data <- list(data)
  }

  # get.error object for pulling UERE-1 variances
  CTMM <- list(axes=axes,error=1)

  # array of residuals to rbind to
  z <- NULL

  # calculate mean and residuals
  for(i in 1:length(data))
  {
    ## regular code
    # UERE-1 error variances
    w <- get.error(data[[i]],CTMM)
    if(attr(w,"flag")==0) { stop("Failed to import DOP column from Movebank file.") }
    # relative weights
    w <- 1/w
    # locations
    dz <- get.telemetry(data[[i]],axes)
    # stationary mean
    mu <- c(w %*% dz)/sum(w)
    # detrend the mean for error/residuals
    dz <- dz - mu
    # UERE-1 standardize residuals
    dz <- sqrt(w) * dz

    z <- rbind(z,dz)

    ## need special code for simultaneous 3D error calibration with HDOP & VDOP having the same UERE !!!
    #
    #
  }

  # total degrees of freedom
  n <- (nrow(z)-length(data)) * length(axes)
  # x,y data would be standardized if UERE==1 back from error=1 in CTMM list
  UERE <- sqrt(sum(z^2)/n)

  # name UERE properly
  names(UERE) <- DOP.match(axes)

  if(diagnostic)
  {
    r <- c(z) / (UERE)
    CI <- chisq.ci(1,DOF=n)

    KDE <- stats::density(r,bw="SJ")
    graphics::hist(r,breaks="scott",freq=FALSE,main="residual distribution",xlab="standardized error",ylim=c(0,max(KDE$y,1/sqrt(2*pi*CI[1]))))
    graphics::lines(KDE)
    DNORM <- Vectorize(function(x){ exp(-x^2/2/CI[2])/sqrt(2*pi*CI[2]) })
    graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red",lwd=2)
    DNORM <- Vectorize(function(x){ exp(-x^2/2/CI[1])/sqrt(2*pi*CI[1]) })
    graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red")
    DNORM <- Vectorize(function(x){ exp(-x^2/2/CI[3])/sqrt(2*pi*CI[3]) })
    graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red")
    # should have used a t-distribution?

    return(sqrt(chisq.ci(UERE^2,DOF=n)))
  }
  else
  { return( UERE ) }
}
