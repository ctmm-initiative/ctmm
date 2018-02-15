# global variables for dop/uere/error functions (special axes)
DOP.LIST <- list(unknown=list(axes=NA,DOP=NA,VAR=NA,COV=NA) ,
                 horizontal=list(axes=c("x","y"),DOP="HDOP",VAR="VAR.xy",COV=c("COV.x.x","COV.x.y","COV.y.y")),
                 vertical=list(axes="z",DOP="VDOP",VAR="VAR.z",COV=NA),
                 speed=list(axes=c("vx","vy"),DOP="SDOP",VAR="VAR.v",COV=c("COV.vx.vx","COV.vx.vy","COV.vy.vy")) )


# match DOP type for DOP.LIST by axes argument
DOP.match <- function(axes)
{
  DOP.LIST <- DOP.LIST[-1] # skip unknown case
  NAMES <- names(DOP.LIST)
  for(i in 1:length(DOP.LIST)) { if(all(axes==DOP.LIST[[i]]$axes)) { return(NAMES[i]) } }
  # match was not found
  warning("axes=",paste(axes,collapse=",")," not of known DOP type.")
  return("unknown")
}


# try to assign one UERE to one data
try.assign.uere <- function(data,UERE,TYPE=names(UERE))
{
  if(is.null(TYPE)) # universal assignment from unnamed UERE
  {
    # try every UERE type assignment if those axes are present in the data
    for(TYPE in names(DOP.LIST)) { if(all(DOP.LIST[[TYPE]]$axes %in% names(data))) {data <- try.assign.uere(data,UERE,TYPE) } }
    attr(data,"info")$UERE <- UERE # overwrite individual UEREs
  }
  else # non-destructive assignment from named UERE
  {
    NAMES <- names(attr(data,"info")$UERE)
    # global variable DOP.LIST
    LIST <- DOP.LIST[[TYPE]]
    axes <- LIST$axes
    DOP <- LIST$DOP
    VAR <- LIST$VAR

    if(is.null(UERE) && (DOP %in% names(data))) # null assignment
    {
      data[[VAR]] <- NULL

      # code to delete a possible element seems complicated in R ???
      if(TYPE %in% NAMES) { attr(data,"info")$UERE <- attr(data,"info")$UERE[-which(NAMES==TYPE)] }
      # else { attr(data,"info")$UERE <- NULL }
    }
    else # numeric assignment
    {
      if(DOP %in% names(data))
      { data[[VAR]] <- (UERE*data[[DOP]])^2/length(axes) }
      else if(!(VAR %in% names(data)))
      {
        data[[VAR]] <- (UERE)^2/length(axes)
        message("DOP values missing. Assuming DOP=1 homoskedastic errors.")
      }

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
# look for every DOP without a VAR and calculate UERE - override option for all
# option to integrate HDOP & VDOP if they share the same UERE
uere <- function(data,override=FALSE,integrate=FALSE,diagnostic=FALSE)
{
  if(class(data)=="telemetry" || class(data)=="data.frame")
  {
    # return present UERE
    UERE <- attr(data,"info")$UERE
    if(!is.null(UERE)) { return(UERE) }

    # promote to list of calibration data and proceed
    data <- list(data)
  }

  # types of data present
  TYPES <- names(DOP.LIST[-1])
  IN <- sapply(TYPES,function(TYPE) { all(DOP.LIST[[TYPE]]$axes %in% names(data[[1]])) })
  TYPES <- TYPES[IN]
  # remove cases where calibration exists
  if(length(TYPES) && !override)
  {
    OUT <- sapply(TYPES,function(TYPE) {DOP.LIST[[TYPE]]$VAR %in% names(data[[1]])})
    TYPES <- TYPES[!OUT]
  }

  UERE <- NULL
  N <- NULL
  for(TYPE in TYPES)
  {
    axes <- DOP.LIST[[TYPE]]$axes

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
    }
    # total degrees of freedom
    N <- c(N,(nrow(z)-length(data)) * length(axes))
    # x,y data would be standardized if UERE==1 back from error=1 in CTMM list
    UERE <- c(UERE,sqrt(sum(z^2)/last(N)))

    names(N)[length(N)] <- TYPE
    names(UERE)[length(UERE)] <- TYPE

    if(diagnostic)
    {
      r <- c(z) / (last(UERE))
      CI <- chisq.ci(1,DOF=N[TYPE])

      KDE <- stats::density(r,bw="SJ")
      graphics::hist(r,breaks="scott",freq=FALSE,main=paste(TYPE,"residual distribution"),xlab="standardized error",ylim=c(0,max(KDE$y,1/sqrt(2*pi*CI[1]))))
      graphics::lines(KDE)
      DNORM <- Vectorize(function(x){ exp(-x^2/2/CI[2])/sqrt(2*pi*CI[2]) })
      graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red",lwd=2)
      DNORM <- Vectorize(function(x){ exp(-x^2/2/CI[1])/sqrt(2*pi*CI[1]) })
      graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red")
      DNORM <- Vectorize(function(x){ exp(-x^2/2/CI[3])/sqrt(2*pi*CI[3]) })
      graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red")
      # should have used a t-distribution?

      ## right-justifying a set of labels: thanks to Uwe Ligges
      LEG <- sprintf("%.2f",rev(sqrt(chisq.ci(UERE[TYPE]^2,DOF=N[TYPE]))))
      LEN <- max(graphics::strwidth(LEG))
      OUT <- graphics::legend("topright",legend=rep(" ",length(LEG)),text.width=LEN,title="UERE",box.col=NA,col='red',bg=grDevices::grey(0.9),lwd=c(1,2,1),xjust=1,yjust=1)
      graphics::text(OUT$rect$left+OUT$rect$w,OUT$text$y,LEG,pos=2)
    }
  }

  if(integrate && all(c('horizontal','vertical') %in% TYPES) && all(c("HDOP","VDOP") %in% names(data[[1]])))
  {
    P <- (UERE['horizontal']/UERE['vertical'])^2
    # F-test
    P <- min( pf(P,N['horizontal'],N['vertical'],lower.tail=TRUE) , pf(P,N['horizontal'],N['vertical'],lower.tail=FALSE) )
    message("p[UERE] = ",P)

    UNITY <- sqrt((2*UERE['horizontal']^2 + UERE['vertical']^2)/3)
    UERE['horizontal'] <- UERE['vertical'] <- UNITY
  }

  return(UERE)
}


## prepare error array, also return a flag #
# 0 : no error
# 1 : constant error parameter fit
# 2 : proportional error parameter fit to DOP value
# 3 : full error no fit (circle)
# 4 : " " (ellipse)
get.error <- function(DATA,CTMM,flag=FALSE,circle=FALSE)
{
  n <- length(DATA$t)
  axes <- CTMM$axes
  COLS <- names(DATA)

  if(CTMM$error)
  {
    TYPE <- DOP.match(axes)
    # DOP.LIST is global variable from uere.R
    TYPE <- DOP.LIST[[TYPE]]
    AXES <- TYPE$axes
    COV <- TYPE$COV
    VAR <- TYPE$VAR
    DOP <- TYPE$DOP

    if(all(COV %in% COLS) && !circle) # calibrated error ellipses - ARGOS
    {
      error <- get.telemetry(DATA,COV[c(1,2,2,3)]) # pull matrix elements
      dim(error) <- c(nrow(error),2,2) # array of matrices
      FLAG <- 4
    }
    else if(VAR %in% COLS) # calibrated error circles - VAR=(UERE*HDOP)^2/2
    {
      error <- DATA[[VAR]]
      FLAG <- 3
    }
    else if(DOP %in% COLS) # fitted errors - HDOP
    {
      error <- (CTMM$error*DATA[[DOP]])^2/length(axes)
      FLAG <- 2
    }
    else # fitted errors - no HDOP
    {
      error <- rep(CTMM$error^2/length(axes),n)
      FLAG <- 1
    }
  } # END error
  else # no error
  {
    FLAG <- 0
    error <- rep(0,n)
  }

  if(flag) { return(FLAG) }
  else
  {
    attr(error,"flag") <- FLAG
    return(error)
  }
}
