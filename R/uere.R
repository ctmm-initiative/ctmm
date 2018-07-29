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
  if(!is.null(TYPE) && TYPE=="error") { TYPE <- NULL } # generic fitted UERE

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


# what DOP types are in the data
get.dop.types <- function(data)
{
  # types of data present
  TYPES <- names(DOP.LIST[-1])
  IN <- sapply(TYPES,function(TYPE) { all(DOP.LIST[[TYPE]]$axes %in% names(data[[1]])) })
  TYPES <- TYPES[IN]

  return(TYPES)
}


# return the UERE from set data or calculate it from calibration data
# look for every DOP without a VAR and calculate UERE - override option for all
# option to integrate HDOP & VDOP if they share the same UERE
uere <- function(data,override=FALSE,integrate=FALSE,trace=FALSE)
{
  data <- listify(data)

  if(!override)
  {
    UERE <- sapply(data,function(d){attr(d,"info")$UERE}) # (dim,tag)
    if(class(UERE)!="list")
    {
      UERE <- array(UERE,c(length(UERE)/length(data),length(data))) # sapply will not return consistent shape array
      UERE <- apply(UERE,1,unique)
      if(class(UERE)!="list" && is.null(dim(UERE))) { return(attr(data[[1]],"info")$UERE) } # unique UERE array to return
    }
  }

  TYPES <- get.dop.types(data) # types of DOP data present
  # remove cases where calibration exists
  if(length(TYPES) && !override)
  {
    OUT <- sapply(TYPES,function(TYPE) {DOP.LIST[[TYPE]]$VAR %in% names(data[[1]]) || all(DOP.LIST[[TYPE]]$COV %in% names(data[[1]]))})
    TYPES <- TYPES[!OUT]
  }

  # calculate residuals with UERE=1
  data <- residuals.calibration(data,TYPES)
  # data <- do.call(rbind,data)
  # estimate UERE from residuals
  UERE <- NULL
  N <- NULL
  for(TYPE in TYPES)
  {
    axes <- DOP.LIST[[TYPE]]$axes
    z <- NULL
    for(i in 1:length(data)) { z <- rbind(z,get.telemetry(data[[i]],axes)) }

    # total degrees of freedom
    N <- c(N,(nrow(z)-length(data)) * length(axes))
    # x,y data would be standardized if UERE==1 back from error=1 in CTMM list
    UERE <- c(UERE,sqrt(sum(z^2)/last(N)))

    names(N)[length(N)] <- TYPE
    names(UERE)[length(UERE)] <- TYPE

    if(trace)
    {
      CI <- sqrt(chisq.ci(UERE[TYPE]^2,DOF=N[TYPE]))
      message("UERE[",TYPE,"] = ", sprintf("%.4f - %.4f",CI[1],CI[3]) )
    }
  }

  # Do x,y,z axes all share the same UERE?
  if(integrate && all(c('horizontal','vertical') %in% TYPES) && all(c("HDOP","VDOP") %in% names(data[[1]])))
  {
    P <- (UERE['horizontal']/UERE['vertical'])^2
    # F-test
    P <- min( stats::pf(P,N['horizontal'],N['vertical'],lower.tail=TRUE) , stats::pf(P,N['horizontal'],N['vertical'],lower.tail=FALSE) )
    message("p[UERE] = ",P)

    UNITY <- sqrt((2*UERE['horizontal']^2 + UERE['vertical']^2)/3)
    UERE['horizontal'] <- UERE['vertical'] <- UNITY
  }

  return(UERE)
}


# calculate residuals of calibration data
# add axes argument for other uses?
residuals.calibration <- function(data,TYPES=get.dop.types(data),...)
{
  # enforce list structure
  data <- listify(data)

  for(TYPE in TYPES)
  {
    axes <- DOP.LIST[[TYPE]]$axes

    # get.error object for pulling UERE-1 variances
    CTMM <- list(axes=axes,error=1)

    # calculate mean and residuals
    for(i in 1:length(data))
    {
      # DOP <- c(DOP,data[[i]][[DOP.LIST[[TYPE]]$DOP]])
      # UERE-1 error variances
      w <- get.error(data[[i]],CTMM)
      if(attr(w,"flag")==0) { stop("Failed to import DOP column from Movebank file.") }
      # now these are the weights
      w <- 1/w
      # locations
      z <- get.telemetry(data[[i]],axes)
      # stationary mean
      mu <- c(w %*% z)/sum(w)
      # detrend the mean for error/residuals
      z <- t(t(z) - mu)
      # x <- rbind(x,dz)
      # UERE-1 standardize residuals
      z <- sqrt(w) * z
      # store back in the object
      data[[i]][,axes] <- z
    }
  }

  return(data)
}


## prepare error array, also return a flag #
# 0 : no error
# 1 : constant error parameter fit
# 2 : proportional error parameter fit to DOP value
# 3 : full error no fit (circle)
# 4 : " " (ellipse)
# circle : force variance   scalar output
# DIM : force covariance matrix output with dim [DIM,DIM]
get.error <- function(DATA,CTMM,flag=FALSE,circle=FALSE,DIM=FALSE)
{
  n <- nrow(DATA)
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

    if(all(COV %in% COLS)) # calibrated error ellipses - ARGOS
    {
      FLAG <- 4
      if(flag) { return(FLAG) }

      error <- get.telemetry(DATA,COV[c(1,2,2,3)]) # pull matrix elements
      dim(error) <- c(nrow(error),2,2) # array of matrices

      # reduce to VAR
      if(circle) { error <- (error[,1,1]+error[,2,2])/2 }
    }
    else if(VAR %in% COLS) # calibrated error circles - VAR=(UERE*HDOP)^2/2
    {
      FLAG <- 3
      if(flag) { return(FLAG) }

      error <- DATA[[VAR]]
    }
    else if(DOP %in% COLS) # fitted errors - HDOP
    {
      FLAG <- 2
      if(flag) { return(FLAG) }

      error <- (CTMM$error*DATA[[DOP]])^2/length(axes)
    }
    else # fitted errors - no HDOP
    {
      FLAG <- 1
      if(flag) { return(FLAG) }

      error <- rep(CTMM$error^2/length(axes),n)
    }
  } # END errors
  else # no error
  {
    FLAG <- 0
    if(flag) { return(FLAG) }

    error <- rep(0,n)
  }

  # upgrade variance scalar to covariance matrix
  if(FLAG<4 && DIM) { error <- outer(error,diag(DIM)) } # [n,d,d] }

  attr(error,"flag") <- FLAG
  return(error)
}
