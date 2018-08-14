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
uere <- function(data,override=FALSE,precision=1/2,trace=FALSE)
{
  data <- listify(data)

  TYPES <- get.dop.types(data) # types of DOP data present

  # loop over axes/types
  UERE <- lapply(TYPES,function(type){uere.type(data,override=override,precision=precision,trace=trace,type=type)})
  names(UERE) <- TYPES

  # integrate attempt
  # Do x,y,z axes all share the same UERE?
  # integrate <- FALSE # never seems to hold
  # if(integrate && all(c('horizontal','vertical') %in% TYPES) && all(c("HDOP","VDOP") %in% names(data[[1]])))
  # {
  #   P <- (UERE['horizontal']/UERE['vertical'])^2
  #   # F-test
  #   P <- min( stats::pf(P,N['horizontal'],N['vertical'],lower.tail=TRUE) , stats::pf(P,N['horizontal'],N['vertical'],lower.tail=FALSE) )
  #   message("p[UERE] = ",P)
  #
  #   UNITY <- sqrt((2*UERE['horizontal']^2 + UERE['vertical']^2)/3)
  #   UERE['horizontal'] <- UERE['vertical'] <- UNITY
  # }

  return(UERE)
}


# uere values for one set of axes
uere.type <- function(data,override=FALSE,trace=FALSE,type='horizontal',precision=1/2,...)
{
  TOL <- .Machine$double.eps^precision

  axes <- DOP.LIST[[type]]$axes
  z <- lapply(1:length(data),function(i){get.telemetry(data[[i]],axes)})

  # location classes
  if("class" %in% names(data[[1]]))
  {
    CLASS <- sapply(data,function(D){levels(D$class)})
    CLASS <- unique(c(CLASS))
  }
  else
  { CLASS <- "all" }

  # extract UEREs
  UERE <- lapply(data,function(D){attr(D,"UERE")[[type]]})
  # test for equivalence
  if(length(UERE)>1)
  {
    SAME <- sapply(UERE[-1],function(U){identical(UERE[[1]],U)})
    if(all(SAME)) { UERE <- UERE[[1]] }
    else { UERE <- NULL }
  }
  else
  { UERE <- UERE[[1]] }

  # null UERE structure
  if(is.null(UERE))
  {
    # null UERE structure given location classes
    UERE <- numeric(length(CLASS))
    names(UERE) <- CLASS
  }

  # what class UEREs will we be fitting
  if(override) { EST <- rep(TRUE,length(UERE)) }
  else { EST <- !UERE }
  names(EST) <- CLASS

  # DOP values
  DOP <- DOP.LIST[[type]]$DOP
  DOP <- lapply(data,function(D){ if(DOP %in% names(D)) { D[[DOP]] } else { rep(1,length(D$t)) } })

  # weights
  w <- lapply(DOP,function(D){length(axes)/D^2})

  # class indicators
  C <- lapply(data,function(D){get.class.mat(D,CLASS)})
  # ML DOF
  DOF.ML <- vapply(C,colSums,UERE) # (class,animals)
  dim(DOF.ML) <- c(length(UERE),length(data))
  DOF.ML <- rowSums(DOF.ML) # (class)

  # initial guess of UEREs
  UERE[EST] <- 10

  # iterative fit
  ERROR <- Inf
  while(ERROR>TOL)
  {
    # precicions
    P <- vapply(1:length(data),function(i){c(w[[i]] %*% C[[i]]) * (1/UERE)},UERE) # (class,animals)
    dim(P) <- c(length(UERE),length(data))
    Pc <- rowSums(P) # (class)
    Pk <- colSums(P) # (animals)

    # REML degrees of freedom
    DOF <- DOF.ML - colSums(t(P)/Pk) # (class)

    # means
    mu <- vapply(1:length(data),function(i){t(w[[i]] * z[[i]]) %*% C[[i]] %*% (1/UERE)},z[[1]][1,]) # (axes,animals)
    dim(mu) <- c(length(axes),length(data))
    mu <- t(mu)/Pk # (animals,axes)

    # updated UERE esitmates
    UERE2 <- vapply(1:length(data),function(i){colSums((t(z[[i]])-mu[i,])^2 %*% (w[[i]] * C[[i]]))},UERE) # (class,animals)
    dim(UERE2) <- c(length(UERE),length(data))
    UERE2 <- rowSums(UERE2) / DOF
    UERE2 <- sqrt(UERE2/length(axes))

    ERROR <- abs(UERE2[EST]-UERE[EST])/UERE2[EST]
    ERROR <- max(ERROR)

    UERE[EST] <- UERE2[EST]

    # stop after 1 run if 1 location class
    if(length(CLASS)<=1) { break }
  }

  attr(UERE,"DOF") <- DOF

  if(trace)
  {
    if(length(CLASS)>1)
    {
      for(i in 1:length(CLASS))
      {
        CI <- sqrt(chisq.ci(UERE[i]^2,DOF=DOF[i]))
        message("UERE[",type,",",CLASS[i],"] = ", sprintf("%.4f - %.4f",CI[1],CI[3]) )
      }
    }
    else
    {
      CI <- sqrt(chisq.ci(UERE^2,DOF=DOF))
      message("UERE[",type,"] = ", sprintf("%.4f - %.4f",CI[1],CI[3]) )
    }
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


# try to assign one UERE to one data
try.assign.uere <- function(data,UERE,TYPE)
{
  if(is.null(TYPE)) # universal assignment from unnamed UERE
  {
    # try every UERE type assignment if those axes are present in the data
    for(TYPE in names(DOP.LIST)) { if(all(DOP.LIST[[TYPE]]$axes %in% names(data))) {data <- try.assign.uere(data,UERE,TYPE) } }
    # attr(data,"UERE")[[TYPE]] <- UERE # overwrite individual UEREs
  }
  else # non-destructive assignment from named UERE
  {
    NAMES <- names(attr(data,"UERE"))
    # global variable DOP.LIST
    LIST <- DOP.LIST[[TYPE]]
    axes <- LIST$axes
    DOP <- LIST$DOP
    VAR <- LIST$VAR

    if(is.null(UERE) && (DOP %in% names(data))) # null assignment
    {
      data[[VAR]] <- NULL

      # delete UERE value
      attr(data,"UERE")[[TYPE]] <- 0
    }
    else # numeric assignment
    {
      # total UERE
      U <- c(get.class.mat(data,names(UERE)) %*% UERE)

      if(DOP %in% names(data))
      { data[[VAR]] <- (U*data[[DOP]])^2/length(axes) }
      else if(!(VAR %in% names(data)))
      {
        data[[VAR]] <- (U)^2/length(axes)
        message("DOP values missing. Assuming DOP=1.")
      }

      # overwrite UERE value
      attr(data,"UERE")[[TYPE]] <- UERE
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
  if(class(UERE)=="numeric") { UERE <- list(horizontal=UERE) }

  # promote to list and revert back if DROP
  DROP <- FALSE
  if(class(data)=="telemetry" || class(data)=="data.frame")
  {
    data <- list(data)
    DROP <- TRUE
  }

  for(i in 1:length(data)) { for(j in 1:max(1,length(UERE))) { data[[i]] <- try.assign.uere(data[[i]],UERE[[j]],TYPE=names(UERE)[j]) } }

  if(DROP) { data <- data[[1]] }

  # demote to data.frame
  return(data)
}


#####
# class indicator matrix
get.class.mat <- function(data,LEVELS=levels(data$class))
{
  if("class" %in% names(data))
  {
    C <- sapply(LEVELS,function(lc){data$class==lc}) # (time,class)
    colnames(C) <- LEVELS
  }
  else
  { C <- cbind(rep(1,length(data$t))) }

  return(C)
}
