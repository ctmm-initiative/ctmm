# global variables for dop/uere/error functions (special axes)
DOP.LIST <- list(unknown=list(axes=NA,geo=NA,DOP=NA,VAR=NA,COV=NA) ,
                 horizontal=list(axes=c("x","y"),geo=c("longitude","latitude"),DOP="HDOP",VAR="VAR.xy",COV=c("COV.x.x","COV.x.y","COV.y.y")),
                 vertical=list(axes="z",geo="z",DOP="VDOP",VAR="VAR.z",COV=NA),
                 speed=list(axes=c("vx","vy"),geo=c("speed","heading"),DOP="SDOP",VAR="VAR.v",COV=c("COV.vx.vx","COV.vx.vy","COV.vy.vy")) )


# is the data calibrated
is.calibrated <- function(data,type="horizontal")
{
  UERE <- attr(data,"UERE")

  # classes in data
  CLASS <- levels(data$class)
  if(!length(CLASS)) { CLASS <- "all" }

  # UERE of classes in data only
  UERE <- UERE[CLASS,type]
  UERE <- !is.na(UERE)
  UERE <- mean(UERE) # fraction calibrated

  return(UERE)
}

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

  data <- listify(data)
  # all kinds of data
  NAMES <- sapply(data,names)
  NAMES <- unique(c(NAMES))

  IN <- sapply(TYPES,function(type){all(DOP.LIST[[type]]$axes %in% NAMES) || all(DOP.LIST[[type]]$geo %in% NAMES)})
  TYPES <- TYPES[IN]

  return(TYPES)
}


# return the UERE from set data
uere <- function(data)
{
  if(class(data)=="list" && length(data)==1)
  { data <- data[[1]] }
  if(class(data)=="telemetry")
  {
    UERE <- attr(data,"UERE")
    return(UERE)
  }

  UERE <- lapply(data,function(d){attr(d,"UERE")})
  SAME <- sapply(UERE[-1],function(U){identical(U,UERE[[1]])})
  if(all(SAME)) { return( UERE[[1]] ) }
  else { return(UERE) }
}


# calculate UERE value from calibration data
# look for every DOP without a VAR and calculate UERE - override option for all
uere.fit <- function(data,precision=1/2)
{
  data <- listify(data)

  TYPES <- get.dop.types(data) # types of DOP data present

  # loop over axes/types
  UERE <- lapply(TYPES,function(type){uere.type(data,precision=precision,trace=trace,type=type)})
  names(UERE) <- TYPES

  # passed by slot
  DOF <- lapply(UERE,function(U){attr(U,"DOF")})
  names(DOF) <- TYPES
  # passed by slot
  AICc <- vapply(UERE,function(U){attr(U,"AICc")},numeric(1))
  names(AICc) <- TYPES
  # passed by slot
  Z2 <- vapply(UERE,function(U){attr(U,"Zsq")},numeric(1))
  names(Z2) <- TYPES

  CLASS <- names(UERE[[1]])
  if(!length(CLASS)) { CLASS <- "all" }

  # RRRRRRRRRR WHY DOES R DROP DIMENSIONS & DIMENSION NAMES ......
  UERE <- matrix( simplify2array(UERE) , length(CLASS) , length(TYPES) )
  dimnames(UERE) <- list(CLASS,TYPES)
  DOF <- matrix( simplify2array(DOF) , length(CLASS) , length(TYPES) )
  dimnames(DOF) <- list(CLASS,TYPES)

  UERE <- new.UERE(UERE,DOF=DOF,AICc=AICc,Zsq=Z2)

  return(UERE)
}


# uere values for one set of axes
uere.type <- function(data,trace=FALSE,type='horizontal',precision=1/2,...)
{
  TOL <- .Machine$double.eps^precision

  axes <- DOP.LIST[[type]]$axes # axes of data type
  DOP <- DOP.LIST[[type]]$DOP # DOP type

  for(i in 1:length(data))
  {
    # toss out Inf DOP values (used for missing data types)
    if(DOP %in% names(data[[i]]))
    {
      IN <- (data[[i]][[DOP]] < Inf)
      data[[i]] <- data[[i]][IN,]
    }
    else # make sure data has some DOP value
    { data[[i]][[DOP]] <- 1 }

    # null class cannot be NULL
    if("class" %nin% names(data[[i]]))
    {
      data[[i]]$class <- "all"
      data[[i]]$class <- as.factor(data[[i]]$class)
    }
  }

  # all location classes
  CLASS <- lapply(data,function(D){levels(D$class)})
  CLASS <- unique(unlist(CLASS))

  # null UERE structure given location classes
  UERE <- rep(NA_real_,length(CLASS))
  names(UERE) <- CLASS

  z <- lapply(1:length(data),function(i){get.telemetry(data[[i]],axes)})

  # make sure axes present in dataset (e.g., only 2D fixes)
  IN <- as.logical(sapply(z,length))
  data <- data[IN]
  z <- z[IN]

  # don't have enough data to estimate any UERE
  if(!length(data) || length(data)>=sum(sapply(data,nrow)))
  {
    attr(UERE,"DOF") <- numeric(length(UERE)) # sampling distribution
    attr(UERE,"AICc") <- NA_real_
    attr(UERE,"Zsq") <- NA_real_

    return(UERE)
  }

  EST <- is.na(UERE)
  names(EST) <- CLASS

  # DOP values # DOP values ensured above
  DOP <- lapply(data,function(D){ D[[DOP]] })

  # weights
  w <- lapply(DOP,function(D){length(axes)/D^2})

  # class indicators
  Ci <- lapply(data,function(D){get.class(D,CLASS)}) # (animal;time)
  Cim <- lapply(data,function(D){get.class.mat(D,CLASS)}) # (animal;time,class)

  # ML DOF
  DOF.ML <- vapply(Cim,colSums,UERE) # (class,animals)
  dim(DOF.ML) <- c(length(UERE),length(data))
  DOF.ML <- rowSums(DOF.ML) # (class)
  names(DOF.ML) <- CLASS

  # initial guess of UEREs
  UERE[EST] <- 10

  # check for missing levels/classes
  BAD <- (DOF.ML<=0)
  if(any(BAD))
  {
    # missing UEREs should not impact calculation
    UERE[BAD] <- Inf
    EST[BAD] <- FALSE
  }

  # iterative fit
  DEBUG <- FALSE
  BREAK <- FALSE
  repeat
  {
    if(type=="speed") # known mean
    {
      # ML/REML degrees of freedom lost
      Kc <- numeric(length(CLASS))

      # known means
      mu <- array(0,c(length(data),length(axes)))
    }
    else # unknown mean
    {
      # precicion weights
      Pck <- vapply(1:length(data),function(i){c(w[[i]] %*% Cim[[i]]) * (1/UERE^2)},UERE) # (class,animals)
      dim(Pck) <- c(length(UERE),length(data))
      Pc <- rowSums(Pck) # (class)
      Pk <- colSums(Pck) # (animals)

      # REML degrees of freedom lost
      Kc <- colSums(t(Pck)/Pk) # (class)

      # means
      CdUERE2 <- lapply(1:length(data),function(i){ 1/UERE[Ci[[i]]]^2 }) # (animals;time,class)
      mu <- vapply(1:length(data),function(i){t(w[[i]] * z[[i]]) %*% CdUERE2[[i]]},z[[1]][1,]) # (axes,animals)
      dim(mu) <- c(length(axes),length(data))
      mu <- t(mu)/Pk # (animals,axes)

      if(DEBUG)
      {
        loglike <- -(1/2)*vapply(1:length(data),function(i){
          LL <- length(axes) * sum( log(UERE[Ci[[i]]]^2/w[[i]]) )
          LL <- LL + sum( (t(z[[i]])-mu[i,])^2 %*% (w[[i]] * CdUERE2[[i]]) )
          LL <- LL + length(axes)*log( w[[i]] %*% CdUERE2[[i]] ) # REML
          return(LL) },numeric(1))
        loglike <- sum(loglike)
        message(type," log(Like) = ",loglike)
      }
    }
    DOF <- DOF.ML - Kc # (class)

    if(BREAK) { break }

    # updated UERE esitmates
    UERE2 <- vapply(1:length(data),function(i){colSums((t(z[[i]])-mu[i,])^2) %*% (w[[i]] * Cim[[i]])},UERE) # (class,animals)
    dim(UERE2) <- c(length(UERE),length(data))
    UERE2 <- rowSums(UERE2) / (length(axes)*DOF)
    UERE2 <- sqrt(UERE2)
    UERE2 <- UERE2[EST]

    ERROR <- abs(UERE2-UERE[EST])/UERE2
    ERROR <- max(ERROR)

    UERE[EST] <- UERE2

    if(ERROR<TOL) { BREAK <- TRUE }
  }

  # missing UEREs should not impact calculation
  if(any(BAD)) { UERE[BAD] <- 1 }

  ### AICc values ###
  AICc <- vapply(1:length(data),function(i){ sum( log( (2*pi)*(UERE[Ci[[i]]]^2 / w[[i]]) ) ) },numeric(1)) # (animals)
  AICc <- length(axes) * sum(AICc)
  if(type=="speed")
  { dof <- DOF.ML }
  else
  { dof <- DOF.ML - colSums((t(Pck)/Pk)^2) } # (class)
  # missing UEREs should not impact calculation
  if(any(BAD))
  {
    Kc[BAD] <- 0
    dof[BAD] <- 0
  }
  AICc <- AICc + length(axes)^2*sum( (DOF.ML+Kc)*dof/(length(axes)*dof-2) )
  if(length(axes)*dof<=2) # does not have a mean value
  {
    AICc <- Inf
    warning("Sample size too small for AICc to exist.")
  }

  ### Z_reduced^2 ###
  Pi <- lapply(1:length(data),function(k){c( w[[k]] / UERE[Ci[[k]]]^2 )}) # (animal;time)
  u2 <- lapply(1:length(data),function(k){c( colMeans((t(z[[k]])-mu[k,])^2) * Pi[[k]] )}) # (animal;time)
  if(type=="speed") # known mean
  {
    alpha <- lapply(1:length(data),function(k){ array(1,length(Pi[[k]])) }) # (animal;time)
    beta <- lapply(1:length(data),function(k){ ( DOF.ML[Ci[[k]]] - 1 ) / DOF[Ci[[k]]] }) # (animal;time)
    dofi <- lapply(1:length(data),function(k){ DOF.ML[Ci[[k]]] - 1 }) # (animal;time)
  }
  else
  {
    alpha <- lapply(1:length(data),function(k){ R <- Pck[Ci[[k]],k] ; ( R - Pi[[k]] )/R }) # (animal;time)
    gamma <- lapply(1:length(data),function(k){ t(Pck[Ci[[k]],,drop=FALSE] - Pi[[k]])/outer(Pk,Pi[[k]],'-') }) # (animal;animal',time)
    beta <- lapply(1:length(data),function(k){ ( DOF.ML[Ci[[k]]] - 1 - colSums(gamma[[k]]) ) / DOF[Ci[[k]]] }) # (animal;time)
    dofi <- lapply(1:length(data),function(k){ DOF.ML[Ci[[k]]] - 1 - colSums(gamma[[k]]^2) }) # (animal;time)
  }
  t2 <- lapply(1:length(data),function(k){ beta[[k]] * u2[[k]] / ( alpha[[k]]^2 - alpha[[k]]*u2[[k]]/DOF[Ci[[k]]] ) }) # (animal,time)
  Z <- log(unlist(t2))/2
  dofi <- unlist(dofi)
  M <- -(dofi-1)/dofi/(2*length(axes))
  VAR <- (dofi+1)/dofi/(2*length(axes))
  ERR <- VAR + M^2
  Z2 <- (Z^2/ERR)

  # fix missing UEREs
  if(any(BAD)) { UERE[BAD] <- NA }

  attr(UERE,"DOF") <- dof # sampling distribution
  attr(UERE,"AICc") <- AICc
  attr(UERE,"Zsq") <- mtmean(Z2) # minimally trimmed mean in case of weird speed=0 estimates

  return(UERE)
}


# summarize uere object
summary.UERE <- function(object,level=0.95,...)
{
  if(class(object)=="list") { return(summary.UERE.list(object,...)) }

  TYPE <- colnames(object)
  N <- sapply(TYPE,function(type){length(DOP.LIST[[type]]$axes)}) # (type)
  DOF <- attr(object,"DOF") # (class,type)
  DOF <- t(t(DOF)*N)

  UERE <- sapply(1:length(object),function(i){chisq.ci(object[i]^2,DOF=DOF[i],level=level)}) #(3,class*type)
  UERE <- sqrt(UERE)
  dim(UERE) <- c(3,dim(object)) # (3,class,type)
  dimnames(UERE)[[1]] <- c("low","ML","high")
  dimnames(UERE)[2:3] <- dimnames(object)
  UERE <- aperm(UERE,c(2,1,3))
  return(UERE)
}


# summarize list of uere objects
summary.UERE.list <- function(object,...)
{
  AICc <- sapply(object,function(U) { attr(U,"AICc") } )
  AICc <- AICc - min(AICc)
  Zsq <- sapply(object,function(U) { attr(U,"Zsq") } )

  TAB <- cbind(AICc,Zsq)
  colnames(TAB) <- c("\u0394AICc","Z[red]\u00B2")
  # Encoding(colnames(TAB)) <- "UTF-8"
  rownames(TAB) <- names(object)

  IND <- order(AICc)
  TAB <- TAB[IND,]

  return(TAB)
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
# some NA value means Inf VAR; all NA values means delete VAR column
try.assign.uere <- function(data,UERE,TYPE="horizontal")
{
  # global variable DOP.LIST
  LIST <- DOP.LIST[[TYPE]]
  axes <- LIST$axes
  DOP <- LIST$DOP
  VAR <- LIST$VAR
  COV <- LIST$COV

  # class indices
  Ci <- get.class(data,names(UERE))
  # partition in to NA-UERE and UERE>0
  NAS <- is.na(UERE)
  POS <- !NAS

  if(all(NAS))
  {
    # delete the VAR/COV columns
    data[[VAR]] <- NULL
    if(any(COV %in% names(data))) { data[COV] <- NULL }
  }
  else
  {
    # UERE[time] > 0
    if(any(NAS)) { UERE[NAS] <- 0 }
    U.POS <- UERE[Ci]

    if(DOP %in% names(data))
    { data[[VAR]] <- (U.POS*data[[DOP]])^2/length(axes) }
    else if(!(VAR %in% names(data)))
    {
      data[[VAR]] <- (U.POS)^2/length(axes)
      message("DOP values missing. Assuming DOP=1.")
    }

    # apply NA-UEREs
    if(any(NAS))
    {
      # indication of NA time
      U.NAS <- NAS[Ci]
      data[[VAR]][U.NAS] <- Inf
    }

    # covariance matrix (dual tagged)
    if(all(COV %in% names(data)))
    {
      data[[COV[1]]] <- data[[VAR]]
      data[[COV[2]]] <- 0
      data[[COV[3]]] <- data[[VAR]]
    }
  }

  return(data)
}


# create null UERE object for @UERE slot of telemetry data
uere.null <- function(data)
{
  if("class" %in% names(data))
  { CLASS <- levels(data$class) }
  else
  { CLASS <- "all" }

  TYPES <- get.dop.types(data)

  UERE <- matrix(NA_real_,length(CLASS),length(TYPES))
  rownames(UERE) <- CLASS
  colnames(UERE) <- TYPES

  UERE <- new.UERE(UERE,DOF=UERE,AICc=UERE[1,])

  return(UERE)
}

# store the UERE
# axes determined by name(value) consistent with DOP.LIST global variable above
# NULL is for universal UERE, else "horizontal", "vertical", "speed"
"uere<-" <- function(data,value)
{
  UERE <- value

  # default ambiguous assignment
  if(class(UERE)=="numeric" || class(UERE)=="integer")
  {
    UERE <- cbind(UERE)
    colnames(UERE) <- "horizontal"
    if(is.null(names(value))) { rownames(UERE) <- "all" }
    UERE <- new.UERE(UERE,DOF=NA*UERE,AICc=NA*UERE[1,])
  }

  DOF <- attr(UERE,"DOF")

  # promote to list and revert back if DROP
  DROP <- FALSE
  if(class(data)=="telemetry" || class(data)=="data.frame")
  {
    data <- list(data)
    DROP <- TRUE
  }

  for(i in 1:length(data))
  {
    # make sure of UERE slot compatibility
    data[[i]] <- droplevels(data[[i]])

    # all class/type from data columns
    UERE.NULL <- uere.null(data[[i]])
    CLASS <- rownames(UERE.NULL)
    TYPES <- colnames(UERE.NULL)

    if(length(UERE))
    {
      # all class/type in @UERE slot
      CLASS <- c(CLASS,rownames(attr(data[[i]],"UERE")))
      TYPES <- c(TYPES,colnames(attr(data[[i]],"UERE")))

      # all class/type in UERE-value
      CLASS <- c(CLASS,rownames(UERE))
      TYPES <- c(TYPES,colnames(UERE))

      CLASS <- unique(CLASS)
      TYPES <- unique(TYPES)
    }

    # setup null UERE
    UERE.NULL <- matrix(NA_real_,length(CLASS),length(TYPES))
    rownames(UERE.NULL) <- CLASS
    colnames(UERE.NULL) <- TYPES
    DOF.NULL <- UERE.NULL
    AICc.NULL <- UERE.NULL[1,]
    names(AICc.NULL) <- TYPES # R drops dimnames...

    # copy over @UERE slot
    if(length(attr(data[[i]],"UERE")) && length(UERE))
    {
      CLASS <- rownames(attr(data[[i]],"UERE"))
      TYPES <- colnames(attr(data[[i]],"UERE"))
      UERE.NULL[CLASS,TYPES] <- attr(data[[i]],"UERE")
      DOF.NULL[CLASS,TYPES] <- attr(attr(data[[i]],"UERE"),"DOF")
      AICc.NULL[TYPES] <- attr(attr(data[[i]],"UERE"),"AICc")
    }

    # copy over UERE-value (overwriting conflicts)
    if(length(UERE))
    {
      CLASS <- rownames(UERE)
      TYPES <- colnames(UERE)
      UERE.NULL[CLASS,TYPES] <- methods::getDataPart(UERE)
      DOF.NULL[CLASS,TYPES] <- attr(UERE,"DOF")
      AICc.NULL[TYPES] <- attr(UERE,"AICc")
    }

    # calibrate data
    TYPES <- get.dop.types(data[[i]])
    for(TYPE in TYPES)
    {
      # R does not preserve dimension names ???
      UERE.ROW <- UERE.NULL[,TYPE]
      names(UERE.ROW) <- rownames(UERE.NULL)
      data[[i]] <- try.assign.uere(data[[i]],UERE.ROW,TYPE=TYPE)
    }

    # store UERE that was applied
    UERE.NULL <- new.UERE(UERE.NULL,DOF=DOF.NULL,AICc=AICc.NULL)
    attr(data[[i]],"UERE") <- UERE.NULL
  }

  if(DROP) { data <- data[[1]] }

  # demote to data.frame
  return(data)
}


#####
# class index at times
# make robust to LEVELS being different from levels
get.class <- function(data,LEVELS=levels(data$class))
{
  if("class" %in% names(data))
  {
    Ci <- as.integer(data$class)

    Levels <- levels(data$class)
    if(any(LEVELS %nin% Levels)) # some LEVELS not in data
    {
      MISS <- LEVELS[LEVELS %nin% Levels]
      Levels <- c(Levels,MISS) # append missing LEVELS

      LEVELS <- sapply(Levels,function(L){ which(LEVELS==L) })
      Ci <- LEVELS[Ci]
    }
  }
  else
  { Ci <- rep(1,length(data$t)) }

  return(Ci)
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
