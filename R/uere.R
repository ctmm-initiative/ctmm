# global variables for dop/uere/error functions (special axes)
DOP.LIST <- list(
  unknown=list(axes=NA,geo=NA,DOP=NA,VAR=NA,COV=NA,COV.geo=NA,units=NA) ,
  horizontal=list(axes=c("x","y"),geo=c("longitude","latitude"),DOP="HDOP",VAR="VAR.xy",COV=c("COV.x.x","COV.x.y","COV.y.y"),COV.geo=c("COV.major","COV.minor","COV.angle"),units="distance"),
  vertical=list(axes="z",geo="z",DOP="VDOP",VAR="VAR.z",COV=NA,COV.geo=NA,units="distance"),
  speed=list(axes=c("vx","vy"),geo=c("speed","heading"),DOP="SDOP",VAR="VAR.v",COV=c("COV.vx.vx","COV.vx.vy","COV.vy.vy"),COV.geo=NA,units="speed"),
  frequency=list(axes='f',geo=NA,DOP=NA,VAR=NA,COV=NA,COV.geo=NA,units='frequency'),
  mass=list(axes='m',geo=NA,DOP=NA,VAR=NA,COV=NA,COV.geo=NA,units='mass')
)


# are the data calibrated
is.calibrated <- function(data,type="horizontal")
{
  if(class(data)[1]=="list") { return( mean( sapply(data,is.calibrated) ) ) }

  DOF <- attr(data,"UERE")$DOF # RMS UERE matrix

  # classes in data
  CLASS <- levels(data$class)
  if(!length(CLASS))
  {
    CLASS <- "all"
    if(is.null(dimnames(DOF))) { dimnames(DOF) <- list(CLASS,type) }
  }

  # UERE of classes in data only
  DOF <- DOF[CLASS,type]
  DOF <- !is.na(DOF) & DOF>0
  DOF <- mean(DOF) # fraction calibrated

  return(DOF)
}

# match DOP type for DOP.LIST by axes argument
DOP.match <- function(axes)
{
  DOP.LIST <- DOP.LIST[-1] # skip unknown case
  NAMES <- names(DOP.LIST)
  for(i in 1:length(DOP.LIST)) { if(all(axes==DOP.LIST[[i]]$axes)) { return(NAMES[i]) } }
  # match was not found
  # warning("axes=",paste(axes,collapse=",")," not of known DOP type.")
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
  NAMES <- unique(unlist(NAMES))

  IN <- sapply(TYPES,function(type){all(DOP.LIST[[type]]$axes %in% NAMES) || all(DOP.LIST[[type]]$geo %in% NAMES)})
  TYPES <- TYPES[IN]

  return(TYPES)
}


# return the RMS UERE object from set data
uere <- function(data)
{
  if(class(data)[1]=="list" && length(data)==1)
  { data <- data[[1]] }

  if(class(data)[1] %in% c("telemetry","variogram"))
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
  LIST <- lapply(TYPES,function(type){uere.type(data,precision=precision,trace=trace,type=type)})

  UERE <- lapply(LIST,function(L){ L$UERE })
  names(UERE) <- TYPES

  DOF <- lapply(LIST,function(L){ L$DOF })
  names(DOF) <- TYPES

  AICc <- sapply(LIST,function(L){ L$AICc })
  names(AICc) <- TYPES

  Zsq <- sapply(LIST,function(L){ L$Zsq })
  names(Zsq) <- TYPES

  VAR.Zsq <- sapply(LIST,function(L){ L$VAR.Zsq })
  names(VAR.Zsq) <- TYPES

  N <- sapply(LIST,function(L){ L$N })
  names(N) <- TYPES

  CLASS <- names(UERE[[1]])
  if(!length(CLASS)) { CLASS <- "all" }

  # RRRRRRRRRR WHY DOES R DROP DIMENSIONS & DIMENSION NAMES ......
  UERE <- matrix( simplify2array(UERE) , length(CLASS) , length(TYPES) )
  dimnames(UERE) <- list(CLASS,TYPES)
  DOF <- matrix( simplify2array(DOF) , length(CLASS) , length(TYPES) )
  dimnames(DOF) <- list(CLASS,TYPES)

  UERE <- list(UERE=UERE,DOF=DOF,AICc=AICc,Zsq=Zsq,VAR.Zsq=VAR.Zsq,N=N)
  UERE <- new.UERE(UERE)

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
  # don't overwrite ARGOS when calibrating GPS
  ARGOS <- grepl('argos',tolower(CLASS))
  if(any(ARGOS) && type=='horizontal') { UERE[ARGOS] <- 1 }

  z <- lapply(1:length(data),function(i){get.telemetry(data[[i]],axes)})

  # make sure axes present in dataset (e.g., only 2D fixes)
  if(type=='speed') { IN <- sapply(z,length)>0 } else { IN <- sapply(z,length)>length(axes) }
  data <- data[IN]
  z <- z[IN]

  # don't have enough data to estimate any UERE
  if(!length(data))
  {
    UERE <- list(UERE=UERE)
    UERE$DOF <- numeric(length(UERE)) # sampling distribution
    UERE$AICc <- NA_real_
    UERE$Zsq <- NA_real_
    UERE$VAR.Zsq <- NA_real_
    UERE$N <- NA_real_

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
      # precision weights
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

    ERROR <- abs(UERE2-UERE[EST])/max(UERE[EST],UERE2)
    ERROR <- nant(ERROR,0)
    ERROR <- max(ERROR)

    UERE[EST] <- UERE2

    if(ERROR<TOL) { BREAK <- TRUE }
  }

  ### AICc values ###
  if(any(BAD)) { UERE[BAD] <- NA } # missing UEREs should not impact calculation
  AICc <- vapply(1:length(data),function(i){ sum( log( (2*pi)*(UERE[Ci[[i]]]^2 / w[[i]]) ), na.rm=TRUE) },numeric(1)) # (animals)
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
  AICc <- AICc + length(axes)^2*sum( ((DOF.ML+Kc)*dof/(length(axes)*dof-2))[EST] )
  if(any(length(axes)*dof[EST]<=2)) # does not have a mean value
  {
    AICc <- Inf
    warning("Sample size too small for AICc to exist.")
  }

  ### Z_reduced^2 ###
  if(any(BAD)) { UERE[BAD] <- Inf } # missing UEREs should not impact calculation
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
  t2 <- lapply(1:length(data),function(k){ clamp( beta[[k]] * u2[[k]] / ( alpha[[k]]^2 - alpha[[k]]*u2[[k]]/DOF[Ci[[k]]] ),0,Inf) }) # (animal,time)
  Z <- log(unlist(t2))/2
  # M <- -clamp(dofi-1,0,Inf)/dofi/(2*length(axes)) # asymptotic only
  # VAR <- (dofi+1)/dofi/(2*length(axes)) # asymptotic only
  d1 <- length(axes)
  d2 <- unlist(dofi)
  if(length(CLASS)==1) { d2 <- d2[1] } # minimize computation if possible
  d2 <- length(axes) * clamp(d2,0,Inf) # finish and fix for tiny samples in a class
  M2 <- function(d2)
  {
    if(d2==0) { return(Inf) }

    d1 <- d1/2
    d2 <- d2/2

    L1 <- -(log(d1)-digamma(d1))/2
    L2 <- -(log(d2)-digamma(d2))/2

    Q1 <- L1^2 + trigamma(d1)/4
    Q2 <- L2^2 + trigamma(d2)/4

    R <- Q1 + Q2 - 2*L1*L2

    return(R)
  }
  V2 <- sapply(d2,M2)
  Z2 <- (Z^2/V2)

  # fix missing UEREs
  if(any(BAD)) { UERE[BAD] <- NA }
  if(any(ARGOS) && type=='horizontal') { dof[ARGOS] <- Inf }

  UERE <- list(UERE=UERE)
  UERE$DOF <- dof # sampling distribution
  UERE$AICc <- AICc
  UERE$Zsq <- mtmean(Z2) # minimally trimmed mean in case of weird speed=0 estimates
  UERE$VAR.Zsq <- mtmean(Z2^2) - UERE$Zsq^2
  UERE$N <- sum(Z2<Inf)

  return(UERE)
}


# summarize uere object
summary.UERE <- function(object,level=0.95,...)
{
  if(class(object)[1]=="list") { return(summary.UERE.list(object,level=level,...)) }

  TYPE <- colnames(object$UERE)
  N <- sapply(TYPE,function(type){length(DOP.LIST[[type]]$axes)}) # (type)
  DOF <- object$DOF # (class,type)
  DOF <- t(t(DOF)*N)

  UERE <- sapply(1:length(object$UERE),function(i){chisq.ci(object$UERE[i]^2,DOF=DOF[i],level=level)}) #(3CI,class*type)
  UERE <- sqrt(UERE)
  dim(UERE) <- c(3,dim(object$UERE)) # (3CI,class,type)
  dimnames(UERE)[[1]] <- NAMES.CI
  dimnames(UERE)[2:3] <- dimnames(object$UERE)
  UERE <- aperm(UERE,c(2,1,3)) # (class,3CI,type)
  return(UERE)
}


# summarize list of uere objects
summary.UERE.list <- function(object,level=0.95,drop=TRUE,CI=FALSE,...)
{
  NAMES <- names(object)
  if(is.null(NAMES)) { NAMES <- 1:length(object) }

  # determine all types of data
  # better be consistent across joint models
  # might not be consistent within joint models' individual models
  TYPES <- NULL
  for(i in 1:length(object))
  {
    if(class(object[[i]])[1]=="UERE") # these should be all the same
    {
      TYPES <- c(TYPES,names(object[[i]]$AICc))
      TYPES <- unique(TYPES)
    }
    else if(class(object[[i]])[1]=="list") # these can differ within list
    {
      for(j in 1:length(object[[i]]))
      {
        TYPES <- c(TYPES,names(object[[i]][[j]]$AICc))
        TYPES <- unique(TYPES)
      }
    }
  }

  # aggregate lists of individual models to compare with joint models
  for(i in 1:length(object))
  {
    if(class(object[[i]])[1]=="list")
    {
      AICc <- rep(0,length(TYPES))
      names(AICc) <- TYPES
      N <- AICc
      Zsq <- AICc
      VAR.Zsq <- AICc

      for(j in 1:length(object[[i]]))
      {
        U <- object[[i]][[j]]
        IN <- TYPES %in% names(U$AICc)
        if(length(IN))
        {
          AICc[IN] <- AICc[IN] + U$AICc[IN]
          N[IN] <- N[IN] + U$N[IN]
          Zsq[IN] <- Zsq[IN] + U$N[IN] * U$Zsq[IN]
          VAR.Zsq[IN] <- VAR.Zsq[IN] + U$N[IN] * (U$VAR.Zsq[IN] + U$Zsq[IN]^2)
        }
      }

      Zsq <- Zsq/N
      VAR.Zsq <- VAR.Zsq/N - Zsq^2

      object[[i]]$AICc <- AICc
      object[[i]]$N <- N
      object[[i]]$Zsq <- Zsq
      object[[i]]$VAR.Zsq <- VAR.Zsq
    }
  }

  AIC <- sapply(object,function(U) { U$AICc } )
  N <- sapply(object,function(U) { U$N } )
  Zsq <- sapply(object,function(U) { U$Zsq } )
  VAR.Zsq <- sapply(object,function(U) { U$VAR.Zsq } )

  DIM <- dim(AIC)
  if(!length(DIM))
  {
    DIM <- c(1,length(AIC))
    dim(AIC) <- DIM
    dim(Zsq) <- DIM
    dim(VAR.Zsq) <- DIM
    dim(N) <- DIM
  }

  TAB <- list()
  for(i in 1:DIM[1])
  {
    AIC[i,] <- AIC[i,] - min(AIC[i,])

    if(CI)
    {
      # is N correct here for 1D and 2D ???
      TAB[[i]] <- sapply(1:DIM[2],function(j){ chisq.ci(Zsq[i,j],VAR=VAR.Zsq[i,j]/N[i,j],level=level) }) # (3CIS,models)
      TAB[[i]] <- cbind(AIC[i,],t(TAB[[i]]))
      colnames(TAB[[i]]) <- c("\u0394AICc","(       ","Z[red]\u00B2","       )")
    }
    else
    {
      TAB[[i]] <- cbind(AIC[i,],Zsq[i,])
      colnames(TAB[[i]]) <- c("\u0394AICc","Z[red]\u00B2")
    }
    # Encoding(colnames(TAB)) <- "UTF-8"


    rownames(TAB[[i]]) <- NAMES

    IND <- order(AIC[i,])
    TAB[[i]] <- TAB[[i]][IND,]
  }
  names(TAB) <- rownames(AIC)

  if(length(TAB)==1 && drop) { return(TAB[[1]]) }
  else { return(TAB) }
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

    # calculate mean and residuals
    for(i in 1:length(data))
    {
      n <- nrow(data[[i]])

      if(is.calibrated(data[[i]])<1) { uere(data[[i]]) <- uere(data[[i]]) } # force calibration for plotting
      UERE <- uere(data[[i]])
      ERROR <- UERE$UERE[,'horizontal']
      names(ERROR) <- rownames(UERE$UERE) # R drops dimnames
      ERROR <- ctmm(error=ERROR,axes=axes)
      ERROR <- get.error(data[[i]],ERROR,calibrate=TRUE)
      ELLIPSE <- attr(ERROR,"ellipse")

      # locations
      z <- get.telemetry(data[[i]],axes)

      if(!ELLIPSE)
      {
        # now these are the weights
        w <- 1/ERROR
        W <- sum(w)

        # stationary mean
        mu <- c(w %*% z)/W

        # detrend the mean for error/residuals
        z <- t(t(z) - mu)
        # x <- rbind(x,dz)

        # UERE-1 standardize residuals
        z <- sqrt(w) * z
      }
      else
      {
        w <- vapply(1:n,function(j){PDsolve(ERROR[j,,])},diag(2)) # [x,x,t]
        dim(w) <- c(2,2,n)
        w <- aperm(w,c(3,1,2)) # [t,x,x]
        W <- apply(w,2:3,sum) # [x,y]

        # stationary mean
        mu <- vapply(1:n,function(j){w[j,,] %*% z[j,]},1:2) # [x,t]
        mu <- apply(mu,1,sum)
        mu <- c(PDsolve(W) %*% mu)

        # detrend the mean for error/residuals
        z <- t(t(z) - mu)

        # standardize/studentize
        w <- vapply(1:n,function(j){sqrtm(ERROR[j,,])},diag(2)) # [x,x,t]
        dim(w) <- c(2,2,n)
        w <- aperm(w,c(3,1,2)) # [t,x,x]

        z <- vapply(1:n,function(j){w[j,,] %*% z[j,]},1:2) # [x,t]
        z <- t(z) # [t,x]
      }

      # store back in the object
      data[[i]][,axes] <- z

      uere(data) <- NULL
    } # end data loop
  } # end type loop

  return(data)
}


## prepare error array, also return a ellipse necessity #
# circle : force variance   scalar output
# DIM : force covariance matrix output with dim [DIM,DIM]
get.error <- function(data,CTMM,circle=FALSE,DIM=FALSE,calibrate=TRUE)
{
  n <- nrow(data)
  axes <- CTMM$axes
  COLS <- names(data)
  ELLIPSE <- FALSE
  error <- rep(0,n)

  if(any(CTMM$error>0))
  {
    TYPE <- DOP.match(axes)
    UERE.DOF <- attr(data,"UERE")$DOF[,TYPE]
    names(UERE.DOF) <- rownames(attr(data,"UERE")$DOF)

    # expand to what classes are in the UERE object
    ERROR <- rep(FALSE,length(UERE.DOF))
    names(ERROR) <- names(UERE.DOF)
    ERROR[names(CTMM$error)] <- CTMM$error

    UERE.FIT <- ERROR & !is.na(UERE.DOF) & UERE.DOF<Inf # will we be fitting any error parameters?
    UERE.FIX <- ERROR & (is.na(UERE.DOF) | UERE.DOF==Inf)
    UERE.PAR <- names(UERE.FIT)[UERE.FIT>0] # names of fitted UERE parameters

    # make sure that non-fitted parameters are logicals
    if(any(!UERE.FIT)) { ERROR[!UERE.FIT] <- as.logical(ERROR[!UERE.FIT]) }

    # reduce to what classes are actually in the data - also fixes bad sorting
    if("class" %in% names(data))
    {
      LEVELS <- levels(data$class)
      ERROR <- ERROR[LEVELS]
      UERE.DOF <- UERE.DOF[LEVELS]
      UERE.FIT <- UERE.FIT[LEVELS]
      UERE.FIX <- UERE.FIX[LEVELS]
      UERE.PAR <- UERE.PAR[UERE.PAR %in% LEVELS]
    }

    CLASS <- get.class.mat(data)
    FIT <- as.logical(c(CLASS %*% UERE.FIT)) # times where errors will be fitted (possibly with prior)
    # FIX <- as.logical(c(CLASS %*% UERE.FIX)) # times where errors are fixed

    # DOP.LIST is global variable from uere.R
    TYPE <- DOP.LIST[[TYPE]]
    AXES <- TYPE$axes
    COV <- TYPE$COV
    VAR <- TYPE$VAR
    DOP <- TYPE$DOP

    ELLIPSE <- all(COV %in% COLS)
    CIRCLE <- VAR %in% COLS

    # location classes with fixed (known) parameters
    if(ELLIPSE) # at least some (fixed) error ellipses - ARGOS / VHF - this should also include other fixed errors (promoted to ellipses)
    {
      ellipse <- get.telemetry(data,COV[c(1,2,2,3)]) # pull matrix elements
      dim(ellipse) <- c(nrow(ellipse),2,2) # array of matrices

      # reduce to VAR
      if(circle)
      {
        error <- (ellipse[,1,1]+ellipse[,2,2])/2
        ELLIPSE <- FALSE
      }
    }
    else if(CIRCLE) # some fixed error circles (no ellipses)
    {
      error <- data[[VAR]] # FITs will overwrite where appropriate
    }

    # location classes with fitted error parameters
    if(any(UERE.FIT))
    {
      if(calibrate) # apply RMS UERE to DOP values
      { CLASS <- c( CLASS %*% ERROR ) }
      else # just calculate relative variance
      { CLASS <- c( CLASS %*% as.logical(ERROR) ) }

      if(DOP %in% COLS) # apply DOP values
      { CLASS <- CLASS * data[[DOP]] }

      CLASS <- CLASS^2/length(axes) # VAR=(RMS[UERE]*HDOP)^2/2 in 2D
      error[FIT] <- CLASS[FIT]
      # FIX was already taken care of in if(ELLIPSE) & if(CIRCLE)
    }

    # promote circular errors to elliptical errors
    if(ELLIPSE)
    {
      # copy over error circles
      if(any(UERE.FIT))
      {
        ellipse[FIT,1,1] <- ellipse[FIT,2,2] <- error[FIT]
        ellipse[FIT,1,2] <- ellipse[FIT,2,1] <- 0
      }

      error <- ellipse
    }
  } # END errors

  # upgrade variance scalar to covariance matrix
  if(!ELLIPSE && !circle && DIM) { error <- outer(error,diag(DIM)) } # [n,d,d] }

  attr(error,"ellipse") <- ELLIPSE
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

  # don't overwrite ARGOS when calibrating GPS
  ARGOS <- grepl('argos',tolower(levels(data$class)))
  if(TYPE=='horizontal' && any(ARGOS))
  {
    ARGOS <- levels(data$class)[ARGOS]
    SAFE <- (data$class != ARGOS)
  }
  else
  { SAFE <- rep(TRUE,length(data$t)) }

  if(all(NAS))
  {
    # delete the VAR/COV columns
    data[[VAR]] <- NULL
    if(any(COV %in% names(data))) { data[COV] <- NULL }
  }
  else if(any(SAFE))
  {
    # UERE[time] > 0
    if(any(NAS)) { UERE[NAS] <- 0 }
    U.POS <- UERE[Ci]

    if(DOP %in% names(data))
    { data[[VAR]][SAFE] <- (U.POS*data[[DOP]])[SAFE]^2/length(axes) }
    else if(!(VAR %in% names(data)))
    {
      data[[VAR]][SAFE] <- (U.POS)[SAFE]^2/length(axes)
      message("DOP values missing. Assuming DOP=1.")
    }

    # apply NA-UEREs
    if(any(NAS))
    {
      # indication of NA time
      U.NAS <- NAS[Ci]
      data[[VAR]][U.NAS & SAFE] <- Inf
    }

    # covariance matrix (dual tagged)
    if(all(COV %in% names(data)))
    {
      data[[COV[1]]][SAFE] <- data[[VAR]][SAFE]
      data[[COV[2]]][SAFE] <- 0
      data[[COV[3]]][SAFE] <- data[[VAR]][SAFE]
    }

    if(TYPE=="horizontal" && all(LIST$COV.geo %in% names(data)))
    {
      data$COV.major[SAFE] <- data[[VAR]][SAFE]
      data$COV.minor[SAFE] <- data[[VAR]][SAFE]
      data$COV.angle[SAFE] <- 0
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

  M <- matrix(1,length(CLASS),length(TYPES))
  rownames(M) <- CLASS
  colnames(M) <- TYPES

  V <- M[1,]
  names(V) <- colnames(V)

  UERE <- list(UERE=1*M,DOF=0*M,AICc=Inf*V,Zsq=Inf*V,VAR.Zsq=Inf*V,N=0*V)
  UERE <- new.UERE(UERE)

  return(UERE)
}


# store the UERE
# axes determined by name(value) consistent with DOP.LIST global variable above
# NULL is for universal UERE, else "horizontal", "vertical", "speed"
"uere<-" <- function(data,value)
{
  UERE <- value

  # default ambiguous assignment - overrides everything
  if(class(UERE)[1]=="numeric" || class(UERE)[1]=="integer")
  {
    # in case of different location classes
    if(class(data)[1]=="list")
    {
      data <- lapply(data,function(d){"uere<-"(d,value)})
      return(data)
    }

    # at some point we might want to check the dimensions of value for different types

    uere(data) <- NULL
    UERE <- uere(data)
    UERE$UERE[,] <- value
    UERE$DOF[] <- Inf
    UERE$AICc[] <- Inf
    UERE$Zsq[] <- Inf
    UERE$VAR.Zsq[] <- Inf
    UERE$N[] <- Inf
    uere(data) <- UERE
    return(data)
  }
  else if(class(UERE)[1]=="ctmm")
  {
    # in case of different location classes
    if(class(data)[1]=="list")
    {
      data <- lapply(data,function(d){"uere<-"(d,value)})
      return(data)
    }

    # only calibrate axis of value
    axes <- value$axes
    AXES <- length(axes)
    TYPE <- DOP.match(axes)
    ERROR <- value$error
    CLASS <- names(ERROR)
    ECLASS <- paste("error",CLASS)

    if(all(ECLASS %in% rownames(value$COV)))
    {
      # sqrt(chi^2) relation
      N <- 2/AXES * ERROR^2 / (2^2*diag(value$COV[ECLASS,ECLASS,drop=FALSE]))
      names(N) <- CLASS
    }
    else
    {
      N <- UERE
      N[] <- Inf
    }

    UERE <- uere(data)
    CLASS <- CLASS[CLASS %in% rownames(UERE$UERE)]
    ERROR <- ERROR[CLASS]
    N <- N[CLASS]

    UERE$UERE[CLASS,TYPE] <- ERROR
    UERE$DOF[CLASS,TYPE] <- N
    UERE$AICc[TYPE] <- Inf
    UERE$Zsq[TYPE] <- Inf
    UERE$VAR.Zsq[TYPE] <- Inf
    UERE$N[TYPE] <- sum(N)
    uere(data) <- UERE
    return(data)
  }

  DOF <- UERE$DOF

  # promote to list and revert back if DROP
  DROP <- FALSE
  if(class(data)[1]=="telemetry" || class(data)[1]=="data.frame")
  {
    data <- list(data)
    DROP <- TRUE
  }

  for(i in 1:length(data))
  {
    # make sure of UERE slot compatibility
    # data[[i]] <- droplevels(data[[i]])
    # why was I doing this?

    # all class/type from data columns
    UERE.NULL <- uere.null(data[[i]])
    CLASS <- rownames(UERE.NULL$UERE)
    TYPES <- colnames(UERE.NULL$UERE)

    if(length(UERE$UERE))
    {
      # all class/type in @UERE slot
      CLASS <- c(CLASS,rownames(attr(data[[i]],"UERE")$UERE))
      TYPES <- c(TYPES,colnames(attr(data[[i]],"UERE")$UERE))

      # all class/type in UERE-value
      CLASS <- c(CLASS,rownames(UERE$UERE))
      TYPES <- c(TYPES,colnames(UERE$UERE))

      CLASS <- unique(CLASS)
      TYPES <- unique(TYPES)
    }

    # setup null UERE
    UERE.NULL <- matrix(1,length(CLASS),length(TYPES))
    rownames(UERE.NULL) <- CLASS
    colnames(UERE.NULL) <- TYPES
    DOF.NULL <- 0*UERE.NULL
    V <- UERE.NULL[1,]
    names(V) <- TYPES # R drops dimnames...
    AICc.NULL <- Inf*V
    Zsq.NULL <- AICc.NULL
    VAR.Zsq.NULL <- AICc.NULL
    N.NULL <- 0*V

    # copy over @UERE slot
    if(length(attr(data[[i]],"UERE")) && length(UERE$UERE))
    {
      CLASS <- rownames(attr(data[[i]],"UERE")$UERE)
      TYPES <- colnames(attr(data[[i]],"UERE")$UERE)
      UERE.NULL[CLASS,TYPES] <- attr(data[[i]],"UERE")$UERE
      DOF.NULL[CLASS,TYPES] <- attr(data[[i]],"UERE")$DOF
      AICc.NULL[TYPES] <- attr(data[[i]],"UERE")$AICc
      Zsq.NULL[TYPES] <- attr(data[[i]],"UERE")$Zsq
      VAR.Zsq.NULL[TYPES] <- attr(data[[i]],"UERE")$VAR.Zsq
      N.NULL[TYPES] <- attr(data[[i]],"UERE")$N
    }

    # copy over UERE-value (overwriting conflicts)
    if(length(UERE))
    {
      CLASS <- rownames(UERE$UERE)
      TYPES <- colnames(UERE$UERE)
      UERE.NULL[CLASS,TYPES] <- UERE$UERE
      DOF.NULL[CLASS,TYPES] <- UERE$DOF
      AICc.NULL[TYPES] <- UERE$AICc
      Zsq.NULL[TYPES] <- UERE$Zsq
      VAR.Zsq.NULL[TYPES] <- UERE$VAR.Zsq
      N.NULL[TYPES] <- UERE$N
    }

    # calibrate data
    TYPES <- get.dop.types(data[[i]])
    for(TYPE in TYPES)
    {
      # R does not preserve dimension names ???
      UERE.ROW <- UERE.NULL[,TYPE]
      names(UERE.ROW) <- rownames(UERE.NULL)
      # don't calibrate data if NULL assignment
      if(!is.null(value)) { data[[i]] <- try.assign.uere(data[[i]],UERE.ROW,TYPE=TYPE) }
      else { data[[i]] <- try.assign.uere(data[[i]],NA*UERE.ROW,TYPE=TYPE) }
    }

    # store UERE that was applied
    UERE.NULL <- list(UERE=UERE.NULL,DOF=DOF.NULL,AICc=AICc.NULL,Zsq=Zsq.NULL,VAR.Zsq=VAR.Zsq.NULL,N=N.NULL)
    UERE.NULL <- new.UERE(UERE.NULL)
    attr(data[[i]],"UERE") <- UERE.NULL
  }

  if(DROP) { data <- data[[1]] }

  # demote to data.frame
  return(data)
}


#####
# class index at times
get.class <- function(data,LEVELS=levels(data$class))
{
  if("class" %in% names(data))
  {
    Ci <- rep(NA_integer_,nrow(data))

    for(i in 1:length(LEVELS))
    {
      SUB <- which(data$class==LEVELS[i])
      if(length(SUB)) { Ci[SUB] <- i }
    }
  }
  else
  { Ci <- rep(1,length(data$t)) }

  return(Ci)
}


#####
# class indicator matrix - [time,class]
get.class.mat <- function(data,LEVELS=levels(data$class))
{
  if("class" %in% names(data))
  {
    C <- sapply(LEVELS,function(lc){data$class==lc}) # (time,class)
    dim(C) <- c(nrow(data),length(LEVELS))
    colnames(C) <- LEVELS
  }
  else
  { C <- cbind(rep(1,length(data$t))) }

  return(C)
}


########
get.UERE.DOF <- function(x)
{
  N <- x$DOF
  if(is.null(N) || is.na(N)) { N <- 0 }
  return(N)
}


############
# get/set dimnames of UERE objects
classnames <- function(object)
{ rownames(object$UERE) }

"classnames<-" <- function(object,value)
{
  rownames(object$UERE) <- value
  rownames(object$DOF) <- value
  return(object)
}

"typenames<-" <- function(object,value)
{
  colnames(object$UERE) <- value
  colnames(object$DOF) <- value
  names(object$AICc) <- value
  names(object$Zsq) <- value
  names(object$VAR.Zsq) <- value
  names(object$N) <- value
  return(object)
}
