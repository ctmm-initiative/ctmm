# raster::extent keeps working differently and S4 methods are impossible to debug

# generic for matrix and data.frame
extent.default <- function(x,...)
{
  x <- as.data.frame(x)
  RANGE <- data.frame(row.names=c('min','max'))

  for(COL in colnames(x))
  {
    if(class(x[[COL]])[1] %in% c("numeric","POSIXct"))
    {
      if(COL=="longitude") # use circular statistics
      { RANGE[[COL]] <- quantile.longitude(x[[COL]],na.rm=TRUE) }
      else
      {
        # RANGE[[COL]] <- stats::quantile(x[[COL]],probs=probs,na.rm=TRUE)
        RANGE[[COL]] <- range(x[[COL]],na.rm=TRUE,finite=TRUE)
      }
    }
  }

  return(RANGE)
}
extent.telemetry <- extent.default
extent.data.frame <- extent.default
extent.matrix <- extent.default

as.matrix.Extent <- function(x,...)
{
  x <- cbind(x=c(x@xmin,x@xmax),y=c(x@ymin,x@ymax))
  rownames(x) <- c('min','max')
  return(x)
}


# what is the range of a list of objects of any type
extent.list <- function(x,...)
{
  # list of element ranges
  RANGE <- lapply(x,function(y){extent(y,...)})

  # # shared columns only
  # COLS <- lapply(RANGE,colnames)
  # COLS <- Reduce(intersect,COLS)
  # RANGE <- lapply(RANGE,function(R){R[,COLS]})

  # all columns
  COLS <- lapply(RANGE,colnames)
  COLS <- unlist(COLS)
  COLS <- unique(COLS)
  for(i in 1:length(x))
  {
    # fill missing columns with NAs
    NIN <- COLS %nin% colnames(RANGE[[i]])
    if(any(NIN)) { RANGE[[i]][,COLS[NIN]] <- NA }
    # sort columns
    RANGE[[i]] <- RANGE[[i]][,COLS,drop=FALSE]
  }

  # concatenate ranges
  RANGE <- Reduce(rbind,RANGE)
  # final range x & y
  RANGE <- extent(RANGE,level=1)
  # level has already been applied
  return(RANGE)
}


# range of Gaussian contours
extent.ctmm <- function(x,level=0.95,level.UD=0.95,complete=FALSE,...)
{
  axes <- x$axes

  level <- max(level)
  level.UD <- max(level.UD)

  if(is.null(x$mu)) { stop("This model has no mean location. Try ctmm.guess.") }

  alpha <- 1 - level
  alpha.UD <- 1 - level.UD
  # standard Gaussian quantile stuff
  if(length(axes)==2)
  { z <- sqrt(-2*log(alpha.UD)) }
  else if(length(axes)==1)
  { z <- stats::qnorm(1-alpha.UD/2) }

  # proportionality constants for outer CIs
  sigma <- x$sigma

  # capture outer contour if present
  const <- 1
  if(is.na(level)) { alpha <- 1 } # will use ML contour
  if("COV" %in% names(x))
  {
    K <- length(x$tau)
    const <- confint.ctmm(x,alpha)["area",3]/ifelse(length(axes)==2,sqrt(det(sigma)),sigma)
  }

  if(complete && length(axes)==2)
  {
    sigma <- const*sigma@.Data
    RANGE <- ellipsograph(x$mu,sigma,level=level.UD,PLOT=FALSE)
    LL <- project(RANGE,from=projection(x))
    colnames(LL) <- GEO
    RANGE <- cbind(RANGE,LL)
    RANGE <- extent(RANGE)
  }
  else # x,y only - faster
  {
    buff <- z*sqrt(const*diag(sigma))

    if(length(axes)==2)
    {
      X <- x$mu[1] + c(-1,1)*buff[1]
      Y <- x$mu[2] + c(-1,1)*buff[2]

      RANGE <- data.frame(x=X,y=Y)
    }
    else
    {
      Z <- x$mu[1] + c(-1,1)*buff

      RANGE <- data.frame(Z)
      colnames(RANGE) <- axes
    }
  }

  row.names(RANGE) <- c("min","max")
  return(RANGE)
}


# range of UD contours
extent.UD <- function(x,level=0.95,level.UD=0.95,complete=FALSE,...)
{
  level <- max(level)
  level.UD <- max(level.UD)

  PROJ <- attr(x,"info")$projection

  if(level.UD==1 || (!is.na(level) && level==1 && !is.null(x$DOF.area)))
  {
    RANGE <- data.frame(x=c(-1,1),y=c(-1,1))*Inf
    row.names(RANGE) <- c("min","max")
    return(RANGE)
  }

  # capture ML contour
  if(is.null(x$DOF.area) || any(x$DOF.area<=.Machine$double.eps) || is.na(level))
  { P <- max(level.UD) }
  else # capture outer contour
  { P <- CI.UD(x,max(level.UD),max(level),P=TRUE)[3] }

  # too coarse fix
  P <- max(P,min(x$CDF))

  if(complete)
  {
    # list of contours (multiple if multiple regions/modes)
    R <- grDevices::contourLines(x=x$r$x,y=x$r$y,z=x$CDF,levels=P)
    R <- lapply(R,function(L){cbind(L$x,L$y)})
    R <- do.call(rbind,R)

    X <- range(R[,1])
    Y <- range(R[,2])

    RANGE <- data.frame(x=X,y=Y)
    row.names(RANGE) <- c("min","max")

    R <- project(R,from=PROJ)
    RANGE$longitude <- quantile.longitude(R[,1],probs=c(0,1)) # circular extent
    RANGE$latitude <- range(R[,2])
  }
  else # faster code
  {
    # do we extend this far?
    MAT <- (x$CDF <= P)

    X <- apply(MAT,1,any)
    Y <- apply(MAT,2,any)
    rm(MAT)

    # now indices
    X <- which(X)
    Y <- which(Y)

    # range indices
    X <- c(X[1],last(X))
    Y <- c(Y[1],last(Y))

    # location ranges
    X <- x$r$x[X]
    Y <- x$r$y[Y]

    RANGE <- data.frame(x=X,y=Y)
    row.names(RANGE) <- c("min","max")
  }

  return(RANGE)
}


# extent of a variogram object
extent.variogram <- function(x,level=0.95,threshold=2,...)
{
  level <- max(level)

  alpha <- 1-level
  max.lag <- last(x$lag)
  ACF <- !is.null(attr(x,"info")$ACF)

  if(!ACF) # SVF plot
  {
    min.SVF <- 0
    # maximum CI on SVF
    max.SVF <-  max(x$SVF * CI.upper(x$DOF,min(alpha)))
    # limit plot range to threshold * max SVF point estimate (otherwise hard to see)
    max.SVF <- min(max.SVF,threshold*max(x$SVF))
  }
  else # ACF plot
  {
    min.SVF <- min(0,min(x$SVF))

    max.SVF <- max(0,max(x$SVF[-1]))
  }

  EXT <- data.frame(x=c(0,max.lag),y=c(min.SVF,max.SVF))
  rownames(EXT) <- c("min","max")

  return(EXT)
}


# intersection of extents
min.extent <- function(...,na.rm=FALSE)
{
  LIST <- list(...)

  if(!length(LIST))
  {
    INF <- rbind(min=-Inf,max=Inf)
    INF <- cbind(INF,INF)
    return(INF)
  }

  MIN <- LIST[[1]]
  MIN[1,1] <- max(sapply(LIST,function(L){L[1,1]}),na.rm=na.rm)
  MIN[1,2] <- max(sapply(LIST,function(L){L[1,2]}),na.rm=na.rm)
  MIN[2,1] <- min(sapply(LIST,function(L){L[2,1]}),na.rm=na.rm)
  MIN[2,2] <- min(sapply(LIST,function(L){L[2,2]}),na.rm=na.rm)

  # fix for no intersection
  MIN[,1] <- range(MIN[,1])
  MIN[,2] <- range(MIN[,2])

  return(MIN)
}

