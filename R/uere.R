# store the UERE
"uere<-" <- function(data,value)
{
  UERE <- value

  DROP <- FALSE
  # promote to list
  if(class(data)=="telemetry" || class(data)=="data.frame")
  {
    data <- list(data)
    DROP <- TRUE
  }

  for(i in 1:length(data))
  {
    attr(data[[i]],"info")$UERE <- UERE

    NAMES <- names(data[[i]])

    if(is.null(UERE)) # delete UERE information
    {
      data[[i]]$HERE <- NULL
      data[[i]]$VERE <- NULL
    }
    else # store UERE information
    {
      if("HDOP" %in% NAMES) { data[[i]]$HERE <- UERE*data[[i]]$HDOP }
      if("VDOP" %in% NAMES) { data[[i]]$VERE <- UERE*data[[i]]$VDOP }
      if(!any(c("HDOP","VDOP") %in% NAMES)){ stop("Failed to import DOP column from Movebank file.") }
    }
  }

  if(DROP) { data <- data[[1]] }

  # demote to data.frame
  return(data)
}

# return the UERE from set data or calculate it from calibration data
uere <- function(data,diagnostic=FALSE)
{
  if(class(data)=="telemetry" || class(data)=="data.frame")
  {
    # return present UERE
    UERE <- attr(data,"info")$UERE
    if(!is.null(UERE)) { return(UERE) }

    # promote to list of calibration data and proceed
    data <- list(data)
  }

  x <- list()
  y <- list()
  # calculate mean and residuals
  for(i in 1:length(data))
  {
    if(is.null(data[[i]]$HDOP))
    { stop("Failed to import GPS.HDOP column from Movebank file.") }

    w <- 1/data[[i]]$HDOP^2
    mu <- c( w %*% data[[i]]$x , w %*% data[[i]]$y )/sum(w)

    # detrend data
    data[[i]]$x <- data[[i]]$x - mu[1]
    data[[i]]$y <- data[[i]]$y - mu[2]

    w <- 1/data[[i]]$HDOP
    y[[i]] <- w * data[[i]]$x
    x[[i]] <- w * data[[i]]$y
  }
  x <- unlist(x)
  y <- unlist(y)

  # total degrees of freedom
  n <- length(x) - length(data)

  # residuals and UERE
  UERE <- sqrt(sum(x^2+y^2)/n)

  if(diagnostic)
  {
    r <- c(x,y) / sqrt(UERE/2)
    CI <- chisq.ci(1,DOF=2*n)

    KDE <- stats::density(r,bw="SJ")
    graphics::hist(r,breaks="scott",freq=FALSE,main="residual distribution",xlab="standardized error",ylim=c(0,max(KDE$y,1/sqrt(2*pi*CI[1]))))
    graphics::lines(KDE)
    DNORM <- Vectorize(function(x){ exp(-x^2/2/CI[2])/sqrt(2*pi*CI[2]) })
    graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red",lwd=2)
    DNORM <- Vectorize(function(x){ exp(-x^2/2/CI[1])/sqrt(2*pi*CI[1]) })
    graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red")
    DNORM <- Vectorize(function(x){ exp(-x^2/2/CI[3])/sqrt(2*pi*CI[3]) })
    graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red")

    return(sqrt(chisq.ci(UERE^2,DOF=2*n)))
  }
  else
  { return( UERE ) }
}
