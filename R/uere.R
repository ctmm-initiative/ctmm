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

  n <- 0
  mu <- list()
  r <- list()
  # calculate mean and residuals
  for(i in 1:length(data))
  {
    if(is.null(data[[i]]$HDOP))
    { stop("Failed to import GPS.HDOP column from Movebank file.") }

    w <- 2/data[[i]]$HDOP^2
    n <- n + 2*length(w)
    mu[[i]] <- c( w %*% data[[i]]$x , w %*% data[[i]]$y )/sum(w)
    r[[i]] <- c( sqrt(w) * (data[[i]]$x - mu[[i]][1]) , sqrt(w) * (data[[i]]$y - mu[[i]][2]) )
  }
  r <- unlist(r)

  # total degrees of freedom
  n <- n - 2*length(data)

  # residuals and UERE
  UERE <- sqrt(sum(r^2/n))
  CI <- sqrt(chisq.ci(UERE^2,DOF=n))

  if(diagnostic)
  {
    graphics::hist(r,breaks="scott",freq=FALSE,main="residual distribution",xlab="normalized residual")
    graphics::lines(stats::density(r,bw="SJ"))
    DNORM <- Vectorize(function(x){ stats::dnorm(x,sd=UERE) })
    graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red",lwd=2)
    DNORM <- Vectorize(function(x){ stats::dnorm(x,sd=CI[1]) })
    graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red")
    DNORM <- Vectorize(function(x){ stats::dnorm(x,sd=CI[3]) })
    graphics::curve(DNORM,from=min(r),to=max(r),n=1001,add=TRUE,col="red")

    return(CI)
  }
  else
  { return( UERE ) }
}
