uere <- function(data,diagnostic=FALSE)
{
  # promote to list of 
  if(class(data)=="telemetry" || class(data)=="data.frame") { data <- list(data)  }
  
  n <- 0
  mu <- list()
  r <- list()
  # calculate mean and residuals
  for(i in 1:length(data))
  {
    if(is.null(data[[i]]$HDOP)) { stop("Failed to import GPS.HDOP column from Movebank file || missing HDOP column in telemetry object.") }
    
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