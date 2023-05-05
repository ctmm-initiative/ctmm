# attach lunar, solar, and seasonal phases to telemetry object
annotate <- function(object,by="all",cores=1,...)
{
  if(class(object)[1]=="list") { return(plapply(object,annotate,cores=cores,fast=FALSE)) }
  # else one at a time below

  if(class(by)[1]=="list")
  {
    # for the future
  }

  #
  if(class(by)[1]=="data.frame")
  {
    # FOR THE FUTURE
  }

  if(by=='all' || 'moon' %in% by)
  {
    ADD <- suncalc::getMoonIllumination(object$timestamp,keep="fraction")$fraction
    object$moonlight <- ADD # 0:1 fraction
  }

  data <- data.frame(date=object$timestamp,lat=object$latitude,lon=object$longitude)
  if(by=='all' || 'sun' %in% by)
  {
    ADD <- suncalc::getSunlightPosition(data=data,keep="altitude")$altitude
    object$sunlight <- sin(ADD) # -1:1 relative flux when positive
  }

  if(by=='all' || 'season' %in% by)
  {
    data$date <- as.Date(object$timestamp)
    ADD <- suncalc::getSunlightTimes(data=data,keep=c("sunrise","sunset"))
    ADD$sunrise <- as.numeric(ADD$sunrise)
    ADD$sunset <- as.numeric(ADD$sunset)
    ADD <- "day" %#% (ADD$sunset - ADD$sunrise) # fraction of mean day that sun is up
    object$season <- 2*ADD - 1 # roughly centered on 0 and possibly -1:1
  }

  # need to find an R package that can calculate solstice & equinox times
  if(by=='all' || 'tropic' %in% by)
  {
    # current year
    data <- format(as.Date(object$timestamp,format="%d/%m/%Y",tz="GMT"),"%Y")
    data <- as.numeric(data)
    # enclosing year ends
    data <- cbind( before=as.POSIXct(paste0(data,"/01/01"),format="%Y/%m/%d",tz="GMT") , after=as.POSIXct(paste0(data+1,"/01/01"),format="%Y/%m/%d",tz="GMT") )
    data <- as.numeric(data)
    dim(data) <- c(nrow(object),2)
    object$tropic <- (object$t-data[,1])/(data[,2]-data[,1])
  }

  if(by %in% c('switch','sundial'))
  {
    STUFF <- get.sundial(object)
    for(col in names(STUFF)) { object[[col]] <- STUFF[[col]] }
  }

  return(object)
}


################
# color by time/phase
color <- function(object,by="time",col.fn=NULL,alpha=1,dt=NULL,cores=1,...)
{
  object <- listify(object)

  # color wheel placements
  if(by=="individual") { index <- color.individual(object,cores=cores,...) }

  # default color functions
  if(is.null(col.fn))
  {
    if(by=="individual") { col.fn <- function(i,alpha){grDevices::hsv(i,alpha=alpha)} }
    else { col.fn <- function(i,alpha){grDevices::rgb(i,0,1-i,alpha)} }
  }

  if(class(object[[1]])[1] %in% c('UD','CTMM')) # color densities by individual only
  {
    if(by!='individual') { stop("Only by='individual' supported by object class.") }

    # apply color function
    index <- col.fn(index,alpha) # index already calculated above
  }
  else # telemetry or data.frame
  {
    BYS <- c("individual","time","sun","moon","season","tropic")
    # check for annotations if necessary
    check <- function(cby,column) { if(by==cby && column %nin% names(object[[1]])) stop("Data are not annotated.") }
    check('sun','sunlight')
    check('moon','moonlight')
    check('season','season')
    check('tropic','tropic')
    if(by %nin% BYS && by %nin% names(object[[1]])) { stop("Data are not annotated.") }

    # determine index
    if(by %in% c("individual","tropic")) # cyclic scale
    {
      if(by=="individual") # numeric index calculated above
      { index <- lapply(1:length(object), function(i){ rep(index[i],length(object[[i]]$t)) }) }
      else if(by=="tropic")
      { index <- lapply(object,function(o){o$tropic}) }
    }
    else # relative scale (normalized)
    {
      if(by=="time")
      { index <- lapply(object,function(o){o$t}) }
      else if(by=="sun")
      { index <- lapply(object,function(o){o$sunlight}) }
      else if(by=="moon")
      { index <- lapply(object,function(o){o$moonlight}) }
      else if(by=="season")
      { index <- lapply(object,function(o){o$season}) }
      else # custom column
      { index <- lapply(object,function(o){o[[by]]}) }

      # scale index to 0:1
      MIN <- sapply(index,function(i){min(i)})
      MIN <- min(MIN)
      MAX <- sapply(index,function(i){max(i)})
      MAX <- max(MAX)
      index <- lapply(index,function(i){(i-MIN)/(MAX-MIN)})
    }

    # calculate alpha channel
    ALPHA <- lapply(object,function(o){ diff(o$t) })
    if(is.null(dt)) { dt <- stats::median( sapply(ALPHA,function(a){ stats::median(a) }) ) }
    ALPHA <- lapply(ALPHA,function(a){ a <- pmin(c(Inf,a),c(a,Inf)); a <- alpha*clamp(a/dt); return(a) })

    # apply color function
    index <- lapply(1:length(index),function(i){ col.fn(index[[i]],ALPHA[[i]]) })

    if(length(index)==1) { index <- index[[1]] }
  }

  names(index) <- names(object)
  return(index)
}


simplify.color <- function(object)
{
  if(class(object)[1]=="list")
  {
    NAMES <- names(object)
    for(i in 1:length(object))
    {
      col <- object[[i]]
      col <- grDevices::col2rgb(col,alpha=TRUE)
      # opacity weighted average
      w <- col['alpha',]
      w <- w/sum(w)
      col <- col[1:3,] %*% w
      col <- grDevices::rgb(red=col[1],green=col[2],blue=col[3],maxColorValue=255)
      object[[i]] <- col
    }
    object <- unlist(object)
    names(object) <- NAMES
  }

  return(object)
}


# multiply alpha
malpha <- function(col,alpha=1)
{
  n <- max(length(col),length(alpha))
  col <- array(col,n)
  alpha <- array(alpha,n)

  col <- grDevices::col2rgb(col,alpha=TRUE)
  col[4,] <- col[4,] * alpha
  col <- grDevices::rgb(col[1,],col[2,],col[3,],col[4,],maxColorValue=255)
  return(col)
}

grad.white <- function(col)
{
  col <- c(grDevices::col2rgb(col))
  col <- c(grDevices::rgb2hsv(col[1],col[2],col[3]))
  white <- (255:0)/255
  col <- grDevices::hsv(col[1],white,col[3])
  return(col)
}


##############################
# COLOR BY INDIVIDUAL
# a greedy algorithm
color.individual <- function(object,cores=1,...)
{
  object <- listify(object)
  n <- length(object)
  if(n==1) { return(0) }
  if(n==2) { return(c(0,1/2)) }
  if(n==3) { return(c(0,1/3,2/3)) }

  CLASS <- class(object[[1]])[1]
  # could update this to robust overlap calculation
  if(CLASS %in% c('telemetry','data.frame'))
  { OVER <- plapply(object,function(o){ ctmm.fit(o,ctmm(isotropic=(nrow(o)<3)),method='ML',COV=FALSE) },cores=cores) }
  else if(CLASS %in% c('ctmm','UD'))
  { OVER <- object }
  OVER <- overlap(OVER,debias=FALSE,COV=FALSE)$CI[,,2]
  diag(OVER) <- 0 # no self interactions

  # who has the worst spatial overlap
  ORDER <- apply(OVER,1,function(o){mean(sort(o,decreasing=TRUE)[1:2])}) # need 2 values to break ties... full sort is wasteful
  ORDER <- order(ORDER,decreasing=TRUE)

  POS <- ORDER[1:2] # current order in color wheel
  N <- length(POS) # individuals currently on color wheel

  # insert individuals one-by-one
  while(N < n)
  {
    MAT <- outer(1:(N+1),1:(N+1),'-') # color distance matrix
    MAT <- abs(MAT)
    MAT <- ifelse(MAT>(N+1)/2,N+1-MAT,MAT) # cyclic property of colors
    MAT <- (2*pi/2)/(N+1) * MAT # map 0:1 to 0:pi for angle overlap
    MAT <- cos(MAT) # color overlap matrix

    MAX <- rep(NA,N) # to store worst values
    for(i in 1:N) # cycle through gaps
    {
      # ? what is the worst combined overlap if we insert individual N+1 in gap i ?
      TRY <- c(POS[1:i],ORDER[N+1]) # attempted sorting
      if(i<N) { TRY <- c(TRY,POS[(i+1):N]) }
      TRY <- OVER[TRY,TRY] # corresponding overlap matrix
      TRY <- TRY * MAT # combined overlap
      MAX[i] <- max(TRY) # worst combined overlap
    }
    # best worst insertion
    MIN <- which.min(MAX)

    # insert individual N+1
    NEW <- c(POS[1:MIN],ORDER[N+1])
    if(MIN<N) { NEW <- c(NEW,POS[(MIN+1):N]) }

    POS <- NEW
    N <- length(POS)
  }

  theta <- (1:N - 1)/N
  theta[POS] <- theta
  return(theta)
}
