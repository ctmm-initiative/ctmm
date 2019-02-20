# attach lunar, solar, and seasonal phases to telemetry object
annotate <- function(object,by="all",cores=1,...)
{
  if(class(object)=="list") { return(plapply(object,annotate,cores=cores,fast=FALSE)) }
  # else one at a time below

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

  return(object)
}


################
# color by time/phase
color <- function(object,by="time",col.fn=NULL,alpha=1,dt=NULL,...)
{
  object <- listify(object)

  BYS <- c("individual","time","sun","moon","season","tropic")
  # check for annotations if necessary
  check <- function(cby,column) { if(by==cby && column %nin% names(object[[1]])) stop("Data are not annotated.") }
  check('sun','sunlight')
  check('moon','moonlight')
  check('season','season')
  check('tropic','tropic')
  if(by %nin% BYS && by %nin% names(object[[1]])) { stop("Data are not annotated.") }

  if(is.null(col.fn))
  {
    if(by=="individual") { col.fn <- function(i,alpha){grDevices::hsv(i,alpha=alpha)} }
    else { col.fn <- function(i,alpha){grDevices::rgb(i,0,1-i,alpha)} }
  }

  # determine index
  if(by %in% c("individual","tropic")) # cyclic scale
  {
    if(by=="individual")
    {
      index <- color.individual(object,...)
      index <- lapply(1:length(object), function(i){ rep(index[i],length(object[[i]]$t)) })
    }
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

  names(index) <- names(object)
  return(index)
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


##############################
# COLOR BY INDIVIDUAL
# precision: fraction of numerical precision to target---error seems larger than it should be (1/2)
# maxit: maximum number of iterations to attempt when equilibrating
color.individual <- function(object,precision=1/4,STEP=1,maxit=Inf,...)
{
  object <- listify(object)
  n <- length(object)
  if(n==1) { return(0) }
  if(n==2) { return(c(0,1/2)) }

  STEP <- 1 # initial equilibration step size (decays harmonically)
  tol <- .Machine$double.eps^precision

  CLASS <- class(object[[1]])
  # could update this to robust overlap calculation
  if(CLASS %in% c('telemetry','data.frame'))
  { OVER <- lapply(object,function(o){ ctmm.fit(o,ctmm(isotropic=(nrow(o)<3)),method='ML') }) }
  else if(CLASS %in% c('ctmm','UD'))
  { OVER <- object }
  OVER <- overlap(OVER,debias=FALSE,COV=FALSE)[,,2]

  # supply small minimum repulsion/overlap for some initial color separation
  repulsion <- min(1/length(object),1/10)
  repulsion <- max(repulsion-min(OVER),0) # don't necessarily need it
  OVER <- (OVER + repulsion)/(1 + repulsion) # O(repulsion^2) change

  # who has the most overlap with everyone
  SOVER <- rowSums(OVER)
  ORDER <- order(SOVER,decreasing=TRUE)

  # initially greedy energy minimization algorithm matching highest overlaps to largest color distances
  r <- array(Inf,c(n,2)) # current locations
  colnames(r) <- c('x','y')
  rownames(r) <- names(object)

  N <- 2 # individuals currently on color wheel
  # initial solution
  r[ORDER[1:2],] <- array(c(1,-1,0,0),c(2,2))

  # find best gap to insert next color
  insert <- function()
  {
    # current angles
    theta <- atan2(r[ORDER[1:N],'y'],r[ORDER[1:N],'x']) # -pi:pi
    theta <- sort(theta)
    # midpoints between angles
    M <- (theta[-1]+theta[-N])/2
    M[N] <- ( (theta[1]+2*pi) + theta[N] )/2
    # midpoint locations to evaluate
    M <- cbind(x=cos(M),y=sin(M))
    # difference vectors (mid,other)
    dX <- outer(M[,'x'],r[ORDER[1:N],'x'],'-')
    dY <- outer(M[,'y'],r[ORDER[1:N],'y'],'-')
    dr <- sqrt(dX^2+dY^2)
    U <- -OVER[ORDER[N+1],ORDER[1:N]]*log(dr) # potential energies of next individual
    U <- rowSums(U) # total potential energies of next individual
    MIN <- which.min(U) # minimum energy midpoint
    return(M[MIN,])
  }

  # adjust spacings to minimize energy
  equilibrate <- function()
  {
    ORDER <- ORDER[1:N]
    r <- r[ORDER,]
    OVER <- OVER[ORDER,ORDER]

    i <- 0
    ERROR <- Inf
    tol <- N*tol
    while(i<=maxit && ERROR>=tol)
    {
      OLD <- diff(sort(atan2(r[,'y'],r[,'x'])))
      dX <- outer(r[,'x'],r[,'x'],'-')
      dY <- outer(r[,'y'],r[,'y'],'-')
      dr2 <- dX^2+dY^2
      diag(dr2) <- Inf # zero self interaction
      FX <- OVER * dX / dr2 # 1/r force on first index from second index
      FY <- OVER * dY / dr2
      FX <- rowSums(FX) # net force on index
      FY <- rowSums(FY) # net force on index
      Fr <- cbind(x=FX,y=FY) # net force vectors
      MAX <- sqrt(max(rowSums(Fr^2))) # max force, used to scale time step
      r <- r + (STEP/(i+1)/MAX)*Fr # run forward in time
      r <- r / sqrt(rowSums(r^2)) # pull back to color wheel
      r[1,] <- c(1,0) # fix most overlapping individual to first color

      ERROR <- diff(sort(atan2(r[,'y'],r[,'x']))) - OLD
      ERROR <- max(abs(ERROR))
      i <- i + 1
    }
    return(r)
  }

  # insert individuals one-by-one
  while(N < n)
  {
    r[ORDER[N+1],] <- insert() # insert next charge into best gap
    N <- N + 1
    # plot(r[ORDER[1:N],],xlim=c(-1,1),ylim=c(-1,1),asp=1,cex=2*SOVER[ORDER[1:N]]/max(SOVER)); points(0,0,col='red')
    r[ORDER[1:N],] <- equilibrate() # equilibrate current charges
    # plot(r[ORDER[1:N],],xlim=c(-1,1),ylim=c(-1,1),asp=1,cex=2*SOVER[ORDER[1:N]]/max(SOVER)); points(0,0,col='red')
  }

  # normalize interaction strengths to equispace colors
  # plot(r[ORDER[1:N],],xlim=c(-1,1),ylim=c(-1,1),asp=1,cex=2*SOVER[ORDER[1:N]]/max(SOVER)); points(0,0,col='red')
  while(min(OVER)<0.5)
  {
    OVER <- (OVER + sqrt(1/N))/(1 + sqrt(1/N)) # O(1/N) change
    r <- equilibrate()
    # plot(r[ORDER[1:N],],xlim=c(-1,1),ylim=c(-1,1),asp=1,cex=2*SOVER[ORDER[1:N]]/max(SOVER)); points(0,0,col='red')
    # title(min(OVER))
  }

  # format results to unit interval
  theta <- atan2(r[,'y'],r[,'x'])/(2*pi)
  theta <- theta - min(theta) # no longer necessary

  # make equispaced (should be close)
  # ORDER <- order(theta)
  # theta <- (ORDER-1)/N

  return(theta)
}
