# periodogram class
new.periodogram <- methods::setClass("periodogram",representation("data.frame",info="list"))

# extend subset method
subset.periodogram <- function(x,...)
{
  info <- attr(x,"info")
  x <- data.frame(x)
  x <- subset.data.frame(x,...)
  x < - droplevels(x) # why is this here?
  new.periodogram(x,info=info)
}

`[.periodogram` <- function(x,...)
{
  info <- attr(x,"info")
  x <- data.frame(x)
  x <- "[.data.frame"(x,...)
  # if(class(x)[1]=="data.frame") { x <- new.periodogram(x,info=info) }
  x <- new.periodogram(x,info=info)
  return(x)
}


# periodogram wrapper
periodogram <- function(data,CTMM=NULL,dt=NULL,res.freq=1,res.time=1,fast=NULL,axes=c("x","y"))
{
  # listify for general purpose code
  data <- listify(data)
  if(!is.null(CTMM)) { CTMM <- listify(CTMM) }

  # detrend the mean
  for(i in 1:length(data))
  {
    # Some warning about S4 objects that I don't understand or care about...
    suppressWarnings(
    if(is.null(CTMM[[i]]))
    { data[[i]] <- set.telemetry(data[[i]],t(t(data[[i]][,axes,drop=FALSE]) - colMeans(data[[i]][,axes,drop=FALSE])),axes=axes) }
    else # use the better result if provided # this is unfinished
    {
      MU <- drift.mean(CTMM[[i]],data[[i]]$t) %*% CTMM[[i]]$mu
      data[[i]] <- set.telemetry(data[[i]],as.matrix(data[[i]][,axes,drop=FALSE]) - MU[,axes],axes=axes)
    }
    )
  }

  # default worst temporal resolution
  dts <- sapply(data,function(d) { stats::median(diff(d$t)) })
  if(is.null(dt)) { dt <- max(dts) }

  # default best sampling period
  Ts <- sapply(data,function(d) { last(d$t)-d$t[1] })
  T <- max(Ts)

  # round off grid number
  n <- round(T*res.freq/(dt/res.time)) + 1

  # intelligently select algorithm
  if(is.null(fast))
  {
    if(n < 10^4)
    { fast <- FALSE }
    else
    { fast <- 2 }
  }

  # fast time pre-gridding
  if(fast)
  {
    # transform times to aligned index + interpolation padding
    data <- lapply(data, function(d) { d$t <- d$t - d$t[1] ; d$t <- (d$t - grid.init(d$t,dt/res.time))/(dt/res.time) ; if(d$t[1]<0) {d$t <- d$t + 1} ; return(d) } )

    # all encompassing grid size, including interpolation buffer
    n <- max(sapply(data,function(d) { last(d$t) })) + 2*res.time
    n <- ceiling(n*res.freq)

    if(fast==2) { n <- composite(n) }
  }

  # frequency resolution
  df <- 1/(2*n*dt/res.time)

  # frequency grid
  f <- seq(from=df,by=df,length.out=n-1)

  # T=T*res.freq+dt/res.time
  if(fast)
  { lsp <- lapply(data,function(d) periodogram.fast(d,n=n,dn=res.time,axes=axes)) }
  else
  { lsp <- lapply(data,function(d) periodogram.slow(d,f=f,axes=axes)) }

  # DOF of one frequency in each LSP
  if(res.time==1)
  { SUM <- n-1 }
  else
  {
    LOW <- (f < 1/(2*dt))
    SUM <- sum(LOW)
  }
  DOF <- (dts/dt)*2*(sapply(data,function(d){length(data$t)})-1)/SUM

  # average the periodograms if there are multiple
  LSP <- rowSums( sapply(1:length(lsp),function(i){DOF[i]*lsp[[i]]$LSP}) )
  SSP <- rowSums( sapply(1:length(lsp),function(i){DOF[i]*lsp[[i]]$SSP}) )

  # DOF of one frequency in ave LSP
  DOF <- sum(DOF)

  LSP <- LSP/DOF
  SSP <- SSP/DOF

  result <- data.frame(f=f,LSP=LSP,SSP=SSP,DOF=DOF)

  if(res.time>1)
  {
    # toss high-frequency garbage from Lagrange interpolation
    result <- result[LOW,]

    # check for leakage from Lagrange interpolation high-frequency garbage
    if(any(is.nan(result$LSP)) || any(is.nan(result$SSP)) || any(result$LSP<0)|| any(result$SSP<0))
    { stop("Lagrange interpolation failed. Try a different value of res.time.") }
  }

  result <- new.periodogram(result,info=mean_info(data))
  attr(result,"info")$axes <- axes


  return(result)
}

# FFT periodogram code
periodogram.fast <- function(data,n=NULL,dn=NULL,axes=c("x","y"))
{
  t <- data$t
  t <- t + dn # starts at first grid element now (preceeded by interpolation buffer)

  SPAN <- (1-dn):dn
  OFFSET <- (dn-1) - SPAN

  # denomenator of all Lagrange polynomial terms
  DEN <- outer(SPAN,SPAN,"-")
  diag(DEN) <- rep(1,length(SPAN))
  DEN <- sapply(1:length(SPAN),function(i){ prod(DEN[i,]) })

  z <- get.telemetry(data,axes)
  COL <- ncol(z)

  W <- numeric(n)
  Z <- array(0,c(n,COL))

  # "spread" data across grid -- Lagrange interpolaton of the harmonics from Press & Rybicki
  for(i in 1:length(t))
  {
    FLOOR <- floor(t[i])
    p <- t[i] - FLOOR

    if(dn==1) # efficient code for simplest case
    {
      q <- 1-p
      LAGRANGE <- c(q,p)
    }
    else
    {
      dt <- p + OFFSET

      LAGRANGE <- prod(dt)/dt/DEN

      # avoid near coincidence singularity
      LAGRANGE[dn] <- prod(dt[-dn])/DEN[dn]
      LAGRANGE[dn+1] <- prod(dt[-(dn+1)])/DEN[dn+1]
    }

    IND <- FLOOR + SPAN
    W[IND] <- W[IND] + LAGRANGE
    Z[IND,] <- Z[IND,] + (LAGRANGE %o% z[i,])
  }

  # double padded fourier transforms
  W <- FFT(pad(W,2*n))
  Z <- FFT(rpad(Z,2*n))

  # double frequency argument
  # this code will fail for large n
  # W2 <- rep(W,4*n)[seq(1,4*n,2)]
  # same thing but without rep
  W2 <- c(W,W)[seq(1,4*n,2)]
  W2 <- Conj(W2)

  # LSP and SSP denominator
  DEN <- W[1]^2 - abs(W2)^2

  LSP <- W[1]*rowSums(abs(Z)^2) - W2*rowSums(Z^2)
  LSP <- Re(LSP/DEN/COL)

  SSP <- W[1]*abs(W)^2 - W2*W^2
  SSP <- Re(SSP/DEN)

  # Stop before Nyquist periodicity
  LSP <- LSP[2:n]
  SSP <- SSP[2:n]

  result <- data.frame(LSP=LSP,SSP=SSP)

  return(result)
}


# slow periodogram code
periodogram.slow <- function(data,f=NULL,axes=c("x","y"))
{
  t <- data$t
  z <- get.telemetry(data,axes)
  COL <- ncol(z)

  # double angle matrix
  theta <- (4 * pi) * (f %o% t)

  # LSP lag shifts
  tau <- atan( rowSums(sin(theta)) / rowSums(cos(theta)) ) / (4*pi*f)

  # lagged angle matrix
  theta <- (2 * pi) * ((f %o% t) - (f * tau))

  # trigonometric matrices
  COS <- cos(theta)
  SIN <- sin(theta)

  LSP <- rowSums((COS %*% z)^2)/rowSums(COS^2) + rowSums((SIN %*% z)^2)/rowSums(SIN^2)
  LSP <- LSP/(2*COL)

  # sampling schedule periodogram
  SSP <- rowSums(COS)^2/rowSums(COS^2) + rowSums(SIN)^2/rowSums(SIN^2)
  SSP <- SSP/2

  result <- data.frame(LSP=LSP,SSP=SSP)

  return(result)
}


# estimate the envelope of an oversampled/oscillatory periodigram
max_periodogram <- function(LSP,df=stats::median(diff(LSP$f)))
{
  #P <- log(P)
  dP <- diff(LSP$P)

  n <- length(dP)
  dP.left <- dP[1:(n-1)]
  dP.right <- dP[2:n]

  # first index of 3-index sequence containing a local maximum
  START <- which((dP.left>=0) & (dP.right<=0))

  # local maximum in periodogram (quadratic approximation)
  dP.left <- dP[START]/df
  dP.right <- dP[START+1]/df
  dP.mean <- (dP.left+dP.right)/2
  dP.curve <- (dP.left-dP.right)/df # this is negative

  df <- - dP.mean/dP.curve
  f <- LSP$f[START+1] + df
  P <- LSP$P[START+1] + dP.mean*df + (1/2)*dP.curve*df^2
  #P <- exp(P)

  result <- data.frame(f=f,P=P)

  return(result)
}


# plot periodograms
plot.periodogram <- function(x,max=FALSE,diagnostic=FALSE,col="black",transparency=0.25,grid=TRUE,...)
{
  DAY <- 1 %#% 'day' # daily periods

  # frequency in 1/seconds
  f <- x$f

  LSP <- data.frame(f=f,P=log(x$LSP))
  if(max) { LSP <- max_periodogram(LSP) }
  # re-scale to max at 0
  LSP$P <- LSP$P - max(LSP$P)

  at=labels=gcol=c()

  # tick specification function
  ticker <- function(time,div,name)
  {
    # fundamental
    at <<- c(at,time)
    labels <<- c(labels,name)
    gcol <<- c(gcol,grDevices::rgb(0.5,0.5,0.5,1))
    # harmonics
    DIV <- 2:div
    at <<- c(at,time/DIV)
    labels <<- c(labels,array(NA,length(DIV)))
    gcol <<- c(gcol,grDevices::rgb(0.5,0.5,0.5,1/DIV))
  }

  # yearly periods
  ticker(1 %#% 'year',12,"year")

  # lunar periods
  ticker(1 %#% 'month',29,"month")

  # diurnal periods
  ticker(DAY,24,"day")

  col <- malpha(col,alpha=((f[1]/f)^transparency))
  plot(1/LSP$f,LSP$P,log="x",xaxt="n",xlab="Period",ylab="Log Spectral Density",col=col,...)

  if(diagnostic)
  {
    LSP <- data.frame(f=f,P=log(x$SSP))
    if(max) { LSP <- max_periodogram(LSP) }
    LSP$P <- LSP$P - max(LSP$P)

    col <- malpha("red",alpha=((f[1]/f)^transparency))
    graphics::points(1/LSP$f,LSP$P,col=col,...)
  }

  graphics::axis(1,at=at,labels=labels) # tck can't be a vector apparently... :(
  if(grid){ graphics::abline(v=at,col=gcol) }
}
#methods::setMethod("plot",signature(x="periodogram",y="missing"), function(x,y,...) plot.periodogram(x,...))
#methods::setMethod("plot",signature(x="periodogram"), function(x,...) plot.periodogram(x,...))
