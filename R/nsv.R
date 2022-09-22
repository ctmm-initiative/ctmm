# non-stationary variogram
nsv <- function(data,CTMM=NULL,cycle=FALSE,axes=c('x','y'),dt=NULL)
{
  GRID <- pregridder(data$t,dt=dt)
  dt <- GRID$dt
  t0 <- GRID$t0
  n <- round((last(data$t)-t0)/dt)
  t <- t0 + (0:n)*dt
  ts <- as.POSIXct(t,origin=EPOCH,tz="UTC")

  # lazy linear interpolation
  if(is.null(CTMM))
  { CTMM <- ctmm.fit(data,ctmm(sigma=1,range=FALSE,tau=Inf,axes=axes)) }
  data <- predict(data,CTMM,t=t)

  VAR <- 0
  for(a in axes) { VAR <- VAR + outer(data[[a]],FUN='-')^2 }
  VAR <- VAR/(2*length(axes))

  VAR <- list(SVF=VAR,t=t)
  class(VAR) <- "NSV"
  return(VAR)
}

plot.nsv <- function(x,...)
{
  SVF <- x$SVF
  t <- x$t

  zlim <- c(0,max(SVF))
  graphics::image(t,t,SVF,zlim=zlim,asp=1,axes=FALSE,...)
  graphics::lines(t[1]*c(1,1),last(t)*c(1,1),...)
  graphics::lines(c(last(x),t[1]),c(last(t),t[1]),...)
}
