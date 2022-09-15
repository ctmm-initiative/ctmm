# location difference vector
# data is a list of two telemetry objects
# CTMM is a list of two ctmm fit objects corresponding to data
difference <- function(data,CTMM,dt=NULL,...)
{
  check.projections(data)
  INFO <- mean.info(data)

  t1 <- max(data[[1]]$t[1],data[[2]]$t[1])
  t2 <- min(last(data[[1]]$t),last(data[[2]]$t))

  if(is.null(dt)) { dt <- min( stats::median(diff(data[[1]]$t)) , stats::median(diff(data[[2]]$t)) ) }
  dt <- dt/10 # search at finer scale
  n <- (t2-t1)/dt
  n <- ceiling(n)
  dt <- (t2-t1)/n
  t <- seq(t1,t2,length.out=n+1)
  n <- length(t)

  # predict over fine grid
  data[[1]] <- predict(data[[1]],CTMM[[1]],t=t)
  data[[2]] <- predict(data[[2]],CTMM[[2]],t=t)

  axes <- CTMM$axes
  for(z in axes) { data[[1]][[z]] <- data[[1]][[z]] - data[[2]][[z]] }

  VAR <- DOP.LIST$horizontal$VAR
  VARS <- c(VAR,DOP.LIST$horizontal$COV)
  for(v in VARS) { data[[1]][[v]] <- data[[1]][[v]] + data[[2]][[v]] }

  data <- data[[1]][,c('t',axes,VARS)]

  # now find the canonical times that minimize error
  # initial guess
  CANON <- rep(FALSE,nrow(data))
  if(data[[VAR]][1]<data[[VAR]][2]) { CANON[1] <- TRUE }
  for(i in 2:(n-1)) { if(data[[VAR]][i]<min(data[[VAR]][c(i-1,i+1)])) { CANON[i] <- TRUE } }
  if(data[[VAR]][n]<data[[VAR]][n-1]) { CANON[n] <- TRUE }

  IND <- which(CANON)
  m <- length(IND)
  DATA <- array(0,c(m,ncol(data)))
  colnames(DATA) <- colnames(data)
  for(j in 1:m)
  {
    i <- IND[j]
    # end points should be exact
    if(i==1 || i==n)
    { DAT <- data[i,] }
    else # in between interpolate quadratically
    {
      DAT <- data[i+(-1):1,]
      # coefficients
      b0 <- DAT[2,]
      DAT <- t( t(DAT) - DAT[2,] )
      b1 <- (DAT[3,]-DAT[1,])/2
      b2 <- (DAT[3,]+DAT[1,])/2
      # minimum variance index
      x <- -b1[VAR]/(2*b2[VAR])
      x <- nant(x,sign(b1[VAR]))
      # interpolate
      DAT <- b0 + b1*x + b2*x^2
      # fix time
      DAT['t'] <- data$t[i] + x*dt
      # !!! UNFINISHED
    }
    DATA[j,] <- DAT
  } # end for
  rm(data)

  # fix repeats for constant VAR...
  SUB <- c(1,diff(DATA$t)) > 0
  DATA <- DATA[SUB,]

  # make this a telemetry object
  DATA <- new.telemetry(DATA,info=INFO)
  # make sure UERE is correct

  return(DATA)
}
