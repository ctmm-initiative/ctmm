# function for examining sampling schedule

dt.plot <- function(data,...)
{
  # sampling intervals
  dt <- listify(data)
  dt <- lapply(dt,function(d){diff(d$t)})
  dt <- unlist(dt)

  #dt <- diff(data$t) # sampling intervals
  dt <- dt[dt>0] # remove for log scale
  dt <- sort(dt) # Beth's idea
  RANGE <- c(dt[1],last(dt))

  col <- grDevices::rgb(1,1,1,0)
  x <- c(1,length(dt))
  y <- c(dt[1],last(dt))
  plot(x,y,log='y',yaxt='n',ylab='Time Intervals',xlab='Sorted Index',col=col,...)
  # title("Time Intervals")

  TIME <- c('sec','min','hour','day','month','year')
  MULT <- list()
  MULT[['sec']] <- c(1,5,15,30)
  MULT[['min']] <- MULT[['sec']]
  MULT[['hour']] <- c(1,6,12)
  MULT[['day']] <- c(1,7,14)
  MULT[['month']] <- c(1,3,6)
  for(i in 1:length(TIME))
  {
    u <- TIME[i]
    t <- 1 %#% u

    graphics::axis(side=2,at=t,label=u,las=2,tick=FALSE,mgp=c(0,0.2,0))

    if(i<length(TIME))
    {
      MAX <- (1 %#% TIME[i+1])/t

      for(j in 1:length(MULT[[i]]))
      {
        DIVS <- round( MAX/MULT[[i]][j] )
        for(k in 1:DIVS)
        { graphics::abline(h=k*MULT[[i]][j]*t,col=grDevices::grey(1-1/DIVS)) }
      }
    }

    graphics::abline(h=t)
  }

  graphics::points(dt,...)
}
