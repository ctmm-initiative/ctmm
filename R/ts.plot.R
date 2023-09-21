ts.plot <- function(data,col='black',ylab="",...)
{
  if('timestamp' %in% names(data))
  {
    t <- data$timestamp
    xlab <- "Time"
  }
  else
  {
    t <- data$t
    t <- t-t[1]
    UNITS <- unit(t,"time")
    # ...
  }

  CI <- data[,NAMES.CI]

  MAX <- max(CI[,'est'])
  MIN <- min(CI[,'est'])
  DIFF <- MAX-MIN
  ylim <- c( max(MIN-DIFF,min(CI[,'low'])) , min(MAX+DIFF,max(CI[,'high'])) )

  graphics::plot(t,CI[,'est'],type="l",col=col,ylim=ylim,xlab=xlab,ylab=ylab,...)
  graphics::points(t,CI[,'low'],type="l",col=malpha(col,alpha=0.1),...)
  graphics::points(t,CI[,'high'],type="l",col=malpha(col,alpha=0.1),...)
  #graphics::polygon(c(t,rev(t)),c(CI[,'low'],CI[,'high']),col=malpha(col,alpha=0.1),border=NA,...)

  # ERROR BAR PLOT
  # plot(t,CI[,2],pch=19,...)
  #SUB <- CI[,3]-CI[,1] > .Machine$double.eps # still does not avoid annoying warning
  #suppressWarnings( graphics::arrows(t[SUB],CI[SUB,1],t[SUB],CI[SUB,3],length=0.05,angle=90,code=3,...) )
}
