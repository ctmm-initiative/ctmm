ts.plot <- function(data,col='black',ylab="",ylim=NULL,...)
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
    xlab <- "Time"
  }

  CI <- data[,NAMES.CI]

  MAX <- max(CI[,'est'])
  MIN <- min(CI[,'est'])
  DIFF <- MAX-MIN
  if(is.null(ylim)) { ylim <- c( max(MIN-DIFF,min(CI[,'low'])) , min(MAX+DIFF,max(CI[,'high'])) ) }

  graphics::plot(range(t),range(CI[,c('low','high')]),col=rgb(1,1,1,0),ylim=ylim,xlab=xlab,ylab=ylab,...)
  graphics::polygon(c(t,rev(t),t[1]),c(CI[,'low'],CI[,'high'],CI[1,'low']),col=malpha(col,alpha=0.25),border=NA,...)
  graphics::points(t,CI[,'est'],type="l",col=col,...)

  # ERROR BAR PLOT
  # plot(t,CI[,2],pch=19,...)
  #SUB <- CI[,3]-CI[,1] > .Machine$double.eps # still does not avoid annoying warning
  #suppressWarnings( graphics::arrows(t[SUB],CI[SUB,1],t[SUB],CI[SUB,3],length=0.05,angle=90,code=3,...) )
}
