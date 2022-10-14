ts.plot <- function(data,...)
{
  t <- data$t
  CI <- data[,NAMES.CI]

  # base plot
  plot(t,CI[,2],pch=19,...)

  SUB <- CI[,3]-CI[,1] > .Machine$double.eps # still does not avoid annoying warning
  suppressWarnings( graphics::arrows(t[SUB],CI[SUB,1],t[SUB],CI[SUB,3],length=0.05,angle=90,code=3,...) )
}
