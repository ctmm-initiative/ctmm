compass <- function(loc=NULL,cex=3,...)
{
  proj <- get0('projection',plot.env)
  scale <- unique( c( get0('x.scale',plot.env), get0('y.scale',plot.env) ) )

  if(is.null(loc))
  {
    loc <- graphics::par('usr')
    dx <- min( diff(loc[1:2]), diff(loc[3:4]) )
    loc <- loc[c(2,4)] - 0.15*dx
    loc <- loc * scale
  }
  loc <- data.frame(x=loc[1],y=loc[2])

  geo <- project(loc,from=proj,to=DATUM)
  colnames(geo) <- c("longitude","latitude")
  loc <- cbind(loc,geo)

  # angle to rotate to north
  srt <- northing(loc,proj,angle=TRUE)

  x <- loc$x/scale
  y <- loc$y/scale

  # graphics::text(x,y,labels="\u2742",adj=c(0.5,0.5),srt=srt,cex=1.5*cex,...) # too thick
  graphics::text(x,y,labels="\u27A2",adj=c(0.5,0.5),srt=srt,cex=cex,...)
  # graphics::text(x,y,labels="N",adj=c(0.5,0.5),srt=srt-90,cex=cex,...) # offset
}
