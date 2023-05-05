transition <- function(data,n=3,filename="transition",height=2160,...)
{
  # rotate to vertical alignment
  projection(data) <- median(data,k=2)
  TEMP <- data
  data$x <- -TEMP$y
  data$y <- data$x
  rm(TEMP)

  EXT <- extent(data)[,c('x','y')]
  ASP <- apply(EXT,2,diff)
  width <- ASP[1]/ASP[2] * height

  COL <- color(data,by='time')

  i <- 0
  FILE <- paste0(filename,"-",i,".png")
  grDevices::png(FILE,width=width,height=height)
  graphics::par(mar=c(0,0,0,0)) # restore after !!!!!!!!!!!!!!
  plot(data,col=COL,bty="n",axes=F,frame.plot=F,xaxt='n',ann=FALSE,yaxt='n',...)
  grDevices::dev.off()

  P <- diff(range(data$t))

  for(i in 1:n)
  {
    col <- grDevices::col2rgb(COL,alpha=TRUE)
    alpha <- col['alpha',]/255
    col <- col[1:3,]
    col <- grDevices::rgb2hsv(col)

    SUB <- data$t[1] + P*(i-1)/n <= data$t & data$t <= data$t[1] + P*i/n
    col['s',!SUB] <- 0
    col <- grDevices::hsv(h=col['h',],s=col['s',],v=col['v',],alpha=alpha)

    FILE <- paste0("../presentations LaTeX/coati",i,".png")
    grDevices::png(FILE,width=width,height=height)
    graphics::par(mar=c(0,0,0,0))
    plot(data,col=col,bty="n",axes=F,frame.plot=F,xaxt='n',ann=FALSE,yaxt='n',...)
    grDevices::dev.off()
  }
}
