# animate telemetry objects
# time label arguments (title?)
# time label units arguments
video <- function(x,ext=extent(x),fps=60,dt=NULL,ghost=0,timestamp=FALSE,file="ctmm.mp4",res=720,col="red",pch=1,cex=NULL,lwd=1,par.list=list(),...)
{
  x <- listify(x)
  n <- length(x)

  col <- format_par(col,x)
  pch <- format_par(pch,x)
  lwd <- format_par(lwd,x)
  #type <- format_par(type,x,all=TRUE)

  # automagic the plot point size
  if(is.null(cex))
  {
    p <- sum(sapply(x, function(d) { length(d$t) } ))
    if(p>1000) { cex <- 1000/p } else { cex <- 1 }
  }
  cex <- format_par(cex,x)

  res <- array(res,2)

  # time step for frames
  if(is.null(dt))
  {
    dt <- sapply(x,function(y){stats::median(diff(y$t))})
    dt <- stats::median(dt)
  }

  TIMES <- lapply(x,function(y){y$t})
  TIMES <- unlist(TIMES)
  t0 <- grid.init(TIMES,dt=dt) # best fraction of time to subtract for grid alignment
  rm(TIMES)

  t1 <- ext$t[1] - t0 # best initial time for grid
  while(ext$t[1] < t1-dt/2) { t1 <- t1 - dt }
  while(ext$t[1] > t1+dt/2) { t1 <- t1 + dt }

  t2 <- t1 + ceiling((ext$t[2]-t1)/dt)*dt
  while(ext$t[2] > t2+dt/2) { t2 <- t2 + dt }
  while(ext$t[2] < t2-dt/2) { t2 <- t2 - dt }

  nmax <- ceiling((t2-t1)/dt) + 1

  TIMES <- seq(t1,t2,length.out=nmax)

  # exp(-ghost.max/ghost) == 1/255
  ghost.max <- -ghost*log(1/255)

  animation::saveVideo({
    if(length(par.list)) { do.call(graphics::par,par.list) }

    I1 <- I2 <- rep(0,n)
    pb <- utils::txtProgressBar(style=3)
    for(j in 1:nmax)
    {
      t <- TIMES[j]
      TIMESTAMP <- rep(NA,n)
      y <- list() # to copy over
      COL <- ALPHA <- PCH <- LWD <- CEX <- list()
      # update current indices
      for(i in 1:n)
      {
        ALPHA[[i]] <- 1

        # increment forward
        if(I1[i]==0) { if(x[[i]]$t[1] <= t+dt/2) { I1[i] <- 1 } }

        if(I2[i]==0) { if(x[[i]]$t[1] <= t+dt/2) { I2[i] <- 1 } }

        if(I1[i]>0 && I1[i]<=nrow(x[[i]]))
        { while(I1[i]<=nrow(x[[i]]) && x[[i]]$t[I1[i]] < t-dt/2-ghost.max) { I1[i] <- I1[i]+1 } }

        if(I2[i]>0 && I2[i]<nrow(x[[i]]))
        { while(I2[i]<nrow(x[[i]]) && x[[i]]$t[I2[i]+1] < t+dt/2) { I2[i] <- I2[i]+1 } }

        # subset
        if(I1[i]>0 && I2[i]>0 && I1[i]<=I2[i] && I1[i]<=nrow(x[[i]]) && I2[i]<=nrow(x[[i]]))
        {
          I <- I1[i]:I2[i]
          if(I2[i]>I1[i]) # multiple times in window
          {
            # take closest time to window center
            I <- I[which.min(abs(t-x[[i]]$t[I]))]
            # interpolation should be done in predict/simulate --- not here
            if(ghost)
            {
              I <- I1[i]:I # include past locations
              ALPHA[[i]] <- exp(-(x[[i]]$t[last(I)]-x[[i]]$t[I])/dt)
            }
          }
          if(timestamp) { TIMESTAMP[i] <- x[[i]]$timestamp[last(I)] }
        }
        else
        { I <- NULL }
        y[[i]] <- x[[i]][I,]

        # possible time-dependent formatting
        COL[[i]] <- pull(col[[i]],I)
        if(length(I)) { COL[[i]] <- malpha(COL[[i]],ALPHA[[i]]) }
        PCH[[i]] <- pull(pch[[i]],I)
        LWD[[i]] <- pull(lwd[[i]],I)
        CEX[[i]] <- pull(cex[[i]],I)
        # TYPE[[i]] <- pull(type[[i]],I)
      } # end y populate

      plot(y,ext=ext,col=COL,pch=PCH,lwd=LWD,cex=CEX,...)
      if(timestamp)
      {
        TIMESTAMP <- stats::median(TIMESTAMP,na.rm=TRUE)
        graphics::title(TIMESTAMP)
      }
      utils::setTxtProgressBar(pb,j/nmax)
    } # end time loop
    close(pb)
  },video.name=file,interval=1/fps,nmax=nmax,ani.width=res[1],ani.height=res[2])
  # ,ani.res=max(SIZE)
}
