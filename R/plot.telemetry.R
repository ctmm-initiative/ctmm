################################
# ZOOM INTO TELEMETRY DATA
# this needs some work
zoom.telemetry <- function(x,fraction=1,...)
{
  manipulate::manipulate(
    { plot(x,fraction=fraction,...) },
    fraction = manipulate::slider(0, 2.0,initial=fraction,step=1/100)
  )
}
methods::setMethod("zoom",signature(x="telemetry"), function(x,fraction=1,...) zoom.telemetry(x,fraction=fraction,...))
methods::setMethod("zoom",signature(x="UD"), function(x,fraction=1,...) zoom.telemetry(x,fraction=fraction,...))


##############
new.plot <- function(data=NULL,CTMM=NULL,UD=NULL,level.UD=0.95,level=0.95,fraction=1,add=FALSE,xlim=NULL,ylim=NULL,...)
{
  dist <- list()
  dist$name <- "meters"
  dist$scale <- 1

  if(!add)
  {
    ext <- NULL

    # bounding locations from data
    if(!is.null(data))
    { ext <- rbind(ext,extent(data)) }

    # bounding locations from UDs
    if(!is.null(UD))
    { ext <- rbind(ext,extent(UD,level=level,level.UD=level.UD)) }

    # bounding locations from Gaussian CTMM
    if(!is.null(CTMM) & !any(is.na(level.UD)))
    { ext <- rbind(ext,extent(CTMM,level=level,level.UD=level.UD)) }

    # combine ranges
    ext <- extent.telemetry(ext)

    # bounding box
    mu <- c(mean(ext$x),mean(ext$y))
    buff <- c(diff(ext$x),diff(ext$y))/2

    # now zoom in/out to some fraction of the grid
    buff <- fraction*buff

    ext$x <- mu[1] + buff[1]*c(-1,1)
    ext$y <- mu[2] + buff[2]*c(-1,1)

    # try to obey xlim/ylim if provided
    if(!is.null(xlim) || !is.null(ylim))
    {
      max.diff <- max(diff(xlim),diff(ylim))*c(-1,1)/2

      if(is.null(ylim))
      { ylim <- mu[2] + max.diff }
      else if(is.null(xlim))
      { xlim <- mu[1] + max.diff }

      ext$x <- xlim
      ext$y <- ylim
    }

    # Get best unit scale
    dist <- unit(unlist(ext),"length")

    xlab <- paste("x ", "(", dist$name, ")", sep="")
    ylab <- paste("y ", "(", dist$name, ")", sep="")

    # residuals have no units
    if(!is.null(data) && !is.null(attr(data[[1]],"info")$residual))
    {
      xlab <- "x"
      ylab <- "y"
    }

    ext <- ext/dist$scale

    # empty base layer plot
    plot(ext, xlab=xlab, ylab=ylab, col=grDevices::rgb(1,1,1,0), asp=1, ...)
    assign("x.dim","length",pos=plot.env)
    assign("y.dim","length",pos=plot.env)
    assign("x.scale",dist$scale,pos=plot.env)
    assign("y.scale",dist$scale,pos=plot.env)
  }

  return(dist)
}
# setup environment
plot.env <- new.env()

#######################################
# PLOT TELEMETRY DATA
#######################################
plot.telemetry <- function(x,CTMM=NULL,UD=NULL,level.UD=0.95,level=0.95,DF="CDF",error=TRUE,velocity=FALSE,col="red",col.level="black",col.DF="blue",col.grid="white",transparency.error=0.25,pch=1,type='p',labels=NULL,fraction=1,add=FALSE,xlim=NULL,ylim=NULL,cex=NULL,lwd=1,lwd.level=1,...)
{
  alpha <- 1-level
  alpha.UD <- 1-level.UD

  # listify everything for generality
  x <- listify(x)
  CTMM <- listify(CTMM)
  UD <- listify(UD)

  # median time step of data
  dt <- lapply(x,function(X){diff(X$t)})
  dt <- stats::median(unlist(dt))

  dist <- new.plot(data=x,CTMM=CTMM,UD=UD,level.UD=level.UD,level=level,fraction=fraction,add=add,xlim=xlim,ylim=ylim,...)

  # plot cm per unit of distance plotted (m or km)
  cmpkm <- 2.54*mean(graphics::par("fin")*diff(graphics::par("plt"))[-2]/diff(graphics::par("usr"))[-2])
  # plot px per unit of distance plotted (m or km)
  pxpkm <- mean(grDevices::dev.size("px")*diff(graphics::par("plt"))[-2]/diff(graphics::par("usr"))[-2])

  #########################
  # PLOT GAUSSIAN CONTOURS AND DENSITY
  if(!is.null(CTMM))
  {
    # contours colour
    col.level <- array(col.level,length(CTMM))
    col.DF <- array(col.DF,length(CTMM))

    for(i in 1:length(CTMM))
    {
      # scale units
      CTMM[[i]] <- unit.ctmm(CTMM[[i]],dist$scale)

      # plot denisty function lazily reusing KDE code
      pdf <- kde(data.frame(CTMM[[i]]$mu[1,,drop=FALSE]),H=methods::getDataPart(CTMM[[i]]$sigma),axes=c("x","y"),res=500)
      plot.df(pdf,DF=DF,col=col.DF[[i]],...)

      # plot ML estimate, regular style
      plot.ctmm(CTMM[[i]],alpha.UD,col=col.level[[i]],lwd=lwd.level,...)

      # plot CIs dashed if present
      if("COV" %in% names(CTMM[[i]]) && !is.na(level))
      {
        # proportionality constants for outer CIs
        const <- confint.ctmm(CTMM[[i]],alpha)["area",]
        const <- const[c(1,3)]/const[2]
        sigma <- CTMM[[i]]$sigma

        for(j in 1:2)
        {
          CTMM[[i]]$sigma <- const[j]*sigma
          plot.ctmm(CTMM[[i]],alpha.UD,col=scales::alpha(col.level[[i]],0.5),lwd=lwd.level/2,...)
        }
      }
    }
  }

  ##########################################
  # PLOT KDE CONTOURS... AND DENSITY
  if(!is.null(UD))
  {
    UD <- lapply(UD,function(ud){ unit.UD(ud,length=dist$scale) })
    plot.UD(UD,level.UD=level.UD,level=level,DF=DF,col.level=col.level,col.DF=col.DF,col.grid=col.grid,labels=labels,fraction=fraction,add=TRUE,xlim=xlim,ylim=ylim,cex=cex,lwd=lwd.level,...)
  }

  #########################
  # PLOT TELEMETRY DATA

  # prepare point characteristics
  prepare.p <- function(pchar,all=FALSE)
  {
    if(!is.list(pchar))
    {
      if(length(x)>1)
      { pchar <- array(pchar,length(x)) }
      else if(!all)
      { pchar <- list(array(pchar,length(x[[1]]$t))) }
    }
    return(pchar)
  }

  col <- prepare.p(col)
  pch <- prepare.p(pch)
  type <- prepare.p(type,all=TRUE)
  error <- rep(error,length(x))

  # automagic the plot point size
  if(is.null(cex))
  {
    p <- sum(sapply(x, function(d) { length(d$t) } ))
    if(p>1000) { cex <- 1000/p } else { cex <- 1 }
  }
  cex <- prepare.p(cex)

  # standard deviations to plot for circle/ellipse
  z <- sqrt(-2*log(alpha.UD))

  # now plot individually
  for(i in 1:length(x))
  {
    r <- x[[i]][,c('x','y')]/dist$scale

    # scaled error info
    ERROR <- get.error(x[[i]],list(error=TRUE,axes=c('x','y'))) / dist$scale^2
    # don't want to throw z in here yet, in case of kernels
    FLAG <- attr(ERROR,"flag") # nothing, circle or ellipse?
    # we aren't plotting if UERE is missing

    # error=FALSE or no UERE
    if(!error[i] || FLAG<=1) # DEFAULT POINTS
    { graphics::points(r, cex=cex[[i]], col=col[[i]], pch=pch[[i]], type=type[[i]], lwd=lwd, ...) }
    else if(error[i]<3) # CIRCLE/ELLIPSE
    {
      # scale radii
      # ERROR <- z^2 * ERROR

      if("COV.major" %in% names(x[[i]]))
      {
        x[[i]][["COV.major"]] <- x[[i]][["COV.major"]]*(z/dist$scale)^2
        x[[i]][["COV.minor"]] <- x[[i]][["COV.minor"]]*(z/dist$scale)^2
      }

      # set coloring to outer rim or solid disc
      if(error[i]==1) # RIM
      {
        # elliptical circumference
        if(all(DOP.LIST$horizontal$COV %in% names(x[[i]])))
        {
          if("COV.major" %in% names(x[[i]])) # eigen-values already present
          {
            A2 <- x[[i]][["COV.major"]]
            B2 <- x[[i]][["COV.minor"]]
          }
          else
          {
            B2 <- vapply(1:(dim(ERROR)[1]),function(i){ eigen(ERROR[i,,],only.values=TRUE)$values },numeric(2)) # (big/small,n)
            A2 <- B2[1,]
            B2 <- B2[2,]
          }

          B2 <- ifelse(B2<A2,B2/A2,1) # prevent 0/0
          alpha <- (4*z) * sqrt(A2) * pracma::ellipke(sqrt(1-B2))$e
          rm(A2,B2)
        }
        else # circular circumference
        { alpha <- (2*pi*z)*sqrt(ERROR) }

        # circumference of a pixel in physical units, spread around circumference
        alpha <- clamp((4/pxpkm) / alpha)^transparency.error

        bg <- NA
        fg <- scales::alpha(col[[i]],alpha)
      }
      else if(error[i]==2) # DISC
      {
        if((all(DOP.LIST$horizontal$COV %in% names(x[[i]]))))
        {
          if("COV.major" %in% names(x[[i]]))
          { alpha <- (pi*z^2) * sqrt(x[[i]][["COV.major"]] * x[[i]][["COV.minor"]]) }
          else
          { alpha <- (pi*z^2) * sqrt(apply(ERROR,1,det)) }
        }
        else
        { alpha <- (pi*z^2) * ERROR }

        # area of a pixel in physical units, spread over disc
        alpha <- clamp((1/pxpkm^2) / alpha)^transparency.error

        bg <- scales::alpha(col[[i]],alpha)
        fg <- scales::alpha(col[[i]],alpha/2)
      }

      # plot circle
      if(FLAG<4)
      {
        # convert to radii
        ERROR <- z * sqrt(ERROR)
        graphics::symbols(x=r$x,y=r$y,circles=ERROR,fg=fg,bg=bg,inches=FALSE,add=TRUE,lwd=lwd,...)
      }
      else if(FLAG==4)
      { for(j in 1:nrow(r)) { ellipsograph(mu=as.numeric(r[j,]),sigma=ERROR[j,,],level=level.UD,fg=fg[j],bg=bg[j],lwd=lwd,...) } }
    } # end circle/ellipse plot
    else if(error[i]==3) # kernels
    {
      # setup grid with correct extent
      EXT <- array( graphics::par("usr") , c(2,2) )
      GRID <- kde.grid(r,ERROR,res=max(grDevices::dev.size("px")),EXT=EXT)
      # calculate kernels
      UD <- kde(r,ERROR,grid=GRID)
      UD <- new.UD(UD,info=list())
      # plot kernels
      plot.UD(UD,level.UD=NA,level=NA,DF='PDF',col.DF=col[[i]],col.level=NA,col.grid=NA,add=TRUE,...)
    } # end kernel plot

    # also plot velocity vectors at dt scale
    if(velocity && all(c("vx","vy") %in% names(x[[i]])))
    {
      dr <- x[[i]][,c("vx","vy")]/dist$scale*dt

      arr.length <- sqrt(rowSums(dr^2))
      arr.length <- 0.1*cmpkm*arr.length

      # make uncertain speeds more transparent
      if(DOP.LIST$speed$VAR %in% names(x[[i]]))
      {
        # magnitude of estimate relative to uncertainty
        ERROR <- get.error(x[[i]],ctmm(axes=c("vx","vy")),circle=TRUE)
        alpha <- sqrt(rowSums(get.telemetry(x[[i]],axes=c("vx","vy"))^2)/ERROR)
        alpha <- clamp(alpha)
        col[[i]] <- scales::alpha(col[[i]],alpha)
      }

      shape::Arrows(x0=r$x, y0=r$y, x1=(r$x+dr$vx), y1=(r$y+dr$vy), col=col[[i]], code=2, segment=T, arr.adj=1, arr.length=arr.length, arr.type="curved", lwd=lwd)
    }
  } # end telemetry loop
}
# SET METHODS FOR PLOT.TELEMETRY
#methods::setMethod("plot",signature(x="telemetry",y="missing"), function(x,y,...) plot.telemetry(x,...))
#methods::setMethod("plot",signature(x="telemetry",y="telemetry"), function(x,y,...) plot.telemetry(list(x,y),...))
#methods::setMethod("plot",signature(x="telemetry",y="ctmm"), function(x,y,...) plot.telemetry(x,model=y,...))
#methods::setMethod("plot",signature(x="telemetry",y="UD"), function(x,y,...) plot.telemetry(x,akde=y,...))
#methods::setMethod("plot",signature(x="telemetry"), function(x,...) plot.telemetry(x,...))


##############
plot.UD <- function(x,level.UD=0.95,level=0.95,DF="CDF",col.level="black",col.DF="blue",col.grid="white",labels=NULL,fraction=1,add=FALSE,xlim=NULL,ylim=NULL,cex=NULL,lwd=1,...)
{
  x <- listify(x)

  dist <- new.plot(UD=x,fraction=fraction,add=add,xlim=xlim,ylim=ylim,level.UD=level.UD,level=level,...)

  # contours colour
  if(length(col.level)==length(level.UD) && length(col.level) != length(x))
  { col.level <- t(array(col.level,c(length(level.UD),length(x)))) }
  else
  { col.level <- array(col.level,c(length(x),length(level.UD))) }
  col.level <- array(col.level,c(length(x),length(level.UD),3))

  # contour labels
  if(is.null(labels)) { labels <- round(100*level.UD) }
  if((length(labels)==length(level.UD) || length(labels)==3*length(level.UD)) && length(labels) != length(x))
  {
    labels <- array(labels,c(length(level.UD),3,length(x)))
    labels <- aperm(labels,c(3,1,2))
  }
  else
  { labels <- array(labels,c(length(x),length(level.UD),3)) }

  col.DF <- array(col.DF,length(x))
  col.grid <- array(col.grid,length(x))

  # UNIT CONVERSIONS
  for(i in 1:length(x))
  {
    # unit conversion
    x[[i]] <- unit.UD(x[[i]],length=dist$scale)

    # ML DENSITY PLOTS
    plot.df(x[[i]],DF=DF,col=col.DF[[i]],...)
  }

  # plot grid
  for(i in 1:length(x))
  {
    if(sum(diag(x[[i]]$H)>0))
    {
      H <- covm(x[[i]]$H)
      theta <- H@par["angle"]
      ecc <- H@par["eccentricity"]
      sigma <- H@par["area"]

      X <- x[[i]]$r$x
      Y <- x[[i]]$r$y

      COS <- cos(theta)
      SIN <- sin(theta)

      # grid spacing
      du <- sqrt(sigma*exp(+ecc/2))
      dv <- sqrt(sigma*exp(-ecc/2))

      # bandwidth axes
      u <- outer(+X*COS,+Y*SIN,"+")
      v <- outer(-X*SIN,+Y*COS,"+")

      # extent of data
      B <- (x[[i]]$PDF > 0)
      u <- u[B]
      v <- v[B]

      ex.u <- range(u)
      ex.v <- range(v)

      mu.u <- mean(ex.u)
      mu.v <- mean(ex.v)

      # grid numbers
      n.u <- diff(ex.u)/du
      n.v <- diff(ex.v)/dv

      n.u <- ceiling(n.u/2)
      n.v <- ceiling(n.v/2)

      # grid nodes
      u <- mu.u + du*(-n.u):n.u
      v <- mu.v + dv*(-n.v):n.v

      # transform back
      X <- outer(u*COS,-v*SIN,"+")
      Y <- outer(u*SIN,+v*COS,"+")

      for(j in 1:length(u)) { graphics::segments(x0=X[j,1],y0=Y[j,1],x1=last(X[j,]),y1=last(Y[j,]),col=col.grid[i],...) }
      for(j in 1:length(v)) { graphics::segments(x0=X[1,j],y0=Y[1,j],x1=last(X[,j]),y1=last(Y[,j]),col=col.grid[i],...) }
    }

    # not sure why this is necessary
    graphics::box(lwd=lwd,...)
  }

  # CONTOURS
  for(i in 1:length(x))
  {
    if(!any(is.na(col.level[i,,])) && !any(is.na(level.UD)))
    {
      # make sure that correct style is used for low,ML,high even in absence of lows and highs
      plot.kde(x[[i]],level=level.UD,labels=labels[i,,2],col=scales::alpha(col.level[i,,2],1),lwd=lwd,...)

      if(!is.na(level) && !is.null(x[[i]]$DOF.area))
      {
        P <- sapply(level.UD, function(l) { CI.UD(x[[i]],l,level,P=TRUE)[-2] } )
        plot.kde(x[[i]],level=P,labels=c(t(labels[i,,c(1,3)])),col=scales::alpha(c(t(col.level[i,,c(1,3)])),0.5),lwd=lwd/2,...)
      }
    }
  }

}

##################################
# plot PDF stored as KDE object
plot.df <- function(kde,DF="CDF",col="blue",...)
{
  col <- scales::alpha(col,(0:255)/255)

  if(DF=="PDF")
  {
    zlim <- c(0,max(kde$PDF))
  }
  else if(DF=="CDF")
  {
    zlim <- c(0,1)
    kde$CDF <- 1 - kde$CDF
  }

  graphics::image(kde$r,z=kde[[DF]],useRaster=TRUE,zlim=zlim,col=col,add=TRUE,...)
}


#############################
# Plot a KDE object's contours
plot.kde <- function(kde,level=0.95,labels=round(level*100),col="black",...)
{
  # record current option
  # MAX <- getOption("max.contour.segments")

  # do something that works
  options(max.contour.segments=.Machine$integer.max)
  drawlabels <- !(labels==FALSE | is.na(labels))
  graphics::contour(kde$r,z=kde$CDF,levels=level,labels=labels,labelcex=1,drawlabels=drawlabels,col=col,add=TRUE,...)

  # reinstate initial option (or default if was NULL--can't set back to NULL???)
  # if(is.null(MAX)) { MAX <- 25000 }
  # options(max.contour.segments=MAX)
}


##############################
# Plot Gaussian ctmm contours
plot.ctmm <- function(model,alpha=0.05,col="blue",bg=NA,...)
{
  mu <- model$mu # mean vector
  sigma <- model$sigma # covariance matrix

  ellipsograph(mu=mu,sigma=sigma,level=1-alpha,fg=col,bg=bg,...)
}

###################
# mu - mean vector
# sigma - covariance matrix
ellipsograph <- function(mu,sigma,level=0.95,fg=graphics::par("col"),bg=NA,...)
{
  Eigen <- eigen(sigma)
  std <- sqrt(Eigen$values)
  vec <- Eigen$vectors

  # confidence level = 1-alpha
  alpha <- 1 - level
  z <- sqrt(-2*log(alpha))

  num <- 100 # number of plotted points
  theta <- 2*pi*(0:num)/(num+1)
  Sin <- sin(theta)
  Cos <- cos(theta)

  x <- mu[1] + z*(Cos*std[1]*vec[1,1] + Sin*std[2]*vec[1,2])
  y <- mu[2] + z*(Cos*std[1]*vec[2,1] + Sin*std[2]*vec[2,2])

  graphics::xspline(x, y=y, shape=-1, open=FALSE, border=fg, col=bg, ...)
}
