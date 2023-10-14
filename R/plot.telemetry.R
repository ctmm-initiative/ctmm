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
new.plot <- function(data=NULL,CTMM=NULL,UD=NULL,R=NULL,col.bg="white",col.R="green",legend=FALSE,level.UD=0.95,level=0.95,units=TRUE,fraction=1,add=FALSE,xlim=NULL,ylim=NULL,ext=NULL,cex=NULL,...)
{
  RESIDUALS <- !is.null(data) && !is.null(attr(data[[1]],"info")$residual)

  if(is.null(cex)) { cex <- graphics::par('cex') }
  if(class(cex)[1]=='list') { cex <- sapply(cex,stats::median) }
  cex <- stats::median(cex)

  dist <- list()
  dist$name <- "meters"
  dist$scale <- 1
  axes <- c("x","y")

  if(!is.null(ext))
  {
    xlim <- ext$x
    ylim <- ext$y
  }

  if(!add)
  {
    ext <- NULL

    # bounding locations from data
    if(!is.null(data))
    { ext <- rbind(ext,extent(data)[,axes]) }

    # bounding locations from UDs
    if(!is.null(UD))
    {
      if(class(UD[[1]])[1]=='RS') # backup extent for RS objects
      {
        ext <- data.frame(x=1:2,y=1:2)
        rownames(ext) <- c('min','max')
        ext['min','x'] <- min(sapply(UD,function(U){U$r$x[1]}))
        ext['max','x'] <- min(sapply(UD,function(U){last(U$r$x)}))
        ext['min','y'] <- min(sapply(UD,function(U){U$r$y[1]}))
        ext['max','y'] <- min(sapply(UD,function(U){last(U$r$y)}))
      }
      else
      { ext <- rbind(ext,extent(UD,level=level,level.UD=level.UD)[,axes]) }
    }

    # bounding locations from Gaussian CTMM
    if(!is.null(CTMM) & !any(is.na(level.UD)))
    { ext <- rbind(ext,extent(CTMM,level=level,level.UD=level.UD)[,axes]) }

    # # bounding locations from standard normal quantiles
    # if(RESIDUALS && !any(is.na(level.UD)))
    # {
    #   Z <- c(1,1) * sqrt(-2*log(1-max(level.UD)))
    #   ext <- rbind(ext,-Z)
    #   ext <- rbind(ext,+Z)
    # }

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
    dist <- unit(unlist(ext),"length",SI=!units)

    xlab <- paste("x ", "(", dist$name, ")", sep="")
    ylab <- paste("y ", "(", dist$name, ")", sep="")

    # residuals have no units
    if(RESIDUALS)
    {
      xlab <- "x"
      ylab <- "y"
    }

    ext <- ext/dist$scale

    # empty base layer plot
    plot(ext, xlab=xlab, ylab=ylab, col=grDevices::rgb(1,1,1,0), asp=1, cex=cex, ...)

    # plot background color
    lim <- graphics::par('usr')
    xlim <- lim[1:2]
    dx <- diff(xlim)
    ylim <- lim[3:4]
    dy <- diff(ylim)
    # dx <- dy <- 0
    graphics::rect(xlim[1]-dx/2,ylim[1]-dy/2,xlim[2]+dx/2,ylim[2]+dy/2,border=col.bg,col=col.bg)
    # this can cover the plot box
    graphics::box()

    # plot information for further layering
    projection <- unique(c(projection(data),projection(CTMM),projection(UD))) # some objects could be NULL
    if(length(projection)>1 && !RESIDUALS) { stop("Multiple projections not yet supported.") }
    assign("projection",projection,pos=plot.env)
    # dimensional type
    assign("x.dim","length",pos=plot.env)
    assign("y.dim","length",pos=plot.env)
    # unit name
    assign("x.units",dist$name,pos=plot.env)
    assign("y.units",dist$name,pos=plot.env)
    # unit conversion
    assign("x.scale",dist$scale,pos=plot.env)
    assign("y.scale",dist$scale,pos=plot.env)

    scale <- dist$scale
  } # end !add
  else # get distance information from environment
  {
    name <- unique( c( get0('x.units',plot.env), get0('y.units',plot.env) ) )
    scale <- unique( c( get0('x.scale',plot.env), get0('y.scale',plot.env) ) )
    if(length(name)==1 && length(scale)==1)
    {
      dist$name <- name
      dist$scale <- scale
    }
    projection <- get('projection',plot.env)
  }

  # PLOT RASTER / SUITABILITY
  if(!is.null(R))
  { plot_R(R,col=col.R,legend=legend,projection=projection,scale=scale) }

  return(dist)
}
# setup environment
plot.env <- new.env()

#######################################
# PLOT TELEMETRY DATA
#######################################
plot.telemetry <- function(x,CTMM=NULL,UD=NULL,col.bg='white',
                           cex=NULL,col="red",lwd=1,pch=1,type='p',error=TRUE,transparency.error=0.25,velocity=FALSE,
                           DF="CDF",col.UD="blue",col.grid="white",labels=NULL,convex=FALSE,level=0.95,level.UD=0.95,col.level="black",lwd.level=1,
                           SP=NULL,border.SP=TRUE,col.SP=NA,
                           R=NULL,col.R="green",legend=FALSE,
                           fraction=1,xlim=NULL,ylim=NULL,ext=NULL,units=TRUE,add=FALSE,...)
{
  alpha <- 1-level
  alpha.UD <- 1-level.UD

  # list-ify everything for generality
  x <- listify(x)
  CTMM <- listify(CTMM)
  UD <- listify(UD)
  R <- listify(R)
  # fix argument order
  if(class(CTMM[[1]])[1]=="UD") { UD <- CTMM; CTMM <- NULL }
  if(class(CTMM[[1]])[1]=="RasterLayer") { R <- CTMM; CTMM <- NULL }

  # catch 3D UDs
  if(length(dim(UD[[1]]$CDF))==3)
  { return(plot3d(data=x,UD=UD,level=level,level.UD=level.UD,xlim=xlim,ylim=ylim,ext=ext,
                  cex=cex,col=col,lwd=lwd,pch=pch,type=type,error=error,transparency.error=transparency.error,velocity=velocity,
                  DF=DF,col.UD=col.UD,col.grid=col.grid,labels=labels,col.level=col.level,lwd.level=lwd.level,
                  SP=SP,border.SP=border.SP,col.SP=col.SP,
                  fraction=fraction,units=units,add=add,...)) }

  # median time step of data
  dt <- lapply(x,function(X){diff(X$t)})
  dt <- stats::median(unlist(dt))

  # are we plotting residuals?
  RESIDUALS <- !is.null(x) && !is.null(attr(x[[1]],"info")$residual)
  # standard normal model
  if(RESIDUALS)
  {
    error <- FALSE
    CTMM <- list(ctmm(sigma=1,mu=c(0,0)))
  }

  dist <- new.plot(data=x,CTMM=CTMM,UD=UD,col.bg=col.bg,R=R,col.R=col.R,legend=legend,level.UD=level.UD,level=level,units=units,fraction=fraction,add=add,xlim=xlim,ylim=ylim,ext=ext,cex=cex,...)

  # plot cm per unit of distance plotted (m or km)
  cmpkm <- 2.54*mean(graphics::par("fin")*diff(graphics::par("plt"))[-2]/diff(graphics::par("usr"))[-2])
  # plot px per unit of distance plotted (m or km)
  pxpkm <- mean(grDevices::dev.size("px")*diff(graphics::par("plt"))[-2]/diff(graphics::par("usr"))[-2])

  #########################3
  # PLOT SHAPEFILES
  if(!is.null(SP)) { plot_SP(SP=SP,border.SP=border.SP,col.SP=col.SP,PROJ=ctmm::projection(x),...) }

  # unlist annotation colors
  col.level <- simplify.color(col.level)
  col.UD <- simplify.color(col.UD)

  #########################
  # PLOT GAUSSIAN CONTOURS AND DENSITY
  if(!is.null(CTMM))
  {
    # contours colour
    col.level <- array(col.level,length(CTMM))
    col.UD <- array(col.UD,length(CTMM))

    for(i in 1:length(CTMM))
    {
      # scale units
      CTMM[[i]] <- unit.ctmm(CTMM[[i]],dist$scale)

      # plot denisty function lazily reusing KDE code
      pdf <- agde(CTMM[[i]],res=500)
      plot_df(pdf,DF=DF,col=col.UD[[i]],...)

      # plot ML estimate, regular style
      plot_ctmm(CTMM[[i]],alpha.UD,col=col.level[[i]],lwd=lwd.level,...)

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
          plot_ctmm(CTMM[[i]],alpha.UD,col=malpha(col.level[[i]],0.5),lwd=lwd.level/2,...)
        }
      }
    }
  }

  ##########################################
  # PLOT KDE CONTOURS... AND DENSITY
  if(!is.null(UD))
  {
    # UD <- lapply(UD,function(ud){ unit.UD(ud,length=dist$scale) }) # now done in plot.UD
    plot.UD(UD,level.UD=level.UD,level=level,DF=DF,col.level=col.level,col.UD=col.UD,col.grid=col.grid,labels=labels,fraction=fraction,add=TRUE,xlim=xlim,ylim=ylim,ext=ext,cex=cex,lwd.level=lwd.level,convex=convex,...)
  }

  #########################
  # PLOT TELEMETRY DATA

  col <- format_par(col,x)
  pch <- format_par(pch,x)
  lwd <- format_par(lwd,x)
  type <- format_par(type,x,all=TRUE)
  error <- rep(error,length(x))

  # automagic the plot point size
  if(is.null(cex) && length(x))
  {
    p <- sum(sapply(x, function(d) { length(d$t) } ))
    if(p>1000) { cex <- 1000/p } else { cex <- 1 }
  }
  cex <- format_par(cex,x)

  # standard deviations to plot for circle/ellipse
  z <- sqrt(-2*log(alpha.UD))

  # now plot individually
  for(i in 1%:%length(x))
  {
    if(!nrow(x[[i]])) { next } # skip empty data

    x[[i]] <- unit.telemetry(x[[i]],length=dist$scale)

    r <- x[[i]][,c('x','y')]

    # scaled error info
    if(error[i] && is.calibrated(x[[i]])<1) { uere(x[[i]]) <- uere(x[[i]]) } # force calibration for plotting

    if(error[i])
    {
      UERE <- uere(x[[i]])
      ERROR <- UERE$UERE[,'horizontal']
      names(ERROR) <- rownames(UERE$UERE) # R drops dimnames
      ERROR <- ctmm(error=ERROR,axes=c('x','y'))
      ERROR <- get.error(x[[i]],ERROR,calibrate=TRUE)
      # don't want to throw z in here yet, in case of kernels
      ELLIPSE <- length(dim(ERROR))>0 # circle or ellipse?
    }

    # we aren't plotting if UERE is missing
    # error=FALSE or no UERE
    if(!error[i] || all(attr(x[[i]],"UERE")$UERE[,'horizontal']==0)) # DEFAULT POINTS
    { graphics::points(r, cex=cex[[i]], col=col[[i]], pch=pch[[i]], type=type[[i]], lwd=lwd[[i]], ...) }
    else if(error[i]<3) # CIRCLE/ELLIPSE
    {
      # scale radii
      # ERROR <- z^2 * ERROR

      if("COV.major" %in% names(x[[i]]))
      {
        x[[i]][["COV.major"]] <- x[[i]][["COV.major"]]*(z)^2
        x[[i]][["COV.minor"]] <- x[[i]][["COV.minor"]]*(z)^2
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
            B2 <- vapply(1:(dim(ERROR)[1]),function(i){ eigen(ERROR[i,,])$values },numeric(2)) # (big/small,n)
            A2 <- clamp(B2[1,],0,Inf)
            B2 <- clamp(B2[2,],0,A2)
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
        fg <- malpha(col[[i]],alpha)
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

        bg <- malpha(col[[i]],alpha)
        fg <- malpha(col[[i]],alpha/2)
      }

      # plot circle
      if(!ELLIPSE)
      {
        # convert to radii
        ERROR <- z * sqrt(ERROR)
        graphics::symbols(x=r$x,y=r$y,circles=ERROR,fg=fg,bg=bg,inches=FALSE,add=TRUE,lwd=lwd[[i]],...)
      }
      else if(ELLIPSE)
      {
        if(length(lwd[[i]])<nrow(r)) { lwd[[i]] <- rep(lwd[[i]],nrow(r)) }
        for(j in 1:nrow(r)) { ellipsograph(mu=as.numeric(r[j,]),sigma=ERROR[j,,],level=level.UD,fg=fg[j],bg=bg[j],lwd=lwd[[i]][j],...) }
      }
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
      plot.UD(UD,level.UD=NA,level=NA,DF='PDF',col.UD=col[[i]],col.level=NA,col.grid=NA,add=TRUE,...)
    } # end kernel plot

    # also plot velocity vectors at dt scale
    if(velocity && all(c("vx","vy") %in% names(x[[i]])))
    {
      dr <- x[[i]][,c("vx","vy")]*dt

      arr.length <- sqrt(rowSums(dr^2))
      arr.length <- 0.1*cmpkm*arr.length

      # make uncertain speeds more transparent
      if(DOP.LIST$speed$VAR %in% names(x[[i]]))
      {
        # magnitude of estimate relative to uncertainty
        ERROR <- get.error(x[[i]],ctmm(axes=c("vx","vy"),error=TRUE),circle=TRUE)
        alpha <- sqrt(rowSums(get.telemetry(x[[i]],axes=c("vx","vy"))^2)/ERROR)
        alpha <- clamp(alpha)
        col[[i]] <- malpha(col[[i]],alpha)
      }

      shape::Arrows(x0=r$x, y0=r$y, x1=(r$x+dr$vx), y1=(r$y+dr$vy), col=col[[i]], code=2, segment=T, arr.adj=1, arr.length=arr.length, arr.type="curved", lwd=lwd[[i]])
    }
  } # end telemetry loop
}
# SET METHODS FOR PLOT.TELEMETRY
#methods::setMethod("plot",signature(x="telemetry",y="missing"), function(x,y,...) plot.telemetry(x,...))
#methods::setMethod("plot",signature(x="telemetry",y="telemetry"), function(x,y,...) plot.telemetry(list(x,y),...))
#methods::setMethod("plot",signature(x="telemetry",y="ctmm"), function(x,y,...) plot.telemetry(x,model=y,...))
#methods::setMethod("plot",signature(x="telemetry",y="UD"), function(x,y,...) plot.telemetry(x,akde=y,...))
#methods::setMethod("plot",signature(x="telemetry"), function(x,...) plot.telemetry(x,...))

plot.ctmm <- function(x,data=NULL,UD=NULL,col.bg='white',
                           cex=NULL,col="red",lwd=1,pch=1,type='p',error=TRUE,transparency.error=0.25,velocity=FALSE,
                           DF="CDF",col.UD="blue",col.grid="white",labels=NULL,convex=FALSE,level=0.95,level.UD=0.95,col.level="black",lwd.level=1,
                           SP=NULL,border.SP=TRUE,col.SP=NA,
                           R=NULL,col.R="green",legend=FALSE,
                           fraction=1,xlim=NULL,ylim=NULL,ext=NULL,units=TRUE,add=FALSE,...)
{ plot.telemetry(x=data,CTMM=x,UD=UD,col.bg=col.bg,cex=cex,col=col,lwd=lwd,pch=pch,type=type,error=error,
                 transparency.error=transparency.error,velocity=velocity,DF=DF,col.UD=col.UD,col.grid=col.grid,
                 labels=labels,convex=convex,level=level,level.UD=level.UD,col.level=col.level,lwd.level=lwd.level,
                 SP=SP,border.SP=border.SP,col.SP=col.SP,R=R,col.R=col.R,legend=legend,fraction=fraction,xlim=xlim,
                 ylim=ylim,ext=ext,units=units,add=add,...) }

# format point characteristics for dataset x
format_par <- function(pchar,x,all=FALSE)
{
  if(!is.list(pchar) && !is.null(pchar))
  {
    if(length(x)>1)
    {
      pchar <- array(pchar,length(x))
      pchar <- as.list(pchar)
    }
    else if(!all)
    { pchar <- list(array(pchar,length(x[[1]]$t))) }
  }
  return(pchar)
}


# get par that may be a constant or variable
pull <- function(pchar,i)
{
  if(length(pchar)>1) { pchar <- pchar[i] }
  return(pchar)
}


# plot raster layer
plot_R <- function(R,col="green",legend=FALSE,projection="",scale=1)
{
  R <- listify(R)

  for(i in 1:length(R))
  {
    PROJ <- raster::projection(R[[i]])

    if(PROJ!=projection) # reproject raster
    {
      message("Reprojecting raster.")

      ext <- graphics::par('usr') * scale
      ext <- raster::extent(ext)

      res <- raster::res(R[[i]])
      res <- min(res)

      if(grepl("+units=km",PROJ)) { res <- res * 1000 }
      if(grepl("long",PROJ) && grepl("lat",PROJ)) { res <- 2*pi*DATA.EARTH$R.EQ/360 * res }

      CRS <- sp::CRS(projection)

      # res <- res/10 # ensure that reprojection looks good

      TEMP <- raster(resolution=res,ext=ext,crs=CRS)

      R[[i]] <- raster::projectRaster(R[[i]],TEMP)
    }

    if(scale==1000) # km scale (not meters)
    { raster::extent(R[[i]]) <- raster::extent(R[[i]])[]/1000 }

    if(length(col)==length(R))
    { COL <- col[i] }
    else
    { COL <- col }
    COL <- malpha(COL,(0:255)/255)

    maxpixels <- raster::ncell(R[[i]])
    raster::plot(R[[i]],col=COL,legend=legend,maxpixels=maxpixels,add=TRUE)
  }
}


# plot shapefiles
plot_SP <- function(SP=NULL,border.SP=TRUE,col.SP=NA,PROJ=NULL,...)
{
  x.scale <- get0('x.scale',plot.env)
  if(x.scale==1000 && grepl("+units=m",PROJ)) # km scale (not meters)
  {
    PROJ <- strsplit(PROJ,"units=m")[[1]]
    PROJ <- paste0(PROJ[1],"units=km",PROJ[2])
  }

  # spTransform is bad
  if(!any(grepl('sf',class(SP)))) { SP <- sf::st_as_sf(SP) }
  SP <- sf::st_transform(SP,crs=sf::st_crs(PROJ))
  SP <- sf::as_Spatial(SP)

  sp::plot(SP,col=col.SP,border=border.SP,add=TRUE)
}


##############
plot.UD <- function(x,col.bg="white",DF="CDF",col.UD="blue",col.grid="white",labels=NULL,convex=FALSE,level=0.95,level.UD=0.95,col.level="black",lwd.level=1,
                    SP=NULL,border.SP=TRUE,col.SP=NA,
                    R=NULL,col.R="green",legend=FALSE,
                    fraction=1,xlim=NULL,ylim=NULL,ext=NULL,units=TRUE,add=FALSE,...)
{
  x <- listify(x)

  # catch 3D UDs
  if(length(dim(x[[1]]$CDF))==3)
  { return(plot3d(UD=x,level=level,level.UD=level.UD,xlim=xlim,ylim=ylim,ext=ext,
                  DF=DF,col.UD=col.UD,col.grid=col.grid,labels=labels,col.level=col.level,lwd.level=lwd.level,
                  SP=SP,border.SP=border.SP,col.SP=col.SP,
                  fraction=fraction,units=units,add=add,...)) }

  if(class(x[[1]])[1]=="RS") { DF <- 'RS' }

  dist <- new.plot(UD=x,R=R,col.bg=col.bg,col.R=col.R,legend=legend,units=units,fraction=fraction,xlim=xlim,ylim=ylim,ext=ext,level.UD=level.UD,level=level,add=add,...)

  # PLOT SHAPEFILES
  if(!is.null(SP)) { plot_SP(SP=SP,border.SP=border.SP,col.SP=col.SP,PROJ=ctmm::projection(x),...) }

  # contours colour
  if(length(col.level)==length(level.UD) && length(col.level) != length(x))
  { col.level <- t(array(col.level,c(length(level.UD),length(x)))) }
  else
  { col.level <- array(col.level,c(length(x),length(level.UD))) }
  col.level <- array(col.level,c(length(x),length(level.UD),3))

  # contour labels
  if(is.null(labels)) { labels <- round(100*level.UD) }

  LABEL.CI <- (length(labels)<3*length(level.UD))

  if((length(labels)==length(level.UD) || length(labels)==3*length(level.UD)) && length(labels) != length(x))
  {
    labels <- array(labels,c(length(level.UD),3,length(x)))
    labels <- aperm(labels,c(3,1,2))
  }
  else
  { labels <- array(labels,c(length(x),length(level.UD),3)) }

  if(LABEL.CI)
  {
    labels[,,1] <- paste(labels[,,1],"(low)")
    labels[,,3] <- paste(labels[,,3],"(high)")
  }

  col.UD <- array(col.UD,length(x))
  col.grid <- array(col.grid,length(x))

  # UNIT CONVERSIONS
  for(i in 1:length(x))
  {
    # unit conversion
    x[[i]] <- unit.UD(x[[i]],length=dist$scale)

    # ML DENSITY PLOTS
    plot_df(x[[i]],DF=DF,col=col.UD[[i]],...)
  }

  if(DF %nin% c("PDF","CDF")) { return(invisible(NULL)) } # NPR

  # plot grid
  for(i in 1:length(x))
  {
    if("H" %in% names(x[[i]]) && sum(diag(x[[i]]$H)>0))
    {
      H <- covm(x[[i]]$H)
      theta <- H@par["angle"]
      sigma <- eigenvalues.covm(H)

      X <- x[[i]]$r$x
      Y <- x[[i]]$r$y

      COS <- cos(theta)
      SIN <- sin(theta)

      # grid spacing
      du <- sqrt(sigma[1])
      dv <- sqrt(sigma[2])

      # bandwidth axes
      u <- outer(+X*COS,+Y*SIN,"+")
      v <- outer(-X*SIN,+Y*COS,"+")

      # extent of data
      B <- (x[[i]]$PDF > 0)
      if(any(B))
      {
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
    }

    # not sure why this is necessary
    graphics::box(lwd=lwd.level,...)
  }

  # CONTOURS
  for(i in 1:length(x))
  {
    if(!any(is.na(col.level[i,,])) && !any(is.na(level.UD)))
    {
      # make sure that correct style is used for low,ML,high even in absence of lows and highs
      plot_kde(x[[i]],level=level.UD,labels=labels[i,,2],col=malpha(col.level[i,,2],1),lwd=lwd.level,convex=convex,...)

      if(!is.na(level) && !is.null(x[[i]]$DOF.area))
      {
        P <- sapply(level.UD, function(l) { CI.UD(x[[i]],l,level,P=TRUE)[-2] } )
        plot_kde(x[[i]],level=P,labels=c(t(labels[i,,c(1,3)])),col=malpha(c(t(col.level[i,,c(1,3)])),0.5),lwd=lwd.level/2,convex=convex,...)
      }
    }
  }

}
plot_RS <- plot.UD

##################################
# plot PDF stored as KDE object
plot_df <- function(kde,DF="CDF",col="blue",...)
{
  alpha <- grDevices::col2rgb(col,alpha=TRUE)[4,]
  alpha <- max(alpha)
  alpha <- min(alpha,254) # overflow bug otherwise
  col <- malpha(col,(0:alpha)/255)

  if(DF %in% "PDF")
  {
    zlim <- c(0,max(kde[[DF]]))
  }
  else if(DF=="CDF")
  {
    zlim <- c(0,1)
    kde[[DF]] <- 1 - kde[[DF]]
  }
  else # NPR
  {
    zlim <- range(kde[[DF]],na.rm=TRUE)
  }

  # imageRaster is faster but is not reliably called from image
  TEST <- try( graphics::image(kde$r,z=kde[[DF]],useRaster=TRUE,zlim=zlim,col=col,add=TRUE,...) )
  if(!is.null(TEST)) { graphics::image(kde$r,z=kde[[DF]],useRaster=FALSE,zlim=zlim,col=col,add=TRUE,...) }
}


#############################
# Plot a KDE object's contours
plot_kde <- function(kde,level=0.95,labels=round(level*100),col="black",convex=FALSE,...)
{
  # record current option
  # MAX <- getOption("max.contour.segments")

  # do something that works
  drawlabels <- !(labels==FALSE | is.na(labels))
  if(!convex)
  {
    options(max.contour.segments=.Machine$integer.max)
    graphics::contour(kde$r,z=kde$CDF,levels=level,labels=labels,labelcex=1,drawlabels=drawlabels,col=col,add=TRUE,...)
  }
  else
  {
    for(i in 1:length(level))
    {
      SP <- convex(kde,level=level[i]) # now spatial polygons
      sp::plot(SP,border=col,add=TRUE,...)
    }
  }

  # reinstate initial option (or default if was NULL--can't set back to NULL???)
  # if(is.null(MAX)) { MAX <- 25000 }
  # options(max.contour.segments=MAX)
}


##############################
# Plot Gaussian ctmm contours
plot_ctmm <- function(model,alpha=0.05,col="blue",bg=NA,...)
{
  mu <- model$mu # mean vector
  sigma <- model$sigma # covariance matrix

  lapply(alpha,function(a){ellipsograph(mu=mu,sigma=sigma,level=1-a,fg=col,bg=bg,...)})
}

###################
# mu - mean vector
# sigma - covariance matrix
ellipsograph <- function(mu,sigma,level=0.95,fg=graphics::par("col"),bg=NA,PLOT=TRUE,...)
{
  Eigen <- eigen(sigma)
  std <- Eigen$values
  std[1] <- clamp(std[1],0,Inf)
  std[2] <- clamp(std[2],0,std[1])
  std <- sqrt(std)
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

  if(PLOT) { graphics::xspline(x, y=y, shape=-1, open=FALSE, border=fg, col=bg, ...) }
  else { return( cbind(x=x,y=y) ) }
}
