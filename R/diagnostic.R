diagnostic <- function(data,UD,RSF,R=list(),variable=NULL,level=0.95,smooth=TRUE,interpolate=TRUE,...)
{
  # how to sample rasters
  interpolate <- rep(interpolate,length(R))
  interpolate <- ifelse(interpolate,"bilinear","simple")
  names(interpolate) <- names(R)

  axes <- RSF$axes
  GEO <- c('longitude','latitude')
  CTMM <- UD@CTMM
  if(smooth && any(CTMM$error>0)) { data[,c(axes,GEO)] <- predict(data,CTMM=CTMM,t=data$t,complete=TRUE)[,c(axes,GEO)] }

  RVARS <- names(R)
  formula <- RSF$formula

  if(is.null(variable))
  {
    if(length(RVARS)==1)
    { variable <- RVARS[1] }
    else
    {
      for(r in RVARS)
      {
        diagnostic(data,UD=UD,RSF=RSF,R=R,variable=r,level=level,smooth=smooth,interpolate=interpolate,...)
        return()
      }
    }
  }

  # rlim=FALSE
  # if(length(rlim)==1)
  # {
  #   if(rlim)
  #   { rlim <- raster::cellStats() }
  #   else
  #   { rlim <- c(-Inf,Inf) }
  # }

  weights <- UD$weights

  # evaluate raster data
  for(r in names(R))
  {
    xy <- get.telemetry(data,GEO)
    PROJ <- raster::projection(R[[r]])
    xy <- project(xy,to=PROJ)
    data[[r]] <- raster::extract(R[[r]],xy,method=interpolate[r])
  }
  rm(xy)

  VARS <- all.vars(formula)
  DVARS <- VARS[ VARS %nin% RVARS ]
  for(D in DVARS) { data[[D]] <- as.numeric(data[[D]]) } # model.matrix will rename otherwise

  DATA <- data.frame(data)
  DATA$x <- DATA$x - RSF$mu[1]
  DATA$y <- DATA$y - RSF$mu[2]
  log.p <- -(DATA$x^2+DATA$y^2)/(2*RSF$sigma@par['major']) # isotropic

  MODEL <- stats::model.matrix(formula,data=DATA)
  # but what terms change with 'variable'?
  DATA[[variable]] <- 0
  RM <- stats::model.matrix(formula,data=DATA)
  RM <- abs(RM-MODEL)
  RM <- apply(RM,2,max)
  RM <- RM>.Machine$double.eps
  RM <- colnames(MODEL)[RM] # these terms all change with 'variable'
  # there isn't an easier way to do this?

  VARS <- names(RSF$beta)
  VARS <- VARS[VARS %nin% RM] # only vars to condition on
  if(length(VARS)) { log.p <- log.p + c(MODEL[,VARS,drop=FALSE] %*% RSF$beta[VARS]) }
  p <- exp(log.p)

  axes <- variable
  weights <- weights * p
  weights <- weights / sum(weights)
  error <- 0.001
  res <- 10
  n <- UD$DOF.H[1]
  w2d <- 1/n
  w2o <- (n-1)/n
  MISE <- function(h) { w2d/sqrt(2*h^2) + w2o/sqrt(2+2*h^2) - 2/sqrt(2+h^2) + 1/sqrt(2) }
  # silverman's rule of thumb
  DIM <- 1
  h <- (4/(DIM+2)/n)^(1/(DIM+4))
  # find 1D bandwidth for same n
  control <- list(precision=.Machine$double.eps^(1/4))
  h <- optimizer(h,MISE,lower=0,control=control)$par
  # find bias for same h
  bias <- 1 - 1/n + h^2
  # integral K^2
  R <- 1/sqrt(4*pi)
  RFIT <- ctmm.fit(data,ctmm(axes=axes))
  H <- h^2 * methods::getDataPart(RFIT$sigma)
  EXT <- extent(RFIT,level=1-error)[,axes,drop=FALSE] # Gaussian extent (includes uncertainty)
  dr <- c(sqrt(H)/res)
  grid <- format.grid(NULL,axes=axes)
  grid <- kde.grid(data,H=H,axes=axes,alpha=error,res=res,dr=dr,grid=grid,EXT.min=EXT)
  KDE <- kde(data,H=H,axes=axes,CTMM=RFIT,bias=bias,W=weights,alpha=error,dr=dr,grid=grid,...)

  alpha <- 1 - level
  z <- stats::qnorm(1-alpha/2)

  RANGE <- range(data[[variable]])
  r <- KDE$r[[1]]
  dr <- KDE$dr
  SUB <- r>=RANGE[1]-dr & r<=RANGE[2]+dr
  r <- r[SUB]
  PDF <- KDE$PDF[SUB]
  log.PDF <- log(PDF)
  ZERO <- stats::approx(r,log.PDF,RANGE[1],rule=2)$y
  log.PDF <- log.PDF - ZERO
  VAR.log.PDF <- c(R/(n*sqrt(H))) / PDF
  SE.log.PDF <- z * sqrt(VAR.log.PDF)

  plot(range(r),range(log.PDF-SE.log.PDF,log.PDF+SE.log.PDF),xlab=variable,ylab="log(\u03BB)",col=grDevices::rgb(1,1,1,0))
  graphics::points(r,log.PDF,type='l',lwd=2)
  graphics::polygon(c(r,rev(r)),c(log.PDF-SE.log.PDF,rev(log.PDF+SE.log.PDF)),border=NA,col=malpha('black',0.25))

  VARS <- names(RSF$beta)
  VARS <- VARS[VARS %in% RM] # only vars to condition on
  if(length(VARS))
  {
    r <- data[[variable]]
    IND <- sort(r,index.return=TRUE)$ix
    r <- r[IND]
    MODEL <- MODEL[IND,]

    # # minimize CI width via normalization
    # MID <- mean(RANGE)
    # for(V in VARS)
    # {
    #   AVE <- stats::approx(r,MODEL[,V],ties=mean)$y
    #   MODEL[,V] <- MODEL[,V] - AVE
    # }

    EST <- MODEL[,VARS,drop=FALSE] %*% RSF$beta[VARS]
    VAR <- sapply(1:nrow(MODEL),function(i){ MODEL[i,VARS,drop=FALSE] %*% RSF$COV[VARS,VARS,drop=FALSE] %*% t(MODEL[i,VARS,drop=FALSE]) })
    SE <- z*sqrt(VAR)

    # start at zero
    ZERO <- stats::approx(r,EST,RANGE[1],rule=2,ties=mean)$y
    EST <- EST - ZERO # make this exact with model matrix and subtract before VAR calculation !!!!!!!!!!!!!!!!!

    graphics::points(r,EST,col='red',type='l',lwd=2)
    graphics::polygon(c(r,rev(r)),c(EST-SE,rev(EST+SE)),border=NA,col=malpha('red',0.25))
  }
}
