intensity <- function(data,UD,RSF,R=list(),variable=NULL,empirical=FALSE,level=0.95,ticks=TRUE,smooth=TRUE,interpolate=TRUE,...)
{
  # rlim=NULL
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
        intensity(data,UD=UD,RSF=RSF,R=R,variable=r,level=level,smooth=smooth,interpolate=interpolate,...)
        return()
      }
    }
  }

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
  w <- weights * p # WHY??????????
  w <- w / sum(w)
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
  RK2 <- 1/sqrt(4*pi)
  RFIT <- ctmm.fit(data,ctmm(axes=axes))
  H <- h^2 * methods::getDataPart(RFIT$sigma)
  EXT <- extent(RFIT,level=1-error)[,axes,drop=FALSE] # Gaussian extent (includes uncertainty)
  dr <- c(sqrt(H)/res)
  # if(!is.null(rlim))
  # {
  #   # fix dr to fit in rlim
  #   RANGE <- diff(rlim)
  #   dr <- RANGE/ceiling(RANGE/dr)
  #   # format extent to include margin
  #   rlim <- cbind(rlim)
  #   colnames(rlim) <- variable
  #   rownames(rlim) <- c('min','max')
  #   grid <- list(extent=rlim)
  # }
  # else
  { grid <- NULL }
  grid <- format_grid(grid,axes=axes)
  grid <- kde.grid(data,H=H,axes=axes,alpha=error,res=res,dr=dr,grid=grid,EXT.min=EXT)
  KDE <- kde(data,H=H,axes=axes,CTMM=RFIT,bias=bias,W=w,alpha=error,dr=dr,grid=grid,...)

  alpha <- 1 - level
  z <- stats::qnorm(1-alpha/2)

  R <- data[[variable]]
  if(length(variable))
  {
    IND <- sort(R,index.return=TRUE)$ix
    R <- R[IND]
    MODEL <- MODEL[IND,]
    weights <- weights[IND]

    VARS <- names(RSF$beta)
    for(V in VARS)
    {
      # derived from normalization and MLE
      AVE <- c(weights %*% MODEL[,V])
      MODEL[,V] <- MODEL[,V] - AVE
    }

    EST <- c(MODEL[,VARS,drop=FALSE] %*% RSF$beta[VARS])
    VAR <- sapply(1:nrow(MODEL),function(i){ MODEL[i,VARS,drop=FALSE] %*% RSF$COV[VARS,VARS,drop=FALSE] %*% t(MODEL[i,VARS,drop=FALSE]) })
    SE <- z*sqrt(VAR)

    # location of minimum uncertainty
    MIN <- which.min(SE)
    EST <- EST - EST[MIN]
  }

  RANGE <- range(data[[variable]])
  r <- KDE$r[[1]]
  dr <- KDE$dr
  SUB <- r>=RANGE[1]-dr & r<=RANGE[2]+dr
  r <- r[SUB]
  PDF <- KDE$PDF[SUB]
  log.PDF <- log(PDF)
  # matches below for linear model
  if(length(variable)) { ZERO <- stats::approx(r,log.PDF,R[MIN],rule=2,ties=mean)$y }
  # the above can fail with no PDF support at R[MIN]
  if(!length(variable) || is.na(ZERO))
  {
    MIN <- which.max(PDF) # min PDF uncertainty
    ZERO <- log.PDF[MIN]
  }
  log.PDF <- log.PDF - ZERO
  VAR.log.PDF <- c(RK2/n/sqrt(H)) / PDF
  SE.log.PDF <- z * sqrt(VAR.log.PDF)

  lRANGE <- NULL
  if(empirical)
  {
    lRANGE <- c(pmax(log.PDF-SE.log.PDF,ifelse(log.PDF<0,2*log.PDF,-log.PDF)),pmin(log.PDF+SE.log.PDF,ifelse(log.PDF>0,2*log.PDF,-log.PDF)))
    lRANGE <- lRANGE[abs(lRANGE)<Inf]
    lRANGE <- range(lRANGE,na.rm=TRUE)
  }
  if(length(variable))
  {
    lRANGE <- range(lRANGE,pmin(EST-SE,ifelse(EST<0,2*EST,-EST)),pmin(EST+SE,ifelse(EST>0,2*EST,-EST)))
    lRANGE <- range(lRANGE,na.rm=TRUE)
  }
  ylab <- paste0("log(\u03BB)")
  plot(RANGE,lRANGE,xlab=variable,ylab=ylab,col=grDevices::rgb(1,1,1,0))

  if(empirical)
  {
    graphics::points(r,log.PDF,type='l',lwd=2)
    graphics::polygon(c(r,rev(r)),c(log.PDF-SE.log.PDF,rev(log.PDF+SE.log.PDF)),border=NA,col=malpha('black',0.25))
  }

  if(length(variable))
  {
    graphics::points(R,EST,col='red',type='l',lwd=2)
    graphics::polygon(c(R,rev(R)),c(EST-SE,rev(EST+SE)),border=NA,col=malpha('red',0.25))
  }

  # ticks at top corresponding to sampled resource values?
  if(ticks) { graphics::axis(side=3,at=R,labels=FALSE,col=malpha('black',0.5)) }
}
