rsf.fit <- function(data,UD,beta=NULL,R=list(),formula=NULL,integrated=TRUE,reference="auto",level.UD=0.99,isotropic=TRUE,debias=TRUE,smooth=TRUE,standardize=TRUE,error=0.01,trace=TRUE,...)
{
  STATIONARY <- TRUE
  CTMM <- UD@CTMM
  axes <- CTMM$axes
  GEO <- c('longitude','latitude')

  # control list for optimizer
  control <- list()
  # pass trace argument (demoted)
  if(trace) { control$trace <- trace-1 }

  # smooth the data, but don't drop
  if(smooth && any(CTMM$error>0))
  { data[,c(axes,GEO)] <- predict(data,CTMM=CTMM,t=data$t,complete=TRUE)[,c(axes,GEO)] }
  n <- nrow(data)

  if(!integrated) # prepare available region polygon
  {
    if(class(level.UD)[1]=="numeric")
    {
      level.UD <- SpatialPolygonsDataFrame.UD(UD,level.UD=level.UD)
      # subset to point estimate contour
      level.UD <- level.UD@polygons[[2]]
    }
    else # project polygon to data's
    { level.UD <- sp::spTransform(level.UD,sp::CRS(projection(data))) }
    AREA <- level.UD@area
  }
  else if(isotropic && !CTMM$isotropic) # keep things simple for now
  {
    message('RSF code is isotropic for the moment.')
    CTMM <- simplify.ctmm(CTMM,'minor')
    CTMM <- ctmm.fit(data,CTMM,trace=trace)
  }

  # for simulation - uncorrelated, error-less
  IID <- CTMM
  IID$tau <- NULL
  IID$omega <- FALSE
  IID$error <- FALSE

  # extract weights
  W <- mean(UD$DOF.area) # +1 for mean not being detrended accounted for
  w <- UD$weights * W

  R <- listify(R)
  TEST <- sapply(R,raster::inMemory)
  TEST <- any(!TEST)
  if(TEST) { message('Some rasters are not loaded in RAM and may be slow to process. See help("raster::readAll").') }

  if(length(R) && is.null(names(R)))
  {
    warning("R is not a named list of rasters.")
    names(R) <- paste("R",1:length(R))
  }

  if(!is.null(formula) && standardize)
  {
    message("Users are responsible for standardizing ",names(R)," when formula argument is supplied.")
    standardize <- FALSE
  }
  else
  {
    RSCALE <- array(1,length(R))
    names(RSCALE) <- names(R)
  }

  # expand any categorical variables into indicator variables
  STUFF <- expand.factors(R=R,formula=formula,reference=reference,data=data)
  R <- STUFF$R
  formula <- STUFF$formula

  ## The basic setup is:
  # RVARS - raster variables
  # DVARS - data annotation variables - overridden by RVARS
  # TERMS - all predictors - some may be RVARS
  # CVARS - composite variables in TERMS that are single RVARS || DVARS
  # OFFSET - offset T/F variables

  RVARS <- names(R)

  if(is.null(formula))
  {
    TERMS <- RVARS
    DVARS <- NULL
    CVARS <- NULL
    OFFSET <- NULL

    if(length(R))
    {
      formula <- paste(TERMS,collapse="+")
      formula <- paste("~",formula)
      formula <- formula(formula)
    }
    else
    { formula <- formula(~0) }
  }
  else
  {
    VARS <- all.vars(formula)
    DVARS <- VARS[ VARS %nin% RVARS ]

    TERMS <- attr(stats::terms(formula),"term.labels")

    CVARS <- TERMS[ TERMS %nin% VARS ] # remaining terms that are not simple variables

    OFFSET <- stats::terms(formula)
    OFFSET <- attr(OFFSET,"variables")[ attr(OFFSET,"offset") ]

    if(attr(stats::terms(formula),"response")[1]>0) { stop("Response variable not yet supported.") }
  }
  VARS <- c(RVARS,CVARS) # vars that need to be recorded per location
  if(length(DVARS)) { STATIONARY <- FALSE }

  INTERCEPT <- DVARS %in% TERMS
  if(any(INTERCEPT))
  {
    TERMS <- TERMS[ TERMS %nin% DVARS ]
    warning(sum(INTERCEPT),"+",attr(stats::terms(formula),"intercept")," intercept terms ignored.")
  }

  # evaluate term on data.frame
  evaluate <- function(term,envir)
  {
    if(length(dim(envir))==2)
    {
      term <- gsub(":","*",term) # multiplication
      # term <- gsub("I(","(",term) # function evaluations # don't think this is necessary
      envir <- data.frame(envir) # matrices can't be environments, but data.frames can't matrix multiply...
      if(!STATIONARY) { envir <- cbind(envir,data) }
      RET <- eval(parse(text=term),envir=envir)
    }
    else # loop over runs
    {
      if(STATIONARY) # time doesn't matter
      {
        DIM <- dim(envir)
        dim(envir) <- c(DIM[1]*DIM[2],DIM[3])
        colnames(envir) <- VARS
        RET <- evaluate(term,envir)
        dim(RET) <- DIM[1:2]
      }
      else # !STATIONARY
      {
        RET <- sapply(1:dim(envir)[1],function(i)
                {
                  ENVIR <- envir[i,,];
                  dim(ENVIR) <- dim(envir)[-1];
                  colnames(ENVIR) <- VARS
                  ENVIR <- cbind(ENVIR,data);
                  evaluate(term,ENVIR)
                }) # [time,track]
        dim(RET) <- dim(envir)[2:1]
        RET <- t(RET) # [track,time]
      }
    }
    return(RET)
  }

  ## prepare raster data ##
  # I would like to save with raw raster objects to save memory
  # but raster::getValuesBlock is strangely slow
  # and rescaling is good for numerics
  PROJ <- X <- Y <- Z <- Z.ind <- list()
  dX <- dY <- dZ <- rep(NA,length(R))
  for(i in 1 %:% length(R))
  {
    PROJ[[i]] <- raster::projection(R[[i]])

    if(standardize)
    {
      # # raster::median is not defined correctly
      # R[[i]] <- R[[i]] - raster::median(R[[i]],na.rm=TRUE)
      # RSCALE[i] <- raster::median(abs(R[[i]]),na.rm=TRUE)
      R[[i]] <- R[[i]] - raster::cellStats(R[[i]],'mean',na.rm=TRUE)
      RSCALE[i] <- raster::cellStats(R[[i]],'sd',na.rm=TRUE)
      R[[i]] <- R[[i]]/RSCALE[i]
    }

    DIM <- dim(R[[i]])
    # if(length(dim(R[[i]]))==2) # RasterLayer (static)
    # {
    #   R[[i]] <- aperm(R[[i]],2:1) # [x,y]
    # }
    if(DIM[3]>1) # RasterStack or RasterBrick
    {
      X[[i]] <- raster::xFromCol(R[[i]],1:DIM[2])
      Y[[i]] <- raster::yFromRow(R[[i]],1:DIM[1])
      Z[[i]] <- raster::getZ(R[[i]])

      dX[i] <- stats::median(diff(X[[i]]))
      dY[i] <- stats::median(diff(Y[[i]]))

      R[[i]] <- raster::as.array(R[[i]])[,,]

      STATIONARY <- FALSE
      R[[i]] <- aperm(R[[i]],c(2,1,3)) # [x,y,z]

      if(class(Z[[i]])[1]=="Date") { Z[[i]] <- as.POSIXct(Z[[i]]) } # numeric is in days
      if(class(Z[[i]])[1]=="POSIXct") { Z[[i]] <- as.numeric(Z[[i]]) }

      dZ[i] <- stats::median(diff(Z[[i]]))

      # index for data times in raster stack
      Z.ind[[i]] <- (data$t - Z[[i]][1])/dZ[i] + 1
    } # end XYZ
  } # end R loop
  if(length(X)) { names(X) <- names(Y) <- names(dX) <- names(dY) <- RVARS[1:length(X)] }
  if(length(Z)) { names(Z) <- names(dZ) <- names(Z.ind) <- RVARS[1:length(Z)] }

  # setup integrated spatial covariates
  if(!integrated)
  {
    SVARS <- NULL
    NSPA <- 0
  }
  else
  {
    mu <- c(CTMM$mu)
    std <- sqrt(CTMM$sigma@par[c('major','minor')])
    names(std) <- names(mu) <- axes

    SVARS <- axes
    NSPA <- 2
    SCALE <- std

    if(isotropic)
    {
      SVARS <- c(SVARS,'rr')
      NSPA <- NSPA + 1
      SCALE <- c(SCALE,prod(std))
    }
    else
    {
      SVARS <- c(SVARS,'xx','yy','xy')
      NSPA <- NSPA + 3
      SCALE <- c(SCALE,std[1]^2,std[2]^2,prod(std))
      theta <- CTMM$sigma@par['angle']
    }

    names(SCALE) <- SVARS
  }
  VARS <- c(VARS,SVARS)
  TERMS <- c(TERMS,SVARS)

  # data and simulation sampled covariates
  DATA <- matrix(0,nrow(data),length(VARS))
  colnames(DATA) <- VARS

  # # generic SCALE code for arbitrary formula terms not possible
  # SCALES <- c(RVARS,DVARS)
  # SCALES <- t(SCALES)
  # SCALES <- as.data.frame(SCALES)
  # SCALE <- rep(0,length(TERMS))
  # for(i in 1:length(TERMS)) { SCALE[i] <- evaluate(TERMS[i],SCALES) }

  # initial estimates
  BETA <- rep(0,length(TERMS))
  names(BETA) <- TERMS
  if(!is.null(beta))
  {
    if(class(beta)[1]=="ctmm") { beta <- beta$beta }

    if(is.null(names(beta))) # integrated terms default to 0
    { BETA <- pad(beta,length(TERMS)) }
    else # use relevant initial estimates
    {
      beta <- beta[names(beta) %in% TERMS]
      BETA[names(beta)] <- beta
    }
  }
  beta <- BETA
  names(BETA) <- TERMS
  # beta <- beta * SCALE

  # store raster covariates
  for(i in 1:length(RVARS))
  {
    r <- RVARS[i]
    # store data-sampled covariates
    if(i==1 || PROJ[[i]]!=PROJ[[i-1]])
    {
      xy <- get.telemetry(data,GEO)
      xy <- project(xy,to=PROJ[[i]])
    }

    DIM <- dim(R[[i]])
    if(DIM[3]==1) # 2D
    {
      DATA[,r] <- raster::extract(R[[i]],xy,method="bilinear")
      # DATA[,RI] <- bint(R[[i]],t(xy))
    }
    else # 3D space-time stack
    {
      # convert to indices
      XY <- xy
      XY[,1] <- (XY[,1] - X[[i]][1])/dX[i] + 1
      XY[,2] <- (XY[,2] - Y[[i]][1])/dY[i] + 1

      XY <- cbind(XY,Z.ind[[i]])
      DATA[,r] <- tint(R[[i]],t(XY))
    }
  }

  # calculate and store composite covariates
  for(v in CVARS) { DATA[,v] <- evaluate(v,DATA) }

  ### store sampled integrated terms
  if(integrated)
  {
    # de-trend mean
    DATA[,'x'] <- data$x - mu['x']
    DATA[,'y'] <- data$y - mu['y']

    # rotate to major-minor axes
    if(!isotropic) { DATA[,axes] <- rotate.vec(DATA[,axes],-theta) }

    # standardize
    DATA[,'x'] <- DATA[,'x']/std['x']
    DATA[,'y'] <- DATA[,'y']/std['y']

    # variance/covariance terms (relative to pilot estimate)
    if(isotropic) # beta is correction to standardized 1/sigma
    {
      DATA[,'rr'] <- -( DATA[,'x']^2 + DATA[,'y']^2 )/2
    }
    else # beta is correction to standardized solve(sigma)
    {
      DATA[,'xx'] <- -DATA[,'x']^2 /2
      DATA[,'yy'] <- -DATA[,'y']^2 /2
      DATA[,'xy'] <- -DATA[,'x']*DATA[,'y']
    }
  } # end if(integrated)

  DAVE <- c(w %*% DATA)
  if(any(is.na(DAVE))) { stop("NA values in sampled rasters.") }
  names(DAVE) <- VARS
  SATA <- array(0,c(0,dim(DATA))) # simulated data [track,time,vars]
  dimnames(SATA)[[3]] <- VARS

  # partition function defined (scaled) - calling from rsf.loglike
  # if(!is.null(NORM))
  # {
  #   loglike <- (DAVE[,TERMS,drop=FALSE] %*% beta) - W*NORM # Z includes scaling
  #   if(length(OFFSET))
  #   {
  #     ONE <- rep(TRUE,nrow(data))
  #     for(o in OFFSET) { ONE <- ONE & evaluate(o,DATA) }
  #     if(any(!ONE)) { loglike <- -Inf }
  #   }
  #   # log-normal terms from importance sampling
  #   if(integrated)
  #   {
  #     SUMNORM <- rowSums(DATA[,axes]^2) # [n]
  #     SUMNORM <- - 1/2*c(w %*% SUMNORM) # + log(1) terms
  #     loglike <- loglike + SUMNORM
  #   }
  #   return(loglike)
  # }

  nloglike <- function(beta,zero=0,verbose=FALSE)
  {
    SAMP <- exp(SATA[,,TERMS,drop=FALSE] %.% beta) # [track,time]
    SAMP[] <- nant(SAMP,0) # treat NA as inaccessible region

    if(length(OFFSET))
    {
      ONE <- rep(TRUE,nrow(SATA))
      for(o in OFFSET) { ONE <- ONE & evaluate(o,SATA) }
      SAMP <- ONE * SAMP
    }

    nll <- -c(DAVE[TERMS] %*% beta)

    if(STATIONARY)
    {
      SAMP <- c(SAMP)
      MEAN <- mean(SAMP)
      log.MEAN <- log(MEAN)

      if(debias || verbose) # numerical error variance (per W^2)
      { VAR.log <- stats::var(SAMP)/length(SAMP) /MEAN^2 }
      if(debias) # MEAN-log != log-MEAN bias for small N (per W)
      { log.MEAN <- log.MEAN + W/2*VAR.log } # +1/2 from Taylor series

      nll <- W *(nll/W + log.MEAN - zero/W)

      if(verbose) { VAR.log <- W^2 * VAR.log }
    }
    else # !STATIONARY
    {
      MEAN <- apply(SAMP,2,mean) # [time]
      log.MEAN <- log(MEAN)

      if(debias || verbose) # numerical error variance (per w^2)
      { VAR.log <- apply(SAMP,2,stats::var)/dim(SAMP)[1] /MEAN^2 }
      if(debias) # MEAN-log != log-MEAN bias for small N (per w)
      { log.MEAN <- log.MEAN + w/2*VAR.log } # +1/2 from Taylor series

      nll <- sum(nll/n + w*log.MEAN - zero/n) # can't divide by w

      if(verbose) { VAR.log <- sum(w^2 * VAR.log) }
    }

    nll <- c(nll)
    if(verbose)
    {
      RET <- list(loglike=-nll,Z=log.MEAN,VAR.loglike=VAR.log)
      return(RET)
    }
    else
    { return(nll) }
  }

  N <- ifelse(STATIONARY,1,8) # starting value, will increase iteratively
  N.OLD <- 0
  loglike <- -Inf
  ERROR.BIG <- TRUE
  STDloglike <- Inf
  while(ERROR.BIG)
  {
    # double the random sample
    if(integrated)
    { SIM <- simulate(IID,t=1:(N*nrow(data)),complete=TRUE) }
    else
    {
      SIM <- sp::spsample(level.UD,n=N*nrow(data),type="random")
      # sp::spsample drops projection information and throws annoying warning when trying to fix
      suppressWarnings( sp::proj4string(SIM) <- sp::CRS(projection(data)) )
      SIM <- sp::spTransform(SIM,sp::CRS(DATUM))
      SIM <- SIM@coords
      colnames(SIM) <- GEO
    }
    SATA <- fbind(SATA,array(0,c(N,nrow(data),length(VARS))))
    dimnames(SATA)[[3]] <- VARS
    SUB <- N.OLD+1:N # new indices

    # store simulation-sampled raster covariates
    for(i in 1:length(RVARS))
    {
      r <- RVARS[i]

      # store data-sampled covariates
      if(i==1 || PROJ[[i]]!=PROJ[[i-1]])
      {
        if(integrated)
        { xy <- get.telemetry(SIM,GEO) }
        else
        { xy <- SIM }

        xy <- project(xy,to=PROJ[[i]])
      }

      DIM <- dim(R[[i]])
      if(DIM[3]==1)
      {
        SATA[SUB,,r] <- raster::extract(R[[i]],xy,method="bilinear")
        #SATA[SUB,,r] <- bint(R[[r]],t(xy))
      }
      else
      {
        XY <- xy
        XY[,1] <- (XY[,1] - X[[i]][1])/dX[i] + 1
        XY[,2] <- (XY[,2] - Y[[i]][1])/dY[i] + 1

        XY <- cbind(XY,rep(data$t,N))
        SATA[SUB,,r] <- tint(R[[i]],t(xy))
      }
    }

    for(v in CVARS) { SATA[SUB,,v] <- evaluate(v,SATA[SUB,,,drop=FALSE]) }

    if(integrated)
    {
      # de-trend
      SATA[SUB,,'x'] <- SIM$x - mu['x']
      SATA[SUB,,'y'] <- SIM$y - mu['y']

      # rotate
      if(!isotropic)
      { SATA[SUB,,axes] <- rotate.vec(SATA[SUB,axes],-theta) }

      # standardize
      SATA[SUB,,'x'] <- SATA[SUB,,'x']/std['x']
      SATA[SUB,,'y'] <- SATA[SUB,,'y']/std['y']

      # variance/covariance terms
      if(isotropic)
      {
        SATA[SUB,,'rr'] <- -( SATA[SUB,,'x']^2 + SATA[SUB,,'y']^2 )/2
      }
      else
      {
        SATA[SUB,,'xx'] <- -SATA[SUB,,'x']^2/2
        SATA[SUB,,'yy'] <- -SATA[SUB,,'y']^2/2
        SATA[SUB,,'xy'] <- -SATA[SUB,,'x']*SATA[SUB,,'y']
      }
    }

    # update SIM count
    N <- nrow(SATA)

    if(trace) { message("Maximizing likelihood @ n=",N) }
    # adaptive precision - error target is for parameters or sqrt(STDloglike)
    precision <- clamp(STDloglike,0,1)/sqrt(2) # expected STD[loglike] this round
    precision <- max(precision,error) # no point in exceeding threshold error
    precision <- precision^2 # improve over Monte-Carlo error
    precision <- log(precision)/log(.Machine$double.eps) # relative to machine error
    precision <- clamp(precision,1/8,1/2)
    control$precision <- precision

    RESULT <- optimizer(beta,nloglike,control=control)
    beta.OLD <- beta
    beta <- RESULT$par
    loglike.OLD <- loglike
    loglike <- -RESULT$value
    COV <- RESULT$covariance

    STUFF <- nloglike(beta,zero=-loglike,verbose=TRUE)
    VAR.loglike <- STUFF$VAR.loglike
    STDloglike <- sqrt(VAR.loglike)
    STDbeta <- (beta-beta.OLD)/sqrt(diag(COV))
    if(trace)
    {
      message(" SD[log(\u2113)] = ",STDloglike)
      message(" dlog(\u2113) = ",loglike-loglike.OLD)
      message(" d\u03B2/SD[\u03B2] = ",paste(STDbeta,collapse=' '))
    }

    # numbers for next iteration
    N.OLD <- N

    # estimated errors for this round
    ERROR.BIG <- c(STDloglike,abs(c(loglike-loglike.OLD,STDbeta))/sqrt(2)) > error
    ERROR.BIG <- any(ERROR.BIG)
  } # while(ERROR.BIG)

  # compute hessian
  if(trace) { message("Calculating Hessian") }
  DIFF <- genD(par=beta,fn=nloglike,zero=-loglike,Richardson=2,mc.cores=1)
  hess <- DIFF$hessian
  grad <- DIFF$gradient
  # more robust covariance calculation than straight inverse
  COV <- cov.loglike(hess,grad)
  dimnames(COV) <- list(TERMS,TERMS)

  if(integrated)
  {
    # fixed log-normal terms from importance sampling
    SUMNORM <- apply(DATA[,axes,drop=FALSE]^2,1,sum) # [time]
    SUMNORM <- - 1/2*c(w %*% SUMNORM) # + log(1) terms
    loglike <- loglike + SUMNORM

    # un-shift precision matrices # mu was zero centered
    if(isotropic)
    {
      beta['rr'] <- 1 + beta['rr']
    }
    else
    {
      beta['xx'] <- 1 + beta['xx']
      beta['yy'] <- 1 + beta['yy']
      beta['xy'] <- 0 + beta['xy']
    }
  }

  # This doesn't account for bias correction
  # partition function (per W)
  # STUFF <- nloglike(beta,verbose=TRUE)
  # NORM <- STUFF$Z

  if(integrated)
  {
    # un-scale mu, sigma, COV[...]
    beta[SVARS] <- beta[SVARS] / SCALE
    COV[SVARS,] <- COV[SVARS,] / SCALE
    COV[,SVARS] <- COV[,SVARS] / SCALE

    # log-likelihood adjustment (unscale)
    loglike <- loglike - W*log(prod(std)) # 1/2 from sqrt
    # partition function adjustment (unscale)
    # NORM <- NORM + log(prod(std)) # per W

    # convert spatial parameters to adjusted mean and covariance estimates
    fn <- function(beta)
    {
      if(isotropic)
      {
        # beta['rr'] is precision
        beta['rr'] <- 1/beta['rr']

        # beta[c('x','y')] is correction to precision*mu
        beta[axes] <- beta['rr']*beta[axes] + c(CTMM$mu)
      }
      else
      {
        sigma <- matrix(beta[c('xx','xy','xy','yy')],2,2)
        sigma <- PDsolve(sigma)

        beta[c('xx','yy','xy')] <- covm(sigma)@par # (major,minor,angle)

        beta[axes] <- c(sigma %*% beta[axes]) + c(CTMM$mu)
      }

      return(beta)
    }

    J <- genD(par=beta,fn=fn,order=1,drop=FALSE)$gradient
    beta <- fn(beta)
    COV <- J %*% COV %*% t(J)

    if(!isotropic) { beta['xy'] <- beta['xy'] + CTMM$sigma@par['angle'] }

    TERMS <- TERMS[ TERMS %nin% SVARS ]
    SVARS <- axes
    if(isotropic)
    { SVARS <- c(SVARS,"major") }
    else
    { SVARS <- c(SVARS,"major","minor","angle") }
    TERMS <- c(TERMS,SVARS)
  } # end if(integrated)
  else
  {
    loglike <- loglike - W*log(AREA)
    # NORM <- NORM + log(AREA)
  }

  names(beta) <- TERMS
  dimnames(COV) <- list(TERMS,TERMS)

  # debias code (variance only)
  if(integrated && debias)
  {
    if(!is.null(CTMM$MLE))
    { DEBIAS <- sqrt(det.covm(CTMM$sigma)/det.covm(CTMM$MLE$sigma)) }
    else
    { DEBIAS <- W/(W-1) }

    DSCALE <- rep(1,length(TERMS))
    names(DSCALE) <- TERMS
    DSCALE['major'] <- DEBIAS
    if(!isotropic) { DSCALE['minor'] <- DEBIAS }

    beta <- beta * DSCALE
    COV <- COV * outer(DSCALE)

    # impact on RSF variables missing
  }

  ## convert to ctmm() object format with extra beta slot & features
  if(integrated)
  {
    mu <- beta[axes]
    COV.mu <- COV[axes,axes]

    beta <- rm.name(beta,axes)
    COV <- rm.name(COV,axes)

    if(isotropic)
    { PAR <- 'major' }
    else
    { PAR <- c('major','minor','angle') }
    sigma <- beta[PAR]
    sigma <- covm(sigma,isotropic=isotropic,axes=axes)
    beta <- rm.name(beta,PAR)
  }
  else
  {
    mu <- CTMM$mu
    sigma <- covm(Inf,isotropic=isotropic,axes=axes)
    COV.mu <- diag(Inf,length(axes))
    dimnames(COV.mu) <- list(axes,axes)
  }

  if(standardize)
  {
    beta <- beta/RSCALE
    COV[RVARS,RVARS] <- COV[RVARS,RVARS] / outer(RSCALE)
  }

  # package results and return
  # turn this into ctmm object
  RSF <- ctmm(axes=axes,mu=mu,COV.mu=COV.mu,beta=beta,sigma=sigma,isotropic=isotropic,COV=COV)
  RSF@info <- CTMM@info
  RSF$loglike <- loglike
  RSF$VAR.loglike <- VAR.loglike
  # RSF$Z <- NORM
  RSF$integrated <- integrated
  RSF$formula <- formula
  if(!integrated) { RSF$level.UD <- level.UD }

  # copy over autocorrelation information in a reasonable way
  if(integrated)
  {
    RSF$tau <- CTMM$tau
    RSF$omega <- CTMM$omega
    RSF$error <- CTMM$error

    VAR <- diag(CTMM$COV)
    COR <- stats::cov2cor(CTMM$COV)
    NEW <- names(VAR)[names(VAR) %nin% TERMS]
    OLD <- names(VAR)[names(VAR) %in% TERMS]
    ALL <- c(colnames(COV),NEW)
    DIM <- length(ALL)

    BLANK <- matrix(0,DIM,DIM)
    dimnames(BLANK) <- list(ALL,ALL)
    # copy over RSF stuff
    BLANK[rownames(COV),colnames(COV)] <- COV

    # adjust old parameters to equal new parameters, while preserving correlations
    VAR[OLD] <- diag(COV)[OLD]
    STD <- sqrt(VAR)
    COV.OLD <- COR * outer(STD)
    BLANK[NEW,NEW] <- COV.OLD[NEW,NEW]
    BLANK[NEW,OLD] <- COV.OLD[NEW,OLD]
    BLANK[OLD,NEW] <- COV.OLD[OLD,NEW]

    RSF$COV <- BLANK
    RSF$features <- rownames(RSF$COV)
  }
  else # !integrated
  {
    if(length(beta)) { dimnames(RSF$COV) <- list(TERMS,TERMS) }
    RSF$features <- TERMS
    RSF$range <- FALSE
  }

  # check if some data fell outside of polygon
  if(!integrated)
  {
    # if any data isn't inside polygon then loglike = -Inf
    level.UD <- level.UD@Polygons[[1]]@coords
    TEST <- sp::point.in.polygon(data$x,data$y,level.UD[,1],level.UD[,2])
    TEST <- sum(which(TEST==0)) # number of exterior points
    if(TEST>0)
    {
      warning(TEST," data points outside of available area.")
      RSF$loglike <- -Inf # P(outside)=0 -> log(P)=-Inf
    }
  }

  # AIC,AICc,BIC
  RSF$method <- CTMM$method # not sure if this is the best idea
  RSF <- ic.ctmm(RSF,W) # don't include mu, sigma
  RSF$range <- TRUE

  return(RSF)
}


# expand raster factors into
expand.factors <- function(R,formula,reference="auto",data=NULL,fixed=FALSE)
{
  NAMES <- names(R)
  for(i in 1:length(R))
  {
    if(raster::is.factor(R[[i]]))
    {
      NAME <- NAMES[i]
      FACT <- R[[i]]
      R <- R[-i]
      NAMES <- NAMES[-i]

      LEVELS <- raster::levels(FACT)[[1]]$ID # assuming one layer

      if(fixed)
      {
        TERMS <- attr(stats::terms(formula),"term.labels")
        # I am not good with regexp
        HEAD <- paste0(NAME,"[")
        reference <- TERMS[grepl(HEAD,TERMS,fixed=TRUE)][1] # *NAME[#/ref]*
        reference <- strsplit(reference,HEAD,fixed=TRUE)[[1]][2] # #/ref]*
        reference <- strsplit(reference,"/",fixed=TRUE)[[1]][2] # ref]*
        reference <- strsplit(reference,"]",fixed=TRUE)[[1]][1] # ref
        reference <- which(LEVELS==reference)
      }
      else if(reference=="auto") # fix base layer
      {
        PROJ <- raster::projection(FACT)
        data <- get.telemetry(data,c("longitude","latitude"))
        colnames(data) <- c('x','y')
        data <- project(data,to=PROJ)
        data <- raster::extract(FACT,data) # don't interpolate
        UNIQUE <- unique(data) # may be missing some levels
        data <- tabulate(match(data,UNIQUE)) # tallies
        reference <- UNIQUE[which.max(data)]
        message("Reference category set to ",reference,".")
      }
      # else the reference category is specified by `reference`

      DIFF <- LEVELS %nin% reference
      FACT <- lapply(LEVELS[DIFF],function(l){FACT==l})
      LEVELS <- paste0(NAME,"[",LEVELS,"/",LEVELS[reference],"]")[DIFF]
      names(FACT) <- LEVELS

      R <- c(R,FACT)

      # expand terms
      if(!fixed && !is.null(formula))
      {
        formula <- as.character(formula)[2]
        formula <- sapply(LEVELS,function(l){gsub(NAME,l,formula,fixed=TRUE)})
        formula <- paste(formula,collapse="+")
        formula <- paste("~",formula)
        formula <- eval(parse(text=formula))
        formula <- simplify.formula(formula)
      }
    }
  }
  if(fixed) { return(R) }
  RETURN <- list(R=R,formula=formula)
}


# prepare raster data for processing
R.prepare <- function(R)
{
  PROJ <- raster::projection(R)

  DIM <- dim(R)
  X <- raster::xFromCol(R,1:DIM[2])
  Y <- raster::yFromRow(R,1:DIM[1])
  Z <- raster::getZ(R)

  dX <- stats::median(diff(X))
  dY <- stats::median(diff(Y))

  R <- raster::as.array(R)[,,] # last dim will be dropped if length-1
  if(length(dim(R))==2) # RasterLayer (static)
  {
    R <- aperm(R,2:1)
    dZ <- NULL
  } # [x,y]
  else if (length(dim(R))==3) # RasterStack or RasterBrick
  {
    STATIONARY <- FALSE
    R <- aperm(R,c(2,1,3)) # [x,y,z]

    if(class(Z)[1]=="Date") { Z <- as.POSIXct(Z) } # numeric is in days
    if(class(Z)[1]=="POSIXct") { Z <- as.numeric(Z) }

    dZ <- stats::median(diff(Z))
  }

  R <- list(R=R,PROJ=PROJ,X=X,Y=Y,Z=Z,dX=dX,dY=dY,dZ=dZ)
  return(R)
}


# sample location xy (in projection proj) from raster R with grid coordinates X,Y in (in projection PROJ)
R.extract <- function(xy,proj,R,X,Y,Z=NULL,PROJ,dX,dY,dZ=NULL)
{
  DIM <- dim(R)

  xy <- project(xy,from=proj,to=PROJ)
  # continuous index
  xy[,1] <- (xy[,1] - X[1])/dX + 1
  xy[,2] <- (xy[,2] - Y[1])/dY + 1

  if(length(DIM)==2)
  {
    E <- bint(R,t(xy))
  }
  else # xyt
  {
    # missing t axis
    # this is not fully coded yet, but it is not used either
    E <- tint(R,t(xy))
  }

  return(E)
}

# Evaluate raster on new spatial grid
R.grid <- function(r,proj,R)
{
  R <- R.prepare(R)
  PROJ <- R$PROJ
  X <- R$X
  Y <- R$Y
  Z <- R$Z
  dX <- R$dX
  dY <- R$dY
  dZ <- R$dZ
  R <- R$R

  DIM <- c(length(r$x),length(r$y))
  xy <- array(0,c(DIM,2))
  xy[,,1] <- r$x
  xy <- aperm(xy,c(2,1,3))
  xy[,,2] <- r$y
  xy <- aperm(xy,c(2,1,3))
  dim(xy) <- c(prod(DIM),2)

  # xy <- array(0,c(prod(DIM),2))
  # for(i in 1:DIM[1]) { for(j in 1:DIM[2]) { xy[i+(j-1)*DIM[1],] <- c(r$x[i],r$y[j]) } }
  colnames(xy) <- c('x','y')

  if(length(dim(R))==2)
  {
    G <- R.extract(xy,proj=proj,R=R,X=X,Y=Y,PROJ=PROJ,dX=dX,dY=dY)
    G <- array(G,DIM)
  }
  else if(length(dim(R))==3)
  {
    G <- array(0,c(DIM,length(Z)))
    for(i in 1:length(Z))
    { G[,,i] <- R.extract(xy,proj,R=R[,,i],X=X,Y=Y,PROJ=PROJ,dX=dX,dY=dY) }
  }

  return(G)
}

# evaluate habitat suitability raster(s)
R.suit <- function(R,CTMM,data=NULL)
{
  DIM <- dim(R[[1]])
  beta <- CTMM$beta

  offset <- stats::terms(CTMM$formula)
  offset <- attr(offset,"variables")[ attr(offset,"offset") ]

  PREP <- FALSE
  S <- 1

  if(length(beta))
  {
    BETA <- names(beta)

    PREP <- !all(BETA %in% names(R))
    if(!PREP)
    {
      # these will never be raster stacks
      S <- vapply(BETA,function(B){beta[B]*R[[B]]},R[[1]])
      dim(S) <- c(prod(DIM),length(beta))
    }
    else # formula required
    {
      R <- lapply(R,c)
      R <- data.frame(R)
      # need to copy over data if time varying formula
      if(!is.null(data)) { for(COL in names(data)) { R[[COL]] <- data[[COL]] } }

      n <- length(beta)

      # fix formula multiplication
      for(i in 1:n) { BETA[i] <- gsub(":","*",BETA[i]) }

      S <- sapply(1:n,function(i){beta[i]*eval(parse(text=BETA[i]),envir=R)})
    }
    S <- rowSums(S)
    S <- array(S,DIM)
    S <- exp(S)
  } # end beta

  if(length(offset))
  {
    if(all(offset %in% names(R)))
    { for(off in offset) { S <- S * R[[off]] } }
    else # formula required
    {
      # didn't prepare data before
      if(!PREP)
      {
        R <- lapply(R,c)
        R <- data.frame(R)
        # need to copy over data if time varying formula
        if(!is.null(data)) { for(COL in names(data)) { R[[COL]] <- data[[COL]] } }
      }

      # fix formula multiplication
      for(i in 1:length(offset)) { offset[i] <- gsub(":","*",offset[i]) }

      O <- sapply(1:length(offset),function(i){eval(parse(text=offset[i]),envir=R)})
      O <- apply(O,1,prod)
      O <- array(O,DIM)

      S <- O*S
    }
  } # end offset

  return(S)
}


# is RSF model stationary or non-stationary
is.stationary <- function(R,CTMM)
{
  STATIONARY <- TRUE

  DIM <- sapply(R,function(r){length(dim(r))})
  DIM <- max(DIM)
  if(DIM==3) { STATIONARY <- FALSE }

  formula <- CTMM$formula
  if(!is.null(formula))
  {
    VARS <- all.vars(formula)
    DVARS <- VARS[ VARS %nin% names(R) ]
    if(length(DVARS)) { STATIONARY <- FALSE }
  }

  return(STATIONARY)
}


# UNFINISHED
# this is for cross validation only
rsf.loglike <- function(data,CTMM,R=list(),smooth=TRUE,...)
{
  isotropic <- CTMM$isotropic
  beta <- CTMM$beta
  Z <- CTMM$Z
  integrated <- CTMM$integrated
  level.UD <- CTMM$level.UD
  formula <- CTMM$formula

  n <- nrow(data)
  UD <- list(CTMM=CTMM,weights=rep(1,n),DOF.area=n)

#  rsf.mcint(data,UD,R=R,formula=formula,integrated=integrated,level.UD=level.UD,isotropic=isotropic,beta=beta,smooth=smooth,NORM=Z,...)
}

# model selection on anisotropy

# model selection on phenomenological spatial parameters

# model selection on covariates
