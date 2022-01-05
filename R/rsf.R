rsf.fit <- function(data,UD,beta=NULL,R=list(),formula=NULL,integrated=TRUE,level.UD=0.99,isotropic=TRUE,debias=TRUE,smooth=TRUE,error=0.01,trace=TRUE,...)
{ rsf.mcint(data,UD,beta=beta,R=R,formula=formula,integrated=integrated,level.UD=level.UD,isotropic=isotropic,debias=debias,smooth=smooth,error=error,trace=trace,...) }

rsf.mcint <- function(data,UD,beta=NULL,R=list(),formula=NULL,integrated=TRUE,level.UD=0.99,isotropic=TRUE,debias=TRUE,smooth=TRUE,error=0.01,trace=TRUE,NORM=NULL,...)
{
  # error <- error^2 # convert from error[beta] to error[loglike]
  STATIONARY <- TRUE
  CTMM <- UD@CTMM
  axes <- CTMM$axes
  GEO <- c('longitude','latitude')

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

  if(length(R) && is.null(names(R)))
  {
    warning("R is not a named list of rasters.")
    names(R) <- paste("R",1:length(R))
  }

  ## The basic setup is:
  # RVARS - raster variables
  # DVARS - data annotation variables - overridden by RVARS
  # TERMS - all predictors - some may be RVARS
  # CVARS - composite variables in TERMS that are not in RVARS || DVARS
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
        colnames(envir) <- TERMS
        RET <- evaluate(term,envir)
        dim(RET) <- DIM[1:2]
      }
      else # !STATIONARY
      {
        RET <- sapply(1:dim(envir)[1],function(i)
                {
                  ENVIR <- envir[i,,];
                  dim(ENVIR) <- dim(envir)[-1];
                  colnames(ENVIR) <- TERMS
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

    DIM <- dim(R[[i]])
    X[[i]] <- raster::xFromCol(R[[i]],1:DIM[2])
    Y[[i]] <- raster::yFromRow(R[[i]],1:DIM[1])
    Z[[i]] <- raster::getZ(R[[i]])

    dX[i] <- stats::median(diff(X[[i]]))
    dY[i] <- stats::median(diff(Y[[i]]))

    R[[i]] <- raster::as.array(R[[i]])[,,] # last dim will be dropped if length-1
    if(length(dim(R[[i]]))==2) # RasterLayer (static)
    {
      # R[[i]] <- R[[i]][,,1]
      R[[i]] <- aperm(R[[i]],2:1) # [x,y]
    }
    else if (length(dim(R[[i]]))==3) # RasterStack or RasterBrick
    {
      STATIONARY <- FALSE
      R[[i]] <- aperm(R[[i]],c(2,1,3)) # [x,y,z]

      if(class(Z[[i]])[1]=="Date") { Z[[i]] <- as.POSIXct(Z[[i]]) } # numeric is in days
      if(class(Z[[i]])[1]=="POSIXct") { Z[[i]] <- as.numeric(Z[[i]]) }

      dZ[i] <- stats::median(diff(Z[[i]]))

      # index for data times in raster stack
      Z.ind[[i]] <- (data$t - Z[[i]][1])/dZ[i] + 1
    }
    else
    { stop("") }

    # if(names(R)[i] %in% TERMS) { R[[i]] <- R[[i]] - mean(R[[i]]) } # could go wrong with effect modifiers

    # RSCALE[i] <- sqrt(stats::var(c(R[[i]])))
    # R[[i]] <- R[[i]]/RSCALE[i]
  }
  names(X) <- names(Y) <- names(dX) <- names(dY) <- RVARS
  if(length(Z)) { names(Z) <- names(dZ) <- names(Z.ind) <- RVARS }
  # names(RSCALE) <- RVARS

  ## prepare annotated data
  # DSCALE <- rep(1,length(DVARS))
  # for(i in 1%:% length(DVARS))
  # {
  #   DSCALE[i] <- sqrt(stats::var(data[[DVARS]]))
  #   data[[DVARS]] <- data[[DVARS]]/DSCALE[i]
  # }
  # names(DSCALE) <- DVARS

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
  for(r in RVARS)
  {
    # store data-sampled covariates
    xy <- get.telemetry(data,GEO)
    xy <- project(xy,to=PROJ[[i]])
    # convert to indices
    xy[,1] <- (xy[,1] - X[[r]][1])/dX[r] + 1
    xy[,2] <- (xy[,2] - Y[[r]][1])/dY[r] + 1

    if(length(dim(R[[r]]))==2)
    { DATA[,r] <- bint(R[[r]],t(xy)) }
    else # 3D
    {
      xy <- cbind(xy,Z.ind[[r]])
      DATA[,r] <- tint(R[[r]],t(xy))
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

    # if(!STATIONARY) { t.rand <- sample.time(N) } # random time indices according to weights

    # store simulation-sampled raster covariates
    # OUT <- 0
    for(r in RVARS)
    {
      if(integrated)
      { xy <- get.telemetry(SIM,GEO) }
      else
      { xy <- SIM }
      xy <- project(xy,to=PROJ[[i]])
      xy[,1] <- (xy[,1] - X[[r]][1])/dX[r] + 1
      xy[,2] <- (xy[,2] - Y[[r]][1])/dY[r] + 1

      # catch truncation error and record !!!
      BAD <- (xy[,1]<0) | (nrow(R[[r]])+1<xy[,1]) | (xy[,2]<0) | (ncol(R[[r]])+1<xy[,2])
      BADS <- sum(BAD)

      BADS <- BADS/(length(SIM))
      if(BADS > error) { warning(100*BADS,"% of random samples fell outside of rasters.") }

      # if(BADS)
      # {
      #   OUT <- OUT + BADS
      #   BAD <- which(BAD)
      #
      #   # remove bad samples from this set
      #   xy <- xy[-BAD,]
      #   if(!STATIONARY) { t.rand <- t.rand[-BAD] }
      #
      #   # remove bad samples from simulation
      #   SIM <- SIM[-BAD,]
      #
      #   # remove bad samples from all rasters
      #   SATA <- SATA[-(N.OLD+BAD),]
      # }

      # SUB <- N.OLD+1%:%nrow(SIM)

      if(length(SUB))
      {
        if(length(dim(R[[r]]))==2)
        { SATA[SUB,,r] <- bint(R[[r]],t(xy)) }
        else
        {
          xy <- cbind(xy,rep(data$t,N))
          SATA[SUB,,r] <- tint(R[[r]],t(xy))
        }
      }
    } # end for(i in 1:length(R))

    # OUT <- OUT/(OUT+nrow(SIM))
    # if(OUT > error) { warning(100*OUT,"% of random samples fell outside of rasters and were discarded.") }
    #
    # SUB <- N.OLD+1:nrow(SIM)

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
    # parscales will default to 1 if beta are 0
    control <- list(precision=precision)
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

  # package results and return
  # turn this into ctmm object
  RSF <- ctmm(axes=axes,mu=mu,COV.mu=COV.mu,beta=beta,sigma=sigma,isotropic=isotropic,COV=COV)
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

  rsf.mcint(data,UD,R=R,formula=formula,integrated=integrated,level.UD=level.UD,isotropic=isotropic,beta=beta,smooth=smooth,NORM=Z,...)
}

# model selection on anisotropy

# model selection on phenomenological spatial parameters

# model selection on covariates
