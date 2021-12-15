rsf.fit <- function(data,UD,R=list(),mod=NULL,int=TRUE,iRSF=TRUE,level.UD=0.99,isotropic=TRUE,beta=NULL,debias=TRUE,smooth=TRUE,error=0.01,trace=TRUE,...)
{ rsf.mcint(data,UD,R=R,int=int,mod=mod,beta=beta,iRSF=iRSF,level.UD=level.UD,isotropic=isotropic,debias=debias,smooth=smooth,error=error,trace=trace,...) }

rsf.mcint <- function(data,UD,R=list(),mod=NULL,int=TRUE,iRSF=TRUE,level.UD=0.99,isotropic=TRUE,beta=NULL,debias=TRUE,smooth=TRUE,error=0.01,trace=TRUE,Z=NULL,...)
{
  # error <- error^2 # convert from error[beta] to error[loglike]

  CTMM <- UD@CTMM
  axes <- CTMM$axes

  # smooth the data
  if(smooth && any(CTMM$error>0))
  { data <- predict(data,CTMM=CTMM,t=data$t,complete=TRUE) }

  if(!iRSF) # prepare available region polygon
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
  NRES <- length(R)

  int <- rep(int,NRES)

  if(length(R) && is.null(names(R)))
  {
    message("R is not a named list of covariates.")
    names(R) <- 1:length(R)
  }
  RNAMES <- names(R)

  # number of spatial covariates
  if(!iRSF)
  { NSPA <- 0 }
  else if(isotropic)
  { NSPA <- 3 }
  else
  { NSPA <- 5 }

  if(!isotropic) { theta <- CTMM$sigma@par['angle'] }

  # get raster data
  # I should probably stick with raster format, but I find working raster objects a pain compared to arrays
  # will need to change this code to conserve memory
  PROJ <- X <- Y <- dX <- dY <- list()
  RSCALE <- rep(1,length(R))
  for(i in 1 %:% length(R))
  {
    PROJ[[i]] <- raster::projection(R[[i]])

    DIM <- dim(R[[i]])
    X[[i]] <- raster::xFromCol(R[[i]],1:DIM[2])
    Y[[i]] <- raster::yFromRow(R[[i]],1:DIM[1])

    dX[[i]] <- stats::median(diff(X[[i]]))
    dY[[i]] <- stats::median(diff(Y[[i]]))

    R[[i]] <- raster::as.array(R[[i]])[,,1] # not considering time yet
    R[[i]] <- aperm(R[[i]],2:1) # [x,y]

    R[[i]] <- R[[i]] - mean(R[[i]])

    RSCALE[i] <- sqrt(stats::var(c(R[[i]])))
    R[[i]] <- R[[i]]/RSCALE[i]
  }

  # data and simulation sampled covariates
  RDAT <- matrix(0,nrow(data),length(R)+NSPA)

  if(iRSF)
  {
    NAMES <- c(names(R),axes)
    if(isotropic) { NAMES <- c(NAMES,'rr') }
    else { NAMES <- c(NAMES,'xx','yy','xy') }
  }
  else
  { NAMES <- names(R) }
  colnames(RDAT) <- NAMES

  # for x-y scaling
  mu <- c(CTMM$mu)
  std <- sqrt(CTMM$sigma@par[c('major','minor')])
  names(std) <- names(mu) <- axes

  if(iRSF)
  {
    SCALE <- c(RSCALE,std)
    if(isotropic) { SCALE <- c(SCALE,prod(std)) }
    else { SCALE <- c(SCALE,std[1]^2,std[2]^2,prod(std)) }
  }
  else
  { SCALE <- RSCALE }
  names(SCALE) <- NAMES

  # initial estimates
  if(is.null(beta)) { beta <- rep(0,ncol(RDAT)) }
  names(beta) <- NAMES
  beta <- beta * SCALE

  # store observed covariates
  for(i in 1 %:% length(R))
  {
    # store data-sampled covariates
    xy <- get.telemetry(data,c('longitude','latitude'))
    xy <- project(xy,to=PROJ[[i]])
    # convert to indices
    xy[,1] <- (xy[,1] - X[[i]][1])/dX[[i]] + 1
    xy[,2] <- (xy[,2] - Y[[i]][1])/dY[[i]] + 1

    RDAT[,i] <- bint(R[[i]],t(xy))
  }

  ### store sampled integrated terms
  if(iRSF)
  {
    # de-trend mean
    RDAT[,'x'] <- data$x - mu['x']
    RDAT[,'y'] <- data$y - mu['y']

    # rotate to major-minor axes
    if(!isotropic) { RDAT[,axes] <- rotate.vec(RDAT[,axes],-theta) }

    # standardize
    RDAT[,'x'] <- RDAT[,'x']/std['x']
    RDAT[,'y'] <- RDAT[,'y']/std['y']

    # variance/covariance terms (relative to pilot estimate)
    if(isotropic) # beta is correction to standardized 1/sigma
    {
      RDAT[,'rr'] <- -( RDAT[,'x']^2 + RDAT[,'y']^2 )/2
    }
    else # beta is correction to standardized solve(sigma)
    {
      RDAT[,'xx'] <- -RDAT[,'x']^2 /2
      RDAT[,'yy'] <- -RDAT[,'y']^2 /2
      RDAT[,'xy'] <- -RDAT[,'x']*RDAT[,'y']
    }
  } # end if(iRSF)

  RAVE <- w %*% RDAT
  names(RAVE) <- colnames(RDAT)
  RSIM <- NULL

  # partition function defined (scaled) - calling from rsf.loglike
  if(!is.null(Z))
  {
    loglike <- (RAVE %*% beta) - W*Z # Z includes scaling
    # log-normal terms from importance sampling
    if(iRSF)
    {
      SUMNORM <- rowSums(RDAT[,axes]^2) # [n]
      SUMNORM <- - 1/2*c(w %*% SUMNORM) # + log(1) terms
      loglike <- loglike + SUMNORM
    }
    return(loglike)
  }

  # log-likelihood sans constant terms
  nloglike <- function(beta,zero=0)
  { c( -(RAVE %*% beta) + W*(log(mean(exp(RSIM %*% beta)))-zero/W) ) }
  # not sure this zero-ing is useful here

  N <- ceiling(W)^2 # starting value, will increase iteratively
  N.OLD <- 0
  loglike <- -Inf
  ERROR.BIG <- TRUE
  STDloglike <- Inf
  while(ERROR.BIG)
  {
    # double the random sample
    if(iRSF)
    { SIM <- simulate(IID,t=1:N,complete=TRUE) }
    else
    {
      SIM <- sp::spsample(level.UD,n=N,type="random")
      # sp::spsample drops projection information and throws annoying warning when trying to fix
      suppressWarnings( sp::proj4string(SIM) <- sp::CRS(projection(data)) )
      SIM <- sp::spTransform(SIM,sp::CRS(DATUM))
      SIM <- SIM@coords
      colnames(SIM) <- c('longitude','latitude')
    }
    RSIM <- rbind(RSIM,matrix(0,N,length(R)+NSPA))
    colnames(RSIM) <- NAMES

    OUT <- 0
    for(i in 1 %:% length(R))
    {
      # store simulation-sampled covariates
      if(iRSF)
      { xy <- get.telemetry(SIM,c('longitude','latitude')) }
      else
      { xy <- SIM }
      xy <- project(xy,to=PROJ[[i]])
      xy[,1] <- (xy[,1] - X[[i]][1])/dX[[i]] + 1
      xy[,2] <- (xy[,2] - Y[[i]][1])/dY[[i]] + 1

      # catch truncation error and record !!!
      BAD <- (xy[,1]<0) | (nrow(R[[i]])+1<xy[,1]) | (xy[,2]<0) | (ncol(R[[i]])+1<xy[,2])
      BADS <- sum(BAD)
      if(BADS)
      {
        OUT <- OUT + BADS
        BAD <- which(BAD)

        # remove bad samples from this set
        xy <- xy[-BAD,]

        # remove bad samples from simulation
        SIM <- SIM[-BAD,]

        # remove bad samples from all rasters
        RSIM <- RSIM[-(N.OLD+BAD),]
      }

      SUB <- N.OLD+1%:%nrow(SIM)
      if(length(SUB)) { RSIM[SUB,i] <- bint(R[[i]],t(xy)) }
    } # end for(i in 1:length(R))

    OUT <- OUT/(OUT+nrow(SIM))
    if(OUT > error) { warning(100*OUT,"% of random samples fell outside of rasters and were discarded.") }

    SUB <- N.OLD+1:nrow(SIM)

    if(iRSF)
    {
      # de-trend
      RSIM[SUB,'x'] <- SIM$x - mu['x']
      RSIM[SUB,'y'] <- SIM$y - mu['y']

      # rotate
      if(!isotropic)
      { RSIM[SUB,axes] <- rotate.vec(RSIM[SUB,axes],-theta) }

      # standardize
      RSIM[SUB,'x'] <- RSIM[SUB,'x']/std['x']
      RSIM[SUB,'y'] <- RSIM[SUB,'y']/std['y']

      # variance/covariance terms
      if(isotropic)
      {
        RSIM[SUB,'rr'] <- -( RSIM[SUB,'x']^2 + RSIM[SUB,'y']^2 )/2
      }
      else
      {
        RSIM[SUB,'xx'] <- -RSIM[SUB,'x']^2/2
        RSIM[SUB,'yy'] <- -RSIM[SUB,'y']^2/2
        RSIM[SUB,'xy'] <- -RSIM[SUB,'x']*RSIM[SUB,'y']
      }
    }

    # update SIM count
    N <- nrow(RSIM)

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

    EXP <- exp( c(RSIM %*% beta) )
    Elambda <- mean(EXP)
    Vlambda <- stats::var(EXP)
    STDloglike <- sqrt( W^2/N * Vlambda/Elambda^2 )
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
  dimnames(COV) <- list(colnames(RDAT),colnames(RDAT))

  if(iRSF)
  {
    # fixed log-normal terms from importance sampling
    SUMNORM <- rowSums(RDAT[,axes]^2) # [n]
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

  # partition function (per W)
  Z <- log(mean(exp(RSIM %*% beta)))

  # un-scale beta & COV
  beta <- beta / SCALE
  COV <- COV / outer(SCALE)

  if(iRSF)
  {
    # log-likelihood adjustment (unscale)
    loglike <- loglike - W*log(prod(std))
    # partition function adjustment (unscale)
    Z <- Z + log(prod(std)) # per W

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

  } # end if(iRSF)
  else
  {
    loglike <- loglike - W*log(AREA)
    Z <- Z + log(AREA)
  }

  if(length(RNAMES)) { RNAMES <- paste('RSF',RNAMES) }
  if(iRSF)
  {
    NAMES <- c(RNAMES,axes)
    if(isotropic)
    { NAMES <- c(NAMES,'major') }
    else
    { NAMES <- c(NAMES,'major','minor','angle')}
  }
  else
  { NAMES <- RNAMES }
  names(beta) <- NAMES
  dimnames(COV) <- list(NAMES,NAMES)

  # debias code (variance only)
  if(iRSF && debias)
  {
    if(!is.null(CTMM$MLE))
    { DEBIAS <- sqrt(det.covm(CTMM$sigma)/det.covm(CTMM$MLE$sigma)) }
    else
    { DEBIAS <- W/(W-1) }

    DSCALE <- rep(1,length(NAMES))
    names(DSCALE) <- NAMES
    DSCALE['major'] <- DEBIAS
    if(!isotropic) { DSCALE['minor'] <- DEBIAS }

    beta <- beta * DSCALE
    COV <- COV * outer(DSCALE)

    # impact on RSF variables missing
  }

  ## convert to ctmm() object format with extra beta slot & features
  if(iRSF)
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
  RSF$VAR.loglike <- STDloglike^2
  RSF$Z <- Z
  RSF$iRSF <- iRSF
  if(!iRSF) { RSF$level.UD <- level.UD }

  # copy over autocorrelation information in a reasonable way
  if(iRSF)
  {
    RSF$tau <- CTMM$tau
    RSF$omega <- CTMM$omega
    RSF$error <- CTMM$error

    VAR <- diag(CTMM$COV)
    COR <- stats::cov2cor(CTMM$COV)
    NEW <- names(VAR)[names(VAR) %nin% NAMES]
    OLD <- names(VAR)[names(VAR) %in% NAMES]
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
    RSF$features <- c(RNAMES,CTMM$features)
  }
  else
  {
    if(length(beta)) { dimnames(RSF$COV) <- list(RNAMES,RNAMES) }
    RSF$features <- RNAMES
    RSF$range <- FALSE
  }

  # check if some data fell outside of polygon
  if(!iRSF)
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


# this is for cross validation only
rsf.loglike <- function(data,CTMM,R=list(),mod=NULL,int=TRUE,smooth=TRUE,...)
{
  isotropic <- CTMM$isotropic
  beta <- CTMM$beta
  Z <- CTMM$Z
  iRSF <- CTMM$iRSF
  level.UD <- CTMM$level.UD

  n <- nrow(data)
  UD <- list(CTMM=CTMM,weights=rep(1,n),DOF.area=n)

  rsf.mcint(data,UD,R=R,mod=mod,int=int,iRSF=iRSF,level.UD=level.UD,isotropic=isotropic,beta=beta,smooth=smooth,Z=Z,...)
}

# model selection on anisotropy

# model selection on phenomenological spatial parameters

# model selection on covariates
