rsf.fit <- function(data,UD,R,beta=NULL,synoptic=TRUE,isotropic=FALSE,error=0.01,trace=FALSE,...)
{
  CTMM <- UD@CTMM
  IID <- CTMM
  IID$tau <- NULL # for simulation - uncorrelated

  # extract weights
  W <- mean(UD$DOF.area) # sum(w)
  w <- UD$weights * W

  R <- listify(R)
  NRES <- length(R)

  if(is.null(names(R)))
  {
    message("R is not a named list of covariates.")
    names(R) <- paste0("R-",1:length(R))
  }
  RNAMES <- names(R)

  # number of spatial covariates
  if(isotropic)
  { NSPA <- 3 }
  else
  { NSPA <- 5 }

  # get raster data
  # I should probably stick with raster format, but I find working raster objects a pain compared to arrays
  PROJ <- X <- Y <- dX <- dY <- list()
  RSCALE <- rep(1,length(R))
  for(i in 1:length(R))
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

    # subset R to relevant area ?
    #
    #
  }

  # data and simulation sampled covariates
  RDAT <- matrix(0,nrow(data),length(R)+NSPA)

  NAMES <- c(names(R),'x','y')
  if(isotropic) { NAMES <- c(NAMES,'rr') }
  else { NAMES <- c(NAMES,'xx','yy','xy') }
  colnames(RDAT) <- NAMES

  # for x-y scaling
  mu <- c(CTMM$mu)
  std <- sqrt(diag(CTMM$sigma))
  names(std) <- names(mu) <- c('x','y')

  SCALE <- c(RSCALE,1/std)
  if(isotropic) { SCALE <- c(SCALE,1/prod(std)) }
  else { SCALE <- c(SCALE,1/std[1]^2,1/std[2]^2,1/prod(std)) }
  names(SCALE) <- NAMES

  # initial estimates
  if(is.null(beta)) { beta <- rep(0,ncol(RDAT)) }
  names(beta) <- NAMES
  beta <- beta * SCALE

  COR <- stats::cov2cor(CTMM$sigma@.Data)
  iCOR <- PDsolve(COR)
  detCOR <- det(COR)

  # store observed covariates
  for(i in 1:length(R))
  {
    # store data-sampled covariates
    xy <- get.telemetry(data,c('longitude','latitude'))
    xy <- project(xy,to=PROJ[[i]])
    # convert to indices
    xy[,1] <- (xy[,1] - X[[i]][1])/dX[[i]] + 1
    xy[,2] <- (xy[,2] - Y[[i]][1])/dY[[i]] + 1

    RDAT[,i] <- bint(R[[i]],t(xy))
  }

  # store sampled synoptic terms
  # mean terms
  RDAT[,'x'] <- (data$x - mu['x'])/std['x']
  RDAT[,'y'] <- (data$y - mu['y'])/std['y']

  # variance/covariance terms (relative to pilot estimate)
  if(isotropic)
  {
    RDAT[,'rr'] <- RDAT[,'x']^2 + RDAT[,'y']^2
  }
  else
  {
    RDAT[,'xx'] <- RDAT[,'x']^2
    RDAT[,'yy'] <- RDAT[,'y']^2
    RDAT[,'xy'] <- RDAT[,'x'] * RDAT[,'y']
  }

  RAVE <- w %*% RDAT
  names(RAVE) <- colnames(RDAT)
  # subtract off pilot estimates after averaging (less computation)
  if(isotropic)
  {
    RAVE['rr'] <- RAVE['rr'] - W*2
  }
  else
  {
    RAVE['xx'] <- RAVE['xx'] - W*iCOR[1,1]
    RAVE['yy'] <- RAVE['yy'] - W*iCOR[2,2]
    RAVE['xy'] <- RAVE['xy'] - W*iCOR[1,2]
  }
  RSIM <- NULL

  # log-likelihood sans constant terms
  nloglike <- function(beta,zero=0)
  { c( -(RAVE %*% beta) + W*(log(mean(exp(RSIM %*% beta)))-zero/W) ) }
  # not sure this zero-ing is useful here

  N <- ceiling(W)^2 # starting value, will increase iteratively
  N.OLD <- 0
  STDloglike <- Inf # iterate until we meet our error tolerance
  while(STDloglike > error)
  {
    # double the random sample
    SIM <- simulate(IID,t=1:N,complete=TRUE)
    RSIM <- rbind(RSIM,matrix(0,N,length(R)+NSPA))
    colnames(RSIM) <- colnames(RDAT)

    OUT <- 0
    for(i in 1:length(R))
    {
      # store simulation-sampled covariates
      xy <- get.telemetry(SIM,c('longitude','latitude'))
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

      SUB <- N.OLD+1:nrow(SIM)
      RSIM[SUB,i] <- bint(R[[i]],t(xy))
    } # end for(i in 1:length(R))
    if(OUT) { warning(OUT/N,"% of random samples fell outside of rasters and were discarded.") }

    RSIM[SUB,'x'] <- (SIM$x - mu['x'])/std['x']
    RSIM[SUB,'y'] <- (SIM$y - mu['y'])/std['y']

    # variance/covariance terms
    if(isotropic)
    {
      RSIM[SUB,'rr'] <- RSIM[SUB,'x']^2 + RSIM[SUB,'y']^2 - 2
    }
    else
    {
      RSIM[SUB,'xx'] <- RSIM[SUB,'x']^2 - iCOR[1,1]
      RSIM[SUB,'yy'] <- RSIM[SUB,'y']^2 - iCOR[2,2]
      RSIM[SUB,'xy'] <- RSIM[SUB,'x']*RSIM[SUB,'y'] - iCOR[1,2]
    }

    # update SIM count
    N <- nrow(RSIM)

    # parscales will default to 1 if beta are 0
    if(trace) { message("Maximizing likelihood(N=",N,").") }
    RESULT <- optimizer(beta,nloglike,...)
    beta <- RESULT$par
    loglike <- -RESULT$value

    EXP <- exp(RSIM %*% beta)
    Elambda <- mean(EXP)
    Vlambda <- stats::var(EXP)
    STDloglike <- sqrt( W^2/N * Vlambda/Elambda^2 )
    if(trace) { message("STD[likelihood(N=",N,")] = ",STDloglike) }

    # numbers for next iteration
    N.OLD <- N
    N <- 2*N
  }

  # compute hessian
  if(trace) { message("Calculating Hessian.") }
  DIFF <- genD(par=beta,fn=nloglike,zero=-loglike,Richardson=2,mc.cores=1)
  hess <- DIFF$hessian
  grad <- DIFF$gradient
  # more robust covariance calculation than straight inverse
  COV <- cov.loglike(hess,grad)
  dimnames(COV) <- list(colnames(RDAT),colnames(RDAT))

  # log-normal terms from importance sampling
  SUMNORM <-  vapply(1:N,function(i){RDAT[i,c('x','y')] %*% iCOR %*% RDAT[i,c('x','y')]},1) # [n]
  SUMNORM <- -W/2*log(detCOR) - 1/2*c(w %*% SUMNORM)
  loglike <- loglike + SUMNORM

  # unscale beta & COV
  beta <- beta / SCALE
  COV <- COV / SCALE
  COV <- COV / SCALE

  # log-likelihood adjustment (unscale)
  loglike <- loglike - W*log(prod(std))

  # convert spatial parameters to adjusted mean and covariance estimates
  #
  #
  #

  # package results and return
  R <- list(beta=beta,COV=COV,loglike=loglike,VAR.loglike=STDloglike^2)
  return(R)
}

# model selection on anisotropy

# model selection on phenomenological spatial parameters

# model selection on covariates
