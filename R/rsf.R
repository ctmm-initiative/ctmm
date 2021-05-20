rsf.fit <- function(data,UD,R,beta=NULL,isotropic=FALSE,error=0.001,...)
{
  CTMM <- UD@CTMM
  # extract weights
  w <- UD$weights * UD$DOF.area
  W <- sum(w)

  R <- listify(R)
  NRES <- length(R)

  # number of spatial covariates
  if(isotropic)
  { NSPA <- 3 }
  else
  { NSPA <- 5 }

  # I should probably stick with raster format, but I find working raster objects a huge pain compared to arrays
  PROJ <- X <- Y <- dX <- dY <- list()
  SCALE <- rep(1,length(R))
  for(i in 1:length(R))
  {
    PROJ[[i]] <- raster::projection(R[[i]])

    DIM <- dim(R[[i]])
    X[[i]] <- raster::xFromCol(R[[i]],1:DIM[2])
    Y[[i]] <- raster::yFromRow(R[[i]],1:DIM[1])

    dX[[i]] <- stats::median(diff(X[[i]]))
    dY[[i]] <- stats::median(diff(Y[[i]]))

    R[[i]] <- as.array(R[[i]])[,,1] # not considering time yet
    R[[i]] <- aperm(R[[i]],2:1) # [x,y]

    SCALE[i] <- sqrt(stats::var(R[[i]]))
    R[[i]] <- R[[i]]/SCALE[i]

    # subset R to relevant area ?
  }

  N <- nrow(data)
  IID <- CTMM
  IID$tau <- NULL
  SIM <- simulate(CTMM,t=1:N,complete=TRUE) # will increase this number

  RDAT <- array(0,c(nrow(data),length(R)+NSPA))
  RSIM <- array(0,c(nrow(SIM),length(R)+NSPA))

  # initial estimates
  if(is.null(beta))
  { beta <- rep(0,ncol(RDAT)) }
  else
  { beta <- beta * SCALE }

  colnames(RDAT)[1:NRES] <- colnames(RSIM)[1:NRES] <- names(R)
  colnames(RDAT)[NRES + 1:2] + colnames(RSIM)[NRES + 1:2] <- c('x','y')
  if(isotropic)
  { colnames(RDAT)[NRES+3] <- colnames(RSIM)[NRES+3] <- 'rr' }
  else
  { colnames(RDAT)[NRES+3:5] <- colnames(RSIM)[NRES+3:5] <- c('xx','yy','xy') }

  # for x-y scaling
  mu <- c(CTMM$mu)
  if(isotropic)
  { std <- sqrt(mean(diag(CTMM$sigma))) }
  else
  { std <- sqrt(diag(CTMM$sigma)) }

  for(i in 1:length(R))
  {
    xy <- get.telemetry(data,c('longitude','latitude'))
    xy <- project(xy,to=PROJ[[i]])
    # convert to indices
    xy[,1] <- (xy[,1] - X[[i]][1])/dX[[i]] + 1
    xy[,2] <- (xy[,2] - Y[[i]][1])/dY[[i]] + 1

    RDAT[,i] <- mint(R[[i]],xy)

    xy <- get.telemetry(SIM,c('longitude','latitude'))
    xy <- project(xy,to=PROJ[[i]])
    xy[,1] <- (xy[,1] - X[[i]][1])/dX[[i]] + 1
    xy[,2] <- (xy[,2] - Y[[i]][1])/dY[[i]] + 1

    # catch truncation error and record !!!!!!!!!!!!!!!!!
    #
    #

    RSIM[,i] <- mint(R[[i]],xy)
  }

  # mean terms
  RDAT[,'x'] <- (data$x - mu[1])/std[1]
  RDAT[,'y'] <- (data$y - mu[2])/std[2]

  RSIM[,'x'] <- (SIM$x - mu[1])/std[1]
  RSIM[,'y'] <- (SIM$y - mu[2])/std[2]

  # variance/covariance terms
  if(isotropic)
  {
    RDAT[,'rr'] <- RDAT[,'x']^2 + RDAT[,'y']^2

    RSIM[,'rr'] <- RSIM[,'x']^2 + RSIM[,'y']^2
  }
  else
  {
    RDAT[,'xx'] <- RDAT[,'x']^2
    RDAT[,'yy'] <- RDAT[,'y']^2
    RDAT[,'xy'] <- RDAT[,'x'] * RDAT[,'y']

    RSIM[,'xx'] <- RSIM[,'x']^2
    RSIM[,'yy'] <- RSIM[,'y']^2
    RSIM[,'xy'] <- RSIM[,'x'] * RSIM[,'y']
  }

  RAVE <- w %*% RDAT

  # log-likelihood
  nloglike <- function(beta)
  { c( -(RAVE %*% beta) + W*log( mean(exp(RSIM %*% beta)) ) ) }

  # parscales will default to 1 if beta are 0
  RESULT <- optimizer(beta,nloglike,...)
  beta <- RESULT$par
  loglike <- -RESULT$value

  EXP <- exp(RSIM %*% beta)
  Elambda <- mean(EXP)
  Vlambda <- stats::var(EXP)
  STDloglike <- sqrt( W^2/N * Vlambda/Elambda^2 )

  # loop within error tolerance
}
