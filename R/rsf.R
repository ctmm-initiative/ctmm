rsf.fit <- function(data,UD,R,error=0.001,...)
{
  # z <- -qnorm(error)
  # EXT <- extent(UD,level=1-error,complete=TRUE)

  FLAT <- projection(data)
  CTMM <- UD@CTMM
  # extract weights
  W <- UD$weights * UD$DOF.area

  R <- listify(R)

  # I should probably stick with raster format, but I find working with this a huge pain compared to arrays
  PROJ <- X <- Y <- dX <- dY <- list()
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

    # subset R to relevant area ?
  }

  N <- nrow(data)
  IID <- CTMM
  IID$tau <- NULL
  SIM <- simulate(CTMM,t=1:N,complete=TRUE)

  RDAT <- RSIM <- list()
  for(i in 1:length(R))
  {
    xy <- get.telemetry(data,c('longitude','latitude'))
    xy <- project(xy,to=PROJ[[i]])
    # convert to indices
    xy[,1] <- (xy[,1] - X[[i]][1])/dX[[i]] + 1
    xy[,2] <- (xy[,2] - Y[[i]][1])/dY[[i]] + 1
    #
    # VECTORIZE bint()
  }
  #
  #

  # record truncation error
  #
  #
}
