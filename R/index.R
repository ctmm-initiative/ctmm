# x is a UD object only (for now)
# R is a single raster between 0 and 1
index <- function(x,R,...)
{
  I <- index.UD(x,R,...)
  error <- x$error
  if(is.null(error)) { error <- 0.001 } # default numerical error in UD for old objects
  # SWITCH TO INTERPOLATION FOR SMOOTHNESS
  if(exp(-I^2) <= 2*error) # calculation is probably not good somewhere around this amount
  { I <- index.ctmm(x@CTMM,R,...) }
  return(I)
}

index.UD <- function(UD,R,...)
{
  R <- R.grid(UD$r,proj=projection(UD),R,...)
  R <- log(sum(R * UD$PDF)) + log(prod(UD$dr))

  R <- sqrt(max(-R,0))
  return(R)
}

index.ctmm <- function(CTMM,R,...)
{
  # raster coordinates
  PROJ <- raster::projection(R)
  DIM <- dim(R)[2:1]
  names(DIM) <- c('x','y')
  X <- raster::xFromCol(R,1:DIM['x'])
  Y <- raster::yFromRow(R,1:DIM['y'])

  # project raster to CTMM coordinates
  xy <- array(0,c(DIM,2)) # [x,y,2]
  dimnames(xy) <- list(NULL,NULL,c('x','y'))
  xy[,,'x'] <- X
  xy <- aperm(xy,c(2,1,3)) # [y,x,2]
  xy[,,'y'] <- Y
  xy <- aperm(xy,c(2,1,3)) # [x,y,2]

  dim(xy) <- c(prod(DIM),2) # [xy,2]
  colnames(xy) <- c('x','y')
  xy <- project(xy,from=PROJ,to=projection(CTMM))

  # cell areas
  dim(xy) <- c(DIM,2) # [x,y,2]
  dimnames(xy) <- list(NULL,NULL,c('x','y'))
  dx <- xy[-1,,'x'] - xy[-DIM['x'],,'x']
  dy <- xy[,-1,'y'] - xy[,-DIM['y'],'y']
  # average left and right intervals
  dX <- array(0,DIM)
  dY <- array(0,DIM)
  dX[-1,] <- dx/2
  dY[,-1] <- dy/2
  dX[-DIM['x'],] <- dX[-DIM['x'],] + dx/2
  dY[,-DIM['y']] <- dY[,-DIM['y']] + dy/2
  rm(dx,dy)
  dA <- abs(dX*dY)
  rm(dX,dY)
  dim(dA) <- prod(DIM)
  dim(xy) <- c(prod(DIM),2) # [xy,2]
  colnames(xy) <- c('x','y')

  # detrend mean
  xy[,'x'] <- xy[,'x'] - CTMM$mu[1,'x']
  xy[,'y'] <- xy[,'y'] - CTMM$mu[1,'y']

  # standardize
  COV <- CTMM$sigma
  isqCOV <- isqrtm.covm(COV)
  xy <- xy %*% isqCOV

  # quadratic loss
  xy <- xy[,'x']^2 + xy[,'y']^2
  xy <- (-1/2) * xy

  # underflow prevention
  R <- raster::as.array(R)
  dim(R) <- dim(R)[1:2]
  R <- t(R)
  MIN <- min( xy[R>0] )
  SHIFT <- max(MIN,0)
  xy <- xy - SHIFT
  xy <- exp(xy)

  dim(R) <- prod(DIM)
  R <- log(sum(R * xy * dA)) - (log(2*pi) + log(det.covm(COV))/2)
  R <- R + SHIFT # add back minimum

  R <- sqrt(max(-R,0))
  return(R)
}
