# interpolate vector by continuous index (vectorized by index)
# vec is a vector, ind is a continuous index
vint <- function(vec,ind,return.ind=FALSE)
{
  n <- length(vec)

  lo <- floor(ind)
  hi <- ceiling(ind)

  # extrapolate
  lo <- pmax(1,lo)
  hi <- pmax(2,hi)
  # extrapolate
  lo <- pmin(n-1,lo)
  hi <- pmin(n,hi)

  if(return.ind) { return(rbind(lo,hi)) }

  # IND <- abs(vec[c(lo,hi)])
  # if(any(IND==Inf))
  # {
  #   IND <- which.max(IND)
  #   return(vec[c(lo,hi)[IND]])
  # }

  # linear interpolation
  vec <- vec[lo] + (vec[hi]-vec[lo])*(ind-lo)

  # fix NaNs (should be +/- Inf)
  INF <- abs(lo)==Inf
  if(any(INF)) { vec[INF] <- lo[INF] }
  INF <- abs(hi)==Inf
  if(any(INF)) { vec[INF] <- hi[INF] }

  return(vec)
}


# same thing as above but with a block-vector mat
mint <- function(mat,ind)
{
  IND <- vint(mat[1,],ind,return.ind=TRUE)
  mat <- mat[,IND[1,]] + (mat[,IND[2,]]-mat[,IND[1,]])*(ind-IND[1,])
  return(mat)
}


# bi-linear interpolation
bint <- function(M,ind)
{
  ind <- cbind(ind) # vectorize
  # rownames(ind) <- c('x','y')

  # index each dimension
  INDx <- vint(M[,1],ind[1,],return.ind=TRUE)
  INDy <- vint(M[1,],ind[2,],return.ind=TRUE)

  BINT <- function(i)
  {
    dX <- c(INDx[2,i]-ind[1,i],ind[1,i]-INDx[1,i]) / (INDx[2,i]-INDx[1,i])
    dX <- nant(dX,1) # in case of 0/0
    dX <- dX/sum(dX) # in case of two NaNs

    dY <- c(INDy[2,i]-ind[2,i],ind[2,i]-INDy[1,i]) / (INDy[2,i]-INDy[1,i])
    dY <- nant(dY,1)
    dY <- dY/sum(dY)

    c( dX %*% M[INDx[,i],INDy[,i]] %*% dY )
  }

  M <- vapply(1:ncol(ind),BINT,0)
  return(M)
}


# tri-linear interpolation
tint <- function(M,ind)
{
  ind <- cbind(ind) # vectorize
  # rownames(ind) <- c('x','y','z')

  # index each dimension
  INDx <- vint(M[,1,1],ind[1,],return.ind=TRUE)
  INDy <- vint(M[1,,1],ind[2,],return.ind=TRUE)
  INDz <- vint(M[1,1,],ind[3,],return.ind=TRUE)

  CINT <- function(i)
  {
    dX <- c(INDx[2,i]-ind[1,i],ind[1,i]-INDx[1,i]) / (INDx[2,i]-INDx[1,i])
    dX <- nant(dX,1)
    dX <- dX/sum(dX) # in case of two NaNs

    dY <- c(INDy[2,i]-ind[2,i],ind[2,i]-INDy[1,i]) / (INDy[2,i]-INDy[1,i])
    dY <- nant(dY,1)
    dY <- dY/sum(dY) # in case of two NaNs

    dZ <- c(INDz[2,i]-ind[3,i],ind[3,i]-INDz[1,i]) / (INDz[2,i]-INDz[1,i])
    dZ <- nant(dZ,1)
    dZ <- dZ/sum(dZ) # in case of two NaNs

    M <- M[INDx[,i],INDy[,i],INDz[,i]] # tensor block
    M <- dX %.% M # tensor contraction of first index
    M <- dY %*% M %*% dZ
    c(M)
  }

  M <- vapply(1:ncol(ind),CINT,0)
  return(M)
}
