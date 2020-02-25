# 2D rotation matrix
rotate <- function(theta)
{
  COS <- cos(theta)
  SIN <- sin(theta)
  R <- rbind( c(COS,-SIN), c(SIN,COS) )
  return(R)
}


rotate.vec <- function(z,theta)
{
  R <- rotate(theta)
  z <- z %*% t(R)
  return(z)
}


rotate.mat <- function(M,theta)
{
  R <- rotate(theta)
  tR <- t(R)
  n <- dim(M)[1]
  M <- vapply(1:n,function(i){R %*% M[i,,] %*% tR},diag(2)) # (2,2,n)
  M <- aperm(M,c(3,1,2))
  return(M)
}


# rotation matrices array for multiple angles
rotates <- function(theta)
{
  R <- vapply(1:length(theta),function(i){rotate(theta[i])},diag(2)) # (2,2,n)
  R <- aperm(R,c(3,1,2)) # (n,2,2)
  return(R)
}


# apply rotation matrices R to vectors z
rotates.vec <- function(z,R)
{
  DIM <- dim(z) # (n,2*M) or (n,2,M)
  dim(z) <- c(DIM[1],2,prod(DIM[-1])/2) # (n,2,M)
  z <- vapply(1:nrow(z),function(i){R[i,,] %*% z[i,,]},array(0,dim(z)[-1])) # (2,M,n)
  z <- aperm(z,c(3,1,2)) # (n,2,M)
  dim(z) <- DIM
  return(z)
}


# apply rotation matrices R to matrices M
rotates.mat <- function(M,R) # could add tR=t(R) argument for small cost savings
{
  M <- vapply(1:dim(M)[1],function(i){R[i,,] %*% M[i,,] %*% t(R[i,,])},diag(2)) # (2,2,n)
  M <- aperm(M,c(3,1,2)) # (n,2,2)
  return(M)
}

# eccentricity squeeze transformation
squeeze <- function(z,smgm)
{
  z[,1] <- z[,1] / smgm
  z[,2] <- z[,2] * smgm
  return(z)
}

# eccentricity transform matrices
squeeze.mat <- function(M,smgm) # (n,2,2)
{
  R <- c(1/smgm,smgm)
  M <- aperm(M,c(2,3,1)) # (2,2,n)
  M <- R * M # (2',2,n)
  M <- aperm(M,c(2,3,1)) # (2,n,2')
  M <- R * M # (2',n,2')
  M <- aperm(M,c(2,3,1)) # (n,2',2')
  return(M)
}


##### det shouldn't fail because R dropped indices
det.numeric <- function(x,...) { x }
determinant.numeric <- function(x,logarithm=TRUE,...)
{
  SIGN <- sign(x)
  if(logarithm)
  { x <- log(abs(x)) }

  RESULT <- list(modulus=x,sign=SIGN)
  attr(RESULT$modulus,"logarithm") <- logarithm

  class(RESULT) <- "det"

  return(det)
}


# dyadic product default
outer <- function(X,Y=X,FUN="*",...) { base::outer(X,Y,FUN=FUN,...) }


# riffle columns of two matrices
riffle <- function(u,v,by=1)
{
  DIM <- dim(u) # (row,col)
  SUB <- 0:(by-1)
  u <- vapply(seq(1,DIM[2],by),function(i){cbind(u[,i+SUB],v[,i+SUB])},array(0,c(DIM[1],2*by))) # (row,2*by,col/by)
  dim(u) <- c(DIM[1],2*DIM[2])
  return(u)
}


# adjoint of matrix
Adj <- function(M) { t(Conj(M)) }

# Hermitian part of matrix
He <- function(M) { (M + Adj(M))/2 }


# map function for real-valued PSD matrices
PDfunc <-function(M,func=function(m){1/m},force=FALSE,pseudo=FALSE,tol=.Machine$double.eps)
{
  DIM <- dim(M)[1]
  if(is.null(DIM)) { DIM <- 1 }
  # tol <- max(tol,.Machine$double.eps)

  if(DIM==1)
  { M <- c(M) }
  else if(DIM==2)
  {
    TR <- (M[1,1] + M[2,2])/2 # half trace
    BIGNUM <- TR^2 > .Machine$double.xmax * .Machine$double.eps
    if(BIGNUM)
    {
      DET <- (M[1,1]/TR)*(M[2,2]/TR) - (M[1,2]/TR)*(M[2,1]/TR)
      DET <- 1 - DET
    }
    else
    {
      DET <- M[1,1]*M[2,2] - M[1,2]*M[2,1]
      DET <- TR^2 - DET
    }
    if(DET<0) # this shouldn't ever happen with Hermitian matrices
    {
      if(!force && !pseudo) { stop("Matrix not positive definite.") }
      else { DET <- 0 }
    }
    else if(M[1,2]==0) # DET==0 ++
    {
      M <- diag(M)
      V <- array(0,c(2,2,2))
      V[1,1,1] <- V[2,2,2] <- 1
    }
    else
    {
      DET <- sqrt(DET) # now root term
      if(BIGNUM) { DET <- DET * TR }

      # hermitian formula
      V <- diag(1/2,2) %o% c(1,1) + ((M-TR*diag(2))/(DET*2)) %o% c(1,-1)

      M <- TR + c(1,-1)*DET
    }
  }
  else if(DIM>2) # arbitrary DIM
  {
    M <- eigen(M)
    V <- M$vectors
    M <- Re(M$values)

    V <- vapply(1:DIM,function(i){Re(V[,i] %o% Conj(V[,i]))},diag(1,DIM))
  }

  if(any(M<0) && !force && !pseudo) { stop("Matrix not positive definite.") }

  # negative eigenvalues indicate size of numerical error
  # MIN <- last(M)
  # if(MIN<0) { tol <- max(tol,2*abs(MIN)) }

  FORCE <- (M < tol) -> PSEUDO
  # PSEUDO <- (abs(M) < tol) # why abs(M)?

  if(any(FORCE) && force) { M[FORCE] <- tol }
  M <- func(M)
  if(any(PSEUDO) && pseudo) { M[PSEUDO] <- 0 }

  if(DIM==1)
  { M <- cbind(M) }
  else
  {
    # add up from smallest contribution to largest contribution
    INDEX <- sort(abs(M),method="quick",index.return=TRUE)$ix
    M <- lapply(INDEX,function(i){nant(M[i]*V[,,i],0)}) # 0/0 -> 0
    M <- Reduce('+',M)
  }

  return(M)
}


# Positive definite solver
PDsolve <- function(M,force=FALSE,pseudo=FALSE,tol=.Machine$double.eps)
{
  DIM <- dim(M)
  if(is.null(DIM))
  {
    M <- as.matrix(M)
    DIM <- dim(M)
  }
  if(DIM[1]==0) { return(M) }

  # check for Inf & invert those to 0 (and vice versa)
  INF <- diag(M)==Inf
  ZERO <- sapply(1:nrow(M),function(i){all(M[i,]==0)})
  if(any(INF) || any(ZERO))
  {
    # 1/Inf == 0
    if(any(INF))
    { M[INF,INF] <- 0 }

    # 1/0 == Inf
    if(any(ZERO))
    {
      M[ZERO,ZERO] <- 0
      diag(M)[ZERO] <- Inf
    }

    # regular inverse of remaining dimensions
    REM <- !(INF|ZERO)
    if(any(REM)) { M[REM,REM] <- PDsolve(M[REM,REM,drop=FALSE]) }

    return(M)
  }

  if(!force && !pseudo)
  {
    if(DIM[1]==1) { return(matrix(1/M,c(1,1))) }
    if(DIM[1]==2)
    {
      DET <- M[1,1]*M[2,2]-M[1,2]*M[2,1]
      if(DET<=0) { return(diag(Inf,2)) } # force positive definite / diagonal
      SWP <- M[1,1] ; M[1,1] <- M[2,2] ; M[2,2] <- SWP
      M[1,2] <- -M[1,2]
      M[2,1] <- -M[2,1]
      return( M/DET )
    }
  }

  # symmetrize
  M <- He(M)

  # rescale
  W <- abs(diag(M))
  W <- sqrt(W)
  W <- W %o% W

  # now a correlation matrix that is easier to invert
  M <- M/W

  # try ordinary inverse
  M.try <- try(qr.solve(M,tol=tol),silent=TRUE)
  # fall back on decomposition
  if( class(M.try)[1] == "matrix")
  { M <- M.try }
  else
  { M <- PDfunc(M,func=function(m){1/m},force=force,pseudo=pseudo,tol=tol) }

  # back to covariance matrix
  M <- M/W

  # symmetrize
  M <- He(M)

  return(M)
}


# sqrtm fails on 1x1 matrix
# it also gives annoying notes if I don't cast it right
# also, my matrices are PSD
sqrtm <- function(M,force=FALSE,pseudo=FALSE)
{
  if(class(M)[1]=="covm")
  {
    par <- attr(M,"par")
    par['area'] <- sqrt(par['area'])
    M <- covm(par,isotropic=attr(M,'isotropic'),axes=dimnames(M)[[1]])
    return(M)
  }

  DIM <- dim(M)[1]
  if(is.null(DIM)) { DIM <- 1 }
  TOL <- DIM[1]*.Machine$double.eps

  if(DIM==1)
  {
    if(M>=0)
    { M <- sqrt(M) }
    else
    {
      if(force || pseudo) { M <- 0 } # round off error
      else { stop("Matrix is not positive definite.") } # complex sqrt
    }
  }
  else if(DIM==2)
  {
    TR <- M[1,1] + M[2,2]
    DET <- M[1,1]*M[2,2] - M[1,2]*M[2,1]

    if(DET<0 || TR^2<4*DET || diag(M)<0)
    {
      if(force || pseudo)
      { M <- PDfunc(M,func=sqrt,force=force,pseudo=pseudo) }
      else
      { stop("Matrix is not positive definite.") }
    }
    else
    {
      S <- sqrt(DET)
      M <- (M + S*diag(2))/sqrt(TR+2*S)
    }
  }
  else
  {
    if(all(diag(M)>=-TOL))
    {
      # R <- Matrix::Matrix(M,sparse=FALSE,doDiag=FALSE) # complains still.... ???
      # R <- as(R,"dpoMatrix")
      R <- expm::sqrtm(M) # reduces back to class "matrix" ?
    }
    else
    { R <- diag(-1,nrow=DIM) }

    if(all(Re(diag(R))>=-TOL && abs(Im(diag(R)))<=TOL))
    {
      M <- Re(R)
      TEST <- (diag(M)<=0)
      if(any(TEST)) { diag(M)[TEST] <- 0 }
    }
    else
    {
      if(force || pseudo)
      { M <- PDfunc(M,func=sqrt,force=force,pseudo=pseudo) }
      else
      { stop("Matrix is not positive definite.") }
    }
  }

  return(M)
}


# condition number
conditionNumber <- function(M)
{
  M <- try(eigen(M,only.values=TRUE)$values)
  if(class(M)[1]=="numeric")
  {
    M <- last(M)/M[1]
    M <- nant(M,Inf) # worst case
    return(M)
  }
  else
  { return(Inf) } # worst case
}


# Positive definite part of matrix
PDclamp <- function(M,lower=0,upper=Inf)
{
  # Inf fix
  INF <- diag(M)==Inf
  if(any(INF))
  {
    if(any(INF))
    {
      M[INF,INF] <- 0
      diag(M)[INF] <- upper # Inf
    }

    # regular clamp of remaining dimensions
    REM <- !INF
    if(any(REM)) { M[REM,REM] <- PDclamp(M[REM,REM,drop=FALSE],lower=lower,upper=upper) }
  }
  else
  {
    # symmetrize
    M <- He(M)

    M <- PDfunc(M,function(m){clamp(m,lower,upper)},pseudo=TRUE)

    # symmetrize
    M <- He(M)
  }

  return(M)
}


# relatively smallest eigen value of matrix
mat.min <- function(M)
{
  if(any(is.na(M)) || any(abs(M)==Inf) || any(diag(M)==0)) { return(0) }
  diag(M) <- abs(diag(M))
  M <- stats::cov2cor(M)
  M <- eigen(M,only.values=TRUE)$values
  M <- last(M)
  return(M)
}


# smallest maximum of matrices (for inner products on normalized vectors)
# if MAX=FALSE, largest minimum of matrices
ext.mat <- function(...,MAX=TRUE)
{
  MATS <- list(...)
  DIM <- dim(MATS[[1]])[1]
  MATS <- lapply(MATS,eigen)

  VAL <- NULL
  VEC <- NULL
  for(i in 1:length(MATS))
  {
    VAL <- c(VAL,MATS[[i]]$values)
    VEC <- cbind(VEC,MATS[[i]]$vectors)
  }
  rm(MATS)

  ORDER <- order(VAL,decreasing=MAX)
  VAL <- VAL[ORDER]
  VEC <- VEC[,ORDER]

  MATS <- list()
  for(i in 1:DIM)
  {
    MATS[[i]] <- VAL[i] * (VEC[,i] %o% VEC[,i])
    # Grahm-Schmidt orthogonalization
    if(i<DIM) { for(j in (i+1):DIM) { VEC[,j] <- VEC[,j] - c(VEC[,i] %*% VEC[,j]) * VEC[,i] } }
  }

  MATS <- Reduce("+",MATS)
  return(MATS)
}

