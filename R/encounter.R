encounter <- function(object,include=NULL,exclude=NULL,...)
{
  UD <- object
  check.projections(UD)

  # by default, we exclude twin encounters
  if(is.null(include) && is.null(exclude)) { exclude <- diag(1,length(UD)) }
  if(is.null(include)) { include <- 1 - exclude }

  axes <- UD[[1]]$axes
  AXES <- length(axes)

  # Gaussian approximation
  CTMM <- lapply(UD,function(U){U@CTMM})
  CTMM <- encounter.ctmm(CTMM,include=include,exclude=exclude)
  DOF.area <- rep(DOF.area(CTMM),AXES)
  names(DOF.area) <- axes

  DIM <- dim(UD[[1]]$PDF)
  PDF <- matrix(0,DIM[1],DIM[2])
  for(i in 1:(length(UD)-1))
  {
    for(j in (i+1):length(UD))
    { PDF <- PDF + include[i,j] * UD[[i]]$PDF * UD[[j]]$PDF }
  }

  dV <- prod(UD[[1]]$dr)
  GAMMA <- sum(PDF) * dV # useful for relative encounter rates!
  PDF <- PDF / GAMMA

  object <- list()
  object$PDF <- PDF
  object$CDF <- pmf2cdf( PDF * dV )
  object$axes <- axes
  object$r <- UD[[1]]$r
  object$dr <- UD[[1]]$dr

  # resolution (add up like covariance?)
  IN <- 0
  H <- matrix(0,AXES,AXES)
  for(i in 1:(length(UD)-1))
  {
    for(j in (i+1):length(UD))
    {
      IN <- IN + include[i,j]
      H <- H + include[i,j] * PDsolve( PDsolve(UD[[i]]$H) + PDsolve(UD[[i]]$H) )
    }
  }
  H <- H/IN
  dimnames(H) <- list(axes,axes)

  object$H <- H
  object$DOF.area <- DOF.area

  info <- mean.info(UD)
  object <- new.UD(object,info=info,type='encounter',CTMM=CTMM)

  return(object)
}

################
# Gaussian encounters
encounter.ctmm <- function(CTMM,include=NULL,exclude=NULL,...)
{
  #check.projections(CTMM)

  # Gaussian / cumulants
  axes <- CTMM[[1]]$axes
  AXES <- length(axes)
  isotropic <- all(sapply(CTMM,function(C){C$isotropic}))

  fn <- function(CTMM)
  {
    IN <- 0
    M1 <- array(0,AXES)
    M2 <- matrix(0,AXES,AXES)

    for(i in 1:(length(CTMM)-1))
    {
      for(j in (i+1):length(CTMM))
      {
        Pi <- PDsolve(CTMM[[i]]$sigma)
        Pj <- PDsolve(CTMM[[j]]$sigma)
        sigma <- PDsolve(Pi + Pj)
        mu <- sigma %*% (Pi %*% CTMM[[i]]$mu[1,] + Pj %*% CTMM[[j]]$mu[1,])
        mu <- c(mu)

        IN <- IN + include[i,j]
        M1 <- M1 + include[i,j] * mu
        M2 <- M2 + include[i,j] * (sigma + outer(mu))
      }
    }

    M1 <- M1/IN
    M2 <- M2/IN

    mu <- M1
    sigma <- M2 - outer(M1)
    sigma <- covm(sigma)@par
    if(isotropic) { sigma <- sigma[1] }

    return(c(mu,sigma))
  }

  STUFF <- gauss.comp(fn,CTMM,COV=TRUE)

  I <- 1:AXES
  mu <- STUFF$MLE[I]
  names(mu) <- axes
  COV.mu <- STUFF$COV[I,I]
  dimnames(COV.mu) <- list(axes,axes)

  I <- (AXES+1):length(STUFF$MLE)
  sigma <- STUFF$MLE[I]
  sigma <- covm(sigma,axes=axes,isotropic=isotropic)
  COV <- STUFF$COV[I,I]
  NAMES <- names(sigma@par)[1:length(I)] # isotropic
  dimnames(COV) <- list(NAMES,NAMES)

  info <- mean.info(CTMM)

  CTMM <- ctmm(mu=mu,sigma=sigma,COV.mu=COV.mu,COV=COV,axes=axes,isotropic=isotropic,info=info)
  return(CTMM)
}
