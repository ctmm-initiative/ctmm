encounter <- function(object,include=NULL,exclude=NULL,debias=FALSE,...)
{
  UD <- object
  check.projections(UD)

  # by default, we exclude twin encounters
  if(is.null(include) && is.null(exclude)) { exclude <- diag(1,length(UD)) }
  if(is.null(include)) { include <- 1 - exclude }

  axes <- UD[[1]]$axes
  AXES <- length(axes)

  CTMM <- lapply(UD,function(U){U@CTMM})

  # Gaussian approximation
  STUFF <- encounter.ctmm(CTMM,include=include,exclude=exclude,debias=debias)
  CTMM <- STUFF$CTMM
  BIAS <- STUFF$BIAS # individual biases
  bias <- STUFF$bias # pairwise biases

  DOF.area <- rep(DOF.area(CTMM),AXES)
  names(DOF.area) <- axes

  # bias correction (precision)
  for(i in 1:length(object))
  {
    UD[[i]]$PMF <- UD[[i]]$PDF * prod(UD[[i]]$dr)
    if(debias) { UD[[i]]$PMF <- debias.volume(UD[[i]]$PMF,1/BIAS[i])$PMF }
  }

  DIM <- dim(UD[[1]]$PDF)
  PMF <- matrix(0,DIM[1],DIM[2])
  for(i in 1:(length(UD)-1))
  {
    for(j in (i+1):length(UD))
    {
      PROD <- UD[[i]]$PMF * UD[[j]]$PMF
      if(debias) # bias correction (variance)
      {
        gamma <- sum(PROD) # don't muck with weights
        PROD <- PROD / gamma
        PROD <- debias.volume(PROD,bias[i,j])$PMF
        PROD <- PROD * gamma # don't muck with weights
      }
      PMF <- PMF + include[i,j] * PROD
    }
  }

  dV <- prod(UD[[1]]$dr)
  GAMMA <- sum(PMF) # useful for relative encounter rates
  PMF <- PMF / GAMMA

  object <- list()
  object$PDF <- PMF / dV
  object$CDF <- pmf2cdf( PMF )
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
encounter.ctmm <- function(CTMM,include=NULL,exclude=NULL,debias=FALSE,...)
{
  # check.projections(CTMM)

  # Gaussian / cumulants
  axes <- CTMM[[1]]$axes
  AXES <- length(axes)
  isotropic <- all(sapply(CTMM,function(C){C$isotropic}))

  # Wishart DOFs
  DOF <- sapply(CTMM,DOF.wishart)
  # precision matrices
  P <- lapply(CTMM,function(M){PDsolve(M$sigma)})
  # pairwise DOFs
  dof <- matrix(0,length(DOF),length(DOF))
  for(i in 1:length(DOF))
  {
    for(j in 1:length(DOF))
    {
      dof[i,j] <- 3/2 * tr(P[[i]]+P[[j]])^2 / ( tr(P[[i]])^2/DOF[i] + tr(P[[j]])^2/DOF[j] )
    }
  }
  # clamp DOFs
  # DOF <- clamp(DOF,AXES+2,Inf)
  # dof <- clamp(dof,AXES+2,Inf)
  # relative biases of inverse-Wishart
  # BIAS <- DOF/(DOF-AXES-1)
  # bias <- dof/(dof-AXES-1)
  # Pade approximants with better behavior
  BIAS <- 1 + (AXES+1)/DOF
  bias <- 1 + (AXES+1)/dof

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
        if(debias) # bias correction (precision)
        {
          Pi <- Pi / BIAS[i]
          Pj <- Pj / BIAS[j]
        }
        Pij <- Pi + Pj
        sigma <- PDsolve(Pij)

        mu <- sigma %*% (Pi %*% CTMM[[i]]$mu[1,] + Pj %*% CTMM[[j]]$mu[1,])
        mu <- c(mu)

        # intrinsic weight from multiplying two densities and renormalizing
        include[i,j] <- include[i,j] / sqrt( det(CTMM[[i]]$sigma) * det(CTMM[[j]]$sigma) * det(Pij) )

        # the interplay between the bias correction steps and weighting has to be the same here as for the AKDEs
        if(debias) { sigma <- sigma / bias[i,j] } # bias correction (variance)

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
  COV.mu <- STUFF$COV[I,I,drop=FALSE]
  dimnames(COV.mu) <- list(axes,axes)

  I <- (AXES+1):length(STUFF$MLE)
  sigma <- STUFF$MLE[I]
  sigma <- covm(sigma,axes=axes,isotropic=isotropic)
  COV <- STUFF$COV[I,I,drop=FALSE]
  NAMES <- names(sigma@par)[1:length(I)] # isotropic
  dimnames(COV) <- list(NAMES,NAMES)

  info <- mean.info(CTMM)

  CTMM <- ctmm(mu=mu,sigma=sigma,COV.mu=COV.mu,COV=COV,axes=axes,isotropic=isotropic,info=info)
  return(list(CTMM=CTMM,BIAS=BIAS,bias=bias))
}
