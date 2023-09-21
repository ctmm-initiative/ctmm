cde <- function(object,include=NULL,exclude=NULL,debias=FALSE,...)
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
  STUFF <- cde.ctmm(CTMM,include=include,exclude=exclude,debias=debias)
  CTMM <- STUFF$CTMM
  BIAS <- STUFF$BIAS # individual biases
  bias <- STUFF$bias # pairwise biases

  DOF.area <- rep(DOF.area(CTMM),AXES)
  names(DOF.area) <- axes

  # bias correction (precision)
  for(i in 1:length(object))
  {
    UD[[i]]$PMF <- UD[[i]]$PDF * prod(UD[[i]]$dr)
    if(debias) { UD[[i]]$PMF <- debias.volume(UD[[i]]$PMF,BIAS[i])$PMF }
  }

  GRID <- grid.union(object) # r,dr of grid union
  DIM <- c(length(GRID$r$x),length(GRID$r$y))
  PMF <- matrix(0,DIM[1],DIM[2]) # initialize CDE PMF

  for(i in 1:(length(UD)-1))
  {
    for(j in (i+1):length(UD))
    {
      # compute overlapping grid subset
      SUB <- grid.intersection(list(GRID,UD[[i]],UD[[j]]))

      if(length(SUB[[1]]$x) && length(SUB[[1]]$y))
      {
        PROD <- UD[[i]]$PMF[SUB[[2]]$x,SUB[[2]]$y] * UD[[j]]$PMF[SUB[[3]]$x,SUB[[3]]$y]
        if(debias) # bias correction (variance)
        {
          gamma <- sum(PROD) # don't muck with weights
          PROD <- PROD / gamma
          PROD <- debias.volume(PROD,bias[i,j])$PMF
          PROD <- PROD * gamma # don't muck with weights
        }
        PMF[SUB[[1]]$x,SUB[[1]]$y] <- PMF[SUB[[1]]$x,SUB[[1]]$y] + include[i,j] * PROD
      }
    }
  }

  dV <- prod(UD[[1]]$dr)
  GAMMA <- sum(PMF) # useful for relative encounter rates
  PMF <- PMF / GAMMA

  object <- GRID
  object$PDF <- PMF / dV
  object$CDF <- pmf2cdf(PMF)
  object$axes <- axes
  object$weight <- GAMMA # store overall weight for future averaging

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

  info <- mean_info(UD)
  object <- new.UD(object,info=info,type='range',variable="encounter",CTMM=CTMM)

  return(object)
}

################
# Gaussian encounters
cde.ctmm <- function(CTMM,include=NULL,exclude=NULL,debias=FALSE,...)
{
  # check.projections(CTMM)

  # Gaussian / cumulants
  axes <- CTMM[[1]]$axes
  AXES <- length(axes)
  isotropic <- sapply(CTMM,function(C){C$isotropic})
  DIM <- ifelse(isotropic,1,AXES)
  WC <- 1 + DIM # inverse-Wishart/chi^2 constant (DOF threshold)
  MULT <- ifelse(isotropic,AXES,1) # DOF multiplier
  WC <- WC/MULT

  CLAMP <- 1 # DOF minimum

  # Wishart DOFs - DOF.wishart() is flaky
  DOF <- sapply(CTMM,DOF.area)
  BIAS <- 1/(1+WC/clamp(DOF,CLAMP,Inf)) # asymptotic

  # finished with this
  isotropic <- all(isotropic)

  # precision matrices
  P <- lapply(CTMM,function(M){PDsolve(M$sigma)})

  # pairwise DOFs (asymptotic)
  dof <- matrix(0,length(DOF),length(DOF))
  # pairwise bias - asymptotic
  bias <- matrix(1,length(DOF),length(DOF))

  for(i in 1:length(DOF))
  {
    for(j in 1:length(DOF))
    {
      Pi <- tr(P[[i]])
      Pj <- tr(P[[j]])
      Pij <- Pi + Pj
      dof[i,j] <- Pij^2 / ( Pi^2/DOF[i] + Pj^2/DOF[j] )
      wc <- (Pi/Pij)*WC[i] + (Pj/Pij)*WC[j] # placeholder interpolation
      bias[i,j] <- 1 + wc/clamp(dof[i,j],2*CLAMP,Inf)
    }
  }

  fn <- function(CTMM)
  {
    IN <- 0
    M1 <- array(0,AXES)
    M2 <- matrix(0,AXES,AXES)

    # precision matrices # have to recalculate these for gradients
    P <- lapply(CTMM,function(M){PDsolve(M$sigma)})

    for(i in 1:(length(CTMM)-1))
    {
      for(j in (i+1):length(CTMM))
      {
        Pi <- P[[i]]
        Pj <- P[[j]]
        if(debias) # bias correction (precision)
        {
          Pi <- Pi * BIAS[i]
          Pj <- Pj * BIAS[j]
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

  info <- mean_info(CTMM)

  CTMM <- ctmm(mu=mu,sigma=sigma,COV.mu=COV.mu,COV=COV,axes=axes,isotropic=isotropic,info=info)
  return(list(CTMM=CTMM,BIAS=BIAS,bias=bias))
}


# relative encounter rates
encounter <- function(object,debias=FALSE,level=0.95,normalize=FALSE,self=TRUE,...)
{
  units <- FALSE

  R <- overlap(object,debias=debias,level=level,method="Rate",...)
  R$CI <- nant(R$CI,0)
  R$DOF <- nant(R$DOF,0)

  if(normalize)
  {
    # calculate mean self-rate
    s <- diag(R$CI[,,'est'])
    dof <- diag(R$DOF)
    s <- meta.chisq(s,dof)$CI["mean","est"]

    R$CI <- R$CI/s
  }

  # fix diagonals # self encounter rate
  if(self)
  {
    diag(R$CI[,,1]) <- diag(R$CI[,,2]) <- diag(R$CI[,,3]) <- Inf
    diag(R$DOF) <- Inf
  }

  R$CI <- R$CI * pi # standardized encounter probability (1-meter)

  return(R)
}
