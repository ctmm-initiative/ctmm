# TODO zero before sum of log.det, quadratic
# TODO multivariate case
ctpm.loglike <- function(data,CTMM,REML=FALSE,profile=TRUE,zero=0,verbose=FALSE)
{
  lag <- attr(data,"lag")
  trait <- c( get.telemetry(data,axes=CTMM$axes) )
  n <- length(trait)
  K <- length(CTMM$tau)

  range <- CTMM$range
  if(!range)
  { REML <- TRUE }
  else # catch bad parameters
  {
    if(length(CTMM$tau) && CTMM$tau[1]==Inf)
    {
      if(profile)
      {
        CTMM$loglike <- -Inf
        return(CTMM)
      }
      else
      {
        return(-Inf)
      }
    }
  }

  if(range) { N <- n } else { N <- n-1 } # condition off marginalized point
  # degrees of freedom
  if(REML) { DOF <- n-1 } else { DOF <- N }
  # REML variance debias factor # ML constant
  VAR.MULT <- N/DOF

  # default identity link function for traits
  link <- get.link(CTMM)
  # full transform is used so that link function can be fitted / selected
  grad <- link$grad(trait)
  trait <- link$fn(trait)

  # include log-like-link in subtracted 'zero'
  zero <- zero - sum(log(abs(grad)))

  # this is just a numerical approximation for matrix inversion
  mu <- mean(trait)
  trait <- trait - mu

  if(!length(CTMM$tau)) # IID analytic solution
  {
    Q <- sum(trait^2)
    log.det <- 0 # IID COR contribution
    COV.mu <- 1/n # modulo variance
  }
  else # requires special preconditioning solver (matrix.R)
  {
    TEST <- CTMM
    TEST$sigma <- covm(1,axes=CTMM$axes,CTMM$isotropic) # unit variance for SVF/ACF below
    COR <- svf.func(TEST)
    if(range) # stationary likelihood
    { COR <- COR$ACF }
    else # likelihood that avoids explicit conditioning
    { COR <- COR$SVF }
    COR <- Vectorize(COR)
    COR <- COR(lag)
    dim(COR) <- dim(lag)

    mu <- mean(trait)

    # data and mean vector
    S <- cbind(trait,rep(1,n))
    # still need 1s for !range Woodbury formulas

    # ACF = var 1o1 - SVF # Woodbury identities ensue
    if(!range) { COR <- -COR }

    S <- PHYsolve(COR,S,CND=!range)
    Q <- S$Q
    log.det <- S$log.det

    # MLE mean
    mu <- mu + Q[1,2]/Q[2,2]

    # modulo VAR[trait] not yet profiled
    COV.mu <- 1/Q[2,2]

    # complete quadratic term
    Q <- Q[1,1] - Q[1,2]*Q[2,1]/Q[2,2]
    # this is also the Woodbury identity for !range
  } # end numerical matrix solver

  # profile the variance / diffusion rate
  if(profile)
  {
    sigma <- Q / DOF
    Q <- 0 # Q - MLE(Q) per N
  }
  else
  {
    sigma <- attr(CTMM$sigma,"par")[1]
    Q <- Q/N/sigma - 1/VAR.MULT # Q - MLE(Q) per N
  }
  COV.mu <- sigma * COV.mu

  # unchanging constants
  LL.CONST <- - 1/2/VAR.MULT - 1/2*log(2*pi)

  loglike <- -1/2*(log.det/N + log(sigma) + Q) # per N
  loglike <- N*(loglike + (LL.CONST-zero/N))

  if(REML || !range) { loglike <- loglike + log(COV.mu)/2 }
  # also from matrix determinant lemma if !range

  loglike <- nant(loglike,-Inf) # Inf - Inf

  if(verbose)
  {
    # assign variables
    if(profile)
    {
      sigma <- covm(sigma,isotropic=TRUE,axes=CTMM$axes)
      CTMM$sigma <- sigma
    }

    CTMM <- ctmm.repair(CTMM,K=K)

    if(range)
    {
      CTMM$mu <- mu
      CTMM$COV.mu <- COV.mu
    }

    CTMM$loglike <- loglike
    #attr(CTMM,"info") <- attr(data,"info")

    CTMM <- ctmm.ctmm(CTMM)
    return(CTMM)
  }
  else
  {
    return(loglike)
  }
}


# almost exact solver for ill-conditioned PD matrices with big uniformly correlated blocks (phylogenetic covariances)
# solves for Q = t(v) * 1/M * v
# solves for log(det(M)) = trace(log(M))
PHYsolve <- function(M,v,CND=FALSE,method="eigen")
{
  n <- nrow(M)
  v <- cbind(v)
  log.det <- 0

  FAIL <- list(Q=diag(Inf,ncol(v)),log.det=Inf)

  if(method=="svd") # this method stays positive longer but turns to garbage without warning
  {
    SVD <- svd(M)

    if(any(abs(Im(SVD$d))>.Machine$double.eps*n) || any(SVD$d<.Machine$double.eps*n))
    { return(FAIL) }

    SVD$d <- abs(SVD$d)

    log.det <- sum(log(SVD$d))

    # average anti-symmetric numerical errors
    U <- Re(SVD$u + SVD$v)/2

    v <- t(U) %*% v
    SVD$d <- 1/SVD$d

    Q <- t(v) %*% (SVD$d * v)
  }
  else if(method=="eigen")
  {
    EIGEN <- eigen(M,TRUE)

    if(any(abs(Im(EIGEN$values))>.Machine$double.eps*n) || last(EIGEN$values)<=.Machine$double.eps*n)
    { return(FAIL) }

    EIGEN$values <- abs(EIGEN$values)

    log.det <- sum(log(EIGEN$values))

    v <- t(EIGEN$vectors) %*% v
    EIGEN$values <- 1/EIGEN$values

    Q <- t(v) %*% (EIGEN$values * v)
  }
  else if(method=="expm")
  {
    # this should be more accurate than eigen() ?
    LOG <- expm::logm(M)
    log.det <- log.det + sum(Re(diag(LOG)))

    # this is supposed to be really good if the logm() is well calculated
    u <- sapply(1:ncol(v),function(i){ expm::expAtv(LOG,v[,i],-1)$eAtv })
    Q <- t(v) %*% u
  }
  else if(method=="solve")
  {
    # UNFINISHED
  }
  else if(method=="Shur")
  {
    # UNFINISHED
  }

  return(list(Q=Q,log.det=log.det))
}
