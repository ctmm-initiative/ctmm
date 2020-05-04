# population-level parameter estimates for normally distributed parameters and parameter uncertainties
meta.normal <- function(MU,SIGMA,debias=TRUE,isotropic=FALSE,precision=1/2)
{
  N <- dim(MU)[1]
  DIM <- dim(MU)[2]

  tol <- .Machine$double.eps^precision
  REML <- debias

  # initial guesses
  mu <- colMeans(MU)
  sigma <- 0
  for(i in 1:N) { sigma <- sigma + outer(MU[i,]-mu) }
  sigma <- sigma/(N-REML)

  ERROR <- Inf
  loglike <- loglike.old <- -Inf
  while(ERROR>tol && loglike>=loglike.old)
  {
    loglike.old <- loglike
    sigma.old <- sigma

    # estimate mu exactly
    P <- array(0,c(N,DIM,DIM))
    mu <- P.mu <- 0
    for(i in 1:N)
    {
      P[i,,] <- PDsolve(sigma + SIGMA[i,,])
      P.mu <- P.mu + P[i,,]
      mu <- mu + c(P[i,,] %*% MU[i,])
    }
    COV.mu <- PDsolve(P.mu)
    mu <- c(COV.mu %*% mu)

    loglike <- REML/2*log(det(COV.mu)) + DIM*(N-REML)/2*log(2*pi)
    # gradient with respect to sigma
    RHS <- 0
    LHS <- P.mu
    for(i in 1:N)
    {
      D <- mu - MU[i,]
      RHS <- RHS + (P[i,,] %*% outer(D) %*% P[i,,])
      if(debias) { LHS <- LHS - (P[i,,] %*% COV.mu %*% P[i,,]) }
      loglike <- loglike - 1/2*log(det(sigma + SIGMA[i,,])) - 1/2*c(D %*% P[i,,] %*% D)
    }

    K <- sqrtm.covm(covm(sigma)) %*% mpow.covm(covm(LHS),-1/2)
    sigma <- K %*% RHS %*% t(K)

    # Standardized error
    ERROR <- sigma - sigma.old # absolute error
    K <- mpow.covm(covm(sigma),-1/2)
    ERROR <- K %*% ERROR %*% K # standardize to make ~1
    ERROR <- ERROR %*% ERROR # square to make positive
    ERROR <- sqrt(mean(diag(ERROR %*% ERROR))) # error in standard deviations
  }

  DUP <- upper.tri(sigma,diag=TRUE)

  # we still need hessian(sigma) for sigma CIs
  log.like <- function(par,REML=TRUE)
  {
    sigma <- array(0,c(DIM,DIM))
    sigma[DUP] <- par
    sigma <- t(sigma)
    sigma[DUP] <- par

    # sum up log-likelihood
    loglike <- REML/2*log(det(COV.mu)) + DIM*(N-REML)/2*log(2*pi)
    for(i in 1:N)
    {
      D <- mu - MU[i,]
      loglike <- loglike - 1/2*log(det(sigma + SIGMA[i,,])) - 1/2*c(D %*% P[i,,] %*% D)
    }
  }

  par <- sigma[DUP]

  parscale <- sqrt( diag(sigma) )
  parscale <- sigma %o% sigma
  parscale <- parscale[DUP]

  lower <- array(-Inf,c(DIM,DIM))
  diag(lower) <- 0

  DIFF <- genD(par,log.like,parscale=parscale,lower=lower,upper=Inf)
  COV.sigma <- PDsolve(-DIFF$HESS)

  loglike <- log.like(par,REML=FALSE)

  # TODO AIC
  # TODO exact AICc

  return(list(mu=mu,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike))
}


###########
# log-transformed parameters
# debias includes bias correction for chi^2 to log(chi^2)
# matrix casts location covariance, diffusion rate, velocity covariance all as distinct matrices for above bias correction
log.ctmm <- function(CTMM,features,debias=TRUE)
{
  par <- get.parameters(CTMM,features)
  COV <- CTMM$COV

  #!!! re-write to diagonalize sigma first
  #!!! then log transform everything
  #!!! then diagonalize everything according to COV
  #!!! then find element that most overlaps zero off-diagonal
  #!!! bias correct all other elements as if chi^2

  # features to log transform
  SUB <- features[features %in% dimnames(COV)[[1]] && features %in% POSITIVE.PARAMETERS && par>0]

  # Jacobian for log transformation
  J.new <- function()
  {
    J <- diag(1,nrow(COV))
    dimnames(J) <- dimnames(COV)
    return(J)
  }

  # log transform positive parameters
  J <- J.new()
  for(s in SUB)
  {
    J[s,s] <- 1/par[s]
    par[s] <- log(par[s])
  }
  # transform covariance (from logarithms)
  COV <- J %*% COV %*% t(J)

  # diagonalize and log-chi^2 debias relevant parameter estimates
  EIGEN <- eigen(COV[SUB,SUB])
  par[SUB] <- t(EIGEN$vectors) %*% par[SUB] # diagonalize parameters
  DOF <- 2/EIGEN$values # log-chi^2 VAR-DOF relation
  par[SUB] <- par[SUB] + (digamma(DOF/2)-log(DOF/2)) # log-chi^2 bias correction
  par[SUB] <- EIGEN$vectors %*% par[SUB] # transform back (still under logarithm)

  # transform from major-minor-angle to xx-yy-xy basis for aggregation
  if(!CTMM$isotropic)
  {
    # transform COV back
    J <- J.new()
    sigma <- exp(par[c('major','minor')])
    J['major','major'] <- sigma['major']
    J['minor','minor'] <- sigma['minor']
    COV <- J %*% COV %*% t(J)

    # now transform the whole matrix forward under logarithm
    J <- J.new()

    theta <- par["angle"]
    u <- c(cos(theta),sin(theta))
    v <- c(-sin(theta),cos(theta))
    NAMES <- c("major","minor","angle")
    DUP <- c(1,4,2) # "xx","yy","xy" # upper triangle

    par[NAMES] <- par['major']*(u%o%u)[DUP] + par['minor']*(v%o%v)[DUP]
    J[NAMES,NAMES] <- cbind( (u%o%u)[DUP]/sigma['major'], (v%o%v)[DUP]/sigma['minor'], (par['major']-par['minor'])*(u%o%v+v%o%u)[DUP] )
    COV <- J %*% COV %*% t(J)
  }

  # fill out missing VAR with Inf after transformation --- prevent NaNs
  TCOV <- diag(Inf,length(features))
  dimnames(TCOV) <- list(features,features)
  NAMES <- dimnames(COV)[[1]]
  TCOV[NAMES,NAMES] <- COV[NAMES,NAMES]
  COV <- TCOV

  return(list(par=par,COV=COV))
}

#####################
# inverse transformation of above
exp.ctmm <- function(sigma,COV,debias=TRUE)
{

}


###########
mean.ctmm <- function(x,sufficient="Wishart",prior="Inverse-Wishart",method="exact",debias=TRUE,precision=1/2,...)
{
  sufficient <- match.arg(summary,c("Wishart","chisq","log-normal"))
  prior <- match.arg(prior,c("Inverse-Wishart","log-normal"))
  method <- match.arg(method,c("exact","Laplace","MCMC"))

  tol <- .Machine$double.eps^precision

  axes <- x[[1]]$axes
  AXES <- length(axes)
  isotropic <- FALSE # for now

  N <- length(x)

  ####################
  # MEAN STUFF
  ####################
  # Gaussian-Gaussian in all cases
  MU <- array(0,c(N,AXES))
  SIGMA <- array(0,c(N,AXES,AXES))
  for(i in 1:N)
  {
    MU[i,] <- x[[i]]$mu
    SIGMA[i,,] <- x[[i]]$COV.mu

    # fill in with zeroes for non-stationary means
    # !!!
    # !!!
  }
  STUFF <- meta.normal(MU,SIGMA,debias=debias)
  mu <- STUFF$mu
  COV.mu <- STUFF$COV.mu
  sigma.mu <- STUFF$sigma
  COV.sigma.mu <- STUFF$COV.sigma


  #####################
  # VARIANCE/COVARIANCE STUFF
  #####################
  features <- x[[1]]$features

  # analyticlly solvable
  if(method %in% c("exact","Laplace") && sufficient=="log-normal" && prior=="log-normal")
  {
    DIM <- length(features)
    MU <- array(0,c(N,DIM))
    SIGMA <- array(0,c(N,DIM,DIM))

    for(i in 1:N)
    {
      STUFF <- log.ctmm(x[[1]])
      log.sigma <- STUFF$log.sigma
      P <- length(log.sigma)
      MU[1:P] <- log.sigma
      J <- diag(1,c(DIM,DIM))
      J[1:P,1:P] <- STUFF$d.log.sigma
    }


    # debias shift


  }
  else if(method=="exact" && sufficient=="Wishart" && prior=="Inverse-Wishart")
  {
    SIGMA <- array(0,c(N,AXES,AXES))
    DOF <- array(0,N)

    for(i in 1:N)
    {
      # extract covariance matrix SIGMA and covariance of covariance matrix COV
      SIGMA[i,,] <- x[[i]]$sigma # Wishart matrix / n
      if(x[[i]]$isotropic) # chi^2
      { PAR <- 'major' }
      else # approximate Wishart DOF (exact if Wishart)
      { PAR <- c('major','minor') }
      EST <- SIGMA[[i]]@par[PAR]
      DOF[i] <- (2/AXES) * c(EST %*% PDsolve(x[[i]]$COV[PAR,PAR]) %*% EST) # average multiple DOFs if not Wishart
    }

    # EM algorithm
    nu <- sum(DOF) # initial estimate that seems reasonble
    S <- Reduce('+',vapply(1:N,function(i){ (nu+DOF[i])*DOF[i]*SIGMA[i,,,drop=FALSE] },diag(1,AXES))) / sum((nu+DOF)*DOF) # initial estimate (weighted average close to perturbative solution)

    count <- Inf
    while(count>2)
    {
      count <- 0
      L <- -Inf
      ERROR <- Inf
      while(ERROR>tol) # iterative weighted average (UNTESTED)
      {
        S <- Reduce('+', vapply(S %*% PDsolve((nu*S+DOF[i]*SIGMA[i,,,drop=FALSE])/(nu+DOF[i])) %*% S,S) )/N
        L.OLD <- L
        L <- N/2*log(det(S)) - sum( vapply(1:N,function(i){ (nu+DOF[i]) * log(det(nu*S+DOF[i]*S[i,,])) },1) )/2
        ERROR <- L-L.OLD
        count <- count + 1
      }

      L <- -Inf
      ERROR <- Inf
      while(ERROR>tol) # Newton Raphson
      {
        L.OLD <- L
        CONST <- N/2*log(det(nu*S)) - Reduce('+',vapply(1:N,function(i){ log(det(nu*S+DOF[i]*SIGMA[i,,,drop=FALSE])) },S))/2
        L <- nu*CONST - N/2*mpsigamma(nu/2,deriv=-1,dim=AXES) + sum(vapply(1:N,function(i){ mpsigamma((nu+DOF[i])/2,deriv=-1,dim=AXES) },1))/2
        L1 <- CONST - N/2*mpsigamma(nu/2,deriv=0,dim=AXES) + sum(vapply(1:N,function(i){ mpsigamma((nu+DOF[i])/2,deriv=0,dim=AXES) },1))/2
        L2 <- - N/2*mpsigamma(nu/2,deriv=1,dim=AXES) + sum(vapply(1:N,function(i){ mpsigamma((nu+DOF[i])/2,deriv=1,dim=AXES) },1))/2
        nu <- clamp(nu-L1/L2,0,Inf)
        ERROR <- L-L.OLD
        count <- count + 1
      }
    } # end alternating EM algorithm

    like <- function(par)
    {
      nu <- par[1]
      S <- covm(par[-1],axes=axes,isotropic=isotropic)

      R <- N/2*log(det(S)) - sum( vapply(1:N,function(i){ (nu+DOF[i]) * log(det(nu*S+DOF[i]*S[i,,])) },1) )/2
      R <- R - N/2*mpsigamma(nu/2,deriv=-1,dim=AXES) + sum(vapply(1:N,function(i){ mpsigamma((nu+DOF[i])/2,deriv=-1,dim=AXES) },1))/2
      return(R)
    }

    ################
    # covariance matrix of hyper-parameter estimates
    par <- nu
    parscale <- 1
    lower <- 0
    upper <- 0
    NAMES <- 'nu'

    S <- covm(S,axes=axes,isotropic=isotropic)
    par <- c(par,S@par[1])
    parscale <- c(parscale,S@par[1])
    lower <- c(lower,0)
    upper <- c(upper,Inf)
    NAMES <- c(NAMES,"major")

    if(!isotropic)
    {
      par <- c(par,S@par[2:3])
      parscale <- c(parscale,S@par[1],pi/2)
      lower <- c(lower,0,-Inf)
      upper <- c(upper,0,Inf)
      NAMES <- c(NAMES,c('minor','angle'))
    }

    LIKE <- like(par)
    DIFF <- genD(par=par,fn=like,zero=LIKE,lower=lower,upper=upper,parscale=parscale,Richardson=2,mc.cores=1)
    hess <- -DIFF$hessian  # negative log likelihood
    grad <- -DIFF$gradient # negative log likelihood

    # more robust covariance calculation than straight inverse
    COV <- cov.loglike(hess,grad)
    dimnames(COV) <- list(NAMES,NAMES)
    COV <- COV[-1,-1] # nu is a nuisance parameter
  }

  ####################
  # AIC/BIC/MSPE...
  ####################

}
