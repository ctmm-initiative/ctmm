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

  return(list(mu=mu,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike,isotropic=isotropic))
}


###########
# log-transformed parameters
# debias includes bias correction for chi^2 to log(chi^2)
# matrix casts location covariance, diffusion rate, velocity covariance all as distinct matrices for above bias correction
log.ctmm <- function(CTMM,features,debias=TRUE)
{
  isotropic <- CTMM$isotropic
  par <- get.parameters(CTMM,features)
  COV <- CTMM$COV

  ### log transform all positive parameters
  # features to log transform
  COV.NAMES <- dimnames(COV)[[1]]
  SUB <- features[(features %in% COV.NAMES) & (features %in% POSITIVE.PARAMETERS) & (par>0)]

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

  # finish logarithm of sigma matrix (and not just eigen values)
  if(!isotropic)
  {
    angle <- par['angle']
    par['angle'] <- 0 # off-diagonal of log(sigma) after rotation by -angle

    J <- J.new()
    J['angle','angle'] <- par['major'] - par['minor']
    COV <- J %*% COV %*% t(J)
  }

  # log chi^2 bias correction
  if(debias)
  {
    # diagonalize and log-chi^2 debias relevant parameter estimates
    EIGEN <- eigen(COV[SUB,SUB])
    dimnames(EIGEN$vectors) <- list(SUB,SUB)
    names(EIGEN$values) <- SUB
    # fix signs
    if(isotropic) { VAR <- "major" } else { VAR <- c("major","minor") }
    # VAR goes in log numerator for chi^2 variates: variance, diffusion, MS speed, ...
    for(i in 1:nrow(EIGEN$vectors)) { if(sum(EIGEN$vectors[i,VAR])<0) { EIGEN$vectors[i,] <- -EIGEN$vectors[i,] } }
    # transform to diagonalized basis with VARs in log numerator
    par[SUB] <- t(EIGEN$vectors) %*% par[SUB] # diagonalize parameters
    DOF <- 2/EIGEN$values # log-chi^2 VAR-DOF relation
    BIAS <- digamma(DOF/2)-log(DOF/2) # negative bias for log(chi^2) variates
    if(!isotropic)
    {
      # some of the eigen parameter is orientation - which is ~Gaussian and not ~log chi^2 (no bias)
      OVER <- abs(EIGEN$vectors['angle',]) # overlap with orientation eigen-parameter
      BIAS <- (1-OVER)*BIAS # first-order correction (could start at second-order?)
    }
    par[SUB] <- par[SUB] + BIAS # log-chi^2 bias correction
    par[SUB] <- EIGEN$vectors %*% par[SUB] # transform back (still under logarithm)
  }

  # un-diagonalize log(sigma)
  if(!CTMM$isotropic)
  {
    u <- c(cos(angle),sin(angle))
    v <- c(-sin(angle),cos(angle))
    NAMES <- c("major","minor","angle") # input
    UP <- c(1,4,2) # "xx","yy","xy" # upper triangle of log(sigma) # output

    par[NAMES] <- par['major']*(u%o%u)[UP] + par['minor']*(v%o%v)[UP]

    J <- J.new()
    J[NAMES,NAMES] <- cbind( (u%o%u)[UP], (v%o%v)[UP], (u%o%v+v%o%u)[UP] )
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

# orthogonal transformation on matrix -> linear transformation on vector
# t(O) %*% M %*% O -> L %*% m
orth2lin <- function(O,sym=TRUE)
{
  N <- dim(O)[1]
  L <- array(0,c(N,N,N,N))

  for(i in 1:N)
  {
    for(j in 1:N)
    {
      IN <- array(0,c(N,N))
      IN[i,j] <- 1
      OUT <- t(O) %*% IN %*% O
      L[,,i,j] <- OUT
    }
  }

  dim(L) <- c(N^2,N^2)
  return(L)
}

#####################
# inverse transformation of above
exp.ctmm <- function(object,debias=TRUE)
{
  mu <- object$mu # mean log parameters (log chi^2)
  COV.mu <- object$COV.mu # uncertainty in mean log parameters
  sigma <- object$sigma # dispersion of mean logs (determines chi^2 DOFs)
  COV.sigma <- object$COV.sigma # uncertainty in dispersion of mean logs

  isotropic <- object$isotropic
  NAMES <- names(mu)
  N <- length(mu)

  J.new <- function()
  {
    J <- diag(1,N)
    dimnames(J) <- list(NAMES,NAMES)
    return(J)
  }

  # diagonalize log(sigma)
  if(!isotropic)
  {
    SIGMA <- matrix(sigma[c("major","angle","angle","minor")],c(2,2))
    EIGEN <- eigen(SIGMA)
    U <- EIGEN$vectors

    # t(U) %*% SIGMA %*% U
    mu[c('major','minor')] <- EIGEN$values
    mu['angle'] <- 0
    angle <- U[,1] # major axis
    angle <- sign(angle[1]) * angle # positive x component
    angle <- atan2(angle[2],angle[1])

    # transformation matrix J: "xx","yy","xy" -> "major","minor","0-off"
    J <- J.new()
    SIG <- c("major","minor","angle")
    J[SIG,SIG] <- orth2lin(U)[c(1,2,4),c(1,2,4)]

    # linear transform on major, minor, angle
    COV.mu <- J %*% COV.mu %*% t(J)
    sigma <- J %*% sigma %*% t(J)

    # same kind of thing but in even more dimensions
    LJ <- orth2lin(t(J))
    COV.sigma <- LJ %*% COV.sigma %*% LJ
  }

  # reverse bias correction
  if(debias)
  {
    EIGEN <- eigen(COV.mu)
    names(EIGEN$values) <- NAMES
    dimnames(EIGEN$vectors) <- list(NAMES,NAMES)
    # fix signs
    if(isotropic) { VAR <- "major" } else { VAR <- c("major","minor") } # these go in the numerator before log
    for(i in 1:nrow(EIGEN$vectors)) { if(sum(EIGEN$vectors[i,VAR])<0) { EIGEN$vectors[i,] <- -EIGEN$vectors[i,] } }
    # transform to diagonalized basis with VARs in log numerator
    mu <- t(EIGEN$vectors) %*% mu
    DOF <- 2/EIGEN$values # log-chi^2 VAR-DOF relation
    BIAS <- digamma(DOF/2)-log(DOF/2) # negative bias for log(chi^2) variates
    if(!isotropic)
    {
      # some of the eigen parameter is orientation - which is ~Gaussian and not ~log chi^2 (no bias)
      OVER <- abs(EIGEN$vectors['angle',]) # overlap with orientation eigen-parameter
      BIAS <- (1-OVER)*BIAS # first-order correction (could start at second-order?)
    }
    mu <- mu - BIAS # log-chi^2 bias addition
    mu <- EIGEN$vectors %*% mu # transform back (still under logarithm)
  }

  # exp transformation
  SUB <- NAMES %in% POSITIVE.PARAMETERS
  mu[SUB] <- exp(mu[SUB])

  J <- J.new()
  for(s in SUB) { J[s,s] <- mu[s] }

  COV.mu <- J %*% COV.mu %*% t(J) # delta method
  sigma <- J %*% sigma %*% t(J) # not exactly sure about this?
  LJ <- orth2lin(t(J))                 # would be true if above is true
  COV.sigma <- LJ %*% COV.sigma %*% t(LJ) # would be true if above is true

  # undiagonalize sigma-location
  if(!isotropic)
  {
    mu['angle'] <- angle

    # transform uncertainty
    J <- J.new()
    J['angle','angle'] <- 1/(mu['major']-mu['minor'])

    COV.mu <- J %*% COV.mu %*% t(J)

    # transform sigma-par
    sigma <- J %*% sigma %*% t(J)
    LJ <- orth2lin(t(J))
    COV.sigma <- LJ %*% COV.sigma %*% t(LJ)
  }

  return(list(mu=mu,COV.mu=COV.mu,sigma=sigma,COV.sigma=COV.sigma))
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
  mu <- STUFF$mu # mean of means
  COV.mu <- STUFF$COV.mu # uncertainty in mean of means estimate
  PCOV.mu <- STUFF$sigma # dispersion of means
  COV.PCOV.mu <- STUFF$COV.sigma # uncertainty in dispersion of means

  #####################
  # VARIANCE/COVARIANCE STUFF
  #####################
  features <- unique( sapply(x,function(y){y$features}) )

  # analyticlly solvable
  if(method %in% c("exact","Laplace") && sufficient=="log-normal" && prior=="log-normal")
  {
    DIM <- length(features)
    MU <- array(0,c(N,DIM))
    SIGMA <- array(0,c(N,DIM,DIM))

    for(i in 1:N)
    {
      STUFF <- log.ctmm(x[[1]],features,debias=debias)
      MU[i,] <- STUFF$par
      SIGMA[i,,] <- STUFF$COV
    }

    # aggregate log parameters
    STUFF <- meta.normal(MU,SIGMA,debias=debias)
    # transform results back
    STUFF <- exp.ctmm(STUFF,debias=debias)
    par <- STUFF$mu
    COV <- STUFF$COV.mu
    PCOV <- STUFF$sigma
    COV.PCOV <- STUFF$COV.sigma
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
