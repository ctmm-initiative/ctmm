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


##########
mean.ctmm.mu <- function(x,debias=TRUE,isotropic=FALSE,...)
{
  axes <- x[[1]]$axes
  AXES <- length(axes)
  N <- length(x)

  # Gaussian-Gaussian in all cases
  MU <- array(0,c(N,AXES))
  SIGMA <- array(0,c(N,AXES,AXES))
  for(i in 1:N)
  {
    MU[i,] <- x[[i]]$mu
    SIGMA[i,,] <- x[[i]]$COV.mu

    # TODO !!!
    # fill in with zeroes for non-stationary means
    # TODO !!!
  }
  R <- meta.normal(MU,SIGMA,debias=debias)
  # R$mu # mean of means
  # R$COV.mu # uncertainty in mean of means estimate
  # rename these
  names(R)[ which(names(R)=="sigma") ] <- "PCOV.mu" # dispersion of means
  names(R)[ which(names(R)=="COV.sigma") ] <- "COV.PCOV.mu" # uncertainty in dispersion of means
  return(R)
}


#############
mean.ctmm.features <- function(x,sufficient="log-normal",prior="log-normal",method="exact",debias=TRUE,precision=1/2,...)
{
  sufficient <- match.arg(summary,c("Wishart","chisq","log-normal"))
  prior <- match.arg(prior,c("Inverse-Wishart","log-normal"))
  method <- match.arg(method,c("exact","Laplace","MCMC"))

  axes <- x[[1]]$axes
  AXES <- length(axes)
  N <- length(x)

  features <- unique( sapply(x,function(y){y$features}) )

  # analyticlly solvable
  if(method %in% c("exact","Laplace") && sufficient=="log-normal" && prior=="log-normal")
  {
    DIM <- length(features)
    MU <- array(0,c(N,DIM))
    SIGMA <- array(0,c(N,DIM,DIM))

    for(i in 1:N)
    {
      R <- log.ctmm(x[[1]],features,debias=debias)
      MU[i,] <- R$par
      SIGMA[i,,] <- R$COV
    }

    # aggregate log parameters
    R <- meta.normal(MU,SIGMA,debias=debias)
    # transform results back
    R <- exp.ctmm(R,debias=debias)
    names(R)[ which(names(R)=="mu") ] <- "par" # mean parameters
    names(R)[ which(names(R)=="COV.mu") ] <- "COV" # mean parameter uncertinaty
    names(R)[ which(names(R)=="sigma") ] <- "PCOV" # parameter dispersion
    names(R)[ which(names(R)=="COV.sigma") ] <- "COV.PCOV" # parameter dispersion uncertainty
    R <- set.parameters(R,R$par)
  }
  else if(method=="exact" && sufficient=="Wishart" && prior=="Inverse-Wishart")
  {
    #
    #
    #
  }
  else
  {
    #
    #
    #
  }

  return(R)
}


###########
mean.ctmm <- function(x,...)
{
  isotropic <- FALSE # for now
  debias <- TRUE

  MU <- mean.ctmm.mu(x,debias=debias,isotropic=isotropic,...)
  CTMM <- mean.ctmm.features(x,debias=debias,isotropic=isotropic)

  # copy features over
  CTMM$mu <- MU$mu
  CTMM$COV.mu <- MU$COV.mu
  CTMM$PCOV.mu <- MU$PCOV.mu
  CTMM$COV.PCOV.mu <- MU$COV.PCOV.mu

  # add log-likelihoods
  CTMM$loglike <- CTMM$loglike + MU$loglike
  CTMM$AIC <- CTMM$AIC + MU$AIC
  CTMM$AICc <- CTMM$AICc + MU$AICc
  CTMM$BIC <- CTMM$BIC + MU$BIC

  # TODO: how do we structure the 3 isotropic types

  # return final result
  return(CTMM)
}


mean.select <- function(x,IC="AICc",...)
{
  # is location mean isotropic
  #
  # is location covariance isotropic
  #
  # is covariance of location covariance isotropic
}
