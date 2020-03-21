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
    #
    #
  }

  # initial estimates
  mu <- apply(MU,2,mean)
  COV.mu <- vapply(1:N,function(i) {D <- MU[i,]-mu; D%o%D} ) # AXES,AXES,n
  COV.mu <- apply(COV.mu,1:2,mean)
  dim(COV.mu) <- c(AXES,AXES)

  #####################
  # VARIANCE/COVARIANCE STUFF
  #####################

  # analyticlly solvable
  if(method=="exact" && sufficient=="Wishart" && prior=="Inverse-Wishart")
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
  else if(method=="exact" && sufficient=="log-normal" && prior=="log-normal")
  {
    #
  }

  ####################
  # AIC/BIC/MSPE...
  ####################

}
