###############################
# Propagator/Green's function and Two-time correlation from Langevin equation for Kalman filter and simulations
langevin <- function(dt,CTMM)
{
  tau <- CTMM$tau
  sigma <- methods::getDataPart(CTMM$sigma)
  K <- continuity(CTMM)

  if(K==0)
  {
    Green <- array(0,c(1,1))
    Sigma <- array(1,c(1,1))
  }
  else if(K==1)
  {
    c0 <- exp(-dt/tau)
    Green <- array(c0,c(1,1))
    Sigma <- array(1-c0^2,c(1,1))
  }
  else if(K==2)
  {
    f <- 1/tau
    Omega2 <- prod(f)
    nu <- abs(diff(f))/2
    TT <- 2*mean(f)/Omega2

    if(dt==Inf) # make this numerically relative in future
    {
      Green <- rbind( c(0,0) , c(0,0) )
      Sigma <- rbind( c(1,0) , c(0,Omega2) )
    }
    else
    {
      # should we use the exponential or hyperbolic representation?
      if(nu*dt>0.8813736) # exponential
      { DAMP <- TRUE }
      else # hyperbolic
      {
        DAMP <- FALSE
        f <- mean(f)
      }

      if(DAMP) # very overdamped
      {
        Exp <- exp(-dt/tau)/diff(tau)
        c0 <- diff(Exp*tau)
        c1 <- -diff(Exp)
        c2 <- diff(Exp/tau)
      }
      else # mildly overdamped to underdamped
      {
        Exp <- exp(-f*dt)

        SincE <- sinch(nu*dt)*Exp
        CosE <- cosh(nu*dt)*Exp

        c0 <- CosE + (f*dt)*SincE
        c1 <- -(Omega2*dt)*SincE
        c2 <- -Omega2*(CosE - (f*dt)*SincE)
      }

      Green <- rbind( c(c0,-c1/Omega2) , c(c1,-c2/Omega2) )
      Sigma <- -TT*c1^2  #off-diagonal term
      Sigma <- rbind( c(1,0)-c(c0^2+c1^2/Omega2,Sigma) , c(0,Omega2)-c(Sigma,c1^2+c2^2/Omega2) )
    }
  }

  # this needs to be generalized to an outer product for complicated models
  if(is.null(dim(sigma))) { Sigma <- sigma * Sigma }
  else
  {
    Sigma <- outer(Sigma,sigma)
    Sigma <- aperm(Sigma,c(1,3,4,2)) # (k,k,d,d) -> (k,d,k,d)
  }

  # Green will need to be generalized?

  # should I have done this here or in the Kalman filter

  return(list(Green=Green, Sigma=Sigma))
}

#############################################################
# Internal Kalman filter/smoother for multiple derivatives, dimensions, trends
# Kalman filter/smoother for matrix P operator and multiple mean functions
# this is a lower level function
# more for generalizability/readability than speed at this point
############################################
# precompute argument: use precomputed data & solutions
# # FALSE: don't assume or return computed glob
# # +1: store a computed glob in the environment
# # -1: use a computed glob from the environment
kalman <- function(z,u,dt,CTMM,error=rep(0,nrow(z)),smooth=FALSE,sample=FALSE,residual=FALSE,precompute=FALSE)
{
  # STUFF THAT CAN BE PRECOMPUTED IF DOING MULTIPLE SIMULATIONS
  if(precompute>=0)
  {
    n <- nrow(z)
    DATA <- 1:ncol(z)

    # glob data and mean functions together for Kalman filter
    z <- cbind(z,u)
    VEC <- ncol(z)

    # indices of mean functions
    MEAN <- (last(DATA)+1):VEC

    tau <- CTMM$tau
    K <- max(1,length(tau))  # dimension of hidden state per spatial dimension

    # observed dimensions
    OBS <- 1

    Id <- diag(OBS)
    IdH <- diag(K)

    # observable state projection operator (will eventially upgrade to use velocity telemetry)
    P <- array(0,c(K,OBS))
    P[1:OBS,] <- 1

    # forecast estimates for zero-mean, unit-variance filters
    zFor <- array(0,c(n,K,VEC))
    sFor <- array(0,c(n,K,K))

    # forcast residuals
    zRes <- array(0,c(n,OBS,VEC))
    sRes <- array(0,c(n,OBS,OBS))

    # concurrent/filtered estimates
    zCon <- array(0,c(n,K,VEC))
    sCon <- array(0,c(n,K,K))

    # propagation information
    Green <- array(0,c(n,K,K))
    Sigma <- array(0,c(n,K,K))

    # Propagators from Langevin equation
    for(i in 1:n)
    {
      # does the time lag change values? Then update the propagators.
      if(i==1 || dt[i] != dt[i-1])
      { Langevin <- langevin(dt=dt[i],CTMM=CTMM) }

      Green[i,,] <- Langevin$Green
      Sigma[i,,] <- Langevin$Sigma
    }

    # first zForcast is properly zeroed
    sFor[1,,] <- Sigma[1,,]

    for(i in 1:n)
    {
      # residual covariance
      sForP <- sFor[i,,] %*% P # why do I need this?
      ERROR <- error[i]*Id
      sRes[i,,] <- ((t(P) %*% sForP) + ERROR)

      # forcast residuals
      zRes[i,,] <- z[i,] - (t(P) %*% zFor[i,,])

      if(all(abs(sRes[i,,])<Inf)){ Gain <- sForP %*% PDsolve(sRes[i,,]) }
      else { Gain <- sForP %*% (0*Id) } # solve() doesn't like this case...

      # concurrent estimates
      zCon[i,,] <- zFor[i,,] + (Gain %*% zRes[i,,])
      # manifestly symmetric form
      # sCon[i,,] <- (sFor[i,,] - (Gain %*% t(sForP)))
      # manifestly positive-definite form (Joseph)
      JOSEPH <- IdH - (Gain %*% t(P))
      sCon[i,,] <- (JOSEPH %*% sFor[i,,] %*% t(JOSEPH))
      if(error[i]<Inf & error[i]>0) { sCon[i,,] <- sCon[i,,] + (Gain %*% ERROR %*% t(Gain)) } # otherwise Gain==0
      # this is supposed to be more stable numerically

      # update forcast estimates for next iteration
      if(i<n)
      {
        #update forcast estimates now
        zFor[i+1,,] <- Green[i+1,,] %*% zCon[i,,]
        sFor[i+1,,] <- ((Green[i+1,,] %*% sCon[i,,] %*% t(Green[i+1,,])) + Sigma[i+1,,])
      }
    }

    # return (standardized) residuals of Kalman filter only
    if(residual) { return(zRes[,1,]/sqrt(sRes[,1,1])) }

    # general quadratic form, tracing over times
    M <- length(MEAN) + length(DATA)
    M <- matrix(0,M,M)
    # hard-code weights for location observation only
    # not doing this all at once to prevent round-off error
    # M <- Adj(zRes[,1,]) %*% (zRes[,1,]/sRes[,1,1])
    M[MEAN,MEAN] <- Adj(zRes[,1,MEAN]) %*% (zRes[,1,MEAN]/sRes[,1,1])

    # estimate mean parameter
    W <- as.matrix(M[MEAN,MEAN])
    iW <- PDsolve(W)

    M[MEAN,DATA] <- Adj(zRes[,1,MEAN]) %*% (zRes[,1,DATA]/sRes[,1,1])
    # don't think I need adjoint cross term?
    # M[DATA,MEAN] <- Adj(zRes[,1,DATA]) %*% (zRes[,1,MEAN]/sRes[,1,1])

    mu <- iW %*% M[MEAN,DATA]

    # returned profiled mean
    if(!smooth && !sample)
    {
      # Finish detrending the effect of a stationary mean
      MU <- zRes[,1,MEAN]
      dim(MU) <- c(n,length(MEAN))
      MU <- MU %*% mu
      dim(MU) <- c(n,1,length(DATA))
      zRes[,1,DATA] <- zRes[,1,DATA,drop=FALSE] - MU
      # there has to be a better way to do this?

      # hard coded for position observations
      M[DATA,DATA] <- Adj(zRes[,1,DATA]) %*% (zRes[,1,DATA]/sRes[,1,1])
      sigma <- M[DATA,DATA]/n

      # log det autocorrelation matrix == trace log autocorrelation matrix
      logdet <- mean(log(sRes)) # this is 1/n times the full term

      return(list(mu=mu,W=W,iW=iW,sigma=sigma,logdet=logdet))
    }

    # delete residuals
    rm(zRes,sRes)

    #####################
    # KALMAN SMOOTHER
    #####################
    # Finish detrending the effect of a stationary mean
    MU <- zFor[,,MEAN]
    dim(MU) <- c(n*K,length(MEAN))
    MU <- MU %*% mu
    dim(MU) <- c(n,K,length(DATA))
    zFor[,,DATA] <- zFor[,,DATA,drop=FALSE] - MU
    # there has to be a better way to do this?
    MU <- zCon[,,MEAN]
    dim(MU) <- c(n*K,length(MEAN))
    MU <- MU %*% mu
    dim(MU) <- c(n,K,length(DATA))
    zCon[,,DATA] <- zCon[,,DATA,drop=FALSE] - MU
    # why does R drop dimensions so randomly?

    # delete u(t)
    zCon <- zCon[,,DATA,drop=FALSE]
    zFor <- zFor[,,DATA,drop=FALSE]
    # drop=FALSE must be here for BM/OU and I don't fully understand why
  }
  # END PRECOMPUTE BLOCK

  ###############
  # CAN START WITH PRECOMPUTE HERE
  if(precompute){ STUFF <- c('n','K','IdH','u','mu','DATA','zCon','sCon','zFor','sFor','Green','Sigma') }
  # STORE PRECOMPUTED STUFF FOR LATER EVALUATIONS || PULL PRECOMPUTED STUFF FOR FAST EVALUATIONS
  if(precompute>0) { for(thing in STUFF) { assign(thing,get(thing),pos=Kalman.env) } }
  else if(precompute<0) { for(thing in STUFF) { assign(thing,get(thing,pos=Kalman.env)) } }

  #################
  # RANDOM SAMPLER: end point
  #################
  if(sample)
  {
    zCon[n,,] <- sapply(DATA,function(d){MASS::mvrnorm(mu=zCon[n,,d],Sigma=sCon[n,,])})
    sCon[n,,] <- array(0,c(K,K))
  }

  # upgrade concurrent estimates to Kriged estimates
  for(i in (n-1):1)
  {
    # DOUBLE CHECK THIS
    L <- sCon[i,,] %*% t(Green[i+1,,]) %*% PDsolve(sFor[i+1,,])

    # overwrite concurrent estimate with smoothed estimate
    # RTS smoother
    zCon[i,,] <- zCon[i,,] + L %*% (zCon[i+1,,]-zFor[i+1,,])
    # manifestly symmetric backsolver
    # sCon[i,,] <- sCon[i,,] + L %*% (sCon[i+1,,]-sFor[i+1,,]) %*% t(L)
    # manifestly PSD backsolver
    JOSEPH <- IdH - (L %*% Green[i+1,,])
    sCon[i,,] <- (JOSEPH %*% sCon[i,,] %*% t(JOSEPH)) + (L %*% (sCon[i+1,,]+Sigma[i+1,,]) %*% t(L))

    #################
    # RANDOM SAMPLER
    #################
    if(sample)
    {
      # covariance matrices can contain small negative eigen values if dt is small
      SQRT <- PDfunc(sCon[i,,],func=function(x){sqrt(abs(x))},pseudo=TRUE)
      zCon[i,,] <- sapply(DATA,function(d){zCon[i,,d] + (SQRT %*% stats::rnorm(K))})
      # old slower method - solves matrix twice
      # sCon[i,,] <- PDclamp(sCon[i,,]) # prevent small numerical errors from tiny dt
      # zCon[i,,] <- sapply(DATA,function(d){MASS::mvrnorm(mu=zCon[i,,d],Sigma=sCon[i,,])})
      sCon[i,,] <- array(0,c(K,K))
    }
  }

  # restore stationary mean to locations only
  zCon[,1,] <- zCon[,1,] + (cbind(u) %*% mu)

  zname <- c("position")
  if(K>1) { zname <- c(zname,"velocity") }
  dimnames(zCon) <- list(NULL,zname,c("x","y")[DATA])
  dimnames(sCon) <- list(NULL,zname,zname)

  # return smoothed states
  # this object is temporary
  state <- list(CTMM=CTMM,Z=zCon,S=sCon)
  class(state) <- "state"

  return(state)
}


# store Kalman filter/smoother objects to pass between evaluations
Kalman.env <- new.env()
