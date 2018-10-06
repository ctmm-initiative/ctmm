###############################
# Propagator/Green's function and Two-time correlation from Langevin equation for Kalman filter and simulations
# random CTMM objects need to be run through get.taus() first, to precompute various parameters
langevin <- function(dt,CTMM)
{
  K <- CTMM$K
  tau <- CTMM$tau
  sigma <- methods::getDataPart(CTMM$sigma)

  if(K<=1 && (length(tau)==0 || tau[1]==0))
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
    Omega2 <- CTMM$Omega2

    if(dt==Inf) # make this numerically relative in future
    {
      Green <- rbind( c(0,0) , c(0,0) )
      Sigma <- rbind( c(1,0) , c(0,Omega2) )
    }
    else # finite time matrices
    {
      TT <- CTMM$TfOmega2
      f <- CTMM$f.nu[1] # mean(f)
      nu <- CTMM$f.nu[2] # nu || omega

      # function representation choice
      if(tau[1]>tau[2] && nu*dt>0.8813736) # exponential functions
      {
        Exp <- exp(-dt/tau)/diff(tau)
        c0 <- diff(Exp*tau)
        c1 <- -diff(Exp)
        c2 <- diff(Exp/tau)
      }
      else # trigonometric and hyperbolic-trigonometric functions
      {
        Exp <- exp(-f*dt)

        if(tau[1]>tau[2]) # hyperbolic-trigonometric
        {
          SincE <- sinch(nu*dt)*Exp
          CosE <- cosh(nu*dt)*Exp
        }
        else # trigonometric
        {
          SincE <- sinc(nu*dt)*Exp
          CosE <- cos(nu*dt)*Exp
        }

        c0 <- CosE + (f*dt)*SincE
        c1 <- -(Omega2*dt)*SincE
        c2 <- -Omega2*(CosE - (f*dt)*SincE)
      } # end function representation

      Green <- rbind( c(c0,-c1/Omega2) , c(c1,-c2/Omega2) )
      Sigma <- -TT*c1^2  #off-diagonal term
      Sigma <- rbind( c(1,0)-c(c0^2+c1^2/Omega2,Sigma) , c(0,Omega2)-c(Sigma,c1^2+c2^2/Omega2) )
    } # end finite time matrices
  }

  # fix the dimension of the filter
  DIM <- dim(sigma)[1]
  if(is.null(DIM)) # 1D filter
  { Sigma <- sigma * Sigma }
  else # 2D filter
  {
    Sigma <- outer(Sigma,sigma) # (k,k,d,d)
    Sigma <- aperm(Sigma,c(1,3,2,4)) # (k,k,d,d) -> (k,d,k,d)
    dim(Sigma) <- c(K*DIM,K*DIM)

    Green <- outer(Green,diag(DIM)) # (k,k,d,d)
    Green <- aperm(Green,c(1,3,2,4)) # (k,d,k,d)
    dim(Green) <- c(K*DIM,K*DIM)
  }

  return(list(Green=Green, Sigma=Sigma))
}

#############################################################
# Internal Kalman filter/smoother for multiple derivatives, dimensions, trends
# Kalman filter/smoother for matrix P operator and multiple mean functions
# this is a lower level function
# more for generalizability/readability than speed at this point
#############################################
# DIM is spatial dimension
# K is dynamic dimension
############################################
# precompute argument: use precomputed data & solutions
# # FALSE: don't assume or return computed glob
# # +1: store a computed glob in the environment
# # -1: use a computed glob from the environment
kalman <- function(z,u,dt,CTMM,error=NULL,smooth=FALSE,sample=FALSE,residual=FALSE,precompute=FALSE)
{
  # STUFF THAT CAN BE PRECOMPUTED IF DOING MULTIPLE SIMULATIONS
  if(precompute>=0)
  {
    n <- nrow(z)

    DIM <- dim(CTMM$sigma)[1] # infer necessary dimension
    if(is.null(DIM)) { DIM <- 1 } # default dimension, scalar sigma
    if(is.null(error)) { error <- array(0,c(n,DIM,DIM)) }

    tau <- CTMM$tau
    K <- max(1,length(tau))  # dimension of hidden state per spatial dimension
    OBS <- 1 # observed dynamical dimensions: OBS <= K

    DATA <- 1:(ncol(z)/DIM) # DATA indices in following array dim

    # glob data and mean functions together for Kalman filter
    z <- cbind(z,u) ; rm(u)
    VEC <- ncol(z)/(OBS*DIM)
    dim(z) <- c(n,OBS*DIM,VEC)

    # indices of mean functions
    if(last(DATA)==VEC) { MEAN <- NULL } else { MEAN <- (last(DATA)+1):VEC }

    Id <- diag(OBS*DIM) # check if still need this !!!
    IdH <- diag(K*DIM) # used for Joseph form

    # observable state projection operator: t(P) %*% FULL -> OBS
    P <- array(0,c(K,DIM,OBS,DIM))
    P[1,,1,] <- diag(1,DIM) # positions recorded 1:1
    dim(P) <- c(K*DIM,OBS*DIM)

    # forecast estimates for zero-mean, unit-variance filters
    zFor <- array(0,c(n,K*DIM,VEC))
    sFor <- array(0,c(n,K*DIM,K*DIM))

    # concurrent/filtered estimates
    zCon <- zFor # (n,K*VEC)
    sCon <- sFor # (n,K*DIM,K*DIM)

    # forcast residuals - not everything is observed
    # zRes <- array(0,c(n,OBS*DIM,VEC))
    zRes <- z ; rm(z)
    sRes <- array(0,c(n,OBS*DIM,OBS*DIM))

    # propagation information
    Green <- array(0,c(n,K*DIM,K*DIM))
    Sigma <- array(0,c(n,K*DIM,K*DIM))

    # Propagators from Langevin equation
    CTMM <- get.taus(CTMM) # pre-compute some stuff for Langevin equation solutions
    for(i in 1:n)
    {
      # does the time lag change values? Then update the propagators.
      if(i==1 || dt[i] != dt[i-1])
      { Langevin <- langevin(dt=dt[i],CTMM=CTMM) }

      Green[i,,] <- Langevin$Green # (K*DIM,K*DIM)
      Sigma[i,,] <- Langevin$Sigma # (K*DIM,K*DIM)
    }

    # first zForcast is properly zeroed
    sFor[1,,] <- Sigma[1,,] # (K*DIM,K*DIM)

    for(i in 1:n)
    {
      # residual covariance
      sForP <- sFor[i,,] %*% P # (K*DIM,OBS*DIM)
      sRes[i,,] <- ((t(P) %*% sForP) + error[i,,]) # (OBS*DIM,OBS*DIM)

      # forcast residuals
      zRes[i,,] <- zRes[i,,] - (t(P) %*% zFor[i,,]) # (OBS*DIM,VEC) # zRes is initially just z

      Gain <- sForP %*% PDsolve(sRes[i,,]) # (K*DIM,OBS*DIM) # updated to invert Inf uncertainty to 0

      # concurrent estimates
      zCon[i,,] <- zFor[i,,] + (Gain %*% zRes[i,,]) # (K*DIM,VEC)
      # manifestly positive-definite form (Joseph)
      JOSEPH <- IdH - (Gain %*% t(P)) # (K*DIM,K*DIM)
      sCon[i,,] <- (JOSEPH %*% sFor[i,,] %*% t(JOSEPH)) # (K*DIM,K*DIM)
      DIAG <- diag(cbind(error[i,,])) # diag() throws an error if you try to secure nrow/ncol
      if(DIAG<Inf && DIAG>0) { sCon[i,,] <- sCon[i,,] + (Gain %*% error[i,,] %*% t(Gain)) } # don't attempt 0*Inf*0
      # this is supposed to be more stable numerically

      # update forcast estimates for next iteration
      if(i<n)
      {
        #update forcast estimates now
        zFor[i+1,,] <- Green[i+1,,] %*% zCon[i,,] # (K*DIM,VEC)
        sFor[i+1,,] <- ((Green[i+1,,] %*% sCon[i,,] %*% t(Green[i+1,,])) + Sigma[i+1,,]) # (K*DIM,K*DIM)
      }
    }

    # returned likelihood & profiled mean or residuals
    if(!smooth && !sample)
    {
      # special code for DIM-1 circle - rotate (u,0) x-trend -> (0,u) y-trend and pretend like we ran y-trend through the Kalman filter
      if(DIM==1 && CTMM$circle)
      {
        u <- zRes[,,MEAN] # (n,OBS*DIM,2*M)
        dim(u) <- c(n,OBS*DIM,2,length(MEAN)/2)
        u <- aperm(u,c(3,4,1,2)) # (2,M,n,OBS*DIM)
        dim(u) <- c(2,length(MEAN)/2*n*OBS*DIM) # spatial dependence first
        # rotate x trends to y trends
        R <- cbind( 0:1 , -(1:0) ) # R : x -> y
        u <- R %*% u
        dim(u) <- c(length(MEAN),n*OBS*DIM)
        u <- aperm(u,c(2,1)) # (n*OBS*DIM,2*M)
        dim(u) <- c(n*OBS*DIM,length(MEAN))
        # riffle y-trend residuals in between x-trend residuals
        dim(zRes) <- c(n*OBS*DIM,VEC)
        zRes <- cbind(zRes[,DATA],riffle(zRes[,MEAN],u,by=2)) # riffle x & y trend terms
        VEC <- VEC + length(MEAN)
        MEAN <- c(MEAN,last(MEAN) + 1:length(MEAN))
        dim(zRes) <- c(n,OBS*DIM,VEC)
        rm(u)
      }

      isRes <- vapply(1:n,function(i){PDsolve(sRes[i,,])},diag(OBS*DIM)) # (OBS*DIM,OBS*DIM,n)
      dim(isRes) <- c(OBS*DIM,OBS*DIM,n) # R arrays are awful
      uisRes <- vapply(1:n,function(i){isRes[,,i] %*% zRes[i,,MEAN]},array(0,c(OBS*DIM,length(MEAN)))) # (OBS*DIM,MEAN,n) - dont need data terms
      dim(uisRes) <- c(OBS*DIM,length(MEAN),n) # R arrays are awful

      ### everything below is hard-coded for OBS==1 ###
      # if DIM==1, none of this does anything interesting #
      dim(zRes) <- c(n*OBS*DIM,VEC)
      # match above for space-time trace
      uisRes <- aperm(uisRes,c(2,3,1)) # c(MEAN,n,DIM)
      dim(uisRes) <- c(length(MEAN),n*OBS*DIM)

      # tracing (inner product) over times (and space if DIM>1)
      if(!(DIM==1 && CTMM$circle)) # (MEAN,MEAN) quadratic terms
      {
        D <- as.matrix(uisRes %*% zRes[,DATA]) # (MEAN,DATA) quadratic terms
        W <- as.matrix(uisRes %*% zRes[,MEAN]) # (MEAN,MEAN)
        iW <- PDsolve(W) # (MEAN,MEAN)

        M.MEAN <- MEAN
        M.DATA <- DATA
      }
      else # don't waste time calculating 0 correlations between x and y (riffled) trends - also DIM-2 spatial-trace calculations
      {
        # DIM-2 indices
        M.DATA <- 1:(length(DATA)/2)
        M.MEAN <- seq(last(M.DATA)+1,VEC/2)

        # include spatial contraction
        dim(uisRes) <- c(2,length(M.MEAN),n*OBS*DIM) # (2,M,n)
        uisRes <- aperm(uisRes,c(2,3,1)) # (M,n,2)
        dim(uisRes) <- c(length(M.MEAN),n*OBS*DIM*2) # (M,n*2)
        dim(zRes) <- c(n*OBS*DIM*2,VEC/2)

        D <- as.matrix(uisRes %*% zRes[,M.DATA]) # (M,D)

        # faster and more accurate calculations
        W <- iW <- matrix(0,length(M.MEAN),length(M.MEAN)) # (M,M)
        SUB1 <- seq(1,length(M.MEAN),2) # x,x terms
        SUB2 <- seq(2,length(M.MEAN),2) # y,y terms
        W[SUB1,SUB1] <- W[SUB2,SUB2] <- (uisRes[SUB1,] %*% zRes[,length(M.DATA)+SUB1]) # x,x==y,y  &&  x,y=0
        iW[SUB1,SUB1] <- iW[SUB2,SUB2] <- PDsolve(W[SUB1,SUB1])
      }
      # if DIM==1, dim(D)==c(M,2), else dim(D)==c(2*M,1) - reversed order
      rm(uisRes)

      mu <- iW %*% D # (MEAN,DATA)
      # if DIM==1, dim(mu)==c(M,2), else dim(mu)==c(2*M,1) - reversed order

      # separate data from trend
      u <- zRes[,M.MEAN,drop=FALSE] # (n*OBS*DIM,MEAN)
      zRes <- zRes[,M.DATA,drop=FALSE] # (n*OBS*DIM,DATA)
      # detrend data
      u <- u %*% mu # (n*OBS*DIM,DATA)
      zRes <- zRes - u # (n*OBS*DIM,DATA)
      rm(u)

      dim(zRes) <- c(n,OBS*DIM,length(DATA)) # go back to original DIM for sRes multiplication

      # return (standardized) residuals of Kalman filter only
      if(residual)
      {
        sRes <- vapply(1:n,function(i){PDfunc(sRes[i,,],function(m){1/sqrt(abs(m))})},diag(OBS*DIM)) # (OBS*DIM,OBS*DIM,n)
        dim(sRes) <- c(OBS*DIM,OBS*DIM,n) # R arrays are awful
        zRes <- vapply(1:n,function(i){sRes[,,i] %*% zRes[i,,]},array(0,c(OBS*DIM,length(DATA)))) # (OBS*DIM,length(DATA),n)
        dim(zRes) <- c(OBS*DIM*length(DATA),n) # flatten
        zRes <- t(zRes) # (n,length(M.DATA))
        return(zRes)
      }

      # reformat mean coefficients for pre-transform structure u %*% mu
      # if DIM==1, dim(mu)==c(M,2), else dim(mu)==c(2*M,1) - reversed order
      if(DIM>1 || CTMM$circle) { mu <- array(mu,c(prod(dim((mu)))/2,2)) }
      # hard coded for 2D !

      # variance/covariance from detrended terms - preventing numerical errors
      # you can't do this simple matrix multiplication in one line because R's ability to handle array dimensions is awful
      zisRes <- vapply(1:n,function(i){isRes[,,i] %*% zRes[i,,]},array(0,c(OBS*DIM,length(DATA)))) # (OBS*DIM,DATA,n)

      # 1x1 sRes multiplication no longer needed - can fully promote DIM-1 circle stuff to DIM-2
      if(DIM==1 && CTMM$circle)
      {
        DIM <- 2
        DATA <- 1:(length(DATA)/2)
      }

      dim(zisRes) <- c(OBS*DIM,length(DATA),n) # R arrays are awful
      zisRes <- aperm(zisRes,c(2,3,1)) # c(DATA,n,DIM)
      dim(zisRes) <- c(length(DATA),n*OBS*DIM)
      # match on right for space-time trace
      dim(zRes) <- c(n*OBS*DIM,length(DATA))
      sigma <- (zisRes %*% zRes)/(n*DIM)

      # log det autocorrelation matrix == trace log autocorrelation matrix
      logdet <- mean(log(apply(sRes,1,det))) # this is 1/n times the full term

      return(list(mu=mu,W=W,iW=iW,sigma=sigma,logdet=logdet))
    }

    # delete residuals
    rm(zRes,sRes)

    #####################
    # KALMAN SMOOTHER
    #####################
    # Finish detrending the effect of a stationary mean and delete u(t) ## we do this in simulate/predict now ##
    # zFor <- zFor[,DATA,drop=FALSE] - (zFor[,MEAN,drop=FALSE] %*% mu) # (n,DATA)
    # zCon <- zCon[,DATA,drop=FALSE] - (zCon[,MEAN,drop=FALSE] %*% mu) # (n,DATA)
    # why does R drop dimensions with such strange peculiarity?

    #################
    # RANDOM SAMPLER
    #################

    # initialize new matrices
    L <- array(0,c(n,K*DIM,K*DIM))
    SQRT <- array(0,c(n,K*DIM,K*DIM)) # standard deviation matrix

    if(sample)
    {
      # this is all we need precomputed for generating random location
      SQRT[n,,] <- PDfunc(sCon[n,,],func=function(x){sqrt(abs(x))},pseudo=TRUE) # (K*DIM,K*DIM)
      # with a known location, collapse the state uncertainty
      sCon[n,,] <- 0
    }

    # upgrade concurrent estimates to Kriged estimates (covariance only)
    for(i in (n-1):1)
    {
      L[i,,] <- sCon[i,,] %*% t(Green[i+1,,]) %*% PDsolve(sFor[i+1,,]) # (K*DIM,K*DIM)
      # manifestly symmetric backsolver
      # sCon[i,,] <- sCon[i,,] + L %*% (sCon[i+1,,]-sFor[i+1,,]) %*% t(L)
      # manifestly PSD backsolver
      JOSEPH <- IdH - (L[i,,] %*% Green[i+1,,]) # (K*DIM,K*DIM)
      sCon[i,,] <- (JOSEPH %*% sCon[i,,] %*% t(JOSEPH)) + (L[i,,] %*% (sCon[i+1,,]+Sigma[i+1,,]) %*% t(L[i,,])) # (K*DIM,K*DIM)

      #################
      # RANDOM SAMPLER
      #################
      if(sample)
      {
        # get what we need to generate random state and collapse uncertainty
        SQRT[i,,] <- PDfunc(sCon[i,,],func=function(x){sqrt(abs(x))},pseudo=TRUE) # (K*DIM,K*DIM)
        sCon[i,,] <- 0
      }
    }

    rm(Green,Sigma,sFor)
  } # END PRECOMPUTE BLOCK

  ###############
  # CAN START WITH PRECOMPUTE HERE
  if(precompute){ STUFF <- c('n','K','DIM','DATA','zCon','sCon','zFor','L','SQRT') }
  # STORE PRECOMPUTED STUFF FOR LATER EVALUATIONS || PULL PRECOMPUTED STUFF FOR FAST EVALUATIONS
  if(precompute>0) { for(thing in STUFF) { assign(thing,get(thing),pos=Kalman.env) } }
  else if(precompute<0) { for(thing in STUFF) { assign(thing,get(thing,pos=Kalman.env)) } }

  #################
  # RANDOM SAMPLER: end point
  #################
  if(sample) { zCon[n,,] <- zCon[n,,] + SQRT[n,,] %*% array(stats::rnorm(K*DIM*length(DATA)),c(K*DIM,length(DATA))) } # (K*DIM,length(DATA))

  # smooth/simulate data finish
  for(i in (n-1):1)
  {
    # RTS smoother
    zCon[i,,] <- zCon[i,,] + L[i,,] %*% (zCon[i+1,,]-zFor[i+1,,]) # (K*DIM)

    # RANDOM SAMPLER
    if(sample) { zCon[i,,] <- zCon[i,,] + SQRT[i,,] %*% array(stats::rnorm(K*DIM*length(DATA)),c(K*DIM,length(DATA))) } # (K*DIM,length(DATA))
  }

  # name dimensions
  znames <- c(outer(c('position','velocity')[1:K],CTMM$axes[1:DIM],paste))
  dimnames(zCon) <- list(NULL,znames,CTMM$axes[DATA])
  dimnames(sCon) <- list(NULL,znames,znames)

  # return smoothed states
  # this object is temporary
  state <- list(CTMM=CTMM,Z=zCon,S=sCon)
  class(state) <- "state"

  return(state)
}


# store Kalman filter/smoother objects to pass between evaluations
Kalman.env <- new.env()
