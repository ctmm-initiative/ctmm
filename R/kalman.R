# numerically stable evaluation of
# 1 - exp(-x)^2
dexp2 <- function(x,Exp=exp(-x)) { ifelse(Exp<0.7071068,1-Exp^2,2*Exp*sinh(x)) }
# 1 - exp(-x)^1
dexp1 <- function(x,Exp=exp(-x)) { ifelse(Exp<0.5,1-Exp,2*sqrt(Exp)*sinh(x/2)) }

###############################
# Propagator/Green's function and Two-time correlation from Langevin equation for Kalman filter and simulations
# random CTMM objects need to be run through get.taus() first, to precompute various parameters
langevin <- function(dt,CTMM,DIM=1)
{
  K <- CTMM$K
  tau <- CTMM$tau
  sigma <- methods::getDataPart(CTMM$sigma)

  if(K<=1) # IID-BM-OU
  {
    # IID limit
    Green <- array(0,c(1,1))
    Sigma <- array(1,c(1,1))

    if(K)
    {
      if(tau[1]==Inf) # BM
      {
        if(dt<Inf) { Green[1,1] <- 1 }
        # absorbing 1/tau into sigma # VAR -> Diffusion
        Sigma[1,1] <- 2*dt
      }
      else if(dt<Inf) # (BM,OU,IID]
      {
        dtau <- dt/tau
        c0 <- exp(-dtau)
        Green[1,1] <- c0
        Sigma[1,1] <- dexp2(dtau,Exp=c0)
      }
    } # >IID
  } # IID-BM-OU
  else if(K==2) # IOU-OUF-OUO
  {
    Omega2 <- CTMM$Omega2
    #IID limit
    Green <- rbind( c(0,0) , c(0,0) )
    Sigma <- rbind( c(1,0) , c(0,Omega2) )

    f <- CTMM$f.nu[1] # mean(f)
    nu <- CTMM$f.nu[2] # nu || omega
    TT <- CTMM$TfOmega2 # 2 f / Omega^2
    fdt <- f*dt

    if(tau[1]==Inf) # IOU
    {
      dtau <- dt/tau[2]
      Exp <- exp(-dtau)
      DExp <- dexp1(dtau,Exp) # 1-exp(dt/tau[2])

      if(dt<Inf)
      {
        Green[1,1] <- 1
        Green[1,2] <- tau[2]*DExp
        Green[2,2] <- Exp
      }

      # remember that sigma is D=sigma/tau[1]
      DExp2 <- DExp^2 # (1-exp(-dt/tau[2]))^2
      Sigma[1,1] <- clamp( 2*dt - tau[2]*(2*DExp+DExp2) ,0,Inf) # does this still underflow?
      Sigma[2,2] <- dexp2(dtau,Exp)/tau[2]
      if(dt<Inf) { Sigma[c(2,3)] <- clamp(DExp2,0, sqrt(Sigma[1,1]*Sigma[2,2]) ) } # how does this get so far off?
      # 0 at dt=Inf
    } # END IOU
    else if(dt<Inf) # (IOU,OUF/OUO,IID]
    {
      # function representation choice
      nudt <- nu*dt
      EXP <- (tau[1]>tau[2] && nudt>0.8813736)
      if(EXP) # exponential functions
      {
        dtau <- dt/tau
        dift <- diff(tau)
        Exp0 <- exp(-dtau)
        Exp <- Exp0/dift
        c0 <- diff(Exp*tau)
        c1 <- -diff(Exp)
        c2 <- diff(Exp/tau)
      }
      else # trigonometric and hyperbolic-trigonometric functions
      {
        Exp <- exp(-fdt)

        if(tau[1]>tau[2]) # hyperbolic-trigonometric
        {
          Sin0 <- sinh(nudt)
          Sinc0 <- sinch(nudt,Sin0)
          Cos0 <- cosh(nudt)
        }
        else # trigonometric
        {
          Sin0 <- sin(nudt)
          Sinc0 <- sinc(nudt,Sin0)
          Cos0 <- cos(nudt)
        }
        SincE <- Sinc0*Exp
        CosE <- Cos0*Exp

        c0 <- CosE + fdt*SincE
        c1 <- -(Omega2*dt)*SincE
        c2 <- -Omega2*(CosE - fdt*SincE)
      } # end function representation

      Green[1,1] <- c0
      Green[2,1] <- c1
      Green[1,2] <- -c1/Omega2
      Green[2,2] <- -c2/Omega2

      # initially canceling terms
      if(EXP)
      {
        dift2 <- dift^2
        T2 <- tau^2
        S1 <- dexp2(dtau[1],Exp0[1])
        S2 <- dexp2(dtau[2],Exp0[2])
        S12 <- 2*tau[1]*tau[2]*dexp1(fdt,Exp0[1]*Exp0[2])
        Sigma[1,1] <- (T2[1]*S1 - S12 + T2[2]*S2)/dift2
        Sigma[2,2] <- (T2[2]*S1 - S12 + T2[1]*S2)/dift2 * Omega2
      }
      else
      {
        CROSS <- fdt*Sinc0*Exp
        OUTER <- Cos0^2*dexp2(fdt,Exp) - CROSS^2
        CROSS <- 2*Cos0*Exp*CROSS
        Sin2 <- Sin0^2

        if(tau[1]>tau[2])
        {
          Sigma[1,1] <- OUTER - Sin2 - CROSS
          Sigma[2,2] <- (OUTER - Sin2 + CROSS) * Omega2
        }
        else
        {
          Sigma[1,1] <- OUTER + Sin2 - CROSS
          Sigma[2,2] <- (OUTER + Sin2 + CROSS) * Omega2
        }
      }

      # initially vanishing terms
      c12 <- c1^2
      Sigma[1,1] <- Sigma[1,1] - c12/Omega2
      Sigma[c(2,3)] <- TT*c12
      Sigma[2,2] <- Sigma[2,2] - c12
    } # end OUF/OUO
  }

  # fix the dimension of the filter
  if(DIM==1) # 1D filter
  { Sigma <- sigma * Sigma }
  else # 2D filter
  {
    if(length(sigma)==1) { sigma <- diag(sigma,DIM) }

    K <- max(1,K)

    Sigma <- outer(Sigma,sigma) # (k,k,d,d)
    Sigma <- aperm(Sigma,c(1,3,2,4)) # (k,k,d,d) -> (k,d,k,d)
    dim(Sigma) <- c(K*DIM,K*DIM)

    # BM/IOU prior fix
    NAN <- is.nan(Sigma)
    if(any(NAN)) { Sigma[NAN] <- 0 }

    Green <- outer(Green,diag(DIM)) # (k,k,d,d)
    Green <- aperm(Green,c(1,3,2,4)) # (k,d,k,d)
    dim(Green) <- c(K*DIM,K*DIM)
  }

  return(list(Green=Green, Sigma=Sigma))
}


Langevin <- function(t,dt=c(Inf,diff(t)),CTMM,DIM=1)
{
  n <- length(dt)
  tau <- CTMM$tau
  K <- max(1,length(tau))  # dimension of hidden state per spatial dimension

  # propagation information
  Green <- array(diag(K*DIM),c(K*DIM,K*DIM,n))
  Green <- aperm(Green,c(3,1,2)) # [n,K*DIM,K*DIM]
  Sigma <- array(0,c(n,K*DIM,K*DIM))

  dynamics <- CTMM$dynamics
  # default stationary process
  if(is.null(dynamics) || dynamics==FALSE || dynamics=="stationary")
  {
    for(i in 1:n)
    {
      # does the time lag change values? Then update the propagators.
      if(i==1 || dt[i] != dt[i-1])
      { LANGEVIN <- langevin(dt=dt[i],CTMM=CTMM,DIM=DIM) }

      Green[i,,] <- LANGEVIN$Green # (K*DIM,K*DIM)
      Sigma[i,,] <- LANGEVIN$Sigma # (K*DIM,K*DIM)
    }
  }
  else if(dynamics=="change.point")
  {
    CP <- CTMM[[dynamics]] # change points
    CS <- get.states(CTMM) # states
    j <- 1
    for(i in 1:n)
    {
      while(t[i]>CP$stop[j]) # does this time increment cross a change point?
      {
        DT <- CP$stop[j] - max(t[i-1],CP$start[j]) # time until change
        LANGEVIN <- langevin(dt=DT,CTMM=CTMM[[CS[j]]],DIM=DIM)
        Green[i,,] <- LANGEVIN$Green %*% Green[i,,]
        Sigma[i,,] <- (LANGEVIN$Green %*% Sigma[i,,] %*% t(LANGEVIN$Green)) + LANGEVIN$Sigma
        dt[i] <- dt[i] - DT
        j <- j + 1
      }

      if(i>1 && CP$start[j]>t[i-1]) # did we cross a change point?
      {
        LANGEVIN <- langevin(dt=dt[i],CTMM=CTMM[[CS[j]]],DIM=DIM)
        Green[i,,] <- LANGEVIN$Green %*% Green[i,,]
        Sigma[i,,] <- (LANGEVIN$Green %*% Sigma[i,,] %*% t(LANGEVIN$Green)) + LANGEVIN$Sigma
      }
      else # we did not cross a change point
      {
        # do we need a fresh calculation?
        if(i==1 || dt[i] != dt[i-1]) { LANGEVIN <- langevin(dt=dt[i],CTMM=CTMM,DIM=DIM) }

        Green[i,,] <- LANGEVIN$Green # (K*DIM,K*DIM)
        Sigma[i,,] <- LANGEVIN$Sigma # (K*DIM,K*DIM)
      }
    } # end for i in 1:n
  } # end change.point dynamics

  R <- list(Green=Green,Sigma=Sigma)
  return(R)
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
kalman <- function(z,u,t=NULL,dt=c(Inf,diff(t)),CTMM,error=NULL,DIM=1,smooth=FALSE,sample=FALSE,residual=FALSE,precompute=FALSE)
{
  # STUFF THAT CAN BE PRECOMPUTED IF DOING MULTIPLE SIMULATIONS
  if(precompute>=0)
  {
    n <- nrow(z)
    if(CTMM$range) { N <- n } else { N <- n-1 } # condition off first point

    if(is.null(error)) { error <- array(0,c(n,DIM,DIM)) }

    CTMM <- get.taus(CTMM,simplify=TRUE) # pre-compute some stuff for Langevin equation solutions
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

    # Propagators from Langevin equation
    LANGEVIN <- Langevin(t=t,dt=dt,CTMM=CTMM,DIM=DIM)
    Green <- LANGEVIN$Green
    Sigma <- LANGEVIN$Sigma

    # first zForcast is properly zeroed
    sFor[1,,] <- Sigma[1,,] # (K*DIM,K*DIM)

    for(i in 1:n)
    {
      # residual covariance
      # 0*Inf -> 0; 0 projection of Inf covariance
      sForP <- nant(sFor[i,,] %*% P,0) # (K*DIM,OBS*DIM)
      sRes[i,,] <- (nant(t(P) %*% sForP,0) + error[i,,]) # (OBS*DIM,OBS*DIM)
      # nant's here might be stop-gap solutions

      # forecast residuals
      zRes[i,,] <- zRes[i,,] - (t(P) %*% zFor[i,,]) # (OBS*DIM,VEC) # zRes is initially just z

      Gain <- sForP %*% PDsolve(sRes[i,,]) # (K*DIM,OBS*DIM) # updated to invert Inf uncertainty to 0
      # 0/0 NaN have Gain of 1 (P) # # Inf*epsilon have Gain of 1 (P)
      INF <- is.nan(Gain) | (Gain==Inf)
      if(any(INF)) { Gain[INF] <- P[INF] }
      #for(j in 1:ncol(Gain)) { if(any(is.nan(Gain[,j])) || any(abs(Gain[,j])==Inf)) { Gain[,j] <- P[,j] } }

      # concurrent estimates
      zCon[i,,] <- zFor[i,,] + (Gain %*% zRes[i,,]) # (K*DIM,VEC)
      # manifestly positive-definite form (Joseph)
      JOSEPH <- IdH - (Gain %*% t(P)) # (K*DIM,K*DIM)
      # this is supposed to be more stable numerically
      sCon[i,,] <- nant(JOSEPH %*% sFor[i,,],0) %*% t(JOSEPH) + nant(Gain %*% error[i,,] %*% t(Gain),0) # (K*DIM,K*DIM) # 0^2/0 NaN have variance of 0
      # sCon[i,,] <- nant(sCon[i,,],0) # 0^2/0 NaN have variance of 0
      # DIAG <- diag(cbind(error[i,,])) # diag() throws an error if you try to secure nrow/ncol
      # if(DIAG<Inf && DIAG>0) { sCon[i,,] <- sCon[i,,] + (Gain %*% error[i,,] %*% t(Gain)) } # don't attempt 0*Inf*0

      # update forecast estimates for next iteration
      if(i<n)
      {
        # update forecast estimates now
        zFor[i+1,,] <- Green[i+1,,] %*% zCon[i,,] # (K*DIM,VEC)
        # work around for when initial location estimates are NULL
        # avoid NaNs from 0*Inf that should be 0
        TCON <- sCon[i,,]
        dim(TCON) <- c(K*DIM,K*DIM)
        INF <- diag(TCON)==Inf
        ANY <- any(INF)
        if(ANY) { diag(TCON)[INF] <- 1i } # will set back to Inf
        TCON <- Green[i+1,,] %*% TCON %*% t(Green[i+1,,]) + Sigma[i+1,,] # (K*DIM,K*DIM)
        # restore Inf variances after avoiding NaNs -- delete finite correlations with infinite variance (impossible?)
        if(ANY)
        {
          INF <- Im(diag(TCON)) > 0
          TCON <- Re(TCON)
          TCON[INF,] <- TCON[,INF] <- 0
          diag(TCON)[INF] <- Inf
        }
        sFor[i+1,,] <- TCON
      }
    }

    # returned likelihood & profiled mean or residuals
    if(!smooth && !sample)
    {
      # !range
      # remove first location from dataset
      # this is mu[1,] with covariance error[1,,]
      # then how does this propagate along other mu terms?

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
      uisRes <- nant(uisRes,0) # 0/0 infinite precision zero shouldn't be anything but zero in weighting

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
        # iW <- PDsolve(W) # (MEAN,MEAN)
        iW <- PDsolve(W,force=TRUE) # (MEAN,MEAN)

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

      mu <- nant(iW %*% D,0) # (MEAN,DATA)
      # if DIM==1, dim(mu)==c(M,2), else dim(mu)==c(2*M,1) - reversed order
      # fix for BM/IOU conditioning off first location (not really a model parameter)
      if(!CTMM$range && length(mu) && CTMM$mean!="zero")
      {
        DIMS <- dim(zRes)
        dim(zRes) <- c(n,OBS,DIM,VEC)
        # overwrite NaN stationary mean value with first location (maximizes likelihood)
        if(DIM==1)
        { mu[1,] <- zRes[1,1,,DATA] }
        else # account for dimension packing in mean
        { mu[1+c(0,nrow(mu)/2),1] <- zRes[1,1,,DATA] }
        dim(zRes) <- DIMS
        # ! TODO: remove stationary mean terms from COV
        # ! TODO: remove stationary mean terms from COV
        # ! TODO: remove stationary mean terms from COV
      }

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
        if(!CTMM$range) { zRes <- zRes[-1,,drop=FALSE] }
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

      sigma <- (zisRes %*% zRes)/(N*DIM)
      # divergence fix
      INF <- diag(sigma)==Inf
      if(any(INF))
      {
        # force positive definite
        sigma[INF,] <- sigma[,INF] <- 0
        # restore infinite variance
        diag(sigma)[INF] <- Inf
      }

      # log det autocorrelation matrix == trace log autocorrelation matrix
      logdet <- apply(sRes,1,det) # just determinants
      # zero variances are improbable, but not balanced out elsewhere in log-like
      logdet <- ifelse(logdet>0,log(abs(logdet)),Inf)
      if(CTMM$range) # this is 1/n times the full term
      { logdet <- mean(logdet) }
      else # this is 1/(n-1) times the full term, after first is dropped
      { logdet <- mean(logdet[-1]) }

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
      # Inf * 0 -> 0
      TL <- nant( sCon[i,,] %*% t(Green[i+1,,]) ,0)
      INV <- PDsolve(sFor[i+1,,],force=TRUE,tol=0)
      # 0/0 & Inf/Inf -> 1
      #NAN <- diag(TL)==0 & diag(INV)==Inf
      # Inf * 1 -> Inf # even though off-diagnals contribute Inf * 0
      # INF <- diag(TL)!=0 & diag(INV)==Inf
      # complete multiplication
      TL <- TL %*% INV # (K*DIM,K*DIM)
      NAN <- is.nan(diag(TL))
      # now take limits
      if(any(NAN))
      {
        TL[NAN,] <- 0
        TL[,NAN] <- 0
        diag(TL)[NAN] <- 1 # Inf/Inf->1
      }
      # if(any(INF))
      # {
      #   TL[INF,] <- 0
      #   TL[,INF] <- 0
      #   diag(TL)[INF] <- Inf # Inf*1->Inf
      # }
      L[i,,] <- TL
      # manifestly symmetric backsolver
      # sCon[i,,] <- sCon[i,,] + L %*% (sCon[i+1,,]-sFor[i+1,,]) %*% t(L)
      # manifestly PSD backsolver
      JOSEPH <- IdH - (L[i,,] %*% Green[i+1,,]) # (K*DIM,K*DIM)
      sCon[i,,] <- nant(JOSEPH %*% sCon[i,,],0) %*% t(JOSEPH) + (L[i,,] %*% (sCon[i+1,,]+Sigma[i+1,,]) %*% t(L[i,,])) # (K*DIM,K*DIM)

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
