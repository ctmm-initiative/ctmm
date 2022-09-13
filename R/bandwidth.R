# fast lag counter
lag.DOF <- function(data,dt=NULL,weights=NULL,lag=NULL,FLOOR=NULL,p=NULL)
{
  t <- data$t
  # intelligently select algorithm
  n <- length(t)

  # important information for later
  w2 <- sum(weights^2)
  # otherwise this gets corrupted by gridding/smearing

  # smear the weights over an evenly spaced time grid
  lag <- gridder(t,z=cbind(weights),W=F,dt=dt,lag=lag,FLOOR=FLOOR,p=p,finish=FALSE)
  weights <- lag$z
  lag <- lag$lag

  n <- length(lag)
  N <- composite(2*n)

  DOF <- FFT(pad(weights,N))
  DOF <- abs(DOF)^2
  # include only positive lags
  DOF <- Re(IFFT(DOF)[1:n])

  # fix initial and total DOF
  DOF[1] <- w2
  DOF[-1] <- DOF[-1]*nant((1-DOF[1])/sum(DOF[-1]),0)

  # add positive and negative lags
  # DOF[-1] <- 2*DOF[-1]
  # correct numerical error / normalize to weight sum
  # DOF <- DOF/sum(DOF)

  return(list(DOF=DOF,lag=lag))
}


##################################
# Bandwidth optimizer
#lag.DOF is an unsupported option for end users
bandwidth <- function(data,CTMM,VMM=NULL,weights=FALSE,fast=TRUE,dt=NULL,error=0.01,precision=1/2,PC="Markov",verbose=FALSE,trace=FALSE,...)
{
  # temporary solution
  if(class(CTMM)[1]=="list" && class(CTMM[[1]])[1]=="UD")
  { return(bandwidth.pop(data,CTMM,precision=precision)) }

  PC <- match.arg(PC,c("Markov","circulant","IID","direct"))
  trace2 <- ifelse(trace,trace-1,FALSE)

  if(length(CTMM$tau)==0 || all(CTMM$tau==0)) { weights <- FALSE }

  n <- length(data$t)

  sigma <- methods::getDataPart(CTMM$sigma)
  if(!is.null(VMM)) { sigmaz <- methods::getDataPart(VMM$sigma) }

  if(weights)
  { WEIGHTS <- TRUE }
  else
  {
    WEIGHTS <- FALSE
    weights <- rep(1/n,n)
  }

  if(is.null(CTMM$tau)) # IID bandwidth optimization (faster algorithm)
  {
    if(is.null(VMM))
    {
      MISE <- function(h) { 1/n/(2*h^2) + (n-1)/n/(2+2*h^2) - 2/(2+h^2) + 1/2 }
    }
    else # assuming this would also be IID if horizontal is
    {
      MISE <- function(h) { 1/n/(2*h^2)^(3/2) + (n-1)/n/(2+2*h^2)^(3/2) - 2/(2+h^2)^(3/2) + 1/2^(3/2) }
    }
  }
  else # autocorrelated bandwidth optimization (slower, more general algorithm)
  {
    if(fast & is.null(dt))
    {
      dt <- diff(data$t)
      dt <- dt[dt>0]
      # dt <- sort(dt)
      DT <- stats::median(dt)
      dt.min <- min(dt) # smallest dt
      dt <- dt.min

      if((DT-dt)/DT > error)
      { dt <- DT/2 } # some allowance for error
      else
      { dt <- DT }

      UNITS <- unit(dt,"time")
      STRING <- paste(c(dt/UNITS$scale,UNITS$name),collapse=" ")

      if(trace)
      { message("Default grid size of ",STRING," chosen for bandwidth(...,fast=TRUE).") }

      if(WEIGHTS && dt.min<=dt/2)
      {
        UNITS <- unit(dt.min,"time")
        STRING2 <- paste(c(dt.min/UNITS$scale,UNITS$name),collapse=" ")

        warning("Minimum time interval of ",STRING2," is much smaller than grid size of ",STRING,".")
      }
    }

    if(!WEIGHTS) # class(weights)[1]=='numeric' || class(weights)[1]=='integer' ||  # fixed weight lag information
    {
      # if(class(weights)[1]=='numeric' || class(weights)[1]=='integer') # fixed weights
      # { weights <- weights/sum(weights) }
      # else # uniform weights
      # { weights <- rep(1/n,n) }

      # for fixed weights we can calculate the pair number straight up
      DOF <- lag.DOF(data,dt=dt,weights=weights)
      lag <- DOF$lag
      DOF <- DOF$DOF
    }
    else if(WEIGHTS & fast) # grid information for FFTs
    {
      DOF <- NULL
      lag <- pregridder(data$t,dt=dt)
      FLOOR <- lag$FLOOR
      p <- lag$p
      lag <- lag$lag
    }
    else if(WEIGHTS & !fast) # exact lag matrix
    {
      DOF <- NULL
      p <- NULL
      FLOOR <- NULL
      lag <- outer(data$t,data$t,'-')
      lag <- abs(lag) # SVF not defined correctly for negative lags yet
    }
    else
    { stop("bandwidth weights argument misspecified.") }

    # standardized SVF
    CTMM$sigma <- covm(diag(1,2))
    svf <- svf.func(CTMM,moment=FALSE)$svf
    g <- Vectorize(svf)

    G <- lag #copy structure if lag is matrix
    G[] <- g(lag) # preserve structure... why is this necessary here?

    if(!is.null(VMM)) # 3D AKDE #
    {
      # standardized SVF
      VMM$sigma <- covm(1,axes=VMM$axes,isotropic=VMM$isotropic)
      svfz <- svf.func(VMM,moment=FALSE)$svf
      gz <- Vectorize(svfz)

      GZ <- lag
      GZ[] <- gz(lag)
    }

    # construct approximate inverse for preconditioner
    if(fast & PC=="Toeplitz")
    {
      LAG <- last(lag) + dt + lag
      IG <- g(LAG) # double time domain
      if(!is.null(VMM)) { stop("PC=Toeplitz not implemented in 3D yet.") }
    }
    else
    { IG <- NULL }
    # right now this is just double lag information on SVF

    # Mean Integrated Square Error modulo a constant
    ##########################################################
    if(!is.null(DOF)) # fixed weights, DOF pre-calculated
    {
      MISE <- function(h)
      {
        if(any(h<=0)) { return(Inf) }
        if(is.null(VMM)) # 2D
        { sum(DOF/(G+h^2))/2 - 2/(2+h^2) + 1/2 }
        else # 3D
        { sum(DOF/(G+h[1]^2)/sqrt(GZ+h[2]^2))/2^(3/2) - 2/(2+h[1]^2)/sqrt(2+h[2]^2) + 1/2^(3/2) }
      }
    }
    else if(is.null(DOF)) # solve for weights
    {
      MISE <- function(h)
      {
        if(any(h<=0)) { return(Inf) }

        # MISE matrix (local copy)
        if(is.null(VMM)) #2D
        { G <- 1/(G+h^2)/2 }
        else # 3D
        { G <- (1/2^(3/2))/(G+h[1]^2)/sqrt(GZ+h[2]^2) }

        # finish approximate inverse matrix
        if(fast & PC=="Toeplitz")
        {
          IG <- (1/2)/(IG+h^2) # second-half lags
          IG <- c(G,IG) # combine with first-half lags

          IG <- FFT(IG) # quasi-diagonalization
          IG <- Re(IG) # symmetrize in time
          IG <- 1/IG # quasi-inversion
          IG <- IFFT(IG) # back to time domain
          IG <- Re(IG)
          IG <- IG[1:length(lag)] # halve time domain -- back to original time domain
        }

        if(PC=="Markov")
        {
          MARKOV <- list()
          if(is.null(VMM)) # 2D
          {
            MARKOV$ERROR <- (1/2)/(1+h^2) # asymptote of G, effective white-noise process
            MARKOV$VAR <- (1/2)/(0+h^2) - MARKOV$ERROR # variance of effective OU process
          }
          else # 3D
          {
            MARKOV$ERROR <- (1/2^(3/2))/(1+h[1]^2)/sqrt(1+h[2]^2) # asymptote of G, effective white-noise process
            MARKOV$VAR <- (1/2^(3/2))/(0+h[1]^2)/sqrt(1+h[2]^2) - MARKOV$ERROR # variance of effective OU process
          }

          TAU <- CTMM$tau[1]
          # effective decay timescale for the matrix entries
          if(fast)
          { TAU <- -1/stats::glm.fit(lag,(G-MARKOV$ERROR)/MARKOV$VAR,family=stats::gaussian(link="log"),start=-1/TAU)$coefficients }
          else # just uses the first row... probably not as good but don't want to use all lags and bias to small lags
          { TAU <- -1/stats::glm.fit(lag[1,],(G[1,]-MARKOV$ERROR)/MARKOV$VAR,family=stats::gaussian(link="log"),start=-1/TAU)$coefficients }

          MARKOV$t <- (data$t-data$t[1])/TAU # relative times for an exponential decay process
        }
        else
        { MARKOV <- NULL }

        error <- .Machine$double.eps^precision
        SOLVE <- PQP.solve(G,FLOOR=FLOOR,p=p,lag=lag,error=error,PC=PC,IG=IG,MARKOV=MARKOV)
        weights <<- SOLVE$P
        if(trace){ message(SOLVE$CHANGES," feasibility assessments @ ",round(SOLVE$STEPS/SOLVE$CHANGES,digits=1)," conjugate gradient steps/assessment") }
        if(is.null(VMM)) # 2D
        { return( SOLVE$MISE - 2/(2+h^2) + 1/2 ) }
        else # 3D
        { return( SOLVE$MISE - 2/(2+h[1]^2)/sqrt(2+h[2]^2) + 1/2^(3/2) ) }

        # else # ODL QUADPROG CODE
        # {
        #   dvec <- rep(0,n)
        #   AmaT <- cbind(rep(1,n),diag(n))
        #   bvec <- c(1,rep(0,n))
        #   SOLVE <- quadprog::solve.QP(G,dvec,AmaT,bvec,meq=1)
        #   weights <<- SOLVE$solution
        #   return( 2*SOLVE$value - 2/(2+h^2) + 1/2 )
        # }
      }
    }
  }

  # delete old MISE information, just incase
  empty.env(PQP.env)

  # optimization arguments
  control <- list(precision=precision/2,trace=trace2)
  # reduce precision relative to QP solver

  if(is.null(VMM) || is.null(CTMM$tau)) # 2D or 1 bandwidth parameter
  {
    # bounds for h
    hlim <- c(1/n^(1/6)/2,sqrt(2))
    # hlim[1] is half Silverman's rule of thumb. difference must be asymptotically small.
    # hlim[2] is the exact solution when n=1

    # can't do better than IID data
    if(!is.null(CTMM$tau))
    {
      # analogous model(s) with no autocorrelation
      IID <- CTMM
      IID$tau <- NULL
      VID <- VMM
      VID$tau <- NULL

      hlim[1] <- bandwidth(data,IID,VMM=VID,precision=precision,verbose=TRUE)$h[1]
    }

    # can't do worse than uniform weights
    # if(WEIGHTS && !is.null(CTMM$tau))
    # { hlim[2] <- bandwidth(data,CTMM,VMM=VMM,weights=FALSE,fast=fast,dt=dt,precision=precision,verbose=TRUE)$h[1] }
    # Not sure how helpful this is?

    MISE <- optimizer(sqrt(prod(hlim)),MISE,lower=hlim[1],upper=hlim[2],control=control)
    h <- MISE$par
    MISE <- MISE$value

    #MISE <- stats::optimize(f=MISE,interval=hlim,tol=error)
    #h <- MISE$minimum
    #MISE <- MISE$objective

    MISE <- MISE / sqrt(det(2*pi*sigma)) # constant

    if(!is.null(VMM)) { h <- c(h,h) }
  }
  else # 3D
  {
    #!! could do an IID evaluation here for minimum h !!#
    #!! could do an n=1 evaluation here for maximum h !!#

    h <- 4/5/n^(1/7) # Use Silverman's rule of thumb as initial guess
    MISE <- optimizer(c(h,h),MISE,lower=c(0,0),control=control)
    h <- MISE$par
    MISE <- MISE$value

    #MISE <- stats::optim(par=c(h,h),fn=MISE,control=list(maxit=.Machine$integer.max,reltol=error))
    #h <- MISE$par
    #MISE <- MISE$value

    MISE <- MISE / sqrt(det(2*pi*sigma)) / sqrt(2*pi*sigmaz) # constant
  }

  # delete used MISE information
  empty.env(PQP.env)

  H <- h^2

  if(is.null(VMM))
  {
    DOF.H <- ( 1/(2*H)^2 - 1/(2+2*H)^2 ) / ( 1/(2+H)^2 - 1/(2+2*H)^2 )

    H <- H*sigma

    axes <- CTMM$axes
  }
  else
  {
    # numerator
    DOF.H <- ( 1/(2*H[1])^2/sqrt(2*H[2]) + 1/(2*H[1])/(2*H[2])^(3/2) - 1/(2+2*H[1])^2/sqrt(2+2*H[2]) - 1/(2+2*H[1])/(2+2*H[2])^(3/2) )
    # denominator
    DOF.H <- DOF.H / ( 1/(2+H[1])^2/sqrt(2+H[2]) + 1/(2+H[1])/sqrt(2+H[2])^(3/2) - 1/(2+2*H[1])^2/sqrt(2+2*H[2]) - 1/(2+2*H[1])/sqrt(2+2*H[2])^(3/2) )

    VH <- H[2]*sigmaz
    HH <- H[1]*sigma

    H <- matrix(0,3,3)
    H[1:2,1:2] <- HH
    H[3,3] <- VH

    axes <- c(CTMM$axes,VMM$axes)
  }

  rownames(H) <- axes
  colnames(H) <- axes

  if(verbose)
  {
    h <- c(h[1],h)
    CTMM$sigma <- covm(sigma,axes=CTMM$axes,isotropic=CTMM$isotropic)
    DOF.area <- rep(DOF.area(CTMM),2)

    if(!is.null(VMM))
    {
      VMM$sigma <- covm(sigmaz,axes=VMM$axes,isotropic=TRUE)
      DOF.area[3] <- DOF.area(VMM)
    }

    if(is.null(CTMM$tau))
    {
      bias <- 1 - 1/n + h^2
      COV <- cbind(bias*CTMM$sigma)
    }
    else
    {
      # weights were optimized and now DOF can be calculated
      if(is.null(DOF) & fast)
      {
        DOF <- lag.DOF(data,dt=dt,weights=weights,lag=lag,FLOOR=FLOOR,p=p)
        lag <- DOF$lag
        DOF <- DOF$DOF
      }

      STUFF <- akde.bias(CTMM,H=H[1:2,1:2],lag=lag,DOF=DOF,weights=weights)
      COV <- STUFF$COV
      bias <- STUFF$bias

      if(!is.null(VMM))
      {
        STUFF <- akde.bias(VMM,H=VH,lag=lag,DOF=DOF,weights=weights)
        bias[3] <- STUFF$bias
        ZERO <- rep(0,3)
        COV <- rbind(cbind(COV,ZERO),ZERO)
        COV[3,3] <- STUFF$COV
        dimnames(COV) <- list(axes,axes)
      }
    }

    names(bias) <- axes
    names(h) <- axes
    names(DOF.area) <- axes

    H <- list(H=H,h=h,bias=bias,COV=COV,DOF.H=DOF.H,DOF.area=DOF.area,weights=weights,MISE=MISE)

    if(fast) { H$dt <- dt } # store dt argument used for calculation

    class(H) <- "bandwidth"
  } # end if verbose

  if(trace) { message("Bandwidth optimization complete.") }
  return(H)
}


# calculate bandwidth and weights the sample of a population
bandwidth.pop <- function(data,UD,weights=FALSE,ref="Gaussian",precision=1/2,...)
{
  ref <- match.arg(ref,c("Gaussian","AKDE"))

  CTMM <- lapply(UD,function(ud){ud@CTMM})
  MEAN <- mean(CTMM)
  MEAN <- mean.pop(MEAN)
  axes <- MEAN$axes
  n <- length(data)

  if(all(weights==FALSE)) # uniform weights
  { w <- rep(1/n,n) }
  else # individual weights
  { w <- weights/sum(weights) }
  W.OPT <- all(weights==TRUE) # optimize weights

  # prepare semi-variance lag matrix/vector
  S <- list() # semi-variance
  DOF <- list()
  DET <- rep(NA,length(data))
  for(i in 1:length(data))
  {
    DET[i] <- det(CTMM[[i]]$sigma)

    dt <- UD[[i]]$dt
    if(is.null(dt)) # exact matrix calculation
    {
      lag <- outer(data[[i]]$t,data[[i]]$t,'-')
      lag <- abs(lag) # SVF not defined correctly for negative lags yet
    }
    else
    {
      # for fixed weights we can calculate the pair number straight up
      STUFF <- lag.DOF(data[[i]],dt=dt,weights=UD[[i]]$weights)
      lag <- STUFF$lag
      DOF[[i]] <- STUFF$DOF
    }

    # standardized SVF
    svf <- CTMM[[i]]
    svf$sigma <- covm(diag(1,2))
    svf <- svf.func(svf,moment=FALSE)$svf
    s <- Vectorize(svf)

    S[[i]] <- lag # copy structure if lag is matrix
    S[[i]][] <- s(lag) # preserve structure... why is this necessary here?
  }

  if(ref=="AKDE")
  {
    Q.R <- rates(UD,normalize=FALSE,self=FALSE)[,,'est'] / rates(CTMM,normalize=FALSE,self=FALSE)[,,'est']
    Q.R <- stats::cov2cor(Q.R)

    # MEAN.UD <- mean(UD)
    # L.R <- rates(c(list(MEAN.UD),UD),normalize=FALSE,self=FALSE)[,,'est'] / rates(c(list(MEAN),CTMM),normalize=FALSE,self=FALSE)[,,'est']
    # L.R <- stats::cov2cor(L.R)
    # L.R <- L.R[1,-1]
  }

  MISE <- function(h,finish=FALSE)
  {
    if(h==0) { return(Inf) }
    H <- h^2

    # QUAD
    Q <- matrix(0,length(data),length(data))

    # QUAD diagonal (slowest term - more optimized)
    for(i in 1:length(data))
    {
      DEN <- 1/( sqrt(DET[i])*2*(S[[i]]+H) )

      if(!is.null(dim(S[[i]])))
      { Q[i,i] <- UD[[i]]$weights %*% DEN %*% UD[[i]]$weights }
      else
      { Q[i,i] <- sum(DOF[[i]] * DEN) }
    }

    # QUAD off-diagonals
    # Gaussian
    for(i in 1:length(data))
    {
      for(j in (i+1)%:%length(data))
      {
        MU <- c( CTMM[[i]]$mu - CTMM[[j]]$mu )
        SIG <- (1+H)*(CTMM[[i]]$sigma + CTMM[[j]]$sigma)
        Q[i,j] <- Q[j,i] <- exp(-(MU %*% PDsolve(SIG) %*% MU)/2) / sqrt(det(SIG))
      }
    }

    # LINEAR
    # Gaussian
    L <- rep(0,length(data))
    for(i in 1:length(data))
    {
      MU <- c( CTMM[[i]]$mu - MEAN$mu )
      SIG <- MEAN$sigma + (1+H)*CTMM[[i]]$sigma
      L[i] <- exp(-(MU %*% PDsolve(SIG) %*% MU)/2) / sqrt(det(SIG))
    }

    # correct non-Gaussian overlap
    if(ref=="AKDE")
    {
      Q <- Q * Q.R
      # L <- L * L.R
    }

    if(W.OPT)
    {
      # # minimize w.Q.w - L.w | 1==1.w
      # # minimize w.Q.w - 2*L.w + 2*lambda*(1-1.w)
      # # Q.w - L - lambda*1 = 0
      # # w = solve(Q).(L+lambda*1)
      # # 1 = 1.solve(Q).L + 1.solve(Q).1 * lambda
      # # lambda = (1-1.solve(Q).L)/(1.solve(Q).1)
      # iQ <- PDsolve(Q)
      # lambda <- (1-sum(iQ%*%L))/sum(iQ)
      # w <- c( iQ %*% (L+lambda) )

      CON <- cbind(rep(1,length(L)),diag(length(L)))
      BOUND <- c(1,rep(0,length(L)))
      QPS <- quadprog::solve.QP(Q,L,CON,BOUND,meq=1)
      mise <- 2*QPS$value
    }
    else # fixed weights
    {
      # MISE + constant modulo constant
      mise <- w %*% Q %*% w - 2*(w %*% L)
      mise <- c(mise)
    }

    if(!finish) { return(mise) }

    mise <- ( mise + 1/sqrt(det(2*MEAN$sigma)) ) / (2*pi)

    if(W.OPT)
    {
      w <- QPS$solution
      # fix small numerical errors
      w[w<0] <- 0
      w <- w/sum(w)
    }

    R <- list(MISE=mise,weights=w)
    return(R)
  }

  control <- list(precision=precision)

  # h relative to population COV
  n <- sapply(data,nrow)
  n <- sum(n)
  hmax <- sqrt(2) * (det(MEAN$sigma)/min(DET))^(1/4)
  hlim <- c(1/n^(1/6)/2,hmax)
  h <- optimizer(sqrt(prod(hlim)),MISE,lower=hlim[1],upper=hlim[2],control=control)
  h <- h$par
  # names(h) <- names(data)
  H <- h^2

  STUFF <- MISE(h,finish=TRUE)
  MISE <- STUFF$MISE
  weights <- STUFF$weights
  names(weights) <- names(data)

  # KDE bias
  M1 <- M2 <- 0
  for(i in 1:length(UD))
  {
    MU <- c( CTMM[[i]]$mu - MEAN$mu )
    M1 <- M1 + weights[i] * MU
    M2 <- M2 + weights[i] * ( UD[[i]]$COV + H*CTMM[[i]]$sigma + outer(MU) )
  }
  COV <- cbind(M2) - outer(M1) # sample covariance

  # variance inflation factor
  bias <- ( det(COV)/det(MEAN$sigma) )^(1/length(axes))

  DOF.area <- DOF.area(MEAN)
  DOF.H <- ( 1/(2*H)^2 - 1/(2+2*H)^2 ) / ( 1/(2+H)^2 - 1/(2+2*H)^2 )

  axes <- MEAN$axes
  bias <- c(bias,bias)
  names(bias) <- axes
  h <- c(h,h)
  names(h) <- axes
  DOF.area <- c(DOF.area,DOF.area)
  names(DOF.area) <- axes

  H <- list(h=h,bias=bias,COV=COV,DOF.area=DOF.area,weights=weights,MISE=MISE,CTMM=MEAN)

  class(H) <- "bandwidth"

  return(H)
}