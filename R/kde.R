# S3 generic
akde <- function(data,CTMM,VMM=NULL,debias=TRUE,smooth=TRUE,error=0.001,res=10,grid=NULL,...) { UseMethod("akde") }

# fast lag counter
lag.DOF <- function(data,dt=NULL,weights=NULL,lag=NULL,FLOOR=NULL,p=NULL)
{
  t <- data$t
  # intelligently select algorithm
  n <- length(t)

  # important infor for later
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
  DOF[-1] <- DOF[-1]*((1-DOF[1])/sum(DOF[-1]))

  # add positive and negative lags
  # DOF[-1] <- 2*DOF[-1]
  # correct numerical error / normalize to weight sum
  # DOF <- DOF/sum(DOF)

  return(list(DOF=DOF,lag=lag))
}


##################################
# Bandwidth optimizer
#lag.DOF is an unsupported option for end users
bandwidth <- function(data,CTMM,VMM=NULL,weights=FALSE,fast=TRUE,dt=NULL,precision=1/2,PC="Markov",verbose=FALSE,trace=FALSE)
{
  PC <- match.arg(PC,c("Markov","circulant","IID","direct"))

  n <- length(data$t)
  error <- 2^(log(.Machine$double.eps,2)*precision)

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
      dt <- min(diff(data$t))
      if(trace)
      {
        UNITS <- unit(dt,"time")
        message("Default grid size of ",dt/UNITS$scale," ",UNITS$name," chosen for bandwidth(...,fast=TRUE).")
      }
    }

    if(!WEIGHTS) # class(weights)=='numeric' || class(weights)=='integer' ||  # fixed weight lag information
    {
      # if(class(weights)=='numeric' || class(weights)=='integer') # fixed weights
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

    MISE <- stats::optimize(f=MISE,interval=hlim,tol=error)
    h <- MISE$minimum
    MISE <- MISE$objective / sqrt(det(2*pi*sigma)) # constant

    if(!is.null(VMM)) { h <- c(h,h) }
  }
  else # 3D
  {
    h <- 4/5/n^(1/7) # User Silverman's rule of thumb as initial guess
    MISE <- stats::optim(par=c(h,h),fn=MISE,control=list(maxit=.Machine$integer.max,reltol=error))
    h <- MISE$par
    MISE <- MISE$value
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
    { bias <- 1 - 1/n + h^2 }
    else
    {
      # weights were optimized and now DOF can be calculated
      if(is.null(DOF) & fast)
      {
        DOF <- lag.DOF(data,dt=dt,weights=weights,lag=lag,FLOOR=FLOOR,p=p)
        lag <- DOF$lag
        DOF <- DOF$DOF
      }

      bias <- akde.bias(CTMM,H=H[1:2,1:2],lag=lag,DOF=DOF,weights=weights)

      if(!is.null(VMM))  { bias[3] <- akde.bias(VMM,H=VH,lag=lag,DOF=DOF,weights=weights) }
    }

    names(bias) <- axes
    names(h) <- axes
    names(DOF.area) <- axes

    H <- list(H=H,h=h,DOF.H=DOF.H,bias=bias,DOF.area=DOF.area,weights=weights,MISE=MISE)
    class(H) <- "bandwidth"
  }

  if(trace) { message("Bandwidth optimization complete.") }
  return(H)
}

# bias of Gaussian Reference function AKDE
# generalize to non-stationary mean
akde.bias <- function(CTMM,H,lag,DOF,weights)
{
  sigma <- methods::getDataPart(CTMM$sigma)

  # weighted correlation
  ACF <- Vectorize( svf.func(CTMM,moment=FALSE)$ACF )
  if(is.null(DOF)) # exact version
  {
    COV <- lag # copy structure
    COV[] <- ACF(lag) # preserve structure... why do I have to do it this way?
    COV <- c(weights %*% COV %*% weights)
  }
  else # fast version
  { COV <- sum(DOF*ACF(lag)) }

  # variance inflation factor
  bias <- ( det(cbind((1-COV)*sigma + H))/det(cbind(sigma)) )^(1/length(CTMM$axes))
  # remove cbind if/when I can get det.numeric working

  # name dimensions of bias
  axes <- CTMM$axes
  bias <- rep(bias,length(axes))
  names(bias) <- axes

  return(bias)
}


###################
# homerange wrapper function
homerange <- function(data,CTMM,method="AKDE",...)
{ akde(data,CTMM,...) }


#######################################
# wrap the kde function for our telemetry data format and CIs.
akde.telemetry <- function(data,CTMM,VMM=NULL,debias=TRUE,smooth=TRUE,error=0.001,res=10,grid=NULL,...)
{
  if(class(CTMM)=="ctmm") # calculate bandwidth etc.
  {
    axes <- CTMM$axes

    # smooth out errors (which also removes duplicate times)
    z <- NULL
    if(!is.null(VMM))
    {
      axis <- VMM$axes
      if(VMM$error && smooth) { z <- predict(VMM,data=data,t=data$t)[,axis] }
      axes <- c(axes,axis)
    }
    if(CTMM$error && smooth)
    {
      data <- predict(CTMM,data=data,t=data$t)
      CTMM$error <- FALSE # smoothed error model (approximate)
    }
    # copy back !!! times and locations
    if(!is.null(VMM))
    {
      data[,axis] <- z
      VMM$error <- FALSE # smoothed error model (approximate)
    }

    # calculate optimal bandwidth and some other information
    KDE <- bandwidth(data=data,CTMM=CTMM,VMM=VMM,verbose=TRUE,...)
  }
  else if(class(CTMM)=="bandwidth") # bandwidth information was precalculated
  {
    KDE <- CTMM
    axes <- names(KDE$h)
  }
  else
  { stop(paste("CTMM argument is of class",class(CTMM))) }

  if(debias) { debias <- KDE$bias }

  # absolute resolution
  dr <- sqrt(diag(KDE$H))/res

  KDE <- c(KDE,kde(data,KDE$H,axes=axes,bias=debias,W=KDE$weights,alpha=error,dr=dr,grid=grid))

  KDE <- new.UD(KDE,info=attr(data,"info"))

  return(KDE)
}


# if a single H matrix is given, make it into an array of H matrices
prepare.H <- function(H,n)
{
  # one matrix given
  if(length(dim(H))==2)
  {
    d <- nrow(H)
    H <- array(H,c(d,d,n))
    H <- aperm(H,c(3,1,2))
  }

  return(H)
}


############################################
# construct a grid for the density function
kde.grid <- function(data,H,axes=c("x","y"),alpha=0.001,res=1,dr=NULL)
{
  R <- get.telemetry(data,axes) # (times,dim)
  n <- nrow(R) # (times)

  H <- prepare.H(H,n) # (times,dim,dim)

  # how far to extend range from data as to ensure alpha significance in total probability
  z <- sqrt(-2*log(alpha))
  dH <- z * apply(H,1,function(h){sqrt(diag(h))}) # (dim,times)
  dH <- t(dH) # (times,dim)

  # now to find the necessary extent of our grid
  EXT <- rbind( apply(R-dH,2,min) , apply(R+dH,2,max) ) # (ext,dim)
  dEXT <- EXT[2,]-EXT[1,]

  # grid center
  mu <- apply(EXT,2,mean)

  if(is.null(dr)) { dr <- dEXT/res }

  res <- dEXT/dr

  # half resolution
  res <- ceiling(res/2+1)

  # grid locations
  R <- lapply(1:length(res),function(i){ (-res[i]:res[i])*dr[i] + mu[i] } ) # (grid,dim)
  names(R) <- axes

  grid <- list(R=R,dH=dH,dr=dr,EXT=EXT,dEXT=dEXT)
  return(grid)
}


##################################
# construct my own kde objects
# was using ks-package but it has some bugs
# alpha is the error goal in my total probability
kde <- function(data,H,axes=c("x","y"),bias=FALSE,W=NULL,alpha=0.001,res=NULL,dr=NULL,grid=NULL)
{
  r <- get.telemetry(data,axes)
  n <- nrow(r)

  # normalize weights
  if(is.null(W)) { W <- rep(1,length(data$x)) }
  W <- W/sum(W)

  # format bandwidth matrix
  H <- prepare.H(H,n)

  if(is.null(grid)) { grid <- kde.grid(data,H=H,axes=axes,alpha=alpha,res=res,dr=dr) }

  R <- grid$R
  # generalize this for future grid option use
  dH <- grid$dH
  dr <- grid$dr

  R0 <- sapply(R,first)
  # corner origin to minimize arithmetic later
  #r <- t(t(r) - R0)
  #R <- lapply(1:length(R0),function(i){ R[[i]] - R0[i] })

  # probability mass function
  PMF <- array(0,sapply(R,length))
  for(i in 1:n)
  {
    # sub-grid lower/upper bound indices
    i1 <- floor((r[i,]-dH[i,]-R0)/dr) + 1
    i2 <- ceiling((r[i,]+dH[i,]-R0)/dr) + 1

    SUB <- lapply(1:length(i1),function(d){ i1[d]:i2[d] })

    # I can't figure out how to do this in one line
    if(length(SUB)==1)
    { PMF[SUB[[1]]] <- PMF[SUB[[1]]] + W[i]*pnorm1(R[[1]][SUB[[1]]]-r[i,1],H[i,,],dr,alpha) }
    else if(length(SUB)==2)
    { PMF[SUB[[1]],SUB[[2]]] <- PMF[SUB[[1]],SUB[[2]]] + W[i]*pnorm2(R[[1]][SUB[[1]]]-r[i,1],R[[2]][SUB[[2]]]-r[i,2],H[i,,],dr,alpha) }
    else if(length(SUB)==3)
    { PMF[SUB[[1]],SUB[[2]],SUB[[3]]] <- PMF[SUB[[1]],SUB[[2]],SUB[[3]]] + W[i]*pnorm3(R[[1]][SUB[[1]]]-r[i,1],R[[2]][SUB[[2]]]-r[i,2],R[[3]][SUB[[3]]]-r[i,3],H[i,,],dr,alpha) }
  }

  if(sum(bias)) # debias area/volume
  {
    if(length(dr)==2) # AREA debias
    {
      # debias the area
      PMF <- debias.volume(PMF,bias=min(bias))
      CDF <- PMF$CDF
      PMF <- PMF$PMF
    }
    else if(length(dr)==3) # VOLUME debias
    {
      # I'm assuming z-bias is smallest
      vbias <- min(bias)
      abias <- min(bias[1:2])/min(bias)

      # volume then area correction
      PMF2 <- debias.volume(PMF,bias=vbias)$PMF
      PMF2 <- debias.area(PMF2,bias=abias)$PMF

      # area then volumen correction
      PMF <- debias.area(PMF,bias=abias)$PMF
      PMF <- debias.volume(PMF,bias=vbias)$PMF

      # average the two orders to cancel lowest-order errors
      PMF <- (PMF+PMF2)/2
      rm(PMF2)

      # retabulate the cdf
      CDF <- pmf2cdf(PMF)
    }
  }
  else # finish off PDF
  { CDF <- pmf2cdf(PMF) }

  dV <- prod(dr)
  result <- list(PDF=PMF/dV,CDF=CDF,r=R,dr=dr)
  return(result)
}

################################
# debias the entire volume of the PDF
debias.volume <- function(PMF,bias=1)
{
  # cdf: cell probability -> probability included in contour
  CDF <- pmf2cdf(PMF,finish=FALSE)
  DIM <- CDF$DIM
  IND <- CDF$IND
  CDF <- CDF$CDF

  # convert from variance bias to volume bias
  bias <- sqrt(bias)^length(DIM)

  # counting volume by dV
  VOL <- 1:length(CDF)

  # evaluate the debiased cdf on the original volume grid
  CDF <- stats::approx(x=VOL/bias,y=CDF,xout=VOL,yleft=0,yright=1)$y

  # recalculate pdf
  PMF <- cdf2pmf(CDF)

  # back in spatial order # back in table form
  PMF[IND] <- PMF
  PMF <- array(PMF,DIM)

  CDF[IND] <- CDF
  CDF <- array(CDF,DIM)

  return(list(PMF=PMF,CDF=CDF))
}

#####################################
# this is really a foliation debias of 2D slice area in 3D space
debias.area <- function(PMF,bias=1)
{
  DIM <- dim(PMF)

  # debias the CDF areas slice by slice
  for(i in 1:DIM[3])
  {
    # slice probability mass
    dP <- sum(PMF[,,i])
    if(dP) # only do this if there is some probability to shape
    {
      # normalize slice
      PMF[,,i] <- PMF[,,i]/dP

      # debias slice areas (equiprobability ok)
      PMF[,,i] <- debias.volume(PMF[,,i],bias=bias)$PMF

      # un-normalize slice
      PMF[,,i] <- PMF[,,i]*dP
    }
  }

  # not using this ATM
  # CDF <- pmf2cdf(PMF)

  return(list(PMF=PMF,CDF=NULL))
}

########################
# cdf: cell probability -> probability included in contour
pmf2cdf <- function(PMF,finish=TRUE)
{
  #cdf <- pdf * dV
  DIM <- dim(PMF)
  PMF <- c(PMF) # flatten table

  PMF <- sort(PMF,decreasing=TRUE,method="quick",index.return=TRUE)

  IND <- PMF[[2]] # sorted indices
  PMF <- PMF[[1]]

  PMF <- cumsum(PMF)

  if(finish)
  {
    # back in spatial order # back in table form
    PMF[IND] <- PMF
    PMF <- array(PMF,DIM)
    return(PMF)
  }
  else { return(list(CDF=PMF,IND=IND,DIM=DIM)) }
}

###################
# assume already sorted, return sorted
cdf2pmf <- function(CDF)
{
  # this method has some numerical noise from discretization/aliasing
  PMF <- diff(c(0,CDF))

  # ties or something is causing numerical errors with this transformation
  # I don't completely understand the issue, but it seems to follow a pattern
  DECAY <- diff(PMF)
  # this quantity should be negative, increasing
  for(i in 1:length(DECAY))
  {
    # when its positive, its usually adjoining a comparable negative value
    if(DECAY[i]>0)
    {
      # find the adjoining error
      j <- which.min(DECAY[i + -1:1]) - 2 + i
      # tie out the errors
      DECAY[i:j] <- mean(DECAY[i:j])
    }
  }
  PMF <- -rev(cumsum(c(0,rev(DECAY))))

  return(PMF)
}


#######################
# robust bi-variate CDF (mean zero assumed)
#######################
pnorm2 <- function(X,Y,sigma,dr,alpha=0.001)
{
  dx <- dr[1]
  dy <- dr[2]
  cdf <- array(0,c(length(X),length(Y)))

  # eigensystem of kernel covariance
  v <- eigen(sigma)
  s <- v$values

  # effective degree of degeneracy at best resolution
  ZERO <- sum(s <= 0 | (min(dx,dy)/2)^2/s > -2*log(alpha))

  # correlation
  S <- sqrt(sigma[1,1]*sigma[2,2])
  if(S>0)
  {
    rho <- sigma[1,2]/S
    # prevent some tiny numerical errors just in case
    rho <- clamp(rho,min=-1,max=1)
  }
  else { rho <- 0 }

  # main switch
  if(ZERO==0 && abs(rho)<1) # no degeneracy
  {
    # relative grid size (worst case)
    z <- sqrt((dx^2+dy^2)/s[2])

    if(z^3/12<=alpha) # midpoint integration
    {
      cdf <- (dx*dy) * Gauss(X,Y,sigma)
    }
    else if(z^5/2880<=alpha) # Simpson integration
    {
      W <- c(1,4,1)
      cdf <- NewtonCotes(X,Y,sigma,W,dx,dy)
    }
    else if(z^7/1935360<=alpha) # Boole integration
    {
      W <- c(7,32,12,32,7)
      cdf <- NewtonCotes(X,Y,sigma,W,dx,dy)
    }
    else # exact calculation
    {
      # offset to corners
      x <- c(X-dx/2,last(X)+dx/2)
      y <- c(Y-dy/2,last(Y)+dy/2)

      # standardized all locations
      x <- (x)/sqrt(sigma[1,1])
      y <- (y)/sqrt(sigma[2,2])

      # dimensions
      n.x <- length(x)
      n.y <- length(y)

      # corner grid of cell probabilities
      CDF <- outer(x,y,function(x,y){pbivnorm::pbivnorm(x,y,rho=rho)})

      # integrate over cell and add
      cdf <- CDF[-1,-1] - CDF[-n.x,-1] - CDF[-1,-n.y] + CDF[-n.x,-n.y]

      # pbivnorm is very fragile
      # if(any(is.nan(cdf))) { stop(" dx=",dx," dy=",dy," sigma=",sigma," x=",x," y=",y," rho=",rho," CDF=",cdf)}
    }
  }
  else if(ZERO==1 || abs(rho)==1) # line degeneracy
  {
    # max variance
    s <- s[1]
    # unit vector of max variance
    v <- v$vectors[,1]

    # crossings along X grid
    x.cell <- c()
    y.cross <- c()
    m.y <- v[2]/v[1]
    if(abs(m.y)<Inf)
    {
      x.cell <- c(X-dx/2,last(X)+dx/2)
      y.cross <- (x.cell)*m.y
    }

    # crossings along Y grid
    y.cell <- c()
    x.cross <- c()
    m.x <- v[1]/v[2]
    if(abs(m.x)<Inf)
    {
      y.cell <- c(Y-dy/2,last(Y)+dy/2)
      x.cross <- (y.cell)*m.x
    }

    # all crossings
    x.cross <- c(x.cell,x.cross)
    y.cross <- c(y.cross,y.cell)

    # standardized location along line
    z.cross <- ((x.cross)*v[1] + (y.cross)*v[2])/sqrt(s)
    z.cross <- sort(z.cross,method="quick")
    z.cross <- unique(z.cross)

    for(i in 1:(length(z.cross)-1))
    {
      # what cell is this line segment in?
      z.mid <- mean(z.cross[i:(i+1)])

      x <- sqrt(s)*v[1]*z.mid
      y <- sqrt(s)*v[2]*z.mid

      r <- abs(x-X) <= dx/2
      c <- abs(y-Y) <= dy/2

      cdf[r,c] <- (stats::pnorm(z.cross[i+1])-stats::pnorm(z.cross[i]))/(sum(r)*sum(c))
    }
  }
  else if(ZERO==2) # point degeneracy
  {
    # the closest point(s)
    r <- abs(X) <= dx/2
    c <- abs(Y) <= dy/2

    # increment the closest point(s)
    cdf[r,c] <- cdf[r,c] + 1/(sum(r)*sum(c))
  }
  else stop("something is wrong in matrix: sigma == ",sigma)

  return(cdf)
}


# This function is not ready for Kriging
pnorm3 <- function(X,Y,Z,sigma,dr,alpha=0.001)
{
  cdf <- prod(dr) * Gauss3(X,Y,Z,sigma)

  return(cdf)
}

# UNFINISHED
pnorm1 <- function(X,sigma,dr,alpha=0.001) { 0 }

#################
# Newton-Cotes integrators
NewtonCotes <- function(X,Y,sigma,W,dx=mean(diff(X)),dy=mean(diff(Y)))
{
  W <- W/sum(W)
  n <- length(W)
  m <- n-1

  # refined grid to have Simpson's rule in between X,Y points
  # changed from to= to length.out= to avoid roundoff error
  x <- seq(from=X[1]-dx/2,by=dx/m,length.out=length(X)*m+1)
  y <- seq(from=Y[1]-dy/2,by=dy/m,length.out=length(Y)*m+1)

  # weight arrays
  w.x <- dx * array(W[-n],length(x))
  w.y <- dy * array(W[-n],length(y))

  # weight table
  W <- (w.x %o% w.y)

  cdf <- W * Gauss(x,y,sigma)

  # coarsen grid
  # index order is (x,y)
  cdf <- vapply(1:length(X)-1,function(i){colSums(cdf[1:n+m*i,])},rep(0,length(y)))
  # index order is (y,X)
  cdf <- vapply(1:length(Y)-1,function(i){colSums(cdf[1:n+m*i,])},rep(0,length(X)))
  # index order is (X,Y)

  return(cdf)
}

#####################
# gaussian pdf
Gauss <- function(X,Y,sigma=NULL,sigma.inv=solve(sigma),sigma.GM=sqrt(det(sigma)))
{
  cdf <- outer(X^2*sigma.inv[1,1],Y^2*sigma.inv[2,2],"+")/2
  cdf <- cdf + (X %o% Y)*sigma.inv[1,2]
  cdf <- exp(-cdf)/(2*pi*sigma.GM)
  return(cdf)
}


#####################
# gaussian pdf
# assumes uncorrelated z-axis
Gauss3 <- function(X,Y,Z,sigma=NULL,sigma.inv=solve(sigma[1:2,1:2]),sigma.GM=sqrt(det(sigma[1:2,1:2])))
{
  cdf <- Gauss(X,Y,sigma=sigma[1:2,1:2],sigma.inv=sigma.inv,sigma.GM=sigma.GM)
  cdf <- cdf %o% (exp(-Z^2/(2*sigma[3,3]))/sqrt(2*pi*sigma[3,3]))

  return(cdf)
}


#####################
# AKDE CIs
CI.UD <- function(object,level.UD=0.95,level=0.95,P=FALSE)
{
  if(is.null(object$DOF.area) && P)
  {
    names(level.UD) <- "ML"
    return(level.UD)
  }

  dV <- prod(object$dr)

  # point estimate
  area <- sum(object$CDF <= level.UD) * dV
  names(area) <- "ML"

  # chi square approximation of uncertainty
  if(!is.null(object$DOF.area))
  {
    area <- chisq.ci(area,DOF=2*object$DOF.area[1],alpha=1-level)
    names(area) <- c("low","ML","high")
  }

  if(!P) { return(area) }

  # probabilities associated with these areas
  P <- round(area / dV)

  # fix lower bound
  P[1] <- max(1,P[1])
  # fix upper bound to not overflow
  P[3] <- min(length(object$CDF),P[3])
  if(P[3]==length(object$CDF)) { warning("Outer contour extends beyond raster.") }

  P <- sort(object$CDF,method="quick")[P]

  # recorrect point estimate level
  P[2] <- level.UD

  names(P) <- c("low","ML","high")
  return(P)
}

#######################
# summarize details of akde object
summary.UD <- function(object,level.UD=0.95,level=0.95,units=TRUE,...)
{
  # do we convert base units?
  if(units) { thresh <- 1 } else { thresh <- Inf }

  area <- CI.UD(object,level.UD,level)

  # pretty units
  unit.info <- unit(area[2],"area",thresh=thresh)
  name <- unit.info$name
  scale <- unit.info$scale

  area <- array(area/scale,c(1,3))
  rownames(area) <- paste("area (",name,")",sep="")

  colnames(area) <- c("low","ML","high")

  SUM <- list()

  SUM$DOF <- c(object$DOF.area[1],object$DOF.H)
  names(SUM$DOF) <- c("area","bandwidth")

  SUM$CI <- area

  return(SUM)
}
#methods::setMethod("summary",signature(object="UD"), function(object,...) summary.UD(object,...))
