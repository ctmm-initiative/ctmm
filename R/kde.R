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

  if(length(CTMM$tau)==0 || all(CTMM$tau==0)) { weights <- FALSE }

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

    MISE <- optimizer(sqrt(prod(hlim)),MISE,lower=hlim[1],upper=hlim[2],control=list(precision=precision))
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

    h <- 4/5/n^(1/7) # User Silverman's rule of thumb as initial guess
    MISE <- optimizer(c(h,h),MISE,lower=c(0,0),control=list(precision=precision))
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


# AKDE single or list
# (C) C.H. Fleming (2016-2019)
# (C) Kevin Winner & C.H. Fleming (2016)
akde <- function(data,CTMM,VMM=NULL,debias=TRUE,weights=FALSE,smooth=TRUE,error=0.001,res=10,grid=NULL,...)
{
  if(length(projection(data))>1) { stop("Data not in single coordinate system.") }
  validate.grid(data,grid)
  # force grids to be compatible
  COMPATIBLE <- length(data)>1 && !is.grid.complete(grid)

  DROP <- class(data)[1] != "list"
  data <- listify(data)
  CTMM <- listify(CTMM)
  VMM <- listify(VMM)

  n <- length(data)
  weights <- array(weights,n)

  # loop over individuals for bandwidth optimization
  CTMM0 <- list()
  KDE <- list()
  DEBIAS <- list()
  for(i in 1:n)
  {
    CTMM0[[i]] <- CTMM[[i]] # original model fit
    if(class(CTMM[[i]])[1]=="ctmm") # calculate bandwidth etc.
    {
      axes <- CTMM[[i]]$axes

      # smooth out errors (which also removes duplicate times)
      z <- NULL
      if(!is.null(VMM[[i]]))
      {
        axis <- VMM[[i]]$axes
        if(VMM[[i]]$error && smooth) { z <- predict(VMM[[i]],data=data[[i]],t=data[[i]]$t)[,axis] }
        axes <- c(axes,axis)
      }
      if(CTMM[[i]]$error && smooth)
      {
        data[[i]] <- predict(CTMM[[i]],data=data[[i]],t=data[[i]]$t)
        CTMM[[i]]$error <- FALSE # smoothed error model (approximate)
      }
      # copy back !!! times and locations
      if(!is.null(VMM[[i]]))
      {
        data[[i]][,axis] <- z
        VMM[[i]]$error <- FALSE # smoothed error model (approximate)
      }

      # calculate optimal bandwidth and some other information
      KDE[[i]] <- bandwidth(data=data[[i]],CTMM=CTMM[[i]],VMM=VMM[[i]],weights[i],verbose=TRUE,...)
    }
    else if(class(CTMM)[1]=="bandwidth") # bandwidth information was precalculated
    {
      KDE[[i]] <- CTMM[[i]]
      axes <- names(KDE[[i]]$h)
    }
    else
    { stop(paste("CTMM argument is of class",class(CTMM)[1])) }

    DEBIAS[[i]] <- ifelse(debias,KDE[[i]]$bias,FALSE)
  } # end loop over individuals

  COL <- length(axes)

  # determine desired (absolute) resolution
  dr <- sapply(1:n,function(i){sqrt(diag(KDE[[i]]$H))/res}) # (axes,individuals)
  dim(dr) <- c(COL,n)
  dr <- apply(dr,1,min)

  # # combine everything for grid calculation
  # H.all <- lapply(1:n,function(i){N=length(data[[i]]$t);K=prepare.H(KDE[[i]]$H,n=N,axes=axes);dim(K)=c(N,COL^2);return(K)})
  # H.all <- do.call("rbind", H.all)
  # dim(H.all) <- c(length(H.all)/COL^2,COL,COL)
  #
  # # take only the necessary columns
  # data.all <- lapply(data,function(d) { d[,c("t",axes)] })
  # data.all <- do.call("rbind", data.all)
  #
  # # complete grid specification
  # EXT <- extent(CTMM,level=1-error) # Gaussian extent (includes uncertainty)
  # grid <- kde.grid(data.all,H=H.all,axes=axes,alpha=error,res=res,dr=dr,grid=grid,EXT.min=EXT)

  # loop over individuals
  for(i in 1:n)
  {
    EXT <- extent(CTMM[[i]],level=1-error) # Gaussian extent (includes uncertainty)
    if(COMPATIBLE) # force grids compatible
    {
      grid$align.to.origin <- TRUE
      grid$dr <- dr
    }
    GRID <- kde.grid(data[[i]],H=KDE[[i]]$H,axes=axes,alpha=error,res=res,dr=dr,grid=grid,EXT.min=EXT) # individual grid

    KDE[[i]] <- c(KDE[[i]],kde(data[[i]],KDE[[i]]$H,axes=axes,bias=DEBIAS[[i]],W=KDE[[i]]$weights,alpha=error,dr=dr,grid=GRID))

    KDE[[i]] <- new.UD(KDE[[i]],info=attr(data[[i]],"info"),type='range',CTMM=ctmm())
    # in case bandwidth is pre-calculated...
    if(class(CTMM0[[i]])[1]=="ctmm") { attr(KDE[[i]],"CTMM") <- CTMM0[[i]] }
  }

  names(KDE) <- names(data)
  if(DROP) { KDE <- KDE[[1]] }
  return(KDE)
}


# if a single H matrix is given, make it into an array of H matrices
# output [n,d,d]
prepare.H <- function(H,n,axes=c('x','y'))
{
  d <- length(axes)

  # one variance given - promote to matrix first - then passes to later stage
  if(length(H)==1)
  { H <- H*diag(d) }
  else if(is.null(dim(H))) # array of variances given
  {
    H <- sapply(H,function(h) h * diag(d),simplify='array') # [d,d,n]
    H <- aperm(H,c(3,1,2))
  }

  # one matrix given
  if(length(dim(H))==2 && all(dim(H)==c(d,d)))
  {
    H <- array(H,c(d,d,n))
    H <- aperm(H,c(3,1,2))
  }

  return(H)
}


##################################
# construct my own KDE objects
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
  H <- prepare.H(H,n,axes=axes)

  # fill in grid information
  grid <- kde.grid(data,H=H,axes=axes,alpha=alpha,res=res,dr=dr,grid=grid)

  R <- grid$r
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

    # constrain to within grid
    i1 <- pmax(i1,1)
    i2 <- pmin(i2,dim(PMF))

    SUB <- lapply(1:length(i1),function(d){ i1[d]:i2[d] })

    # I can't figure out how to do this in one line
    if(length(SUB)==1) # 1D
    { PMF[SUB[[1]]] <- PMF[SUB[[1]]] + W[i]*pnorm1(R[[1]][SUB[[1]]]-r[i,1],H[i,,],dr,alpha) }
    else if(length(SUB)==2) # 2D
    { PMF[SUB[[1]],SUB[[2]]] <- PMF[SUB[[1]],SUB[[2]]] + W[i]*pnorm2(R[[1]][SUB[[1]]]-r[i,1],R[[2]][SUB[[2]]]-r[i,2],H[i,,],dr,alpha) }
    else if(length(SUB)==3) # 3D
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

      # area then volume correction
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
  result <- list(PDF=PMF/dV,CDF=CDF,axes=axes,r=R,dr=dr)
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

  # extrapolation issue with bias<1
  # introduce a NULL cell with 0 probability
  VOL <- c(0,VOL)
  CDF <- c(0,CDF)

  # evaluate the debiased cdf on the original volume grid
  CDF <- stats::approx(x=VOL/bias,y=CDF,xout=VOL,yleft=0,yright=1)$y

  # drop null cell
  CDF <- CDF[-1]

  # fix tiny numerical errors from ?roundoff?
  CDF <- sort(CDF)

  # recalculate pdf
  PMF <- cdf2pmf(CDF)
  # fix inconsistency from yright=1
  CDF <- cumsum(PMF)

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

  PMF <- sort(PMF,decreasing=TRUE,index.return=TRUE)

  IND <- PMF$ix # sorted indices
  PMF <- PMF$x  # sorted values

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
  # should now be positive from re-sort
  # but its not necessarily decreasing like its supposed to, because of small numerical errors
  for(i in 2:length(PMF)) { PMF[i] <- min(PMF[i-1:0]) }

  # # OLD SOLUTION
  # # ties or something is causing numerical errors with this transformation
  # # I don't completely understand the issue, but it seems to follow a pattern
  # DECAY <- diff(PMF)
  # # this quantity should be negative, increasing to zero
  # for(i in 1:length(DECAY))
  # {
  #   # when its positive, its usually adjoining a comparable negative value
  #   if(DECAY[i]>0)
  #   {
  #     # find the adjoining error
  #     j <- which.min(DECAY[i + -1:1]) - 2 + i
  #     # tie out the errors
  #     DECAY[i:j] <- mean(DECAY[i:j])
  #   }
  # }
  # PMF <- -rev(cumsum(c(0,rev(DECAY))))

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
      # what cell(s) is the line segment in?
      z.mid <- mean(z.cross[i:(i+1)])

      x <- sqrt(s)*v[1]*z.mid
      y <- sqrt(s)*v[2]*z.mid

      r <- abs(x-X)
      c <- abs(y-Y)
      # the closest point(s)
      r <- r==min(r)
      c <- c==min(c)

      cdf[r,c] <- (stats::pnorm(z.cross[i+1])-stats::pnorm(z.cross[i]))/(sum(r)*sum(c))
    }
  }
  else if(ZERO==2) # point degeneracy
  {
    r <- abs(X)
    c <- abs(Y)
    # the closest point(s)
    r <- r==min(r)
    c <- c==min(c)

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
# Newton-Cotes integrators (2D)
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
    names(level.UD) <- NAMES.CI[2] # point estimate
    return(level.UD)
  }

  dV <- prod(object$dr)

  # point estimate
  area <- sum(object$CDF <= level.UD) * dV
  names(area) <- NAMES.CI[2] # point estimate

  # chi square approximation of uncertainty
  if(!is.null(object$DOF.area))
  {
    area <- chisq.ci(area,DOF=2*object$DOF.area[1],alpha=1-level)
    names(area) <- NAMES.CI
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

  names(P) <- NAMES.CI
  return(P)
}

#######################
# summarize details of akde object
summary.UD <- function(object,level=0.95,level.UD=0.95,units=TRUE,...)
{
  type <- attr(object,'type')
  if(type!='range') { stop(type," area is not generally meaningful, biologically.") }

  area <- CI.UD(object,level.UD,level)
  if(length(area)==1) { stop("Object is not a range distribution.") }

  # pretty units
  # do we convert base units?
  unit.info <- unit(area[2],"area",SI=!units)
  name <- unit.info$name
  scale <- unit.info$scale

  area <- array(area/scale,c(1,3))
  rownames(area) <- paste("area (",name,")",sep="")

  colnames(area) <- NAMES.CI

  SUM <- list()

  SUM$DOF <- c(object$DOF.area[1],object$DOF.H)
  names(SUM$DOF) <- c("area","bandwidth")

  SUM$CI <- area

  return(SUM)
}
#methods::setMethod("summary",signature(object="UD"), function(object,...) summary.UD(object,...))
