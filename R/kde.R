# bias of Gaussian Reference function AKDE
# generalize to non-stationary mean
akde.bias <- function(CTMM,H,lag,DOF,weights)
{
  axes <- CTMM$axes
  sigma <- methods::getDataPart(CTMM$sigma)

  # weighted correlation
  ACF <- Vectorize( svf.func(CTMM,moment=FALSE)$ACF )
  if(is.null(DOF)) # exact version
  {
    MM <- lag # copy structure
    MM[] <- ACF(lag) # preserve structure... why do I have to do it this way?
    MM <- c(weights %*% MM %*% weights)
  }
  else # fast version
  { MM <- sum(DOF*ACF(lag)) }
  VAR <- 1-MM # sample variance (relative)
  COV <- cbind(VAR*sigma) # sample covariance
  dimnames(COV) <- list(axes,axes)

  # variance inflation factor
  bias <- ( det(COV+H)/det(cbind(sigma)) )^(1/length(axes))
  # remove cbind if/when I can get det.numeric working

  # name dimensions of bias
  bias <- rep(bias,length(axes))
  names(bias) <- axes

  R <- list(bias=bias,COV=COV)
  return(R)
}


# population AKDE
pkde <- function(data,UD,kernel="individual",weights=FALSE,ref="Gaussian",...)
{
  CLASS <- class(UD[[1]])[1]
  if(CLASS!="UD") { stop("UD class is ",CLASS) }
  akde(data,CTMM=UD,kernel=kernel,weights=weights,ref=ref,...)
}


# AKDE single or list
# (C) C.H. Fleming (2016-2022)
# (C) Kevin Winner & C.H. Fleming (2016)
akde <- function(data,CTMM,VMM=NULL,R=list(),SP=NULL,SP.in=TRUE,variable="utilization",debias=TRUE,weights=FALSE,smooth=TRUE,error=0.001,res=10,grid=NULL,...)
{
  if(variable!="utilization")
  { stop("variable=",variable," not yet supported by akde(). See npr() or revisitation().") }

  if(length(projection(data))>1) { stop("Data not in single coordinate system.") }
  validate.grid(data,grid)

  # if called by pakde
  if(class(CTMM)[1]=="list" && class(CTMM[[1]])[1]=="UD")
  {
    UD <- CTMM
    # extract models
    CTMM <- lapply(UD,function(ud){ud@CTMM})
  }
  else
  { UD <- NULL }

  DROP <- class(data)[1] != "list"
  data <- listify(data)
  CTMM <- listify(CTMM)
  VMM <- listify(VMM)

  if(length(data)!=length(CTMM))
  { stop("length(data)==",length(data),", but length(CTMM)==",length(CTMM)) }

  DOF <- sapply(CTMM,DOF.area)
  SUB <- DOF<error
  if(any(SUB))
  {
    warning("Fit object returned. DOF[area] = ",paste(DOF[SUB],collapse="; "))
    SUB <- !SUB
    if(any(SUB)) { CTMM[SUB] <- akde(data[SUB],CTMM[SUB],VMM=VMM[SUB],R=R,SP=SP,SP.in=SP.in,variable=variable,debias=debias,weights=weights,smooth=smooth,error=error,res=res,grid=grid,...) }
    return(CTMM)
  }

  # force grids to be compatible
  COMPATIBLE <- length(data)>1 && !is.grid.complete(grid)

  axes <- CTMM[[1]]$axes
  AXES <- length(axes)

  n <- length(data)
  if(class(weights)[1]!="list")
  {
    # assume numerical or boolean / by individual or by time
    if(length(weights)==1 || length(weights)==n) # by individual
    { weights <- as.list(array(weights,n)) }
    else # by time
    { weights <- list(weights) }
  }

  # loop over individuals for bandwidth optimization
  CTMM0 <- VMM0 <- list()
  KDE <- list()
  DEBIAS <- list()
  for(i in 1:n)
  {
    CTMM0[[i]] <- CTMM[[i]] # original model fit
    if(class(CTMM[[i]])[1]=="ctmm") # calculate bandwidth etc.
    {
      axes <- CTMM[[i]]$axes

      # smooth out errors (which also removes duplicate times)
      if(!is.null(VMM[[i]]))
      {
        axis <- VMM[[i]]$axes
        if(any(VMM[[i]]$error>0) && smooth)
        { z <- predict(VMM[[i]],data=data[[i]],t=data[[i]]$t)[[axis]] } # [,axis] fails?
        else
        { z <- data[[i]][[axis]] }
        axes <- c(axes,axis)
      }
      if(any(CTMM[[i]]$error>0) && smooth)
      {
        data[[i]] <- predict(CTMM[[i]],data=data[[i]],t=data[[i]]$t)
        CTMM[[i]]$error <- FALSE # smoothed error model (approximate)
      }
      # copy back !!! times and locations
      if(!is.null(VMM[[i]]))
      {
        data[[i]][,axis] <- z
        VMM0[[i]] <- VMM[[i]] # original model fit
        VMM[[i]]$error <- FALSE # smoothed error model (approximate)
      }

      # calculate individual optimal bandwidth and some other information
      if(is.null(UD)) { KDE[[i]] <- bandwidth(data=data[[i]],CTMM=CTMM[[i]],VMM=VMM[[i]],weights=weights[[i]],verbose=TRUE,error=error,...) }
    }
    else if(class(CTMM)[1]=="bandwidth") # bandwidth information was precalculated
    {
      KDE[[i]] <- CTMM[[i]]
      axes <- names(KDE[[i]]$h)
    }
    else
    { stop(paste("CTMM argument is of class",class(CTMM)[1])) }

    if(is.null(UD))
    {
      if(debias)
      { DEBIAS[[i]] <- KDE[[i]]$bias }
      else
      { DEBIAS[[i]] <- FALSE }
    }
  } # end loop over individuals

  grid <- format_grid(grid,axes=axes)
  COL <- length(axes)

  if(!is.null(UD)) # format population KDE as an individual # then set resolution
  {
    i <- 1 # one population
    # weights and bandwidth
    KDE[[i]] <- bandwidth.pop(data,UD,weights=weights,...)
    kernel <- list(...)$kernel
    if(is.null(kernel)) { kernel <- "individual" }

    DOF <- DOF.area(KDE[[i]]$CTMM)
    if(DOF<error)
    {
      warning("Population fit object returned. DOF[area] = ",DOF)
      return(KDE[[i]]$CTMM)
    }

    # bias information
    if(debias)
    { DEBIAS[[i]] <- KDE[[i]]$bias }
    else
    { DEBIAS[[i]] <- FALSE }
    # bandwidth [d,d,n]
    H <- KDE[[i]]$h^2

    # assemble bandwidths
    KDE[[i]]$H <- list()
    for(j in 1:n)
    {
      if(kernel=="individual")
      { KDE[[i]]$H[[j]] <- prepare.H(H*CTMM[[j]]$sigma,nrow(data[[j]])) }
      else if(kernel=="population")
      { KDE[[i]]$H[[j]] <- prepare.H(H*KDE[[i]]$CTMM$sigma,nrow(data[[j]])) }
      DIM <- dim(KDE[[i]]$H[[j]])
      dim(KDE[[i]]$H[[j]]) <- c(DIM[1],AXES^2)
    }
    KDE[[i]]$H <- do.call(rbind,KDE[[i]]$H)
    DIM <- dim(KDE[[i]]$H)
    dim(KDE[[i]]$H) <- c(DIM[1],AXES,AXES)

    # assemble weights # weights * weights
    weights <- KDE[[i]]$weights
    KDE[[i]]$weights <- lapply(1:n,function(j){KDE[[i]]$weights[j]*UD[[j]]$weights})
    KDE[[i]]$weights <- unlist(KDE[[i]]$weights)

    # assemble data
    data <- lapply(data,function(d){d[,c('t',axes)]})
    data <- do.call(rbind,data)
    data <- list(data)
    dr <- sapply(1:n,function(j){sqrt(min(1,H)*diag(CTMM[[j]]$sigma))/res})
    dim(dr) <- c(COL,n)

    # population GRF
    CTMM0 <- CTMM <- list(KDE[[i]]$CTMM)

    n <- 1 # one population
    DROP <- TRUE
  }
  else
  {
    # determine desired (absolute) resolution from smallest of all individual bandwidths
    dr <- sapply(1:n,function(i){sqrt(pmin(diag(KDE[[i]]$H),diag(CTMM0[[i]]$sigma)))/res}) # (axes,individuals)
    dim(dr) <- c(COL,n)
  }
  dr <- apply(dr,1,grid$dr.fn) # minimum by default

  if(COMPATIBLE) # force grids compatible
  {
    grid$align.to.origin <- TRUE
    if("dr" %nin% names(grid)) { grid$dr <- dr }
  }

  # loop over individuals
  for(i in 1:n)
  {
    if(is.null(VMM))
    { EXT <- CTMM[[i]] }
    else
    { EXT <- list(horizontal=CTMM[[i]],vertical=VMM[[i]]) }
    EXT <- extent(EXT,level=1-error)[,axes,drop=FALSE] # Gaussian extent (includes uncertainty)
    GRID <- kde.grid(data[[i]],H=KDE[[i]]$H,axes=axes,alpha=error,res=res,dr=dr,grid=grid,EXT.min=EXT) # individual grid

    KDE[[i]] <- c(KDE[[i]],kde(data[[i]],H=KDE[[i]]$H,axes=axes,CTMM=CTMM0[[i]],SP=SP,SP.in=SP.in,RASTER=R,bias=DEBIAS[[i]],W=KDE[[i]]$weights,alpha=error,dr=dr,grid=GRID,...))

    if(!is.null(UD))
    {
      # KDE[[i]]$H <- KDE[[i]]$h^2
      KDE[[i]]$H <- NULL
      KDE[[i]]$weights <- weights
    }

    KDE[[i]] <- new.UD(KDE[[i]],info=attr(data[[i]],"info"),type='range',variable="utilization",CTMM=ctmm())
    # in case bandwidth is pre-calculated...
    if(class(CTMM0[[i]])[1]=="ctmm") { attr(KDE[[i]],"CTMM") <- CTMM0[[i]] }
    if(!is.null(VMM)) { KDE[[i]]$VMM <- VMM0[[i]] }
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
  { H <- c(H)*diag(d) }
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
kde <- function(data,H,axes=c("x","y"),CTMM=list(),SP=NULL,SP.in=TRUE,RASTER=list(),bias=FALSE,W=NULL,alpha=0.001,res=NULL,dr=NULL,grid=NULL,variable=NA,normalize=TRUE,trace=FALSE,grad=FALSE,...)
{
  DIM <- length(axes)

  if(!is.na(variable))
  {
    if(variable %in% c("revisitation"))
    { VARIABLE <- "speed" }
    else # data column better have the name given
    { VARIABLE <- variable }
  }

  # format bandwidth matrix
  H <- prepare.H(H,nrow(data),axes=axes)

  # normalize weights
  if(is.null(W))
  { W <- rep(1,length(data$x)) }
  else # unweighted times can be skipped
  {
    SUB <- (W/max(W))>.Machine$double.eps
    W <- W[SUB]
    data <- data[SUB,]
    H <- H[SUB,,,drop=FALSE]
  }
  if(normalize) { W <- W/sum(W) }

  n <- nrow(data)
  r <- get.telemetry(data,axes)

  # fill in grid information
  grid <- kde.grid(data,H=H,axes=axes,alpha=alpha,res=res,dr=dr,grid=grid)

  R <- grid$r
  # generalize this for future grid option use
  dH <- grid$dH
  dr <- grid$dr
  dV <- prod(dr)

  R0 <- sapply(R,first)
  # corner origin to minimize arithmetic later
  #r <- t(t(r) - R0)
  #R <- lapply(1:length(R0),function(i){ R[[i]] - R0[i] })

  # stationary versus non-stationary suitability
  if(length(RASTER))
  {
    STUFF <- expand.factors(RASTER,CTMM$formula,fixed=TRUE)
    RASTER <- STUFF$R
    data <- STUFF$data

    proj <- CTMM@info$projection
    # calculate RASTERs on spatial grid
    RASTER <- lapply(RASTER,function(r){R.grid(grid$r,proj=proj,r)})
    # this needs to be moved up for multiple individuals?

    STATIONARY <- is.stationary(CTMM,RASTER)

    # calculate one suitability layer for all times
    if(STATIONARY) { RASTER <- R.suit(RASTER,CTMM) }
    # otherwise we calculate one suitability per time/kernel
  }

  if(length(SP) && DIM==2)
  {
    proj <- CTMM@info$projection
    # this is inaccurate
    # SP <- sp::spTransform(SP,proj)
    if(!any(grepl('sf',class(SP)))) { SP <- sf::st_as_sf(SP) }
    SP <- sf::st_transform(SP,crs=sf::st_crs(proj))
    # sf methods are slow
    SP <- sf::as_Spatial(SP)

    # create raster template
    dx <- grid$dr[1]
    dy <- grid$dr[2]

    xmn <- grid$r$x[1]-dx/2
    xmx <- last(grid$r$x)+dx/2

    ymn <- grid$r$y[1]-dy/2
    ymx <- last(grid$r$y)+dy/2

    RSP <- matrix(0,length(grid$r$y),length(grid$r$x))
    RSP <- raster::raster(RSP,xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx,crs=proj)

    # rasterize SP -> RSP
    RSP <- raster::rasterize(SP,RSP,background=NA)
    RSP <- raster::as.matrix(RSP)
    RSP <- t(RSP)[,dim(RSP)[1]:1]
    RSP <- !is.na(RSP)

    if(!SP.in) { RSP <- !RSP }

    # incorporate into RASTER
    if(length(RASTER))
    {
      if(STATIONARY)
      { RASTER <- nant(RASTER,0) * RSP }
      else
      { RASTER <- lapply(RASTER,function(r){nant(r,0) * RSP}) }
    }
    else
    {
      RASTER <- RSP
      STATIONARY <- TRUE
    }
  }
  else
  { RSP <- NULL }

  # probability mass function
  PMF <- array(0,sapply(R,length))
  # Nadaraya-Watson numerator (for regressions)
  if(!is.na(variable)) { NUM <- PMF }
  # gradient information
  if(grad)
  {
    GRAD <- array(0,c(dim(PMF),2))
    HESS <- array(0,c(dim(PMF),2,2))
  }

  SUB <- lapply(1:length(dim(PMF)),function(d){1:dim(PMF)[d]})

  if(trace) { pb <- utils::txtProgressBar(style=3) } # time loops
  for(i in 1:n)
  {
    # sub-grid lower/upper bound indices
    i1 <- floor((r[i,]-dH[i,]-R0)/dr) + 1
    i2 <- ceiling((r[i,]+dH[i,]-R0)/dr) + 1

    # constrain to within grid
    i1 <- pmax(i1,1)
    i2 <- pmin(i2,dim(PMF))

    CHECK <- i2>=i1
    if(any(!CHECK)) { stop("Grid incompatible with data.") }

    SUB <- lapply(1:length(i1),function(d){ i1[d]:i2[d] })

    # I can't figure out how to do these cases in one line?
    if(DIM==1) # 1D
    {
      dPMF <- pnorm1(R[[1]][SUB[[1]]]-r[i,1],H[i,,],dr,alpha)

      # assume grid has exact bounds
      if(length(SP)) { dPMF <- dPMF / sum(dPMF) }

      PMF[SUB[[1]]] <- PMF[SUB[[1]]] + W[i]*dPMF
    }
    else if(DIM==2) # 2D
    {
      # extract suitability sub-grid
      if(length(RASTER))
      {
        if(STATIONARY) # one suitability raster for all kernels
        { R.SUB <- RASTER[ SUB[[1]], SUB[[2]] ] }
        else
        {
          # time interpolate if necessary
          R.SUB <- list()
          for(j in 1:length(RASTER))
          {
            if(length(dim(RASTER))==2)
            { R.SUB[[j]] <- RASTER[[j]][ SUB[[1]], SUB[[2]] ] }
            else
            {
              T.INT <- (data$t[i]-RASTER[[j]]$Z[1])/RASTER[[j]]$dZ + 1
              T.SUB <- c(floor(T.INT),ceiling(T.INT))
              if(T.SUB[1]==T.SUB[2])
              { R.SUB[[j]] <- RASTER[[j]][ SUB[[1]], SUB[[2]], T.SUB[1] ] }
              else
              {
                T.CON <- c(T.SUB[2]-T.INT,T.INT-T.SUB[1])
                R.SUB[[j]] <- T.CON[1]*RASTER[[j]][ SUB[[1]], SUB[[2]], T.SUB[1] ] + T.CON[2]*RASTER[[j]][ SUB[[1]], SUB[[2]], T.SUB[2] ]
              }
            } # end time interpolation
          } # end raster stack preparation
          names(R.SUB) <- names(RASTER)
          R.SUB <- R.suit(R.SUB,CTMM,data[i,])
        } # end non-stationary suitability calculation
      } # end suitability modulus

      X <- R[[1]][SUB[[1]]]-r[i,1]
      Y <- R[[2]][SUB[[2]]]-r[i,2]
      dPMF <- pnorm2(X,Y,H[i,,],dr,alpha)

      # apply suitability and re-normalize
      if(length(RASTER))
      {
        SUM <- sum(dPMF) # preserve truncation error
        dPMF <- dPMF * R.SUB
        dPMF <- nant(SUM/sum(dPMF),0) * dPMF
      }

      dPMF <- W[i]*dPMF

      PMF[SUB[[1]],SUB[[2]]] <- PMF[SUB[[1]],SUB[[2]]] + dPMF

      if(!is.na(variable))
      { NUM[SUB[[1]],SUB[[2]]] <- NUM[SUB[[1]],SUB[[2]]] + data[[VARIABLE]][i]*dPMF }

      if(grad)
      {
        iH <- PDsolve(H[i,,])
        DIM <- c(length(X),length(Y))
        G <- array(0,c(DIM,2))
        G[,,1] <- array(X,DIM)
        G[,,2] <- t( array(Y,rev(DIM)) )
        G <- G %.% iH

        dG <- G
        dG[,,1] <- dPMF * dG[,,1]
        dG[,,2] <- dPMF * dG[,,2]

        GRAD[SUB[[1]],SUB[[2]],] <- GRAD[SUB[[1]],SUB[[2]],] - dG

        GG <- array(G,c(DIM,2,2))
        GG[,,1,1] <- dPMF * G[,,1] * G[,,1]
        GG[,,1,2] <- dPMF * G[,,1] * G[,,2]
        GG[,,2,1] <- dPMF * G[,,2] * G[,,1]
        GG[,,2,2] <- dPMF * G[,,2] * G[,,2]

        HESS[SUB[[1]],SUB[[2]],,] <- HESS[SUB[[1]],SUB[[2]],,] - (dPMF %o% iH) + GG
      }
    }
    else if(DIM==3) # 3D
    { PMF[SUB[[1]],SUB[[2]],SUB[[3]]] <- PMF[SUB[[1]],SUB[[2]],SUB[[3]]] + W[i]*pnorm3(R[[1]][SUB[[1]]]-r[i,1],R[[2]][SUB[[2]]]-r[i,2],R[[3]][SUB[[3]]]-r[i,3],H[i,,],dr,alpha) }

    if(trace) { utils::setTxtProgressBar(pb,i/n) }
  } # end time loop

  if(!is.na(variable)) # revisitation is treated separately
  { return(NUM/PMF) } # E[variable|data]

  if(sum(bias)) # debias area/volume
  {
    if(DIM<=2) # AREA or WIDTH debias
    {
      # debias the area (2D) or width (1D)
      PMF <- debias.volume(PMF,bias=min(bias))
      CDF <- PMF$CDF
      PMF <- PMF$PMF
    }
    else if(DIM==3) # VOLUME debias
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

  result <- list(PDF=PMF/dV,CDF=CDF,axes=axes,r=R,dr=dr)
  if(grad) { result$grad <- GRAD; result$hess <- HESS }
  if(trace) { close(pb) }
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
  # !!! EXPAND INTO DIFFERENT RESOLUTION INTEGRATORS

  cdf <- prod(dr) * Gauss3(X,Y,Z,sigma)

  return(cdf)
}

# UNFINISHED
pnorm1 <- function(X,sigma,dr,alpha=0.001)
{
  # !!! EXPAND INTO DIFFERENT RESOLUTION INTEGRATORS

  cdf <- dr * Gauss1(X,sigma)

  return(cdf)
}


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
# 2D Gaussian pdf
Gauss <- function(X,Y,sigma=NULL,sigma.inv=solve(sigma),sigma.GM=sqrt(det(sigma)))
{
  cdf <- outer(X^2*sigma.inv[1,1],Y^2*sigma.inv[2,2],"+")/2
  cdf <- cdf + (X %o% Y)*sigma.inv[1,2]
  cdf <- exp(-cdf)/(2*pi*sigma.GM)
  return(cdf)
}


#####################
# 3D Gaussian pdf
# assumes uncorrelated z-axis
Gauss3 <- function(X,Y,Z,sigma=NULL,sigma.inv=solve(sigma[1:2,1:2]),sigma.GM=sqrt(det(sigma[1:2,1:2])))
{
  cdf <- Gauss(X,Y,sigma=sigma[1:2,1:2],sigma.inv=sigma.inv,sigma.GM=sigma.GM)
  cdf <- cdf %o% Gauss1(Z,sigma[3,3])

  return(cdf)
}


#####################
# 1D Gaussian pdf
Gauss1 <- function(X,sigma=NULL)
{ exp(-X^2/(2*sigma))/sqrt(2*pi*sigma) }


#####################
# AKDE CIs
CI.UD <- function(object,level.UD=0.95,level=0.95,P=FALSE,convex=FALSE)
{
  if(is.null(object$DOF.area) && P)
  {
    names(level.UD) <- NAMES.CI[2] # point estimate
    return(level.UD)
  }

  # reverse interpolation for CDF
  interpolate <- function(y,val)
  {
    x <- last(which(y < val))
    if(!length(x))
    { x <- 0 }
    else if(x==length(y))
    { x <- length(y) }
    else # interpolate
    {
      Y <- y[x + 0:1] - val
      beta <- diff(Y)
      # 0 == Y[1] + beta * dx
      dx <- -Y[1]/beta
      x <- x + dx
    }
    return(x)
  }

  SORT <- sort(object$CDF,method="quick")
  dV <- prod(object$dr)

  if(convex)
  {
    # warning("sp polygon areas are unreliable.")
    SP <- convex(object,level=level.UD,convex=convex)
    # This estimate seems to be garbage in many cases
    area <- sum( sapply(SP@polygons, function(POLY){POLY@area}) )
  }
  else
  {
    # simple estimate
    # area <- sum(object$CDF <= level.UD) * dV

    # better estimate
    area <- interpolate(SORT,level.UD) * dV
  }
  names(area) <- NAMES.CI[2] # point estimate

  # chi square approximation of uncertainty
  if(!is.null(object$DOF.area)) { area <- chisq.ci(area,DOF=2*object$DOF.area[1],alpha=1-level) }

  if(!P) { return(area) }

  # probabilities associated with these areas
  IND <- area / dV # fractional index
  P <- vint(SORT,IND) # interpolated probability

  # fix lower bound
  P <- pmax(P,0)

  # recorrect point estimate level
  P[2] <- level.UD

  # fix upper bound to not overflow
  P <- pmin(P,1)
  if(P[3]==1) { warning("Outer contour extends beyond raster.") }

  names(P) <- NAMES.CI
  return(P)
}

#######################
# summarize details of akde object
summary.UD <- function(object,convex=FALSE,level=0.95,level.UD=0.95,units=TRUE,...)
{
  type <- attr(object,'type')
  if(type %nin% c('range','revisitation')) { stop(type," area is not generally meaningful, biologically.") }

  area <- CI.UD(object,level.UD,level,convex=convex)
  if(length(area)==1) { stop("Object is not a range distribution.") }

  DOF <- c(object$DOF.area[1],object$DOF.H)
  R <- summary_UD_format(area,DOF=DOF,units=units)
  return(R)
}
#methods::setMethod("summary",signature(object="UD"), function(object,...) summary.UD(object,...))

summary_UD_format <- function(CI,DOF,units=TRUE)
{
  # pretty units
  # do we convert base units?
  unit.info <- unit(CI[2],"area",SI=!units)
  name <- unit.info$name
  scale <- unit.info$scale

  CI <- array(CI/scale,c(1,3))
  rownames(CI) <- paste("area (",name,")",sep="")

  colnames(CI) <- NAMES.CI

  SUM <- list()

  SUM$DOF <- DOF
  if(length(SUM$DOF)==1) { SUM$DOF[2] <- NA }
  names(SUM$DOF) <- c("area","bandwidth")

  SUM$CI <- CI
  class(SUM) <- "area"

  return(SUM)
}


extract <- function(r,UD,DF="CDF",...)
{
  if(length(dim(r)))
  { r <- t(r) }
  else
  { r <- cbind(r) }

  # continuous index
  r[1,] <- (r[1,]-UD$r$x[1])/UD$dr['x'] + 1
  r[2,] <- (r[2,]-UD$r$y[1])/UD$dr['y'] + 1

  bint(UD[[DF]],r)
}
