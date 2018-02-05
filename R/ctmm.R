new.covm <- methods::setClass("covm", representation("matrix",par="numeric",isotropic="logical"))

#######################################
# convenience wrapper for new.ctmm
ctmm <- function(tau=NULL,isotropic=FALSE,range=TRUE,circle=FALSE,error=FALSE,axes=c("x","y"),...)
{
  tau.names <- c("position","velocity","acceleration")
  List <- list(...)
  List <- List[!sapply(List,is.null)]

  info <- List$info
  List$info <- NULL
  if(is.null(info)) { info=list() }

  List$error <- error
  List$circle <- circle
  List$axes <- axes

  # put covariance into universal format
  if(length(axes)==1) { isotropic <- TRUE }
  if(!is.null(List$sigma)) { List$sigma <- covm(List$sigma,isotropic=isotropic,axes=axes) }
  List$isotropic <- isotropic

  # label tau elements
  K <- length(tau)
  CPF <- List$CPF
  if(is.null(CPF)) { CPF <- FALSE }
  if(K && !CPF)
  {
    tau <- sort(tau,decreasing=TRUE)
    names(tau) <- tau.names[1:K]
    # tau <- tau[tau>0]
  }
  if(!length(tau)) { tau <- NULL } # NULL, NA, integer(0), ...
  List$CPF <- CPF
  if(CPF && length(tau)) { names(tau) <- c("period","decay") }
  List$tau <- tau

  # label parameter covariance matrix
  # cov.names <- c("area",names(tau))
  # if(circle) { cov.names <- c(cov.names,"circle") }
  # List$COV can resolve to List$COV.mu apparently
  # if(!is.null(List[["COV"]])) { dimnames(List$COV) <- list(cov.names,cov.names) }

  if(length(tau))
  {
    if(tau[1]==Inf)
    { range <- FALSE }
    else
    { range <- TRUE }
  }
  List$range <- range

  if(!range & circle) { stop("Inconsistent model options: range=FALSE, circle=TRUE.") }

  # default mean function
  if(is.null(List$mean)) { List$mean <- "stationary" }

  if(!is.null(List$mu))
  {
    # format simple-style stationary mean properly
    List$mu <- rbind(List$mu)
    colnames(List$mu) <- axes
  }

  # FIX THIS
  #if(!is.null(List$COV.mu)) { dimnames(List$COV.mu) <- list(axes,axes) }

  # supply default parameters / check sanity / label dimensions
  # drift <- get(List$mean)
  # List <- drift@clean(List)

  result <- new.ctmm(List,info=info)

  return(result)
}
#print.ctmm <- function(x,...) { print.listof(x,...) }

ctmm.ctmm <- function(CTMM)
{
  List <- methods::getDataPart(CTMM)
  names(List) <- names(CTMM) # bug workaround
  List$info <- attr(CTMM,"info")
  CTMM <- do.call(ctmm,List)
  return(CTMM)
}

# 2D covariance matrix universal format
covm <- function(pars,isotropic=FALSE,axes=c("x","y"))
{
  if(is.null(pars))
  { return(NULL) }
  else if(class(pars)=="covm")
  { return(pars) }

  if(length(axes)==1)
  {
    sigma <- array(pars,c(1,1))
    pars <- sigma[1,1]

    names(pars) <- c("area")
  }
  else if(length(axes)==2)
  {
    if(length(pars)==1)
    {
      pars <- c(pars,0,0)
      sigma <- diag(pars[1],2)
    }
    else if(length(pars)==3)
    { sigma <- sigma.construct(pars) }
    else if(length(pars)==4)
    {
      sigma <- pars
      pars <- sigma.destruct(sigma,isotropic=isotropic)
    }

    # isotropic error check
    if(isotropic)
    {
      pars <- c(mean(diag(sigma)),0,0)
      sigma <- diag(pars[1],2)
    }

    names(pars) <- c("area","eccentricity","angle")
  }

  dimnames(sigma) <- list(axes,axes)

  new.covm(sigma,par=pars,isotropic=isotropic)
}

# construct covariance matrix from 1-3 parameters
sigma.construct <- function(pars)
{
  GM <- pars[1]
  if(length(pars)==1)
  {
    e <- 0
    theta <- 0
  }
  else
  {
    e <- pars[2]
    theta <- pars[3]
  }

  u <- c(cos(theta),sin(theta))
  v <- c(-sin(theta),cos(theta))
  e <- exp(e/2)

  sigma <- GM * ( (u%o%u)*e + (v%o%v)/e )

  return(sigma)
}

# reduce covariance matrix to 1-3 parameters
sigma.destruct <- function(sigma,isotropic=FALSE) # last arg not implemented
{
  stuff <- eigen(sigma)

  e <- stuff$values
  GM <- sqrt(prod(e))
  e <- log(e[1]/e[2])

  if(e==0)
  { theta <- 0 }
  else
  {
    theta <- stuff$vectors[,1]
    theta <- atan(theta[2]/theta[1])
  }

  return(c(GM,e,theta))
}

# return the COV matrix for covm par representation
COV.covm <- function(sigma,n,k=1)
{
  isotropic <- sigma@isotropic
  par <- sigma@par
  sigma <- methods::getDataPart(sigma)

  A <- par["area"]
  ecc <- par["eccentricity"]
  theta <- par["angle"]

  DOF.mu <- n
  COV.mu <- sigma/n

  if(isotropic || length(sigma)==1)
  {
    COV <- cbind( A^2/n )
    dimnames(COV) <- list("area","area")
  }
  else
  {
    # orient eccentricity
    ecc <- sign(sigma[1,1]-sigma[2,2]) * ecc

    # covariance matrix for c( sigma_xx , sigma_yy , sigma_xy )
    n <- n-k
    COV <- diag(0,3)
    COV[1:2,1:2] <- 2/n * sigma^2
    COV[3,] <- c( 2*sigma[1,1]*sigma[1,2] , 2*sigma[2,2]*sigma[1,2] , sigma[1,1]*sigma[2,2]+sigma[1,2]^2 )/n
    COV[,3] <- COV[3,]

    # gradient matrix d sigma / d par
    grad <- c( sigma[1,1] , sigma[2,2] , sigma[1,2] )/A
    grad <- cbind( grad, A/2*c( sinh(ecc/2)+cosh(ecc/2)*cos(2*theta) , sinh(ecc/2)-cosh(ecc/2)*cos(2*theta) , cosh(ecc/2)*sin(2*theta) ) )
    grad <- cbind( grad, c( -2*sigma[1,2] , 2*sigma[1,2] , A*sinh(ecc/2)*2*cos(2*theta) ) )

    # gradient matrix d par / d sigma via inverse function theorem
    grad <- solve(grad)

    COV <- (grad) %*% COV %*% t(grad)
    COV <- He(COV) # asymmetric errors
    dimnames(COV) <- list(names(par),names(par))
  }

  return(list(COV=COV,COV.mu=COV.mu,DOF.mu=DOF.mu))
}


# blank generic function
pars <- function(...) { UseMethod("pars") }

# return the canonical parameters of a covariance matrix
pars.covm <- function(COVM)
{
  if(COVM@isotropic)
  { return(COVM@par[1]) }
  else
  { return(COVM@par) }
}

# returns the canonical parameters of a tau vector
pars.tauv <- function(tau,tauc=tau)
{
  if(length(tauc)==0)
  { return(NULL) }
  else if(tauc[1] < Inf)
  { return(tau[1:length(tauc)]) }
  else if(length(tauc)==1)
  { return(NULL) }
  else
  { return(tau[-1]) }
}

############################
# coarce infinite parameters into finite parameters appropriate for numerics
###########################
ctmm.prepare <- function(data,CTMM,REML=FALSE,precompute=TRUE)
{
  K <- length(CTMM$tau)  # dimension of hidden state per spatial dimension
  axes <- CTMM$axes
  range <- TRUE

  if(K>0)
  {
    # numerical limit for rangeless processes
    if(CTMM$tau[1]==Inf)
    {
      range <- FALSE
      # aim for OU/OUF decay that is half way to machine epsilon
      CTMM$tau[1] <- log(2^((.Machine$double.digits-1)/2))*(last(data$t)-data$t[1])
      # CTMM$tau[1] <- (last(data$t)-data$t[1]) / (.Machine$double.eps)^(1/4)

      # diffusion -> variance
      if(!is.null(CTMM$sigma))
      {
        TAU <- CTMM$tau[1]
        CTMM$sigma <- CTMM$sigma*TAU
        CTMM$sigma@par[1] <- CTMM$sigma@par[1]*TAU
      }
    }

    # continuity reduction
    CTMM$tau = CTMM$tau[CTMM$tau!=0]
    K <- length(CTMM$tau)
  }
  # I am this lazy
  # if(K==0) { K <- 1 ; CTMM$tau <- 0 }

  CTMM$range <- range

  # evaluate mean function for this data set if no vector is provided
  if(precompute && is.null(CTMM$mean.vec))
  {
    drift <- get(CTMM$mean)
    U <- drift(data$t,CTMM)
    CTMM$mean.vec <- U

    if(REML)
    {
      UU <- t(U) %*% U
      CTMM$REML.loglike <- (length(axes)/2)*(log(det(UU)) + ncol(U)*log(2*pi))
    }
  }

  return(CTMM)
}

# undo the above
ctmm.repair <- function(CTMM,K=length(CTMM$tau))
{
  # repair dropped zero timescales
  if(K && length(CTMM$tau)) { CTMM$tau <- replace(numeric(K),1:length(CTMM$tau),CTMM$tau) }
  else if(K) { CTMM$tau <- numeric(K) }

  if(!CTMM$range)
  {
    K <- length(CTMM$tau)

    # variance -> diffusion
    TAU <- CTMM$tau[1]
    CTMM$sigma <- CTMM$sigma/TAU
    CTMM$sigma@par[1] <- CTMM$sigma@par[1]/TAU
    CTMM$tau[1] <- Inf

    ## delete garbate estimates
    # CTMM$mu <- NULL
    # CTMM$COV.mu <- NULL
    # CTMM$DOF.mu <- NULL
  }

  # erase evaluated mean vector from ctmm.prepare
  CTMM$mean.vec <- NULL

  return(CTMM)
}

## prepare error array, also return a flag #
# 0 : no error
# 1 : constant error parameter fit
# 2 : proportional error parameter fit to DOP value
# 3 : full error no fit
get.error <- function(DATA,CTMM,flag=FALSE)
{
  n <- length(DATA$t)
  axes <- CTMM$axes
  COLS <- names(DATA)

  if(CTMM$error)
  {
    TYPE <- DOP.match(axes)
    # DOP.LIST is global variable from uere.R
    TYPE <- DOP.LIST[[TYPE]]
    AXES <- TYPE$axes
    ERE <- TYPE$ERE
    DOP <- TYPE$DOP

    if(ERE %in% COLS) # calibrated errors - HERE
    {
      error <- DATA[[ERE]]^2/length(axes)
      FLAG <- 3
    }
    else if(DOP %in% COLS) # fitted errors - HDOP
    {
      error <- (CTMM$error*DATA[[DOP]])^2/length(axes)
      FLAG <- 2
    }
    else # fitted errors - no HDOP
    {
      error <- rep(CTMM$error^2/length(axes),n)
      FLAG <- 1
    }
  } # END error
  else # no error
  {
    FLAG <- 0
    error <- rep(0,n)
  }

  if(flag) { return(FLAG) }
  else
  {
    attr(error,"flag") <- FLAG
    return(error)
  }
}

# degree of continuity in the model
continuity <- function(CTMM)
{
  K <- sum(CTMM$tau > 0)
  return(K)
}
