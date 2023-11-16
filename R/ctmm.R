#######################################
# convenience wrapper for new.ctmm
ctmm <- function(tau=NULL,omega=FALSE,isotropic=FALSE,range=TRUE,circle=FALSE,error=FALSE,axes=c("x","y"),...)
{
  tau.names <- c("position","velocity","acceleration")
  List <- list(...)
  List <- List[!sapply(List,is.null)]

  info <- List$info
  List$info <- NULL
  if(is.null(info)) { info=list() }

  List$error <- error
  if(any(error>0)) { List$errors <- TRUE } # Boolean
  List$circle <- circle
  List$axes <- axes

  # put covariance into universal format
  if(length(axes)==1) { isotropic <- TRUE }
  if(!is.null(List$sigma)) { List$sigma <- covm(List$sigma,isotropic=isotropic[1],axes=axes) }
  names(isotropic) <- c("sigma","mu")[1:length(isotropic)]
  List$isotropic <- isotropic

  # label tau elements
  K <- length(tau)
  if(K)
  {
    tau <- sort(tau,decreasing=TRUE)
    names(tau) <- tau.names[1:K]
    # tau <- tau[tau>0]
  }
  if(!length(tau)) { tau <- NULL } # NULL, NA, integer(0), ...
  List$tau <- tau

  List$omega <- as.numeric(omega)
  # requires matching decay timescales
  if(omega) { List$tau <- c(1,1)/mean(1/List$tau) }

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

  if(!range && length(tau)>1 && tau[1]==tau[2])
  { stop("Ballistic motion not yet supported.") }

  if(!range && circle) { stop("Inconsistent model options: range=FALSE, circle=TRUE.") }

  # default mean function
  if(is.null(List$mean)) { List$mean <- "stationary" }

  if(!is.null(List$mu))
  {
    # format simple-style stationary mean properly
    List$mu <- rbind(List$mu)
    colnames(List$mu) <- axes
  }

  # default (output) axes link
  if(is.null(List$link)) { List$link <- "identity" }

  # default (input) time-link
  if(is.null(List$timelink)) { List$timelink <- "identity" }

  # default dynamics model
  if(is.null(List$dynamics)) { List$dynamics <- "stationary" }

  # FIX THIS
  #if(!is.null(List$COV.mu)) { dimnames(List$COV.mu) <- list(axes,axes) }

  # supply default parameters / check sanity / label dimensions
  #
  #

  result <- new.ctmm(List,info=info)

  return(result)
}


ctmm.ctmm <- function(CTMM)
{
  List <- methods::getDataPart(CTMM)
  names(List) <- names(CTMM) # bug workaround
  List$info <- attr(CTMM,"info")
  CTMM <- do.call(ctmm,List)
  return(CTMM)
}


# dynamical states
get.states <- function(CTMM)
{
  S <- CTMM$dynamics
  if(length(S))
  {
    S <- CTMM[[S]]
    S <- names(S)
    S <- unique(S)
  }
  return(S)
}


# compute all tau/omega related quantities that we need
get.taus <- function(CTMM,zeroes=FALSE,simplify=FALSE)
{
  STATES <- get.states(CTMM)
  for(S in STATES) { CTMM[[S]] <- get.taus(CTMM,zeroes=zeroes,simplify=simplify) }

  CTMM$tau <- sort(CTMM$tau,decreasing=TRUE)
  if(simplify) { CTMM$tau <- CTMM$tau[CTMM$tau>0] }
  K <- if(zeroes) { length(CTMM$tau) } else { continuity(CTMM) }
  CTMM$K <- K

  if(!simplify)
  { VARS <- dimnames(CTMM$COV)[[1]] } # variances estimated (can be NULL)
  else
  { VARS <- NULL } # don't care what model was attempted

  # can't use range boolean because of approximations in ctmm.loglike
  if(K==1 && CTMM$tau[1]<Inf) # OU
  {
    CTMM$tau.names <- "tau position"
    CTMM$TAU <- CTMM$tau[1]
  }
  else if(K>1 && CTMM$tau[1]==Inf) # IOU
  {
    CTMM$tau.names <- "tau velocity"
    CTMM$TAU <- CTMM$tau[2]
    CTMM$Omega2 <- 1/CTMM$tau[2] # tau[1]*Omega^2 (everything is renormalized via tau[1] for IOU)
    CTMM$J.Omega2 <- -1/CTMM$tau[2]^2

    CTMM$f <- 1/CTMM$tau
    CTMM$f.nu <- c( mean(CTMM$f) , +diff(CTMM$f)/2 ) # (f,nu)
    CTMM$TfOmega2 <- 2*CTMM$f.nu[1]/CTMM$Omega2 # 2f/Omega^2 with renormalized Omega^2 -> tau[1]*Omega^2
  }
  else if(K>1 && all(c("tau position","tau velocity","omega") %in% VARS)) # population model with tau.p, tau.v, omega
  {
    CTMM$tau.names <- c('tau position','tau velocity','omega')
    CTMM$f <- 1/CTMM$tau # make sure this isn't used outside of this function
    CTMM$f.nu <- c( CTMM$f , CTMM$omega ) # (f,omega)
    CTMM$Omega2 <- prod(CTMM$f) + CTMM$omega^2
    # CTMM$TfOmega2 <- 2*CTMM$f.nu[1]/CTMM$Omega2 # what would this be for here?
    # decay & oscillation period
    CTMM$TAU <- c(1,1,2*pi)/CTMM$f.nu

    # d(f,omega)/d(tau,omega) # needed for SVC CIs
    CTMM$J.nu.tau <- diag(c(-CTMM$f^2,1))

    # d(tau,T)/d(tau,omega) # needed for summary
    CTMM$J.TAU.tau <- diag(c(1,1,-2*pi/CTMM$omega^2))

    CTMM$J.Omega2 <- c(CTMM$f[2:1],2*CTMM$omega) %*% CTMM$J.nu.tau # dOmega^2/d(f,omega) d(f,omega)/d(tau,tau,omega)
  }
  else if(K>1 && (CTMM$tau[1]>CTMM$tau[2] || all(c("tau position","tau velocity") %in% VARS))) # overdamped
  {
    # the canonical parameterization in this regime is (tau1,tau2) as this includes (0,0) and extends to (Inf,Inf) with rectangular boundary conditions
    CTMM$tau.names <- c('tau position','tau velocity')
    CTMM$TAU <- CTMM$tau
    CTMM$f <- 1/CTMM$tau
    CTMM$f.nu <- c( mean(CTMM$f) , +diff(CTMM$f)/2 ) # (f,nu)
    CTMM$Omega2 <- prod(CTMM$f) # Omega^2
    CTMM$TfOmega2 <- 2*CTMM$f.nu[1]/CTMM$Omega2 # 2f/Omega^2 term

    ### Jacobians ###
    # d(f1,f2)/d(tau1,tau2) # needed for SVF CIs
    CTMM$J.f.tau <- -diag(CTMM$f^2)
    # d(tau1,tau2)/d(f1,f2)
    CTMM$J.tau.f <- -diag(CTMM$tau^2)

    # d(f,nu)/d(tau1,tau2) # needed for diff(f)==0 test
    CTMM$J.nu.tau <- rbind( -c(1,1)*CTMM$f^2/2 , -c(-1,+1)*CTMM$f^2/2 )

    # dOmega^2/d(tau1,tau2) # needed for speed CIs
    CTMM$J.Omega2 <- -CTMM$Omega2/CTMM$tau
  }
  else if(K>1 && CTMM$tau[1]==CTMM$tau[2] && !CTMM$omega && "omega"%nin%VARS) # critically damped
  {
    # the canonical parameterization in this regime is (tau) as this includes (0) and extends to (Inf)
    CTMM$tau.names <- c('tau')
    CTMM$TAU <- CTMM$tau[1]
    CTMM$f <- 1/CTMM$tau
    CTMM$f.nu <- c( CTMM$f[1] , 0 )
    CTMM$Omega2 <- prod(CTMM$f)
    CTMM$TfOmega2 <- 2*CTMM$f.nu[1]/CTMM$Omega2

    CTMM$J.tau.f <- -CTMM$tau[1]^2
    CTMM$J.f.tau <- -CTMM$f^2

    CTMM$J.Omega2 <- -2/CTMM$tau[1]^3
  }
  else if(K>1 && (CTMM$omega || "omega"%in%VARS)) # underdamped
  {
    # the canonical parameterization in this regime is utlimately (tau,omega) as this includes (0,0) and extends to (Inf,Inf)
    CTMM$tau.names <- c('tau','omega')
    CTMM$f <- 1/CTMM$tau # make sure this isn't used outside of this function
    CTMM$f.nu <- c( mean(CTMM$f) , CTMM$omega ) # (f,omega)
    CTMM$Omega2 <- sum(CTMM$f.nu^2)
    CTMM$TfOmega2 <- 2*CTMM$f.nu[1]/CTMM$Omega2
    # decay & oscillation period
    CTMM$TAU <- c(1,2*pi)/CTMM$f.nu

    # d(f,omega)/d(tau,omega) # needed for SVC CIs
    CTMM$J.nu.tau <- diag(c(-CTMM$f[1]^2,1))

    # d(tau,T)/d(tau,omega) # needed for summary
    CTMM$J.TAU.tau <- diag(c(1,-2*pi/CTMM$omega^2))

    CTMM$J.Omega2 <- 2*CTMM$f.nu %*% CTMM$J.nu.tau # dOmega^2/d(f,omega) d(f,omega)/d(tau1,tau2)
  }
  else # IID || BM
  { CTMM$tau.names <- NULL }

  return(CTMM)
}


# returns the canonical parameters of a tau vector
pars_tauv <- function(tau,tauc=tau)
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


###########################
ctmm.prepare <- function(data,CTMM,precompute=TRUE,tau=TRUE,DIM=length(CTMM$axes),calibrate=FALSE,verbose=TRUE,...)
{
  axes <- CTMM$axes

  # prepare timescale related info
  if(tau)
  {
    K <- length(CTMM$tau)  # dimension of hidden state per spatial dimension
    # range <- TRUE

    if(K>0)
    {
      # continuity reduction
      CTMM$tau = CTMM$tau[CTMM$tau>0]
      K <- length(CTMM$tau)
    }

    CTMM$K <- K
  }

  # error parameters
  TYPE <- DOP.match(axes)
  if(TYPE!="unknown")
  {
    UERE <- attr(data,"UERE")$UERE[,TYPE]
    names(UERE) <- rownames(attr(data,"UERE")$UERE)
  }
  else
  {
    UERE <- 0
    names(UERE) <- "all"
  }

  if(is.null(names(CTMM$error))) # fix manually specified error parameters
  {
    if(identical(CTMM$error,TRUE)) # lazy set fix
    { CTMM$error <- UERE }
    else # could be 10 meters for multiple classes
    { CTMM$error <- array(CTMM$error,length(UERE)) }
    names(CTMM$error) <- names(UERE)
  } # otherwise leave this alone, as could be a proper subset of data$class levels
  if(any(CTMM$error>0)) { CTMM$errors <- TRUE }

  # evaluate mean function for this data set if no vector is provided or vector can change
  if(precompute && (is.null(CTMM$mean.vec) || length(drift.pars(CTMM))))
  {
    U <- drift.mean(CTMM,data$t,verbose=verbose)
    CTMM$mean.vec <- U

    range <- CTMM$range
    AXES <- length(CTMM$axes)
    if(!range && ncol(U)==1)
    {
      CTMM$UU <- matrix(0,1,1)
      CTMM$REML.loglike <- 0  # extra term for REML likelihood (no stationary contributions here)
    }
    else
    {
      UU <- t(U) %*% U
      if(!range) { UU[1,1] <- 0 } # zero stationary contribution (remove below)
      CTMM$UU <- UU
      CTMM$REML.loglike <- AXES/2*PDlogdet(UU[-1,-1]) # extra term for REML likelihood
    }
  }

  # evaluate error stuff
  if(precompute && (is.null(CTMM$error.mat) || is.null(CTMM$class.mat)))
  {
    CTMM$class.mat <- get.class.mat(data)

    # construct error matrix, if UERE is unknown construct error matrix @ RMS UERE=1 for relative variances
    # ERROR <- CTMM
    # ERROR$error <- as.logical(ERROR$error)
    CTMM$error.mat <- get.error(data,CTMM,DIM=DIM,calibrate=calibrate,...) # store error matrix (modulo UERE if fit)
    # this is more of an error structure matrix (modulo variance)
  } # end precomputed matrices

  return(CTMM)
}

# undo the above
ctmm.repair <- function(CTMM,K=length(CTMM$tau))
{
  # repair dropped zero timescales
  if(K && length(CTMM$tau)) { CTMM$tau <- replace(numeric(K),1:length(CTMM$tau),CTMM$tau) }
  else if(K) { CTMM$tau <- numeric(K) }

  if(FALSE && !CTMM$range) # no longer needed
  {
    K <- length(CTMM$tau) # why do I need this now?
  }
  CTMM$K <- NULL

  # erase evaluated mean vector from ctmm.prepare
  CTMM$class.mat <- NULL
  CTMM$mean.vec <- NULL
  CTMM$UU <- NULL
  CTMM$REML.loglike <- NULL # deleted in ctmm.fit
  CTMM$error.mat <- NULL

  return(CTMM)
}


# degree of continuity in the model
continuity <- function(CTMM)
{
  K <- sum(CTMM$tau > 0)
  return(K)
}
