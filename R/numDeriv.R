###############
# calculate first and second derivatives efficiently
# assumes to start with approximant of maximum (par) and inverse hessian (covariance)
####################
genD <- function(par,fn,zero=FALSE,lower=-Inf,upper=Inf,step=NULL,precision=1/2,parscale=NULL,mc.cores=detectCores(),Richardson=2,order=2,drop=TRUE,control=list(),...)
{
  NAMES <- names(par)

  if(is.null(parscale)) { parscale <- pmin(abs(par),abs(par-lower),abs(upper-par)) }
  if(any(parscale==0)) { parscale[parscale==0] <- 1 }

  if(is.null(step)) { step <- sqrt(2*.Machine$double.eps^precision) }

  DEBUG <- control$DEBUG
  if(is.null(DEBUG)) { DEBUG <- FALSE }

  # convert to natural units - prevents confusion with zero tolerance (bacterium)
  par <- par/parscale
  func <- function(par)
  {
    RETURN <- try(fn(parscale*par))
    if(class(RETURN)[1]!="numeric")
    {
      # debugging bad likelihoods
      if(DEBUG)
      {
        debug(fn)
        try(fn(parscale*par))
        undebug(fn)
      }
      RETURN <- Inf
    }
    return(RETURN)
  }

  ### nuDeriv specific code ###

  d <- 0.01
  zero.tol <- sqrt(.Machine$double.eps/7e-7)
  n <- length(par) # input dimension

  # calculate hessian and gradient simultaneously
  # output dimension is assumed 1
  if(order==2)
  {
    method.args <- list(eps=1e-4,d=d,zero.tol=zero.tol,r=Richardson,v=2,show.details=FALSE)
    STUFF <- numDeriv::genD(func,par,method.args=method.args,...)
    D <- STUFF$D

    grad <- D[1:n]
    D <- D[-(1:n)]

    # Bates and Watts ordering is not R ordering
    hess <- lower.tri(diag(n),diag=TRUE)
    SUM <- 0
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        SUM <- SUM + hess[i,j]
        hess[i,j] <- hess[i,j]*SUM
      }
    }
    hess <- hess + t(hess) - diag(diag(hess))
    hess[] <- D[hess]
    # seems like there should have been a more consise way to do that
  }
  else
  { hess <- array(NA,c(n,n)) } # wrong for N>1

  # gradient directions: NA is symmetric
  side <- rep(NA,n)

  # do we need to calculate one-sided gradients because of a boundary?
  # parameters near zero on zero boundary
  TEST <- abs(par) <= zero.tol

  # bounded below by zero
  LO <- TEST & lower==0
  if(any(LO)) { side[LO] <- 1 }

  # bounded above by zero
  HI <- TEST & upper==0
  if(any(HI)) { side[HI] <- -1 }

  # parameters away from zero
  TEST <- !TEST

  LO <- TEST & par-d<=lower
  if(any(LO)) { side[LO] <- 1 }

  HI <- TEST & upper<=par+d
  if(any(HI)) { side[HI] <- -1 }

  # avoiding this gives us a small speed up
  if(any(!is.na(side)) || order==1)
  {
    method.args <- list(eps=1e-4,d=0.0001,zero.tol=zero.tol,r=Richardson,v=2,show.details=FALSE)
    grad <- numDeriv::jacobian(func,par,side=side,method.args=method.args,...)
  }

  # convert back from natural units
  if(is.null(dim(grad)) || drop && nrow(grad)==1)
  {
    grad <- c(grad)/parscale
    hess <- t(hess / parscale) / parscale
  }
  else
  { grad <- t( t(grad) / parscale ) }

  # Inf fix
  hess <- nant(hess,0) # seems reasonable for inverse covariance
  grad <- nant(grad,0) # not sure about this

  if(length(dim(grad)))
  {
    # TODO rownames
    colnames(grad) <- NAMES
  }
  else
  {
    names(grad) <- NAMES
    dimnames(hess) <- list(NAMES,NAMES)
  }

  return(list(gradient=grad,hessian=hess))
}












################################
# multi-core second-order derivatives
genD.mcDeriv <- function(par,fn,zero=0,lower=-Inf,upper=Inf,PERIOD=F,step=NULL,covariance=NULL,cores=detectCores(),cheap=FALSE)
{
  DIM <- length(par)

  # what points are on the boundary, so that we have to do one-sided derivatives
  LO <- (par <= lower)
  UP <- (par >= upper)
  DBOX <- UP - LO # encoded with directions to boundary
  BOX <- as.logical(DBOX)

  # delete off boxed cross correlations
  if(any(BOX) && any(!BOX))
  {
    covariance[BOX,!BOX] <- 0
    covariance[!BOX,BOX] <- 0
  }
  stddev <- sqrtm(covariance)
  stdize <- PDsolve(stddev)

  # initialize gradient and hessian
  grad <- array(0,DIM)
  hess <- array(0,c(DIM,DIM))

  # central evaluation point
  P <- cbind(par)
  # standardized evaluation points
  S <- cbind(rep(0,DIM))

  # evaluation points perpendicular to boundaries
  if(any(BOX))
  {
    STD <- sqrt(diag(covariance)[BOX])

    DIM.BOX <- sum(BOX)
    # boxed directions
    DIR.BOX <- diag(DIM)[,BOX,drop=FALSE]

    # pick points away from the boundary for one-sided derivatives
    # 2*step (reflected)
    Q <- -2*diag(DBOX)[,BOX,drop=FALSE]
    Q <- (step * stddev) %*% Q
    Q <- line.boxer(Q,par,lower=lower,upper=upper,period=PERIOD)

    # 1*step (regular)
    P <- cbind(P,Q,(par+Q)/2)

    Q <- stdize %*% (Q-par)
    S <- cbind(S,Q,Q/2)
  }

  # regular evaluation points
  if(any(!BOX))
  {
    DIM.FREE <- sum(!BOX)
    DIR.FREE <- diag(DIM)[,!BOX,drop=FALSE]

    # x-x diagonal terms
    Q <- (step * stddev) %*% DIR.FREE
    Q <- cbind(Q,-Q)
    Q <- line.boxer(Q,par,lower=lower,upper=upper,period=PERIOD)
    P <- cbind(P,Q)

    Q <- stdize %*% (Q-par)
    S <- cbind(S,Q)

    # cross terms
    if(sum(!BOX)>1)
    {
      DIM.CROSS <- (DIM.FREE^2-DIM.FREE)/2

      # x-y cross terms
      DIR.MINUS <- lapply(1:(DIM.FREE-1),function(i){ DIR.FREE[,i] - DIR.FREE[,-(1:i)] })
      DIR.MINUS <- do.call(cbind,DIR.MINUS)/sqrt(2)

      Q <- (step * stddev) %*% DIR.MINUS
      Q <- cbind(Q,-Q)
      Q <- line.boxer(Q,par,lower=lower,upper=upper,period=PERIOD)
      P <- cbind(P,Q)

      Q <- stdize %*% (Q-par)
      S <- cbind(S,Q)

      # x+y cross terms
      DIR.PLUS <- lapply(1:(DIM.FREE-1),function(i){ DIR.FREE[,i] + DIR.FREE[,-(1:i)] })
      DIR.PLUS <- do.call(cbind,DIR.PLUS)/sqrt(2)

      Q <- (step * stddev) %*% DIR.PLUS
      Q <- cbind(Q,-Q)
      Q <- line.boxer(Q,par,lower=lower,upper=upper,period=PERIOD)
      P <- cbind(P,Q)

      Q <- stdize %*% (Q-par)
      S <- cbind(S,Q)
    }
  }

  if(zero)
  { func <- function(p) { fn(p,zero=zero) } }
  else # doesn't require zero argument
  { func <- function(p) { fn(p) } }

  FN <- unlist(plapply(split(P,col(P)),func,cores=cores))

  MIN <- which.min(FN)
  par.best <- P[,MIN]
  fn.best <- FN[MIN]

  # central evaluation point
  fn.par <- FN[1]
  FN <- FN[-1]
  S <- S[,-1,drop=FALSE]

  # evaluation points perpendicular to boundaries
  if(any(BOX))
  {
    P1 <- S[,1:DIM.BOX,drop=FALSE]
    F1 <- FN[1:DIM.BOX]
    S <- S[,-(1:DIM.BOX),drop=FALSE]
    FN <- FN[-(1:DIM.BOX)]

    P2 <- S[,1:DIM.BOX,drop=FALSE]
    F2 <- FN[1:DIM.BOX]
    S <- S[,-(1:DIM.BOX),drop=FALSE]
    FN <- FN[-(1:DIM.BOX)]

    DIFF <- QuadSolve(rep(0,DIM),P1,P2,DIR.BOX,fn.par,F1,F2)
    grad[BOX] <- DIFF$GRAD
    diag(hess)[BOX] <- DIFF$hessian
  }

  # regular evaluation points
  if(any(!BOX))
  {
    P1 <- S[,1:DIM.FREE,drop=FALSE]
    F1 <- FN[1:DIM.FREE]
    S <- S[,-(1:DIM.FREE),drop=FALSE]
    FN <- FN[-(1:DIM.FREE)]

    P2 <- S[,1:DIM.FREE,drop=FALSE]
    F2 <- FN[1:DIM.FREE]
    S <- S[,-(1:DIM.FREE),drop=FALSE]
    FN <- FN[-(1:DIM.FREE)]

    DIFF <- QuadSolve(rep(0,DIM),P1,P2,DIR.FREE,fn.par,F1,F2)
    #DEBUG <<- list(grad=grad,hess=hess,BOX=BOX,DIFF=DIFF)
    grad[!BOX] <- DIFF$GRAD
    diag(hess)[!BOX] <- DIFF$hessian

    # cross terms
    if(sum(!BOX)>1)
    {
      # x-y cross term
      P1 <- S[,1:DIM.CROSS,drop=FALSE]
      F1 <- FN[1:DIM.CROSS]
      S <- S[,-(1:DIM.CROSS),drop=FALSE]
      FN <- FN[-(1:DIM.CROSS)]

      P2 <- S[,1:DIM.CROSS,drop=FALSE]
      F2 <- FN[1:DIM.CROSS]
      S <- S[,-(1:DIM.CROSS),drop=FALSE]
      FN <- FN[-(1:DIM.CROSS)]

      H.MINUS <- QuadSolve(rep(0,DIM),P1,P2,DIR.MINUS,fn.par,F1,F2)$hessian

      # x+y cross term
      P1 <- S[,1:DIM.CROSS,drop=FALSE]
      F1 <- FN[1:DIM.CROSS]
      S <- S[,-(1:DIM.CROSS),drop=FALSE]
      FN <- FN[-(1:DIM.CROSS)]

      P2 <- S[,1:DIM.CROSS,drop=FALSE]
      F2 <- FN[1:DIM.CROSS]
      S <- S[,-(1:DIM.CROSS),drop=FALSE]
      FN <- FN[-(1:DIM.CROSS)]

      H.PLUS <- QuadSolve(rep(0,DIM),P1,P2,DIR.PLUS,fn.par,F1,F2)$hessian

      H.CROSS <- (H.PLUS-H.MINUS)/2

      hess[upper.tri(hess,diag=FALSE)] <- H.CROSS
      hess <- hess + t(hess) - diag(diag(hess))
    }
  }

  condition <- abs(diag(hess))
  condition <- max(condition)/min(condition)

  hess <- stdize %*% hess %*% t(stdize)
  grad <- stdize %*% grad

  RETURN <- list()
  RETURN$gradient <- grad
  RETURN$hessian <- hess
  RETURN$condition <- condition
  RETURN$par.best <- par.best
  RETURN$fn.best <- fn.best

  return(RETURN)
}

