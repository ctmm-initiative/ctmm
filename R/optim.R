# convenience wrapper for numDeriv::genD on scalar functions
###########################################################
genD.numDeriv <- function(par,fn,zero=0,lower=-Inf,upper=Inf,step=NULL,r=2,covariance=NULL,...)
{
  # # square root of covariance
  # stddev <- sqrtm(covariance)
  # parscale <- diag(stddev)
  # # approximate standardization
  # stdize <- PDsolve(stddev)
  # 
  # lower <- lower/parscale
  # upper <- upper/parscale
  # 
  # # fn in terms of approximately standardized variable and possibly zeroed
  # if(zero)
  # {
  #   fn.std <- function(p) { fn(par + c(stddev %*% p),zero=zero) }
  #   fn.scl <- function(p) { fn(par + c(parscale * p),zero=zero) } # preserves boundaries
  # }
  # else # doesn't require zero argument
  # { 
  #   fn.std <- function(p) { fn(par + c(stddev %*% p)) }
  #   fn.scl <- function(p) { fn(par + c(parscale * p)) } # preserves boundaries
  # }
  # 
  # # numDeriv arguments to use
  # # eps <- max(step,sqrt(32*r*.Machine$double.eps))
  # # 32 checked by hand at precision=1... strangely magic number...
  # eps <- 1e-4 # default
  # d <- sqrt(eps) # relative step fraction
  # # can't make these too small for some reason
  # zero.tol <- max(sqrt(2*.Machine$double.eps),2*d) # fixes step to eps (absolute)
  # method.args <- list(r=r,eps=eps,zero.tol=zero.tol,d=d)

  # differentiate standardized function
  # D <- numDeriv::genD(fn.std,rep(0,n),method.args=method.args,...)$D

  # default numDeriv arguments except r=2... seems very sensitive otherwise
  d <- 0.01
  zero.tol <- sqrt(.Machine$double.eps/7e-7)
  D <- numDeriv::genD(fn,par,method.args=list(eps=1e-4,d=d,zero.tol=zero.tol,r=r,v=2,show.details=FALSE),...)$D
    
  n <- length(par)
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
  
  # garbage for now
  condition <- diag(hess)
  condition <- min(sign(condition))*max(abs(condition),1)/min(abs(condition),1)
  
  # transform back from pre-conditioning
  # hess <- stdize %*% hess %*% t(stdize)
  
  # avoiding this gives us a small speed up
  if(any(!is.na(side)))
  {
    # grad <- numDeriv::grad(fn.scl,rep(0,n),side=side,method.args=method.args,...)
    grad <- numDeriv::grad(fn,par,side=side,method.args=list(eps=1e-4,d=0.0001,zero.tol=zero.tol,r=r,v=2,show.details=FALSE),...)
    # grad <- grad/parscale
  }
  # else
  # { grad <- stdize %*% grad }
  
  return(list(gradient=grad,hessian=hess,condition=condition))
}


################################
# multi-core second-order derivatives
genD.mcDeriv <- function(par,fn,zero=0,lower=-Inf,upper=Inf,step=NULL,covariance=NULL,mc.cores=parallel::detectCores(),cheap=FALSE)
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
    Q <- line.boxer(Q,par,lower=lower,upper=upper)
    
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
    Q <- line.boxer(Q,par,lower=lower,upper=upper)
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
      Q <- line.boxer(Q,par,lower=lower,upper=upper)
      P <- cbind(P,Q)
      
      Q <- stdize %*% (Q-par)
      S <- cbind(S,Q)
      
      # x+y cross terms
      DIR.PLUS <- lapply(1:(DIM.FREE-1),function(i){ DIR.FREE[,i] + DIR.FREE[,-(1:i)] })
      DIR.PLUS <- do.call(cbind,DIR.PLUS)/sqrt(2)
      
      Q <- (step * stddev) %*% DIR.PLUS
      Q <- cbind(Q,-Q)
      Q <- line.boxer(Q,par,lower=lower,upper=upper)
      P <- cbind(P,Q)
      
      Q <- stdize %*% (Q-par)
      S <- cbind(S,Q)
    }
  }
  
  if(zero)
  { func <- function(p) { fn(p,zero=zero) } }
  else # doesn't require zero argument
  { func <- function(p) { fn(p) } }
  
  FN <- unlist(parallelsugar::mclapply(split(P,col(P)),func,mc.cores=mc.cores))
  
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


###############
# calculate first and second derivatives efficiently
# assumes to start with approximant of maximum (par) and inverse hessian (covariance)
####################
genD <- function(par,fn,zero=FALSE,lower=-Inf,upper=Inf,step=NULL,precision=1/2,covariance=NULL,mc.cores=parallel::detectCores(),Richardson=2)
{
  DIM <- length(par)
 
  if(is.null(covariance))
  {
    parscale <- pmin(abs(par),abs(par-lower),abs(upper-par))
    if(any(parscale==0)) { parscale[parscale==0] <- 1 }
    covariance <- diag(parscale^2)
  }
  
  if(is.null(step)) { step <- sqrt(2*.Machine$double.eps^precision) }
  
  if(Richardson==1) # parallelized, but no Richardson extrapolation (2nd order)
  { RETURN <- genD.mcDeriv(par,fn,zero=zero,lower=lower,upper=upper,step=step,covariance=covariance,mc.cores=mc.cores) }
  else # Richardson extrapolation, but no parallelization
  { RETURN <- genD.numDeriv(par,fn,zero=zero,lower=lower,upper=upper,step=step,covariance=covariance,r=Richardson) }
  
  return(RETURN)
}


###################################
# make a wrapper that applies optim then afterwards numDeriv, possibly on a boundary with one-sided derivatives if necessary
Optimizer <- function(par,fn,...,method="Nelder-Mead",lower=-Inf,upper=Inf,control=list())
{
  method <- match.arg(method,c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent","pNewton"))
  
  precision <- maxit <- NULL
  default <- list(precision=1/2,maxit=.Machine$integer.max,parscale=pmin(abs(par),abs(par-lower),abs(upper-par)))
  control <- replace(default,names(control),control)
  # check does not like attach
  NAMES <- names(control)
  for(i in 1:length(control)) { assign(NAMES[i],control[[i]]) }
  
  if(any(parscale==0)) { parscale[parscale==0] <- 1 }
  
  if(method=="pNewton") # use mc.optim
  { RESULT <- mc.optim(par=par,fn=fn,lower=lower,upper=upper,control=control) }
  else if(length(par)==1)
  {
    # try optimize (can fail with high skew & kurtosis)
    tol <- 3*sqrt(2*.Machine$double.eps^precision) # digit precision
    ATTEMPT <- stats::optimize(f=fn,...,lower=max(lower,-10*abs(par)),upper=min(upper,10*abs(par)),tol=tol)
    RESULT <- rbind(c(ATTEMPT$minimum,ATTEMPT$objective))

    # log scale backup that can't capture zero boundary
    ndigit <- -log((.Machine$double.eps)^precision,10)
    steptol <- sqrt(2*.Machine$double.eps^precision)
    ATTEMPT <- stats::nlm(function(p){f=fn(par*exp(p))},p=0,...,ndigit=ndigit,gradtol=0,steptol=steptol,stepmax=log(10),iterlim=maxit)
    RESULT <- rbind(RESULT,c(par*exp(ATTEMPT$estimate),ATTEMPT$minimum))
    
    # choose the better estimate
    MIN <- which.min(RESULT[,2])
    RESULT <- list(par=RESULT[MIN,1],value=RESULT[MIN,2])
  }
  else     # use optim
  {
    control$precision <- NULL
    if(!is.null(control$covariance)) { control$parscale <- sqrt(diag(control$covariance)) }
    control$covariance <- NULL
    if(method=="Nelder-Mead")
    {
      lower <- -Inf
      upper <- Inf
      control$reltol <- .Machine$double.eps^precision
    }
    else
    { 
      control$factr <- .Machine$double.eps^precision
    }
    RESULT <- stats::optim(par=par,fn=fn,method=method,lower=lower,upper=upper,control=control,...)
  }
  
  return(RESULT)
}


######################
# apply box constraints to travel in a straight line from p0 by dp
######################
line.boxer <- function(dp,p0=dp[,1],lower=-Inf,upper=Inf)
{
  DIM <- dim(dp)
  if(is.null(dim(dp))) { dp <- cbind(dp) }
  
  single.boxer <- function(dp)
  {
    # where would we go without a boundary
    p <- p0 + dp
    
    # did we hit a boundary?
    LO <- (p <= lower)
    UP <- (p >= upper)
    
    # are we trying to push through that boundary?
    if(any(LO)) { LO <- LO & (dp[LO]<0) }
    if(any(UP)) { UP <- UP & (dp[UP]>0) }
    
    # stop when we hit the boundary  
    # time until we hit the first new boundary
    if(any(LO)) { t.lo <- (lower-p0)[LO]/dp[LO] } else { t.lo <- 1 }
    if(any(UP)) { t.up <- (upper-p0)[UP]/dp[UP] } else { t.up <- 1 }
    t <- min(t.lo,t.up)
    
    # stop at first boundary
    p <- p0 + t*dp
    
    return(p)
  }
  
  p <- apply(dp,2,function(d){single.boxer(d)})
  dim(p) <- DIM
  
  return(p)
}


###########################################
# calculate local derivatives (along directions DIR) from 3 points, with (par,fn.par) the best (global variables)
###########################################
QuadSolve <- function(P0,P1,P2,DIR,F0,F1,F2)
{
  # convert back to displacements
  P1 <- P1 - P0
  P2 <- P2 - P0
  # signed magnitudes of displacement vectors
  # P1, P2, DIR columns will always be parallel
  M1 <- colSums(DIR * P1)
  M2 <- colSums(DIR * P2)
  
  F1 <- F1-F0
  F2 <- F2-F0
  
  G1 <- F1/M1
  G2 <- F2/M2
  # Hessian estimates, not necessarily assuming that par is between P1 & P2
  hessian <- (G2-G1)/(M2-M1)*2
  
  # gradient estimates, also not assuming that par is between P1 & P2
  GRAD <- (G2/M2-G1/M1)/(1/M2-1/M1)
  
  return(list(GRAD=GRAD,hessian=hessian)) 
}


#################################
# simplex & quasi-Newton optimizers for multiple CPU cores
# TODO: replace search lines with steepest ascent search curves?
# TODO: filling up mc queue if p>2n+1
##################################
mc.optim <- function(par,fn,...,lower=-Inf,upper=Inf,control=list())
{
  # check complains about visible bindings
  fnscale <- parscale <- maxit <- precision <- trace <- mc.cores <- hessian <- covariance <- zero <- NULL
  # fix default control arguments
  default <- list(fnscale=1,parscale=pmin(abs(par),abs(par-lower),abs(upper-par)),maxit=100,trace=FALSE,precision=1/2,mc.cores=parallel::detectCores(),hessian=NULL,covariance=NULL,zero=FALSE)
  control <- replace(default,names(control),control)
  # check does not like attach
  NAMES <- names(control)
  for(i in 1:length(control)) { assign(NAMES[i],control[[i]]) }
  
  if(any(parscale==0)) { parscale[parscale==0] <- 1 }
  
  # error tolerance for estimate, not objective function
  # half this for optimal parameters
  ERROR.TOL <- .Machine$double.eps^precision
  # the smallest relative parameter step we can take that will change the objective function by >0 numerically
  STEP.TOL <- sqrt(2*.Machine$double.eps^precision)
  STEP.MIN <- sqrt(2*.Machine$double.eps)
  
  ZERO <- zero # ZERO will be boolean, zero will be constant
  zero <- 0
  # what we will actually be evaluating
  par <- par/parscale
  if(ZERO) # this objective function has the ability to approximately zero its objective value
  { func <- function(par,...) { fn(par*parscale,zero=zero*fnscale,...)/fnscale } }
  else
  { func <- function(par,...) { fn(par*parscale,...)/fnscale } }
  
  DIM <- length(par)
  lower <- array(lower,DIM)/parscale
  upper <- array(upper,DIM)/parscale
  
  # start with the canonical directions
  DIR <- diag(1,DIM) # rows of column vectors
  
  # cost of full 2nd-order Hessian calculation relative to diagonal 2nd-order Hessian calculation
  COST.HESS <- ceiling((1+2*DIM+(DIM^2-DIM)/2*4)/mc.cores)/ceiling((1+2*DIM)/mc.cores)
  # cost of Newton-Raphson iteration relative to line search finisher
  COST.LINE <- ceiling((1+2*DIM)/mc.cores)/ceiling(2/mc.cores)
  
  # do we have better information than parscale?
  # these are all post-fnscale parameters
  if(!is.null(covariance))
  {
    covariance <- t(t(covariance/parscale)/parscale)
    if(is.null(hessian)) { hessian <- PDsolve(covariance) }
  }
  else if(!is.null(hessian))
  {
    hessian <- t(t(hessian*parscale)*parscale)
    if(is.null(covariance)) { covariance <- PDsolve(hessian) }
  }
  else # use parscale
  {
    covariance <- diag(DIM) # inverse hessian
    hessian <- diag(DIM) # hesssian
  }
  
  ######################
  # minimally rotate orthonormal basis DIR to include vec
  ######################
  Gram.Schmidt <- function(DIR,vec)
  {
    # normalize par.diff
    vec <- vec / sqrt(sum(vec^2))
    
    # calculate overlaps
    OVER <- colSums(vec*DIR)
    # replace dimension of largest overlap with par.diff
    MAX <- which.max(abs(OVER))
    DIR[,MAX] <- vec
    OVER[MAX] <- 1 # not necessary, but true
    
    # place vec is in its most canonical slot
    SWAP <- which.max(abs(vec))
    DIR[,c(MAX,SWAP)] <- DIR[,c(SWAP,MAX)]
    OVER[c(MAX,SWAP)] <- OVER[c(SWAP,MAX)]
    MAX <- SWAP # vec is now here
    
    # subtract projection from all but vec
    DIR[,-MAX] <- DIR[,-MAX] - (vec %o% OVER[-MAX])
    
    # re-normalization may be necessary
    MAG <- sqrt(colSums(DIR^2))
    DIR <- t(t(DIR)/MAG)
    
    return(DIR)
  }
  
  
  ##################
  # am I on a boundary and if so, align DIR
  ######################
  is.boxed <- function(par)
  {
    # what points are on the boundary, so that we have to do one-sided derivatives
    LO <- (par <= lower)
    UP <- (par >= upper)
    DBOX <- UP - LO # encoded with directions to boundary
    BOX <- as.logical(DBOX) 
    
    # rotate coordinates so that BOXed dimensions are fixed/isolated in DIR and DIR[,BOX] points strictly towards boundary
    for(i in which(BOX))
    {
      dir <- numeric(DIM)
      dir[i] <- DBOX[i]
      DIR <<- Gram.Schmidt(DIR,dir)
    }

    return(BOX)
  }
  
  # initializing stuff checked in loop
  counts <- 0
  ERROR <- sqrt(DIM) # error (relative to std dev) of all dimensions summed in quadrature
  ERROR.RATE <- 1 # estimated rate of convergence
  LINE.DO <- TRUE # else do line searches in between Newton-Raphson iterations (always initially)
  LINE.DID <- TRUE # did we do a line search previously
  par.target <- par # where to evaluate around next
  fn.par <- Inf # current best objective value fn(par)
  condition <- Inf
  par.diff.old <- rep(Inf,DIM)
  hessian.L <- NULL # line-search hessian
  ######################
  # MAIN LOOP
  ######################
  while(counts < maxit)
  {
    # adjust zero shift
    if(ZERO && fn.par<Inf)
    {
      zero <- zero + fn.par
      # fn.par <- 0 # subject to roundoff error
    }
    
    ################################
    # STEP 1: O(DIM) differentiation
    ################################
    # we differentiate in DIM current directions DIR from point par
    counts <- counts + 1
    
    BOX <- is.boxed(par.target)
    
    # transform Hessian to current coordinates
    hessian <- t(DIR) %*% hessian %*% DIR
    covariance <- t(DIR) %*% covariance %*% DIR
    # standard deviations along current axes 
    STD <- sqrt(abs(diag(covariance)))
    # can update this element from line search... will only help for derivative spacing
    if(!is.null(hessian.L))
    { 
      STD[BEST] <- 1/sqrt(abs(hessian.L))
      hessian.L <- NULL
    }
    
    # sample initial points surrounding the initial point for numerical differentiation
    ERROR.DIM <- sqrt(2/DIM)*ERROR.RATE*ERROR # estimated current RMS error (linear & relative to standard deviation) per dimension
    STEP <- max(sqrt(ERROR.DIM*STEP.TOL),STEP.MIN) # put the differentiation step between the error estimate and the numerical tolerance, with a lower limit
    par.step <- STEP * t(STD*t(DIR))
    P1 <- -par.step # away from boundaries
    P2 <- +par.step # towards boundaries
    # don't sample points across boundary, fold back instead - correct one-sided derivatives implemented
    if(any(BOX)) { P2[,BOX] <- -par.step[,BOX]*2 }
    # we could potentially get cornered doing this... but seems unlikely !!!
    
    # apply these displacements within box constraints
    P1 <- line.boxer(P1,par.target,lower=lower,upper=upper)
    P2 <- line.boxer(P2,par.target,lower=lower,upper=upper)
    # columns are now boxed coordinates
    
    # mc evaluate all points
    P <- cbind(par,par.target,P1,P2)
    # par must be re-evaluated againt to prevent zero shift roundoff error
    fn.queue <- unlist(parallelsugar::mclapply(split(P,col(P)),func,mc.cores=mc.cores))
    # separate back into parts
    fn.par <- fn.queue[1]
    fn.target <- fn.queue[2]
    F1 <- fn.queue[2+1:DIM]
    F2 <- fn.queue[2+DIM+1:DIM]
    
    # calculate axial derivatives to second order
    DIFF <- QuadSolve(par.target,P1,P2,DIR,fn.target,F1,F2)
    GRAD <- DIFF$GRAD
    H <- DIFF$hessian
    
    # will need line search if not normal looking
    if(any(H<=.Machine$double.eps)) { LINE.DO <- TRUE }
    
    # if the gradient is null in one direction, then we've likely profiled that axis to machine precision
    TEST <- abs(F1-fn.target)<=.Machine$double.eps | abs(F2-fn.target)<=.Machine$double.eps
    if(any(TEST))
    {
      GRAD[TEST] <- 0
      H[TEST] <- diag(hessian)[TEST]
    }
    # fix profiled dimensions to change nothing, else will get 0/0
    
    # Make sure the Hessian esitmate didn't blow up
    TEST <- abs(H)<=.Machine$double.eps
    if(any(TEST)) { H[TEST] <- diag(hessian)[TEST] }
    # stick with old estimate
    
    condition <- H/diag(hessian)
    condition <- min(sign(condition))*max(abs(condition),1)/min(abs(condition),1)
    
    # curvature correction factor: new curvature / old curvature (current coordinates)
    FACT <- sqrt(abs(H/diag(hessian))) # correction factor diagonal
    # update curvatures while preserving correlations
    hessian <- t(t(FACT*hessian)*FACT)
    # update covariances the same way as Hessian (prevents requirement of matrix inversion)
    covariance <- t(t(covariance/FACT)/FACT)
    
    # transform gradient to canonical coordinates
    GRAD <- c(DIR %*% GRAD)
    
    # transform Hessian back to canonical coordinates
    hessian <- DIR %*% hessian %*% t(DIR)
    covariance <- DIR %*% covariance %*% t(DIR)
    
    # Newton-Raphson search step from par.target
    par.diff <- -c(covariance %*% GRAD)
    
    #  this only happens when line search was skipped prior... so start doing line searches again
    if(fn.par <= fn.target && !LINE.DID) { LINE.DO <- TRUE }
    
    # pick all-time best evaluated location as our new starting point in line search
    MIN <- which.min(fn.queue)
    # new best par
    par <- P[,MIN]
    fn.par <- fn.queue[MIN]
    # choose best starting point for search step
    if(MIN!=2) # 2 is the middle point par.target, from which search step was calculated
    {
      par.diff2 <- P[,MIN]-par.target # step from old start to new (evaluated) start
      par.diff <- par.diff - par.diff2 # step from new (evaluated) par to same Newton-Raphson estimate
      GRAD <- GRAD + c(hessian %*% par.diff2) # gradient at the new par
    }
    
    # don't search past boundary, but along the boundary
    BOX <- is.boxed(par) # par could have been updated
    if(any(BOX & (par.diff%*%DIR)>0))
    {
      par.diff[BOX] <- 0
      # boundary restricted search
      if(any(!BOX)) { par.diff[!BOX] <- -PDsolve(hessian[!BOX,!BOX]) %*% GRAD[!BOX] }
    }
    
    # test for stopping condition here
    ERROR.OLD <- sqrt(c(par.diff.old %*% hessian %*% par.diff.old)) # old error, but with new hessian estimate
    par.diff -> par.diff.old # store for next time
    ERROR <- sqrt(c(par.diff %*% hessian %*% par.diff)) # relative error
    # linear convergence rate estimate (underestimates super-linear error)
    if(!is.nan(ERROR.OLD) && ERROR.OLD<Inf) { ERROR.RATE <- min(1,ERROR/ERROR.OLD,na.rm=TRUE) }
    else { ERROR.RATE <- 1 }
    # ERROR <- min(1,ERROR) # should I do this
    if(ERROR < STEP.TOL) { break } # smallest numerical change to objective function near minimum
    if(ZERO && max(abs(fn.queue[1:2]))<ERROR.TOL) { break } # absolute tolerance for zeroed objective functions
    if(!ZERO && abs((fn.par-fn.queue[1])/fn.par)<ERROR.TOL) { break } # relative tolerance for non-zeroed objective functions
    
    # calculate best DIR and orthonormalize the remaining DIRs
    BEST <- which.max(abs(par.diff)) # future index of par.diff direction
    DIR <- Gram.Schmidt(DIR,par.diff)
    
    # where we aim to evaluate next
    par.target <- line.boxer(par.diff,p0=par,lower=lower,upper=upper)
    
    # DEBUG <- list()
    # DEBUG$xy <- cbind(par.target,P) * parscale / c(60^2*24,60)
    # DEBUG$z <- c(NA,fn.queue)
    # DEBUG$cex <- DEBUG$z - min(DEBUG$z,na.rm=TRUE)
    # DEBUG$cex <- 1+2*DEBUG$cex/max(DEBUG$cex[DEBUG$cex<Inf],na.rm=TRUE)
    # DEBUG$cex[1] <- 1
    # DEBUG$col <- DEBUG$z - fn.target
    # DEBUG$col <- sapply(1:length(DEBUG$col),function(i){ col <- DEBUG$col[i] ; if(is.na(col)) { "black" } else if(col>0 && col<Inf) { "red" } else if(col<0) { "blue" } else { "black" } })
    # DEBUG$pch <- c(4,2,4,rep(1,length(DEBUG$col)-3))
    # plot(DEBUG$xy[1,],DEBUG$xy[2,],cex=DEBUG$cex,col=DEBUG$col,xlab="Days",ylab="Minutes",pch=DEBUG$pch)
    # title("Newton-Raphson step")
    # lines(DEBUG$xy[1,2:3],DEBUG$xy[2,2:3],col="grey")
    # lines(DEBUG$xy[1,c(1,MIN+1)],DEBUG$xy[2,c(1,MIN+1)],col="grey")
    
    if(trace) { message(sprintf("%f Newton-Raphson step (condition=%f)",zero+fn.par,condition)) }
    LINE.DID <- FALSE
    ##################
    # LINE SEARCH LOOP
    ##################
    if(LINE.DO)
    {
      trace.string <- "Forward"
      
      LINE.DID <- TRUE
      LINE.DO <- FALSE
      LINE.SORTED <- TRUE
      
      # generate a linear sequence of points from old par to the other side of new par
      P <- line.boxer(2*par.diff,p0=par,lower=lower,upper=upper)
      par.diff <- P - par # twice the old par.diff with no boundary (reflection=1)
      M <- sqrt(sum(par.diff^2)) # total search magnitude
      SEQ <- seq(0,M,length.out=max(2,mc.cores)+1)[-1]
      P <- (DIR[,BEST] %o% SEQ) + par
      
      # store for later check
      par.start <- par
      fn.start <- fn.par
      
      # initialize (all) storage results
      par.all <- cbind(par)
      fn.all <- fn.par
      
      # start iteration loop
      while(counts < maxit)
      {
        counts <- counts + 1
        
        # most expensive part
        # evaluate objective function at new P and store to fn.queue
        fn.queue <- unlist(parallelsugar::mclapply(split(P,col(P)),func,mc.cores=mc.cores))
        
        # combine with older results
        par.all <- cbind(par.all,P)
        fn.all <- c(fn.all,fn.queue)
        
        # sort along DIR[,BEST]
        if(!LINE.SORTED)
        {
          SORT <- colSums((par.all-par)*DIR[,BEST])
          SORT <- sort(SORT,method="quick",index.return=TRUE)$ix
          par.all <- par.all[,SORT,drop=FALSE] # ARGH!
          fn.all <- fn.all[SORT]
          
          LINE.SORTED <- TRUE
        }
        
        # new best estimate
        MIN <- which.min(fn.all)
        par <- par.all[,MIN]
        fn.par <- fn.all[MIN]
        
        # DEBUG2 <- list()
        # DEBUG2$xy <- cbind(DEBUG$xy,par.all * parscale / c(60^2*24,60))
        # DEBUG2$z <- c(DEBUG$z,fn.all)
        # DEBUG2$cex <- DEBUG2$z - min(DEBUG2$z,na.rm=TRUE)
        # DEBUG2$cex <- 1+2*DEBUG2$cex/max(DEBUG2$cex[DEBUG2$cex<Inf],na.rm=TRUE)
        # DEBUG2$cex[1] <- 1
        # DEBUG2$col <- DEBUG2$z - fn.target
        # DEBUG2$col <- sapply(1:length(DEBUG2$col),function(i){ col <- DEBUG2$col[i] ; if(is.na(col)) { "black" } else if(col>0 && col<Inf) { "red" } else if(col<0) { "blue" } else { "black" } })
        # DEBUG2$pch <- c(4,2,4,rep(1,length(DEBUG2$col)-3))
        # plot(DEBUG2$xy[1,],DEBUG2$xy[2,],cex=DEBUG2$cex,col=DEBUG2$col,xlab="Days",ylab="Minutes",pch=DEBUG2$pch)
        # title(sprintf("%s line search",trace.string))
        # points(par[1]*parscale[1]/(60^2*24) , par[2]*parscale[2]/60, col="orange")

        if(trace) { message(sprintf("%f %s search",zero+fn.par,trace.string)) }
        
        # numerical degeneracy check
        if((MIN<length(fn.all) && fn.all[MIN]==fn.all[MIN+1]) || (MIN>1 && fn.all[MIN]==fn.all[MIN-1]) || min(fn.all[MIN+c(1,-1)],na.rm=TRUE)-fn.all[MIN]<=ERROR.TOL)
        { 
          par.target <- par
          break 
        }

        ################################################
        # do we need to keep going to capture the minimum?
        ###################################################
        if(all(par==par.start) || MIN==1) # we went too far, ... steep cliff after minimum, or need to back up a touch
        {
          LINE.DO <- TRUE
          
          # left boundary relfect
          if(MIN>1) { i <- -1 } else { i <- 2 }
          P1 <- cbind(par.all[,MIN+i])
          F1 <- fn.all[MIN+i]
          
          # right boundary relfect
          if(MIN<length(fn.all)) { i <- 1 } else { i <- -2 }
          P2 <- cbind(par.all[,MIN+i])
          F2 <- fn.all[MIN+i]
          
          # calculate gradient and curvature
          DIFF <- QuadSolve(par,P1,P2,DIR[,BEST,drop=FALSE],fn.par,F1,F2)
          GRAD <- DIFF$GRAD
          hessian.L <- DIFF$hessian
          # will update search direction hessian with this during next iteration
          par.diff <- -(GRAD/abs(hessian.L))
          
          # boundary limits
          M1 <- min((-(lower-par)/DIR[,BEST])[-DIR[,BEST]<0],(-(upper-par)/DIR[,BEST])[-DIR[,BEST]>0],Inf)
          M2 <- min(((lower-par)/DIR[,BEST])[DIR[,BEST]<0],((upper-par)/DIR[,BEST])[DIR[,BEST]>0],Inf)
          
          if(MIN==1 && par.diff<0)
          {
            trace.string <- "Reverse"
            
            # geometrically expanding reverse expansion 
            if(M1>abs(par.diff)) # far away from boundary
            {
              b <- min(2,(M1/abs(par.diff))^(1/max(1,mc.cores-1)))
              SEQ <- par.diff * b^(0:max(1,mc.cores-1))
            }
            else # even grid to boundary, because boundary is close
            { SEQ <- -seq(0,M1,length.out=mc.cores+1)[-1] }
          }
          else
          {
            trace.string <- "Shrink"
            
            # adjacent limits (half)
            if(MIN>1) { M1 <- min(M1,abs(DIR[,BEST]%*%(par-par.all[,MIN-1]))/2) }
            if(MIN<length(fn.all)) { M2 <- min(M2,abs(DIR[,BEST]%*%(par.all[,MIN+1]-par))/2) }
            
            # Newton-Raphson limit
            M1 <- min(M1[M1<Inf],abs(par.diff))
            M2 <- min(M2[M2<Inf],abs(par.diff))
            
            # probably won't ever happen...
            if(M1==Inf) { M1 <- M2 }
            
            # generate geometrically tightening sequence around par
            SEQ <- (1/2)^(1:ceiling(mc.cores/2)-1)
            SEQ <- c(-SEQ*M1,SEQ*M2)
          }
          
          # combine for evaluation
          P <- (DIR[,BEST] %o% SEQ) + par
          
          LINE.SORTED <- FALSE
          # goto evaluate iteration step
          next
        }
        else if(MIN==length(fn.all)) # we didn't go far enough or we hit a boundary
        {
          trace.string <- "Expansion"

          LINE.DO <- TRUE
          
          # if we hit a boundary, then we can stop at the boundary
          BOX.NEW <- (par <= lower) | (par >= upper)
          if(sum(BOX.NEW)>sum(BOX))
          { 
            par.target <- par
            break 
          }
          
          # Distance to first boundary that we will hit going in direction par.diff||DIR[,BEST]>0
          M.BOX <- min(((upper-par)/DIR[,BEST])[DIR[,BEST]>0],((lower-par)/DIR[,BEST])[DIR[,BEST]<0])
          
          # we didn't hit a boundary yet, but we will eventually hit a boundary, so let's just do that
          if(M.BOX<Inf)
          {
            # generate a linear sequence of points that terminate at the eventual boundary
            SEQ <- seq(0,M.BOX,length.out=max(2,mc.cores)+1)[-1]
            P <- (DIR[,BEST] %o% SEQ) + par
            
            # goto evaluate iteration step
            next
          }
          else # we didn't hit a boundary and we never will, because there is no boundary
          {
            # first let's update the Newton search step
            P1 <- cbind(par.all[,MIN-2])
            F1 <- fn.all[MIN-2]
            P2 <- cbind(par.all[,MIN-1])
            F2 <- fn.all[MIN-1]
            
            DIFF <- QuadSolve(par,P1,P2,DIR[,BEST,drop=FALSE],fn.par,F1,F2)
            GRAD <- DIFF$GRAD
            H <- DIFF$hessian
            
            if(H>0) { M <- max(M,-GRAD/H) }
            # make sure we are accelerating
            M <- max(M,last(diff(SEQ)))
            
            # generate an exponential sequence of points that terminate at the infinite boundary but drop the last point
            # start with a uniform sequence
            SEQ <- seq(0,1,length.out=max(2,mc.cores)+2)
            SEQ <- SEQ[-c(1,length(SEQ))]
            # log transform
            SEQ <- -log(1-SEQ)/diff(SEQ[1:2]) * M
            P <- (DIR[,BEST] %o% SEQ) + par
            
            # goto evaluate iteration step
            next
          }
        }
        else # update is safely encapsulated
        { break }
      }
      # end iteration loop
      ######################################
      # we improved our estimate and can stop the line search
      
      #####################
      # LINE SEARCH FINISHER
      ######################
      
      # interpolate to an even better point
      #####################################
      # take best triplet
      
      # left boundary relfect
      if(MIN>1) { i <- -1 } else { i <- 2 }
      P1 <- cbind(par.all[,MIN+i])
      F1 <- fn.all[MIN+i]

      # right boundary relfect
      if(MIN<length(fn.all)) { i <- 1 } else { i <- -2 }
      P2 <- cbind(par.all[,MIN+i])
      F2 <- fn.all[MIN+i]
      
      # estimate better location if we can (if we are in the middle we can)
      if(1<MIN && MIN<length(fn.all))
      { 
        # calculate gradient and curvature
        DIFF <- QuadSolve(par,P1,P2,DIR[,BEST,drop=FALSE],fn.par,F1,F2)
        GRAD <- DIFF$GRAD
        hessian.L <- DIFF$hessian
        # will update search direction hessian with this during next iteration
        
        par.diff <- -(GRAD/hessian.L)
        # need to differentiate new par.target from evaluated par to save
        par.target.LINE <- par + par.diff*DIR[,BEST]
        # this will always be boxed
        
        # how far off was our aim? fraction of line search
        # ERROR.LINE <- sqrt(((par.target-par.target.LINE)%*%DIR[,BEST])^2/((par.all[,1]-par.all[,length(fn.all)])%*%DIR[,BEST])^2)
        # what was our target (closest)
        MIN.start <- which.min( colSums((par.target-par.all)^2) )
        # should we do a line search next time?
        if(fn.start<fn.all[MIN.start]) { LINE.DO <- TRUE } #  || ERROR.LINE >= 1/4
        
        # what is the next target of Newton-Raphson
        # if(MIN==MIN.start) { par.target <- par } # trust Newton-Raphson more
        # else 
        { par.target <- par.target.LINE } # trust line search more
      }
      else
      {
        par.target <- par
        LINE.DO <- TRUE
      }
      # points(par.target[1]*parscale[1]/(60^2*24) , par.target[2]*parscale[2]/60, col="orange",pch=4)
    } # go back to differentiation
  } # end main loop
  

  if(counts<maxit) { convergence <- 0} else { convergence <- 1 }
  
  # return stuff in similar format to optim
  RETURN <- list()
  RETURN$par <- par*parscale
  RETURN$value <- (fn.par+zero)*fnscale
  RETURN$counts <- counts
  RETURN$convergence <- convergence
  RETURN$hessian <- t(t(hessian/parscale)/parscale)
  RETURN$covariance <- t(t(covariance*parscale)*parscale)
  RETURN$lower <- lower*parscale
  RETURN$upper <- upper*parscale
      
  return(RETURN)
}
