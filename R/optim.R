# convenience wrapper for numDeriv::genD on scalar functions
###########################################################
genD <- function(par,fn,lower=-Inf,upper=Inf,precision=1/2,parscale=pmin(abs(par),abs(par-lower),abs(upper-par)),...)
{
  if(any(parscale==0)) { parscale[parscale==0] <- 1 }
  
  # arguments to use
  r <- 2 # minimum that works... want to use method="simple", but its very limited
  eps <- max(sqrt(2*r*.Machine$double.eps^precision),sqrt(32*r*.Machine$double.eps))
  # 32 checked by hand at precision=1... strangely magic number...
  d <- eps # relative step fraction
  zero.tol <- max(sqrt(.Machine$double.eps),2*d) # fixes step to eps (absolute)
  method.args <- list(r=r,eps=eps,zero.tol=zero.tol,d=d)

  # scale with parscale
  par <- par/parscale
  lower <- lower/parscale
  upper <- upper/parscale
  func <- function(p) { fn(parscale * p) }
  
  n <- length(par)
  D <- numDeriv::genD(func,par,method.args=list(),...)$D
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
  
  # avoiding this gives us a small speed up
  if(any(!is.na(side))) { grad <- numDeriv::grad(func,par,side=side,method.args=method.args,...) }
  
  # unscale with parscale
  grad <- grad/parscale
  hess <- t(t(hess/parscale)/parscale)
  
  return(list(gradient=grad,hessian=hess))
}


###################################
# make a wrapper that applies optim then afterwards numDeriv, possibly on a boundary with one-sided derivatives if necessary
Optimizer <- function(par,fn,...,method="Nelder-Mead",lower=-Inf,upper=Inf,control=list())
{
  method <- match.arg(method,c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent","pNewton"))
  
  precision <- maxit <- NULL
  default <- list(precision=1/2,maxit=.Machine$integer.max)
  control <- replace(default,names(control),control)
  # check does not like attach
  NAMES <- names(control)
  for(i in 1:length(control)) { assign(NAMES[i],control[[i]]) }
  
  
  if(method=="pNewton") # use poptim
  { RESULT <- poptim(par=par,fn=fn,lower=lower,upper=upper,control=control) }
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
    if(method=="Nelder-Mead")
    {
      lower <- -Inf
      upper <- Inf
      control$reltol <- .Machine$double.eps^precision
    }
    else if(method=="L-BFGS-B")
    { control$factr <- .Machine$double.eps^precision }
    RESULT <- stats::optim(par=par,fn=fn,method=method,lower=lower,upper=upper,control=control,...)
  }
  
  return(RESULT)
}


#################################
# simplex & quasi-Newton optimizers for multiple CPU cores
# TODO: replace search lines with steepest ascent search curves?
# TODO: filling up mc queue if p>2n+1
##################################
poptim <- function(par,fn,...,lower=-Inf,upper=Inf,control=list())
{
  # check complains about visible bindings
  fnscale <- parscale <- maxit <- precision <- trace <- mc.cores <- hessian <- covariance <- NULL
  # fix default control arguments
  default <- list(fnscale=1,parscale=pmin(abs(par),abs(par-lower),abs(upper-par)),maxit=100,trace=FALSE,precision=1/2,mc.cores=parallel::detectCores(),hessian=NULL,covariance=NULL)
  control <- replace(default,names(control),control)
  # check does not like attach
  NAMES <- names(control)
  for(i in 1:length(control)) { assign(NAMES[i],control[[i]]) }
  
  if(any(parscale==0)) { parscale[parscale==0] <- 1 }
  
  # error tolerance for estimate, not objective function
  # half this for optimal parameters
  ERROR.TOL <- .Machine$double.eps^(precision/2)
  # the smallest relative parameter step we can take that will change the objective function by >0 numerically
  STEP.MIN <- sqrt(2*.Machine$double.eps)
  
  # what we will actually be evaluating
  par <- par/parscale
  func <- function(par) { fn(par*parscale,...)/fnscale }
  
  DIM <- length(par)
  lower <- array(lower,DIM)/parscale
  upper <- array(upper,DIM)/parscale
  
  # start with the canonical directions
  DIR <- diag(1,DIM) # rows of column vectors
  
  # cost of Newton-Raphson iteration relative to line search finisher
  COST <- ceiling((1+2*DIM)/mc.cores)/max(1,2/mc.cores)
  
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
  hessian.L <- NULL # line-search hessian
  
  ######################
  # apply box constraints to travel in a straight line from (global) par by dp
  ######################
  line.boxer <- function(dp,p0=par.target,slide=FALSE)
  {
    # if we don't hit a boundary we go here
    p <- p0 + dp
    
    # did we hit a boundary?
    LO <- (p <= lower)
    UP <- (p >= upper)

    # are we trying to push through that boundary?
    if(any(LO)) { LO <- LO & (dp[LO]<0) }
    if(any(UP)) { UP <- UP & (dp[UP]>0) }
  
    # stop when we hit the boundary  
    if(!slide)
    {
      # time until we hit the first new boundary
      if(any(LO)) { t.lo <- (lower-p0)[LO]/dp[LO] } else { t.lo <- 1 }
      if(any(UP)) { t.up <- (upper-p0)[UP]/dp[UP] } else { t.up <- 1 }
      t <- min(t.lo,t.up)
      
      # stop at first boundary
      p <- p0 + t*dp
    }
    else # slide along the boundary (project back) --- bad for differentiation, ?good? for iteration
    {
      if(any(LO)) { p[LO] <- lower[LO] }
      if(any(UP)) { p[UP] <- upper[UP] }
    }
      
    return(p)
  }
  
  
  ######################
  # minimally rotate orthonormal basis DIR to include P
  ######################
  Gram.Schmidt <- function(DIR,par.diff)
  {
    # normalize par.diff
    par.diff <- par.diff / sqrt(sum(par.diff^2))
    
    # calculate overlaps
    OVER <- colSums(par.diff*DIR)
    # replace dimension of largest overlap with par.diff
    MAX <- which.max(abs(OVER))
    DIR[,MAX] <- par.diff
    OVER[MAX] <- 1 # not necessary, but clearer
    # OVER[,MAX] <- par.diff
    # make sure par.diff is in its most canonical slot
    SWAP <- which.max(abs(par.diff))
    DIR[,c(MAX,SWAP)] <- DIR[,c(SWAP,MAX)]
    OVER[c(MAX,SWAP)] <- OVER[c(SWAP,MAX)]
    
    # subtract projection from all but par.diff direction
    MAX <- SWAP # par.diff (MAX) was moved to SWAP
    DIR[,-MAX] <- DIR[,-MAX] - (OVER[-MAX] %o% par.diff)
    
    # re-normalization is necessary
    MAG <- sqrt(colSums(DIR^2))
    DIR <- t(t(DIR)/MAG)
    
    return(DIR)
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
  
  counts <- 0
  ERROR <- sqrt(DIM) # error (relative to std dev) of all dimensions summed in quadrature
  ERROR.RATE <- 1 # estimated rate of convergence
  LINE.DO <- TRUE # else do line searches in between Newton-Raphson iterations (always initially)
  LINE.DID <- TRUE # did we do a line search previously
  par.target <- par # where to evaluate around next
  fn.par <- Inf # current best objective value fn(par)
  ######################
  # MAIN LOOP
  ######################
  while(counts < maxit)
  {
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
    ERROR.DIM <- ERROR.RATE*ERROR/sqrt(DIM) # estimated current RMS error (linear & relative to standard deviation) per dimension
    STEP <- max(sqrt(ERROR.DIM*STEP.MIN),STEP.MIN) # put the differentiation step between the error estimate and the numerical minimum that can change the objective function
    par.diff <- STEP * t(STD*t(DIR))
    P1 <- -par.diff # away from boundaries
    P2 <- +par.diff # towards boundaries
    # don't sample points across boundary, fold back instead - correct one-sided derivatives implemented
    if(any(BOX)) { P2[,BOX] <- -par.diff[,BOX]*2 }
    
    # apply these displacements within box constraints
    P1 <- apply(P1,2,line.boxer)
    P2 <- apply(P2,2,line.boxer)
    # columns are now boxed coordinates
    
    # mc evaluate all points
    P <- cbind(par.target,P1,P2)
    fn.queue <- unlist(parallelsugar::mclapply(split(P,col(P)),func,mc.cores=mc.cores))
    # separate back into parts
    F1 <- fn.queue[2:(DIM+1)]
    F2 <- fn.queue[(DIM+2):(1+2*DIM)]
    
    # calculate axial derivatives to second order
    DIFF <- QuadSolve(P[,1],P1,P2,DIR,fn.queue[1],F1,F2)
    GRAD <- DIFF$GRAD
    H <- DIFF$hessian
    
    # if the gradient is null in one direction, then we've fully profiled that axis to machine precision
    TEST <- abs(F1-fn.queue[1])<=.Machine$double.eps | abs(F2-fn.queue[1])<=.Machine$double.eps
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
    
    # will need line search if not normal looking
    if(any(H<=.Machine$double.eps)) { LINE.DO <- TRUE }
    
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
    
    # Newton-Raphson search step from old par
    par.diff <- -c(covariance %*% GRAD)
    
    #  this only happens when line search was skipped prior... so start doing line searches again
    if(fn.queue[1] >= fn.par && !LINE.DID) { LINE.DO <- TRUE }
    # include values in selection
    P <- cbind(par,P)
    fn.queue <- c(fn.par,fn.queue)
    
    # pick all-time best evaluated location as our new starting point in line search
    MIN <- which.min(fn.queue)
    if(MIN>1)
    {
      par.diff2 <- P[,MIN]-par.target # step from old start to new (evaluated) start
      par <- P[,MIN] # new par
      fn.par <- fn.queue[MIN]
      par.diff <- par.diff - par.diff2 # step from new (evaluated) par to same Newton-Raphson estimate
      GRAD <- GRAD + (hessian %*% par.diff2) # gradient at the new par
    }
    
    # don't search past boundary, but along the boundary
    BOX <- is.boxed(par) # par could have been updated
    if(any(BOX & (par.diff%*%DIR)>0))
    {
      par.diff[BOX] <- 0
      # boundary restricted search
      par.diff[!BOX] <- -PDsolve(hessian[!BOX,!BOX]) %*% GRAD[!BOX]
    }
    
    # test for stopping condition here
    ERROR.OLD <- ERROR # old error
    ERROR <- sqrt(c(par.diff %*% hessian %*% par.diff)) # relative error 
    ERROR.RATE <- min(1,ERROR/ERROR.OLD) # lienar convergence rate estimate (underestimates super-linear error)
    # ERROR <- min(1,ERROR) # should I do this
    if(ERROR/sqrt(2) < ERROR.TOL) { break } # smallest numerical change to objective function near minimum
    
    # calculate best DIR and orthonormalize the remaining DIRs
    BEST <- which.max(abs(par.diff)) # future index of par.diff direction
    DIR <- Gram.Schmidt(DIR,par.diff)
    
    # where we aim to evaluate next
    par.target <- line.boxer(par.diff,p0=par,slide=FALSE)
    # slide=TRUE is aggressive
    
    LINE.DID <- FALSE
    ##################
    # LINE SEARCH LOOP
    ##################
    if(LINE.DO)
    {
      LINE.DID <- TRUE
      LINE.DO <- FALSE
      LINE.SORTED <- TRUE
      
      # generate a linear sequence of points from old par to the other side of new par
      P <- line.boxer(2*par.diff,p0=par)
      par.diff <- P - par # twice the old par.diff with no boundary (reflection=1)
      M <- sqrt(sum(par.diff^2)) # total search magnitude
      SEQ <- seq(0,M,length.out=max(2,mc.cores)+1)[-1]
      P <- (DIR[,BEST] %o% SEQ) + par
      
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
          par.all <- par.all[,SORT]
          fn.all <- fn.all[SORT]
          
          LINE.SORTED <- TRUE
        }
        
        # store for later check
        par.start <- par
        # new best estimate
        MIN <- which.min(fn.all)
        par <- par.all[,MIN]
        fn.par <- fn.all[MIN]
        
        # numerical degeneracy check
        if((MIN<length(fn.all) && fn.all[MIN]==fn.all[MIN+1]) || (MIN>1 && fn.all[MIN]==fn.all[MIN-1]))
        { 
          par.target <- par
          break 
        }
        
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

          # what was the shortest (normal) step previously, then for next subdivision
          M <- sort(diff(SEQ),method='quick')[2] * max(2,mc.cores)/(max(2,mc.cores)+1) # fixed subdivision factor
          
          # backstep==1
          par.diff2 <- line.boxer(-par.diff*DIR[,BEST],p0=par) - par
          M1 <- (par.diff2 %*% DIR[,BEST])
          
          # what is target step size, including reflection==1
          par.diff2 <- line.boxer(2*par.diff*DIR[,BEST],p0=par) - par
          M2 <- (par.diff2 %*% DIR[,BEST])
          
          # smallest steps
          M1 <- min(abs(M),abs(M2))*sign(M1)
          M2 <- min(abs(M),abs(M1))*sign(M2)
          
          # generate sequence around par
          SEQ <- seq(M1,M2,length.out=max(2,mc.cores))

          # combine for evaluation
          P <- (DIR[,BEST] %o% SEQ) + par
          
          LINE.SORTED <- FALSE
          # goto evaluate iteration step
          next
        }
        else if(MIN==length(fn.all)) # we didn't go far enough or we hit a boundary
        {
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
      
      # numerical degeneracy check 2/2
      if((MIN<length(fn.all) && fn.all[MIN]==fn.all[MIN+1]) || (MIN>1 && fn.all[MIN]==fn.all[MIN-1]))
      { 
        par.target <- par
        break
      }
      
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
        par.target.L <- line.boxer(par.diff*DIR[,BEST],p0=par) # need to differentiate new par.target from evaluated par to save !!!
        
        # how far off was our aim?
        ERROR.L <- sqrt(((par.target-par.target.L)%*%DIR[,BEST])^2/((par.all[,1]-par.all[,length(fn.all)])%*%DIR[,BEST])^2)
        # what was our target (closest)
        MIN <- which.min( (DIR[,BEST] %*% (par.target-par.all))^2 )
        # should we do a line search next time?
        if(ERROR.L >= 1/4 || fn.all[1]<fn.all[MIN]) { LINE.DO <- TRUE }

        # copy over line search iteration
        par.target <- par.target.L
      }
      else
      {
        par.target <- par
        LINE.DO <- TRUE
      }
    } # go back to differentiation
  } # end main loop
  

  if(counts<maxit) { convergence <- 0} else { convergence <- 1 }
  
  # return stuff in similar format to optim
  RETURN <- list()
  RETURN$par <- par*parscale
  RETURN$value <- fn.par*fnscale
  RETURN$counts <- counts
  RETURN$convergence <- convergence
  RETURN$hessian <- t(t(hessian/parscale)/parscale)
  RETURN$covariance <- t(t(covariance*parscale)*parscale)
  RETURN$lower <- lower*parscale
  RETURN$upper <- upper*parscale
      
  return(RETURN)
}