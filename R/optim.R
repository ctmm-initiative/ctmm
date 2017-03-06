
#################################
# simplex & quasi-Newton optimizers for multiple CPU cores
mcoptim <- function(par,fn,method="Newton",lower=-Inf,upper=Inf,control=list())
{
  method <- match.arg(method,c("pNewton","pSimplex"))
  # method simplex not yet supported
  
  # check complains about visible bindings
  fnscale <- parscale <- ndeps <- maxit <- reltol <- trace <- fast <- mc.cores <- hessian <- covariance <- NULL
  # fix default control arguments
  default <- list(fnscale=1,parscale=pmin(abs(par),abs(par-lower),abs(upper-par)),ndeps=1,maxit=100,reltol=sqrt(.Machine$double.eps),trace=FALSE,fast=TRUE,mc.cores=parallel::detectCores(),hessian=NULL,covariance=NULL)
  control <- replace(default,names(control),control)
  # check does not like attach
  NAMES <- names(control)
  for(i in 1:length(control)) { assign(NAMES[i],control[[i]]) }
  
  if(any(parscale==0)) { parscale[parscale==0] <- 1 }
  
  # what we will actually be evaluating
  par <- par/parscale
  func <- function(par) { fn(par*parscale)/fnscale }
  
  DIM <- length(par)
  lower <- array(lower,DIM)/parscale
  upper <- array(upper,DIM)/parscale
  
  # prototype Hessian and inverse Hessian matrices relative to parscale
  COV <- diag(DIM) # inverse hessian
  HESS <- diag(DIM) # hesssian
  GRAD <- array(0,DIM) # gradient

  # do we have better information than parscale?
  # these are all post-fnscale parameters
  if(!is.null(covariance))
  {
    COV <- t(t(covariance/parscale)/parscale)
    if(is.null(hessian)) { HESS <- PDsolve(COV) }
  }
  if(!is.null(hessian))
  {
    HESS <- t(t(hessian*parscale)*parscale)
    if(is.null(covariance))
    {
      COV <- PDsolve(HESS)
    }
  }

  # standard deviations
  STD <- sqrt(diag(COV))
  H <- NULL # future hessian updates
  
  # start with the canonical directions
  DIR <- diag(1,DIM) # rows of column vectors

  # cost of Newton-Raphson iteration relative to line search finisher
  COST <- ceiling((1+2*DIM)/mc.cores)/max(1,2/mc.cores)
  
  ######################
  # apply box constraints to travel in a straight line from (global) par by dp
  ######################
  line.boxer <- function(dp)
  {
    # if we don't hit a boundary we go here
    p <- par + dp
    
    # did we hit a boundary?
    LO <- (p <= lower)
    UP <- (p >= upper)

    # are we trying to push through that boundary?
    LO <- LO && (dp[LO]<0)
    UP <- UP && (dp[UP]>0)
    
    # time until we hit the first new boundary
    if(any(LO)) { t.lo <- (lower-par)[LO]/dp[LO] } else { t.lo <- 1 }
    if(any(UP)) { t.up <- (upper-par)[UP]/dp[UP] } else { t.up <- 1 }
    t <- min(t.lo,t.up)
    
    # stop at first boundary
    p <- par + t*dp
      
    return(p)
  }
  
  
  ######################
  # minimally rotate orthonormal basis DIR to include P
  ######################
  Gram.Schmidt <- function(DIR,dP)
  {
    # normalize dP
    dP <- dP / sqrt(sum(dP^2))
    
    # calculate overlaps
    OVER <- colSums(dP*DIR)
    # replace dimension of largest overlap with dP
    MAX <- which.max(abs(OVER))
    DIR[,MAX] <- dP
    OVER[MAX] <- 1 # not necessary, but clearer
    # OVER[,MAX] <- dP
    # make sure dP is in its most canonical slot
    SWAP <- which.max(abs(dP))
    DIR[,c(MAX,SWAP)] <- DIR[,c(SWAP,MAX)]
    OVER[c(MAX,SWAP)] <- OVER[c(SWAP,MAX)]
    
    # subtract projection from all but dP direction
    MAX <- SWAP # dP (MAX) was moved to SWAP
    DIR[,-MAX] <- DIR[,-MAX] - (OVER[-MAX] %o% dP)
    
    # re-normalization is necessary
    MAG <- sqrt(colSums(DIR^2))
    DIR <- t(t(DIR)/MAG)
    
    return(DIR)
  }
  
  
  ###########################################
  # calculate local derivatives (along directions DIR) from 3 points, with (par,fp) the best (global variables)
  ###########################################
  QuadSolve <- function(P1,P2,DIR,F1,F2)
  {
    # convert back to displacements
    P1 <- P1 - par
    P2 <- P2 - par
    # signed magnitudes of displacement vectors
    # P1, P2, DIR columns will always be parallel
    M1 <- colSums(DIR * P1)
    M2 <- colSums(DIR * P2)
    
    F1 <- F1-fp
    F2 <- F2-fp
    
    G1 <- F1/M1
    G2 <- F2/M2
    # Hessian estimates, not necessarily assuming that par is between P1 & P2
    HESS <- (G2-G1)/(M2-M1)*2
    
    # gradient estimates, also not assuming that par is between P1 & P2
    GRAD <- (G2/M2-G1/M1)/(1/M2-1/M1)
    
    return(list(GRAD=GRAD,HESS=HESS)) 
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
      dir <- array(0,DIM)
      dir[i] <- DBOX[i]
      DIR <<- Gram.Schmidt(DIR,dir)
    }

    return(BOX)
  }
  
  
  counts <- 0
  ERROR <- DIM # error (relative to std dev) of all dimensions summed in quadrature
  RATE <- 1 # estimated rate of convergence
  SUCCESS <- FALSE # else do line searches in between Newton-Raphson iterations (always initially)
  fp <- Inf
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
    
    BOX <- is.boxed(par)
    
    # transform Hessian to current coordinates
    HESS <- t(DIR) %*% HESS %*% DIR
    COV <- t(DIR) %*% COV %*% DIR
    # standard deviations along current axes 
    STD <- sqrt(abs(diag(COV)))
    # can update this element from line search... will only help initially
    if(!is.null(HESS.L))
    { 
      STD[BEST] <- 1/sqrt(abs(HESS.L))
      HESS.L <- NULL
    }
    
    # sample initial points surrounding the initial point for numerical differentiation
    ERROR.DIM <- RATE*ERROR^2/DIM # estimated current error (quadratic & relative to standard deviation) per dimension
    dP <- (ndeps*ERROR.DIM) * t(STD*t(DIR))
    P1 <- -dP # away from boundaries
    P2 <- +dP # towards boundaries
    # don't sample points across boundary, fold back instead - correct one-sided derivatives implemented
    if(any(BOX)) { P2[,BOX] <- -dP[,BOX]/2 }
    
    # apply these displacements within box constraints
    P1 <- apply(P1,2,line.boxer)
    P2 <- apply(P2,2,line.boxer)
    # columns are now boxed coordinates
    
    # mc evaluate all points
    P <- cbind(par,P1,P2)
    FN <- unlist(parallelsugar::mclapply(split(P,col(P)),func,mc.cores=mc.cores))
    # separate back into parts
    F1 <- FN[2:(DIM+1)]
    F2 <- FN[(DIM+2):(1+2*DIM)]
    
    # need to keep old par earlier incase of failure... this only happens when line search was skipped prior
    OLD.par <- par
    OLD.fp <- fp
    fp <- FN[1] # fn(par)
    if(fp >= OLD.fp) { SUCCESS <- FALSE }
    
    # calculate axial derivatives to second order
    DIFF <- QuadSolve(P1,P2,DIR,F1,F2)
    GRAD <- DIFF$GRAD
    H <- DIFF$HESS
    
    # if the gradient is null in one direction, then we've fully profiled that axis to machine precision
    JUNK <- (F1==fp) | (F2==fp)
    if(any(JUNK))
    {
      GRAD[JUNK] <- 0
      H[JUNK] <- HESS[JUNK,JUNK]
    }
    # fix profiled dimensions to change nothing, else will get 0/0
    
    # curvature correction factor: new curvature / old curvature (current coordinates)
    FACT <- sqrt(abs(H/diag(HESS))) # correction factor diagonal
    COR <- array(FACT,c(DIM,DIM))
    COR <- t(COR) * COR # correction factor matrix (for Hadamard product)
    # update curvatures while preserving correlations
    HESS <- COR * HESS
    # update covariances the same way as Hessian (prevents requirement of matrix inversion)
    COR <- array(1/FACT,c(DIM,DIM))
    COR <- t(COR) * COR # inverse correction factor matrix
    COV <- COR * COV
    
    # transform gradient to canonical coordinates
    GRAD <- c(DIR %*% GRAD)
    
    # transform Hessian back to canonical coordinates
    HESS <- DIR %*% HESS %*% t(DIR)
    COV <- DIR %*% COV %*% t(DIR)
    
    # Newton-Raphson search step from old par
    dP <- -c(COV %*% GRAD)
    
    # include OLD.par & OLD.fp in this consideration
    P <- cbind(P,OLD.par)
    FN <- c(FN,OLD.fp)
    # pick all-time best evaluated location as our new starting point in line search
    MIN <- which.min(FN)
    dP2 <- P[,MIN]-par # step from old par to new (evaluated) par
    par <- P[,MIN] # new par
    fp <- FN[MIN]
    dP <- dP - dP2 # step from new (evaluated) par to same Newton-Raphson estimate
    GRAD <- GRAD + (HESS %*% dP2) # gradient at the new par
    
    # don't search past boundary, but along the boundary
    BOX <- is.boxed(par) # par could have been updated
    if(any(BOX & dP[BOX]>0))
    {
      dP[BOX] <- 0
      dP[-BOX] <- -PDsolve(HESS[-BOX,-BOX],GRAD[-BOX])
    }
    
    # standard deviations along canonical axes 
    STD <- sqrt(abs(diag(COV)))
    # test for stopping condition here
    OLD.ERROR <- ERROR # old error
    ERROR <- sqrt((dP %*% HESS %*% dP)) # relative error
    if(ERROR<OLD.ERROR) { RATE <- ERROR/OLD.ERROR^2 } else { RATE <- 1 } # Newton-Raphson quadratic convergence rate
    ERROR <- min(1,ERROR) # this could happen?
    if(ERROR < reltol) { break } # we are finished here
    
    # !!! need to throw the old par in here incase !SUCCESS
    # calculate best DIR and orthonormalize the remaining DIRs
    BEST <- which.max(abs(dP)) # future index of dP direction
    DIR <- Gram.Schmidt(DIR,dP)
    
    ##################
    # LINE SEARCH LOOP
    ##################
    if(!SUCCESS)
    {
      SUCCESS <- TRUE
      
      # generate a linear sequence of points from old par to the other side of new par
      START <- par # store for later check !!!
      TARGET <- par + dP # store for later comparison
      P <- line.boxer(2*dP)
      dP <- P - par # twice the old dP with no boundary (reflection=1) # !!! replace dP with dP2 here # !!! or at least dP with DIR[,BEST] below
      M <- sqrt(sum(dP^2)) # total search magnitude
      SEQ <- seq(0,M,length.out=max(2,mc.cores)+1)[-1]
      P <- (DIR[,BEST] %o% SEQ) + par
      
      # initialize (all) storage results
      P.ALL <- cbind(par)
      F.ALL <- fp
      
      # start iteration loop
      repeat
      {
        counts <- counts + 1
        
        # most expensive part
        # evaluate objective function at new P and store to FN
        FN <- unlist(parallelsugar::mclapply(split(P,col(P)),func,mc.cores=mc.cores))
        
        # combine with older results
        P.ALL <- cbind(P.ALL,P)
        F.ALL <- c(F.ALL,FN)
        
        # sort along DIR[,BEST]
        SORT <- colSums((P.ALL-par)*DIR[,BEST])
        SORT <- sort(SORT,method="quick",index.return=TRUE)$ix
        P.ALL <- P.ALL[,SORT]
        F.ALL <- F.ALL[SORT]
        
        # new best estimate
        MIN <- which.min(F.ALL)
        par <- P.ALL[,MIN]
        fp <- F.ALL[MIN]
        
        # do we need to keep going to capture the minimum?
        
        if(MIN==1) # we went too far in first step !!! cliff problem --- iterate until different from START
        {
          SUCCESS <- FALSE
          
          # do a sub line search, including backwards reflection
          # what was the shortest step previously
          M <- min(diff(SEQ))
          # tabulate forward steps with half of cores (assuming multiple of 2)
          # change this to check START somehow !!! maybe need another condition
          SEQ <- seq(0,M,length.out=max(1,mc.cores/2)+2)
          SEQ <- SEQ[-c(1,length(SEQ))]
          # tabulate backward steps other half of cores
          SEQ <- c(rev(seq(0,-M,length.out=max(1,mc.cores/2)+1)[-1]),SEQ)
          
          # combine for evaluation
          P <- (DIR[,BEST] %o% SEQ) + par
          
          # goto evaluate iteration step
          next
          
          # if !fast, this might not be an improvement...
          # need a conditional or fast by default
        }
        else if(MIN==length(F.ALL)) # we didn't go far enough or we hit a boundary
        {
          SUCCESS <- FALSE
          
          # if we hit a boundary, then we can stop at the boundary
          NEW.BOX <- (par <= lower) | (par >= upper)
          if(sum(NEW.BOX)>sum(BOX)) { break } # goto hessian stuff
          
          # Distance to first boundary that we will hit going in direction dP
          M.BOX <- min((upper-par)[dP>0]/dP[dP>0],(lower-par)[dP<0]/dP[dP<0])
          
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
            P1 <- cbind(P.ALL[,MIN-2])
            F1 <- F.ALL[MIN-2]
            P2 <- cbind(P.ALL[,MIN-1])
            F2 <- F.ALL[MIN-1]
            
            DIFF <- QuadSolve(P1,P2,DIR[,BEST,drop=FALSE],F1,F2)
            GRAD <- DIFF$GRAD
            H <- DIFF$HESS
            
            if(H>0) { M <- -GRAD/H }
            else { M <- max(M,last(diff(SEQ))) } # start with the last step size to keep accelerating if necessary
            
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
      # we improved our estimate and can stop the line search
      
      # interpolate to an even better point
      if(fast)
      {
        # take best triplet
        if(MIN>1) # first point after boundary
        {
          P1 <- cbind(P.ALL[,MIN-1])
          F1 <- F.ALL[MIN-1]
        }
        else # first point reflected... I don't think this can ever happen
        {
          P1 <- cbind(P.ALL[,MIN+2])
          F1 <- F.ALL[MIN+2]
        }
        
        if(MIN<length(F.ALL)) # second point before boundary
        {
          P2 <- cbind(P.ALL[,MIN+1])
          F2 <- F.ALL[MIN+1]
        }
        else # second point reflected off boundary
        {
          P2 <- cbind(P.ALL[,MIN-2])
          F2 <- F.ALL[MIN-2]
        }
        
        # calculate gradient and curvature
        DIFF <- QuadSolve(P1,P2,DIR[,BEST,drop=FALSE],F1,F2)
        GRAD <- DIFF$GRAD
        HESS.L <- DIFF$HESS
        # will update search direction hessian with this during next iteration
        
        # estimate better location if we can
        if(HESS.L>0 && MIN<length(F.ALL) && 1<MIN)
        { 
          dP <- -(GRAD/HESS.L)
          par <- line.boxer(dP*DIR[,BEST]) # need to differentiate new TARGET from evaluated par to save !!!
          
          # what is the relative error along the search line
          ERROR.L <- sqrt(((par-TARGET)%*%DIR[,BEST])^2/HESS.L)
          RATE.L <- ERROR.L/ERROR^2
          # is line search relatively efficient?
          if((RATE.L*ERROR^2)^COST <= RATE*ERROR^2/sqrt(DIM)) { SUCCESS <- FALSE }
        }
        else
        { SUCCESS <- FALSE }
      }
    } # go back to differentiation
  } # end main loop
  

  if(counts<maxit) { convergence <- 0} else { convergence <- 1 }
  
  # return stuff in similar format to optim
  RETURN <- list()
  RETURN$par <- par*parscale
  RETURN$value <- FN*fnscale
  RETURN$counts <- counts
  RETURN$convergence <- convergence
  RETURN$hessian <- t(t(HESS/parscale)/parscale)
  RETURN$covariance <- t(t(COV*parscale)*parscale)
  RETURN$lower <- lower*parscale
  RETURN$upper <- upper*parscale
      
  return(RETURN)
}