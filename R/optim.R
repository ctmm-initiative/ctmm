
#################################
# simplex & quasi-Newton optimizers for multiple CPU cores
mcoptim <- function(par,fn,method="Newton",lower=-Inf,upper=Inf,control=list(),hessian=FALSE)
{
  method <- match.arg(method,c("Newton","Simplex"))
  # method simplex not yet supported
  
  # check complains about visible bindings
  fnscale <- parscale <- ndeps <- maxit <- reltol <- trace <- fast <- mc.cores <- hessian <- covariance <- NULL
  # fix default control arguments !!! add a covariance/hessian argument that takes preference over parscale
  default <- list(fnscale=1,parscale=pmin(abs(par),abs(par-lower),abs(upper-par)),ndeps=1e-3,maxit=100,reltol=sqrt(.Machine$double.eps),trace=FALSE,fast=TRUE,mc.cores=parallel::detectCores(),hessian=NULL,covariance=NULL)
  control <- replace(default,names(control),control)
  # check does not like attach
  NAMES <- names(control)
  lapply(1:length(control),function(i){ assign(NAMES[i],control[[i]]) })
  
  if(any(parscale==0)) { parscale[parscale==0] <- 1 }
  
  # what we will actually be evaluating
  par <- par/parscale
  func <- function(par) { fn(par*parscale)/fnscale }
  
  DIM <- length(par)
  lower <- array(lower,DIM)/parscale
  upper <- array(upper,DIM)/parscale
  
  # prototype Hessian and inverse Hessian matrices relative to parscale
  STD <- array(1,DIM) # diagonal standard deivations
  COV <- diag(DIM) # inverse hessian
  HESS <- diag(DIM) # hesssian
  GRAD <- array(0,DIM) # gradient

  # do we have better information than parscale?
  # these are all post-fnscale parameters
  if(!is.null(covariance))
  {
    COV <- t(t(covariance/parscale)/parscale)
    STD <- sqrt(diag(COV))
    if(is.null(hessian)) { HESS <- PDsolve(COV) }
  }
  if(!is.null(hessian))
  {
    HESS <- t(t(hessian*parscale)*parscale)
    if(is.null(covariance))
    {
      COV <- PDsolve(HESS)
      STD <- sqrt(diag(COV))
    }
  }

  # start with the canonical directions
  DIR <- diag(1,DIM) # rows of column vectors
  # relative error initially assumed along directions
  ERR <- 1
  
  
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
  Gram.Schmidt <- function(DIR,P)
  {
    # normalize P
    P <- P / sqrt(sum(P^2))
    
    # calculate overlaps
    OVER <- colSums(P*DIR)
    # replace dimension of largest overlap with P
    MAX <- which.max(abs(OVER))
    DIR[MAX,] <- P
    OVER[MAX] <- 1
    # OVER[MAX,] <- P
    # make sure P is in its most canonical slot
    SWAP <- which.max(abs(P))
    DIR[c(MAX,SWAP),] <- DIR[c(SWAP,MAX),]
    OVER[c(MAX,SWAP)] <- OVER[c(SWAP,MAX)]
    
    # subtract projection from all but P direction
    MAX <- SWAP
    DIR[-MAX,] <- DIR[-MAX,] - (OVER[-MAX] %o% P)
    
    # re-normalize just in case
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
  
  # H <- 1 # search direciton hessian
  counts <- 0
  ######################
  # MAIN LOOP
  ######################
  while(ERROR >= reltol)
  {
    counts <- counts + 1
    
    ################################
    # STEP 1: O(DIM) differentiation
    ################################
    # we differentiate in DIM current directions DIR from point par
    
    # what points are on the boundary, so that we have to do one-sided derivatives
    LO <- (par <= lower)
    UP <- (par >= upper)
    BOX < UP - LO
    
    # rotate coordinates so that BOXed dimensions are fixed/isolated in DIR and DIR[,BOX] points strictly towards boundary
    for(i in which(BOX))
    {
      dir <- array(0,DIM)
      dir[i] <- BOX[i]
      DIR <- Gram.Schmidt(DIR,dir)
    }
    # now fix BOX to be TRUE/FALSE
    BOX <- as.logical(BOX)
    
    # sample initial points surrounding the initial point for numerical differentiation
    dP <- (ndeps*ERR)*t(STD*t(DIR))
    P1 <- -dP + par # away from boundaries
    P2 <- +dP + par # towards boundaries
    # don't sample points across boundary, fold back instead - correct one-sided derivatives implemented
    if(any(BOX)) { P2[,BOX] <- -dP[,BOX]/2 + par }
    
    # apply these displacements within box constraints
    P1 <- apply(P1,2,line.boxer)
    P2 <- apply(P2,2,line.boxer)
    # columns are now boxed coordinates
    
    # mc evaluate all points
    P <- c(list(par),apply(P1,2,list),apply(P2,2,list))
    FN <- unlist(parallelsugar::mclapply(P,func,mc.cores=mc.cores))
    # separate back into parts
    F1 <- FN[2:(DIM+1)]
    F2 <- FN[(DIM+2):(1+2*DIM)]
    fp <- FN[1] # fn(par)
    
    # calculate axial derivatives to second order
    DIFF <- QuadSolve(P1,P2,DIR,F1,F2)
    GRAD <- DIFF$GRAD
    H <- DIFF$HESS
    
    # transform Hessian to current coordinates
    HESS <- t(DIR) %*% HESS %*% DIR
    COV <- t(DIR) %*% COV %*% DIR
    # curvature correction factor: new curvature / old curvature (current coordinates)
    H <- sqrt(abs(H/diag(HESS)))
    COR <- array(H,c(DIM,DIM))
    COR <- t(COR) * COR
    # update curvatures while preserving correlations
    HESS <- COR * HESS
    # update covariances the same way as Hessian (prevents requirement of matrix inversion)
    COR <- array(1/H,c(DIM,DIM))
    COR <- t(COR) * COR
    COV <- COR * COV
    
    # standard deviations along DIR axes
    STD <- sqrt(abs(diag(COV)))
    
    # transform gradient to canonical coordinates
    GRAD <- c(DIR %*% GRAD)
    
    # transform Hessian back to canonical coordinates
    HESS <- DIR %*% HESS %*% t(DIR)
    COV <- DIR %*% COV %*% t(DIR)
    
    # Newton-Raphson search step
    dP <- -(COV %*% GRAD)
    # don't search past boundary, but along the boundary
    if(BOX && GRAD[BOX]>0)
    {
      dP[BOX] <- 0
      dP[-BOX] <- -PDsolve(HESS[-BOX,-BOX],GRAD[-BOX])
    }
    
    # test for stopping condition here
    ERROR <- max(dP/STD)
    if(ERROR < reltol) { break } # we are finished here
    
    # calculate best DIR and orthonormalize the remaining DIRs
    BEST <- which.max(abs(dP)) # future index of dP direction
    DIR <- Gram.Schmidt(DIR,dP)
    
    ##################
    # LINE SEARCH LOOP
    ##################
    counts <- counts + 1
    
    # generate a linear sequence of points from old par to the other side of new par
    P <- line.boxer(2*dP)
    dP <- P - par
    M <- sqrt(sum(dP^2)) # total search magnitude
    SEQ <- seq(0,M,length.out=min(2,mc.cores)+1)[-1]
    P <- (DIR[,BEST] %o% SEQ) + par
    
    # initialize (all) storage results
    P.ALL <- cbind(par)
    F.ALL <- fp
    
    # start iteration loop
    MIN <- 1
    while(TRUE)
    {
      # most expensive part
      # evaluate objective function at new P and store to FN
      P <- apply(P,2,list)
      FN <- unlist(parallelsugar::mclapply(P,func,mc.cores=mc.cores))
      
      # combine with older results
      P.ALL <- cbind(P.ALL,P)
      F.ALL <- c(F.ALL,FN)
      
      # sort along DIR[,BEST]
      SORT <- colSums(P.ALL*DIR[,BEST])
      SORT <- sort(SORT,method="quick",index.return=TRUE)$ix
      P.ALL <- P.ALL[,SORT]
      F.ALL <- F.ALL[SORT]
      
      # new best estimate
      MIN <- which.min(F.ALL)
      par <- P.ALL[,MIN]
      fp <- F.ALL[MIN]
      
      # do we need to keep going to capture the minimum?
      
      if(MIN==1) # we went too far in first step
      {
        # do a sub line search, including backwards reflection
        # what was the shortest step previously
        M <- SEQ[1] * M
        # tabulate forward steps with half of cores (assuming multiple of 2)
        SEQ <- seq(0,M,length.out=min(1,mc.cores/2)+2)
        SEQ <- SEQ[-c(1,length(SEQ))]
        # tabulate backward steps other half of cores
        SEQ <- c(rev(seq(0,-M,length.out=min(1,mc.cores/2)+1)[-1]),SEQ)
        
        # combine for evaluation
        P <- (DIR[,BEST] %o% SEQ) + par
        
        # goto evaluate iteration step
        next
        
        # if !fast, this might not be an improvement...
        # need a conditional or fast by default
      }
      else if(MIN==length(F.ALL)) # we didn't go far enough or we hit a boundary
      {
        # if we hit a boundary, then we can stop at the boundary
        NEW.BOX <- (par <= lower) || (par >= upper)
        if(sum(NEW.BOX)>sum(BOX)) { break } # goto hessian stuff
        
        # Distance to first boundary that we will hit going in direction dP
        M.BOX <- min((upper-par)[dP>0]/dP[dP>0],(lower-par)[dP<0]/dP[dP<0])
        
        # we didn't hit a boundary yet, but we will eventually hit a boundary, so let's just do that
        if(M.BOX<Inf)
        {
          # generate a linear sequence of points that terminate at the eventual boundary
          SEQ <- seq(0,M.BOX,length.out=min(2,mc.cores)+1)[-1]
          P <- (DIR[,BEST] %o% SEQ) + par
          
          # goto evaluate iteration step
          next
        }
        else # we didn't hit a boundary and we never will, because there is no boundary
        {
          # generate an exponential sequence of points that terminate at the infinite boundary but drop the last point
          # start with a uniform sequence
          SEQ <- seq(0,1,length.out=min(2,mc.cores)+2)
          SEQ <- t[-c(1,length(SEQ))]
          # log transform
          SEQ <- M * log(1-SEQ)
          P <- (DIR[,BEST] %o% SEQ) + par
          
          # update M bigly to keep accelerating next time around (if necessary)
          M <- last(diff(SEQ)) * M
          
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
        P1 <- cbind(P[MIN-1])
        F1 <- F.ALL[MIN-1]
      }
      else # first point reflected... I don't think this can ever happen
      {
        P1 <- cbind(P[MIN+2])
        F1 <- F.ALL[MIN+2]
      }
      
      if(MIN<length(F.ALL)) # second point before boundary
      {
        P2 <- cbind(P[MIN+1])
        F2 <- F.ALL[MIN+1]
      }
      else # second point reflected off boundary
      {
        P2 <- cbind(P[MIN-2])
        F2 <- F.ALL[MIN-2]
      }
      
      # calculate gradient and curvature
      DIFF <- QuadSolve(P1,P2,DIR[,BEST,drop=FALSE],F1,F2)
      GRAD <- DIFF$GRAD
      H <- DIFF$HESS
      # will update search direction hessian with this during next iteration
      
      # estimate better location if we can
      if(H>0 && MIN<length(F.ALL)) { par <- line.boxer(par - (GRAD/H)*DIR[,BEST]) }
    }
    
    # go back to differentiation
  }  

  if(counts==maxit) { convergence <- 1} else { convergence <- 0 }
  
  # return stuff in similar format to optim
  RETURN <- list()
  RETURN$par <- par
  RETURN$value <- FN*fnscale
  RETURN$counts <- counts
  RETURN$convergence <- convergence
  RETURN$hessian <- t(t(HESS/parscale)/parscale)
  RETURN$covariance <- t(t(COV*parscale)*parscale)
  RETURN$lower <- lower
  RETURN$upper <- upper
      
  return(RETURN)
}