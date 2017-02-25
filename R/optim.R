
#################################
# simplex & quasi-Newton optimizers for multiple CPU cores
mcoptim <- function(par,fn,gr=NULL,method="Newton",lower=-Inf,upper=Inf,control=list(),hessian=FALSE)
{
  method <- match.arg(method,c("Newton","Simplex"))
  
  # check complains about visible bindings
  fnscale <- parscale <- ndeps <- maxit <- abstol <- reltol <- trace <- mc.cores <- NULL
  # fix default control arguments
  default <- list(fnscale=1,parscale=pmin(abs(par),abs(par-lower),abs(upper-par)),ndeps=1e-3,maxit=100,abstol=0,reltol=sqrt(.Machine$double.eps),trace=FALSE,mc.cores=parallel::detectCores())
  control <- replace(default,names(control),control)
  # check does not like attach
  NAMES <- names(control)
  lapply(1:length(control),function(i){ assign(NAMES[i],control[[i]]) })
  
  if(any(parscale==0)) { parscale[parscale==0] <- 1 }
  
  # what we will actually be evaluating
  par <- par/parscale
  func <- function(par) { fn(parscale*par)/fnscale }
  
  DIM <- length(par)
  parscale <- array(parscale,DIM)
  lower <- array(lower,DIM)/parscale
  upper <- array(upper,DIM)/parscale
  
  # prototype Hessian and inverse Hessian matrices relative to parscale
  STD <- array(1,DIM) # diagonal standard deivations
  COV <- diag(DIM) # inverse hessian
  HESS <- diag(DIM) # hesssian
  GRAD <- array(0,DIM) # gradient
  
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
    # calculate overlaps
    OVER <- colSums(P*DIR)
    # replace dimension of largest overlap with P
    MAX <- which.max(abs(OVER))
    # OVER[MAX,] <- P
    # make sure P is in its most canonical slot
    SWAP <- which.max(abs(P))
    Q <- OVER[SWAP,]
    OVER[SWAP,] <- P
    OVER[MAX,] <- Q
    
    # subtract projection from all but P direction
    DIR[-SWAP,] <- DIR[-SWAP,] - (OVER[-SWAP] %o% P)
    
    # re-normalize just in case
    MAG <- sqrt(colSums(DIR^2))
    DIR <- t(t(DIR)/MAG)
    
    return(DIR)
  }
  
  ######################
  # MAIN LOOP
  ######################
  
  ################################
  # STEP 1: O(DIM) differentiation
  ################################
  # we differentiate in DIM current directions DIR from point par
  
  # what points are on the boundary, so that we have to do one-sided derivatives
  LO <- (par <= lower)
  UP <- (par >= upper)
  BOX < UP - LO

  # rotate coordinates so that BOXed dimensions are fixed/isolated in DIR and DIR[,BOX] points strictly away from boundary
  for(i in which(BOX))
  {
    box.dir <- array(0,DIM)
    box.dir[i] <- -BOX[i]
    DIR <- Gram.Schmidt(DIR,box.dir)
  }
  
  # sample initial points surrounding the initial point for numerical differentiation
  dP <- (ndeps*ERR)*t(STD*t(DIR))
  P1 <- +dP + par # away from boundaries
  P2 <- -dP + par # towards boundaries
  # don't sample points across boundary, fold back instead - correct one-sided derivatives implemented
  if(any(BOX)) { P2[,BOX] <- +dP[,BOX]/2 + par }
  
  # apply these displacements within box constraints
  P1 <- apply(P1,2,line.boxer)
  P2 <- apply(P2,2,line.boxer)
  # columns are now boxed coordinates

  # mc evaluate all points
  P <- c(list(par),apply(P1,2,list),apply(P2,2,list))
  FN <- unlist(parallelsugar::mclapply(P,func,mc.cores=mc.cores))
  # separate into parts
  F1 <- FN[2:(DIM+1)]
  F2 <- FN[(DIM+2):(1+2*DIM)]
  FN <- FN[1]
  
  # convert back to displacements
  P1 <- P1 - par
  P2 <- P2 - par
  # signed magnitudes of displacement vectors
  M1 <- colSums(DIR * P1)
  M2 <- colSums(DIR * P2)
  
  F1 <- F1-FN
  F2 <- F2-FN
  G1 <- F1/M1
  G2 <- F2/M2
  # Hessian estimates, not necessarily assuming that par is between P1 & P2
  H <- (G2-G1)/(M2-M1)*2
  
  # gradient estimates
  GRAD <- (G2/M2-G1/M1)/(1/M2-1/M1)
  
  # transform Hessian to current coordinates
  HESS <- t(DIR) %*% HESS %*% DIR
  COV <- t(DIR) %*% COV %*% DIR
  # curvature correction factor: new curvature / old curvature
  COR <- sqrt(abs(H/diag(HESS))) # does this flip the correlations? !!!
  COR <- array(COR,c(DIM,DIM))
  COR <- t(COR) * COR
  # update curvatures while preserving correlations
  HESS <- COR * HESS
  # update covariances the same way as Hessian (prevents requirement of matrix inversion)
  COR <- sqrt(abs(diag(HESS)/H))
  COR <- array(COR,c(DIM,DIM))
  COR <- t(COR) * COR
  COV <- COR * COV

  # standard deviations along DIR axes
  STD <- sqrt(abs(diag(COV)))
  
  # transform gradient to canonical coordinates
  GRAD <- c(DIR %*% GRAD)
  
  # transform Hessian back to canonical coordinates
  HESS <- D %*% HESS %*% t(D)
  COV <- D %*% COV %*% t(D)
  
  # Newton-Raphson search step
  dP <- -(COV %*% GRAD)
  # don't search past boundary, but along the boundary
  if(BOX && GRAD[BOX,BOX]<0)
  {
    dP[BOX] <- 0
    dP[-BOX] <- -qr.solve(HESS[-BOX,-BOX],GRAD[-BOX])
  }
  BEST <- P + dP
  
  # generate a sequence of points from old par to the other side of new par
  P <- line.boxer(2*dP)
  dP <- P - par
  M <- sqrt(sum(dP^2))
  d <- dP/M
  M <- seq(0,M,length.out=min(3,mc.cores+1))[-1] * M
  dP <- (d %o% M)
  P <- dP + par
  
  # evaluate objective function
  P.ALL <- cbind(par)
  F.ALL <- FN
  
  P <- apply(P,2,list)
  FN <- unlist(parallelsugar::mclapply(P,func,mc.cores=mc.cores))
  
  # pick best points including original par
  P.ALL <- cbind(P.ALL,P)
  F.ALL <- c(F.ALL,FN)
  
  # do we need to keep going to capture the minimum?
  MIN <- which.min(F.ALL)
  if(MIN==1) # we went too far
  {
    
  }
  else if(MIN==length(F.ALL)) # we didn't go far enough
  {
    
  }
  else # we improved our esitmate
  {
    
  }
  
  # update relative ERRor
  
  # test for stopping condition

  # shrink steps
  
  # accelerated sequence of M geometrically?
  
  # construct new coordinate system
  
  # sort old directions by overlap with search direction
  
  # Gram–Schmidt from least overlap to most overlap, replacing most overlap with search direction
  
  # go back to differentiating
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}