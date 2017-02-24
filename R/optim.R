
#################################
# simplex & quasi-Newton optimizers for multiple CPU cores
mcoptim <- function(par,fn,gr=NULL,method="Newton",lower=-Inf,upper=Inf,control=list(),hessian=FALSE)
{
  method <- match.arg(method,c("Newton","Simplex"))
  
  # fix default control arguments
  default <- list(fnscale=1,parscale=abs(par),ndeps=1e-3,maxit=100,abstol=0,reltol=sqrt(.Machine$double.eps),trace=FALSE,mc.cores=parallel::detectCores())
  control <- replace(default,names(control),control)
  
  # CRAN DOES NOT LIKE ATTACH :(
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
  
  # prototype Hessian and inverse Hessian matrices
  STD <- parscale # diagonal standard deivations
  COV <- diag(parscale^2) # inverse hessian
  HESS <- diag(1/parscale^2) # hesssian
  GRAD <- array(0,DIM) # gradient
  ERR <- 1 # relative error
  
  ######################
  # apply box constraints to travel in a straight line from (global) par by dp
  ######################
  line.boxer <- function(dp)
  {
    # if we don't hit a boundary we go here
    p <- par + dp
    
    # did we hit a boundary?
    UP <- (p >= upper)
    LO <- (p <= lower)
    if(any(UP) || any(LO)) { BOX <- TRUE } else { BOX <- FALSE }
    
    # time until we hit the first boundary
    if(any(UP)) { t.up <- (upper-par)[UP]/dp[UP] } else { t.up <- 1 }
    if(any(LO)) { t.lo <- (lower-par)[LO]/dp[LO] } else { t.lo <- 1 }
    t <- min(t.lo,t.up)
    
    # stop at first boundary
    p <- par + t*dp
      
    return(p)
  }
  
  ##################
  # two-sided symmetric derivatives
  ##################
  # invent some initial points surrounding the initial guess (by parscale)
  P1 <- array(0,c(DIM,DIM))
  P2 <- P1
  for(i in 1:DIM)
  {
    # set of orthogonal displacements
    P1[i,i] <- + ndeps*STD*ERR
    P2[i,i] <- - ndeps*STD*ERR
    # if par is on a boundary then we need to generalize this so that both P1 & P2 are in the same direction !!!
  }
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
  # magnitudes of displacement vectors
  # M1 <- sqrt(colSums(P1^2))
  M2 <- sqrt(colSums(P2^2))
  # directions of displacement vectors (defining the positive directon)
  D <- t(t(P2)/M2)
  # R applies addition and multiplication to different margins by default... weird
  # consistent magnitudes in the same direction (general purpose code)
  M1 <- colSums(P1*D)
  
  F1 <- F1-FN
  F2 <- F2-FN
  G1 <- F1/M1
  G2 <- F2/M2
  # Hessian estimates in the direction of D, not necessarily assuming that par is between P1 & P2
  H <- (G2-G1)/(M2-M1)*2
  # make sure we are always going downhill

  # gradient estimates
  GRAD <- (G2/M2-G1/M1)/(1/M2-1/M1)

  # transform gradient to canonical coordinates
  GRAD <- c(D %*% GRAD)
  
  # transform Hessian to current coordinates
  HESS <- t(D) %*% HESS %*% D
  COV <- t(D) %*% COV %*% D
  # curvature correction factor: new curvature / old curvature
  COR <- sqrt(H/diag(HESS))
  COR <- array(COR,c(DIM,DIM))
  COR <- t(COR) * COR
  # update curvatures while preserving correlations
  HESS <- COR * HESS
  # update covariances the same way as Hessian (prevents requirement of matrix inversion)
  COR <- sqrt(diag(HESS)/H)
  COR <- array(COR,c(DIM,DIM))
  COR <- t(COR) * COR
  COV <- COR * COV

  # standard deviations
  
  # transform Hessian back to canonical coordinates
  HESS <- D %*% HESS %*% t(D)
  COV <- D %*% COV %*% t(D)
  
  # Newton-Raphson search step
  dP <- -(COV %*% GRAD)
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # take best points avoiding degeneracy
  P <- P[1:(DIM+1)]
  FN <- FN[1:(DIM+1)]
  for(i in 1:DIM)
  {
    if(F1[i]<F2[i])
    {
      P[[i+1]] <- P1[[i]]
      FN[i+1] <- F1[i]
    }
    else
    {
      P[[i+1]] <- P2[[i]]
      FN[i+1] <- F2[i]
    }
  }
  rm(P1,P2,F1,F2)
  
  # sort by objective function
  sorter <- function()
  {
    FN <<- sort(FN,index.return=TRUE)
    P <<- P[FN$ix]
    FN <<- FN$x
  }
  sorter()
  
  # decompose vector in to magnitude and direction
  decompose <- function(par)
  {
    MAG <- sqrt(sum(par^2))
    DIR <- par/MAG
    return(list(MAG=MAG,DIR=DIR))
  }
  
  # 1/vector for gradient estimation
  reciprocate <- function(par)
  {
    par <- decompose(par)
    return(par$DIR/par$MAG)
  }
  
  # Newton iterator that can handle negative curvatures appropriately
  boxsolver <- function(HESS,GRAD,P)
  {
    P - qr.solve(HESS,GRAD)
  }
  
  ################
  # main loop
  TOL <- Inf
  counts <- 0
  NAMES <- c("centroid","reflection","expansion","contraction")
  while(TOL >= abstol && TOL/FN[1] >= reltol && counts < maxit)
  {
    # new points to evaluate
    PN <- vector("list",4)
    # centroid ######## seems like we can use this guy more ###########
    PN[[1]] <- Reduce('+',P[1:DIM])/DIM
    # reflection
    PN[[2]] <- PN[[1]] + alpha*(PN[[1]]-P[[DIM+1]])
    # expansion
    PN[[3]] <- PN[[1]] + gamma*(PN[[2]]-PN[[1]])
    # contraction
    PN[[4]] <- PN[[1]] + rho*(P[[DIM+1]]-PN[[1]])
    
    # fill up remaining cue with approximate Newton search
    # approximate gradient from simplex
    GRAD <- lapply(2:(DIM+1),function(i){ (FN[i]-FN[1])*reciprocate(P[[i]]-P[[1]]) })
    GRAD <- Reduce('+',GRAD)/DIM
    # solve for the Newton step accouting for negative curvature == boundary solution
    GRAD <- boxsolver(HESS,GRAD,P[[1]])
    PN[[5]] <- GRAD$P
    # sequence along the Newton search that we will sample to improve gradient/hessian in that direction
    SEQ <- seq(0,2,length.out=(max(mc.cores,6)-5)+2)
    SEQ <- SEQ[2:length(SEQ)]
    SEQ <- SEQ[SEQ!=1]
    for(i in 6:max(mc.cores,6))
    { PN[[i]] <- P[[1]] + SEQ[i-5]*GRAD$dP }

    # evaluate new points in mc
    PN <- lapply(PN,function(par){boxer(par)})
    EN <- unlist(parallelsugar::mclapply(PN,fn,mc.cores=mc.cores))/fnscale

    # update Hessian/Covariance with NM line search result
    MAX <- which.max(EN[1:4]) # worst point to toss can't be in the middle
    if(MAX==1 || MAX==4)
    {
      dP <- decompose(PN[[1]]-P[[DIM+1]])
      DIR <- dP$DIR
      # Nelder-Mead search lengths
      MAG <- c(-rho,0,alpha,alpha*gamma) * dP$MAG
      MAG <- MAG[-MAX]
      # gradients along this direction
      GRAD <- diff(EN[1:4][-MAX])/diff(MAG)
      # Hessian along this direction
      H <- diff(GRAD)/(MAG[3]-MAG[1])*2
      # remove old variance & fix new variance
      HESS <- HESS + (H - (DIR %*% HESS %*% DIR))*(DIR %o% DIR)
    }
    
    # update Hessian/Covariance with Newton line search result
    #
    
    # to avoid degeneracy, we can only take the best of P[1:DIM] and centroid
    # why isn't this usually done?
    
    # to avoid degeneracy, we can only take the best of reflection, expansion, contraction
    BEST <- sort(EN[2:4],index.return=TRUE)$ix[1] + 1
    
    # to avoid degeneracy, we can only take the 2 best of P[[1]] and Newton direction
    
    # shrinkage step if nothing is found
    if(EN[BEST] >= FN[DIM+1])
    {
      NAME <- "shrink"
      
      P[2:(DIM+1)] <- lapply(P[2:(DIM+1)],function(par){P[[1]] + sigma*(par-P[[1]])})
      FN[2:(DIM+1)] <- unlist(parallelsugar::mclapply(P[2:(DIM+1)],fn,mc.cores=mc.cores))/fnscale
    }
    else
    {
      NAME <- NAMES[BEST]
      
      P <- c(P,PN[BEST])
      FN <- c(FN,EN[BEST])
    }
    sorter()
    P <- P[1:(DIM+1)]
    FN <- FN[1:(DIM+1)]
    
    TOL <- FN[DIM+1]-FN[1]
    counts <- counts + 1
    
    if(trace) { message("step ",counts,"\t",NAME,"\t",TOL) }
  }
  # now we can include the centroid
  P <- c(P,PN[1])
  FN <- c(FN,EN[1])
  sorter()
  
  RETURN <- list()
  RETURN$par <- P[1]
  RETURN$value <- FN[1]
  RETURN$counts <- counts
  if(counts<maxit) { RETURN$convergence <- 0 }
  else { RETURN$convergence <- 1 }
  # return 10 if degeneracy in P[1:DIM]
  
  return(RETURN)
}