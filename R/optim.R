
#################################
# Nelder-Mead-ish & quasi-Newton parallel optimizer for multiple CPU cores
optim_pNMBFGS <- function(par,fn,lower=-Inf,upper=Inf,control=list(),hessian=FALSE)
{
  # CRAN DOES NOT LIKE ATTACH :(
  if(length(control))
  {
    NAMES <- names(control)
    for(i in 1:length(control)) { assign(NAMES[i],control[[i]]) }
  }
  
  if(is.null(control$fnscale)) { fnscale <- 1 }
  if(is.null(control$parscale)) { parscale <- abs(par) }
  if(any(parscale==0)) { parscale[parscale==0] <- 1 }
  if(is.null(control$maxit)) { maxit <- 500 }
  if(is.null(control$abstol)) { abstol <- 0 }
  if(is.null(control$reltol)) { reltol <- sqrt(.Machine$double.eps) }
  if(is.null(control$alpha)) { alpha <- 1 }
  if(is.null(control$gamma)) { gamma <- 2 }
  if(is.null(control$rho)) { rho <- 1/2 }
  if(is.null(control$sigma)) { sigma <- 1/2 }
  if(is.null(control$trace)) { trace <- FALSE }
  if(is.null(control$mc.cores)) { mc.cores <- parallel::detectCores() }
  
  DIM <- length(par)
  lower <- array(lower,DIM)
  upper <- array(upper,DIM)
  parscale <- array(parscale,DIM)
  
  # apply box constraints to parameters
  boxer <- function(par)
  {
    par <- pmax(lower,par)
    par <- pmin(upper,par)
    return(par)
  }
  
  # invent some initial points surrounding the initial guess
  P1 <- lapply(1:DIM,function(i){par})
  P2 <- P1
  for(i in 1:DIM)
  {
    P1[[i]][i] <- P1[[i]][i] - parscale[i]/2
    P2[[i]][i] <- P2[[i]][i] + parscale[i]/2
  }
  # apply box constraints
  P1 <- lapply(P1,function(par){boxer(par)})
  P2 <- lapply(P2,function(par){boxer(par)})
  
  # cat all points to evaluate
  P <- c(list(par),P1,P2)
  # mc evaluate
  FN <- unlist(parallelsugar::mclapply(P,fn,mc.cores=mc.cores))/fnscale
  # separate into parts
  F1 <- FN[2:(DIM+1)]
  F2 <- FN[(DIM+2):(1+2*DIM)]
  
  # estimate initial Hessian matrix from this sample (diagonal)
  HESS <- diag(1,nrow=DIM)
  for(i in 1:DIM) { HESS[i,i] <- (F2[i]-2*FN[1]+F1[i])/(parscale[i]/2)^2 }

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