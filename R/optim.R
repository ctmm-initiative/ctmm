###################################
# make a wrapper that applies optim then afterwards numDeriv, possibly on a boundary with one-sided derivatives if necessary
Optimizer <- function(par,fn,...,method="Nelder-Mead",lower=-Inf,upper=Inf,period=F,control=list())
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
  { RESULT <- mc.optim(par=par,fn=fn,lower=lower,upper=upper,period=period,control=control) }
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
    control$zero <- NULL
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
line.boxer <- function(dp,p0=dp[,1],lower=-Inf,upper=Inf,period=F,period.max=1/2)
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
    if(any(LO)) { LO <- LO & (dp<0) }
    if(any(UP)) { UP <- UP & (dp>0) }
    
    # stop when we hit the boundary  
    # time until we hit the first new boundary
    if(any(LO)) { t.lo <- (lower-p0)[LO]/dp[LO] } else { t.lo <- 1 }
    if(any(UP)) { t.up <- (upper-p0)[UP]/dp[UP] } else { t.up <- 1 }
    t <- min(t.lo,t.up)
    
    # don't go more than period.max fraction of a period in one step
    PERIOD <- as.logical(period)
    if(any(PERIOD)) { t <- min(t,abs(period/dp)[PERIOD]*period.max) }
    
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
  DIR <- cbind(DIR)
  M1 <- colSums(DIR * P1)
  M2 <- colSums(DIR * P2)
  
  F1 <- F1-F0
  F2 <- F2-F0
  
  G1 <- F1/M1
  G2 <- F2/M2
  # Hessian estimates, not necessarily assuming that par is between P1 & P2
  hessian <- (G2-G1)/(M2-M1)*2
  
  # gradient estimates, also not assuming that par is between P1 & P2
  gradient <- (G2/M2-G1/M1)/(1/M2-1/M1)
  
  return(list(gradient=gradient,hessian=hessian)) 
}

# does this data look roughly quadratic near the minimum?
QuadTest <- function(x,y,thresh=0.5)
{
  MIN <- which.min(y)
  DIFF <- QuadSolve(x[MIN],x[MIN-1],x[MIN+1],1,y[MIN],y[MIN-1],y[MIN+1])
  # predicted minimum
  x0 <- x[MIN] - DIFF$gradient/DIFF$hessian
  y0 <- y[MIN] - DIFF$hessian/2*(x[MIN]-x0)^2
  # quadratic estimate
  y.func <- function(x) { y0 + DIFF$hessian/2*(x-x0)^2 }
  r <- y.func(x[MIN+c(-2,2)])
  # relative residuals
  r <- (y[MIN+c(-2,2)]-r)/(r-y0)
  return(all(abs(r)<thresh))
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
  
  # subtract projection from all but vec
  DIR[,-MAX] <- DIR[,-MAX] - (vec %o% OVER[-MAX])
  
  # re-normalization (may be necessary?)
  MAG <- sqrt(colSums(DIR^2))
  DIR <- t(t(DIR)/MAG)
  
  return(DIR)
}


# correlation-preserving rank-1 update
rank1update <- function(H.LINE,LINE,hessian,covariance)
{
  LINE <- LINE/sqrt(sum(LINE^2))
  
  # Hessian diagonal element correction factor
  H0 <- c(LINE %*% hessian %*% LINE)
  FACT <- abs(H.LINE / H0)
  FACT <- sqrt(FACT)
  
  # rank-1 Hessian update - no matrix-matrix multiplication
  A <- FACT-1
  B <- c(hessian %*% LINE)
  O <- outer(LINE)
  hessian <- hessian + A^2*H0*O
  B <- outer(LINE,B)
  B <- B + t(B)
  hessian <- hessian + A*B
  
  # rank-1 covariance update - no matrix-matrix multiplication
  A <- 1/FACT-1
  B <- c(covariance %*% LINE)
  covariance <- covariance + A^2*c(LINE%*%B)*O
  B <- outer(LINE,B)
  B <- B + t(B)
  covariance <- covariance + A*B
  
  return(list(hessian=hessian,covariance=covariance))
}

# best number of calculations to make with min count and mc.cores
mc.min <- function(min,mc.cores=detectCores())
{
  x <- ceiling(min/mc.cores)
  x <- x * mc.cores
  return(x)
}

#################################
# quasi-Newton optimizers for multiple CPU cores
# TODO: replace search lines with steepest ascent search curves?
# TODO: filling up mc queue if p>2n+1
##################################
mc.optim <- function(par,fn,...,lower=-Inf,upper=Inf,period=F,control=list())
{
  DEBUG <- FALSE
  # check complains about visible bindings
  fnscale <- parscale <- maxit <- precision <- trace <- mc.cores <- hessian <- covariance <- zero <- NULL
  # fix default control arguments
  default <- list(fnscale=1,parscale=pmin(abs(par),abs(par-lower),abs(upper-par)),maxit=100,trace=FALSE,precision=1/2,mc.cores=detectCores(),hessian=NULL,covariance=NULL,zero=FALSE)
  control <- replace(default,names(control),control)
  # check does not like attach
  NAMES <- names(control)
  for(i in 1:length(control)) { assign(NAMES[i],control[[i]]) }
  
  if(any(parscale==0)) { parscale[parscale==0] <- 1 }
  
  # machine tolerances for various calculations
  TOL.HESS <- .Machine$double.eps^.25
  TOL.GRAD <- .Machine$double.eps^.5
  TOL.LIKE <- .Machine$double.eps
  # goal errors
  TOL.GOAL <- .Machine$double.eps^precision
  
  DIM <- length(par)
  period <- array(period,DIM)
  PERIOD <- as.logical(period)
  PERIOD <- FALSE # THIS STUFF IS NOT WORKING YET - ALSO NEEDS TO DISABLE LINE SEARCH CODE
  # fix parscale for periodic variables
  if(any(PERIOD))
  { 
    # make pi the period for local tangent transformation
    perscale <- period[PERIOD]/pi
    par[PERIOD] <- par[PERIOD]/perscale
    parscale[PERIOD] <- 1 # effective parscale fixed to 45 degrees
    # tangent limits post transformation
    lower[PERIOD] <- -Inf
    upper[PERIOD] <- Inf
  }

  par <- par/parscale
  lower <- array(lower,DIM)/parscale
  upper <- array(upper,DIM)/parscale
  period <- period/parscale
  
  # start with the canonical directions
  DIR <- diag(1,DIM) # rows of column vectors
  # DIR.theta <- rep(0,(DIM^2-DIM)/2) # rotation angles corresponding to orthogonal transformation DIR
  
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
  
  # initialize local tangent transformation
  if(any(PERIOD))
  {
    # tangent transformation centered on best angle
    angle <- par[PERIOD]
    par[PERIOD] <- 0
    # convert to pi-period scale
    covariance[PERIOD,PERIOD] <- t(t(covariance[PERIOD,PERIOD]/perscale)/perscale)
    hessian[PERIOD,PERIOD] <- t(t(hessian[PERIOD,PERIOD]*perscale)*perscale)
  }
  
  ZERO <- zero # ZERO will be boolean, zero will be constant
  zero <- 0
  TOL.ZERO <- 0 # error from zeroing
  
  # what we will actually be evaluating
  func <- function(par,...)
  { 
    # transform back to 
    if(any(PERIOD)) { par[PERIOD] <- angle + atan(par[PERIOD]) }
    # this objective function has the ability to approximately zero its objective value
    if(ZERO) { RETURN <- fn(par*parscale,zero=zero*fnscale,...)/fnscale }
    else { RETURN <- fn(par*parscale,...)/fnscale } # ordinary objective function
    return(RETURN)
  }
  
  # recenter tangent transformation on periodic par
  reset.tangent <- function()
  {
    angle.old <- angle
    # transform back to angle from par
    angle <<- angle + atan(par[PERIOD])
    # derivative of tangent wrt angle
    FACT <- 1 + par[PERIOD]^2
    # first derivative transformation
    gradient[PERIOD] <<- FACT * gradient[PERIOD]
    gradient.old[PERIOD] <<- FACT * gradient.old[PERIOD]
    # second derivative transformation
    hessian[PERIOD,PERIOD] <<- t(t(hessian[PERIOD,PERIOD]*FACT)*FACT)
    # inverse transformation
    covariance[PERIOD,PERIOD] <<- t(t(covariance[PERIOD,PERIOD]/FACT)/FACT)

    # best angle is new center of tangent transformation
    par[PERIOD] <<- 0
    par.target[PERIOD] <<- tan((atan(par.target[PERIOD]) + angle.old - angle))
    par.target.old[PERIOD] <<- tan((atan(par.target.old[PERIOD]) + angle.old - angle))
  }
  
  ##################
  # am I on a boundary and if so, align DIR
  ######################
  is.boxed <- function(par,fix.dir=FALSE)
  {
    # what points are on the boundary, so that we have to do one-sided derivatives
    LO <- (par <= lower)
    UP <- (par >= upper)
    DBOX <- UP - LO # encoded with directions to boundary
    BOX <- as.logical(DBOX) 
    
    # rotate coordinates so that BOXed dimensions are fixed/isolated in DIR and DIR[,BOX] points strictly towards boundary
    if(fix.dir)
    {
      for(i in which(BOX))
      {
        dir <- numeric(DIM)
        dir[i] <- DBOX[i]
        DIR <<- Gram.Schmidt(DIR,dir)
      }
    }

    return(BOX)
  }
  
  numderiv.diff <- function(p0,DIR,TOL.DIFF)
  {
    # different calculations for the numerical step length
    # standard deviation for each axis
    # STD <- sqrt(abs(diag(covariance)))
    # parameter scale not along current axes
    SCL <- pmax(abs(par),1)
    # squared parameter scale along axes - analog to STD
    SCL <- t(DIR) %*% diag(SCL^2,nrow=DIM) %*% DIR
    # single parameter scale for each axis
    SCL <- sqrt(abs(diag(SCL)))
    # sample initial points around the center for numerical differentiation
    STEP <- TOL.HESS*SCL
    par.step <- t(STEP*t(DIR))
    P1 <- -par.step # away from boundaries
    P2 <- +par.step # towards boundaries
    # don't sample points across boundary, fold back instead - correct one-sided derivatives implemented
    if(any(BOX)) { P2[,BOX] <- -par.step[,BOX]*2 }
    # we could potentially get cornered doing this... but seems unlikely !!!
    
    # apply these displacements within box constraints
    P1 <- line.boxer(P1,p0=p0,lower=lower,upper=upper,period=period)
    P2 <- line.boxer(P2,p0=p0,lower=lower,upper=upper,period=period)
    # columns are now boxed coordinates
    
    return(list(P1=P1,P2=P2))
  }
  
  INIT <- list()
  stage.progress <- function()
  {
    if(stage==1)
    {
      
    }
    
    stage <<- stage + 1
    
    return(INIT)
  }
  
  # initializing stuff checked in loop
  stage <- 1 # Newton, Gradient, Pattern
  counts <- 0
  LINE.DO <- TRUE # else do line searches in between Newton-Raphson iterations (always initially)
  LINE.DID <- FALSE # did we do a line search previously
  LINE.BETTER <- FALSE
  par.target <- par # where to evaluate around next
  par.target.old <- par
  fn.par <- Inf # current best objective value fn(par)
  condition <- Inf
  gradient.old <- rep(0,DIM)
  hessian.LINE <- NULL # line-search hessian
  ######################
  # MAIN LOOP
  ######################
  while(counts < maxit)
  {
    # adjust zero shift
    if(ZERO && fn.par<Inf) { zero <- zero + fn.par }
    
    BOX <- is.boxed(par.target,fix.dir=TRUE)
    
    ################################
    # QUASI-NEWTON RAPHSON
    ###############################
    if(stage==1)
    {
      TOL.STAGE <- max(TOL.GOAL,TOL.HESS)
      ## update hessian from line-search result (Rank-1 update)
      if(LINE.DID && hessian.LINE>0)
      {
        DIFF <- rank1update(hessian.LINE,DIR.STEP,hessian,covariance)
        hessian <- DIFF$hessian
        covariance <- DIFF$covariance
      }
      
      # recenter the tangent transformation on par for periodic parameters
      if(counts>0 && any(PERIOD)) { reset.tangent() }
      
      ## lots of choices here for the basis
      # we could stick with the canonical basis
      if(counts>0 && any(!BOX))
      {
        ## eigen basis
        # I could make this O(n^2) if we decide on this choice
        # DIR[!BOX,!BOX] <- eigen(hessian[!BOX,!BOX])$vectors
        
        ## random basis
        # not sure if I can make this better than O(n^3)
        n <- sum(!BOX)
        dir <- matrix(0,n,n)
        # random angles near zero
        dir[upper.tri(dir)] <- stats::runif((n^2-n)/2,min=-pi,max=pi)
        # anti-symmetric matrix generator
        dir <- dir - t(dir)
        # orthogonal matrix
        DIR[!BOX,!BOX] <- expm::expm(dir,trySym=FALSE)
      }
      
      # transform Hessian to current coordinates
      hessian <- t(DIR) %*% hessian %*% DIR
      covariance <- t(DIR) %*% covariance %*% DIR
      
      # sample for numeric differentiation
      DIFF <- numderiv.diff(par.target,DIR,TOL.HESS)
      # mc evaluate all points
      P <- cbind(par,par.target,DIFF$P1,DIFF$P2)
      # par must be re-evaluated againt to prevent zero shift roundoff error
      counts.diff <- ceiling(ncol(P)/mc.cores)
      counts <- counts + counts.diff
      fn.queue <- unlist(mclapply(split(P,col(P)),func,mc.cores=mc.cores))
      # separate back into parts
      fn.par <- fn.queue[1]
      fn.target <- fn.queue[2]
      F1 <- fn.queue[2+1:DIM]
      F2 <- fn.queue[2+DIM+1:DIM]
      
      if(counts>counts.diff && ZERO) { TOL.ZERO <- abs(fn.par) }
      
      # is the minimum encapsulated?
      encapsulated <- (fn.target<F1) & (fn.target<F2)
      
      # calculate axial derivatives to second order
      DIFF <- QuadSolve(par.target,DIFF$P1,DIFF$P2,DIR,fn.target,F1,F2)
      gradient <- DIFF$gradient
      hessian.diag <- DIFF$hessian

      # ERROR CHECKING
      #######################################
      
      # will need line search if not normal looking
      TEST <- is.nan(hessian.diag)
      if(any(TEST)) { hessian.diag[TEST] <- 0 }
      if(any(hessian.diag<=.Machine$double.eps)) { LINE.DO <- TRUE }
      
      # if the gradient is null in one direction, then we've likely profiled that axis to machine precision
      # fix profiled dimensions to change nothing, else will get 0/0
      TEST <- TEST | abs(F1-fn.target)<=.Machine$double.eps | abs(F2-fn.target)<=.Machine$double.eps
      if(any(TEST))
      {
        gradient[TEST] <- 0
        hessian.diag[TEST] <- diag(hessian)[TEST]
      }
      
      # stick with old estimate
      TEST <- abs(hessian.diag)<=.Machine$double.eps
      if(any(TEST)) { hessian.diag[TEST] <- diag(hessian)[TEST] }
      
      condition <- hessian.diag/diag(hessian)
      MIN <- which.min(abs(condition))
      MAX <- which.max(abs(condition))
      condition <- condition[c(MIN,MAX)]
      condition <- sign(condition) * sqrt(abs(condition))
      
      # transform gradient back to canonical coordinates
      gradient <- c(DIR %*% gradient)

      # update hessian from Newton-Raphson result (rank-1)
      # we do this before the better local derivative update
      # also this provides an objective check on NR step quality
      ########################################
      if(counts>counts.diff)
      {
        # search direction
        DIFF <- par.target-par.target.old
        MAG <- sqrt(sum(DIFF^2))
        DIFF <- DIFF/MAG
        
        # hessian along line between centers where derivatives were calculated
        gradient.diff <- gradient-gradient.old
        hessian.DIFF <- c( DIFF%*%gradient.diff ) / MAG
        if(is.nan(hessian.DIFF)) { hessian.DIFF <- 0 }
        
        # if we just did a line search, then we can compare the two to see if NR would have given a bad Hessian update
        if(LINE.DID)
        {
          # relative comparison between curvature estimates along search direction
          TEST <- abs(log(abs(hessian.LINE/hessian.DIFF)))
          # NR with a factor of 2 curvature estimate?
          if(is.nan(TEST) || TEST >= log(2)) { LINE.DO <- TRUE }
        }
        
        # will need line search if not normal looking
        if(hessian.DIFF<=.Machine$double.eps) { LINE.DO <- TRUE }
        
        if(!LINE.DID && fn.target<fn.par && abs(hessian.DIFF)>.Machine$double.eps)
        {
          # transform DIFF to same frame as hessian/covariance
          DIFF <- rank1update(hessian.DIFF,c(t(DIR)%*%DIFF),hessian,covariance)
          hessian <- DIFF$hessian
          covariance <- DIFF$covariance
        }
      }
      gradient.old <- gradient
      par.target.old <- par.target
      
      # now that we've updated the hessian from the search result...
      if(counts>counts.diff && fn.par<=fn.target && !LINE.DID)
      {
        # Newton-Raphson failed, so don't update the Hessian with local information from here (its worse than where we were)
        LINE.DO <- TRUE
        
        # transform Hessian back to canonical coordinates
        hessian <- DIR %*% hessian %*% t(DIR)
        covariance <- DIR %*% covariance %*% t(DIR)
        
        if(trace) { message(sprintf("%f Newton-Raphson overstep (encapsulated=%d/%d; condition=%f:%f)",zero+fn.par,sum(encapsulated),DIM,condition[1],condition[2])) }
        LINE.TYPE <- "Enclosure"
      }
      else
      {
        # DIAGONAL UPDATE (Rank-DIM)
        ####################################
        # curvature correction factor: new curvature / old curvature (current coordinates)
        FACT <- sqrt(abs(hessian.diag/diag(hessian))) # correction factor diagonal
        # update curvatures while preserving correlations
        hessian <- t(t(hessian*FACT)*FACT)
        # update covariances the same way as Hessian (prevents requirement of matrix inversion)
        covariance <- t(t(covariance/FACT)/FACT)
        
        # transform Hessian back to canonical coordinates
        hessian <- DIR %*% hessian %*% t(DIR)
        covariance <- DIR %*% covariance %*% t(DIR)
        
        # Newton-Raphson search step from par.target
        par.diff <- -c(covariance %*% gradient)
        
        # new best par
        par <- par.target
        fn.par <- fn.target

        # don't search past boundary, but along the boundary
        BOX <- is.boxed(par,fix.dir=TRUE) # par could have been updated
        if(any(BOX & (par.diff%*%DIR)>0))
        {
          par.diff[BOX] <- 0
          # boundary restricted search
          if(any(!BOX)) { par.diff[!BOX] <- -PDsolve(hessian[!BOX,!BOX]) %*% gradient[!BOX] }
        }
        
        # new search direction
        DIR.STEP <- par.diff / sqrt(sum(par.diff^2))
        
        # where we aim to evaluate next
        par.target <- line.boxer(par.diff,p0=par,lower=lower,upper=upper,period=period)

        if(trace) { message(sprintf("%f Newton-Raphson step (encapsulated=%d/%d; condition=%f:%f)",zero+fn.par,sum(encapsulated),DIM,condition[1],condition[2])) }
        LINE.TYPE <- "Prospection"
      }
      
      # NR stopping conditions here
      # alternative stopping conditions in line search
      if(!LINE.DO && !LINE.DID)
      {
        TEST <- (ZERO && abs(fn.target.old-fn.par) < TOL.STAGE + TOL.ZERO)
        TEST <- TEST || (!ZERO && abs((fn.target.old-fn.par)/fn.par) < TOL.STAGE)
        if(TEST) { INIT <- stage.progress() }
      }
      fn.target.old <- fn.target
    }
    
    if(stage==2)
    {
      #######################
      # GRADIENT DESCENT
      #######################
      TOL.STAGE <- max(TOL.GOAL,TOL.GRAD)
      
      LINE.DO <- FALSE
      
      # stop conditions
      stage <- stage + 1
    }
    
    if(stage==3)
    {
      ###################
      # PATTERN SEARCH
      ###################
      TOL.STAGE <- max(TOL.GOAL,TOL.LIKE)
      
      LINE.DO <- FALSE
      
      # stop conditions
      break
    }
      
    LINE.DID <- FALSE
    ##################
    # LINE SEARCH LOOP
    ##################
    if(LINE.DO)
    {
      LINE.DID <- TRUE
      LINE.DO <- FALSE
      LINE.SORTED <- TRUE
      
      # start of line search
      par.start <- par
      fn.start <- fn.par
      
      # initialize (all) storage results
      par.all <- cbind(par)
      fn.all <- fn.par
      
      par.diff <- par.target - par
      if(sum(par.diff^2) <= .Machine$double.eps) { break } # this actually happened once
      
      if(LINE.TYPE=="Prospection")
      {
        # generate a linear sequence of points from par to the other side of par.target
        P <- line.boxer(2*par.diff,p0=par,lower=lower,upper=upper,period=period)
        par.diff <- P - par # twice the old par.diff with no boundary (reflection=1)
        M <- sqrt(sum(par.diff^2)) # total search magnitude
        SEQ <- seq(0,M,length.out=mc.min(4,mc.cores)+1)[-1]
      }
      else if(LINE.TYPE=="Enclosure")
      {
        # generate a linear sequence of points between par and par.target (already evaluated)
        M <- sqrt(sum(par.diff^2))
        SEQ <- seq(0,M,length.out=mc.min(3,mc.cores)+2)[-1]
        SEQ <- SEQ[-length(SEQ)]
        
        par.all <- cbind(par.all,par.target)
        fn.all <- cbind(fn.all,fn.target)
        LINE.SORTED <- FALSE
      }
      
      # new points to evaluate
      P <- (DIR.STEP %o% SEQ) + par
      
      # start iteration loop
      while(counts < maxit)
      {
        counts <- counts + ceiling(ncol(P)/mc.cores)
        
        # most expensive part
        # evaluate objective function at new P and store to fn.queue
        fn.queue <- unlist(mclapply(split(P,col(P)),func,mc.cores=mc.cores))
        
        # combine with older results
        par.all <- cbind(par.all,P)
        fn.all <- c(fn.all,fn.queue)
        
        # sort along DIR.STEP
        if(!LINE.SORTED)
        {
          SORT <- c(DIR.STEP %*% (par.all-par))
          SORT <- sort(SORT,method="quick",index.return=TRUE)$ix
          par.all <- par.all[,SORT,drop=FALSE] # ARGH!
          fn.all <- fn.all[SORT]
          
          LINE.SORTED <- TRUE
        }
        
        # new best estimate
        MIN <- which.min(fn.all)
        par <- par.all[,MIN]
        fn.par <- fn.all[MIN]
        
        if(trace) { message(sprintf("%f %s search",zero+fn.par,LINE.TYPE)) }
        
        if(DEBUG)
        {
          graphics::plot(DIR.STEP %*% (par.all-par),fn.all-fn.par,xlab="Distance from MIN",ylab="Value over MIN")
          graphics::title(sprintf("%s search",LINE.TYPE))
        }
        
        # numerical degeneracy check
        TEST <- FALSE
        TEST <- TEST || (MIN<length(fn.all) && ZERO && fn.all[MIN+1]-fn.all[MIN]<=TOL.STAGE + TOL.ZERO)
        TEST <- TEST || (MIN<length(fn.all) && !ZERO && (fn.all[MIN+1]-fn.all[MIN])/fn.all[MIN]<=TOL.STAGE)
        TEST <- TEST || (MIN>1 && ZERO && fn.all[MIN-1]-fn.all[MIN]<=TOL.STAGE + TOL.ZERO)
        TEST <- TEST || (MIN>1 && !ZERO && (fn.all[MIN-1]-fn.all[MIN])/fn.all[MIN]<=TOL.STAGE)
        if(TEST) { break }
        
        #####################################
        # we didn't go far enough or we hit a boundary
        if(MIN==length(fn.all))
        {
          # if we hit a boundary, then we can stop at the boundary
          if(sum(is.boxed(par))>sum(BOX))
          { 
            par.target <- par
            break 
          }
          
          LINE.TYPE <- "Expansion"
          LINE.DO <- TRUE
          
          # Distance to first boundary that we will hit going in direction DIR.STEP>0
          M.BOX <- min(((upper-par)/DIR.STEP)[DIR.STEP>0],((lower-par)/DIR.STEP)[DIR.STEP<0])
          
          # we didn't hit a boundary yet, but we will eventually hit a boundary, so let's just do that
          if(M.BOX<Inf)
          {
            # generate a linear sequence of points that terminate at the eventual boundary
            SEQ <- seq(0,M.BOX,length.out=mc.cores+1)[-1]
            P <- (DIR.STEP %o% SEQ) + par
            
            # goto evaluate iteration step
            next
          }
          else # we didn't hit a boundary and we never will, because there is no boundary
          {
            # make sure we are accelerating
            M <- max(M,last(diff(SEQ)))
            
            # geometric sequence
            SEQ <- 2^(1:mc.cores) * M
            
            P <- (DIR.STEP %o% SEQ) + par
            
            # goto evaluate iteration step
            next
          }
        } # end expansion search setup
        
        ################################
        # our grid isn't refined enough to check for local quadratic behavior
        TEST <- (MIN<=2) || (MIN>=length(fn.all)-1)
        # is the function locally quadratic?
        if(!TEST) { TEST <- !QuadTest(c(DIR.STEP %*% par.all[,MIN+(-2):2]),fn.all[MIN+(-2):2]) }
        # in either case, we need refinement
        if(TEST)
        {
          LINE.TYPE <- "Refinement"
          LINE.DO <- TRUE
          
          # boundary limits
          M1 <- min(((par-lower)/DIR.STEP)[DIR.STEP>0],((upper-par)/-DIR.STEP)[DIR.STEP<0],Inf)
          M2 <- min(((par-lower)/-DIR.STEP)[DIR.STEP<0],((upper-par)/DIR.STEP)[DIR.STEP>0],Inf)
          
          # adjacence limits
          if(MIN>1) { M1 <- min(M1,DIR.STEP%*%(par-par.all[,MIN-1])) }
          if(MIN<length(fn.all)) { M2 <- min(M2,DIR.STEP%*%(par.all[,MIN+1]-par)) }
          
          # reflect if unknown
          if(M1==Inf) { M1 <- M2 }
          if(M2==Inf) { M2 <- M1 }
          
          # refining aims to make the grid even when filling in gaps adjacent to MIN
          # how many points to refine on each side of MIN
          n <- mc.min(2,mc.cores)
          # left solution if points were continuous
          n1 <- ((n+1)*M1-M2)/(M1+M2)
          # nearest discrete solutions
          n1 <- c(floor(n1),ceiling(n1))
          n2 <- n-n1
          # difference in gap sizes
          dev <- abs(M1/(n1+1) - M2/(n2+1))
          # discrete solution
          n1 <- n1[which.min(dev)]
          n2 <- n-n1
          
          SEQ <- seq(0,M1,length.out=n1+2)[-c(1,n1+2)]
          SEQ <- c(-rev(SEQ),seq(0,M2,length.out=n2+2)[-c(1,n2+2)])
          
          # combine for evaluation
          P <- (DIR.STEP %o% SEQ) + par
          
          LINE.SORTED <- FALSE
          next
        }
        
        # seems like this shouldn't happen?
        if(all(par==par.start))
        {
          LINE.TYPE <- "Shrink"
          LINE.DO <- TRUE
          
          # boundary limits
          M1 <- min(((par-lower)/DIR.STEP)[DIR.STEP>0],((upper-par)/-DIR.STEP)[DIR.STEP<0],Inf)
          M2 <- min(((par-lower)/-DIR.STEP)[DIR.STEP<0],((upper-par)/DIR.STEP)[DIR.STEP>0],Inf)
          
          # adjacence limits
          if(MIN>1) { M1 <- min(M1,(DIR.STEP%*%(par-par.all[,MIN-1]))/2) }
          if(MIN<length(fn.all)) { M2 <- min(M2,(DIR.STEP%*%(par.all[,MIN+1]-par))/2) }
          
          # reflect if unknown
          if(M1==Inf) { M1 <- M2 }
          if(M2==Inf) { M2 <- M1 }
          
          # generate geometrically tightening sequence around par
          SEQ <- (1/2)^(1:ceiling(mc.cores/2)-1)
          SEQ <- c(-M1*SEQ,M2*rev(SEQ))
          
          # combine for evaluation
          P <- (DIR.STEP %o% SEQ) + par
          LINE.SORTED <- FALSE
          # goto evaluate iteration step
          next
        }
        
        # everything must be good, line search finished
        break
      }
      # end iteration loop
      ######################################
      # we improved our estimate and can stop the line search
      
      #####################
      # LINE SEARCH FINISHER
      ######################
      # interpolate to an even better point

      i <- ifelse(MIN==1,MIN+2,MIN-1)
      P1 <- par.all[,i]
      F1 <- fn.all[i]
      
      i <- ifelse(MIN==length(fn.all),MIN-2,MIN+1)
      P2 <- par.all[,i]
      F2 <- fn.all[i]
      
      # calculate gradient and curvature
      DIFF <- QuadSolve(par,P1,P2,DIR.STEP,fn.par,F1,F2)
      gradient.LINE <- DIFF$gradient
      hessian.LINE <- abs(DIFF$hessian)
      # will update search direction hessian with this during next iteration
      
      # estimate better location if we didn't hit a boundary
      if(1<MIN && MIN<length(fn.all)) # MIN==1 shouldn't be possible here?
      { 
        # what was our target (closest if boundary)
        MIN.start <- which.min( colSums((par.target-par.all)^2) )
        # should we do a line search next time? Would Newton-Raphson have outright failed?
        if(fn.start<fn.all[MIN.start]) { LINE.DO <- TRUE }
        
        par.diff <- -(gradient.LINE/hessian.LINE)
        # who peformed better: Newton or line search?
        LINE.BETTER <- (MIN != MIN.start)
        if(TRUE) { par.target <- par + par.diff*DIR.STEP } # non-monotonic
        # what we predict this will be valued at
        fn.predict <- fn.par - hessian.LINE*par.diff^2/2
        
        TEST <- (ZERO && abs(fn.target.old-fn.predict) < TOL.STAGE + TOL.ZERO)
        TEST <- TEST || (!ZERO && abs((fn.target.old-fn.predict)/fn.predict) < TOL.STAGE)
        if(TEST) { INIT <- stage.progress() }
      }
      else
      {
        hessian.LINE <- 0 # won't update
        par.target <- par
        LINE.DO <- TRUE
        
        TEST <- (ZERO && abs(fn.target.old-fn.par) < TOL.STAGE + TOL.ZERO)
        TEST <- TEST || (!ZERO && abs((fn.target.old-fn.par)/fn.par) < TOL.STAGE)
        if(TEST) { INIT <- stage.progress() }
      }
    } # go back to differentiation
  } # end main loop
  

  if(counts<maxit) { convergence <- 0} else { convergence <- 1 }
  if(trace) { message(sprintf("%s in %d parallel function evaluations.",ifelse(convergence,"No convergence","Convergence"),counts)) }
  
  # transform back from tangent
  if(any(PERIOD))
  {
    # recenter on best par
    reset.tangent()
    # back to original angle scale
    par[PERIOD] <- angle * perscale
    covariance[PERIOD,PERIOD] <- t(t(covariance[PERIOD,PERIOD]*perscale)*perscale)
    hessian[PERIOD,PERIOD] <- t(t(hessian[PERIOD,PERIOD]/perscale)/perscale)
  }
  
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
