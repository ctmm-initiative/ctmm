################################
# Return hidden state estimates or simulations
################################
smoother <- function(DATA,CTMM,...)
{
  CTMM <- ctmm.prepare(DATA,CTMM)
  
  t <- DATA$t
  
  dt <- c(Inf,diff(t))
  n <- length(t)
  
  u <- array(1,n)
  
  isotropic <- CTMM$isotropic
  sigma <- CTMM$sigma
  area <- sigma@par[1]
  ecc <- sigma@par[2]
  theta <- sigma@par[3]

  K <- length(CTMM$tau)
  
  circle <- CTMM$circle
  if(circle) { circle <- 2*pi/circle }
  if(abs(circle) == Inf) { circle <- FALSE }
  
  R <- cbind(DATA$x,DATA$y)
  # stationary mean function
  u <- array(1,n)
  
  V <- array(0,dim=c(n,2))
  COV <- array(0,dim=c(n,2,2))
  
  # do we need to orient the data along the major an minor axes of sigma
  ROTATE <- !isotropic && (CTMM$error || circle)
  if(ROTATE)
  { 
    M <- rotate(-theta)
    R <- R %*% t(M)
    mu <- M %*% as.numeric(CTMM$mu)
  }
  
  if(circle)
  {
    # proportional standardization from ellipse to circle
    if(ecc)
    {
      R[,1] <- R[,1] * exp(-ecc/4)
      R[,2] <- R[,2] * exp(+ecc/4)
    }
    R <- cbind(R[,1] + 1i*R[,2])
    
    # corotating frame
    M <- exp(-1i*circle*(t-t[1]))
    R <- M * R
    u <- M * u
    
    CTMM$sigma <- area
    KALMAN <- kalman(R,u,dt=dt,CTMM=CTMM,error=DATA$error,...)
    
    R <- KALMAN$Z[,"position",1]
    if(K>1) { V <- KALMAN$Z[,"velocity",1] }
    COV[,1,1] <- KALMAN$S[,"position","position"]
    COV[,2,2] <- KALMAN$S[,"position","position"]
    
    # spin back + velocity of frame
    M <- Conj(M)
    R <- M * R
    mu <- mu[1] + 1i*mu[2]
    V <- (M * V) + 1i*circle*(R-mu)
    
    R <- cbind(Re(R),Im(R))
    V <- cbind(Re(V),Im(V))
    # unstandardize
    if(ecc)
    {
      R[,1] <- R[,1] * exp(+ecc/4)
      R[,2] <- R[,2] * exp(-ecc/4)

      V[,1] <- V[,1] * exp(+ecc/4)
      V[,2] <- V[,2] * exp(-ecc/4)
    
      COV[,1,1] <- COV[,1,1] * exp(+ecc/2)
      COV[,2,2] <- COV[,2,2] * exp(-ecc/2)
    }
  }
  else if(CTMM$isotropic || !CTMM$error)
  {
    CTMM$sigma <- area
    KALMAN <- kalman(R,u,dt=dt,CTMM=CTMM,error=DATA$error,...)
    
    R <- KALMAN$Z[,"position",]
    if(K>1) { V <- KALMAN$Z[,"velocity",] }
    COV <- KALMAN$S[,"position","position"] %o% covm(c(1,ecc,theta),isotropic)
  }
  else
  {
    #diagonalize data and then run two 1D Kalman filters with separate means
    # major axis likelihood
    CTMM$sigma <- area * exp(+ecc/2)
    KALMAN <- kalman(cbind(R[,1]),u,dt=dt,CTMM=CTMM,error=DATA$error,...)
    
    R[,1] <- KALMAN$Z[,"position",1]
    if(K>1) { V[,1] <- KALMAN$Z[,"velocity",1] }
    COV[,1,1] <- KALMAN$S[,"position","position"]

    # minor axis likelihood
    CTMM$sigma <- area * exp(-ecc/2)
    KALMAN <- kalman(cbind(R[,2]),u,dt=dt,CTMM=CTMM,error=DATA$error,...)
    
    R[,2] <- KALMAN$Z[,"position",1]
    if(K>1) { V[,2] <- KALMAN$Z[,"velocity",1] }
    COV[,2,2] <- KALMAN$S[,"position","position"]
  }
  
  if(ROTATE)
  {
    # transform results back
    M <- rotate(+theta)
    
    # there MUST be an easier way to do a simple inner product but %*% fails
    COV <- aperm(COV,perm=c(2,1,3))
    dim(COV) <- c(2,2*n)
    COV <- M %*% COV
    dim(COV) <- c(2*n,2)
    COV <- COV %*% t(M)
    dim(COV) <- c(2,n,2)
    COV <- aperm(COV,perm=c(2,1,3))
    
    R <- t(M %*% t(R))
    V <- t(M %*% t(V))
  }
  
  RETURN <- list(t=t,R=R,COV=COV)
  if(K>1) { RETURN$V <- V }
  
  return(RETURN)
}


########################################
# fill in data gaps with missing observations
########################################
fill.data <- function(data,CTMM=ctmm(tau=Inf),verbose=FALSE,t=NULL,dt=NULL,res=1,cor.min=0,dt.max=NULL)
{
  # prepare data error
  data$error <- get.error(data,CTMM)
  
  # FIX THE TIME GRID TO AVOID TINY DT
  if(is.null(t))
  {
    t <- data$t
    DT <- diff(t)
    
    # target resolution
    if(is.null(dt)){ dt <- stats::median(DT)/res }
    
    # maximum gap to bridge
    if(is.null(dt.max)) { dt.max <- -log(cor.min)*CTMM$tau[1] }
    
    # this regularization is not perfectly regular, but holds up to sampling drift in caribou data
    t.grid <- c()  # full (even-ish) grid
    dt.grid <- c() # local (numeric) sampling resolution 
    t.new <- c()   # new times in this even grid
    for(i in 1:length(DT))
    {
      n.sub <- round(DT[i]/dt)+1
      n.sub <- max(n.sub,2) # fix for crazy small time-steps
      t.sub <- seq(from=t[i],to=t[i+1],length.out=n.sub)
      dt.sub <- DT[i]/(n.sub-1)
      
      # skip low correlation times in gaps
      INCLUDE <- (t.sub-t[i])<=dt.max/2 | (t[i+1]-t.sub)<=dt.max/2
      t.sub <- t.sub[INCLUDE]
      n.sub <- length(t.sub)

      dt.sub <- rep(dt.sub,n.sub)
      t.grid <- c(t.grid,t.sub)
      dt.grid <- c(dt.grid,dt.sub)
      
      t.new <- c(t.new,t.sub[c(-1,-n.sub)])
    }
    
    # half weight repeated endpoints in grid
    w.grid <- dt.grid
    REPEAT <- where(diff(t.grid)==0)
    w.grid[REPEAT] <- w.grid[REPEAT]/2
    w.grid[REPEAT+1] <- w.grid[REPEAT+1]/2
  }
  else # use a pre-specified time grid
  {
    t.new <- t[!(t %in% data$t)]
  }

  # empty observation row for these times
  blank <- data[1,]
  blank$error <- Inf
  blank <- blank[rep(1,length(t.new)),]
  blank$t <- t.new
  
  # attach empty measurements to data
  data <- rbind(data,blank)
  
  # sort times
  data <- data[sort.list(data$t,na.last=NA,method="quick"),]
  # this is now our fake data set to feed into the kalman smoother
  
  if(verbose)
  { data <- list(data=data,t.grid=t.grid,dt.grid=dt.grid,w.grid=w.grid) }

  return(data)
}

                      
#################################
# Kriged Kernel Density Estimate
# H is your additional smoothing bandwidth matrix (zero by default)
# resolution is the number of kriged locations per median step
# cor.min is roughly the correlation required between locations to bridge them
# dt.max is (alternatively) the maximum gap allowed between locations to bridge them
#################################
occurrence <- function(data,CTMM,H=0,res.time=10,res.space=10,grid=NULL,cor.min=0.5,dt.max=NULL)
{
  dt <- stats::median(diff(data$t))
  
  if(length(H)==1) { H <- diag(H,2) }
  
  info <- attr(data,"info")
  CTMM <- ctmm.prepare(data,CTMM)
  error <- get.error(data,CTMM)
  MIN.ERR <- min(error)
  
  # format data to be relatively evenly spaced with missing observations
  data <- fill.data(data,CTMM,verbose=TRUE,res=res.time,cor.min=cor.min,dt.max=dt.max)
  t.grid <- data$t.grid
  dt.grid <- data$dt.grid
  w.grid <- data$w.grid
  data <- data$data
  # run through smoother to get   
  state <- smoother(data,CTMM,smooth=TRUE)

  # evenly sampled subset: data points (bridge ends) may be counted twice and weighted half
  GRID <- c( where(data$t %in% t.grid) , where(data$t %in% t.grid[where(diff(t.grid)==0)]) )
  GRID <- sort.int(GRID,method="quick")

  # GRID <- data$t %in% t.grid
  
  # t <- state$t[GRID]
  R <- state$R[GRID,]
  COV <- state$COV[GRID,,]
  n <- length(R[,1])
  
  # continuous velocities will give us more information to use
  if(length(CTMM$tau)>1)
  { V <- state$V[GRID,] }
  else # null velocity data otherwise
  { V <- array(0,c(n,2)) }

  # fake data
  data <- data.frame(x=R[,1],y=R[,2])

  # some covariances are slightly negative due to roundoff error
  # so here I clamp the negative singular values to zero
  COV <- lapply(1:n,function(i){ PDclamp(COV[i,,]) })
  COV <- unlist(COV)
  dim(COV) <- c(2,2,n)
  COV <- aperm(COV,perm=c(3,1,2))
  
  # uncertainties/bandwidths for this data
  h <- H # extra smoothing
  H <- array(0,c(n,2,2))
  for(i in 1:n)
  {
    # total covariance is bandwidth + krige uncertainty + kinetic numerical correction
    H[i,,] <- h + COV[i,,] + dt.grid[i]^2/12*(V[i,] %o% V[i,])
    # maybe publish the last part as a technical note
  }
  # there is probably a faster way to do that
  
  # estimate size of data blob
  dr <- diag(CTMM$sigma)
  if(CTMM$range)
  {
    if(length(CTMM$tau)==1) #OU
    { dr <- dt/4 * dr/CTMM$tau[1] }
    else if(length(CTMM$tau)==2) #OUF
    { dr <- dt^2/24 * dr/prod(CTMM$tau)}
  }
  else # these sigmas are var/tau[1] diffusion limits
  {
    if(length(CTMM$tau)==1) #BM
    { dr <- dt/4 * dr }
    else if(length(CTMM$tau)==2) #IOU
    { dr <- dt^2/24 * dr/CTMM$tau[2] }
  }
  if(CTMM$error){ dr <- dr + MIN.ERR }
  dr <- sqrt(dr)
  
  # using the same data format as AKDE, but with only the ML estimate (alpha=1)
  KDE <- kde(data,H=H,W=w.grid,dr=dr/res.space)
  KDE$H <- diag(0,2)
  KDE <- new.UD(KDE,info=info)
  return(KDE)
}


##############################################
# SIMULATE DATA over time array t
simulate.ctmm <- function(object,nsim=1,seed=NULL,data=NULL,t=NULL,dt=NULL,res=1,...)
{
  if(!is.null(seed)){ set.seed(seed) }
  
  info <- attr(object,"info")

  # Gaussian simulation not conditioned off of any data
  if(is.null(data))
  {
    object <- ctmm.prepare(list(t=t),object)
    
    tau <- object$tau
    if(is.null(tau)) { tau = 0 }
    K <- length(tau)
    
    mu <- object$mu
    if(is.null(mu)) { mu <- array(0,c(1,2)) }
    sigma <- object$sigma
    if(is.null(sigma)) { sigma <- diag(1,2) }
    
    # I cannot figure out how to make this "Note" go away!
    # sigma <- Matrix::Matrix(sigma,sparse=FALSE,doDiag=FALSE)
    suppressWarnings(Lambda <- expm::sqrtm(sigma))
    
    K <- length(tau)
    
    n <- length(t)
    dt <- c(Inf,diff(t)) # time lags
    
    # where we will store the data
    x <- rep(0,times=n)
    y <- rep(0,times=n)
    
    # initial hidden state, for standardized process
    Hx <- rep(0,times=K)
    Hy <- rep(0,times=K)
    
    object$sigma <- 1
    for(i in 1:n)
    {
      # tabulate propagators if necessary
      if((i==1)||(dt[i]!=dt[i-1]))
      {
        Langevin <- langevin(dt=dt[i],CTMM=object)
        Green <- Langevin$Green
        Sigma <- Langevin$Sigma
      }
      
      # standardized process
      Hx <- MASS::mvrnorm(n=1,mu=as.vector(Green %*% Hx),Sigma=Sigma)
      Hy <- MASS::mvrnorm(n=1,mu=as.vector(Green %*% Hy),Sigma=Sigma)
      
      x[i] <- Hx[1]
      y[i] <- Hy[1]
    }
    
    # rotate process if necessary
    circle <- object$circle
    if(circle)
    {
      circle <- 2*pi/circle
      r <- exp(1i*circle*(t-t[1]))
      r <- r * (x + 1i*y)
      x <- Re(r)
      y <- Im(r)
    }
    
    r <- cbind(x,y)
    r <- r %*% Lambda
    
    # calculate mean function
    mu <- object$mean.vec %*% mu
    r <- r + mu
    
    data <- data.frame(t=t,x=r[,1],y=r[,2])
  }
  else # condition off of the data
  {
    object <- ctmm.prepare(data,object)
    data <- fill.data(data,CTMM=object,t=t,dt=dt,res=res)
    object$error <- TRUE # avoids unit variance algorithm
    data <- smoother(data,object,sample=TRUE)
    data <- data.frame(t=data$t,x=data$R[,1],y=data$R[,2])
  }
  
  data <- new.telemetry(data,info=info)
  return(data)
}
#methods::setMethod("simulate",signature(object="ctmm"), function(object,...) simulate.ctmm(object,...))
