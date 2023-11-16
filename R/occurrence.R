#################################
# Kriged Kernel Density Estimate
# H is your additional smoothing bandwidth matrix (zero by default)
# resolution is the number of kriged locations per median step
# cor.min is roughly the correlation required between locations to bridge them
# dt.max is (alternatively) the maximum gap allowed between locations to bridge them
#################################

# wrapper for multiple individuals
occurrence <- function(data,CTMM,R=list(),SP=NULL,SP.in=TRUE,H=0,variable="utilization",res.time=10,res.space=10,grid=NULL,cor.min=0.05,dt.max=NULL,buffer=TRUE,...)
{
  if(length(projection(data))>1) { stop("Data not in single coordinate system.") }
  validate.grid(data,grid)

  DROP <- class(data)[1] != "list"
  data <- listify(data)
  CTMM <- listify(CTMM)
  n <- length(data)

  # force grids to be compatible
  COMPATIBLE <- length(data)>1 && !is.grid.complete(grid)

  axes <- CTMM[[1]]$axes

  grid <- format_grid(grid,axes=axes)
  COL <- length(axes)

  KDE <- list()
  for(i in 1:n)
  { KDE[[i]] <- currence(data[[i]],CTMM[[i]],H=H,variable=variable,res.time=res.time,res.space=res.space,grid=grid,cor.min=cor.min,dt.max=dt.max,buffer=buffer,...) }
  W <- sapply(KDE,function(k){sum(k$W)})

  # determine desired (absolute) resolution
  dr <- sapply(1:n,function(i){KDE[[i]]$dr})
  dim(dr) <- c(COL,n)
  dr <- apply(dr,1,grid$dr.fn)

  if(COMPATIBLE) # force grids compatible
  {
    grid$align.to.origin <- TRUE
    if("dr" %nin% names(grid)) { grid$dr <- dr }
  }

  # finish distribution calculation
  for(i in 1:n)
  {
    # using the same data format as AKDE, but with only the ML estimate (alpha=1)
    KDE[[i]]  <- kde(KDE[[i]]$data,CTMM=CTMM[[i]],H=KDE[[i]]$H,W=KDE[[i]]$W,dr=dr,grid=grid,SP=SP,SP.in=SP.in,RASTER=R)
    KDE[[i]]$H <- diag(0,2)
    dimnames(KDE[[i]]$H) <- list(axes,axes)
    KDE[[i]]$W <- W[i]
    KDE[[i]] <- new.UD(KDE[[i]],info=attr(data[[i]],"info"),type='occurrence',CTMM=CTMM[[i]])
  }
  names(KDE) <- names(data)
  if(DROP) { KDE <- KDE[[1]] }

  return(KDE)
}

# occurrence for single indviduals
currence <- function(data,CTMM,H=0,variable="utilization",res.time=10,res.space=10,grid=NULL,cor.min=0.05,dt.max=NULL,buffer=TRUE,...)
{
  if(length(CTMM$tau)<2)
  {
    if(variable=='revisitation')
    { stop("Revisitation estimation requires a continuous-velocity model.") }
    else if(variable=='speed')
    { stop("Speed estimation requires a continuous-velocity model.") }
  }

  validate.grid(data,grid)

  axes <- CTMM$axes
  CTMM0 <- CTMM
  dt <- stats::median(diff(data$t))
  if(is.null(dt.max)) { dt.max <- dt }

  if(length(H)==1) { H <- diag(H,2) }

  info <- attr(data,"info")
  SIGMA <- CTMM$sigma # diffusion matrix for later
  CTMM <- ctmm.prepare(data,CTMM,precompute=FALSE) # not the final t for calculating u
  error <- get.error(data,CTMM,circle=TRUE)
  MIN.ERR <- min(error) # Fix something here?

  # format data to be relatively evenly spaced with missing observations
  data <- fill.data(data,CTMM,verbose=TRUE,res=res.time,cor.min=cor.min,dt.max=dt.max,buffer=buffer)
  t.grid <- data$t.grid
  dt.grid <- data$dt.grid
  w.grid <- data$w.grid
  data <- data$data

  # calculate trend
  drift <- drift.mean(CTMM0,data$t) %*% CTMM0$mu

  # detrend for smoothing - retrend later
  z <- get.telemetry(data,axes=axes)
  data[,axes] <- z - drift

  # smooth mean-zero data # run through smoother to get
  state <- smoother(data,CTMM,smooth=TRUE)

  # skip repeated timestamps in data (full information retained)
  KEEP <- !data$skip
  data <- data[KEEP,]
  state$R <- state$R[KEEP,,drop=FALSE]
  state$COV <- state$COV[KEEP,,,drop=FALSE]
  if(length(CTMM$tau)>1) { state$V <- state$V[KEEP,,drop=FALSE] }
  drift <- drift[KEEP,,drop=FALSE]

  # detrend for simulation - retrend later
  state$R <- state$R + drift

  # evenly sampled subset: data points (bridge ends) may be counted twice and weighted half
  GRID <- c( which(data$t %in% t.grid) , which(data$t %in% t.grid[which(diff(t.grid)==0)]) )
  GRID <- sort.int(GRID)
  # GRID <- data$t %in% t.grid

  # t <- state$t[GRID]
  R <- state$R[GRID,]
  COV <- state$COV[GRID,,]
  n <- length(R[,1])

  # continuous velocities will give us more information to use
  if(length(CTMM$tau)>1 && dt<CTMM$tau[2])
  { V <- state$V[GRID,] }
  else # null velocity data otherwise
  { V <- array(0,c(n,2)) }

  # fake data
  data <- data.frame(x=R[,1],y=R[,2])

  if(variable %in% c("revisitation","speed"))
  {
    data$vx <- state$V[GRID,'vx']
    data$vy <- state$V[GRID,'vy']
    data$COV.vx.vx <- state$VCOV[GRID,'vx','vx']
    data$COV.vx.vy <- state$VCOV[GRID,'vx','vy']
    data$COV.vy.vy <- state$VCOV[GRID,'vy','vy']
    w.grid <- w.grid * speeds_fast(data,append=TRUE)$speed
  }

  # some covariances are slightly negative due to roundoff error
  # so here I clamp the negative singular values to zero
  COV <- vapply(1:n,function(i){ PDclamp(COV[i,,]) },COV[1,,]) # (d,d,n)
  COV <- aperm(COV,perm=c(3,1,2)) # (n,d,d)

  # uncertainties/bandwidths for this data
  H <- vapply(1:n,function(i){ H + COV[i,,] + dt.grid[i]^2/12*(V[i,] %o% V[i,]) },H) # (2,2,n)
  H <- aperm(H,c(3,1,2)) # (n,2,2)

  # estimate size of data blob
  dr <- diag(SIGMA)
  if(CTMM$range && length(CTMM$tau)) { dr <- dr/CTMM$tau[1] }
  # prefactors from mid-bridge variance
  if(length(CTMM$tau)==1) #BM/OU
  { dr <- dt/4 * dr }
  else if(length(CTMM$tau)==2) #IOU/OUF
  { dr <- dt^2/24 * dr/CTMM$tau[2] }

  if(any(CTMM$error>0)){ dr <- dr + MIN.ERR }
  dr <- sqrt(dr)

  # return list of stuff to work with
  STUFF <- list(data=data,H=H,W=w.grid,dr=dr/res.space)
}
