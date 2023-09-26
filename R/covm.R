new.covm <- methods::setClass("covm", representation("matrix",par="numeric",isotropic="logical"))

# 2D covariance matrix universal format
# major - variance along major axis
# minor - variance along minor axis
# angle - angle of major axis
covm <- function(pars,isotropic=FALSE,axes=c("x","y"))
{
  if(is.null(pars)) { return(NULL) }

  if(class(pars)[1]=="covm") { pars <- pars@par }

  if(length(axes)==1)
  {
    sigma <- array(pars,c(1,1))
    pars <- sigma[1,1]

    names(pars) <- c("major")
  }
  else if(length(axes)==2)
  {
    if(length(pars)==1)
    {
      pars <- c(pars,pars,0)
      sigma <- diag(pars[1],2)
    }
    else if(length(pars)==3)
    {
      major <- max(pars[1:2]) # fix if out of order
      minor <- min(pars[1:2])
      theta <- ifelse(pars[1]==major,pars[3],pars[3]+pi/2) # rotate if out of order
      pars <- c(major,minor,theta)
      sigma <- sigma.construct(pars)
    }
    else if(length(pars)==4)
    {
      sigma <- pars
      if(class(pars)[1]=="covm") { pars <- attr(pars,'par') }
      else { pars <- sigma.destruct(sigma,isotropic=isotropic) }
    }

    # isotropic error check
    if(isotropic)
    {
      pars <- mean(pars[1:2])
      pars <- c(pars,pars,0)
      sigma <- diag(pars[1],2)
    }

    names(pars) <- c("major","minor","angle")
  }

  dimnames(sigma) <- list(axes,axes)

  new.covm(sigma,par=pars,isotropic=isotropic)
}


# construct covariance matrix from 1-3 parameters
sigma.construct <- function(pars)
{
  if(length(pars)==1) # isotropic
  {
    major <- pars[1]
    minor <- pars[1]
    theta <- 0

    sigma <- diag(major,2)
  }
  else
  {
    major <- pars[1]
    minor <- pars[2]
    theta <- pars[3]

    u <- c(cos(theta),sin(theta))
    v <- c(-sin(theta),cos(theta))

    sigma <- major*(u%o%u) + minor*(v%o%v)
  }

  return(sigma)
}


# reduce covariance matrix to 1-3 parameters
sigma.destruct <- function(sigma,isotropic=FALSE) # last arg not implemented
{
  INF <- diag(sigma)==Inf
  # catch infinite variances
  if(any(INF))
  {
    # can only represent diagonal matrix in that form
    major <- max(diag(sigma))
    minor <- min(diag(sigma))
    theta <- 0
  }
  else # finite variances
  {
    stuff <- eigen(sigma)
    stuff$values <- clamp(stuff$values,0,Inf)

    major <- stuff$values[1]
    minor <- stuff$values[2]

    if(major==minor)
    { theta <- 0 }
    else
    {
      theta <- stuff$vectors[,1]
      theta <- atan(theta[2]/theta[1])
    }
  }

  stuff <- c(major,minor,theta)
  names(stuff) <- c("major","minor","angle")
  return(stuff)
}


# eigen values of covariance matrix
eigenvalues.covm <- function(sigma)
{
  if(ncol(sigma)==1) { return(sigma@par['major']) }

  sigma <- attr(sigma,'par')[c('major','minor')]

  sigma <- sort(sigma,decreasing=TRUE)

  return(sigma)
}


# total variance or average variance
var.covm <- function(sigma,ave=FALSE)
{
  sigma <- eigenvalues.covm(sigma)

  if(ave) { sigma <- mean(sigma,na.rm=TRUE) }
  else { sigma <- sum(sigma,na.rm=TRUE) }

  return(sigma)
}


# determinant
det.covm <- function(sigma,ave=FALSE)
{
  sigma <- eigenvalues.covm(sigma)
  DIM <- length(sigma)

  sigma <- prod(sigma)

  if(ave) { sigma <- sigma^(1/DIM) }

  return(sigma)
}


# get geometric area/volume
area.covm <- function(sigma) { return( det.covm(sigma,ave=TRUE) ) }


# rescale covm
# default scales to variance-1
scale.covm <- function(sigma,value=1/var.covm(sigma,ave=TRUE))
{
  sigma <- sigma * value

  PARS <- 'major'
  if(length(sigma)==4) { PARS[2] <- 'minor' } # 2D
  attr(sigma,'par')[PARS] <- attr(sigma,'par')[PARS] * value

  return(sigma)
}


# squeeze factor and squeezability
squeezable.covm <- function(CTMM)
{
  AXES <- length(CTMM$axes)
  vars <- eigenvalues.covm(CTMM$sigma)

  # ratio of major axis to (geometric) mean axis (meters/meters)
  fact <- (max(vars)/min(vars))^(1/4)
  # can squeeze data to match variances? or extreme eccentricity
  able <- !is.nan(fact) && 4*abs(log(fact))<log(1/.Machine$double.eps)

  return(list(fact=fact,able=able))
}


# rotate covm by theta
rotate.covm <- function(sigma,theta=-sigma@par['angle'])
{
  if(length(sigma)==1) { return(sigma) }

  isotropic <- sigma@isotropic
  axes <- colnames(sigma)

  sigma <- attr(sigma,"par")
  sigma["angle"] <- sigma['angle'] + theta
  sigma <- covm(sigma,isotropic=isotropic,axes=axes)

  return(sigma)
}


# squeeze covm by factor smgm
squeeze.covm <- function(sigma,smgm=NULL,circle=FALSE)
{
  if(length(sigma)==1) { return(sigma) }

  isotropic <- sigma@isotropic
  axes <- colnames(sigma)

  if(is.null(smgm)) { circle <- TRUE }

  sigma <- attr(sigma,"par")
  if(circle)
  { sigma['major'] <- sigma['minor'] <- sqrt(sigma['major']*sigma['minor']) }
  else
  { sigma[c("major","minor")] <- c(1/smgm,smgm)^2 * sigma[c("major","minor")] }
  sigma <- covm(sigma,isotropic=isotropic,axes=axes)

  return(sigma)
}


# invert covariance matrix
solve_covm <- function(sigma,pseudo=FALSE)
{
  isotropic <- sigma@isotropic
  axes <- colnames(sigma)

  sigma <- attr(sigma,"par")
  PARS <- 1:min(2,length(sigma)) # major, (minor)
  sigma[PARS] <- rev(1/sigma[PARS]) # order matters
  if(length(sigma)>1) { sigma['angle'] <- sigma['angle'] + pi/2 } # reverse ordered

  sigma <- covm(sigma,isotropic=isotropic,axes=axes)

  return(sigma)
}


# square root of covariance matrix
sqrtm.covm <- function(sigma)
{
  isotropic <- sigma@isotropic
  axes <- colnames(sigma)

  sigma <- attr(sigma,"par")
  PARS <- 1:min(2,length(sigma)) # major, (minor)
  sigma[PARS] <- sqrt(sigma[PARS])

  sigma <- covm(sigma,isotropic=isotropic,axes=axes)

  return(sigma)
}

# matrix power
fn.covm <- function(sigma,fn)
{
  isotropic <- sigma@isotropic
  axes <- colnames(sigma)

  sigma <- attr(sigma,"par")
  PARS <- 1:min(2,length(sigma)) # major, (minor)
  sigma[PARS] <- fn(sigma[PARS])

  sigma <- covm(sigma,isotropic=isotropic,axes=axes)

  return(sigma)
}

# matrix power
mpow.covm <- function(sigma,pow) { fn.covm(sigma,function(s){s^pow}) }

# matrix logarithm
log_covm <- function(sigma,pow) { fn.covm(sigma,log) }

# matrix exponental
exp_covm <- function(sigma,pow) { fn.covm(sigma,exp) }


####### calculate variance and variance-covariance from major/minor information
# assumes that variance/covariance parameters come first in COV
axes2var <- function(CTMM,MEAN=TRUE)
{
  PAR <- c('major','minor','angle')
  COV <- CTMM$COV

  if(any(c('minor','angle') %nin% rownames(COV)))
  {
    NAMES <- rownames(COV)
    NAMES[NAMES=='major'] <- 'variance'
    dimnames(COV) <- list(NAMES,NAMES)
    return(COV)
  }

  OLD <- rownames(COV)
  OTHER <- OLD[OLD %nin% PAR]

  if(CTMM$isotropic[1])
  {
    NEW <- c("variance",OTHER)

    if(!MEAN)
    {
      COV['major',] <- 2 * COV['major',]
      COV[,'major'] <- 2 * COV[,'major']
    }
  }
  else
  {
    NEW <- c("variance",OTHER)
    grad <- matrix(0,length(NEW),length(OLD))
    rownames(grad) <- NEW
    colnames(grad) <- OLD
    if(length(OTHER)==1)
    { grad[OTHER,OTHER] <- 1 } # annoying that below isn't general
    else if(length(OTHER)>1)
    { diag(grad[OTHER,OTHER]) <- 1 }

    # convert major,minor/major uncertainty into mean-variance uncertainty
    grad['variance','major'] <- grad['variance','minor'] <- 1
    if(MEAN) { grad <- grad/2 } # average x-y variance

    COV <- grad %*% COV %*% t(grad)

    # backup for infinite covariances
    for(i in 1:nrow(COV))
    {
      if(any( is.nan(COV[i,]) | is.nan(COV[,i]) ))
      {
        COV[i,] <- COV[,i] <- 0
        COV[i,i] <- Inf
      }
    }
  }
  dimnames(COV) <- list(NEW,NEW)

  return(COV)
}


# gradient matrix d sigma / d par from par
J.sigma.par <- function(par)
{
  if(length(dim(par))==2) { par <- par@par }

  major <- par["major"]
  minor <- par["minor"]
  theta <- par["angle"]

  # s_xx = major*cos(theta)^2 + minor*sin(theta)^2
  # s_yy = major*sin(theta)^2 + minor*cos(theta)^2
  # s_xy = major*cos(theta)*sin(theta) - minor*sin(theta)*cos(theta)
  #      = (major-minor)*sin(2*theta)/2

  # gradient matrix d sigma / d par
  grad <- diag(3)
  dimnames(grad) <- list(c('xx','yy','xy'),names(par))
  # d sigma / d major
  grad[,'major'] <- c( cos(theta)^2, sin(theta)^2, +sin(2*theta)/2 )
  # d sigma / d minor
  grad[,'minor'] <- c( sin(theta)^2, cos(theta)^2, -sin(2*theta)/2 )
  # d sigma / d theta
  grad['xx','angle'] <- (minor-major)*sin(2*theta)
  grad['yy','angle'] <- (major-minor)*sin(2*theta)
  grad['xy','angle'] <- (major-minor)*cos(2*theta)

  return(grad)
}

# get linear COV[sigma] from ctmm object
sigma.COV <- function(CTMM)
{
  if(CTMM$isotropic[1])
  {
    VAR <- CTMM$COV['major','major']
    COV <- diag(c(VAR,VAR,0))
    rownames(COV) <- colnames(COV) <- c('xx','yy','xy')
    COV['xx','yy'] <- COV['yy','xx'] <- VAR
  }
  else
  {
    P <- c('major','minor','angle')
    COV <- CTMM$COV[P,P]
    J <- J.sigma.par(CTMM$sigma@par)
    COV <- J %*% COV %*% t(J)
  }

  return(COV)
}

# gradient matrix d par / d sigma from sigma
J.par.sigma <- function(sigma)
{
  if(length(dim(sigma))==2) { sigma <- sigma[c(1,4,2)] }
  names(sigma) <- c("xx","yy","xy")

  grad <- diag(3)
  dimnames(grad) <- list(c('major','minor','angle'),c('xx','yy','xy'))

  # these formulas come from Mathematica
  D <- (sigma['xx']-sigma['yy'])^2 + 4*sigma['xy']^2
  sqrt.D <- sqrt(D)
  R1 <- nant( (sigma['xx']-sigma['yy']) / sqrt.D , 0 )
  R2 <- nant( 4*sigma['xy'] / sqrt.D , 0 )

  grad['major',] <- c(1,1,0)/2 + c(+R1,-R1,+R2)/2
  grad['minor',] <- c(1,1,0)/2 + c(-R1,+R1,-R2)/2
  grad['angle',] <- nant( c(-sigma['xy'],sigma['xy'],sigma['xx']-sigma['yy'])/D , 0)

  return(grad)

  # s_xx = major*cos(theta)^2 + minor*sin(theta)^2
  # s_yy = major*sin(theta)^2 + minor*cos(theta)^2
  # s_xy = major*cos(theta)*sin(theta) - minor*sin(theta)*cos(theta)
  #      = (major-minor)*sin(2*theta)/2

  # s_xx + s_yy        = major + minor
  # s_xx*s_yy - s_xy^2 =

  par.fn <- function(s)
  {
    par <- matrix(s[c("xx","xy","xy","yy")],2,2)
    covm(par)@par
  }

  parscale <- sigma
  DIAG <- c('xx','yy')
  parscale['xy'] <- max( sqrt(prod(abs(sigma[DIAG]))) , mean(abs(sigma[DIAG])), abs(sigma['xy']) )
  lower <- c(0,0,-Inf)

  J <- genD(sigma,par.fn,lower=lower,parscale=parscale,order=1)
  J <- J$grad
  dimnames(J) <- list(c('major','minor','angle'),names(sigma))
  return(J)

  # solver method - could fail
  J <- J.sigma.par(sigma)
  J <- solve(J)
  return(J)

  # analytic calculation - unfinished, requires a limit corner case solution
  TR <- sum(sigma[c("xx","yy")])
  DET <- sigma["xx"]*sigma["yy"] - sigma["xy"]^2
}

# return the COV matrix for covm par representation
COV.covm <- function(sigma,n,k=1,REML=TRUE)
{
  isotropic <- attr(sigma,'isotropic')
  par <- attr(sigma,'par')
  sigma <- methods::getDataPart(sigma)
  DIM <- sqrt(length(sigma))

  DOF.mu <- n
  COV.mu <- sigma/n

  if(REML) { n <- n-k }

  if(isotropic || DIM==1)
  {
    COV <- cbind( 2*par["major"]^2/(n*DIM) )
    dimnames(COV) <- list("major","major")
  }
  else # 2D
  {
    # covariance matrix for c( sigma_xx , sigma_yy , sigma_xy )
    COV <- diag(0,3)
    NAMES <- c('xx','yy','xy')
    dimnames(COV) <- list(NAMES,NAMES)
    COV[c('xx','yy'),c('xx','yy')] <- 2/n * sigma^2
    COV['xy',] <- c( 2*sigma[1,1]*sigma[1,2] , 2*sigma[2,2]*sigma[1,2] , sigma[1,1]*sigma[2,2]+sigma[1,2]^2 )/n
    COV[,'xy'] <- COV['xy',]

    # gradient matrix d sigma / d par
    grad <- J.sigma.par(par)
    # gradient matrix d par / d sigma via inverse function theorem
    grad <- PDsolve(grad,sym=FALSE)

    COV <- (grad) %*% COV %*% t(grad)
    COV <- nant(COV,0) # 0/0 for inactive

    COV <- He(COV) # asymmetric errors
    dimnames(COV) <- list(names(par),names(par))
  }

  return(list(COV=COV,COV.mu=COV.mu,DOF.mu=DOF.mu))
}


# return the canonical parameters of a covariance matrix
pars_covm <- function(COVM)
{
  if(COVM@isotropic)
  { return(COVM@par[1]) }
  else
  { return(COVM@par) }
}
