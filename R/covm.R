new.covm <- methods::setClass("covm", representation("matrix",par="numeric",isotropic="logical"))

# 2D covariance matrix universal format
# major - variance along major axis
# minor - ratio of minor/major variances (1/condition number)
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
      pars <- c(pars,1,0)
      sigma <- diag(pars[1],2)
    }
    else if(length(pars)==3)
    { sigma <- sigma.construct(pars) }
    else if(length(pars)==4)
    {
      sigma <- pars
      if(class(pars)[1]=="covm") { pars <- attr(pars,'par') }
      else { pars <- sigma.destruct(sigma,isotropic=isotropic) }
    }

    # isotropic error check
    if(isotropic)
    {
      pars <- c(pars[1]*(1+pars[2])/2,1,0)
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
    minor <- 1
    theta <- 0
  }
  else
  {
    major <- pars[1]
    minor <- pars[2]
    theta <- pars[3]
  }

  u <- c(cos(theta),sin(theta))
  v <- c(-sin(theta),cos(theta))

  sigma <- major*( (u%o%u) + minor*(v%o%v) )

  return(sigma)
}


# reduce covariance matrix to 1-3 parameters
sigma.destruct <- function(sigma,isotropic=FALSE) # last arg not implemented
{
  INF <- diag(sigma)==Inf
  # catch infinite variances
  if(any(INF))
  {
    # can only represent diagonal infinite matrix in that form
    major <- Inf
    minor <- 1
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

    if(major<=0) { minor <- 1 }
    else { minor <- minor/major }
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
  sigma['minor'] <- sigma['minor']*sigma['major']

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
  attr(sigma,'par')['major'] <- attr(sigma,'par')['major'] * value
  return(sigma)
}


# squeeze factor and squeezability
squeezable.covm <- function(CTMM)
{
  AXES <- length(CTMM$axes)
  minor <- CTMM$sigma@par['minor'] # NA in 1D

  # ratio of major axis to (geometric) mean axis (meters/meters)
  fact <- (1/minor)^(1/4)
  # extreme eccentricity --- cannot squeeze data to match variances
  able <- AXES==2 && !is.nan(fact) && 4*abs(log(fact))<log(1/.Machine$double.eps)

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
  {
    sigma['major'] <- sigma['major'] * sqrt(sigma['minor'])
    sigma['minor'] <- 1
  }
  else
  {
    sigma['minor'] <- sigma['minor'] * sigma['major']
    sigma[c("major","minor")] <- c(1/smgm,smgm)^2 * sigma[c("major","minor")]

    major <- max(sigma[1:2])
    if(major>0) { minor <- min(sigma[1:2])/major }

    sigma[c("major","minor")] <- c(major,minor)
  }
  sigma <- covm(sigma,isotropic=isotropic,axes=axes)

  return(sigma)
}


# invert covariance matrix
solve.covm <- function(sigma,pseudo=FALSE)
{
  isotropic <- sigma@isotropic
  axes <- colnames(sigma)

  sigma <- attr(sigma,"par")
  PARS <- 1:min(2,length(sigma)) # major, (minor)
  sigma[PARS] <- 1/sigma[PARS]

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


####### calculate variance and variance-covariance from major/minor information
axes2var <- function(CTMM,MEAN=TRUE)
{
  COV <- CTMM$COV
  NAMES <- rownames(COV)

  if(CTMM$isotropic)
  {
    NAMES <- c("variance",NAMES[-1])

    if(!MEAN)
    {
      COV['major',] <- 2 * COV['major',]
      COV[,'major'] <- 2 * COV[,'major']
    }
  }
  else
  {
    NAMES <- c("variance",NAMES[-(1:3)])

    sigma <- CTMM$sigma@par
    major <- sigma['major']
    minor <- sigma['minor']

    # convert major,minor/major uncertainty into mean-variance uncertainty
    grad <- c(1+minor,major,0)  # total x-y variance
    if(MEAN) { grad <- grad/2 } # average x-y variance

    P <- nrow(COV)
    if(P>3)
    {
      grad <- rbind( grad , array(0,c(P-3,3)) )
      grad <- cbind( grad , rbind( rep(0,P-3) , diag(1,P-3) ) )
    }
    else
    { grad <- rbind(grad) }

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
  dimnames(COV) <- list(NAMES,NAMES)

  return(COV)
}


# gradient matrix d sigma / d par
J.sigma.par <- function(par)
{
  major <- par["major"]
  minor <- par["minor"]
  theta <- par["angle"]

  # s_xx = major*cos(theta)^2 + major*minor*sin(theta)^2
  # s_yy = major*sin(theta)^2 + major*minor*cos(theta)^2
  # s_xy = major*cos(theta)*sin(theta) - major*minor*sin(theta)*cos(theta)
  #      = major*(1-minor)*sin(2*theta)/2

  # gradient matrix d sigma / d par
  grad <- diag(3)
  names(grad) <- names(par)
  # d sigma / d major
  grad[,1] <- c( cos(theta)^2+minor*sin(theta)^2, sin(theta)^2+minor*cos(theta)^2, +(1-minor)*sin(2*theta)/2 )
  # d sigma / d minor
  grad[,2] <- c( major*sin(theta)^2, major*cos(theta)^2, -major*sin(2*theta)/2 )
  # d sigma / d theta
  grad[1,3] <- major*(-1 + minor)*sin(2*theta)
  grad[2,3] <- major*(+1 - minor)*sin(2*theta)
  grad[3,3] <- major*(+1 - minor)*cos(2*theta)

  return(grad)
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
    COV[1:2,1:2] <- 2/n * sigma^2
    COV[3,] <- c( 2*sigma[1,1]*sigma[1,2] , 2*sigma[2,2]*sigma[1,2] , sigma[1,1]*sigma[2,2]+sigma[1,2]^2 )/n
    COV[,3] <- COV[3,]

    # gradient matrix d sigma / d par
    grad <- J.sigma.par(par)
    # gradient matrix d par / d sigma via inverse function theorem
    grad <- PDsolve(grad)

    COV <- (grad) %*% COV %*% t(grad)
    COV <- nant(COV,0) # 0/0 for inactive

    COV <- He(COV) # asymmetric errors
    dimnames(COV) <- list(names(par),names(par))
  }

  return(list(COV=COV,COV.mu=COV.mu,DOF.mu=DOF.mu))
}


# return the canonical parameters of a covariance matrix
pars.covm <- function(COVM)
{
  if(COVM@isotropic)
  { return(COVM@par[1]) }
  else
  { return(COVM@par) }
}
