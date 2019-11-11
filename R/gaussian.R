gauss.comp <- function(fn,CTMM,COV=TRUE,...)
{
  CTMM <- listify(CTMM)
  axes <- CTMM[[1]]$axes
  AXES <- length(axes)

  # initialize
  par <- NULL
  parscale <- NULL
  lower <- NULL
  upper <- NULL

  # get parameter vector information
  for(i in 1:length(CTMM))
  {
    par <- c(par, CTMM[[i]]$mu[1,] )
    parscale <- c(parscale, sqrt(diag(CTMM[[i]]$sigma)) )
    lower <- c(lower, c(-Inf,-Inf) )
    upper <- c(upper, c(Inf,Inf) )

    if(CTMM[[i]]$isotropic)
    {
      par <- c(par,CTMM[[i]]$sigma@par[1])
      parscale <- c(parscale, CTMM[[i]]$sigma@par[1] )
      lower <- c(lower, 0 )
      upper <- c(upper, Inf )
    }
    else
    {
      par <- c(par,CTMM[[i]]$sigma@par)
      parscale <- c(parscale, c(CTMM[[i]]$sigma@par[1],1,pi/2) )
      lower <- c(lower, c(0,0,-Inf) )
      upper <- c(upper, c(Inf,Inf,Inf) )
    }
  }

  # evaluate function at arbitrary parameter vector
  FN <- function(PAR)
  {
    j <- 0
    for(i in 1:length(CTMM))
    {
      CTMM[[i]]$mu[1,] <- PAR[j + 1:AXES]
      j <- j + AXES
      isotropic <- CTMM[[i]]$isotropic
      if(isotropic)
      {
        CTMM[[i]]$sigma <- covm(PAR[j + 1],isotropic=isotropic,axes=axes)
        j <- j + 1
      }
      else
      {
        CTMM[[i]]$sigma <- covm(PAR[j + 1:3],isotropic=isotropic,axes=axes)
        j <- j + 3
      }
    }
    return(fn(CTMM))
  }

  MLE <- fn(CTMM)
  if(COV)
  {
    COV <- diag(0,length(par))

    grad <- genD(par,FN,lower=lower,upper=upper,parscale=parscale,mc.cores=1,order=1,drop=FALSE,...)$gradient
    # rows=FN, cols=par

    # fill COV matrix
    j <- 0
    for(i in 1:length(CTMM))
    {
      I <- j + 1:AXES
      if(length(dim(CTMM[[i]]$COV.mu))>AXES)
      { COV[I,I] <- CTMM[[i]]$COV.mu[,0,0,] }
      else
      { COV[I,I] <- CTMM[[i]]$COV.mu }
      j <- j + AXES
      isotropic <- CTMM[[i]]$isotropic
      if(isotropic)
      {
        I <- j + 1
        PAR <- names(CTMM[[i]]$sigma@par)[1]
        COV[I,I] <- CTMM[[i]]$COV[PAR,PAR]
        j <- j + 1
      }
      else
      {
        I <- j + 1:3
        PAR <- names(CTMM[[i]]$sigma@par)[1:3]
        COV[I,I] <- CTMM[[i]]$COV[PAR,PAR]
        j <- j + 3
      }
    }

    COV <- grad %*% COV %*% t(grad)
  }
  else
  { COV <- diag(0,length(MLE)) }

  return(list(MLE=MLE,COV=COV))
}
