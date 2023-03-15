sdm.fit <- function(data,R=list(),formula=NULL,area=NULL,reference="auto",standardize=TRUE,integrator="MonteCarlo",error=0.01,max.mem="1 Gb",interpolate=TRUE,trace=TRUE,...)
{
  UD <- sdm.UD(data)

  rsf.fit(data,UD=UD,R=R,fomula=formula,integrated=is.null(area),level.UD=area,reference=reference,smooth=FALSE,standardize=standardize,integrator=integrator,error=error,max.mem=max.mem,interpolate=interpolate,trace=trace,...)
}

sdm.select <- function(data,R=list(),formula=NULL,area=NULL,verbose=FALSE,IC="AICc",trace=TRUE,...)
{
  UD <- sdm.UD(data)

  rsf.select(data,UD=UD,R=R,formula=formula,integrated=is.null(area),level.UD=area,verbose=verbose,IC=IC,trace=trace,smooth=FALSE,...)
}

sdm.UD <- function(data,res=10)
{
  # GUESS <- ctmm(tau=NULL,isotropic=TRUE)
  # CTMM <- ctmm.select(data,GUESS)
  CTMM <- ctmm.fit(data,ctmm(isotropic=TRUE,tau=NULL))

  UD <- list()
  UD$weights <- rep(1,nrow(data))
  UD$DOF.area <- nrow(data) - 1
  UD$dr <- sqrt( c(x=CTMM$sigma['x','x'],y=CTMM$sigma['y','y']) ) / res
  UD <- new.UD(UD,CTMM=CTMM)

  return(UD)
}


sdm.integrate <- function(biased=NULL,bias=NULL,unbiased=NULL)
{
  PAR <- NULL
  copy.beta <- function(beta)
  {
    BETA <- numeric(length(PAR))
    names(BETA) <- PAR
    BETA[names(beta)] <- beta
    return(BETA)
  }

  copy.PRE <- function(pre)
  {
    PRE <- matrix(0,length(PAR),length(PAR))
    dimnames(PRE) <- list(PAR,PAR)
    PRE[rownames(pre),colnames(pre)] <- pre
    return(PRE)
  }

  if(!is.null(biased) && !is.null(bias))
  {
    biased.beta <- biased$beta
    PAR <- names(biased.beta)
    biased.COV <- biased$COV[PAR,PAR]
    # precision matrix
    biased.PRE <- PDsolve(biased.COV)

    bias.beta <- bias$beta
    PAR <- names(bias.beta)
    bias.COV <- bias$COV[PAR,PAR]
    bias.PRE <- PDsolve(bias.COV)

    PAR <- unique(c(names(biased.beta),names(bias.beta)))
    biased.beta <- copy.beta(biased.beta)
    bias.beta <- copy.beta(bias.beta)
    biased.PRE <- copy.PRE(biased.PRE)
    bias.PRE <- copy.PRE(bias.PRE)

    debiased.beta <- biased.beta - bias.beta
    debiased.PRE <- biased.PRE + bias.PRE

    # just in case this is everything
    beta <- debiased.beta
    COV <- PDsolve(debiased.PRE)
  }
  else
  {
    debiased.PRE <- NULL
    debiased.beta <- NULL
  }

  if(!is.null(unbiased))
  {
    if(class(unbiased)[1]=="list")
    {
      unbiased <- mean(unbiased)
      # mean pop?
    }

    unbiased.beta <- unbiased$beta
    PAR <- names(unbiased.beta)
    unbiased.COV <- unbiased$COV[PAR,PAR]
    unbiased.PRE <- PDsolve(unbiased.COV)

    PAR <- unique(c(names(debiased.beta),names(bias.beta)))
    debiased.beta <- copy.beta(debiased.beta)
    unbiased.beta <- copy.beta(unbiased.beta)
    debiased.PRE <- copy.PRE(debiased.PRE)
    unbiased.PRE <- copy.PRE(unbiased.PRE)

    PRE <- debiased.PRE + unbiased.PRE
    COV <- PDsolve(PRE)
    beta <- c( debiased.PRE %*% COV %*% debiased.beta + unbiased.PRE %*% COV %*% unbiased.beta )
    names(beta) <- PAR
  }

  # package into ctmm format
  M <- ctmm(beta=beta,COV=COV)
  # finish mu, sigma later

  return(M)
}
