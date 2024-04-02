# R refuses to call log.list for list objects outside of ctmm
Log <- function(CTMM,EST=NULL,variable="area",debias=TRUE,...)
{
  CTMM <- listify(CTMM)
  NAMES <- names(CTMM)
  CTMM <- log.list(CTMM,variable=variable,debias=debias,...)

  EST <- listify(EST)
  if(length(EST))
  {
    y <- log.list(EST,variable=variable,debias=debias,...)
    # merge with log(CTMM)
    # !!!!!!!!!!!!!!!!!!!!!!
    # !!!!!!!!!!!!!!!!!!!!!!
    # !!!!!!!!!!!!!!!!!!!!!!
  }

  CTMM <- data.frame(log=CTMM$log,VAR.log=CTMM$VAR.log)
  rownames(CTMM) <- NAMES

  INF <- CTMM$VAR.log==Inf
  CTMM$log[INF] <- 0

  return(CTMM)
}


mlog_ctmms <- function(x,variable="area",debias=TRUE,level.UD=0.95,...)
{
  for(i in 1:length(x))
  {
    object <- x[[i]]
    # propagate error uncertainty
    UERE.FIT <- object$error>0 & ( paste("error",names(object$error)) %in% dimnames(object)[[1]] )
    STUFF <- id.parameters(object,profile=FALSE,linear=FALSE,UERE.FIT=UERE.FIT)
    NAMES <- STUFF$NAMES
    parscale <- STUFF$parscale
    lower <- STUFF$lower
    upper <- STUFF$upper
    PAR <- get.parameters(object,NAMES)

    par <- dVAR <- array(0,length(variable)) # log(p)
    J <- matrix(0,length(variable),length(PAR)) # Jacobian matrix

    for(j in 1:length(variable))
    {
      if(variable[j]=="speed") # mean speed
      {
        fn <- function(p)
        {
          object <- set.parameters(object,p)
          speed_deterministic(object)
        }

        par[i] <- fn(PAR)
        J[i,] <- genD(par=PAR,fn=fn,lower=lower,upper=upper,parscale=parscale,Richardson=2,mc.cores=1)$grad
      }
      else if(variable[i]=="diffusion") # max diffusion rate
      {
        STUFF <- diffusion(object)
        par[i] <- STUFF$D
        J[i,names(STUFF$J)] <- STUFF$J
      }
      else if(variable[i]=="kinetic")
      {
        STUFF <- rand.speed(object)
        par[i] <- STUFF$MLE
        J[i,names(STUFF$J)] <- STUFF$J

        STUFF <- drift.speed(object)
        par[i] <- par[i] + STUFF$EST
        dVAR[i] <- STUFF$VAR
      }
    } # end for variables

    COV <- object$COV[PAR,PAR,drop=FALSE]
    COV <- J %*% COV %*% t(J)

    # additional variance from other parameters (speed)
    diag(COV) <- diag(COV) + dVAR

    # log transformation
    COV <- t(COV / par) / par
    par <- log(par)

    # debiasing
    if(debias)
    {
      # log chi^2 bias correction
      VAR <- diag(COV) # don't include features that shouldn't be there
      SUB <- features[VAR[features]<Inf]
      if(debias && length(SUB))
      {
        COR <- diag(1,nrow(COV))
        dimnames(COR) <- dimnames(COV)
        GOOD <- VAR>.Machine$double.eps
        if(any(GOOD)) { COR[GOOD,GOOD] <- stats::cov2cor(COV[GOOD,GOOD,drop=FALSE]) }
        # COR <- nant(COR,0)

        # diagonalize and log-chi^2 debias relevant parameter estimates
        EIGEN <- eigen(COV[SUB,SUB])
        dimnames(EIGEN$vectors) <- list(SUB,SUB)
        names(EIGEN$values) <- SUB

        # fix signs
        if(object$isotropic[1]) { PARS <- "major" } else { PARS <- c("major","minor") }
        # VAR goes in log numerator for unbiased chi^2 estimates: variance, diffusion, MS speed, ...
        for(i in 1:ncol(EIGEN$vectors)) { if(sum(EIGEN$vectors[PARS,i])<0) { EIGEN$vectors[,i] <- -EIGEN$vectors[,i] } }

        # transform to diagonalized basis with VARs in log numerator
        par[SUB] <- t(EIGEN$vectors) %*% par[SUB] # diagonalize parameters
        DOF <- 2/EIGEN$values # log-chi^2 VAR-DOF relation
        BIAS <- log_chi2_bias(DOF) # negative bias for log(chi^2) variates
        BIAS <- pmax(BIAS,log_chi2_bias(1)) # clamp to 1 DOF
        par[SUB] <- par[SUB] - BIAS # E[log(chi^2)] bias correction
        par[SUB] <- c(EIGEN$vectors %*% par[SUB]) # transform back (still under logarithm)

        # log-gamma variance (better than delta method)
        EIGEN$values <- trigamma(DOF/2)
        EIGEN <- EIGEN$vectors %*% (EIGEN$values * t(EIGEN$vectors))

        COR[SUB,SUB] <- stats::cov2cor(EIGEN)
        # diag(COR) <- nant(diag(COR),1)
        # COR <- nant(COR,0)
        VAR[SUB] <- diag(EIGEN)
        COV <- COR * outer(sqrt(VAR))
        COV <- nant(COV,0)
      }
    }

    x[[i]] <- list(par=par,COV=COV)
  }# end for objects

  # format consistently in arrays
  # format Infs
  # !!!
  # !!!
  # !!!

  return(x)
}


log_ctmms <- function(x,variable="area",debias=TRUE,level.UD=0.95,...)
{
  if(length(variable)>1) { return(mlog_ctmms(x,variable=variable,debias=debias,level.UD=level.UD,...)) }

  x <- listify(x)
  x <- import.variable(x,variable=variable,level.UD=level.UD)
  # list(ID=ID,AREA=AREA,DOF=DOF,variable=variable)
  y <- list()
  y$log <- log(x$AREA)

  # 2D correction made in import.variables
  if(!debias)
  { y$VAR.log <- 2/x$DOF }
  else
  {
    y$VAR.log <- trigamma(x$DOF/2)

    BIAS <- log_chi2_bias(x$DOF)
    y$log <- y$log - BIAS
  }

  return(y)
}


log_area <- function(x,variable="area",debias=TRUE,...)
{ log_ctmms(x,variable=variable,debias=debias,...) }


log_UD <- function(x,variable="area",debias=TRUE,level.UD=0.95,...)
{
  x <- listify(x)
  x <- lapply(x,function(y){summary(y,level.UD=level.UD,units=FALSE)})
  x <- log_area(x,debias=debias,...)
  return(x)
}


# speed
log_speed <- function(x,variable="speed",debias=TRUE,...)
{
  x <- listify(x)
  x <- import.variable(x,variable="speed",chi=TRUE)
  # list(ID=ID,AREA=AREA,DOF=DOF,variable=variable)
  y <- list()
  y$log <- log(x$AREA)

  # 2D correction made in import.variables
  if(!debias)
  { y$VAR.log <- 2/x$DOF/4 }
  else
  {
    y$VAR.log <- trigamma(x$DOF/2)
    y$VAR.log <- y$VAR.log/4

    BIAS <- digamma(x$DOF/2) - log(x$DOF/2)
    BIAS <- BIAS/2
    BIAS <- nant(BIAS,0)
    y$log <- y$log - BIAS
  }

  return(y)
}


Exp <- function(est,VAR.est=0,VAR=0,VAR.VAR=0,variable="area",debias=TRUE,level=0.95,units=TRUE,...)
{
  # convert from log-chi to log-chi^2
  R <- exp_log(est=est,VAR.est=VAR.est,VAR=VAR,VAR.VAR=VAR.VAR,...)
  est <- R$mu
  VAR <- R$VAR

  if(variable=="speed")
  { DOF <- chi.dof(est,est^2+VAR) }
  else # chi^2
  { DOF <- 2*est^2/VAR }

  if(debias)
  {
    BIAS <- digamma(DOF/2) - log(DOF/2)
    BIAS <- nant(BIAS,0)
    if(variable=="speed") { BIAS <- BIAS/2 }
    est <- est + BIAS
  }

  CI <- chisq.ci(est,VAR=VAR,level=level)

  # apply units and name
  CI <- summary_UD_format(CI,DOF/2,units=units)

  return(CI)
}


exp_log <- function(est,VAR.est=0,VAR=0,VAR.VAR=0,...)
{
  mu <- exp(est + VAR/2)
  # grad <- c(1,1/2) %o% mu
  VAR.mu <- (mu)^2*VAR.est + (mu/2)^2*VAR.VAR

  R <- list(mu=mu,VAR=VAR.mu)
  return(R)
}
