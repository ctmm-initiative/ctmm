NAMES.POP <- c("pop-low ","pop-ave ","pop-high")

# meta-analysis of chi^2 random variables with inverse-chi^2 prior
meta.chisq <- function(s,dof,level=0.95,level.pop=0.95,IC="AIC",boot=FALSE,error=0.01,debias=TRUE,robust=TRUE,precision=1/2,...)
{
  ROBUST <- TRUE # robust statistics on bootstrap # 'robust' is for population median versus mean
  #tol <- .Machine$double.eps^precision
  n <- length(s)

  # MLE fit and CI return
  DOF <- S <- NULL # S'
  fit <- function(s,complete=TRUE,par=NULL)
  {
    # test for DOF==Inf MLE
    S <- sum(dof*s)/sum(dof)
    INF <- sum( dof^2*(S-s)^2 - 2*dof*S^2 )<=0

    # AIC/BIC constants
    if(complete && !is.na(IC))
    {
      if(IC=="AIC")
      { dIC <- 2 }
      else if(IC=="BIC")
      { dIC <- log(n) }
      else if(IC=="LOOCV")
      { dIC <- 0 }
    }

    CI <- t(array(c(0,0,Inf),c(3,3))) # initialize confidence intervals (0,Inf)

    # trivial model selection / model fitting
    if(INF)
    {
      DOF <- Inf
    }
    else # non-trivial model selection for DOF<Inf
    {
      SUB <- 1:n
      # negative log-likelihood of DOF after profiling S'
      nloglike <- function(par,zero=0)
      {
        n <- length((1:n)[SUB])
        S <- par[1]
        DOF <- 1/par[2]

        if(DOF<=0 || S==0 || S==Inf) { return(Inf) }

        # beta-prime with chi^2 limit -- special functions in series.R
        zero <- zero + sum( ( (dof/2-1)*log(dof*s) + log(dof) )[SUB] )
        -sum( ( -lbetaplog(dof/2,DOF/2,s/S) - dof*(log(2*S)/2) - (dof*s)/(2*S)*log1pxdx((dof*s)/(DOF*S)) + zero/n )[SUB] )
      }

      # log-likelihood at DOF=Inf
      if(complete && !is.na(IC))
      {
        if(IC!="LOOCV")
        { LL0 <- -nloglike(c(S,0)) }
        else if(IC=="LOOCV")
        {
          LL0 <- 0
          for(i in 1:n)
          {
            # leave-one-out MLE
            SUB <- -i
            S <- sum( (dof*s)[SUB] ) / sum( dof[SUB] )
            # (CV) log-likelihood of left out
            SUB <- i
            LL0 <- LL0 - nloglike(c(S,0))
          }
          SUB <- 1:n
        }
      }

      # initial guesses for small N
      if(is.null(par))
      {
        S <- 1/mean(1/s)
        DOF <- 2/(S^2*stats::var(1/s))
        par <- c(S,1/DOF)
      }
      FIT <- optimizer(par,nloglike,lower=0,upper=Inf)
      par <- FIT$par
      LLN <- -FIT$value # log-likelihood at DOF MLE

      if(complete && !is.na(IC))
      {
        if(IC=="AIC")
        { dIC <- dIC - 2*LLN + 2*LL0 }
        else if(IC=="BIC")
        { dIC <- dIC - 2*LLN + 2*LL0 }
        else if(IC=="LOOCV")
        {
          LLN <- 0
          for(i in 1:n)
          {
            # leave-one-out MLE
            SUB <- -i
            par.i <- optimizer(par,nloglike,lower=0,upper=Inf)$par
            # CV log-likelihood of left out
            SUB <- i
            LLN <- LLN - nloglike(par.i)
          }
          SUB <- 1:n
          dIC <- -2*LLN + 2*LL0
          # IC <- "LOOCVIC"
        }
      }

      if(complete && !is.na(IC) && dIC>=0)
      {
        INF <- TRUE
        S <- sum(dof*s)/sum(dof)
        DOF <- Inf
      }
      else
      {
        S <- par[1]
        DOF <- 1/par[2]
      }
    } # end model fitting and selection

    # model selection cannot support population variation
    if(complete && INF) # Inf-N CIs
    {
      if(!is.na(IC))
      { warning("Population variation not supported (\u0394",IC,"=",dIC,").") }
      else
      { warning("Population variation not supported.") }

      S <- sum(dof*s)/sum(dof)
      dof <- sum(dof) # precision of point estimate
      DOF <- Inf # precision of population

      # mean area is the only area
      CI[1,] <- CI[2,] <- CI[3,] <- chisq.ci(S,DOF=dof,level=level,robust=robust)
      CI.DOF <- array(dof,3)
    }
    else if(complete) # finite-N CIs
    {
      H <- array(0,c(2,2))
      # easy derivatives
      dimnames(H) <- list(c("SDOF","DOF"),c("SDOF","DOF"))
      H['DOF','DOF'] <- -sum( psigamma(DOF/2,1) - psigamma((DOF+dof)/2,1) )/4
      H['DOF','SDOF'] <- H['SDOF','DOF'] <- sum( dof*s/(dof*s+DOF*S) )/(2*DOF*S)
      H['SDOF','SDOF'] <- -sum( 1/(DOF*S^2) - (DOF+dof)/(dof*s+DOF*S)^2 )/2

      COV <- cov.loglike(-H)
      J <- diag(1,2)
      # variables of interest
      dimnames(J) <- list(c("S","DOF"),c("SDOF","DOF"))
      J['S','SDOF'] <- 1/DOF
      J['DOF','SDOF'] <- 1/S
      COV <- J %*% COV %*% t(J)

      # delta method for population quantile uncertainty
      # population point estimates are inverse chi^2
      Q <- function(par=c(S,DOF)) { rev( par[1]/chisq.ci(1,DOF=par[2],level=level.pop,robust=robust) ) }
      GRAD <- genD(c(S,DOF),Q,lower=c(0,0),order=1)$gradient # [Q,(S,DOF)]
      VAR <- diag(GRAD %*% COV %*% t(GRAD)) # [Q,Q]
      # point estimates
      CI[,2] <- Q()
      # population low,point,high CIs
      for(i in 1:3) { CI[i,] <- chisq.ci(CI[i,2],COV=VAR[i],level=level) } # point estimate already robust from above
      CI.DOF <- 2*CI[,2]^2/VAR

      # replace (median) point estimate with inverse-chi^2 mean estimate
      if(!robust)
      {
        Q <- function(par=c(S,DOF)) { nant(par[2]/max(par[2]-2,0),1)*par[1] }
        GRAD <- genD(c(S,DOF),Q,lower=c(0,0),order=1)$gradient
        VAR <- c(GRAD %*% COV %*% GRAD)
        # mean area
        CI[2,] <- chisq.ci(Q(),COV=VAR,level=level)
        CI.DOF[2] <- 2*CI[2,2]^2/VAR
      }
    }
    else # point estimates only
    {
      # point estimates only
      CI[,2] <- rev( S/chisq.ci(1,DOF=DOF,level=level.pop,robust=robust) )
      # mean value of inverse-chi^2 for point estimate instead of harmonic mean value
      if(!robust) { CI[2,2] <- nant(DOF/max(DOF-2,0),1)*S }
    }

    if(complete)
    {
      rownames(CI) <- NAMES.POP
      colnames(CI) <- NAMES.CI

      R <- list(CI=CI,S=S,DOF=DOF,CI.DOF=CI.DOF) # S'
      return(R)
    }
    else
    {
      # just return the point estimates
      CI <- CI[,2]
      names(CI) <- NAMES.POP
      return(CI)
    }
  } # end fit function

  # MLE
  STUFF <- fit(s)
  CI <- STUFF$CI
  S <- STUFF$S
  DOF <- STUFF$DOF

  # there is no point in bootstrapping if DOF==Inf
  ## bias is in the opposite direction to give population variance
  ## and MLE is MVU

  if(!boot || DOF==Inf) # MLE
  { CI <- cbind(CI, DOF=STUFF$CI.DOF) }
  else # parametric bootstrap
  {
    S1 <- S2 <- 0 # sums
    AVE <- VAR <- 0 # cumulants: mean, variance
    ENSEMBLE <- NULL

    INF <- c(0,Inf,Inf,0)
    INF <- rbind(INF,INF,INF)
    rownames(INF) <- NAMES.POP
    colnames(INF) <- c(NAMES.CI,"DOF")

    N <- 0
    ERROR <- Inf
    pb <- utils::txtProgressBar(style=3)
    while(ERROR>=error || N<20)
    {
      b <- sapply(1:length(dof),function(i){stats::rbeta(1,dof[i]/2,DOF/2)}) # beta
      b <- b/(1-b) # beta'
      b <- (DOF*S)*(b/dof) # scaled beta'

      #DEBUG <<- list(fit=fit,b=b,dof=dof)
      SAMPLE <- fit(b,complete=FALSE,par=c(S,1/DOF))
      ENSEMBLE <- cbind(ENSEMBLE,SAMPLE)
      N <- ncol(ENSEMBLE)

      if(!ROBUST)
      {
        if(any(SAMPLE==Inf))
        {
          warning("Sampling distribution does not always resolve mean. Try robust=TRUE.")
          return(INF)
        }
        S1 <- S1 + SAMPLE
        S2 <- S2 + SAMPLE^2
      }
      else # keep ensemble sorted with insert sort
      {
        ENSEMBLE <- apply(ENSEMBLE,1,function(s){sort(s,method="radix")}) # [N,3]
        dim(ENSEMBLE) <- c(N,3) # R drops dimensions :(
        ENSEMBLE <- t(ENSEMBLE) # [3,N]
      }

      if(N>1)
      {
        if(!ROBUST) # mean variance
        {
          AVE <- S1/N # mean
          VAR <- S2/N - AVE^2 # variance
          VAR <- N/(N-1) * VAR # debias
          # max relative error
          ERROR <- max( sqrt(VAR/N)/AVE )
        }
        else # median & CI
        {
          # medians
          AVE <- mint(ENSEMBLE,(N+1)/2)
          # standard errors on the medians
          Q1 <- mint(ENSEMBLE,(N+1-sqrt(N))/2)
          Q2 <- mint(ENSEMBLE,(N+1+sqrt(N))/2)
          ERROR <- max(1-Q1/AVE,Q2/AVE-1)

          # correct for Inf AVE
          if(is.nan(ERROR)) { ERROR <- Inf }

          if(N>length(AVE)/error^2)
          {
            warning("Expectation values did not converge after ",N," iterations.")
            break
          }
        }
      }

      utils::setTxtProgressBar(pb,clamp(min(N/20,(error/ERROR)^2)))
    }
    # TODO close progressbar

    if(!ROBUST)
    {
      # more numerically precise cumulant calculation
      AVE <- apply(ENSEMBLE,1,mean)
      VAR <- apply(ENSEMBLE,1,stats::var)
      DOF <- 2*AVE^2/VAR # chi^2 DOF for CIs
      BIAS <- AVE/CI[,2] # first-order bias captured by multiplicative bias
      if(debias) { CI <- CI/BIAS } # first-order bias removed by this
      CI <- sapply(1:3,function(i){chisq.ci(CI[i,2],DOF=DOF[i],level=level,robust=robust)})
      CI <- t(CI)
    }
    else
    {
      alpha <- (1-level)/2
      Q <- apply(ENSEMBLE,1,function(en){stats::quantile(en,c(alpha,0.5,1-alpha))}) # [3,p]
      Q <- t(Q) # [p,3]

      BIAS <- Q[,2]/CI[,2] # first-order bias in the medians captured by this
      CI <- Q
      if(debias) { CI <- CI/BIAS } # first-order bias removed by this
      DOF <- sapply(1:nrow(ENSEMBLE),function(i){chisq.dof(Q[i,2],stats::IQR(ENSEMBLE[i,]))})
    }
    rownames(CI) <- NAMES.POP
    colnames(CI) <- NAMES.CI
    CI <- cbind(CI,DOF=DOF)
    close(pb)
  }

  # mean not estimated
  if(!robust && DOF<=2)
  {
    CI[2,1] <- 0 # low CI
    CI[2,4] <- 0 # DOF
    warning("Population mean not convergent. Consider robust=TRUE.")
  }

  return(CI)
}


# shrinkage estimator
shrink.chisq <- function(s,dof,S,DOF,...)
{
  #
  #
  #
}


# inverse of special function u = f(n) = log(n/2) - digamma(n/2) ~ 1/n
inv.special.F <- function(u,precision=1/2,...)
{
  tol <- .Machine$double.eps^precision
  ERROR <- Inf

  n <- 1/u/2 # really n/2
  while(ERROR>tol)
  {
    f <- log(n) - digamma(n)
    if(f==0) { return(Inf) }
    d <- 1/n - trigamma(n)
    n <- n + (u-f)/d

    ERROR <- abs(f-u)/u
  }
  n <- n*2 # n

  return(n)
}


meta <- function(x,level=0.95,level.UD=0.95,level.pop=0.95,robust=FALSE,IC="AIC",boot=FALSE,error=0.01,debias=TRUE,units=TRUE,plot=TRUE,...)
{
  meta.area(x=x,level=level,level.UD=level.UD,level.pop=level.pop,IC=IC,boot=boot,error=error,debias=debias,robust=robust,units=units,plot=plot,...)
}

# wrapper: meta-analysis of CTMM areas
# TODO range=FALSE ???
meta.area <- function(x,level=0.95,level.UD=0.95,level.pop=0.95,robust=FALSE,IC="AIC",boot=FALSE,error=0.01,debias=TRUE,units=TRUE,plot=TRUE,...)
{
  N <- length(x)
  CLASS <- class(x[[1]])

  # TODO do range or diffusion, but not both
  # TODO do range or diffusion, but not both
  # TODO do range or diffusion, but not both

  AREA <- DOF <- ID <- array(0,N)
  for(i in 1:N)
  {
    ID[i] <- x[[i]]@info$identity
    if(CLASS=="ctmm")
    {
      AREA[i] <- area.covm(x[[i]]$sigma)
      DOF[i] <- DOF.area(x[[i]])
    }
    else # UD or summary(UD)
    {
      if(CLASS=="UD") { x[[i]] <- summary(x[[i]],level.UD=level.UD,units=FALSE) }
      # now summary(UD) list
      DOF[i] <- x[[i]]$DOF['area']
      # convert to SI units
      UNITS <- rownames(x[[i]]$CI) # "area (units)"
      UNITS <- substr(UNITS,nchar("area (")+1,nchar(UNITS)-1)
      AREA[i] <- x[[i]]$CI[2] %#% UNITS
    }
  }

  # 2-dimensions
  DOF <- 2*DOF

  # level.UD coverage (e.g., 95% home ranges rather than straight variance)
  if(CLASS=="ctmm") { AREA <- -2*log(1-level.UD)*pi * AREA }

  # inverse-chi^2 population distribution
  CI <- meta.chisq(AREA,DOF,level=level,level.pop=level.pop,IC=IC,boot=boot,error=error,debias=debias,robust=robust)

  if(plot)
  {
    ID[N+1] <- ifelse(robust,"median","mean")
    AREA[N+1] <- CI[2,2]
    DOF[N+1] <- CI[2,4]

    # basic forest plot
    IND <- (N+1):1
    PLOT <- sapply(1:(N+1),function(i){chisq.ci(AREA[i],DOF=DOF[i],level=level)})
    PLOT <- t(PLOT) # [N+1,3]

    UNITS <- unit(PLOT,"area",SI=!units)
    PLOT <- PLOT/UNITS$scale

    xlab <- paste0(100*level.UD,"% Area (",UNITS$name,")")
    plot(range(PLOT),c(1,N+1),col=grDevices::rgb(1,1,1,0),xlab=xlab,ylab=NA,yaxt="n",...)
    graphics::axis(2,at=IND[1:N],labels=ID[1:N],las=2)
    graphics::axis(2,at=1,labels=ID[N+1],las=2,font=2)
    graphics::abline(v=PLOT[N+1,2],col=grDevices::rgb(0.5,0.5,0.5,0.5))
    graphics::points(PLOT[,2],IND,pch=16)
    suppressWarnings( graphics::arrows(PLOT[,1],IND,PLOT[,3],IND,length=0.05,angle=90,code=3,...) )
  }

  UNITS <- unit(CI[,1:3],"area",SI=!units)
  CI[,1:3] <- CI[,1:3]/UNITS$scale
  rownames(CI) <- paste0(rownames(CI)," area (",UNITS$name,")")

  # 2-dimensions
  CI[,4] <- CI[,4]/2

  return(CI[3:1,])
}


# wrapper: meta-analysis of CTMM speeds
meta.speed <- function(x,level=0.95,level.UD=0.95,level.pop=0.95,units=TRUE,...)
{
  N <- length(x)

  MSS <- DOF <- array(0,N)
  for(i in 1:N)
  {
    MSS[i] <- speed(x[[i]],prior=FALSE,units=FALSE)[2]^2 # mean square speed ~ chi^2
    DOF[i] <- DOF.speed(x[[i]])
  }

  # TODO substitute infinite speed (zero DOF) with 1s so to not produce NaNs
  #
  #

  CI <- meta.chisq(MSS,DOF,level=level,level.pop=level.pop,...)
  CI[,1:3] <- sqrt(CI[,1:3]) # root mean square speed

  UNITS <- unit(CI[,1:3],"speed",SI=!units)
  CI[,1:3] <- CI[,1:3]/UNITS$scale
  rownames(CI) <- paste0(rownames(CI)," speed (",UNITS$name,")")

  return(CI[3:1,])
}







# wrapper: meta-analysis of location covariances
meta.wishart.sigma <- function(x,...)
{
  N <- length(x)
  AXES <- length(x[[1]]$axes)

  SIGMA <- array(0,c(N,AXES,AXES))
  DOF <- array(0,N)
  for(i in 1:N)
  {
    # extract covariance matrix SIGMA and covariance of covariance matrix COV
    SIGMA[i,,] <- x[[i]]$sigma # Wishart matrix / n
    if(x[[i]]$isotropic) # chi^2
    { PAR <- 'major' }
    else # approximate Wishart DOF (exact if Wishart)
    { PAR <- c('major','minor') }
    EST <- SIGMA[[i]]@par[PAR]
    DOF[i] <- (2/AXES) * c(EST %*% PDsolve(x[[i]]$COV[PAR,PAR]) %*% EST) # average multiple DOFs if not Wishart
  }

  R <- meta.wishart(SIGMA,DOF,...)
  return(R)
}


# wrapper: meta-analysis of diffusion rates
meta.wishart.diff <- function(x,...)
{
  # TODO
}

# wrapper: meta-analysis of mean square velocities
meta.wishart.diff <- function(x,...)
{
  # TODO
}


####################
# population-level parameter estimates for Wishart - Inverse Wishart parameters
meta.wishart <- function(SIGMA,DOF,isotropic=FALSE,precision=1/2)
{
  tol <- .Machine$double.eps^precision
  N <- length(DOF)
  axes <- dimnames(SIGMA)[[2]]
  AXES <- length(axes)

  # EM algorithm
  nu <- sum(DOF) # initial estimate that seems reasonble
  S <- Reduce('+',vapply(1:N,function(i){ (nu+DOF[i])*DOF[i]*SIGMA[i,,,drop=FALSE] },diag(1,AXES))) / sum((nu+DOF)*DOF) # initial estimate (weighted average close to perturbative solution)

  count <- Inf
  while(count>2)
  {
    count <- 0
    L <- -Inf
    ERROR <- Inf
    while(ERROR>tol) # iterative weighted average (UNTESTED)
    {
      S <- Reduce('+', vapply(1:N,function(i) { S %*% PDsolve((nu*S+DOF[i]*SIGMA[i,,,drop=FALSE])/(nu+DOF[i])) %*% S },S) )/N
      L.OLD <- L
      L <- N/2*log(det(S)) - sum( vapply(1:N,function(i){ (nu+DOF[i]) * log(det(nu*S+DOF[i]*S[i,,])) },1) )/2
      ERROR <- L-L.OLD
      count <- count + 1
    }

    L <- -Inf
    ERROR <- Inf
    while(ERROR>tol) # Newton Raphson
    {
      L.OLD <- L
      CONST <- N/2*log(det(nu*S)) - Reduce('+',vapply(1:N,function(i){ log(det(nu*S+DOF[i]*SIGMA[i,,,drop=FALSE])) },S))/2
      L <- nu*CONST - N/2*mpsigamma(nu/2,deriv=-1,dim=AXES) + sum(vapply(1:N,function(i){ mpsigamma((nu+DOF[i])/2,deriv=-1,dim=AXES) },1))/2
      L1 <- CONST - N/2*mpsigamma(nu/2,deriv=0,dim=AXES) + sum(vapply(1:N,function(i){ mpsigamma((nu+DOF[i])/2,deriv=0,dim=AXES) },1))/2
      L2 <- - N/2*mpsigamma(nu/2,deriv=1,dim=AXES) + sum(vapply(1:N,function(i){ mpsigamma((nu+DOF[i])/2,deriv=1,dim=AXES) },1))/2
      nu <- clamp(nu-L1/L2,0,Inf)
      ERROR <- L-L.OLD
      count <- count + 1
    }
  } # end alternating EM algorithm

  like <- function(par)
  {
    ## data
    # SIGMA
    # DOF

    ## prior parameters - Inverse Wishart for SIGMA
    nu <- par[1]
    S <- covm(par[-1],axes=axes,isotropic=isotropic)

    ## posterior parameters
    sigma <- DOF*SIGMA + nu*array(S,c(AXES,AXES,N))
    dof <- DOF + nu # posterior t-dof
    sigma <- vapply(1:N,function(i){ sigma[,,i] / sigma[i] },S) # [AXES,AXES,N]
    sigma <- aperm(sigma,c(3,1,2)) # [N,AXES,AXES] # posterior t-sigma

    ## posterior likelihood - t
    log.det.S <- log(det(S))
    R <- sum( vapply(1:N,function(i){ (nu+DOF[i]) * ( log.det.S - log(det(nu*S+DOF[i]*S[i,,])) ) },1) )/2
    R <- R - N/2*mpsigamma(nu/2,deriv=-1,dim=AXES) + sum(vapply(1:N,function(i){ mpsigamma((nu+DOF[i])/2,deriv=-1,dim=AXES) },1))/2
    return(R)
  }

  ################
  # covariance matrix of hyper-parameter estimates
  par <- nu
  parscale <- 1
  lower <- 0
  upper <- 0
  NAMES <- 'nu'

  S <- covm(S,axes=axes,isotropic=isotropic)
  par <- c(par,S@par[1])
  parscale <- c(parscale,S@par[1])
  lower <- c(lower,0)
  upper <- c(upper,Inf)
  NAMES <- c(NAMES,"major")

  if(!isotropic)
  {
    par <- c(par,S@par[2:3])
    parscale <- c(parscale,S@par[1],pi/2)
    lower <- c(lower,0,-Inf)
    upper <- c(upper,0,Inf)
    NAMES <- c(NAMES,c('minor','angle'))
  }

  LIKE <- like(par)
  DIFF <- genD(par=par,fn=like,zero=LIKE,lower=lower,upper=upper,parscale=parscale,Richardson=2,mc.cores=1)
  hess <- -DIFF$hessian  # negative log likelihood
  grad <- -DIFF$gradient # negative log likelihood

  # more robust covariance calculation than straight inverse
  COV <- cov.loglike(hess,grad)
  dimnames(COV) <- list(NAMES,NAMES)
  COV <- COV[-1,-1] # nu is a nuisance parameter

}


# population-level parameter estimates for normally distributed parameters and parameter uncertainties
meta.normal <- function(MU,SIGMA,debias=TRUE,isotropic=FALSE,precision=1/2)
{
  N <- dim(MU)[1]
  DIM <- dim(MU)[2]

  tol <- .Machine$double.eps^precision
  REML <- debias

  # initial guesses
  mu <- colMeans(MU)
  sigma <- 0
  for(i in 1:N) { sigma <- sigma + outer(MU[i,]-mu) }
  sigma <- sigma/(N-REML)

  ERROR <- Inf
  loglike <- loglike.old <- -Inf
  while(ERROR>tol && loglike>=loglike.old)
  {
    loglike.old <- loglike
    sigma.old <- sigma

    # estimate mu exactly
    P <- array(0,c(N,DIM,DIM))
    mu <- P.mu <- 0
    for(i in 1:N)
    {
      P[i,,] <- PDsolve(sigma + SIGMA[i,,])
      P.mu <- P.mu + P[i,,]
      mu <- mu + c(P[i,,] %*% MU[i,])
    }
    COV.mu <- PDsolve(P.mu)
    mu <- c(COV.mu %*% mu)

    loglike <- REML/2*log(det(COV.mu)) + DIM*(N-REML)/2*log(2*pi)
    # gradient with respect to sigma
    RHS <- 0
    LHS <- P.mu
    for(i in 1:N)
    {
      D <- mu - MU[i,]
      RHS <- RHS + (P[i,,] %*% outer(D) %*% P[i,,])
      if(debias) { LHS <- LHS - (P[i,,] %*% COV.mu %*% P[i,,]) }
      loglike <- loglike - 1/2*log(det(sigma + SIGMA[i,,])) - 1/2*c(D %*% P[i,,] %*% D)
    }

    K <- sqrtm.covm(covm(sigma)) %*% mpow.covm(covm(LHS),-1/2)
    sigma <- K %*% RHS %*% t(K)

    # Standardized error
    ERROR <- sigma - sigma.old # absolute error
    K <- mpow.covm(covm(sigma),-1/2)
    ERROR <- K %*% ERROR %*% K # standardize to make ~1
    ERROR <- ERROR %*% ERROR # square to make positive
    ERROR <- sqrt(mean(diag(ERROR %*% ERROR))) # error in standard deviations
  }

  DUP <- upper.tri(sigma,diag=TRUE)

  # we still need hessian(sigma) for sigma CIs
  log.like <- function(par,REML=TRUE)
  {
    sigma <- array(0,c(DIM,DIM))
    sigma[DUP] <- par
    sigma <- t(sigma)
    sigma[DUP] <- par

    # sum up log-likelihood
    loglike <- REML/2*log(det(COV.mu)) + DIM*(N-REML)/2*log(2*pi)
    for(i in 1:N)
    {
      D <- mu - MU[i,]
      loglike <- loglike - 1/2*log(det(sigma + SIGMA[i,,])) - 1/2*c(D %*% P[i,,] %*% D)
    }
  }

  par <- sigma[DUP]

  parscale <- sqrt( diag(sigma) )
  parscale <- sigma %o% sigma
  parscale <- parscale[DUP]

  lower <- array(-Inf,c(DIM,DIM))
  diag(lower) <- 0

  DIFF <- genD(par,log.like,parscale=parscale,lower=lower,upper=Inf)
  COV.sigma <- PDsolve(-DIFF$HESS)

  loglike <- log.like(par,REML=FALSE)

  # AIC
  n <- N
  q <- DIM
  k <- 1
  nu <- (q^2+q)/2
  K <- q*k + nu

  AIC <- 2*K - 2*loglike
  # this is exact if !isotropic
  AICc <- q*(n-k)*2*K/(q*n-K-nu) - 2*loglike
  BIC <- K*log(N) - 2*loglike

  return(list(mu=mu,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike,AIC=AIC,AICc=AICc,BIC=BIC,isotropic=isotropic))
}


###########
# log-transformed parameters
# debias includes bias correction for chi^2 to log(chi^2)
# matrix casts location covariance, diffusion rate, velocity covariance all as distinct matrices for above bias correction
log.ctmm <- function(CTMM,features,debias=TRUE)
{
  isotropic <- CTMM$isotropic
  par <- get.parameters(CTMM,features)
  COV <- CTMM$COV

  ### log transform all positive parameters
  # features to log transform
  COV.NAMES <- dimnames(COV)[[1]]
  SUB <- features[(features %in% COV.NAMES) & (features %in% POSITIVE.PARAMETERS) & (par>0)]

  # Jacobian for log transformation
  J.new <- function()
  {
    J <- diag(1,nrow(COV))
    dimnames(J) <- dimnames(COV)
    return(J)
  }

  # log transform positive parameters
  J <- J.new()
  for(s in SUB)
  {
    J[s,s] <- 1/par[s]
    par[s] <- log(par[s])
  }
  # transform covariance (from logarithms)
  COV <- J %*% COV %*% t(J)

  # finish logarithm of sigma matrix (and not just eigen values)
  if(!isotropic)
  {
    angle <- par['angle']
    par['angle'] <- 0 # off-diagonal of log(sigma) after rotation by -angle

    J <- J.new()
    J['angle','angle'] <- par['major'] - par['minor']
    COV <- J %*% COV %*% t(J)
  }

  # log chi^2 bias correction
  if(debias)
  {
    # diagonalize and log-chi^2 debias relevant parameter estimates
    EIGEN <- eigen(COV[SUB,SUB])
    dimnames(EIGEN$vectors) <- list(SUB,SUB)
    names(EIGEN$values) <- SUB
    # fix signs
    if(isotropic) { VAR <- "major" } else { VAR <- c("major","minor") }
    # VAR goes in log numerator for chi^2 variates: variance, diffusion, MS speed, ...
    for(i in 1:nrow(EIGEN$vectors)) { if(sum(EIGEN$vectors[i,VAR])<0) { EIGEN$vectors[i,] <- -EIGEN$vectors[i,] } }
    # transform to diagonalized basis with VARs in log numerator
    par[SUB] <- t(EIGEN$vectors) %*% par[SUB] # diagonalize parameters
    DOF <- 2/EIGEN$values # log-chi^2 VAR-DOF relation
    BIAS <- digamma(DOF/2)-log(DOF/2) # negative bias for log(chi^2) variates
    if(!isotropic)
    {
      # some of the eigen parameter is orientation - which is ~Gaussian and not ~log chi^2 (no bias)
      OVER <- abs(EIGEN$vectors['angle',]) # overlap with orientation eigen-parameter
      BIAS <- (1-OVER)*BIAS # first-order correction (could start at second-order?)
    }
    par[SUB] <- par[SUB] + BIAS # log-chi^2 bias correction
    par[SUB] <- EIGEN$vectors %*% par[SUB] # transform back (still under logarithm)
  }

  # un-diagonalize log(sigma)
  if(!CTMM$isotropic)
  {
    u <- c(cos(angle),sin(angle))
    v <- c(-sin(angle),cos(angle))
    NAMES <- c("major","minor","angle") # input
    UP <- c(1,4,2) # "xx","yy","xy" # upper triangle of log(sigma) # output

    par[NAMES] <- par['major']*(u%o%u)[UP] + par['minor']*(v%o%v)[UP]

    J <- J.new()
    J[NAMES,NAMES] <- cbind( (u%o%u)[UP], (v%o%v)[UP], (u%o%v+v%o%u)[UP] )
    COV <- J %*% COV %*% t(J)
  }

  # fill out missing VAR with Inf after transformation --- prevent NaNs
  TCOV <- diag(Inf,length(features))
  dimnames(TCOV) <- list(features,features)
  NAMES <- dimnames(COV)[[1]]
  TCOV[NAMES,NAMES] <- COV[NAMES,NAMES]
  COV <- TCOV

  return(list(par=par,COV=COV))
}

# orthogonal transformation on matrix -> linear transformation on vector
# t(O) %*% M %*% O -> L %*% m
orth2lin <- function(O,sym=TRUE)
{
  N <- dim(O)[1]
  L <- array(0,c(N,N,N,N))

  for(i in 1:N)
  {
    for(j in 1:N)
    {
      IN <- array(0,c(N,N))
      IN[i,j] <- 1
      OUT <- t(O) %*% IN %*% O
      L[,,i,j] <- OUT
    }
  }

  dim(L) <- c(N^2,N^2)
  return(L)
}

#####################
# inverse transformation of above
exp.ctmm <- function(object,debias=TRUE)
{
  mu <- object$mu # mean log parameters (log chi^2)
  COV.mu <- object$COV.mu # uncertainty in mean log parameters
  sigma <- object$sigma # dispersion of mean logs (determines chi^2 DOFs)
  COV.sigma <- object$COV.sigma # uncertainty in dispersion of mean logs

  isotropic <- object$isotropic
  NAMES <- names(mu)
  N <- length(mu)

  J.new <- function()
  {
    J <- diag(1,N)
    dimnames(J) <- list(NAMES,NAMES)
    return(J)
  }

  # diagonalize log(sigma)
  if(!isotropic)
  {
    SIGMA <- matrix(sigma[c("major","angle","angle","minor")],c(2,2))
    EIGEN <- eigen(SIGMA)
    U <- EIGEN$vectors

    # t(U) %*% SIGMA %*% U
    mu[c('major','minor')] <- EIGEN$values
    mu['angle'] <- 0
    angle <- U[,1] # major axis
    angle <- sign(angle[1]) * angle # positive x component
    angle <- atan2(angle[2],angle[1])

    # transformation matrix J: "xx","yy","xy" -> "major","minor","0-off"
    J <- J.new()
    SIG <- c("major","minor","angle")
    J[SIG,SIG] <- orth2lin(U)[c(1,2,4),c(1,2,4)]

    # linear transform on major, minor, angle
    COV.mu <- J %*% COV.mu %*% t(J)
    sigma <- J %*% sigma %*% t(J)

    # same kind of thing but in even more dimensions
    LJ <- orth2lin(t(J))
    COV.sigma <- LJ %*% COV.sigma %*% LJ
  }

  # reverse bias correction
  if(debias)
  {
    EIGEN <- eigen(COV.mu)
    names(EIGEN$values) <- NAMES
    dimnames(EIGEN$vectors) <- list(NAMES,NAMES)
    # fix signs
    if(isotropic) { VAR <- "major" } else { VAR <- c("major","minor") } # these go in the numerator before log
    for(i in 1:nrow(EIGEN$vectors)) { if(sum(EIGEN$vectors[i,VAR])<0) { EIGEN$vectors[i,] <- -EIGEN$vectors[i,] } }
    # transform to diagonalized basis with VARs in log numerator
    mu <- t(EIGEN$vectors) %*% mu
    DOF <- 2/EIGEN$values # log-chi^2 VAR-DOF relation
    BIAS <- digamma(DOF/2)-log(DOF/2) # negative bias for log(chi^2) variates
    if(!isotropic)
    {
      # some of the eigen parameter is orientation - which is ~Gaussian and not ~log chi^2 (no bias)
      OVER <- abs(EIGEN$vectors['angle',]) # overlap with orientation eigen-parameter
      BIAS <- (1-OVER)*BIAS # first-order correction (could start at second-order?)
    }
    mu <- mu - BIAS # log-chi^2 bias addition
    mu <- EIGEN$vectors %*% mu # transform back (still under logarithm)
  }

  # exp transformation
  SUB <- NAMES %in% POSITIVE.PARAMETERS
  mu[SUB] <- exp(mu[SUB])

  J <- J.new()
  for(s in SUB) { J[s,s] <- mu[s] }

  COV.mu <- J %*% COV.mu %*% t(J) # delta method
  sigma <- J %*% sigma %*% t(J) # not exactly sure about this?
  LJ <- orth2lin(t(J))                 # would be true if above is true
  COV.sigma <- LJ %*% COV.sigma %*% t(LJ) # would be true if above is true

  # undiagonalize sigma-location
  if(!isotropic)
  {
    mu['angle'] <- angle

    # transform uncertainty
    J <- J.new()
    J['angle','angle'] <- 1/(mu['major']-mu['minor'])

    COV.mu <- J %*% COV.mu %*% t(J)

    # transform sigma-par
    sigma <- J %*% sigma %*% t(J)
    LJ <- orth2lin(t(J))
    COV.sigma <- LJ %*% COV.sigma %*% t(LJ)
  }

  return(list(mu=mu,COV.mu=COV.mu,sigma=sigma,COV.sigma=COV.sigma))
}


##########
mean.ctmm.mu <- function(x,debias=TRUE,isotropic=FALSE,...)
{
  axes <- x[[1]]$axes
  AXES <- length(axes)
  N <- length(x)

  # Gaussian-Gaussian in all cases
  MU <- array(0,c(N,AXES))
  SIGMA <- array(0,c(N,AXES,AXES))
  for(i in 1:N)
  {
    MU[i,] <- x[[i]]$mu
    SIGMA[i,,] <- x[[i]]$COV.mu

    # TODO !!!
    # fill in with zeroes for non-stationary means
    # TODO !!!
  }
  R <- meta.normal(MU,SIGMA,debias=debias)
  # R$mu # mean of means
  # R$COV.mu # uncertainty in mean of means estimate
  # rename these
  names(R)[ which(names(R)=="sigma") ] <- "PCOV.mu" # dispersion of means
  names(R)[ which(names(R)=="COV.sigma") ] <- "COV.PCOV.mu" # uncertainty in dispersion of means
  return(R)
}


#############
mean.ctmm.features <- function(x,sufficient="log-normal",prior="log-normal",method="exact",debias=TRUE,precision=1/2,...)
{
  sufficient <- match.arg(summary,c("Wishart","chisq","log-normal"))
  prior <- match.arg(prior,c("Inverse-Wishart","log-normal"))
  method <- match.arg(method,c("exact","Laplace","MCMC"))

  axes <- x[[1]]$axes
  AXES <- length(axes)
  N <- length(x)

  features <- unique( sapply(x,function(y){y$features}) )

  # analyticlly solvable
  if(method %in% c("exact","Laplace") && sufficient=="log-normal" && prior=="log-normal")
  {
    DIM <- length(features)
    MU <- array(0,c(N,DIM))
    SIGMA <- array(0,c(N,DIM,DIM))

    for(i in 1:N)
    {
      R <- log.ctmm(x[[1]],features,debias=debias)
      MU[i,] <- R$par
      SIGMA[i,,] <- R$COV
    }

    # aggregate log parameters
    R <- meta.normal(MU,SIGMA,debias=debias)
    # transform results back
    R <- exp.ctmm(R,debias=debias)
    names(R)[ which(names(R)=="mu") ] <- "par" # mean parameters
    names(R)[ which(names(R)=="COV.mu") ] <- "COV" # mean parameter uncertinaty
    names(R)[ which(names(R)=="sigma") ] <- "PCOV" # parameter dispersion
    names(R)[ which(names(R)=="COV.sigma") ] <- "COV.PCOV" # parameter dispersion uncertainty
    R <- set.parameters(R,R$par)
  }
  else if(method=="exact" && sufficient=="Wishart" && prior=="Inverse-Wishart")
  {
    #
    #
    #
  }
  else
  {
    #
    #
    #
  }

  return(R)
}


###########
mean.ctmm <- function(x,...)
{
  isotropic <- FALSE # for now

  MU <- mean.ctmm.mu(x,debias=debias,isotropic=isotropic,...)
  CTMM <- mean.ctmm.features(x,debias=debias,isotropic=isotropic)

  # copy features over
  CTMM$mu <- MU$mu
  CTMM$COV.mu <- MU$COV.mu
  CTMM$PCOV.mu <- MU$PCOV.mu
  CTMM$COV.PCOV.mu <- MU$COV.PCOV.mu

  # add log-likelihoods
  CTMM$loglike <- CTMM$loglike + MU$loglike
  CTMM$AIC <- CTMM$AIC + MU$AIC
  CTMM$AICc <- CTMM$AICc + MU$AICc
  CTMM$BIC <- CTMM$BIC + MU$BIC

  # TODO: how do we structure the 3 isotropic types

  # return final result
  return(CTMM)
}


mean.select <- function(x,IC="AICc",...)
{
  # is location mean isotropic
  #
  # is location covariance isotropic
  #
  # is covariance of location covariance isotropic
}
