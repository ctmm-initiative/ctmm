summary.meta <- function(object,IC="AICc",...)
{
  if(class(object)[1]=="list")
  { return(summary.meta.list(object,IC=IC,...)) }
  else
  { return(summary.meta.single(object,...)) }
}

summary.meta.list <- function(object,IC="AICc",...)
{

}

summary.meta.single <- function(object,...)
{

}


# meta-analysis of chi^2 random variables with inverse-Gaussian prior
meta.chisq <- function(s,dof,level=0.95,level.pop=0.95,IC="AICc",method='mle',boot=FALSE,iterate=FALSE,error=0.01,debias=TRUE,precision=1/2,CI.FN="chisq",...)
{
  # discard null estimates
  ZERO <- dof<=.Machine$double.eps
  s <- s[!ZERO]
  dof <- dof[!ZERO]

  # tol <- .Machine$double.eps^precision
  n <- length(s)
  # for model ICs
  dIC <- Inf

  # Bessel's correction
  BC <- n/(n-debias)

  # % chi^2 versus IG in sampling distribution of mu (very approximate)
  DRATIO <- 0 # (0,1) -> (chi^2,IG)

  # kappa for debiased standard deviation for CIs
  k.std <- function(par)
  {
    mu <- par[1]
    k <- par[2] # must already be debiased
    if(k==0) { return(0) }
    theta <- n/(k*mu) # mu sampling distribution theta
    BIAS <- sqrt(2/pi*theta)*besselK(theta,1,expon.scaled=TRUE)
    BIAS <- BIAS * exp(lgamma((n+1)/2) - lgamma(n/2)) / sqrt(n/2)
    k <- k / BIAS^2
    return(k)
  }

  # this is MVU for both chi^2 and IG
  # this is 1st order debiased in between
  inverse.mean <- function(mu,Vm2,MIN=2)
  {
    dof <- 2/Vm2
    dof <- max(dof,MIN)
    1/mu * (1-debias*2/dof)
  }
  # Vm2 = VAR[mu]/mu^2

  ################
  # MLE fit and CI return
  fit.mle <- function(s,par=NULL)
  {
    MODELS <- 2
    complete <- is.null(par)
    NAME <- rep(NA_character_,MODELS)
    PAR <- t(array(c(NA_real_,0,1),c(3,MODELS)))
    LL <- rep(-Inf,MODELS)

    ####################
    # parameter estimation

    # SUB <- 1:n
    # negative log-likelihood
    DOF.LO <- dof < 1/.Machine$double.eps
    N.LO <- sum(DOF.LO)
    DOF.HI <- !DOF.LO
    N.HI <- sum(DOF.HI)
    nloglike <- function(par,zero=0)
    {
      nll <- 0

      mu <- par[1] # mean
      k <- ifelse(length(par)>=2,par[2],0) # IG var/mu^3
      theta <- 1/(k*mu) # delta at Inf
      rho <- ifelse(length(par)>=3,par[3],1) # GIG shape (IG at 1)

      if(k<0 || k==Inf || mu<=0 || mu==Inf) { return(Inf) }

      if(theta==Inf) # chi^2-delta == chi^2
      {
        zero <- zero + sum( -log(s) + dof/2*log(dof/2)  - lgamma(dof/2) )

        r <- s/mu
        nll <- - sum( +dof/2*log(r) - dof/2*r + zero/n )

        return(nll)
      }

      # moderate to large dof
      if(any(DOF.LO)) # chi^2-IG
      {
        I <- DOF.LO
        dz <- sum( dof[I]/2*log(dof[I]/2) - log(s[I]) - lgamma(dof[I]/2) )

        alpha <- dof[I]*s[I]/mu
        beta <- alpha/theta
        beta <- 1/2*log1p(beta) # log(sqrt(1+dof*s*k))
        nll <- nll - sum( -dof[I]/2*log(mu/s[I]) - (dof[I]+rho)/2*beta + lKK(rho,dof[I],theta,alpha) + (zero/n+dz/N.LO) )
      }

      # extremely large dof case
      if(any(DOF.HI)) # delta-IG -- IG
      {
        I <- DOF.HI
        dz <- sum( -log(s[I]) - log(2) )

        r <- s[I]/mu
        nll <- nll - sum( -theta/2*(r+1/r-2) - BesselK(theta,rho/2,expon.scaled=TRUE,log=TRUE) - rho/2*log(r) + (zero/n+dz/N.HI) )
      }

      return(nll)
    }

    # n.min <- 1-2 (AIC-AICc)
    if(complete)
    {
      NAME[1] <- "Dirac-\u03B4"
      w <- dof/sum(dof)
      mu <- sum(w*s) # exact solution
      mu <- nant(mu,mean(s)) # if any dof==Inf
      PAR[1,1] <- mu
      LL[1] <- -nloglike(mu)
    }

    # n.min <- 2-4 (AIC-AICc)
    NAME[2] <- "inverse-Gaussian"
    if(complete && n>=2) # par==NULL
    {
      w <- 1-exp(-dof) # turn off low dof from mean
      w <- w/sum(w)
      mu <- sum(w*s)
      k <- sum(w/s) - 1/mu
      par <- c(mu,k)
    }
    if((complete && n>=2) || length(par)==2)
    {
      FIT <- optimizer(par,nloglike,lower=c(0,0),upper=c(Inf,Inf))
      PAR[2,1:2] <- FIT$par
      LL[2] <- -FIT$value # log-likelihood at MLE
      # Bessel's correction in log-likelihood
      if(debias) { LL[2] <- -nloglike(c(1,BC)*FIT$par) }
      # but don't debias PAR until after COV calculation
    }

    # # unknown shape parameter rho (generalized inverse Gaussian)
    # # n.min <- 3-6
    # NAME[3] <- "generalized IG"
    # if(complete && n>=3)
    # {
    #   eta <- FIT$par[1]
    #   theta <- 1/FIT$par[2]
    #   rho <- 1
    #   par <- c(eta,1/theta,rho)
    # }
    # if((complete && n>=3) || length(par)==3)
    # {
    #   FIT <- optimizer(par,nloglike,lower=c(0,0,-Inf),upper=c(Inf,Inf,Inf))
    #   PAR[3,] <- FIT$par
    #   LL[3] <- -FIT$value # log-likelihood at MLE
    # }

    if(complete && !is.na(IC))
    {
      ICS <- c("AIC","AICc","BIC")
      dIC <- matrix(0,MODELS,length(ICS))
      colnames(dIC) <- ICS

      k <- c(1,1,1)[1:length(LL)] # mean parameters
      nu <- c(0,1,2)[1:length(LL)] # variance parameters
      K <- k + nu # total parameters
      dIC[,"AIC"]  <- 2*K - 2*LL
      dIC[,"AICc"] <- 2*K*nant(ifelse(debias,n-k,n)/pmax(n-K-nu,0),1) - 2*LL
      dIC[,"BIC"]  <- log(n)*K - 2*LL

      dIC <- cbind(dIC[,IC])
      rownames(dIC) <- NAME
      colnames(dIC) <- paste0("\u0394",IC)
      IND <- sort(c(dIC),index.return=TRUE)$ix
      dIC <- dIC[IND,,drop=FALSE]
      dIC <<- dIC # store to main function
      print(dIC - min(dIC))
    }
    else # unknown DOF (no model selection)
    { IND <- 2  }

    ##################
    # propagate uncertainty
    CI <- t(array(c(0,0,Inf),c(3,4))) # initialize confidence intervals (0,Inf)

    IND <- IND[1] # best model
    PAR <- PAR[IND,] # best model parameter estimates
    par <- PAR[1:IND] # canonical subset of parameters

    if(IND==1) # Dirac delta -- DOF==Inf
    {
      dof <- sum(dof) # precision of point estimate

      # mean area is the only area
      CI[1,] <- chisq.ci(PAR[1],DOF=dof,level=level)
      CI.VAR <- 2*PAR[1]^2/dof

      # inverse-mean area
      CI[2,] <- 1/CI[1,] # not used
      CI[2,2] <- CI[2,2] * dof/max(dof-debias,1) # inverse-chi^2 mean bias correction
      CI.VAR[2] <- 2*CI[2,2]^2/max(dof-3*debias,1)

      # CoV^2 (RVAR)
      CI[3,] <- c(0,0,Inf)
      CI.VAR[3] <- Inf

      # CoV   (RSTD)
      CI[4,] <- c(0,0,Inf)
      CI.VAR[4] <- Inf
    }
    else if(IND==2) # inverse-Gaussian model
    {
      STUFF <- genD(par,nloglike,lower=c(0,0),order=2)
      COV <- cov.loglike(STUFF$hessian,STUFF$gradient,WARN=is.na(IC))

      if(debias) # Bessel's correction to point estimate of k=1/lambda and COV
      {
        par[2] <- BC * par[2]
        COV[2,] <- BC * COV[2,] # COV[BC*kappa] = BC COV[kappa] BC
        COV[,2] <- BC * COV[,2]
        COV[1,] <- sqrt(BC) * COV[1,] # VAR[mu] ~ k/n
        COV[,1] <- sqrt(BC) * COV[,1]
        PAR[1:IND] <- par
      }

      CI.VAR <- COV[1,1]

      # approximate ratio of chi^2 to IG behavior
      DRATIO <<- DD.IG.ratio(par,CI.VAR[1],n=n)

      CI[1,] <- chisq.IG.ci(par[1],VAR=CI.VAR[1],w=DRATIO,level=level,precision=precision)

      # mean inverse area of IG random variable # IG bias corrected
      CI[2,] <- 1/CI[1,] # these numbers aren't really used
      CI[2,2] <- inverse.mean(par[1],Vm2=COV[1,1]/par[1]^2) # debiased inverse
      # debiased relative variance at extremes
      CI.VAR[2] <- 1/(par[1]^2/COV[1,1] - 2*(1-DRATIO)) + DRATIO*COV[1,1]/par[1]*COV[2,2]/par[2]
      CI.VAR[2] <- CI.VAR[2]/par[1]^2

      # CoV^2 (RVAR)
      # CI[5,2] <- par[1]*par[2]
      GRAD <- par[2:1]
      CI.VAR[3] <- GRAD %*% COV %*% GRAD
      CI[3,] <- chisq.ci(par[1]*par[2],VAR=CI.VAR[3],level=level)

      # CoV (RSTD)
      CI[4,] <- sqrt(CI[3,])
      CI.VAR[4] <- CI.VAR[3]/CI[3,2]/4
      if(debias && 0<CI[3,2]) # sqrt bias
      {
        BIAS <- 1 + CI.VAR[3]/CI[3,2]^2/4
        CI[4,] <- CI[4,] * sqrt(BIAS)
        CI.VAR[4] <- CI.VAR[4] * BIAS
      }
    }
    # else if(IND==3) # GIG distribution
    # {
    #   STUFF <- genD(par,nloglike,lower=c(0,0,-Inf),order=2)
    #   COV <- cov.loglike(STUFF$hessian,STUFF$gradient)
    #
    #   Q <- function(par) { GIG.CI(par[1],1/par[2],par[3],level=level.pop,precision=precision) }
    #   GRAD <- genD(par,Q,lower=c(0,0,-Inf),order=1)$gradient
    #   VAR <- diag(GRAD %*% COV %*% t(GRAD))
    #
    #   # same approximation as above... not sure what the GIG sampling distribution is, though... probably not IG or GIG
    #   CI[,2] <- Q(par)
    #   for(i in 1:3) { CI[i,] <- ci.chisq.GIG(CI[i],VAR[i],eta=PAR[IND,1],theta=1/PAR[IND,2],rho=PAR[IND,3]) }
    #   CI <- t(CI) # [pop,est]
    #
    #   CI.DOF <- 2*CI[,2]^2/VAR # DOFs for F-test
    # }

    rownames(CI) <- c("mean","inverse-mean","CoV\u00B2 (RVAR)","CoV  (RSTD)")
    colnames(CI) <- NAMES.CI

    if(!complete) { return(CI[,2]) }

    R <- list(CI=CI,CI.VAR=CI.VAR,mu=PAR[1],k=PAR[2],rho=PAR[3])
    return(R)
  } # end IG/GIG fitting function

  ###############
  fit.blue <- function(s,par=NULL)
  {
    complete <- is.null(par)
    VAR <- 2*s^2/dof
    BFIT <- meta.normal(s,VAR,debias=debias)

    # IG parameters
    mu <- BFIT$mu
    k <- c(BFIT$sigma)/mu^3
    rho <- 1

    # return point estimates
    if(!complete)
    {
      CI <- mu
      CI[2] <- inverse.mean(mu,c(BFIT$sigma)/mu^2)  # 1/mu debiased to first order
      CI[3] <- mu*k
      CI[4] <- sqrt(CI[3])
      return(CI)
    }

    # print model selection table
    if(!is.na(IC))
    {
      BFIT0 <- meta.normal(s,VAR,debias=debias,VARS=FALSE) # no population variance

      dIC <- rbind(BFIT0[[IC]],BFIT[[IC]])
      rownames(dIC) <- c("Dirac-\u03B4","inverse-Gaussian")
      colnames(dIC) <- paste0("\u0394",IC)
      IND <- sort(c(dIC),index.return=TRUE)$ix
      dIC <- dIC[IND,,drop=FALSE]
      dIC <<- dIC # store to main function
      print(dIC - min(dIC))
    }
    else
    { IND <- 2 }

    CI <- t(array(c(0,0,Inf),c(3,4))) # initialize confidence intervals (0,Inf)

    IND <- IND[1] # best model
    if(IND==1) # exact chi^2 result
    {
      mu <- c(BFIT$mu)
      k <- 0
      CI.VAR <- BFIT$COV.mu
      dof <- 2*mu^2/CI.VAR[1] # sampling dof (not population)
      # mean area is the only area
      CI[1,] <- chisq.ci(mu,DOF=dof,level=level)
      CI[2,] <- 1/CI[1,] # not used
      CI[2,2] <- CI[2,2] * max(dof-debias,0)/dof
      CI.VAR[2] <- 2*CI[2,2]^2 / max(dof-debias*3,0)

      # CoV^2 (RVAR)
      CI[3,] <- c(0,0,Inf)
      CI.VAR[3] <- 0

      # CoV   (RSTD)
      CI[4,] <- c(0,0,Inf)
      CI.VAR[4] <- 0
    }
    else if(IND==2) # population variance (IG relations)
    {
      # k == VAR/mean^3
      GRAD <- rbind(c(1,0),c(-3*k/mu,1/mu^3)) # d(mu,k)/d(mean,var)
      COV <- diag(c(BFIT$COV.mu,BFIT$COV.sigma),2)
      COV <- GRAD %*% COV %*% t(GRAD) # (mu,k) COV

      diag(COV) <- nant(diag(COV),Inf)
      COV <- nant(COV,0)

      # delta method for population quantile uncertainty
      # population point estimates are IG
      par <- c(mu,k)
      CI.VAR <- COV[1,1]

      # approximate ratio of chi^2 to IG behavior
      DRATIO <<- DD.IG.ratio(par,CI.VAR[1],n=n)

      CI[1,] <- chisq.IG.ci(par[1],CI.VAR[1],w=DRATIO,level=level,precision=precision)

      CI[2,] <- 1/CI[1,] # numbers not used
      Vm2 <- c(BFIT$COV.mu)/mu^2 # keep fixed
      CI[2,2] <- inverse.mean(mu,Vm2) # 1/mu debiased to first order
      GRAD <- -1/mu^2/(1+debias*Vm2) # gradient evaluated analytically
      CI.VAR[2] <- GRAD^2 * c(BFIT$COV.mu)

      # CoV^2 (RVAR)
      # CI[5,2] <- par[1]*par[2]
      GRAD <- par[2:1]
      CI.VAR[3] <- GRAD %*% COV %*% GRAD
      CI[3,] <- chisq.ci(par[1]*par[2],VAR=CI.VAR[3],level=level)
      # CoV (RSTD)
      CI[4,] <- sqrt(CI[3,])
      CI.VAR[4] <- CI.VAR[3]/CI[3,2]/4
    } # end population variance

    rownames(CI) <- c("mean","inverse-mean","CoV\u00B2 (RVAR)","CoV  (RSTD)")
    colnames(CI) <- NAMES.CI

    R <- list(CI=CI,CI.VAR=CI.VAR,mu=mu,k=k,rho=1)
    return(R)
  } # end normal fit function

  if(method=="mle")
  { fit <- fit.mle }
  else if(method=='blue')
  { fit <- fit.blue }

  # MLE
  STUFF <- fit(s)
  CI <- STUFF$CI
  CI.VAR <- STUFF$CI.VAR
  mu <- STUFF$mu
  k <- STUFF$k
  rho <- STUFF$rho
  par <- c(mu,k) # HARDCODED for IG model

  # there is no point in bootstrapping if DOF==Inf (use exact MVU result instead)
  ## bias is in the opposite direction to give population variance
  ## and MLE is MVU
  if(boot && k) # parametric bootstrap
  {
    S1 <- S2 <- 0 # sums
    AVE <- VAR <- 0 # cumulants: mean, variance
    ENSEMBLE <- NULL

    INF <- c(0,Inf,Inf)
    INF <- rbind(INF,INF,INF,INF)
    rownames(INF) <- c("mean","inverse-mean","CoV\u00B2 (RVAR)","CoV  (RSTD)")
    colnames(INF) <- NAMES.CI

    iERROR <- Inf
    while(iERROR>=error)
    {
      N <- 0
      ERROR <- Inf
      pb <- utils::txtProgressBar(style=3)
      while(ERROR>=error || N<20)
      {
        # sample from prior, then sample from chi^2
        b <- statmod::rinvgauss(n,mu,1/k)
        b <- b * stats::rchisq(n,dof)/dof

        SAMPLE <- fit(b,par)
        ENSEMBLE <- cbind(ENSEMBLE,SAMPLE)
        N <- ncol(ENSEMBLE)

        # Will only happen for GIG distribution
        if(any(SAMPLE==Inf))
        {
          warning("Sampling distribution does not always resolve mean.")
          return(INF)
        }
        S1 <- S1 + SAMPLE
        S2 <- S2 + SAMPLE^2

        if(N>1)
        {
          AVE <- S1/N # mean
          VAR <- S2/N - AVE^2 # variance
          VAR <- N/(N-debias) * VAR # debias
          # max relative error
          ERROR <- max( sqrt(VAR/N)/AVE )
        }

        utils::setTxtProgressBar(pb,clamp(min(N/20,(error/ERROR)^2)))
      } # end single bootstrap loop

      # more numerically precise cumulant calculation
      AVE <- apply(ENSEMBLE,1,mean)
      CI.VAR <- apply(ENSEMBLE,1,stats::var)
      alpha <- 1-level
      Q1 <- apply(ENSEMBLE,1,function(x){stats::quantile(x,probs=alpha/2)})
      Q2 <- apply(ENSEMBLE,1,function(x){stats::quantile(x,probs=1-alpha/2)})
      if(debias) # first-order bias removed by this, while respecting boundary, and keeping CIs proportional
      { CI <- cbind(Q1,AVE,Q2) * (CI[,2]/AVE)^2 } # 1/mu is a little different because we already assigned true mu, but not sure how to leverage that
      else # keep CIs centered on biased estimate
      { CI <- cbind(Q1,CI[,2],Q2) }
      close(pb)

      # iterate debias
      # mu.pred <- AVE[2]
      # debias mu
      # debias k

      if(!iterate) { break }
    } # end bootstrap iteration
  }

  # # mean not estimated
  # if(!robust && method=="MLE")
  # {
  #   CI[2,1] <- 0 # low CI
  #   CI[2,4] <- 0 # DOF
  #   warning("Population mean not convergent. Consider robust=TRUE.")
  # }

  rownames(CI) <- c("mean","inverse-mean","CoV\u00B2 (RVAR)","CoV  (RSTD)")
  colnames(CI) <- NAMES.CI

  if(CI.FN=="beta")
  {
    for(i in c(1,3,4))
    { CI[i,1:3] <- 100*sqrt(beta.ci(CI[i,2],CI.VAR[i],level=level)) }
    CI[2,1:3] <- 1/(100*sqrt(beta.ci(1/CI[2,2],CI.VAR[2]/CI[2,2]^4,level=level)))
    CI.VAR <- CI.VAR/2^2/CI[,2]
  }

  R <- list(CI=CI,VAR=CI.VAR,dIC=dIC)
  return(R)
}


# shrinkage estimator
shrink.chisq <- function(s,dof,S,DOF,...)
{
  #
  #
  #
}


meta <- function(x,variable="area",level=0.95,level.UD=0.95,method="MLE",IC="AICc",boot=FALSE,error=0.01,debias=TRUE,verbose=FALSE,units=TRUE,plot=TRUE,sort=FALSE,mean=TRUE,col="black",...)
{
  method <- canonical.name(method)
  method <- match.arg(method,c("mle","blue"))

  variable <- canonical.name(variable)
  variable <- match.arg(variable,c("area","diffusion","speed","tauposition","tauvelocity","distance","periodicity","cyclicity"))

  meta.uni(x=x,variable=variable,level=level,level.UD=level.UD,IC=IC,boot=boot,error=error,debias=debias,method=method,verbose=verbose,units=units,plot=plot,sort=sort,mean=mean,col=col,...)
}


############
import.variable <- function(x,variable="area",level.UD=0.95,chi=FALSE)
{
  x <- name.list(x)
  N <- length(x)
  ID <- names(x) # may be null

  # TODO do range or diffusion, but not both
  # TODO do range or diffusion, but not both
  # TODO do range or diffusion, but not both

  AREA <- DOF <- array(0,N)
  for(i in 1:N)
  {
    CLASS <- class(x[[i]])
    if(length(ID) < i) { ID[i] <- attr(x[[i]],"info")$identity }

    if(CLASS=="ctmm")
    {
      if(variable=="area")
      {
        AREA[i] <- area.covm(x[[i]]$sigma)
        DOF[i] <- DOF.area(x[[i]]) * 2 # 2D
      }
      else if(variable=="diffusion")
      {
        STUFF <- diffusion(x[[i]],finish=FALSE)
        AREA[i] <- STUFF$D
        DOF[i] <- STUFF$DOF
      }
      else if(variable=="speed")
      {
        STUFF <- speed(x[[i]],units=FALSE)
        AREA[i] <- STUFF$CI[2]
        DOF[i] <- STUFF$DOF * 2 # 2D
      }
      else if(variable=="kinetic")
      {
        STUFF <- summary(x[[i]],units=FALSE)
        DOF[i] <- STUFF$DOF['speed'] * 2 # 2D
        if(DOF[i]>0)
        { AREA[i] <- STUFF$CI['speed (meters/second)','est']^2 }
        else
        { AREA[i] <- Inf }
      }
      else if(variable=="tauposition")
      {
        E <- x[[i]]$tau[1]
        if(is.null(E) || is.na(E))
        {
          AREA[i] <- 0
          DOF[i] <- 0
        }
        else
        {
          AREA[i] <- E
          P <- c("tau position","tau")
          P <- P[ P %in% rownames(x[[i]]$COV) ]
          DOF[i] <- 2*E^2/x[[i]]$COV[P,P]
        }
      }
      else if(variable=="tauvelocity")
      {
        E <- x[[i]]$tau[2]
        if(is.null(E) || is.na(E))
        {
          AREA[i] <- 0
          DOF[i] <- 0
        }
        else
        {
          AREA[i] <- E
          P <- c("tau velocity","tau")
          P <- P[ P %in% rownames(x[[i]]$COV) ]
          DOF[i] <- 2*E^2/x[[i]]$COV[P,P]
        }
      }
      else if(variable=="periodicity")
      {
        R <- periodic.variances(x[[i]])$R
        AREA[i] <- R$MLE
        DOF[i] <- nant( 2*R$MLE^2/R$VAR, 0 )
      }
      else if(variable=="cyclicity")
      {
        R <- periodic.variances(x[[i]])$V
        AREA[i] <- R$MLE
        DOF[i] <- 2*R$MLE^2/R$VAR
      }
    }
    else # UD or summary(UD) or speed()
    {
      if(CLASS=="overlap")
      {
        variable <- "distance"
        DOF[i] <- x[[i]]$DOF[1,2]
        AREA[i] <- -log(x[[i]]$CI[1,2,'est'])
      }
      if(CLASS=="distance")
      {
        variable <- "distance"
        DOF[i] <- x[[i]]$DOF[1,2]
        AREA[i] <- x[[i]]$CI[1,2,'est']
      }
      if(CLASS=="speed" || variable=="speed")
      {
        variable <- "speed"
        DOF[i] <- x[[i]]$DOF * 2 # 2D
        UNITS <- rownames(x[[i]]$CI) # "speed (units)"
        UNITS <- substr(UNITS,nchar("speed (")+1,nchar(UNITS)-1)
        AREA[i] <- x[[i]]$CI[2] %#% UNITS
      }
      else if(variable=="area")
      {
        if(CLASS=="UD") { x[[i]] <- summary(x[[i]],level.UD=level.UD,units=FALSE) }
        # now summary(UD) list
        DOF[i] <- x[[i]]$DOF['area'] * 2 # 2D
        # convert to SI units
        UNITS <- rownames(x[[i]]$CI) # "area (units)"
        UNITS <- substr(UNITS,nchar("area (")+1,nchar(UNITS)-1)
        AREA[i] <- x[[i]]$CI[2] %#% UNITS
      }
    }
  }

  # level.UD coverage (e.g., 95% home ranges rather than straight variance)
  if(variable=="area" && CLASS=="ctmm")
  { AREA <- -2*log(1-level.UD)*pi * AREA }
  else if(variable=="speed" && !chi) # chi DOF changes when chi is approximated as chi^2
  {
    VAR <- chi.var(DOF) # modulo E[X]^2
    DOF <- 2/VAR
  }

  R <- list(ID=ID,AREA=AREA,DOF=DOF,variable=variable)
  return(R)
}


# wrapper: meta-analysis of CTMM areas
# TODO range=FALSE ???
meta.uni <- function(x,variable="area",level=0.95,level.UD=0.95,level.pop=0.95,method="MLE",IC="AICc",boot=FALSE,error=0.01,debias=TRUE,verbose=FALSE,units=TRUE,plot=TRUE,sort=FALSE,mean=TRUE,col="black",...)
{
  N <- length(x)

  # N group comparisons (list of lists that are not summaries)
  SUBPOP <- class(x)[1]=='list' && class(x[[1]])[1]=='list' && !( length(names(x[[1]]))==2 && all(names(x[[1]])==c('DOF','CI')) )

  if(SUBPOP) { CLASS <- class(x[[1]][[1]])[1] }
  else { CLASS <- class(x[[1]])[1] }

  # fix variable argument if necessary
  if(CLASS %in% c("UD","area"))
  { variable <- "area" }
  else if(CLASS=="speed")
  { variable <- "speed" }
  else if(CLASS=="overlap")
  { variable <- "overlap" }
  else if(CLASS=="distance")
  { variable <- "distance" }

  if(variable %in% c('periodicity','cyclicity'))
  { CI.FN <- "beta" }
  else
  { CI.FN <- "chisq" }

  if(SUBPOP)
  {
    ID <- names(x)
    ID[N+1] <- "mean"

    # analyze all N groups separately
    RESULTS <- AREA <- DOF <- list()
    for(i in 1:N)
    {
      STUFF <- import.variable(x[[i]],variable,level.UD=level.UD)
      AREA[[i]] <- STUFF$AREA
      DOF[[i]] <- STUFF$DOF

      message(paste("* Sub-population",ID[i]))
      RESULTS[[i]] <- meta.chisq(AREA[[i]],DOF[[i]],level=level,level.pop=level.pop,IC=IC,boot=boot,error=error,debias=debias,method=method,verbose=TRUE,units=FALSE,CI.FN=CI.FN)
    }
    message("* Joint population")
    RESULTS[[N+1]] <- meta.chisq(unlist(AREA),unlist(DOF),level=level,level.pop=level.pop,IC=IC,boot=boot,error=error,debias=debias,method=method,verbose=TRUE,units=FALSE,CI.FN=CI.FN)
    names(RESULTS) <- ID

    message("* Joint population versus sub-populations (best models)")
    dIC <- sum( sapply(RESULTS[1:N],function(R){R$dIC[1,]}) )
    dIC[2] <- RESULTS[[N+1]]$dIC[1,]
    dIC <- cbind(dIC)
    rownames(dIC) <- c("Sub-population","Joint population")
    colnames(dIC) <- paste0("\u0394",IC)
    IND <- sort(c(dIC),index.return=TRUE)$ix
    dIC <- dIC[IND,,drop=FALSE]
    print(dIC - min(dIC))

    # forest plot object
    PLOT <- sapply(RESULTS,function(R){R$CI[1,1:3]}) # [3,N+1]
  }
  else
  {
    STUFF <- import.variable(x,variable=variable,level.UD=level.UD)
    AREA <- STUFF$AREA
    DOF <- STUFF$DOF
    ID <- STUFF$ID

    # inverse-Gaussian population distribution
    CI <- meta.chisq(AREA,DOF,level=level,level.pop=level.pop,IC=IC,boot=boot,error=error,debias=debias,method=method,CI.FN=CI.FN)
    CI.VAR <- CI$VAR
    CI <- CI$CI

    # basic forest plot
    if(CI.FN=="chisq")
    { PLOT <- sapply(1:(N+1),function(i){chisq.ci(AREA[i],DOF=DOF[i],level=level)}) } # [3,N+1]
    else if(CI.FN=="beta")
    { PLOT <- sapply(1:(N+1),function(i){100*sqrt(beta.ci(AREA[i],2*AREA[i]^2/DOF[i],level=level))}) }

    ID[N+1] <- "mean"
    PLOT[,N+1] <- CI[1,1:3] # overwrite chi^2 CI with better
  }

  if(variable=="tauposition")
  {
    VAR.UNITS <- "time"
    VAR.NAME <- "Position timescale"
  }
  else if(variable=="tauvelocity")
  {
    VAR.UNITS <- "time"
    VAR.NAME <- "Velocity timescale"
  }
  else if(variable %in% c("overlap","distance"))
  {
    VAR.UNITS <- "dissimilarity"
    VAR.NAME <- "Distance"
  }
  else
  {
    VAR.UNITS <- variable
    VAR.NAME <- capitalize(variable)
  }

  if(plot)
  {
    PLOT <- PLOT[,1:(N+mean)] # drop mean if FALSE
    IND <- (N+mean):1

    M <- length(col)
    col <- array(col,N+mean)
    if(M<N+1 && mean) { col[N+1] <- "black" } # default if mean not specified

    UNITS <- unit(PLOT[,N+1],VAR.UNITS,SI=!units)
    PLOT <- PLOT/UNITS$scale

    xlab <- VAR.NAME
    if(length(UNITS$name)) { xlab <- paste0(xlab," (",UNITS$name,")") }
    if(variable=="area") { xlab <- paste0(100*level.UD,"% ",xlab) }
    # base layer plot
    RANGE <- range(PLOT[PLOT<Inf])
    RANGE[2] <- min(RANGE[2],10*PLOT[3,N+1])
    plot(RANGE,c(1,N+mean),col=grDevices::rgb(1,1,1,0),xlab=xlab,ylab=NA,yaxt="n",...)

    # 2nd attempt to fix long labels # still not working, but better than nothing
    CEX.AXIS <- graphics::par("cex.axis")
    CEX.NEW <- CEX.AXIS
    # fix too wide
    MAX <- max(graphics::strwidth(ID,units="inches"))
    OFFSET <- 1.5*graphics::strheight("A",units="inches") # 3x default tick length is still not enough? (2x should be about perfect?)
    CEX.NEW[2] <- (graphics::par("mai")[2]-OFFSET)/(graphics::par("cex")*MAX)
    # fix too tall
    CEX.NEW[3] <- graphics::par('pin')[1] / sum(graphics::strheight(ID,units="inches"))
    # fix too tall & too wide
    CEX.NEW <- min(CEX.NEW)
    # zero is not valid
    CEX.NEW <- max(CEX.NEW,0.01)
    graphics::par(cex.axis=CEX.NEW)
    # will reset when done with axis()

    if(sort)
    {
      SORT <- sort(PLOT[2,1:N],decreasing=sort,index.return=TRUE)$ix
      ID[1:N] <- ID[SORT]
      PLOT[,1:N] <- PLOT[,SORT]
      if(length(col)>1) { col[1:N] <- col[SORT] }
    }

    # colored axes
    for(i in 1:N) { graphics::axis(2,at=IND[i],labels=ID[i],las=2,col.axis=col[i]) }
    if(mean)
    {
      graphics::axis(2,at=1,labels=ID[N+1],las=2,font=2,col.axis=col[N+1])
      # mean vertical line
      graphics::abline(v=PLOT[2,N+1],col=malpha(col[N+1],1/2))
    }
    # error bars
    graphics::points(PLOT[2,],IND,pch=16,col=col)
    suppressWarnings( graphics::arrows(PLOT[1,],IND,PLOT[3,],IND,length=0.05,angle=90,code=3,col=col) )

    # reset cex.axis
    graphics::par(cex.axis=CEX.AXIS)
  }

  if(SUBPOP) # ratios CIs
  {
    RESULTS <- RESULTS[1:N]
    ID <- ID[1:N]

    NUM <- sapply(RESULTS,function(R){R$CI[1,2]})
    NVAR <- sapply(RESULTS,function(R){R$VAR[1]})
    DEN <- sapply(RESULTS,function(R){R$CI[2,2]})
    DVAR <- sapply(RESULTS,function(R){R$VAR[2]})

    dof <- pmax(2*NUM^2/NVAR,1)

    CI <- array(1,c(N,N,3))
    PV <- array(1,c(N,N))
    dimnames(CI) <- list(paste0(ID,"/"),paste0("/",ID),NAMES.CI)
    dimnames(PV) <- list(paste0(ID,"/"),paste0("/",ID))
    for(i in 1:N)
    {
      for(j in (1:N)[-i]) # diagonal == 1/1
      {
        CI[i,j,] <- F.CI(NUM[i],NVAR[i],DEN[j],DVAR[j],level=level)
        PV[i,j] <- stats::pf(NUM[i]/NUM[j],dof[i],dof[j],lower.tail=NUM[i]<NUM[j])
      }
    }

    if(verbose)
    {
      UNITS <- sapply(RESULTS,function(R){R$CI[1,]})
      UNITS <- unit(UNITS,variable,SI=!units,concise=TRUE)

      for(i in 1:N)
      {
        RESULTS[[i]] <- RESULTS[[i]]$CI[c(1,3,4),]
        RESULTS[[i]][1,] <- RESULTS[[i]][1,]/UNITS$scale
        if(length(UNITS$name)) { rownames(RESULTS[[i]])[1] <- paste0(rownames(RESULTS[[i]])[1]," (",UNITS$name,")") }
      }

      RESULTS[[N+1]] <- PV
      names(RESULTS)[N+1] <- "p-value"

      RESULTS[[N+2]] <- CI
      names(RESULTS)[N+2] <- "mean ratio"

      CI <- RESULTS
    }
  }
  else # population levels CIs
  {
    CI <- CI[c(1,3,4),] # mean and CoV^2

    UNITS <- unit(CI[1,],VAR.UNITS,SI=!units,concise=TRUE)
    CI[1,] <- CI[1,]/UNITS$scale

    if(length(UNITS$name)) { rownames(CI)[1] <- paste0(rownames(CI)[1]," (",UNITS$name,")") }
    #rownames(CI)[2] <- "CoV\u00B2 (RVAR)"
    #rownames(CI)[3] <- "CoV  (RSTD)"
  }

  return(CI)
}
