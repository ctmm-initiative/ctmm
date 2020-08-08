NAMES.POP <- c("pop-low ","pop-ave ","pop-high")

# meta-analysis of chi^2 random variables with inverse-Gaussian prior
meta.chisq <- function(s,dof,level=0.95,level.pop=0.95,IC="AICc",method='exact',boot=FALSE,iterate=FALSE,error=0.01,debias=TRUE,precision=1/2,...)
{
  # tol <- .Machine$double.eps^precision
  n <- length(s)
  # for model ICs
  dIC <- Inf

  # CI that are chi^2 when VAR[M]>>VAR.pop
  ci.chisq.GIG <- function(M,VAR,mu,k=0,rho=1)
  {
    # population moments
    # k == 1/lambda
    theta <- 1/(k*mu)
    M.pop <- mu * (besselK(theta,rho/2-1,expon.scaled=TRUE)/besselK(theta,rho/2,expon.scaled=TRUE))
    VAR.pop <- mu^2 * (besselK(theta,rho/2-2,expon.scaled=TRUE)/besselK(theta,rho/2,expon.scaled=TRUE)) - M.pop^2

    # chi^2 VAR.pop==0
    # IG when VAR==VAR.pop/n
    w <- VAR.pop/VAR # (0,n) : (chi^2,IG)
    w <- clamp(w,0,n) / n # (0,1) : (chi^2,IG)
    w <- c(1-w,w)

    # approximation that bridges the two with correct mean and variance
    CI <- w[1]*chisq.ci(M,VAR=VAR,level=level) + w[2]*IG.CI(M,VAR=VAR,level=level,precision=precision)
    # TODO make the second term some kind of GIG?
    # TODO GIG stuff
    return(CI)
  }

  # Bessel's correction
  BC <- n/(n-1)

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
    nloglike <- function(par,zero=0)
    {
      # subset for CV
      # s <- s[SUB]
      # dof <- dof[SUB]
      # n <- length(s)

      mu <- par[1] # mean
      k <- ifelse(length(par)>=2,par[2],0) # IG var/mu^3
      theta <- 1/(k*mu)
      rho <- ifelse(length(par)>=3,par[3],1) # GIG shape

      if(k<0 || k==Inf || mu<=0 || mu==Inf) { return(Inf) }

      zero <- zero + sum( dof/2*log(dof/2) +(dof/2-1)*log(s) - lgamma(dof/2) )

      alpha <- dof*s/mu
      beta <- alpha/theta
      beta <- 1/2*log1p(beta) # log(sqrt(1+dof*s*k))
      -sum( -dof/2*log(mu) - (dof+rho)/2*beta + lKK(rho,dof,theta,alpha) + zero/n )
    }

    # n.min <- 1-2 (AIC-AICc)
    if(complete)
    {
      NAME[1] <- "Dirac-\u03B4"
      mu <- sum(dof*s)/sum(dof) # exact solution
      PAR[1,1] <- mu
      LL[1] <- -nloglike(mu)
    }

    # n.min <- 2-4 (AIC-AICc)
    NAME[2] <- "inverse-Gaussian"
    if(complete && n>=2) # par==NULL
    {
      mu <- mean(s)
      k <- mean(1/s) - 1/mu
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

    ###################
    # return point estimates for bootstrap - no model selection nor uncertainty CIs needed
    if(!complete)
    {
      CI <- rep(NA_real_,3)
      if(length(par)==2)
      {
        mu <- PAR[2,1]
        k <- PAR[2,2]
        if(debias) # kappa for debiased standard deviation
        {
          k <- BC * k
          # DEBUG <<- list(mu=mu,k=k,k.std=k.std,IG.CI=IG.CI)
          k <- k.std(PAR[2,])
        }
        CI <- IG.CI(mu,k=k,level=level.pop,precision=precision) }
      # else if(length(par)==3)
      # { CI <- GIG.CI(eta,theta,rho,level=level.pop,precision=precision) }
      names(CI) <- NAMES.POP
      return(CI)
    }

    if(!is.na(IC))
    {
      ICS <- c("AIC","AICc","BIC")
      dIC <- matrix(0,MODELS,length(ICS))
      colnames(dIC) <- ICS

      k <- c(1,1,1)[1:length(LL)] # mean parameters
      nu <- c(0,1,2)[1:length(LL)] # variance parameters
      K <- k + nu # total parameters
      dIC[,"AIC"]  <- 2*K - 2*LL
      dIC[,"AICc"] <- 2*K*ifelse(debias,n-k,n)/pmax(n-K-nu,0) - 2*LL
      dIC[,"BIC"]  <- log(n)*K - 2*LL

      dIC <- cbind(dIC[,IC])
      rownames(dIC) <- NAME
      colnames(dIC) <- paste0("\u0394",IC)
      IND <- sort(c(dIC),index.return=TRUE)$ix
      dIC <- dIC[IND,,drop=FALSE]
      dIC <<- dIC # store to main function
      print(dIC - min(dIC))

      # some old code that would need to be re-written
      # if(IC=="LOOCV")
      # {
      #   LL0 <- 0
      #   for(i in 1:n)
      #   {
      #     # leave-one-out MLE
      #     SUB <- -i
      #     S <- sum( (dof*s)[SUB] ) / sum( dof[SUB] )
      #     # (CV) log-likelihood of left out
      #     SUB <- i
      #     LL0 <- LL0 - nloglike(c(S,0))
      #   }
      # SUB <- 1:n
      # LLN <- 0
      # for(i in 1:n)
      # {
      #   # leave-one-out MLE
      #   SUB <- -i
      #   par.i <- optimizer(par,nloglike,lower=0,upper=Inf)$par
      #   # CV log-likelihood of left out
      #   SUB <- i
      #   LLN <- LLN - nloglike(par.i)
      # }
      # SUB <- 1:n
      # dIC <- -2*LLN + 2*LL0
    }
    else # unknown DOF (no model selection)
    { IND <- 2  }

    ##################
    # propagate uncertainty
    CI <- t(array(c(0,0,Inf),c(3,3))) # initialize confidence intervals (0,Inf)

    IND <- IND[1] # best model
    PAR <- PAR[IND,] # best model parameter estimates
    par <- PAR[1:IND] # canonical subset of parameters

    if(IND==1) # Dirac delta -- DOF==Inf
    {
      dof <- sum(dof) # precision of point estimate
      # mean area is the only area
      CI[1,] <- CI[2,] <- CI[3,] <- chisq.ci(PAR[1],DOF=dof,level=level)
      CI.DOF <- array(dof,3)
    }
    else if(IND==2) # inverse-Gaussian model
    {
      STUFF <- genD(par,nloglike,lower=c(0,0),order=2)
      COV <- cov.loglike(STUFF$hessian,STUFF$gradient)

      if(debias) # Bessel's correction to point estimate and COV
      {
        par[2] <- BC * par[2]
        COV[2,] <- BC * COV[2,] # COV[BC*kappa] = BC COV[kappa] BC
        COV[,2] <- BC * COV[,2]
        COV[1,] <- sqrt(BC) * COV[1,] # VAR[mu] = VAR[x]/n
        COV[,1] <- sqrt(BC) * COV[,1]
        PAR[1:IND] <- par
      }

      Q <- function(par) { IG.CI(par[1],k=ifelse(debias,k.std(par),par[2]),level=level.pop,precision=precision) }
      GRAD <- genD(par,Q,lower=c(0,0),order=1)$gradient
      VAR <- diag(GRAD %*% COV %*% t(GRAD))

      # if population VAR is small, then CIs are approximately chi^2
      # if population VAR is large, then CIs are approximately IG
      # approximation that bridges the two with correct mean and variance
      CI[,2] <- Q(par)
      for(i in 1:3){ CI[i,] <- ci.chisq.GIG(CI[i,2],VAR=VAR[i],mu=PAR[1],k=PAR[2],rho=PAR[3]) }

      CI.DOF <- 2*CI[,2]^2/VAR # DOFs for F-test
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

    rownames(CI) <- NAMES.POP
    colnames(CI) <- NAMES.CI

    R <- list(CI=CI,CI.DOF=CI.DOF,mu=PAR[1],k=PAR[2],rho=PAR[3])
    return(R)
  } # end IG/GIG fitting function

  ###############
  fit.blue <- function(s,par=NULL)
  {
    complete <- is.null(par)
    VAR <- 2*s^2/dof
    STUFF.V <- meta.normal(s,VAR,debias=debias)
    # list(mu=mu,sigma=sigma,COV.mu=COV.mu,COV.sigma=COV.sigma,loglike=loglike,AIC=AIC,AICc=AICc,BIC=BIC)

    # IG parameters
    mu <- STUFF.V$mu
    k <- c(STUFF.V$sigma)/mu^3
    rho <- 1

    # return point estimates
    if(!complete)
    {
      CI <- IG.CI(mu,VAR=c(STUFF.V$sigma),level=level.pop)
      names(CI) <- NAMES.POP
      return(CI)
    }

    # print model selection table
    if(!is.na(IC))
    {
      STUFF.Z <- meta.normal(s,VAR,debias=debias,VARS=FALSE) # no population variance

      dIC <- rbind(STUFF.Z[[IC]],STUFF.V[[IC]])
      rownames(dIC) <- c("Dirac-\u03B4","inverse-Gaussian")
      colnames(dIC) <- paste0("\u0394",IC)
      IND <- sort(c(dIC),index.return=TRUE)$ix
      dIC <- dIC[IND,,drop=FALSE]
      dIC <<- dIC # store to main function
      print(dIC - min(dIC))
    }
    else
    { IND <- 2 }

    CI <- t(array(c(0,0,Inf),c(3,3))) # initialize confidence intervals (0,Inf)

    IND <- IND[1] # best model
    if(IND==1) # INF complete
    {
      mu <- c(STUFF.Z$mu)
      k <- 0
      VAR <- c(STUFF.Z$COV.mu)
      dof <- 2*mu^2/VAR # sampling dof (not population)
      # mean area is the only area
      CI[1,] <- CI[2,] <- CI[3,] <- chisq.ci(mu,DOF=dof,level=level)
      CI.DOF <- array(dof,3)
    }
    else if(IND==2) # population variance (IG relations)
    {
      # k == VAR/mean^3
      GRAD <- rbind(c(1,0),c(-3*k/mu,1/mu^3)) # d(mu,k)/d(mean,var)
      COV <- diag(c(STUFF.V$COV.mu,STUFF.V$COV.sigma),2)
      COV <- GRAD %*% COV %*% t(GRAD) # (mu,k) COV

      # delta method for population quantile uncertainty
      # population point estimates are IG
      par <- c(mu,k)
      Q <- function(par) { IG.CI(par[1],k=par[2],level=level.pop,precision=precision) }
      GRAD <- genD(par,Q,lower=c(0,0),order=1)$gradient # [Q,(mu,k)]
      VAR <- diag(GRAD %*% COV %*% t(GRAD)) # [Q,Q]
      # point estimates
      CI[,2] <- Q(par)
      # population low,point,high CIs
      for(i in 1:3) { CI[i,] <- ci.chisq.GIG(CI[i,2],VAR[i],mu=mu,k=k) } # point estimate already robust from above
      CI.DOF <- 2*CI[,2]^2/VAR
    } # end population variance

    rownames(CI) <- NAMES.POP
    colnames(CI) <- NAMES.CI

    R <- list(CI=CI,CI.DOF=CI.DOF,mu=mu,k=k,rho=1)
    return(R)
  } # end normal fit function

  if(method=="mle")
  { fit <- fit.mle }
  else if(method=='blue')
  { fit <- fit.blue }

  # MLE
  STUFF <- fit(s)
  CI <- STUFF$CI
  mu <- STUFF$mu
  k <- STUFF$k
  rho <- STUFF$rho
  par <- c(mu,k) # HARDCODED for IG model

  # there is no point in bootstrapping if DOF==Inf
  ## bias is in the opposite direction to give population variance
  ## and MLE is MVU

  if(!boot || k==0) # MLE
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
          VAR <- N/(N-1) * VAR # debias
          # max relative error
          ERROR <- max( sqrt(VAR/N)/AVE )
        }

        utils::setTxtProgressBar(pb,clamp(min(N/20,(error/ERROR)^2)))
      } # end single bootstrap loop

      # more numerically precise cumulant calculation
      AVE <- apply(ENSEMBLE,1,mean)
      VAR <- apply(ENSEMBLE,1,stats::var)
      alpha <- 1-level
      Q1 <- apply(ENSEMBLE,1,function(x){stats::quantile(x,probs=alpha/2)})
      Q2 <- apply(ENSEMBLE,1,function(x){stats::quantile(x,probs=1-alpha/2)})
      DOF <- 2*AVE^2/VAR
      if(debias) # first-order bias removed by this, while respecting boundary, and keeping CIs proportional
      { CI <- cbind(Q1,AVE,Q2) * (CI[,2]/AVE)^2 }
      else # keep CIs centered on biased estimate
      { CI <- t(sapply(1:3,function(i){ci.chisq.GIG(CI[i,2],VAR=2*CI[i,2]^2/DOF[i],mu=mu,k=k)})) }
      rownames(CI) <- NAMES.POP
      colnames(CI) <- NAMES.CI
      CI <- cbind(CI,DOF=DOF)
      close(pb)

      # iterate debias
      # mu.pred <- AVE[2]
      # debias eta
      # debias eta/sqrt(theta)

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

  R <- list(CI=CI,dIC=dIC)
  return(R)
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


meta <- function(x,level=0.95,level.UD=0.95,level.pop=0.95,method="MLE",IC="AICc",boot=FALSE,error=0.01,debias=TRUE,units=TRUE,plot=TRUE,...)
{
  method <- tolower(method)
  method <- match.arg(method,c("mle","blue"))

  meta.area(x=x,level=level,level.UD=level.UD,level.pop=level.pop,IC=IC,boot=boot,error=error,debias=debias,method=method,units=units,plot=plot,...)
}


############
import.area <- function(x,level.UD=0.95)
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

  R <- list(ID=ID,AREA=AREA,DOF=DOF)
  return(R)
}


# wrapper: meta-analysis of CTMM areas
# TODO range=FALSE ???
meta.area <- function(x,level=0.95,level.UD=0.95,level.pop=0.95,method="MLE",IC="AICc",boot=FALSE,error=0.01,debias=TRUE,units=TRUE,plot=TRUE,...)
{
  N <- length(x)

  # N group comparisons (list of lists that are not summaries)
  SUBPOP <- class(x)=='list' && class(x[[1]])=='list' && !all(names(x[[1]])==c('DOF','CI'))
  if(SUBPOP)
  {
    ID <- names(x)
    ID[N+1] <- "Joint population"

    # analyze all N groups separately
    RESULTS <- AREA <- DOF <- ID <- list()
    for(i in 1:N)
    {
      STUFF <- import.area(x[[i]],level.UD=level.UD)
      AREA[[i]] <- STUFF$AREA
      DOF[[i]] <- STUFF$DOF

      print(paste("Sub-population",ID[i]))
      RESULTS[[i]] <- meta.chisq(AREA[[i]],DOF[[i]],level=level,level.pop=level.pop,IC=IC,boot=boot,error=error,debias=debias,method=method,verbose=TRUE)
    }
    print("Joint population")
    RESULTS[[N+1]] <- meta.chisq(unlist(AREA),unlist(DOF),level=level,level.pop=level.pop,IC=IC,boot=boot,error=error,debias=debias,method=method,verbose=TRUE)

    print("Joint population versus sub-populations")
    dIC <- sum( sapply(RESULTS[[1:N]],function(R){R$dIC[1,]}) )
    dIC[2] <- RESULTS[[N+1]]$dIC[1,]
    dIC <- cbind(dIC)
    colnames(dIC) <- c("sub-population","joint population")
    rownames(dIC) <- paste0("\u0394",IC)
    IND <- sort(c(dIC),index.return=TRUE)$ix
    dIC <- dIC[IND,,drop=FALSE]
    print(dIC - min(dIC))

    # forest plot object
    PLOT <- sapply(RESULTS,function(R){R$CI[2,1:3]}) # [3,N+1]

    # effect size matrix return
    #
    #
  }
  else
  {
    STUFF <- import.area(x,level.UD=level.UD)
    AREA <- STUFF$AREA
    DOF <- STUFF$DOF
    ID <- STUFF$ID

    # inverse-chi^2 population distribution
    CI <- meta.chisq(AREA,DOF,level=level,level.pop=level.pop,IC=IC,boot=boot,error=error,debias=debias,method=method)$CI

    # basic forest plot
    PLOT <- sapply(1:(N+1),function(i){chisq.ci(AREA[i],DOF=DOF[i],level=level)}) # [3,N+1]

    ID[N+1] <- "mean"
    PLOT[,N+1] <- CI[2,1:3] # overwrite chi^2 CI with better
  }

  if(plot)
  {
    IND <- (N+1):1

    UNITS <- unit(PLOT,"area",SI=!units)
    PLOT <- PLOT/UNITS$scale

    xlab <- paste0(100*level.UD,"% Area (",UNITS$name,")")
    plot(range(PLOT),c(1,N+1),col=grDevices::rgb(1,1,1,0),xlab=xlab,ylab=NA,yaxt="n",...)
    graphics::axis(2,at=IND[1:N],labels=ID[1:N],las=2)
    graphics::axis(2,at=1,labels=ID[N+1],las=2,font=2)
    graphics::abline(v=PLOT[2,N+1],col=grDevices::rgb(0.5,0.5,0.5,0.5))
    graphics::points(PLOT[2,],IND,pch=16)
    suppressWarnings( graphics::arrows(PLOT[1,],IND,PLOT[3,],IND,length=0.05,angle=90,code=3) )
  }

  if(SUBPOP)
  {
    EST <- sapply(RESULTS[1:N],function(R){R$CI[2,2]})
    DOF <- sapply(RESULTS[1:N],function(R){R$CI[2,4]})
    CI <- outer(EST,EST,'/')
    rownames(CI) <- ID[1:N]
    colnames(CI) <- ID[1:N]
  }
  else
  {
    CI <- CI[,1:3]
    UNITS <- unit(CI,"area",SI=!units)
    CI <- CI/UNITS$scale
    rownames(CI) <- paste0(rownames(CI)," area (",UNITS$name,")")
  }

  return(CI)
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
