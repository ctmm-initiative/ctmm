cluster <- function(x,level=0.95,level.UD=0.95,debias=TRUE,IC="BIC",units=TRUE,plot=TRUE,sort=FALSE,...)
{
  cluster.area(x=x,level=level,level.UD=level.UD,IC=IC,debias=debias,units=units,plot=plot,sort=sort,...)
}


# wrapper: meta-analysis of CTMM areas
# TODO range=FALSE ???
cluster.area <- function(x,level=0.95,level.UD=0.95,IC="BIC",debias=TRUE,units=TRUE,plot=TRUE,sort=FALSE,...)
{
  N <- length(x)
  ID <- names(x)
  STUFF <- import.variable(x,level.UD=level.UD)
  AREA <- STUFF$AREA
  DOF <- STUFF$DOF

  # inverse-chi^2 population distribution
  STUFF <- cluster.chisq(AREA,DOF,level=level,IC=IC,debias=debias)
  CI <- STUFF$CI # c("mu1","CoV1","mu2","CoV2","P1","P2","mu2/mu1")
  P <- STUFF$P # posterior probability of falling into first class
  names(P) <- ID

  AREA <- sapply(1:N,function(i){chisq.ci(AREA[i],DOF=DOF[i],level=level)})
  AREA <- t(AREA)

  # colors # black=normative # red=extreme
  COL <- grDevices::rgb(1-P,0,0)

  if(plot)
  {
    # first population at bottom
    ID <- c(ID,expression(mu[1]))
    AREA <- rbind(AREA,CI[1,])
    COL <- c(COL,grDevices::rgb(0,0,0))

    # second population at top
    ID <- c(expression(mu[2]),ID)
    AREA <- rbind(CI[3,],AREA)
    COL <- c(grDevices::rgb(1,0,0),COL)

    if(sort)
    {
      SORT <- sort(AREA[1+1:N,2],decreasing=sort,index.return=TRUE)$ix
      ID[1+1:N] <- ID[1+SORT]
      AREA[1+1:N,] <- AREA[1+SORT,]
      COL[1+1:N] <- COL[1+SORT]
    }

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

    # basic forest plot
    IND <- (N+2):1

    UNITS <- unit(AREA,"area",SI=!units)
    AREA <- AREA/UNITS$scale

    xlab <- paste0(100*level.UD,"% Area (",UNITS$name,")")
    # base layer plot
    RANGE <- range(AREA)
    RANGE[2] <- min(RANGE[2],10*max(AREA[1,3],AREA[N+2,3]))
    plot(RANGE,c(1,N+2),col=grDevices::rgb(1,1,1,0),xlab=xlab,ylab=NA,yaxt="n",...)

    for(i in 2:(N+1)) { graphics::axis(2,at=IND[i],labels=ID[i],las=2,col.axis=COL[i],lwd=0) }
    graphics::axis(2,at=IND[N+2],labels=ID[N+2],las=2,font=2,col.axis=COL[N+2],lwd=0)
    graphics::axis(2,at=IND[1],labels=ID[1],las=2,font=2,col.axis=COL[1],lwd=0)
    graphics::axis(2,at=IND,labels=FALSE)
    graphics::abline(v=AREA[1,2],col=grDevices::rgb(1,0.0,0.0,0.5))
    graphics::abline(v=AREA[N+2,2],col=grDevices::rgb(0.0,0.0,0.0,0.5))
    graphics::points(AREA[,2],IND,pch=16,col=COL)
    suppressWarnings( graphics::arrows(AREA[,1],IND,AREA[,3],IND,length=0.05,angle=90,code=3,col=COL) )

    # reset cex.axis
    graphics::par(cex.axis=CEX.AXIS)
  }

  UNITS <- unit(CI[c(1,3),],"area",SI=!units,concise=TRUE)
  CI[c(1,3),] <- CI[c(1,3),]/UNITS$scale
  rownames(CI)[c(1,3)] <- paste0(rownames(CI)[c(1,3)]," (",UNITS$name,")")

  return(list("P"=P,CI=CI))
}


# cluster-analysis of chi^2 random variables with mixture of 2 inverse-chi^2 priors
cluster.chisq <- function(s,dof,level=0.95,IC="BIC",debias=TRUE,precision=1/2,...)
{
  # discard null estimates
  ZERO <- dof<=.Machine$double.eps
  s <- s[!ZERO]
  dof <- dof[!ZERO]

  #tol <- .Machine$double.eps^precision
  n <- length(s)
  SORT <- sort(s,index.return=TRUE)$ix

  # MLE fit and CI return
  SUB <- 1:n # relevant subset

  # negative log-likelihood
  # S mean
  # k = VAR/S^3 inverse-Gaussian parameter
  nloglike <- function(w1=1,S1,K1=0,w2=1-w1,S2=S1,K2=K1,zero=0)
  {
    if(K1<0 || K1==Inf || S1<=0 || S1==Inf || K2<0 || K2==Inf || S2<=0 || S2==Inf) { return(Inf) }

    if(w1>0 && (K1<0 || S1==0 || S1==Inf)) { return(Inf) }
    if(w2>0 && (K2<0 || S2==0 || S2==Inf)) { return(Inf) }

    # subset to relevant set
    s <- s[SUB]
    dof <- dof[SUB]
    n <- length(s)

    rho <- 1

    # common factors # common to DD and IG population distributions
    zero <- zero + sum( dof/2*log(dof/2) + (dof/2-1)*log(s) - lgamma(dof/2) )

    if(K1>.Machine$double.eps) # inverse-Gaussian
    {
      theta <- 1/(K1*S1)
      alpha <- dof*s/S1
      beta <- alpha/theta
      beta <- 1/2*log1p(beta) # log(sqrt(1+dof*s*k))
      loglike1 <- -dof/2*log(S1) - (dof+rho)/2*beta + lKK(rho,dof,theta,alpha)
    }
    else # Dirac delta
    { loglike1 <- -dof/2*(log(S1)+s/S1) }

    # no mixture
    if(w1==1)
    {
      loglike1 <- loglike1 + zero/n
      loglike1 <- sum(loglike1)
      return(-loglike1)
    }

    if(K2>.Machine$double.eps) # IG
    {
      theta <- 1/(K2*S2)
      alpha <- dof*s/S2
      beta <- alpha/theta
      beta <- 1/2*log1p(beta) # log(sqrt(1+dof*s*k))
      loglike2 <- -dof/2*log(S2) - (dof+rho)/2*beta + lKK(rho,dof,theta,alpha)
    }
    else # DD
    { loglike2 <- -dof/2*(log(S2)+s/S2) }

    logw1 <- log(w1)
    logw2 <- log(w2)

    # sort for numerical stability in log-likelihood representation (which is required to prevent underflow)
    w.min <- w.max <- loglike.min <- loglike.max <- logw.min <- logw.max <- array(0,n)
    for(i in 1:n)
    {
      if(loglike1[i]+logw1 < loglike2[i]+logw2)
      {
        w.min[i] <- w1
        w.max[i] <- w2
        loglike.min[i] <- loglike1[i]
        loglike.max[i] <- loglike2[i]
        logw.min[i] <- logw1
        logw.max[i] <- logw2
      }
      else
      {
        w.min[i] <- w2
        w.max[i] <- w1
        loglike.min[i] <- loglike2[i]
        loglike.max[i] <- loglike1[i]
        logw.min[i] <- logw2
        logw.max[i] <- logw1
      }
    }

    # representation insensitive to w.min~0 & loglike.min<<loglike.max
    loglike <- logw.max + loglike.max + log(1+nant((w.min/w.max)*exp(loglike.min-loglike.max),0)) + zero/n
    -sum(loglike)
  }

  # lists of stuff organized by models
  L <- list()
  L$dirac <- list()       # 1
  L$gauss <- list()       # 2 # inverse Gaussian
  L$gauss.gauss <- list() #   # inverse Gaussian mixture with equal CoV

  L$dirac$K <- 1 # number of parameters
  L$dirac$COST <- function(par) { nloglike(S1=par) } # negative log-likelihood
  L$dirac$FIT <- function(par=NULL) # fit MLE
  {
    par <- sum(dof[SUB]*s[SUB])/sum(dof[SUB])
    value <- nloglike(S1=par)
    R <- list(par=par,value=value)
    return(R)
  }

  L$gauss$K <- 2
  L$gauss$COST <- function(par,zero=0) { nloglike(S1=par[1],K1=par[2],zero=zero) }
  L$gauss$GUESS <- function()
  {
    w <- 1-exp(-dof[SUB]) # turn off tiny dof
    w <- w/sum(w)
    mu <- sum(w*s[SUB])
    k <- sum(w/s[SUB]) - 1/mu # VAR/mu^3
    c(mu,k)
  }
  L$gauss$FIT <- function(par=NULL)
  {
    # if guess is not specified
    if(is.null(par)) { par <- L$gauss$GUESS() }
    optimizer(par,L$gauss$COST,lower=c(0,0),upper=c(Inf,Inf))
  }

  L$gauss.gauss$K <- 4

  ## greedy partition fit ## used for good initial guess for mixture model
  part <- function(JOINT=NULL,LEFT=JOINT,RIGHT=JOINT)
  {
    NLL <- array(NA,n) # negative log-likelihoods
    MIN <- L[[LEFT]]$K
    MAX <- n - L[[RIGHT]]$K
    i.min <- i <- MAX # current partition at i + 1/2
    MID <- round((MIN+MAX)/2)

    if(!is.null(JOINT))
    {
      # zero-ing is not exact here because I'm using individual likelihoods lazily --- this result is just used for guess though
      # par = (S1,S2,CoV2) # COV2 = k*mu
      COST <- function(par,zero=0)
      {
        SUB <<- SORT[1:i]
        FRAC <- length(SUB)/n
        NLL <-       nloglike(S1=par[1],K1=par[3]/par[1],zero=zero*FRAC)
        SUB <<- SORT[(i+1):n]
        FRAC <- length(SUB)/n
        NLL <- NLL + nloglike(S1=par[2],K1=par[3]/par[2],zero=zero*FRAC)
        return(NLL)
      }
    }

    NLL.MIN <- Inf
    for(i in c(MID:MAX,MID:MIN)) # work from inside out for better parameter fitting
    {
      if(i==MID) { par <- par1 <- par2 <- NULL } # reset initial guess at mid-point # recycle as you go out

      if(is.null(JOINT)) # this code handles 4 different cases
      {
        # left fit
        SUB <<- SORT[1:i]
        FIT <- L[[LEFT]]$FIT(par1)
        par1 <- FIT$par
        NLL[i] <- FIT$value

        # right fit
        SUB <<- SORT[(i+1):n]
        FIT <- L[[RIGHT]]$FIT(par2)
        par2 <- FIT$par
        NLL[i] <- NLL[i] + FIT$value
      }
      else # this is hard-coded to the joint-N model
      {
        # guessed parameters (ave k=VAR/mu^3)
        if(is.null(par))
        {
          SUB <<- SORT[1:i]
          par1 <- L$gauss$GUESS()
          SUB <<- SORT[(i+1):n]
          par2 <- L$gauss$GUESS()
          CoV2 <- mean( prod(par1), prod(par2) )
          par <- c(par1[1],par2[1],CoV2)
        }
        FIT <- optimizer(par,COST,lower=c(0,0,0),upper=c(Inf,Inf,Inf))
        par <- FIT$par
        NLL[i] <- FIT$value
      }

      if(NLL[i] < NLL.MIN) # store best fit
      {
        i.min <- i
        par1.min <- par1
        par2.min <- par2
        par.min <- par
        NLL.MIN <- NLL[i]
      }
    } # end for
    SUB <<- 1:n

    # index -> fraction
    p <- (i.min+0.5)/n
    # all parameters (for guess)
    if(is.null(JOINT))
    { par <- c(p,par1.min,par2.min) }
    else
    { par <- c(p,par.min) }

    return(par)
  } # end part()

  ################
  # model selection
  M <- 7 # total number of candidate models
  NAME <- array("",M)
  COST <- MLE <- PAR <- BC <- PARS <- LOWER <- UPPER <- list()
  LL <- K <- Kc <- array(0,M)
  if(!is.na(IC))
  {
    # (1) Dirac-delta ## exactly solvable
    m <- 1 # model index
    NAME[m] <- "Dirac-\u03B4"
    PARS[[m]] <- c("S1")
    COST[[1]] <- L$dirac$COST
    LOWER[[m]] <- c(0)
    UPPER[[m]] <- c(Inf)
    FIT <- L$dirac$FIT()
    LL[m] <- -FIT$value # log-likelihood at MLE
    MLE[[m]] <- PAR[[m]] <- FIT$par
    names(MLE[[m]]) <- names(PAR[[m]]) <- PARS[[m]]
    BC[[m]] <- rep(1,length(PARS[[m]])); names(BC[[m]]) <- PARS[[m]]

    # (2) inverse-Gaussian # unimodal likelihood
    m <- 2
    NAME[m] <- "inverse-Gaussian"
    PARS[[m]] <- c("S1","K1") # V==1/DOF
    COST[[m]] <- L$gauss$COST
    LOWER[[m]] <- c(0,0)
    UPPER[[m]] <- c(Inf,Inf)
    FIT <- L$gauss$FIT()
    LL[m] <- -FIT$value
    MLE[[m]] <- PAR[[m]]<- FIT$par
    names(MLE[[m]]) <- names(PAR[[m]]) <- PARS[[m]]
    BC[[m]] <- rep(1,length(PARS[[m]])); names(BC[[m]]) <- PARS[[m]]

    if(debias) # debias k=VAR/mu^3 parameter
    {
      N <- max( length(s) ,2)
      BC[[m]]["K1"] <-  N/(N-1)
      PAR[[m]] <- BC[[m]] * PAR[[m]]
      LL[m] <- -COST[[m]](PAR[[m]])
    }

    # (3) Dirac-delta & Dirac-delta
    m <- 3
    NAME[m] <- "Dirac-\u03B4 + Dirac-\u03B4"
    PARS[[m]] <- c("P1","S1","S2")
    par <- part(LEFT="dirac",RIGHT="dirac") # guess == partitioned solution
    COST[[m]] <- function(par,zero=0) { nloglike(w1=par[1],S1=par[2],S2=par[3],zero=zero) }
    LOWER[[m]] <- c(0,0,0)
    UPPER[[m]] <- c(1,Inf,Inf)
    FIT <- optimizer(par,COST[[m]],lower=LOWER[[m]],upper=UPPER[[m]])
    LL[m] <- -FIT$value
    MLE[[m]] <- PAR[[m]] <- FIT$par
    names(MLE[[m]]) <- names(PAR[[m]]) <- PARS[[m]]
    BC[[m]] <- rep(1,length(PARS[[m]])); names(BC[[m]]) <- PARS[[m]]

    # (4) Dirac-delta & inverse-Gaussian
    m <- 4
    NAME[m] <- "Dirac-\u03B4 + inverse-Gaussian"
    PARS[[m]] <- c("P1","S1","S2","K2")
    par <- part(LEFT="dirac",RIGHT="gauss") # guess == partitioned solution
    COST[[m]] <- function(par,zero=0) { nloglike(w1=par[1],S1=par[2],S2=par[3],K2=par[4],zero=zero) }
    LOWER[[m]] <- c(0,0,0,0)
    UPPER[[m]] <- c(1,Inf,Inf,Inf)
    FIT <- optimizer(par,COST[[m]],lower=LOWER[[m]],upper=UPPER[[m]])
    LL[m] <- -FIT$value
    MLE[[m]] <- PAR[[m]] <- FIT$par
    names(MLE[[m]]) <- names(PAR[[m]]) <- PARS[[m]]
    BC[[m]] <- rep(1,length(PARS[[m]])); names(BC[[m]]) <- PARS[[m]]

    if(debias) # debias k=VAR/mu^3 parameter by membership
    {
      N <- max( length(s) * (1-PAR[[m]]["P1"]) , 2)
      BC[[m]]["K2"] <-  N/(N-1)
      PAR[[m]] <- BC[[m]] * PAR[[m]]
      LL[m] <- -COST[[m]](PAR[[m]])
    }

    # (5) inverse-Gaussian & Dirac-delta
    m <- 5
    NAME[m] <- "inverse-Gaussian + Dirac-\u03B4"
    PARS[[m]] <- c("P1","S1","K1","S2")
    par <- part(LEFT="gauss",RIGHT="dirac") # guess == partitioned solution
    COST[[m]] <- function(par,zero=0) { nloglike(w1=par[1],S1=par[2],K1=par[3],S2=par[4],K2=0,zero=zero) }
    LOWER[[m]] <- c(0,0,0,0)
    UPPER[[m]] <- c(1,Inf,Inf,Inf)
    FIT <- optimizer(par,COST[[m]],lower=LOWER[[m]],upper=UPPER[[m]])
    LL[m] <- -FIT$value
    MLE[[m]] <- PAR[[m]] <- FIT$par
    names(MLE[[m]]) <- names(PAR[[m]]) <- PARS[[m]]
    BC[[m]] <- rep(1,length(PARS[[m]])); names(BC[[m]]) <- PARS[[m]]

    if(debias) # debias k=VAR/mu^3 parameter by membership
    {
      N <- max( length(s) * PAR[[m]]["P1"] , 2)
      BC[[m]]["K1"] <-  N/(N-1)
      PAR[[m]] <- BC[[m]] * PAR[[m]]
      LL[m] <- -COST[[m]](PAR[[m]])
    }

    # (6) inverse-Gaussian & inverse-Gaussian (same relative variance)
    m <- 6
    NAME[m] <- "inverse-Gaussian(CoV) + inverse-Gaussian(CoV)"
    PARS[[m]] <- c("P1","S1","S2","RV")
    par <- part(JOINT="gauss") # guess == partitioned solution
    COST[[m]] <- function(par,zero=0) { nloglike(w1=par[1],S1=par[2],K1=par[4]/par[2],S2=par[3],K2=par[4]/par[3],zero=zero) }
    LOWER[[m]] <- c(0,0,0,0)
    UPPER[[m]] <- c(1,Inf,Inf,Inf)
    FIT <- optimizer(par,COST[[m]],lower=LOWER[[m]],upper=UPPER[[m]])
    LL[m] <- -FIT$value
    MLE[[m]] <- PAR[[m]] <- FIT$par
    names(MLE[[m]]) <- names(PAR[[m]]) <- PARS[[m]]
    BC[[m]] <- rep(1,length(PARS[[m]])); names(BC[[m]]) <- PARS[[m]]

    if(debias) # debias k=VAR/mu^3-ish parameter... does this make sense?
    {
      N <- max( length(s) ,2)
      BC[[m]]["RV"] <- N/(N-1) # numerator and denominator weighted separately # half as big as model with two CoV
      PAR[[m]] <- BC[[m]] * PAR[[m]]
      LL[m] <- -COST[[m]](PAR[[m]])
    }
  }

  # (7) inverse-Gaussian & inverse-Gaussian (independent parameters)
  m <- 7
  NAME[m] <- "inverse-Gaussian(CoV\u2081) + inverse-Gaussian(CoV\u2082)"
  PARS[[m]] <- c("P1","S1","K1","S2","K2")
  par <- part(LEFT="gauss",RIGHT="gauss") # guess == partitioned solution
  COST[[m]] <- function(par,zero=0) { nloglike(w1=par[1],S1=par[2],K1=par[3],S2=par[4],K2=par[5],zero=zero) }
  LOWER[[m]] <- c(0,0,0,0,0)
  UPPER[[m]] <- c(1,Inf,Inf,Inf,Inf)
  FIT <- optimizer(par,COST[[m]],lower=LOWER[[m]],upper=UPPER[[m]])
  LL[m] <- -FIT$value
  MLE[[m]] <- PAR[[m]] <- FIT$par
  names(MLE[[m]]) <- names(PAR[[m]]) <- PARS[[m]]
  BC[[m]] <- rep(1,length(PARS[[m]])); names(BC[[m]]) <- PARS[[m]]

  if(debias) # debias k=VAR/mu^3 parameters by membership
  {
    N <- max( length(s) * PAR[[m]]["P1"] ,2)
    BC[[m]]["K1"] <- N/(N-1)
    N <- max( length(s) * (1-PAR[[m]]["P1"]) ,2)
    BC[[m]]["K2"] <- N/(N-1)
    PAR[[m]] <- BC[[m]] * PAR[[m]]
    LL[m] <- -COST[[m]](PAR[[m]])
  }

  if(is.na(IC)) { IND <- M } else { IND <- 1:M }
  for(m in IND)
  {
    K[m] <- length(PARS[[m]])  # number of parameters
    K.mu <- sum(grepl("S[0-9]",PARS[[m]])) # mean parameters
    K.nu <- K[m] - K.mu
    Kc[m] <- K[m]*ifelse(debias,n-K.mu,n)/pmax(n-K[m]-K.nu,0)
  }

  if(!is.na(IC))
  {
    ## select best performing model according to IC
    if(IC=="AIC")
    { dIC <- 2*K - 2*LL }
    else if(IC=="AICc")
    { dIC <- 2*Kc - 2*LL }
    else if(IC=="BIC")
    { dIC <- log(n)*K - 2*LL }

    dIC <- cbind(dIC)
    rownames(dIC) <- NAME
    colnames(dIC) <- paste0("\u0394",IC)
    MIN <- which.min(dIC) # need for later
    dIC <- dIC - min(dIC)
    IND <- sort(c(dIC),index.return=TRUE)$ix
    dIC <- dIC[IND,,drop=FALSE]
    print(dIC)
  }
  else
  { MIN <- 7 }

  par <- PAR[[MIN]] # selected parameters
  PARS <- PARS[[MIN]]
  names(par) <- PARS
  BC <- BC[[MIN]]

  COST <- COST[[MIN]]
  LOWER <- LOWER[[MIN]]
  UPPER <- UPPER[[MIN]]

  SUB <- 1:n
  if(MIN==1) # chi^2 case
  {
    COV <- 2*par["S1"]^2/sum(dof)
    COV <- array(COV,c(1,1))
  }
  else # numeric cases
  {
    COV <- genD(MLE[[MIN]],COST,lower=LOWER,upper=UPPER)
    COV <- cov.loglike(COV$hessian,COV$gradient)
  }
  dimnames(COV) <- list(PARS,PARS)

  if(debias) { COV <- BC * t(BC * COV) }

  ### RETURN PARAMETERS ###
  CI <- array(NA_real_,c(7,3))
  # rownames(CI) <- c("mu1","CoV1","mu2","CoV2","P1","P2","mu2/mu1")
  rownames(CI) <- c("\u03BC\u2081","CoV\u2081","\u03BC\u2082","CoV\u2082","P\u2081","P\u2082","\u03BC\u2082/\u03BC\u2081")
  colnames(CI) <- NAMES.CI
  INF <- c(0,0,Inf)

  # swap if sub-population means out of order
  if("S2" %in% PARS && par["S1"]>par["S2"])
  {
    RARS <- PARS
    RARS["P1"] <- 1 - PARS["P1"]
    RARS[PARS=="S1"] <- "S2"
    RARS[PARS=="S2"] <- "S1"
    if("K1" %in% PARS) { RARS[PARS=="K1"] <- "K2" }
    if("K2" %in% PARS) { RARS[PARS=="K2"] <- "K1" }
    PARS <- RARS

    names(par) <- PARS
    dimnames(COV) <- list(PARS,PARS)
  }

  S1 <- par["S1"]

  ### first sub-population ###
  if(all(c("K1","RV") %nin% PARS)) # Dirac-delta on first sub-population
  {
    CI[1,] <- chisq.ci(par["S1"],VAR=COV["S1","S1"],level=level)
    CI[2,] <- INF

    K1 <- 0
  }
  else # inverse-Gaussian on first sub-population
  {
    # CoV^2 (RVAR)
    if("K1" %in% PARS) # unstructured
    {
      # CoV^2 (RVAR) : par[1]*par[2]
      GRAD <- par[c("K1","S1")]
      VAR <- c(GRAD %*% COV[c("S1","K1"),c("S1","K1")] %*% GRAD)
      CI[2,] <- chisq.ci(par["S1"]*par["K1"],VAR=VAR,level=level)

      K1 <- par["K1"]
    }
    else if("RV" %in% PARS) # constrained relative variance
    {
      # CoV^2 (RVAR)
      VAR <- COV["RV","RV"]
      CI[2,] <- chisq.ci(par["RV"],VAR=VAR,level=level)

      K1 <- par["RV"]/par["S1"]
      K2 <- par["RV"]/par["S2"]
    }

    # CoV (RSTD)
    CI[2,] <- sqrt(CI[2,])
    if(debias && 0<CI[2,2]) # sqrt bias
    {
      BIAS <- 1 + VAR/CI[2,2]^2/4
      CI[2,] <- CI[2,] * sqrt(BIAS)
    }

    if("RV" %in% PARS) { CI[4,] <- CI[2,] } # copy CoV

    # mu_1
    N <- ifelse("P1" %in% PARS,par["P1"],1) * n
    DD.IG <- DD.IG.ratio(c(S1,K1),COV["S1","S1"],N)
    CI[1,] <- chisq.IG.ci(par["S1"],COV["S1","S1"],DD.IG,level=level,precision=precision)
  } # end inverse-Gaussian on first sub-population

  if("S2" %nin% PARS) # one sub-population only
  {
    # copy sub-population
    CI[3:4,] <- CI[1:2,]
    CI[5,] <- c(0,1,1)
    CI[6,] <- c(0,0,1)
    CI[7,] <- c(0,1,Inf)

    P <- rep(1,length(s))
  }
  else if("S2" %in% PARS) # two sub-populations
  {
    S2 <- par["S2"]

    if(all(c("K2","RV") %nin% PARS)) # Dirac-delta on first sub-population
    {
      # mu_2
      CI[3,] <- chisq.ci(par["S2"],VAR=COV["S2","S2"],level=level)
      CI[4,] <- INF

      K2 <- 0
    }
    else # IG
    {
      if("K2" %in% PARS) # unstructured IG
      {
        # CoV^2 (RVAR) : par[1]*par[2]
        GRAD <- par[c("K2","S2")]
        VAR <- c(GRAD %*% COV[c("S2","K2"),c("S2","K2")] %*% GRAD)
        CI[4,] <- chisq.ci(par["S2"]*par["K2"],VAR=VAR,level=level)

        # CoV (RSTD)
        CI[4,] <- sqrt(CI[4,])
        if(debias && 0<CI[4,2]) # sqrt bias
        {
          BIAS <- 1 + VAR/CI[4,2]^2/4
          CI[4,] <- CI[4,] * sqrt(BIAS)
        }

        K2 <- par["K2"]
      }
      # constrained RV already finished and copied over

      # mu_2
      DD.IG <- DD.IG.ratio(c(S2,K2),COV["S2","S2"],(1-par["P1"])*n)
      CI[3,] <- chisq.IG.ci(par["S2"],COV["S2","S2"],DD.IG,level=level,precision=precision)
    }

    CI[5,] <- beta.ci(par["P1"],VAR=COV["P1","P1"],level=level)
    CI[6,] <- beta.ci(1-par["P1"],VAR=COV["P1","P1"],level=level)
    P1 <- par["P1"]

    # this is MVU for both chi^2 and IG
    # this is 1st order debiased in between
    BIAS <- max(1-debias*COV["S1","S1"]/S1^2,0)
    RATIO <- par["S2"]/par["S1"] * BIAS
    GRAD <- c(-par["S2"]/par["S1"]^2,1/par["S1"])
    GRAD[1] <- GRAD[1] / BIAS # scale S1 parameter and its variance the same (approximate, but not exact correction)
    # numerator and denominator are correlated...
    VAR.RATIO <- c(GRAD %*% COV[c("S1","S2"),c("S1","S2")] %*% GRAD)
    # not sure of better choice
    CI[7,] <- IG.ci(RATIO,VAR=VAR.RATIO,level=level,precision=precision)

    P <- rep(1,length(s))
    for(i in 1:n)
    {
      SUB <- i
      P[i] <- exp( -nloglike(w1=P1,S1=S1,K1=K1,w2=0,S2=S2,K2=K2) + nloglike(w1=P1,S1=S1,K1=K1,w2=1-P1,S2=S2,K2=K2) )
    }
    SUB <- 1:length(s)
  }

  STUFF <- list()
  STUFF$CI <- CI
  STUFF$P <- P
  return(STUFF)
}
