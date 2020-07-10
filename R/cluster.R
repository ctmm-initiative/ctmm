cluster <- function(x,level=0.95,level.UD=0.95,level.pop=0.95,robust=FALSE,IC="BIC",boot=FALSE,error=0.01,debias=TRUE,units=TRUE,plot=TRUE,...)
{
  cluster.area(x=x,level=level,level.UD=level.UD,level.pop=level.pop,IC=IC,boot=boot,error=error,debias=debias,robust=robust,units=units,plot=plot,...)
}


# wrapper: meta-analysis of CTMM areas
# TODO range=FALSE ???
cluster.area <- function(x,level=0.95,level.UD=0.95,robust=FALSE,IC="BIC",units=TRUE,classify=TRUE,plot=TRUE,...)
{
  N <- length(x)
  STUFF <- import.area(x,level.UD=level.UD)
  AREA <- STUFF$AREA
  DOF <- STUFF$DOF
  ID <- STUFF$ID

  # inverse-chi^2 population distribution
  STUFF <- cluster.chisq(AREA,DOF,level=level,IC=IC,robust=robust)
  CI <- STUFF$CI
  P <- STUFF$P # probability of falling into second class
  names(P) <- ID

  AREA <- sapply(1:N,function(i){chisq.ci(AREA[i],DOF=DOF[i],level=level)})
  AREA <- t(AREA)

  # colors # black=normative # red=extreme
  COL <- grDevices::rgb(1-P,0,0)

  if(plot)
  {
    # first population at bottom
    ID <- c(ID,ifelse(robust,expression(median[1]),expression(mean[1])))
    AREA <- rbind(AREA,CI[1,])
    COL <- c(COL,grDevices::rgb(0,0,0))

    # second population at top
    ID <- c(ifelse(robust,expression(median[2]),expression(mean[2])),ID)
    AREA <- rbind(CI[2,],AREA)
    COL <- c(grDevices::rgb(1,0,0),COL)

    # basic forest plot
    IND <- (N+2):1

    UNITS <- unit(AREA,"area",SI=!units)
    PLOT <- AREA/UNITS$scale

    xlab <- paste0(100*level.UD,"% Area (",UNITS$name,")")
    plot(range(PLOT),c(1,N+2),col=grDevices::rgb(1,1,1,0),xlab=xlab,ylab=NA,yaxt="n",...)
    for(i in 2:(N+1)) { graphics::axis(2,at=IND[i],labels=ID[i],las=2,col.axis=COL[i],lwd=0) }
    graphics::axis(2,at=IND[N+2],labels=ID[N+2],las=2,font=2,col.axis=COL[N+2],lwd=0)
    graphics::axis(2,at=IND[1],labels=ID[1],las=2,font=2,col.axis=COL[1],lwd=0)
    graphics::axis(2,at=IND,labels=FALSE)
    graphics::abline(v=PLOT[N+2,2],col=grDevices::rgb(0.0,0.0,0.0,0.5))
    graphics::abline(v=PLOT[1,2],col=grDevices::rgb(1,0.0,0.0,0.5))
    graphics::points(PLOT[,2],IND,pch=16,col=COL)
    suppressWarnings( graphics::arrows(PLOT[,1],IND,PLOT[,3],IND,length=0.05,angle=90,code=3,col=COL) )
  }

  UNITS <- unit(CI[1:2,],"area",SI=!units)
  CI[1:2,] <- CI[1:2,]/UNITS$scale
  rownames(CI)[1:2] <- paste0(rownames(CI)[1:2]," area (",UNITS$name,")")

  if(classify)
  { return(list("P"=P,CI=CI)) }
  else
  { return(CI) }
}


# cluster-analysis of chi^2 random variables with mixture of 2 inverse-chi^2 priors
cluster.chisq <- function(s,dof,level=0.95,level.pop=0.95,IC="BIC",boot=FALSE,error=0.01,debias=TRUE,robust=TRUE,precision=1/2,...)
{
  #tol <- .Machine$double.eps^precision
  n <- length(s)
  SORT <- sort(s,index.return=TRUE)$ix

  # MLE fit and CI return
  SUB <- 1:n
  # negative log-likelihood
  nloglike <- function(w1=1,S1,DOF1=Inf,w2=1-w1,S2=S1,DOF2=DOF1,zero=0)
  {
    # subset for CV
    s <- s[SUB]
    dof <- dof[SUB]
    n <- length(s)

    if(w1>0 && (DOF1<=0 || S1==0 || S1==Inf)) { return(Inf) }
    if(w2>0 && (DOF2<=0 || S2==0 || S2==Inf)) { return(Inf) }

    zero <- zero + sum( (dof/2-1)*log(dof*s) + log(dof) ) # common factors

    loglike1 <- -lbetaplog(dof/2,DOF1/2,s/S1) - dof*(log(2*S1)/2) - (dof*s)/(2*S1)*log1pxdx((dof*s)/(DOF1*S1))
    if(w1==1)
    {
      loglike1 <- loglike1 + zero/n
      loglike1 <- sum(loglike1)
      return(-loglike1)
    }
    loglike2 <- -lbetaplog(dof/2,DOF2/2,s/S2) - dof*(log(2*S2)/2) - (dof*s)/(2*S2)*log1pxdx((dof*s)/(DOF2*S2))

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
    loglike <- logw.max + loglike.max + log(1+(w.min/w.max)*exp(loglike.min-loglike.max)) + zero/n
    -sum(loglike)
  }

  # lists of stuff organized by models
  L <- list()
  L$dirac <- list()       # 1
  L$chisq <- list()       # 2
  L$chisq.chisq <- list() # 6

  L$dirac$K <- 1 # number of parameters
  L$dirac$COST <- function(par) { nloglike(S1=par) } # negative log-likelihood
  L$dirac$FIT <- function(par=NULL) # fit MLE
  {
    par <- sum(dof[SUB]*s[SUB])/sum(dof[SUB])
    value <- nloglike(S1=par)
    R <- list(par=par,value=value)
    return(R)
  }

  L$chisq$K <- 2
  L$chisq$COST <- function(par,zero=0) { nloglike(S1=par[1],DOF1=1/par[2],zero=zero) }
  L$chisq$GUESS <- function()
  {
    S <- 1/mean(1/s[SUB])
    DOF <- 2/(S^2*stats::var(1/s[SUB]))
    c(S,1/DOF)
  }
  L$chisq$FIT <- function(par=NULL)
  {
    # if guess is not specified
    if(is.null(par)) { par <- L$chisq$GUESS() }
    optimizer(par,L$chisq$COST,lower=c(0,0),upper=c(Inf,Inf))
  }

  L$chisq.chisq$K <- 3

  ## greedy partition fit ## used for good initial guess for mixture model
  part <- function(JOINT=NULL,LEFT=JOINT,RIGHT=JOINT)
  {
    NLL <- NULL # negative log-likelihood
    i <- ceiling(n/2) -> IND # current partition at i + 1/2
    MIN <- L[[LEFT]]$K
    MAX <- n - L[[RIGHT]]$K - 1

    if(!is.null(JOINT))
    {
      # zero-ing is not exact here because I'm using individual likelihoods lazily --- this result is just used for guess though
      COST <- function(par,zero=0)
      {
        SUB <<- SORT[1:i]
        FRAC <- length(SUB)/n
        NLL <- nloglike(S1=par[1],DOF1=1/par[2],zero=zero*FRAC)
        SUB <<- SORT[(i+1):n]
        FRAC <- length(SUB)/n
        NLL <- NLL + nloglike(S1=par[3],DOF1=1/par[2],zero=zero*FRAC)
        return(NLL)
      }
    }

    # initial fit
    j <- 1 # current list index
    IND <- NULL # all indices condidered thus far
    par <- par1 <- par2 <- PAR <- PAR1 <- PAR2 <- NULL # parameters
    while(!length(IND) || (IND[j]==min(IND) && IND[j]>MIN) || (IND[j]==max(IND) && IND[j]<MAX))
    {
      IND <- c(IND,i)
      j <- length(IND)

      if(is.null(JOINT)) # this code handles 4 different cases
      {
        # left fit
        SUB <<- SORT[1:i]
        FIT <- L[[LEFT]]$FIT(par1)
        PAR1 <- rbind(PAR1,FIT$par)
        NLL[j] <- FIT$value

        # right fit
        SUB <<- SORT[(i+1):n]
        FIT <- L[[RIGHT]]$FIT(par2)
        PAR2 <- rbind(PAR2,FIT$par)
        NLL[j] <- NLL[j] + FIT$value

        j <- which.min(NLL) # current best index
        par1 <- PAR1[j,] # current best parameters
        par2 <- PAR2[j,]
      }
      else # this is hard-coded to the joint-N model
      {
        # guessed parameters (ave 1/DOF)
        if(is.null(par))
        {
          SUB <<- SORT[1:i]
          par1 <- L$chisq$GUESS()
          SUB <<- SORT[(i+1):n]
          par2 <- L$chisq$GUESS()
          par <- c(par1[1],mean(par1[2],par2[2]),par2[1])
        }
        FIT <- optimizer(par,COST,lower=c(0,0,0),upper=c(Inf,Inf,Inf))
        PAR <- rbind(PAR,FIT$par)
        NLL[j] <- FIT$value

        j <- which.min(NLL) # current best index
        par <- PAR[j,] # current best parameters
      }

      # next attempt # bias right/high
      if(IND[j]==max(IND)) { i <- IND[j] + 1 } # move partition higher
      else if(IND[j]==min(IND)) { i <- IND[j] - 1 } # move partition lower
      else { break } # local minima reached
    }
    SUB <<- 1:n

    # all parameters (for guess)
    if(is.null(JOINT)) { par <- c(IND[j]/n,par1,par2) }
    else { par <- c(IND[j]/n,par) }
    return(par)
  }

  ################
  # model selection
  M <- 7 # total number of candidate models
  NAME <- array("",M)
  COST <- PAR <- PARS <- LOWER <- UPPER <- list()
  LL <- K <- array(NA_real_,M)
  if(TRUE && !is.na(IC))
  {
    # (1) Dirac delta ## exactly solvable
    m <- 1 # model index
    NAME[m] <- "Dirac-\u03B4"
    PARS[[m]] <- c("S1")
    COST[[1]] <- L$dirac$COST
    FIT <- L$dirac$FIT()
    LL[m] <- -FIT$value # log-likelihood at MLE
    PAR[[m]] <- FIT$par
    K[m] <- 1 # number of parameters

    # (2) inverse chi^2 # unimodal likelihood
    m <- 2
    NAME[m] <- "inverse-\u03C7\u00B2(N)"
    PARS[[m]] <- c("S1","V1") # V==1/DOF
    COST[[2]] <- L$chisq$COST
    FIT <- L$chisq$FIT()
    LL[m] <- -FIT$value
    PAR[[m]] <- FIT$par
    K[m] <- 2

    # (3) Dirac delta & Dirac delta
    m <- 3
    NAME[m] <- "Dirac-\u03B4 + Dirac-\u03B4"
    PARS[[m]] <- c("P1","S1","S2")
    par <- part(LEFT="dirac",RIGHT="dirac") # guess == partitioned solution
    COST[[m]] <- function(par,zero=0) { nloglike(w1=par[1],S1=par[2],S2=par[3],zero=zero) }
    LOWER[[m]] <- c(0,0,0)
    UPPER[[m]] <- c(1,Inf,Inf)
    FIT <- optimizer(par,COST[[m]],lower=LOWER[[m]],upper=UPPER[[m]])
    LL[m] <- -FIT$value
    PAR[[m]] <- FIT$par
    K[m] <- 3

    # (4) Dirac delta & inverse chi^2
    m <- 4
    NAME[m] <- "Dirac-\u03B4 + inverse-\u03C7\u00B2(N)"
    PARS[[m]] <- c("P1","S1","S2","V2")
    par <- part(LEFT="dirac",RIGHT="chisq") # guess == partitioned solution
    COST[[m]] <- function(par,zero=0) { nloglike(w1=par[1],S1=par[2],S2=par[3],DOF2=1/par[4],zero=zero) }
    LOWER[[m]] <- c(0,0,0,0)
    UPPER[[m]] <- c(1,Inf,Inf,Inf)
    FIT <- optimizer(par,COST[[m]],lower=LOWER[[m]],upper=UPPER[[m]])
    LL[m] <- -FIT$value
    PAR[[m]] <- FIT$par
    K[m] <- 4

    # (5) inverse chi^2 & Dirac delta
    m <- 5
    NAME[m] <- "inverse-\u03C7\u00B2(N) + Dirac-\u03B4"
    PARS[[m]] <- c("P1","S1","V2","S2")
    par <- part(LEFT="chisq",RIGHT="dirac") # guess == partitioned solution
    COST[[m]] <- function(par,zero=0) { nloglike(w1=par[1],S1=par[2],DOF1=1/par[3],S2=par[4],DOF2=Inf,zero=zero) }
    LOWER[[m]] <- c(0,0,0,0)
    UPPER[[m]] <- c(1,Inf,Inf,Inf)
    FIT <- optimizer(par,COST[[m]],lower=LOWER[[m]],upper=UPPER[[m]])
    LL[m] <- -FIT$value
    PAR[[m]] <- FIT$par
    K[m] <- 4

    # (6) inver chi^2 & inverse chi^2 (same relative variance)
    m <- 6
    NAME[m] <- "inverse-\u03C7\u00B2(N) + inverse-\u03C7\u00B2(N)"
    PARS[[m]] <- c("P1","S1","V","S2")
    par <- part(JOINT="chisq") # guess == partitioned solution
    COST[[m]] <- function(par,zero=0) { nloglike(w1=par[1],S1=par[2],DOF1=1/par[3],S2=par[4],DOF2=1/par[3],zero=zero) }
    LOWER[[m]] <- c(0,0,0,0)
    UPPER[[m]] <- c(1,Inf,Inf,Inf)
    FIT <- optimizer(par,COST[[m]],lower=LOWER[[m]],upper=UPPER[[m]])
    LL[m] <- -FIT$value
    PAR[[m]] <- FIT$par
    K[m] <- 4

    # (7) inver chi^2 & inverse chi^2
    m <- 7
    NAME[m] <- "inverse-\u03C7\u00B2(N\u2081) + inverse-\u03C7\u00B2(N\u2082)"
    PARS[[m]] <- c("P1","S1","V1","S2","V2")
    par <- part(LEFT="chisq",RIGHT="chisq") # guess == partitioned solution
    COST[[m]] <- function(par,zero=0) { nloglike(w1=par[1],S1=par[2],DOF1=1/par[3],S2=par[4],DOF2=1/par[5],zero=zero) }
    LOWER[[m]] <- c(0,0,0,0,0)
    UPPER[[m]] <- c(1,Inf,Inf,Inf,Inf)
    FIT <- optimizer(par,COST[[m]],lower=LOWER[[m]],upper=UPPER[[m]])
    LL[m] <- -FIT$value
    PAR[[m]] <- FIT$par
    K[m] <- 5
  }

  ## select best performing model according to IC
  if(IC=="AIC")
  { dIC <- 2*K - 2*LL }
  else if(IC=="BIC")
  { dIC <- log(n)*K - 2*LL }
  else if(IC=="LOOCV")
  {
    # TODO
  }

  if(!is.na(IC))
  {
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

  COST <- COST[[MIN]]
  LOWER <- LOWER[[MIN]]
  UPPER <- UPPER[[MIN]]

  SUB <- 1:n
  COV <- genD(par,COST,lower=LOWER,upper=UPPER)
  COV <- cov.loglike(COV$hessian,COV$gradient)
  dimnames(COV) <- list(PARS,PARS)

  # DOF<=2 warnings
  V <- c("V1","V2","V")
  V <- V[V %in% PARS]
  if(!robust && length(V) && any(1/par[V]<=2))
  { warning("Population mean not convergent. Consider robust=TRUE.") }

  pop.ave <- function(S,V)
  {
    Q <- function(par)
    {
      q <- par[S]
      if(length(V))
      {
        if(!robust) { q <- q /max(1-2*par[V],0) } # inverse-chi^2 bias # N/(N-2) == 1/(1-2/N)
        else { q <- q * 1/stats::qchisq(0.5,1/par[V])/par[V] } # inverse-chi^2(N) median N/Q
      }
      return(q)
    }
    MLE <- Q(par)
    GRAD <- genD(par,Q,lower=LOWER,upper=UPPER,order=1)$gradient
    VAR <- c( GRAD %*% COV %*% GRAD )
    chisq.ci(MLE,COV=VAR,level=level)
  }

  # population-1 average
  V <- c("V1","V")
  V <- V[V %in% PARS]
  CI.1 <- pop.ave("S1",V)

  if(MIN>2)
  {
    # population-2 average
    V <- c("V2","V")
    V <- V[V %in% PARS]
    CI.2 <- pop.ave("S2",V)

    ## % of population in lower component
    MLE <- par["P1"]
    VAR <- COV["P1","P1"]
    CI.w <- beta.ci(MLE,VAR,level=level)

    # are sub-population means out of order !!!

    ## calculate F-like statistic
    if(!robust) # ratio of mean areas
    {
      V <- c("V2","V") # numerator (S2) is inverse-chi^2
      V <- V[V %in% PARS]
      Q <- function(par)
      {
        q <- par["S2"]/par["S1"]
        if(length(V)) { q <- q * 1/max(1-2*par[V],0) }
        return(q)
      }
    }
    else # ratio of median areas
    {
      Q <- function(par)
      {
        q <- par["S2"]/par["S1"]
        if("V2" %in% PARS) { q <- q * stats::qchisq(0.5,1/par["V2"])/par["V2"] }
        if("V1" %in% PARS) { q <- q * stats::qchisq(0.5,1/par["V1"])*par["V1"] }
        return(q)
      }
    }
    MLE <- Q(par)
    GRAD <- genD(par,Q,lower=LOWER,upper=UPPER,order=1)$gradient
    VAR <- c( GRAD %*% COV %*% GRAD )
    CI.F <- lognorm.ci(MLE,COV=VAR,level=level)

    ## classifications -- membership probabilities
    w1 <- par["P1"]
    S1 <- par["S1"]
    S2 <- par["S2"]
    DOF1 <- 1/par["V1"]
    DOF2 <- 1/par["V2"]
    if(is.na(DOF1)) { DOF1 <- 1/par["V"] }
    if(is.na(DOF2)) { DOF2 <- 1/par["V"] }
    if(is.na(DOF1)) { DOF1 <- Inf }
    if(is.na(DOF2)) { DOF2 <- Inf }
    P <- array(1,n)
    for(i in 1:n)
    {
      SUB <- i
      P[i] <- exp( -nloglike(w1=w1,S1=S1,DOF1=DOF1,w2=0,S2=S2,DOF2=DOF2) + nloglike(w1=w1,S1=S1,DOF1=DOF1,w2=1-w1,S2=S2,DOF2=DOF2) )
    }
  }
  else
  {
    CI.w <- c(0,0,0)
    CI.F <- c(1,1,1)
    CI.2 <- CI.1
    P <- array(1,n)
  }

  CI <- rbind(CI.1,CI.2,CI.w,CI.F)
  rownames(CI) <- c("A\u2081","A\u2082","P\u2081","A\u2082/A\u2081")
  colnames(CI) <- NAMES.CI

  STUFF <- list()
  STUFF$CI <- CI
  STUFF$P <- P
  return(STUFF)
}
