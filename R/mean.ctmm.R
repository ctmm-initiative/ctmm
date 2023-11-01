# M %*% S %*% t(M) -> L %*% S[TRI]
quad2lin <- function(M,diag=FALSE)
{
  M <- rbind(M) # [n,m]
  n <- nrow(M)
  m <- ncol(M)

  TRI <- diag(1,m)
  TRI <- upper.tri(TRI,diag=TRUE)
  TRI <- which(TRI)
  p <- length(TRI)

  if(diag)
  { J <- matrix(0,n,p) }
  else
  { J <- matrix(0,p,p) }

  for(i in 1:p)
  {
    # E <- array(0,dim(M))
    E <- array(0,c(m,m))
    E[TRI[i]] <- 1
    E <- E + t(E)
    E <- E/max(E)

    if(diag) # n!=m and only care about n variances
    {
      # O(n^2) slow method like below
      # E <- M %*% E %*% t(M) # [n,n]
      # E <- diag(E)

      # O(n) fast method
      E <- sapply(1:n,function(i){M[i,] %*% E %*% M[i,]})
    }
    else # n==m and want full covariance
    {
      E <- M %*% E %*% t(M) # [n,n]
      E <- E[TRI]
    }
    J[,i] <- E
  }

  return(J)
}

#############
mean.features <- function(x,debias=TRUE,weights=NULL,trace=FALSE,IC="AICc",select="all",formula=FALSE,...)
{
  if(is.null(weights))
  { weights <- rep(1,length(x)) }

  N <- length(x)
  axes <- x[[1]]$axes

  ISO <- sapply(x,function(y){y$isotropic})
  ISO <- all(ISO)

  # only want biological parameters
  for(i in 1:length(x))
  {
    FEATURES <- x[[i]]$features
    ERRORS <- grepl('error',FEATURES)
    if(any(ERRORS))
    {
      x[[i]]$error <- FALSE
      x[[i]]$features <- FEATURES[!ERRORS]
      #FEATURES <- rownames(x[[i]]$COV)
      #ERRORS <- grepl('error',FEATURES) # how can this differ?
      x[[i]]$COV <- x[[i]]$COV[!ERRORS,!ERRORS,drop=FALSE]
    }
  }

  # all possible features
  FEATURES <- lapply(x,function(y){y$features})
  FEATURES <- unlist(FEATURES)
  FEATURES <- unique(FEATURES)

  TAU.EXPAND <- all(c("tau","tau position") %in% FEATURES)
  if(TAU.EXPAND)
  {
    FEATURES <- FEATURES[FEATURES!="tau"]
    if("tau velocity" %nin% FEATURES)
    { FEATURES <- c(FEATURES,"tau velocity") }
  }

  NFEAT <- FEATURES[FEATURES %nin% c('major','minor','angle')]

  M <- length(FEATURES)
  # all possible RSF beta
  BETA <- lapply(x,function(y){names(y$beta)})
  BETA <- unlist(BETA)
  BETA <- unique(BETA)

  MEANS <- VARS <- array(TRUE,M)
  names(MEANS) <- names(VARS) <- FEATURES

  MU <- array(0,c(N,M))
  SIGMA <- array(0,c(N,M,M))
  INF <- diag(Inf,M)
  colnames(MU) <- FEATURES
  dimnames(SIGMA) <- list(NULL,FEATURES,FEATURES)

  # temporary code for linear functional response
  if(class(formula)[1]=="formula")
  { formula <- all.vars(formula) } # now character
  if(class(formula)[1]=="character")
  { formula <- BETA[BETA %in% formula] }
  else
  { formula <- array(formula,length(BETA)) }

  if(length(formula)) # get predictors for model-II regression - linear response function
  {
    TERMS <- formula
    X <- array(0,c(N,length(TERMS)))
    colnames(X) <- TERMS
    SX <- array(0,c(N,length(TERMS),length(TERMS)))
    dimnames(SX) <- list(NULL,TERMS,TERMS)
  }
  else
  { X <- SX <- FALSE }

  for(i in 1:N)
  {
    if(length(formula) && length(x[[i]]$used.mean))
    {
      NAMES <- names(x[[i]]$used.mean)
      X[i,NAMES] <- x[[i]]$used.mean
      SX[i,NAMES,NAMES] <- x[[i]]$used.cov
    }

    x[[i]] <- log_ctmm(x[[i]],debias=debias)
    features <- names(x[[i]]$par)

    if(TAU.EXPAND && "tau" %in% features)
    {
      # expand features
      features[features=="tau"] <- "tau position"
      features <- c(features,"tau velocity")

      # expand point estimates
      NAMES <- names(x[[i]]$par)
      names(x[[i]]$par)[NAMES=="tau"] <- "tau position"
      x[[i]]$par["tau velocity"] <- x[[i]]$par["tau position"]

      # expand covariance
      NAMES <- rownames(x[[i]]$COV)
      NAMES[NAMES=="tau"] <- "tau position"
      dimnames(x[[i]]$COV) <- list(NAMES,NAMES)

      ROW <- x[[i]]$COV['tau position',]
      x[[i]]$COV <- rbind(x[[i]]$COV,'tau velocity'=ROW)
      ROW['tau velocity'] <- ROW['tau position']
      x[[i]]$COV <- cbind(x[[i]]$COV,'tau velocity'=ROW)
    }

    MU[i,features] <- x[[i]]$par
    SIGMA[i,,] <- INF
    SIGMA[i,features,features] <- x[[i]]$COV

    P <- c("major","minor")
    if(!ISO && x[[i]]$isotropic)
    {
      MU[i,'minor'] <- MU[i,'major']
      # perfectly correlated observation
      SIGMA[i,P,P] <- matrix(SIGMA[i,'major','major'],2,2)
      # copy over other correlations as well
      SIGMA[i,'minor',NFEAT] <- SIGMA[i,'major',NFEAT]
      SIGMA[i,NFEAT,'minor'] <- SIGMA[i,NFEAT,'major']
    }
  } # for(i in 1:N)

  ###################
  J <- diag(1,ncol(MU))
  dimnames(J) <- list(FEATURES,FEATURES)

  # everything is already log transformed - diagonalize as much as possible
  if("minor" %in% FEATURES)
  {
    J["major",c("major","minor")] <- c(1/2,1/2) # average variance
    J["minor",c("major","minor")] <- c(1,-1) # difference / eccentricity

    if("tau" %in% FEATURES)
    { J["tau",c("major","minor","tau")] <- c(1/2,1/2,-1) } # D

    if("tau position" %in% FEATURES)
    { J["tau position",c("major","minor","tau position")] <- c(1/2,1/2,-1) } # D
  }
  else
  {
    if("tau" %in% FEATURES)
    { J["tau",c("major","tau")] <- c(1,-1) } # D

    if("tau position" %in% FEATURES)
    { J["tau position",c("major","tau position")] <- c(1,-1) } # D
  }

  # if("tau velocity" %in% FEATURES)
  # { J["tau velocity",c("major","tau position","tau velocity")] <- c(1,-1,-1) } # MSV correlates with D

  # if("omega" %in% FEATURES) # MSV - delta-MSV # sigma cancels out
  # {
  #   if("tau" %in% FEATURES) # OUf
  #   { J["omega",c("tau","omega")] <- c(1,-1) }
  #   else if(all(c("tau position","tau velocity") %in% FEATURES)) # OUF
  #   { J["omega",c("tau position","tau velocity","omega")] <- c(1,1,-2) }
  #   else if("tau velocity" %in% FEATURES) # IOU
  #   { J["omega",c("tau velocity","omega")] <- c(1,-1) }
  #   else if("tau position" %in% FEATURES) # OU
  #   { J["omega",c("tau position","omega")] <- c(1,-1) }
  # }

  # missing data / infinite uncertainties
  INF <- apply(SIGMA,1,function(M){diag(M)==Inf}) # [par,ind]
  dim(INF) <- dim(SIGMA)[2:1]
  INF <- t(INF) # [ind,par]
  colnames(INF) <- colnames(SIGMA)

  tJ <- t(J)
  MU[INF] <- 0 # zero out infinite uncertainties (point estimates could be infinite after log transform)
  MU <- MU %*% tJ

  SIGMA <- SIGMA %.% tJ
  for(i in 1:nrow(SIGMA)) { for(j in 1:ncol(SIGMA)) { if(INF[i,j]) { SIGMA[i,j,] <- SIGMA[i,,j] <- 0; SIGMA[i,j,j] <- Inf } } }
  SIGMA <- aperm(SIGMA,c(1,3,2))
  SIGMA <- SIGMA %.% tJ
  for(i in 1:nrow(SIGMA)) { for(j in 1:ncol(SIGMA)) { if(INF[i,j]) { SIGMA[i,j,] <- SIGMA[i,,j] <- 0; SIGMA[i,j,j] <- Inf } } }
  dimnames(SIGMA) <- list(NULL,FEATURES,FEATURES)

  # can't make this infinite before or will mess up above calculations
  P <- c("minor","angle")
  for(i in 1:length(x))
  {
    if(!ISO && x[[i]]$isotropic)
    {
      INF[i,P] <- TRUE
      SIGMA[i,P,] <- 0
      SIGMA[i,,P] <- 0
      SIGMA[i,P,P] <- diag(Inf,2)
    }
  }

  # number of estimated features (amount of data)
  DIM <- length(FEATURES)

  make.names <- function(MEANS,VARS,OFF=FALSE)
  {
    MEANS <- rbind(MEANS)
    VARS <- rbind(VARS)
    OFF <- OFF && sum(VARS)>1 # no correlations present
    OFF <- ifelse(OFF,"COV","VAR")

    MEANS <- sapply(1:nrow(MEANS),function(i){ paste(FEATURES[MEANS[i,]],collapse=",") })
    VARS <- sapply(1:nrow(VARS),function(i){ paste(FEATURES[VARS[i,]],collapse=",") })

    NAMES <- paste0("E[",MEANS,"] ",OFF,"[",VARS,"]")

    return(NAMES)
  }

  # number of variance-covariance parameters modeled
  num.pars <- function(FIT,mean=FALSE)
  {
    V <- FIT$VARS
    V <- V[upper.tri(V,diag=TRUE)]
    V <- sum(V)
    if(mean) { V <- V + sum(FIT$INT) }
    return(V)
  }

  # data can potentially support these parameters only
  MEAN.MAX <- colSums(!INF)>0
  VAR.MAX <- colSums(!INF)>1
  # these mean parameters can't be turned off
  MEAN.MIN <- rep(TRUE,DIM)
  names(MEAN.MIN) <- FEATURES
  if("minor" %in% FEATURES) { MEAN.MIN[c("minor","angle")] <- FALSE } # can be mean zero
  # MEAN.MIN[BETA] <- FALSE # can be mean-zero

  ##########
  ## build up phase
  ## minimal terms
  OLD <- list()
  MEANS <- MEAN.MIN
  VARS <- rep(FALSE,DIM)
  names(VARS) <- FEATURES

  S <- make.names(MEANS,VARS)
  if(trace) { message("Fitting covariance model ",S) }
  OLD[[S]] <- meta.normal(MU,SIGMA,X=X,SX=SX,debias=debias,weights=weights,VARS=diag(VARS,DIM),MEANS=MEANS,WARN=FALSE,...)

  ## add terms one by one
  GUESS <- OLD[[1]]
  while(any(!VARS) || any(!MEANS))
  {
    NEW <- list()

    ########################
    ## add variance terms ##
    # which terms are presently off
    try <- which(!VARS & VAR.MAX & MEANS) # only add vars to means
    if(length(try))
    {
      TRYS <- array(VARS,c(DIM,length(try)))
      TRYS <- t(TRYS) # [off->on,all]
      colnames(TRYS) <- FEATURES
      for(i in 1:nrow(TRYS)) { TRYS[i,try[i]] <- TRUE }
      if("minor" %in% FEATURES) { TRYS[,"minor"] <- TRYS[,"angle"] <- TRYS[,"minor"] | TRYS[,"angle"] }
      TRYS <- unique(TRYS)

      # models that we will try by turning one term on
      NAMES <- make.names(MEANS,TRYS)
      # don't try what's been done
      SUB <- NAMES %nin% names(OLD)
      SUB <- which(SUB)

      TRYS <- TRYS[SUB,,drop=FALSE]
      NAMES <- NAMES[SUB]

      # fit all
      for(i in 1%:%length(NAMES))
      {
        S <- NAMES[i]
        if(trace) { message("Fitting covariance model ",S) }
        NEW[[S]] <- meta.normal(MU,SIGMA,X=X,SX=SX,debias=debias,weights=weights,VARS=diag(TRYS[i,],DIM),MEANS=MEANS,GUESS=GUESS,WARN=FALSE,...)
      }
    }

    ####################
    ## add mean terms ##
    # which terms are presently off
    try <- which(!MEANS & MEAN.MAX)
    if(length(try))
    {
      TRYS <- array(MEANS,c(DIM,length(try)))
      TRYS <- t(TRYS) # [off->on,all]
      colnames(TRYS) <- FEATURES
      for(i in 1:nrow(TRYS)) { TRYS[i,try[i]] <- TRUE }
      if("minor" %in% FEATURES) { TRYS[,"minor"] <- TRYS[,"angle"] <- TRYS[,"minor"] | TRYS[,"angle"] }
      TRYS <- unique(TRYS)

      # models that we will try by turning one term on
      NAMES <- make.names(TRYS,VARS)
      # don't try what's been done
      SUB <- NAMES %nin% names(OLD)
      SUB <- which(SUB)

      TRYS <- TRYS[SUB,,drop=FALSE]
      NAMES <- NAMES[SUB]

      # fit all
      for(i in 1%:%length(NAMES))
      {
        S <- NAMES[i]
        if(trace) { message("Fitting covariance model ",S) }
        NEW[[S]] <- meta.normal(MU,SIGMA,X=X,SX=SX,debias=debias,weights=weights,VARS=diag(VARS,DIM),MEANS=TRYS[i,],GUESS=GUESS,WARN=FALSE,...)
      }
    }

    ###############
    # wrap up
    if(!length(NEW)) { break } # no new models to fit

    ICS <- sapply(NEW,function(m){m[[IC]]})
    BEST <- which.min(ICS)
    BEST <- NEW[[BEST]]
    GUESS <- BEST
    MEANS <- BEST$INT
    VARS <- diag(BEST$VARS)

    OLD <- c(OLD,NEW)

    if(BEST[[IC]]==Inf) { break }
  } # finish build up variances

  ###########
  ## pair down terms phase
  # start with all variances or best (if IC==Inf)
  NAME.ALL <- make.names(MEAN.MAX,VAR.MAX)
  if(NAME.ALL %in% names(OLD))
  { BEST <- OLD[[NAME.ALL]] }
  else # limited data cannot support all variances
  {
    ICS <- sapply(OLD,function(m){m[[IC]]})
    BEST <- which.min(ICS)
    BEST <- OLD[[BEST]]
  }
  GUESS <- BEST
  MEANS <- BEST$INT
  VARS <- diag(BEST$VARS)

  ##############
  # pair down variances phase
  while(sum(VARS)>1 || sum(MEANS)>1)
  {
    NEW <- list()
    NAMES1 <- NAMES2 <- NULL
    NEW1 <- NEW2 <- NULL

    #################
    ### var terms ###
    # which terms are presently off
    try <- which(VARS)
    if(length(try))
    {
      TRYS <- array(VARS,c(DIM,length(try)))
      TRYS <- t(TRYS) # [off->on,all]
      colnames(TRYS) <- FEATURES
      for(i in 1:nrow(TRYS)) { TRYS[i,try[i]] <- FALSE }
      if("minor" %in% FEATURES) { TRYS[,"minor"] <- TRYS[,"angle"] <- TRYS[,"minor"] & TRYS[,"angle"] }
      TRYS <- unique(TRYS)
      # models that we will try by turning one term off
      NAMES <- NAMES1 <- make.names(MEANS,TRYS)
      # paired down models that we haven't attempted yet
      SUB <- NAMES %nin% names(OLD)
      NEW.NAMES <- NEW1 <- NAMES[SUB]
      TRYS <- TRYS[SUB,,drop=FALSE]

      # fit all
      for(i in 1%:%length(NEW.NAMES))
      {
        S <- NEW.NAMES[i]
        if(trace) { message("Fitting covariance model ",S) }
        NEW[[S]] <- meta.normal(MU,SIGMA,X=X,SX=SX,debias=debias,weights=weights,VARS=diag(TRYS[i,],DIM),MEANS=MEANS,GUESS=GUESS,WARN=FALSE,...)
      }
    }

    #################
    ### mean terms ###
    # which terms are presently off
    try <- which(MEANS & !MEAN.MIN & !VARS)
    if(length(try))
    {
      TRYS <- array(MEANS,c(DIM,length(try)))
      TRYS <- t(TRYS) # [off->on,all]
      colnames(TRYS) <- FEATURES
      for(i in 1:nrow(TRYS)) { TRYS[i,try[i]] <- FALSE }
      if("minor" %in% FEATURES) { TRYS[,"minor"] <- TRYS[,"angle"] <- TRYS[,"minor"] & TRYS[,"angle"] }
      TRYS <- unique(TRYS)
      # models that we will try by turning one term off
      NAMES <- NAMES2 <- make.names(TRYS,VARS)
      # paired down models that we haven't attempted yet
      SUB <- NAMES %nin% names(OLD)
      NEW.NAMES <- NEW2 <- NAMES[SUB]
      TRYS <- TRYS[SUB,,drop=FALSE]

      # fit all
      for(i in 1%:%length(NEW.NAMES))
      {
        S <- NEW.NAMES[i]
        if(trace) { message("Fitting covariance model ",S) }
        NEW[[S]] <- meta.normal(MU,SIGMA,X=X,SX=SX,debias=debias,weights=weights,VARS=diag(VARS,DIM),MEANS=TRYS[i,],GUESS=GUESS,WARN=FALSE,...)
      }
    }

    ###############
    ## wrap up ##
    OLD <- c(OLD,NEW)
    NEW <- OLD[c(NAMES1,NAMES2)] # also consider models that we've already fitted (with the right number of parameters)
    if(!length(c(NEW1,NEW2))) { break } # stopped trying new models

    ICS <- sapply(NEW,function(m){m[[IC]]})
    BEST <- which.min(ICS)
    BEST <- NEW[[BEST]]
    GUESS <- BEST
    MEANS <- BEST$INT
    VARS <- diag(BEST$VARS)
  } # finish pair down variances

  # take best variance model
  ICS <- sapply(OLD,function(m){m[[IC]]})
  I <- which.min(ICS)
  BEST <- OLD[[I]]
  NP <- num.pars(BEST,mean=TRUE)
  NV <- num.pars(BEST,mean=FALSE)

  # now see if correlations are supported
  if(NV>1 && sum(!INF)>=NP)
  {
    GUESS <- BEST
    MEANS <- BEST$INT
    VARS <- diag(BEST$VARS)
    S <- make.names(MEANS,VARS,OFF=TRUE)
    if(S %nin% names(OLD))
    {
      if(trace) { message("Fitting covariance model ",S) }
      OLD[[S]] <- meta.normal(MU,SIGMA,X=X,SX=SX,debias=debias,weights=weights,VARS=VARS,MEANS=MEANS,GUESS=GUESS,WARN=FALSE,...)
    }
  } # end correlation fit

  ICS <- sapply(OLD,function(m){m[[IC]]})
  names(ICS) <- names(OLD)
  I <- sort(ICS,index.return=TRUE)$ix
  ICS <- ICS[I]
  ICS <- ICS - ICS[1]
  I <- I[1]

  if(trace)
  {
    # report model selection
    ICS <- data.frame(ICS)
    colnames(ICS) <- paste0("\u0394",IC)
    message("* Model selection for autocovariance distribution.")
    print(ICS)
  }

  R <- OLD[[I]]
  MEANS <- R$INT
  VARS <- diag(R$VARS)
  R$name <- names(OLD)[I]
  R$isotropic <- !("minor" %in% FEATURES && MEANS["minor"]) # just the mean
  R$axes <- axes

  R$mu <- c(R$beta,R$mu)
  names(R)[ which(names(R)=="mu") ] <- "par" # population mean of features
  names(R)[ which(names(R)=="COV.mu") ] <- "COV" # uncertainty in population mean
  names(R)[ which(names(R)=="sigma") ] <- "POV" # population dispersion of features (under log)
  names(R)[ which(names(R)=="COV.sigma") ] <- "COV.POV" # uncertainty in population dispersion of features (under log)

  # transform from diagonal-ish basis for covariance model
  if(any(VARS))
  {
    JJ <- quad2lin(J)
    dimnames(JJ) <- dimnames(R$COV.POV)
    R$COV.POV <- JJ %*% R$COV.POV %*% t(JJ)

    EXTRA <- length(formula)
    if(EXTRA)
    {
      BASE <- diag(EXTRA+nrow(J))
      dimnames(BASE) <- dimnames(R$COV)
      SUB <- EXTRA+1:nrow(J)
      BASE[SUB,SUB] <- J
      J <- BASE
    }

    J <- solve(J)
    tJ <- t(J)
    R$par <- (R$par %*% tJ)[1,]
    R$COV <- J %*% R$COV %*% tJ
    R$POV <- J %*% R$POV %*% tJ
    # SUB <- FEATURES[variance]
    # J <- J[SUB,SUB]
  }

  # set beta structure (empty)
  BETA <- BETA[BETA %in% FEATURES]

  # fix names
  NEW <- names(R$beta)
  if(length(NEW))
  {
    # implement linear response
    for(i in which(grepl('/',NEW,fixed=TRUE)))
    {
      Q <- strsplit(NEW[i],'/',fixed=TRUE)[[1]]
      Q <- paste0(Q[1],":",Q[2])

      names(R$beta)[i] <- Q
      R$beta[i] <- R$beta[i] * 2
      # this comes from integration:
      # from beta(0) + beta'(0)*x          (slope meta-analysis)
      #   to beta(0)*x + 1/2*beta'(0)*x^2  (RSF)

      I <- which(names(R$par)==NEW[i])
      names(R$par)[i] <- Q
      R$par[i] <- R$par[i] * 2

      I <- which(rownames(R$COV)==NEW[i])
      rownames(R$COV)[I] <- colnames(R$COV)[I] <- Q
      R$COV[I,] <- R$COV[,I] <- R$COV[I,] * 2

      I <- rownames(R$POV)==NEW[i]
      rownames(R$POV)[I] <- colnames(R$POV)[I] <- Q
      R$POV[I,] <- R$POV[,I] <- R$POV[I,] * 2

      for(I in which(grepl(NEW[i],rownames(R$COV.POV),fixed=TRUE)))
      {
        S <- strsplit(rownames(R$COV.POV)[I],"-",fixed=TRUE)[[1]]

        if(S[1]==NEW[i])
        { S[1] <- Q }

        if(S[2]==NEW[i])
        { S[2] <- Q }

        S <- paste(S[1],"-",S[2])
        rownames(R$COV.POV)[I] <- colnames(R$COV.POV)[I] <- S
        R$COV.POV[I,] <-  R$COV.POV[,I] <- R$COV.POV[I,] * 4
      }
    } # end for NEW

    NEW <- names(R$beta) # new slope of slopes

    # add new population variances
    R$POV <- mpad(R$POV,diff=length(NEW),side=-1,padname=NEW)
    # add their uncertainties
    ADD <- outer(NEW,NEW,function(a,b){paste0(a,'-',b)})
    ADD <- ADD[ upper.tri(ADD,diag=TRUE) ] # unique terms
    ADD <- c( ADD , outer(NEW,FEATURES,function(a,b){paste0(a,'-',b)}) )
    R$COV.POV <- mpad(R$COV.POV,diff=length(ADD),side=-1,padname=ADD)

    BETA <- c(NEW,BETA)
    FEATURES <- c(NEW,FEATURES)

    # resort
    SORT <- outer(FEATURES,FEATURES,function(a,b){paste0(a,'-',b)})
    SORT <- SORT[ upper.tri(SORT,diag=TRUE) ]
    R$COV.POV <- R$COV.POV[SORT,SORT]
  } # end if(length(BETA))

  beta <- rep(0,length(BETA))
  names(beta) <- BETA
  R$beta <- beta

  if(length(BETA))
  { R$formula <- stats::as.formula(paste("~",paste(names(beta),collapse=" + "))) }

  R$features <- FEATURES

  # information for fitting that we no longer use
  NEW <- rep(TRUE,length(NEW))
  MEANS <- c(NEW,MEANS)
  VARS <- c(NEW,VARS)
  NIN <- !MEANS & !VARS
  if(any(NIN))
  {
    SUB <- !NIN
    NIN <- FEATURES[NIN]
    FEATURES <- FEATURES[SUB]
    R$par <- R$par[FEATURES]
    VARS <- VARS[SUB]
    # R$par <- NULL
    R$COV <- R$COV[FEATURES,FEATURES]
    R$POV <- R$POV[FEATURES,FEATURES]
    # COV.POV should be okay except when missing POV
    FEAT <- sapply(NIN,function(nin){grepl(nin,rownames(R$COV.POV),fixed=TRUE)})
    FEAT <- !apply(FEAT,1,any)
    R$COV.POV <- R$COV.POV[FEAT,FEAT]
  }

  # reverse log transformation
  # R <- set.parameters(R,par,linear.cov=TRUE)
  R <- exp_ctmm(R,debias=debias,variance=VARS)

  return(R)
}


##########
mean.mu <- function(x,debias=TRUE,weights=NULL,trace=FALSE,IC="AICc",...)
{
  if(is.null(weights))
  { weights <- rep(1,length(x)) }

  axes <- x[[1]]$axes
  AXES <- length(axes)
  N <- length(x)

  # Gaussian-Gaussian in all cases
  MU <- array(0,c(N,AXES))
  SIGMA <- array(0,c(N,AXES,AXES))

  colnames(MU) <- axes
  dimnames(SIGMA) <- list(NULL,axes,axes)

  for(i in 1:N)
  {
    MU[i,] <- x[[i]]$mu
    SIGMA[i,,] <- x[[i]]$COV.mu

    # TODO !!!
    # fill in with zeroes for non-stationary means
    # TODO !!!
  }

  # list of candidate models
  MM <- list()

  # no variance in mean location
  S <- "Dirac-\u03B4(\u03BC)"
  if(trace) { message("Fitting location-mean model ",S) }
  MM[[S]] <- meta.normal(MU,SIGMA,VARS=FALSE,isotropic=TRUE,debias=debias,weights=weights,WARN=FALSE,...)
  GUESS <- MM[[S]]

  # symmetric mean distribution
  if(2*N >= 2+1)
  {
    S <- "isotropic-\u03BC"
    if(trace) { message("Fitting location-mean model ",S) }
    MM[[S]] <- meta.normal(MU,SIGMA,isotropic=TRUE,debias=debias,weights=weights,GUESS=GUESS,WARN=FALSE,...)
    GUESS <- MM[[S]]

    # general distribution
    if(2*N >= 2+3)
    {
      S <- "anisotropic-\u03BC"
      if(trace) { message("Fitting location-mean model ",S) }
      MM[[S]] <- meta.normal(MU,SIGMA,debias=debias,weights=weights,GUESS=GUESS,WARN=FALSE,...)
    }
  }

  ICS <- sapply(MM,function(m){m[[IC]]})
  names(ICS) <- names(MM)
  i <- which.min(ICS)
  i <- names(MM)[i]
  MM <- MM[[i]]
  MM$name <- i

  # report model selection
  if(trace)
  {
    i <- sort(ICS,index.return=TRUE)$ix
    ICS <- ICS[i] # sorted
    ICS <- ICS - ICS[1] # start at zero
    ICS <- data.frame(ICS)
    colnames(ICS) <- paste0("\u0394",IC)
    message("* Model selection for location-mean \u03BC distribution.")
    print(ICS)
  }

  # R$mu # population mean of mean locations
  # R$COV.mu # uncertainty in mean of means estimate
  names(MM)[ which(names(MM)=="sigma") ] <- "POV.mu" # dispersion of means
  names(MM)[ which(names(MM)=="COV.sigma") ] <- "COV.POV.mu" # uncertainty in dispersion of means

  return(MM)
}


###########
mean.ctmm <- function(x,formula=FALSE,weights=NULL,sample=TRUE,debias=TRUE,IC="AIC",trace=TRUE,...)
{
  select <- "all"

  # if(length(x)==1)
  # {
  #   warning("Only one individual in list.")
  #   return(x)
  # }

  n <- length(x)
  info <- mean_info(x)
  axes <- x[[1]]$axes

  if(is.null(weights))
  { weights <- rep(1,n) }
  else
  { weights <- weights/max(weights) }
  names(weights) <- names(x)

  isotropic <- c(TRUE,TRUE)
  names(isotropic) <- c("sigma","mu")

  #######################
  # average of mean locations
  if(sample)
  {
    MM <- mean.mu(x,debias=debias,weights=weights,IC=IC,trace=trace,...)
    isotropic['mu'] <- MM$isotropic

    FM <- mean.features(x,debias=debias,weights=weights,select=select,IC=IC,trace=trace,formula=formula,...)
    isotropic['sigma'] <- FM$isotropic

    x <- copy(from=FM,to=MM)
    x$isotropic <- isotropic

    # STORE EVERYTHING
    x$name <- paste(MM$name,FM$name)
    x$AIC <- MM$AIC + FM$AIC
    x$AICc <- MM$AICc + FM$AICc
    x$BIC <- MM$BIC + FM$BIC
  } # end sample
  else # !sample
  {
    # COV[mu^mu] diagonal-term
    cov.diag <- function(cov)
    {
      COV <- cov^2
      VEC <- c( diag(cov)*cov[1,2] , cov[1,1]*cov[2,2]+cov[1,2]^2 )
      COV <- rbind(COV,VEC[1:2])
      COV <- cbind(COV,VEC)
      COV <- 2*COV
      return(COV)
    }

    cov.off <- function(cov1,cov2)
    {
      COV <- 4 * cov1 * cov2
      VEC <- c(2*cov1[1,1]*cov2[1,2]+2*cov1[1,2]*cov2[1,1],2*cov1[2,2]*cov2[1,2]+2*cov1[1,2]*cov2[2,2],cov1[1,1]*cov2[2,2]+cov1[2,2]*cov2[1,1]+2*cov1[1,2]*cov2[1,2])
      COV <- rbind(COV,VEC[1:2])
      COV <- cbind(COV,VEC)
      return(COV)
    }

    MU <- COV.MU <- COV <- COV.COV <- 0
    w <- weights/sum(weights)
    for(i in 1:n)
    {
      MU <- MU + w[i]*x[[i]]$mu
      COV.MU <- COV.MU + w[i]^2*x[[i]]$COV.mu
      COV <- COV + w[i]*( x[[i]]$sigma + outer(c(x[[i]]$mu)) ) # M2
      COV.COV <- COV.COV + w[i]^2*sigma.COV(x[[i]]) + (w[i]-w[i]^2)^2*cov.diag(x[[i]]$COV.mu)
      for(j in (i+1)%:%n) { COV.COV <- COV.COV - w[i]^2*w[j]^2*cov.off(x[[i]]$COV.mu,x[[j]]$COV.mu) }
    }
    COV <- COV - outer(c(MU)) # M2 - M1^2
    rownames(COV.MU) <- colnames(COV.MU) <- axes

    isotropic['mu'] <- FALSE
    if(!all(sapply(x,function(y){y$isotropic[1]}))) { isotropic['sigma'] <- FALSE }

    sigma <- covm(COV,axes=axes,isotropic=isotropic['sigma'])
    if(isotropic['sigma'])
    {
      COV <- COV.COV[1,1,drop=FALSE]
      rownames(COV) <- colnames(COV) <- 'major'
      features <- 'major'
    }
    else
    {
      J <- J.par.sigma(sigma[c(1,4,2)])
      COV <- J %*% COV.COV %*% t(J)
      features <- c('major','minor','angle')
    }

    x <- ctmm(mu=MU,COV.mu=COV.MU,sigma=sigma,COV=COV,features=features)
    # TODO non-sigma features
  } # end !sample

  x$info <- info
  x <- do.call(ctmm,x)

  x$isotropic <- isotropic
  x$weights <- weights

  # return final result
  # x <- new.ctmm(x,info)
  return(x)
}


## population stationary distribution
mean_pop <- function(CTMM)
{
  # spread (x,y)
  sigma <- CTMM$POV.mu + CTMM$sigma
  # uncertainty of spread (x-x,x-y,y-y)
  if(CTMM$isotropic['mu'])
  { COV <- c(CTMM$COV.POV.mu) * diag(c(1,1,0)) }
  else
  {
    P <- c(1,3,2) # xx,yy,xy
    COV <- CTMM$COV.POV.mu[P,P]
  }

  if(CTMM$isotropic['sigma'])
  {
    P <- 'major'
    COV <- COV + c(CTMM$COV[P,P]+CTMM$POV[P,P])*diag(c(1,1,0))
  }
  else
  {
    P <- c('major','minor','angle')
    J <- J.sigma.par(CTMM$sigma)
    COV <- COV + J%*%(CTMM$COV[P,P]+CTMM$POV[P,P])%*%t(J)
  }

  isotropic <- all(CTMM$isotropic)
  sigma <- covm(sigma,isotropic=isotropic)

  if(isotropic)
  {
    P <- 'major'
    COV <- COV[1,1,drop=FALSE]
  }
  else
  {
    P <- c('major','minor','angle')
    J <- J.par.sigma(sigma)
    COV <- J %*% COV %*% t(J)
  }
  dimnames(COV) <- list(P,P)

  CTMM$isotropic <- isotropic
  CTMM$sigma <- sigma
  CTMM$COV <- COV
  CTMM$POV <- CTMM$COV.POV <- CTMM$COV.POV.mu <- NULL
  CTMM$features <- rownames(COV)
  CTMM$tau <- NULL
  CTMM$circle <- FALSE

  # BETA ARE TOTALLY UNFINISHED
  # will need eventually...

  return(CTMM)
}

