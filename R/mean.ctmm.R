###########
# log-transformed parameters
# debias includes bias correction for chi^2 to log(chi^2)
# matrix casts location covariance, diffusion rate, velocity covariance all as distinct matrices for above bias correction
log.ctmm <- function(CTMM,debias=FALSE,...)
{
  if(class(CTMM)[1]=="list") { return(log.ctmms(CTMM,debias=debias,...)) }

  SIGMA <- c("major","minor","angle")
  isotropic <- CTMM$isotropic
  axes <- CTMM$axes
  features <- CTMM$features
  par <- get.parameters(CTMM,features,linear.cov=FALSE)
  sigma <- CTMM$sigma
  COV <- CTMM$COV

  ### log transform sigma
  # log transform major-minor in sigma and COV
  PARS <- SIGMA[SIGMA %in% features]
  PARS <- PARS[PARS %in% POSITIVE.PARAMETERS]

  COV[PARS,] <- COV[PARS,,drop=FALSE] / par[PARS]
  COV[,PARS] <- t( t(COV[,PARS,drop=FALSE]) / par[PARS] )

  sigma <- log.covm(sigma)

  # convert log(eigen) to log(xy) in COV
  PARS <- SIGMA[SIGMA %in% features]
  if(isotropic)
  { par[PARS] <- sigma@par[PARS] }
  if(!isotropic)
  {
    par[PARS] <- sigma[c(1,4,2)]  # log 'xx', 'yy', 'xy'

    J <- J.sigma.par(sigma@par)
    COV[PARS,] <- J %*% COV[PARS,,drop=FALSE]
    COV[,PARS] <- COV[,PARS,drop=FALSE] %*% t(J)
  }

  ### log transform all positive parameters that haven't been logged
  # features to log transform
  FEAT.ALL <- features[ features %in% POSITIVE.PARAMETERS | grepl("error",features) ]
  FEAT <- FEAT.ALL[FEAT.ALL %nin% SIGMA]

  if(length(FEAT))
  {
    COV[FEAT,] <- COV[FEAT,,drop=FALSE] / par[FEAT]
    COV[,FEAT] <- t( t(COV[,FEAT,drop=FALSE]) / par[FEAT] )
    par[FEAT] <- log(par[FEAT])
  }

  # log chi^2 bias correction
  VAR <- diag(COV) # don't include features that shouldn't be there
  SUB <- features[features %in% c(FEAT.ALL,"angle") & VAR[features]<Inf]
  if(debias && length(SUB))
  {
    COR <- diag(1,nrow(COV))
    dimnames(COR) <- dimnames(COV)
    GOOD <- VAR>.Machine$double.eps
    COR[GOOD,GOOD] <- stats::cov2cor(COV[GOOD,GOOD])
    # COR <- nant(COR,0)

    # diagonalize and log-chi^2 debias relevant parameter estimates
    EIGEN <- eigen(COV[SUB,SUB])
    dimnames(EIGEN$vectors) <- list(SUB,SUB)
    names(EIGEN$values) <- SUB

    # fix signs
    if(isotropic) { PARS <- "major" } else { PARS <- c("major","minor") }
    # VAR goes in log numerator for chi^2 variates: variance, diffusion, MS speed, ...
    for(i in 1:ncol(EIGEN$vectors)) { if(sum(EIGEN$vectors[PARS,i])<0) { EIGEN$vectors[,i] <- -EIGEN$vectors[,i] } }

    # transform to diagonalized basis with VARs in log numerator
    par[SUB] <- t(EIGEN$vectors) %*% par[SUB] # diagonalize parameters
    DOF <- 2/EIGEN$values # log-chi^2 VAR-DOF relation
    BIAS <- digamma(DOF/2)-log(DOF/2) # negative bias for log(chi^2) variates
    BIAS <- nant(BIAS,0)
    BIAS <- pmax(BIAS,digamma(1/2)-log(1/2)) # clamp to 1 DOF
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

  RETURN <- list(par=par,COV=COV,isotropic=isotropic)
  return(RETURN)
}


#####################
# inverse transformation of above
exp.ctmm <- function(CTMM,debias=FALSE,variance=TRUE)
{
  SIGMA <- c("major","minor","angle")
  isotropic <- CTMM$isotropic
  axes <- CTMM$axes
  features <- CTMM$features
  variance <- array(variance,length(features))
  names(variance) <- features
  if("par" %in% names(CTMM))
  { par <- CTMM$par }
  else
  { par <- get.parameters(CTMM,features,linear.cov=TRUE) }
  COV <- CTMM$COV
  POV <- CTMM$POV
  COV.POV <- CTMM$COV.POV
  JP <- diag(1,nrow(POV))
  dimnames(JP) <- dimnames(COV)

  # log chi^2 bias correction
  FEAT.ALL <- features[ features %in% POSITIVE.PARAMETERS | grepl("error",features) ]
  SUB <- features[features %in% c(FEAT.ALL,"angle")]
  if(debias && length(SUB))
  {
    ##################
    ### diagonalize and log-chi^2 debias point estimates
    EIGEN <- (COV+POV)[SUB,SUB]
    EIGEN <- eigen(EIGEN)
    dimnames(EIGEN$vectors) <- list(SUB,SUB)
    names(EIGEN$values) <- SUB
    EIGEN$values <- clamp(EIGEN$values,0,Inf)

    # fix signs
    if(isotropic) { PARS <- "major" } else { PARS <- c("major","minor") }
    # VAR goes in log numerator for chi^2 variates: variance, diffusion, MS speed, ...
    for(i in 1:nrow(EIGEN$vectors)) { if(sum(EIGEN$vectors[PARS,i])<0) { EIGEN$vectors[,i] <- -EIGEN$vectors[,i] } }

    # transform to diagonalized basis with VARs in log numerator
    par[SUB] <- t(EIGEN$vectors) %*% par[SUB] # diagonalize parameters

    # log-gamma variance (better than delta method)
    DOF <- 2*itrigamma(EIGEN$values)
    BIAS <- digamma(DOF/2)-log(DOF/2) # negative bias for log(chi^2) variates
    BIAS <- nant(BIAS,0)
    BIAS <- pmax(BIAS,digamma(1/2)-log(1/2)) # clamp to 1 DOF
    par[SUB] <- par[SUB] + BIAS # E[log-chi^2] bias correction
    par[SUB] <- c(EIGEN$vectors %*% par[SUB]) # transform back (still under logarithm)

    ################
    ### diagonalize and log-chi^2 debias COV
    VAR <- diag(COV)
    COR <- diag(1,nrow(COV))
    dimnames(COR) <- dimnames(COV)
    GOOD <- VAR>.Machine$double.eps
    COR[GOOD,GOOD] <- stats::cov2cor(COV[GOOD,GOOD])

    EIGEN <- COV[SUB,SUB]
    EIGEN <- eigen(EIGEN)
    dimnames(EIGEN$vectors) <- list(SUB,SUB)
    names(EIGEN$values) <- SUB
    EIGEN$values <- clamp(EIGEN$values,0,Inf)

    # fix signs
    # VAR goes in log numerator for chi^2 variates: variance, diffusion, MS speed, ...
    for(i in 1:nrow(EIGEN$vectors)) { if(sum(EIGEN$vectors[PARS,i])<0) { EIGEN$vectors[,i] <- -EIGEN$vectors[,i] } }

    # log-gamma variance (better than delta method)
    DOF <- 2*itrigamma(EIGEN$values)
    EIGEN$values <- 2/DOF
    EIGEN <- EIGEN$vectors %*% (EIGEN$values * t(EIGEN$vectors))

    COR[SUB,SUB] <- stats::cov2cor(EIGEN)
    VAR[SUB] <- diag(EIGEN)
    COV <- COR * outer(sqrt(VAR))

    #################
    ### diagonalize and log-chi^2 debias POV
    VAR <- diag(POV)
    if(all(VAR>.Machine$double.eps))
    {
      COR <- stats::cov2cor(POV)

      EIGEN <- POV[SUB,SUB]
      EIGEN <- eigen(EIGEN)
      dimnames(EIGEN$vectors) <- list(SUB,SUB)
      names(EIGEN$values) <- SUB
      EIGEN$values <- clamp(EIGEN$values,0,Inf)

      # fix signs
      # VAR goes in log numerator for chi^2 variates: variance, diffusion, MS speed, ...
      for(i in 1:nrow(EIGEN$vectors)) { if(sum(EIGEN$vectors[PARS,i])<0) { EIGEN$vectors[,i] <- -EIGEN$vectors[,i] } }

      # log-gamma variance (better than delta method)
      SCALE <- EIGEN$values

      DOF <- 2*itrigamma(EIGEN$values)
      EIGEN$values <- 2/DOF
      POV.SUB <- EIGEN$vectors %*% (EIGEN$values * t(EIGEN$vectors))
      COR[SUB,SUB] <- stats::cov2cor(POV.SUB)
      VAR[SUB] <- diag(POV.SUB)
      POV <- COR * outer(sqrt(VAR))

      SCALE <- sqrt(EIGEN$values/SCALE)
      SCALE <- nant(SCALE,1)
      JP[SUB,SUB] <- EIGEN$vectors %*% diag(SCALE) %*% t(EIGEN$vectors)
    }
  }

  ### exp transform all positive parameters except sigma
  # features to log transform
  FEAT <- features[features %in% POSITIVE.PARAMETERS]
  FEAT <- FEAT[FEAT %nin% SIGMA]

  if(length(FEAT))
  {
    par[FEAT] <- exp(par[FEAT])

    COV[FEAT,] <- COV[FEAT,,drop=FALSE] * par[FEAT]
    COV[,FEAT] <- t( t(COV[,FEAT,drop=FALSE]) * par[FEAT] )

    POV[FEAT,] <- POV[FEAT,,drop=FALSE] * par[FEAT]
    POV[,FEAT] <- t( t(POV[,FEAT,drop=FALSE]) * par[FEAT] )

    JP[FEAT,] <- JP[FEAT,] * par[FEAT]
  }

  ### exp transform sigma
  PARS <- SIGMA[SIGMA %in% features]
  sigma <- par[PARS]

  # convert COV[log(xy)] to COV[log(eigen)] (and POV)
  if(!isotropic)
  {
    J <- J.par.sigma(sigma)

    COV[PARS,] <- J %*% COV[PARS,,drop=FALSE]
    COV[,PARS] <- COV[,PARS,drop=FALSE] %*% t(J)

    POV[PARS,] <- J %*% POV[PARS,,drop=FALSE]
    POV[,PARS] <- POV[,PARS,drop=FALSE] %*% t(J)

    JP[PARS,] <- J %*% JP[PARS,,drop=FALSE]

    sigma <- matrix(sigma[c('major','angle','angle','minor')],2,2)
  }
  # else sigma is just 'major'

  # exp of eigen
  sigma <- covm(sigma,axes=axes,isotropic=isotropic)
  sigma <- exp.covm(sigma)

  # copy over
  par[PARS] <- sigma@par[PARS]

  PARS <- PARS[PARS %in% POSITIVE.PARAMETERS]

  COV[PARS,] <- COV[PARS,,drop=FALSE] * par[PARS]
  COV[,PARS] <- t( t(COV[,PARS,drop=FALSE]) * par[PARS] )

  POV[PARS,] <- POV[PARS,,drop=FALSE] * par[PARS]
  POV[,PARS] <- t( t(POV[,PARS,drop=FALSE]) * par[PARS] )

  JP[PARS,] <- JP[PARS,,drop=FALSE] * par[PARS]

  # finally transform COV.POV
  if(any(variance))
  {
    JP <- quad2lin(JP)
    NAMES <- dimnames(COV.POV)
    COV.POV <- JP %*% COV.POV %*% t(JP)
    dimnames(COV.POV) <- NAMES
  }
  else
  { COV.POV <- NULL }

  CTMM$sigma <- sigma
  CTMM$COV <- COV
  CTMM$POV <- POV
  CTMM$COV.POV <- COV.POV

  CTMM$par <- NULL
  CTMM <- set.parameters(CTMM,par)

  return(CTMM)
}

# M %*% X %*% t(M) -> L %*% X[TRI]
quad2lin <- function(M,diag=FALSE)
{
  M <- rbind(M)
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
    E <- array(0,dim(M))
    E[TRI[i]] <- 1
    E <- E + t(E)
    E <- E/max(E)
    E <- M %*% E %*% t(M) # [n,n]
    if(diag) # n!=m and only care about n variances
    { E <- diag(E) }
    else # n==m and want full covariance
    { E <- E[TRI] }
    J[,i] <- E
  }

  return(J)
}

#############
mean.features <- function(x,debias=TRUE,weights=NULL,trace=FALSE,IC="AICc",...)
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

  for(i in 1:N)
  {
    x[[i]] <- log.ctmm(x[[i]],debias=debias)
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
  }

  if("tau" %in% FEATURES)
  { J["tau",c("major","tau")] <- c(1,-1) } # D

  if("tau position" %in% FEATURES)
  { J["tau position",c("major","tau position")] <- c(1,-1) } # D

  # if("tau velocity" %in% FEATURES)
  # { J["tau velocity",c("major","tau position","tau velocity")] <- c(1,-1,-1) } # MSV correlates with D

  if("omega" %in% FEATURES) # MSV - delta-MSV # sigma cancels out
  {
    if("tau" %in% FEATURES) # OUf
    { J["omega",c("tau","omega")] <- c(1,-1) }
    else if(all(c("tau position","tau velocity") %in% FEATURES)) # OUF
    { J["omega",c("tau position","tau velocity","omega")] <- c(1,1,-2) }
    else if("tau velocity" %in% FEATURES) # IOU
    { J["omega",c("tau velocity","omega")] <- c(1,-1) }
    else if("tau position" %in% FEATURES) # OU
    { J["omega",c("tau position","omega")] <- c(1,-1) }
  }

  tJ <- t(J)
  MU <- MU %*% tJ

  # missing data / infinite uncertainties
  INF <- t(apply(SIGMA,1,function(M){diag(M)==Inf})) # [ind,par]

  SIGMA <- SIGMA %.% tJ
  for(i in 1:nrow(SIGMA)) { for(j in 1:ncol(SIGMA)) { if(INF[i,j]) { SIGMA[i,j,] <- SIGMA[i,,j] <- 0; SIGMA[i,j,j] <- Inf } } }
  SIGMA <- aperm(SIGMA,c(1,3,2))
  SIGMA <- SIGMA %.% tJ
  for(i in 1:nrow(SIGMA)) { for(j in 1:ncol(SIGMA)) { if(INF[i,j]) { SIGMA[i,j,] <- SIGMA[i,,j] <- 0; SIGMA[i,j,j] <- Inf } } }

  # number of estimated features (amount of data)
  DIM <- length(FEATURES)

  MEANS <- list()
  MEANS[["isotropic"]] <- array(TRUE,DIM)
  names(MEANS[["isotropic"]]) <- FEATURES
  if(!ISO)
  {
    MEANS[["anisotropic"]] <- MEANS[["isotropic"]]
    MEANS[["anisotropic"]][c('minor','angle')] <- FALSE
  }

  NAME.ZERO <- "Dirac-\u03B4(\u03A3)"
  NAME.ALL <- "diag(\u03A3)"
  make.names <- function(TRYS)
  {
    NAMES <- sapply(1:nrow(TRYS),function(i){paste0(FEATURES[TRYS[i,]],collapse=",")})
    for(i in 1:length(NAMES))
    {
      if(all(!TRYS[i,])) # Dirac-delta
      { NAMES[i] <- NAME.ZERO }
      else if(DIM==1 && TRYS[i,1]) # full matrix
      { NAMES[i] <- "\u03A3" }
      else if(all(TRYS[i,])) # full diagonal
      { NAMES[i] <- NAME.ALL }
      else if(!grepl(",",NAMES[i])) # part diagonal, but all correlations
      { NAMES[i] <- paste0("\u03A3[",NAMES[i],"]") }
      else # part diagonal, no correlations
      { NAMES[i] <- paste0("diag(\u03A3)[",NAMES[i],"]") }
    }

    return(NAMES)
  }

  FM <- list()
  for(m in names(MEANS))
  {
    ms <- paste0(m,"-\u03C3")
    VAR.MAX <- colSums(!INF)>1 | !MEANS[[m]]
    OLD <- list()

    ##########
    ## build up variance phase
    ## no variance
    S <- NAME.ZERO
    if(trace) { message("Fitting autocovariance model ",S," ",ms) }
    OLD[[S]] <- meta.normal(MU,SIGMA,debias=debias,weights=weights,VARS=FALSE,MEANS=MEANS[[m]],...)

    ## add variances one by one
    VARS <- array(FALSE,DIM)
    GUESS <- OLD[[1]]
    while(any(!VARS))
    {
      # which terms are presently off
      try <- which(!VARS)
      TRYS <- array(VARS,c(DIM,length(try)))
      TRYS <- t(TRYS) # [off->on,all]
      colnames(TRYS) <- FEATURES
      for(i in 1:nrow(TRYS))
      {
        TRYS[i,try[i]] <- TRUE
        TRYS[i,] <- TRYS[i,] & VAR.MAX
      }
      if("minor" %in% FEATURES) { TRYS[,"minor"] <- TRYS[,"angle"] <- TRYS[,"minor"] | TRYS[,"angle"] }
      TRYS <- unique(TRYS)

      # models that we will try by turning one term on
      NAMES <- make.names(TRYS)
      # don't try what's been done
      SUB <- NAMES %nin% names(OLD)
      SUB <- which(SUB)

      if(!length(SUB)) { break } # no new models to fit

      TRYS <- TRYS[SUB,,drop=FALSE]
      NAMES <- NAMES[SUB]

      # fit all
      NEW <- list()
      for(i in 1%:%length(NAMES))
      {
        S <- NAMES[i]
        if(trace) { message("Fitting autocovariance model ",S," ",ms) }
        VARS <- diag(TRYS[i,],length(FEATURES))
        NEW[[S]] <- meta.normal(MU,SIGMA,debias=debias,weights=weights,VARS=VARS,MEANS=MEANS[[m]],GUESS=GUESS,...)
      }

      ICS <- sapply(NEW,function(m){m[[IC]]})
      BEST <- which.min(ICS)
      BEST <- NEW[[BEST]]
      GUESS <- BEST
      VARS <- diag(BEST$VARS)

      OLD <- c(OLD,NEW)

      if(BEST[[IC]]==Inf) { break }
    } # finish build up variances

    ###########
    ## pair down variance phase
    # start with all variances or best (if IC==Inf)
    if(NAME.ALL %in% names(OLD))
    { BEST <- OLD[[NAME.ALL]] }
    else # limited data cannot support all variances
    {
      ICS <- sapply(OLD,function(m){m[[IC]]})
      BEST <- which.min(ICS)
      BEST <- OLD[[BEST]]
    }
    GUESS <- BEST
    VARS <- diag(BEST$VARS)

    ##############
    # pair down variances phase
    while(sum(VARS)>1)
    {
      # which terms are presently off
      try <- which(VARS)
      TRYS <- array(VARS,c(DIM,length(try)))
      TRYS <- t(TRYS) # [off->on,all]
      colnames(TRYS) <- FEATURES
      for(i in 1:nrow(TRYS)) { TRYS[i,try[i]] <- FALSE }
      if("minor" %in% FEATURES)
      {
        TRYS[,"minor"] <- TRYS[,"angle"] <- TRYS[,"minor"] & TRYS[,"angle"]
        TRYS <- unique(TRYS)
      }
      # models that we will try by turning one term off
      NAMES <- make.names(TRYS)
      # paired down models that we haven't attempted yet
      NEW.NAMES <- NAMES[NAMES %nin% names(OLD)]

      # fit all
      NEW <- list()
      for(i in 1%:%length(NEW.NAMES))
      {
        S <- NEW.NAMES[i]
        if(trace) { message("Fitting autocovariance model ",S," ",ms) }
        VARS <- diag(TRYS[i,],length(FEATURES))
        NEW[[S]] <- meta.normal(MU,SIGMA,debias=debias,weights=weights,VARS=VARS,MEANS=MEANS[[m]],GUESS=GUESS,...)
      }

      OLD <- c(OLD,NEW)
      NEW <- OLD[NAMES] # also consider models that we've already fitted (with the right number of parameters)

      ICS <- sapply(NEW,function(m){m[[IC]]})
      BEST <- which.min(ICS)
      BEST <- NEW[[BEST]]
      GUESS <- BEST
      VARS <- diag(BEST$VARS)
    } # finish pair down variances

    # take best variance model
    ICS <- sapply(OLD,function(m){m[[IC]]})
    BEST <- which.min(ICS)
    BEST <- OLD[[BEST]]
    GUESS <- BEST
    VARS <- diag(BEST$VARS)

    # now see if correlations are supported
    V <- sum(VARS)
    if(V>1 && sum(!INF)-sum(MEANS[[m]]) >= V*(V+1)/2)
    {
      if(V<DIM)
      {
        S <- paste0(FEATURES[VARS],collapse=",")
        S <- paste0("\u03A3[",S,"]")
      }
      else
      { S <- "\u03A3" }

      if(S %nin% names(OLD))
      {
        if(trace) { message("Fitting autocovariance model ",S," ",ms) }
        OLD[[S]] <- meta.normal(MU,SIGMA,debias=debias,weights=weights,VARS=VARS,MEANS=MEANS[[m]],GUESS=GUESS,...)
      }
    } # end correlation fit

    names(OLD) <- paste(names(OLD),ms)
    FM <- c(FM,OLD)
  } # for(m in names(MEANS))

  ICS <- sapply(FM,function(m){m[[IC]]})
  names(ICS) <- names(FM)
  NAME <- which.min(ICS)
  NAME <- names(ICS)[NAME]
  R <- FM[[NAME]]
  R$name <- NAME
  R$isotropic <- !grepl("anisotropic",R$name)
  VARS <- diag(R$VARS)

  if(trace)
  {
    # report model selection
    i <- sort(ICS,index.return=TRUE)$ix
    ICS <- ICS[i] # sorted
    ICS <- ICS - ICS[1] # start at zero
    ICS <- data.frame(ICS)
    colnames(ICS) <- paste0("\u0394",IC)
    message("* Model selection for autocovariance distribution.")
    print(ICS)
  }

  R$axes <- axes
  names(R)[ which(names(R)=="mu") ] <- "par" # population mean of features
  names(R)[ which(names(R)=="COV.mu") ] <- "COV" # uncertainty in population mean
  names(R)[ which(names(R)=="sigma") ] <- "POV" # population dispersion of features (under log)
  names(R)[ which(names(R)=="COV.sigma") ] <- "COV.POV" # uncertainty in population dispersion of features (under log)

  # transform from diagonal-ish basis for covariance model
  if(any(VARS))
  {
    J <- solve(J)
    tJ <- t(J)
    R$par <- (R$par %*% tJ)[1,]
    R$COV <- J %*% R$COV %*% tJ
    R$POV <- J %*% R$POV %*% tJ
    # SUB <- FEATURES[variance]
    # J <- J[SUB,SUB]
    J <- quad2lin(J)
    dimnames(J) <- dimnames(R$COV.POV)
    R$COV.POV <- J %*% R$COV.POV %*% t(J)
  }

  # information for fitting that we no longer use
  if(R$isotropic && "minor" %in% FEATURES)
  {
    SUB <- FEATURES %nin% c('minor','angle')
    VARS <- VARS[SUB]
    FEATURES <- FEATURES[SUB]
    R$par <- R$par[FEATURES]
    # R$par <- NULL
    R$COV <- R$COV[FEATURES,FEATURES]
    R$POV <- R$POV[FEATURES,FEATURES]
    # COV.POV should be okay except when missing POV
    FEAT <- grepl('minor',rownames(R$COV.POV),fixed=TRUE) | grepl('angle',rownames(R$COV.POV),fixed=TRUE)
    FEAT <- !FEAT
    R$COV.POV <- R$COV.POV[FEAT,FEAT]
  }
  R$features <- FEATURES

  # set beta structure (empty)
  beta <- array(0,length(BETA))
  names(beta) <- BETA
  if(length(BETA)) { R$beta <- beta }

  # reverse log transformation
  # R <- set.parameters(R,par,linear.cov=TRUE)
  R <- exp.ctmm(R,debias=debias,variance=VARS)

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
  MM[[S]] <- meta.normal(MU,SIGMA,VARS=FALSE,isotropic=TRUE,debias=debias,weights=weights,...)

  # symmetric mean distribution
  if(2*N >= 2+1)
  {
    S <- "isotropic-\u03BC"
    if(trace) { message("Fitting location-mean model ",S) }
    MM[[S]] <- meta.normal(MU,SIGMA,isotropic=TRUE,debias=debias,weights=weights,...)
    GUESS <- MM[[S]]

    # general distribution
    if(2*N >= 2+3)
    {
      S <- "anisotropic-\u03BC"
      if(trace) { message("Fitting location-mean model ",S) }
      MM[[S]] <- meta.normal(MU,SIGMA,debias=debias,weights=weights,GUESS=GUESS,...)
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
mean.ctmm <- function(x,weights=NULL,sample=TRUE,debias=TRUE,IC="AICc",trace=TRUE,...)
{
  # if(length(x)==1)
  # {
  #   warning("Only one individual in list.")
  #   return(x)
  # }

  n <- length(x)
  info <- mean.info(x)
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

    FM <- mean.features(x,debias=debias,weights=weights,IC=IC,trace=trace,...)
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
mean.pop <- function(CTMM)
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
    J <- J.sigma.par(CTMM$sigma@par)
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
    J <- J.par.sigma(sigma[c(1,4,2)])
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

