###########
# log-transformed parameters
# debias includes bias correction for chi^2 to log(chi^2)
# matrix casts location covariance, diffusion rate, velocity covariance all as distinct matrices for above bias correction
log.ctmm <- function(CTMM,debias=FALSE)
{
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
    COR <- stats::cov2cor(COV)
    COR <- nant(COR,0)

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
exp.ctmm <- function(CTMM,debias=FALSE)
{
  SIGMA <- c("major","minor","angle")
  isotropic <- CTMM$isotropic
  axes <- CTMM$axes
  features <- CTMM$features
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
    par[SUB] <- par[SUB] + BIAS # E[log-chi^2] bias correction
    par[SUB] <- c(EIGEN$vectors %*% par[SUB]) # transform back (still under logarithm)

    ################
    ### diagonalize and log-chi^2 debias COV
    VAR <- diag(COV)
    COR <- stats::cov2cor(COV)

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
    if(all(VAR>0))
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
  JP <- quad2lin(JP)
  NAMES <- dimnames(COV.POV)
  COV.POV <- JP %*% COV.POV %*% t(JP)
  dimnames(COV.POV) <- NAMES

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
mean.features <- function(x,debias=TRUE,isotropic=FALSE,variance=TRUE,weights=NULL,...)
{
  if(is.null(weights))
  { weights <- rep(1,length(x)) }

  N <- length(x)
  axes <- x[[1]]$axes

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

  MEANS <- VARS <- rep(TRUE,M)
  names(MEANS) <- names(VARS) <- FEATURES

  if(isotropic && "minor" %in% FEATURES)
  { MEANS[c("minor","angle")] <- VARS[c("minor","angle")] <- c(FALSE,FALSE) }

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
    if(!isotropic && x[[i]]$isotropic)
    {
      MU[i,'minor'] <- MU[i,'major']
      # perfectly correlated observation
      SIGMA[i,P,P] <- matrix(SIGMA[i,'major','major'],2,2)
      # copy over other correlations as well
      SIGMA[i,'minor',NFEAT] <- SIGMA[i,'major',NFEAT]
      SIGMA[i,NFEAT,'minor'] <- SIGMA[i,NFEAT,'major']
    }
    else if(isotropic && !x[[i]]$isotropic)
    {
      MU[i,P] <- c(mean(MU[i,P]),-diff(MU[i,P]))
      J <- rbind( c(1/2,1/2) , c(1,-1) )
      SIGMA[i,P,] <- J %*% SIGMA[i,P,]
      SIGMA[i,,P] <- SIGMA[i,,P] %*% t(J)
    }
  } # for(i in 1:N)

  if(!variance) { VARS <- FALSE }
  R <- meta.normal(MU,SIGMA,MEANS=MEANS,VARS=VARS,debias=debias,weights=weights)
  R$isotropic <- isotropic
  R$axes <- axes
  names(R)[ which(names(R)=="mu") ] <- "par" # population mean of features
  names(R)[ which(names(R)=="COV.mu") ] <- "COV" # uncertainty in population mean
  names(R)[ which(names(R)=="sigma") ] <- "POV" # population dispersion of features (under log)
  names(R)[ which(names(R)=="COV.sigma") ] <- "COV.POV" # uncertainty in population dispersion of features (under log)

  # information for fitting that we no longer use
  if(isotropic && "minor" %in% FEATURES)
  {
    FEATURES <- FEATURES[FEATURES %nin% c('minor','angle')]
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

  # set beta structure
  beta <- rep(0,length(BETA))
  names(beta) <- BETA
  if(length(BETA)) { R$beta <- beta }

  # reverse log transformation
  # R <- set.parameters(R,par,linear.cov=TRUE)
  R <- exp.ctmm(R,debias=debias)

  return(R)
}


##########
mean.mu <- function(x,debias=TRUE,isotropic=FALSE,variance=TRUE,weights=NULL,...)
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
  R <- meta.normal(MU,SIGMA,isotropic=isotropic,VARS=variance,debias=debias,weights=weights)
  # R$mu # population mean of mean locations
  # R$COV.mu # uncertainty in mean of means estimate
  names(R)[ which(names(R)=="sigma") ] <- "POV.mu" # dispersion of means
  names(R)[ which(names(R)=="COV.sigma") ] <- "COV.POV.mu" # uncertainty in dispersion of means

  return(R)
}


###########
mean.ctmm <- function(x,weights=NULL,debias=TRUE,IC="AICc",...)
{
  # if(length(x)==1)
  # {
  #   warning("Only one individual in list.")
  #   return(x)
  # }

  if(is.null(weights))
  { weights <- rep(1,length(x)) }
  else
  { weights <- weights/max(weights) }
  names(weights) <- names(x)

  isotropic <- c(TRUE,TRUE)
  names(isotropic) <- c("sigma","mu")

  #######################
  # average of mean locations
  MU <- mean.mu(x,debias=debias,isotropic=TRUE,weights=weights,...)

  # anisotropic mu ?
  if(length(x)>2)
  {
    AN <- mean.mu(x,debias=debias,isotropic=FALSE,weights=weights,...)
    if(AN[[IC]]<MU[[IC]])
    {
      isotropic["mu"] <- FALSE
      MU <- AN
    }
  }

  # zero variance in mu ?
  if(isotropic["mu"])
  {
    ZV <- mean.mu(x,debias=debias,isotropic=TRUE,variance=FALSE,weights=weights,...)
    if(ZV[[IC]]<=MU[[IC]])
    {
      ZV$POV.mu <- 0 * MU$POV.mu
      ZV$COV.POV.mu <- 0 * MU$COV.POV.mu
      MU <- ZV
    }
  }

  ########################
  # average of features
  ISO <- sapply(x,function(y){y$isotropic})
  ISO <- all(ISO)
  FU <- mean.features(x,debias=debias,isotropic=TRUE,weights=weights,...)

  # anisotropic sigma ?
  if(length(x)>2 && !ISO)
  {
    AN <- mean.features(x,debias=debias,isotropic=FALSE,weights=weights,...)
    if(AN[[IC]]<FU[[IC]])
    {
      isotropic["sigma"] <- FALSE
      FU <- AN
    }
  }

  # zero variance in sigma
  if(isotropic["sigma"])
  {
    ZV <- mean.features(x,debias=debias,isotropic=TRUE,variance=FALSE,weights=weights,...)
    if(ZV[[IC]]<=FU[[IC]])
    {
      ZV$POV.sigma <- 0 * FU$POV.sigma
      ZV$COV.POV.sigma <- 0 * FU$COV.POV.sigma
      FU <- ZV
    }
  }

  info <- mean.info(x)
  x <- copy(from=FU,to=MU)
  x$info <- info
  x <- do.call(ctmm,x)

  # STORE EVERYTHING
  x$isotropic <- isotropic

  x$AIC <- MU$AIC + FU$AIC
  x$AICc <- MU$AICc + FU$AICc
  x$BIC <- MU$BIC + FU$BIC

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

