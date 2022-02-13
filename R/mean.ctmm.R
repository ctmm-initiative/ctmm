###########
# log-transformed parameters
# debias includes bias correction for chi^2 to log(chi^2)
# matrix casts location covariance, diffusion rate, velocity covariance all as distinct matrices for above bias correction
log.ctmm <- function(CTMM,debias=FALSE)
{
  isotropic <- CTMM$isotropic
  axes <- CTMM$axes
  features <- CTMM$features
  par <- get.parameters(CTMM,features,linear.cov=TRUE)
  sigma <- attr(CTMM$sigma,"par")
  COV <- CTMM$COV

  ### log transform sigma
  # log transform major-minor
  if(isotropic)
  {
    PARS <- "major"
    sigma <- sigma[PARS]
  }
  else
  { PARS <- c("major","minor") }
  par[PARS] <- sigma[PARS] <- log(sigma[PARS])
  COV[PARS,] <- COV[PARS,,drop=FALSE] / sigma[PARS]
  COV[,PARS] <- COV[,PARS,drop=FALSE] / sigma[PARS]

  # convert log(eigen) to log(xy)
  sigma <- covm(sigma,isotropic=isotropic,axes=axes)
  if(!isotropic)
  {
    par[PARS] <- sigma[c(1,4,2)]  # log 'xx', 'yy', 'xy'

    PARS <- c("major","minor","angle")
    J <- J.sigma.par(sigma)
    COV[PARS,] <- J %*% COV[PARS,,drop=FALSE]
    COV[,PARS] <- COV[,PARS,drop=FALSE] %*% t(J)
  }

  ### log transform all positive parameters that haven't been logged
  # features to log transform
  FEAT <- features[features %in% POSITIVE.PARAMETERS]
  FEAT <- FEAT[FEAT %nin% names(CTMM$sigma@par)]

  if(length(FEAT))
  {
    COV[FEAT,] <- COV[FEAT,] / par[FEAT]
    COV[,FEAT] <- COV[,FEAT] / par[FEAT]
    par[FEAT] <- log(par[FEAT])
  }

  # log chi^2 bias correction
  if(debias)
  {
    # diagonalize and log-chi^2 debias relevant parameter estimates
    EIGEN <- eigen(COV)
    dimnames(EIGEN$vectors) <- list(features,features)
    names(EIGEN$values) <- features

    # fix signs
    if(isotropic) { PARS <- "major" } else { PARS <- c("major","minor") }
    # VAR goes in log numerator for chi^2 variates: variance, diffusion, MS speed, ...
    for(i in 1:nrow(EIGEN$vectors)) { if(sum(EIGEN$vectors[i,PARS])<0) { EIGEN$vectors[i,] <- -EIGEN$vectors[i,] } }

    # transform to diagonalized basis with VARs in log numerator
    par <- t(EIGEN$vectors) %*% par # diagonalize parameters
    DOF <- 2/EIGEN$values # log-chi^2 VAR-DOF relation
    BIAS <- digamma(DOF/2)-log(DOF/2) # negative bias for log(chi^2) variates
    par <- par + BIAS # E[log-chi^2] bias correction
    par <- EIGEN$vectors %*% par # transform back (still under logarithm)
  }

  RETURN <- list(par=par,COV=COV)
  return(RETURN)
}


#####################
# inverse transformation of above
exp.ctmm <- function(CTMM,debias=FALSE)
{
  isotropic <- CTMM$isotropic
  axes <- CTMM$axes
  features <- CTMM$features
  par <- get.parameters(CTMM,features,linear.cov=TRUE)
  COV <- CTMM$COV
  POV <- CTMM$log.PCOV

  # log chi^2 bias correction
  if(debias)
  {
    # diagonalize and log-chi^2 debias relevant parameter estimates
    EIGEN <- eigen(COV+POV)
    dimnames(EIGEN$vectors) <- list(features,features)
    names(EIGEN$values) <- features

    # fix signs
    if(isotropic) { PARS <- "major" } else { PARS <- c("major","minor") }
    # VAR goes in log numerator for chi^2 variates: variance, diffusion, MS speed, ...
    for(i in 1:nrow(EIGEN$vectors)) { if(sum(EIGEN$vectors[i,PARS])<0) { EIGEN$vectors[i,] <- -EIGEN$vectors[i,] } }

    # transform to diagonalized basis with VARs in log numerator
    par <- t(EIGEN$vectors) %*% par # diagonalize parameters
    DOF <- 2/EIGEN$values # log-chi^2 VAR-DOF relation
    BIAS <- digamma(DOF/2)-log(DOF/2) # negative bias for log(chi^2) variates
    par <- par - BIAS # E[log-chi^2] bias correction
    par <- EIGEN$vectors %*% par # transform back (still under logarithm)
  }

  ### exp transform all positive parameters except sigma
  # features to log transform
  FEAT <- features[features %in% POSITIVE.PARAMETERS]
  FEAT <- FEAT[FEAT %nin% names(CTMM$sigma@par)]

  if(length(FEAT))
  {
    par[FEAT] <- exp(par[FEAT])
    COV[FEAT,] <- COV[FEAT,] * par[FEAT]
    COV[,FEAT] <- COV[,FEAT] * par[FEAT]
  }

  ### exp transform sigma
  sigma <- attr(CTMM$sigma,"par")

  # convert COV[log(xy)] to COV[log(eigen)]
  if(!isotropic)
  {
    PARS <- names(sigma)
    J <- J.par.sigma(sigma)
    COV[PARS,] <- J %*% COV[PARS,,drop=FALSE]
    COV[,PARS] <- COV[,PARS,drop=FALSE] %*% t(J)

    PARS <- c("major","minor")
  }
  else
  { PARS <- "major" }

  # now log (major,minor,angle)
  sigma[PARS] <- exp(sigma[PARS])
  COV[PARS,] <- COV[PARS,] * sigma[PARS]
  COV[,PARS] <- COV[,PARS] * sigma[PARS]

  sigma <- covm(sigma,axes=axes,isotropic=isotropic)
  CTMM$sigma <- sigma
  CTMM$COV <- COV

  RETURN <- list(par=par,COV=COV)
  return(RETURN)
}


#############
mean.features <- function(x,debias=TRUE,isotropic=FALSE,...)
{
  N <- length(x)

  # all possible features
  FEATURES <- lapply(x,function(y){y$features})
  FEATURES <- unlist(FEATURES)
  FEATURES <- unique(FEATURES)
  M <- length(FEATURES)

  MEANS <- VARS <- rep(TRUE,M)
  names(VARS) <- names(VARS) <- FEATURES

  if(isotropic && "minor" %in% FEATURES)
  { MEANS[c("minor","angle")] <- VARS[c("minor","angle")] <- c(FALSE,FALSE) }

  MU <- array(0,c(N,M))
  SIGMA <- array(0,c(N,M,M))
  INF <- diag(Inf,M)
  colnames(MU) <- FEATURES
  dimnames(SIGMA) <- list(NULL,FEATURES,FEATURES)

  for(i in 1:N)
  {
    x[[i]] <- log.ctmm(x[[i]])
    features <- x[[i]]$features
    MU[i,features] <- x[[i]]$par
    SIGMA[i,,] <- INF
    SIGMA[i,features,features] <- x[[i]]$COV

    P <- c("major","minor")
    if(!isotropic && x[[i]]$isotropic)
    {
      MU[i,'minor'] <- MU[i,'major']
      SIGMA[i,P] <- matrix(SIGMA[i,'major','major'],2,2)
    }
    else if(isotropic && !x[[i]]$isotropic)
    {
      MU[i,P] <- c(mean(MU[i,P]),diff(MU[i,P]))
      J <- rbind( c(1/2,1/2) , c(1,-1) )
      SIGMA[i,P,] <- J %*% SIGMA[i,P,]
      SIGMA[i,,P] <- SIGMA[i,,P] %*% t(J)
    }
  }

  R <- meta.normal(MU,SIGMA,MEANS=MEANS,VARS=VARS,debias=debias)
  names(R)[ which(names(R)=="mu") ] <- "par" # population mean of features
  names(R)[ which(names(R)=="COV.mu") ] <- "COV" # uncertainty in population mean
  names(R)[ which(names(R)=="sigma") ] <- "log.PCOV" # population dispersion of features (under log)
  names(R)[ which(names(R)=="COV.sigma") ] <- "COV.log.PCOV" # uncertainty in population dispersion of features (under log)

  # information for fitting that we no longer use
  if(isotropic && "minor" %in% FEATURES)
  {
    FEATURES <- FEATURES[FEATURES %nin% c('minor','angle')]
    par <- par[FEATURES]
    R$par <- NULL
    R$COV <- R$COV[FEATURES,FEATURES]
    R$log.PCOV <- R$log.PCOV[FEATURES,FEATURES]
    # COV.log.PCOV should be okay
  }

  # reverse log transformation
  R <- set.parameters(R,par,linear.cov=TRUE)
  R <- exp.ctmm(R)

  return(R)
}


##########
mean.mu <- function(x,debias=TRUE,isotropic=FALSE,...)
{
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
  R <- meta.normal(MU,SIGMA,isotropic=isotropic,debias=debias)
  # R$mu # population mean of mean locations
  # R$COV.mu # uncertainty in mean of means estimate
  names(R)[ which(names(R)=="sigma") ] <- "PCOV.mu" # dispersion of means
  names(R)[ which(names(R)=="COV.sigma") ] <- "COV.PCOV.mu" # uncertainty in dispersion of means

  return(R)
}


###########
mean.ctmm <- function(x,debias=TRUE,IC="AICc",...)
{
  isotropic <- c(TRUE,TRUE)
  names(isotropic) <- c("sigma","mu")

  # average of mean locations
  MU <- mean.mu(x,debias=debias,isotropic=TRUE,...)
  if(length(x)>2)
  {
    AN <- mean.mu(x,debias=debias,isotropic=FALSE,...)
    if(AN[[IC]]<MU[[IC]])
    {
      isotropic["mu"] <- FALSE
      MU <- AN
    }
  }

  # average of features
  ISO <- sapply(x,function(y){y$isotropic})
  ISO <- all(ISO)
  FU <- mean.features(x,debias=debias,isotropic=TRUE,...)
  if(!ISO)
  {
    AN <- mean.features(x,debias=debias,isotropic=FALSE,...)
    if(AN[[IC]]<FU[[IC]])
    {
      isotropic["sigma"] <- FALSE
      FU <- AN
    }
  }

  info <- mean.info(x)
  x <- FU

  # STORE EVERYTHING
  x$isotropic <- isotropic

  x$AIC <- MU$AIC + FU$AIC
  x$AICc <- MU$AICc + FU$AICc
  x$BIC <- MU$BIC + FU$BIC

  x$mu <- MU$mu
  x$COV.mu <- MU$COV.mu

  # return final result
  x <- new.ctmm(x,info)
  return(x)
}

# TODO !!!!!!!!!!!!!!!!!!!!!
# function to calculate the population range from the above
