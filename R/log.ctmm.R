###########
# log-transformed parameters
# debias includes bias correction for chi^2 to log(chi^2)
# matrix casts location covariance, diffusion rate, velocity covariance all as distinct matrices for above bias correction
log_ctmm <- function(CTMM,debias=FALSE,...)
{
  if(class(CTMM)[1]=="list") { return(log_ctmms(CTMM,debias=debias,...)) }

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

  sigma <- log_covm(sigma)

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
    if(any(GOOD)) { COR[GOOD,GOOD] <- stats::cov2cor(COV[GOOD,GOOD,drop=FALSE]) }
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

  RETURN <- list(par=par,COV=COV,isotropic=isotropic)
  return(RETURN)
}


#####################
# inverse transformation of above
exp_ctmm <- function(CTMM,debias=FALSE,variance=TRUE)
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
    BIAS <- log_chi2_bias(DOF) # negative bias for log(chi^2) variates
    BIAS <- pmax(BIAS,log_chi2_bias(1)) # clamp to 1 DOF
    par[SUB] <- par[SUB] + BIAS # E[log-chi^2] bias correction
    par[SUB] <- c(EIGEN$vectors %*% par[SUB]) # transform back (still under logarithm)

    ################
    ### diagonalize and log-chi^2 debias COV
    VAR <- diag(COV)
    COR <- diag(1,nrow(COV))
    dimnames(COR) <- dimnames(COV)
    GOOD <- VAR>.Machine$double.eps
    if(any(GOOD)) { COR[GOOD,GOOD] <- stats::cov2cor(COV[GOOD,GOOD,drop=FALSE]) }

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
      JP[SUB,SUB] <- EIGEN$vectors %*% diag(SCALE,nrow=length(SCALE)) %*% t(EIGEN$vectors)
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
  sigma <- exp_covm(sigma)

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
