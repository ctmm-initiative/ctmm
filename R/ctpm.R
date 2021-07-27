ctpm.solve <- function(lag,CTMM,log.det=0,PC=FALSE,INV=!PC)
{
  TEST <- CTMM
  TEST$sigma <- covm(1,axes=CTMM$axes,CTMM$isotropic) # unit variance for SVF/ACF below
  COR <- svf.func(TEST)
  if(range) # stationary likelihood
  { COR <- COR$ACF }
  else # likelihood that avoids explicit conditioning
  { COR <- COR$SVF }
  COR <- Vectorize(COR)
  COR <- COR(lag)
  dim(COR) <- dim(lag)
  if(!range) { COR <- -COR } # ACF = COV - SVF

  N <- nrow(lag)

  # 1st: construct an approximate covariance matrix that we can solve exactly
  lambda <- COR
  diag(lambda) <- 0
  lambda <- sum(lambda)/(N^2-N) # average off-diagonal element
  # COR ~ Lambda = (1-lambda) I + lambda 1o1
  # where 1 is the vector of 1s and 1o1 is outer(1,1)
  # this matrix has two eigen-values
  # Lambda * 1 = (1+(N-1)*lambda) * 1
  # Lambda * V = (1-lambda) * V
  # where V is any vector orthogonal 1, which satisfies sum(V)=0
  lambda <- c( 1 + (N-1)*lambda , 1 - lambda )

  log.det <- rep(log(lambda[2]),N)
  log.det[1] <- log(lambda[1])

  ONE <- rep(1/sqrt(N),N) # eigen-vector 1 normalized
  # inverse square-root of Lambda
  ONE <- outer(ONE)

  lambda <- 1/sqrt(lambda)
  pc <- lambda[1]*ONE + lambda[2]*(diag(N)-ONE)

  COR <- pc %*% COR %*% pc

  EIGEN <- eigen(COR,symmetric=TRUE)
  log.det <- log.det + EIGEN$values

  if(INV) # inverse: COR^(-1)
  {
    INV <- Reduce('+',lapply(1:N,function(i){EIGEN$values[i] * outer(EIGEN$vectors[,i])}))
    INV <- pc %*% INV %*% pc
  }
  else if(PC) # preconditioner: COR^(-1/2)
  {
    PC <- Reduce('+',lapply(1:N,function(i){1/sqrt(EIGEN$values[i]) * outer(EIGEN$vectors[,i])}))

    lambda <- sqrt(lambda) # now 1/quarter-root
    pc <- lambda[1]*ONE + lambda[2]*(diag(N)-ONE)

    PC <- pc %*% PC %*% pc
  }

  RETURN <- list(log.det=log.det,INV=INV,PC=PC)
}


# TODO link function adjustment
# TODO zero before sum of log.det, quadratic
# TODO multivariate case
ctpm.loglike <- function(data,CTMM,link=identity,REML=FALSE,profile=TRUE,zero=0,verbose=FALSE)
{
  lag <- attr(data,"lag")
  trait <- c( get.telemetry(data,axes=CTMM$axes) )
  n <- length(trait)
  K <- length(CTMM$tau)

  range <- CTMM$range
  if(!range)
  { REML <- TRUE }
  else # catch bad parameters
  {
    if(length(CTMM$tau) && CTMM$tau[1]==Inf)
    { return(-Inf) }
  }

  VAR.MULT <- n/(n-REML)

  # default identity link function for traits
  if(class(link)=="function")
  {
    if(identical(link,identity))
    { link <- list(fn=identity,grad=function(x){rep(1,length(x))}) }
    else if(identical(link,log))
    { link <- list(fn=log,grad=function(x){1/x}) }
  }
  # full transform is used so that link function can be fitted / selected
  grad <- link$grad(trait)
  trait <- link$fn(trait)

  if(!length(CTMM$tau)) # IID analytic solution
  {
    mu <- mean(trait)
    trait <- trait - mu
    Q <- trait^2
    COV.mu <- 1/n
  }
  else # requires preconditioned solver
  {
    if(length(CTMM$tau)>1) # OU preconditioner
    {
      TEMP <- CTMM
      TEMP$tau <- TEMP$tau[1]
      STUFF <- ctpm.solve(lag,TEMP,PC=TRUE)

      PC <- STUFF$PC
      log.det <- STUFF$log.det


    }

    STUFF <- ctpm.solve(lag,CTMM,INV=TRUE)
    log.det <- STUFF$log.det
    INV <- STUFF$INV



    TEST <- CTMM
    TEST$sigma <- covm(1,axes=CTMM$axes,CTMM$isotropic) # unit variance for SVF/ACF below
    COR <- svf.func(TEST)
    if(range) # stationary likelihood
    { COR <- COR$ACF }
    else # likelihood that avoids explicit conditioning
    { COR <- COR$SVF }
    COR <- Vectorize(COR)
    COR <- COR(lag)
    dim(COR) <- dim(lag)
    if(!range) { COR <- -COR } # ACF = COV - SVF

    # THIS DIDN'T WORK
    # # OU preconditioner for OUF
    # if(length(CTMM$tau)>1)
    # {
    #
    #   TEST$tau <- TEST$tau[1]
    #   PC <- svf.func(TEST)
    #   if(range) # stationary likelihood
    #   { PC <- PC$ACF }
    #   else # likelihood that avoids explicit conditioning
    #   { PC <- PC$SVF }
    #   PC <- Vectorize(PC)
    #   PC <- PC(lag)
    #   dim(PC) <- dim(lag)
    #   if(!range) { PC <- -PC } # ACF = COV - SVF
    #
    #   EIGEN <- eigen(PC,symmetric=TRUE)
    #   # SVF case needs to drop the large negative eigenvalue, which pairs to infinite variance term dropped
    #   if(!range)
    #   {
    #     EIGEN$values <- EIGEN$values[-n]
    #     EIGEN$vectors <- EIGEN$vectors[,-n]
    #   }
    #   N <- length(EIGEN$values)
    #
    #   if(tail(EIGEN$values,1)<=0) { return(-Inf) }
    #
    #   log.det <- log(EIGEN$values)
    #
    #   EIGEN$values <- 1/sqrt(EIGEN$values)
    #   PC <- Reduce('+',lapply(1:N,function(i){EIGEN$values[i] * outer(EIGEN$vectors[,i])}))
    #
    #   COR <- PC %*% COR %*% PC # shrinks matrix condition number
    # }
    # else
    # {
    #   log.det <- 0
    # }

    # how to do this for BM/IOU
    #
    #

    N <- nrow(COR)

    # 1st: construct an approximate covariance matrix that we can solve exactly
    lambda <- COR
    diag(lambda) <- 0
    lambda <- sum(lambda)/(N^2-N) # average off-diagonal element
    # COR ~ Lambda = (1-lambda) I + lambda 1o1
    # where 1 is the vector of 1s and 1o1 is outer(1,1)
    # this matrix has two eigen-values
    # Lambda * 1 = (1+(N-1)*lambda) * 1
    # Lambda * V = (1-lambda) * V
    # where V is any vector orthogonal 1, which satisfies sum(V)=0
    lambda <- c( 1 + (N-1)*lambda , 1 - lambda )

    log.det <- rep(log(lambda[2]),N)
    log.det[1] <- log(lambda[1])

    ONE <- rep(1/sqrt(N),N) # eigen-vector 1 normalized
    # inverse square-root of Lambda
    PC <- outer(ONE)
    PC <- PC/sqrt(lambda[1]) + (diag(N)-PC)/sqrt(lambda[2])

    COR <- PC %*% COR %*% PC

    EIGEN <- eigen(COR,symmetric=TRUE)

    # SVF case needs to drop the large negative eigenvalue, which pairs to infinite variance term dropped
    # if(!range)
    # {
    #   EIGEN$values <- EIGEN$values[-n]
    #   EIGEN$vectors <- EIGEN$vectors[,-n]
    # }

    if(tail(EIGEN$values,1)<=0) { return(-Inf) }

    if(range) # profile mean
    {
      iCOR <- Reduce('+',lapply(1:N,function(i){EIGEN$values[i] * outer(EIGEN$vectors[,i])}))
      iCOR <- PC %*% iCOR %*% PC
      W <- colSums(iCOR)
      COV.mu <- 1/sum(W)
      W <- W * COV.mu
      mu <- c(W %*% trait)
    }
    else # ignored projection (from infinite variance limit)
    {
      mu <- mean(trait)
      COV.mu <- Inf
    }

    # detrend mean
    trait <- c(PC %*% (trait - mu))

    # slow calculation
    # Q <- c(trait %*% iCOR %*% trait)
    # obtuser but faster
    Q <- (trait %*% EIGEN$vectors)^2/EIGEN$values
    # not yet summed
    log.det <- log.det + log(EIGEN$values)
    # not yet summed
    log.det <- mean(log.det) # per N
  }

  # profile the variance / diffusion rate
  if(profile)
  {
    sigma <- sum(Q) / (n-REML)
    Q <- 0 # - MLE per N
  }
  else
  {
    sigma <- attr(CTMM$sigma,"par")[1]
    Q <- mean(Q)/sigma - 1/VAR.MULT # - MLE per N
  }
  COV.mu <- sigma * COV.mu

  # unchanging constants
  LL.CONST <- - 1/2/VAR.MULT - 1/2*log(2*pi)

  loglike <- -Q/2 - log(sigma)/2 - log.det/2 # per N
  loglike <- N*(loglike + (LL.CONST-zero/N))

  if(range && REML) { loglike <- loglike + log(COV.mu)/2 }

  if(verbose)
  {
    # assign variables
    if(profile)
    {
      sigma <- covm(sigma,isotropic=TRUE,axes=CTMM$axes)
      CTMM$sigma <- sigma
    }

    CTMM <- ctmm.repair(CTMM,K=K)

    # mu <- drift@shift(mu,mu.center) # translate back to origin from center
    CTMM$mu <- mu
    CTMM$COV.mu <- COV.mu

    CTMM$loglike <- loglike
    #attr(CTMM,"info") <- attr(data,"info")

    CTMM <- ctmm.ctmm(CTMM)
    return(CTMM)
  }
  else
  {
    return(loglike)
  }
}
