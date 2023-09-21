# overlap <- function(object,...) UseMethod("overlap") #S3 generic

# forwarding function for list of a particular datatype
overlap <- function(object,method="Bhattacharyya",level=0.95,debias=TRUE,...)
{
  object <- name.list(object)
  check.projections(object)

  CLASS <- class(object[[1]])[1]

  if(CLASS=="ctmm")
  { OverlapFun <- overlap.ctmm }
  # else if(CLASS=="telemetry")
  # {
  #   # Generate aligned UDs
  #   object <- akde(object,...)
  #   OverlapFun <- overlap.UD
  # }
  else if(CLASS=="UD")
  { OverlapFun <- overlap.UD }
  else { stop(CLASS," object class not supported by overlap.") }

  n <- length(object)
  DOF <- array(Inf,c(n,n))
  OVER <- array(1,c(n,n,3))
  # tabulate overlaps
  for(i in 1:n)
  {
    if(method=="Rate")
    { START <- i } # include diagonal for normalization
    else
    { START <- i+1 }

    for(j in START%:%n)
    {
      R <- OverlapFun(object[c(i,j)],level=level,debias=debias,method=method,...)
      DOF[i,j] <- DOF[j,i] <- R$DOF
      OVER[i,j,] <- OVER[j,i,] <- R$CI
    }
  }

  dimnames(OVER) <- list(names(object),names(object),NAMES.CI)
  dimnames(DOF) <- list(names(object),names(object))

  R <- list(DOF=DOF,CI=OVER)
  class(R) <- "overlap"

  return(R)
  # utils::getS3method("overlap",CLASS)(object,...)
}


# approximate-exact Wishart DOFs
DOF.wishart <- function(CTMM)
{
  AXES <- length(CTMM$axes)
  # extract covariance matrix SIGMA and covariance of covariance matrix COV
  SIGMA <- CTMM$sigma # Wishart matrix / n
  if(CTMM$isotropic) # chi^2
  { PAR <- 'major' }
  else # approximate Wishart DOF (exact if Wishart)
  { PAR <- c('major','minor') }
  EST <- SIGMA@par[PAR]
  DOF <- CTMM[['COV']][PAR,PAR]
  if(length(DOF)==length(PAR)^2)
  { DOF <- (2/AXES) * c(EST %*% PDsolve(DOF) %*% EST) } # average multiple DOFs if not Wishart
  if(length(DOF)==0) { DOF <- 0 }
  return(DOF)
}


# n/(n-dim-1) for n>=dim+2
# leading order in 1/n for n<=dim+2
# matched to first derivative
soft.clamp <- function(n,DIM)
{
  if(n >= DIM+2) { return(n) }

  A <- -2*DIM - 5*DIM^2 - 4*DIM^3 - DIM^4
  B <- 4 + 16*DIM + 25*DIM^2 + 19*DIM^3 + 7*DIM^4 + DIM^5

  # same to leading order, with matching first derivative
  BIAS <- 1 + (DIM+1)/n + A/n^2 + B/n^3

  # n that would give the above bias with n/(n-dim-1) formula
  (DIM+1)*BIAS/(BIAS-1)
}


#####################
overlap.ctmm <- function(object,level=0.95,debias=TRUE,COV=TRUE,method="Bhattacharyya",distance=FALSE,sqrt=FALSE,...)
{
  CTMM1 <- object[[1]]
  CTMM2 <- object[[2]]
  DIM <- length(CTMM1$axes)

  Dfunc <- get(paste0(method,'D')) # functions in distance.R
  STUFF <- gauss.comp(Dfunc,object,COV=COV)
  MLE <- c(STUFF$MLE)
  VAR <- c(STUFF$COV)
  # this quantity is roughly chi-square
  DOF <- nant(2*MLE^2/VAR,1/VAR)

  # approximate debiasing, correct for IID, equal covariance, REML
  ########################
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]
  COV.mu <- CTMM1$COV.mu + CTMM2$COV.mu

  if(method=="Euclidean") { sigma <- diag(1,DIM) }
  else { sigma <- (CTMM1$sigma + CTMM2$sigma)/2 }

  # trace variances
  s0 <- mean(diag(sigma))
  s1 <- mean(diag(CTMM1$sigma))
  s2 <- mean(diag(CTMM2$sigma))

  # approximate average Wishart DOFs
  # n1 <- DOF.wishart(CTMM1)
  # n2 <- DOF.wishart(CTMM2)
  # the above can be flaky
  n1 <- DOF.area(CTMM1)
  n2 <- DOF.area(CTMM2)

  # hard clamp before soft clamp
  n1 <- clamp(n1,1,Inf)
  n2 <- clamp(n2,1,Inf)

  # using mean variance - additive & rotationally invariant
  n0 <- 4 * s0^2 / (s1^2/n1 + s2^2/n2)
  n0 <- nant(n0,0)
  n0 <- clamp(n0,2,Inf)

  # clamp the DOF not to diverge <=DIM+1
  n0 <- soft.clamp(n0,DIM)
  n1 <- soft.clamp(n1,DIM)
  n2 <- soft.clamp(n2,DIM)

  # expectation value of log det Wishart
  ElogW <- function(s,n,add=TRUE) { add*PDlogdet(s) + mpsigamma(n/2,dim=DIM) - DIM*log(n/2) }

  # inverse Wishart expectation value pre-factor
  BIAS <- nant( n0/(n0-DIM-1) ,1)
  if(method=="Euclidean")
  { BIAS <- 0 } # don't include first term
  else if(method=="Rate")
  { BIAS <- BIAS - 1 } # subtractive rather than multiplicative treatment

  # mean terms
  BIAS <- sum(diag((BIAS*outer(mu) + COV.mu) %*% PDsolve(sigma)))
  if(method=="Bhattacharyya")
  {
    BIAS <- BIAS/8
    # AMGM covariance terms
    BIAS <- BIAS + max( ElogW(sigma,n0)/2 - ElogW(CTMM1$sigma,n1)/4 - ElogW(CTMM2$sigma,n2)/4 , 0)
    # this is actually the expectation value
  }
  else if(method=="Encounter") # encounter-probability overlap measure
  {
    BIAS <- BIAS/4
    # AMGM covariance terms
    BIAS <- BIAS + max( ElogW(sigma,n0)/2 - ElogW(CTMM1$sigma,n1)/4 - ElogW(CTMM2$sigma,n2)/4 , 0)
    # this is actually the expectation value
  }
  else if(method=="Rate") # encounter probability
  {
    BIAS <- BIAS/4
    # covariance terms
    BIAS <- BIAS + max( ElogW(sigma,n0,FALSE)/2 , ElogW(CTMM1$sigma,n1,FALSE)/4 + ElogW(CTMM2$sigma,n2,FALSE)/4 )
  }

  # relative bias instead of absolute bias
  if(method!="Rate") { BIAS <- nant(BIAS/MLE,MLE) }
  # would subtract off estimate to get absolute bias
  if(distance) { BIAS <- 1 + BIAS } # didn't include first term

  # error corrections
  BIAS <- as.numeric(BIAS)
  if(MLE==0) { BIAS <- 1 }
  #####################

  if(level)
  {
    if(method!="Rate")
    {
      if(debias) { MLE <- MLE/BIAS }

      DOF <- 2*MLE^2/VAR
      CI <- chisq.ci(MLE,DOF=DOF,alpha=1-level)

      if(distance) # return distance
      {
        if(sqrt) # sqrt and debias
        {
          CI <- sqrt(CI)
          if(debias) { CI <- CI/chi.bias(max(DOF,1)) }
        }

        return(CI)
      }

      # transform from (square) distance to overlap measure
      CI <- exp(-rev(CI))
    }
    else # method=="Rate" # can be negative
    {
      if(debias) { MLE <- MLE - BIAS }

      MLE <- exp(-MLE) # this is more chi^2
      VAR <- MLE^2 * VAR

      DOF <- 2*MLE^2/VAR
      CI <- chisq.ci(MLE,DOF=DOF,alpha=1-level)
    }

    names(CI) <- NAMES.CI
    R <- list(DOF=DOF,CI=CI)
    return(R)
  }
  else # return BD ingredients
  { return(list(MLE=MLE,VAR=VAR,DOF=DOF,BIAS=BIAS)) }
}


#####################
#overlap density function
overlap.UD <- function(object,level=0.95,debias=TRUE,method="Bhattacharyya",...)
{
  CTMM <- list(attr(object[[1]],"CTMM"),attr(object[[2]],"CTMM"))
  type <- c(attr(object[[1]],"type"),attr(object[[2]],"type"))
  type <- type[type!="range"]
  if(length(type)) { stop(type," overlap is not generally meaningful, biologically.") }

  # check resolution and subset to overlapping grid
  object <- same.grids(object)
  # can now be null mass

  dr <- object[[1]]$dr
  dA <- prod(dr)

  OVER <- object[[1]]$PDF * object[[2]]$PDF
  if(!is.null(OVER) && method=="Bhattacharyya") { OVER <- sqrt(OVER) }

  # overlap point estimate
  OVER <- sum(OVER)*dA

  if(!is.null(OVER) && method=="Encounter")
  {
    OVER <- OVER / sqrt(sum(object[[1]]$PDF^2)*dA*sum(object[[2]]$PDF^2)*dA)
    # no shared support after subsetting
    OVER <- nant(OVER,0)
  }

  if(!is.null(CTMM))
  {
    # calculate Gaussian overlap distance^2 variance, bias, etc.
    CI <- overlap.ctmm(CTMM,method=method,level=FALSE)

    # Bhattacharyya distances
    if(OVER>0) { D <- -log(OVER) }
    else { D <- CI$MLE }

    if(method=="Rate") # D can be negative
    {
      # additive bias
      if(debias)
      {
        D <- D - CI$BIAS
        D <- nant(D,Inf) # Inf - Inf
        OVER <- exp(-D)
      }

      VAR <- CI$VAR
      # CI <- norm.ci(D,VAR=VAR,alpha=1-level)
      # CI <- exp(-rev(CI))

      VAR <- VAR*OVER^2 # VAR[D] -> VAR[OVER]
      DOF <- 2*OVER^2/VAR

      # if(CI[3]<Inf) { CI <- chisq.ci(OVER,DOF=DOF,alpha=1-level) }
      CI <- chisq.ci(OVER,DOF=DOF,alpha=1-level)
    }
    else # normalized overlap # can't be negative
    {
      # relative debias
      if(debias){ D <- D/CI$BIAS }

      # calculate new distance^2 with KDE point estimate
      DOF <- 2*D^2/CI$VAR
      CI <- chisq.ci(D,DOF=DOF,alpha=1-level)

      # transform from (square) distance to overlap measure
      CI <- exp(-rev(CI))
    }
  }
  else
  {
    DOF <- NA
    CI <- c(NA,OVER,NA)
  }

  names(CI) <- NAMES.CI

  R <- list(DOF=DOF,CI=CI)
  return(R)
}
