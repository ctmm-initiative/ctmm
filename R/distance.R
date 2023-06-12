####################
# Bhattacharyya distance between stationary Gaussian distributions
BhattacharyyaD <- function(CTMM)
{
  CTMM1 <- CTMM[[1]]
  CTMM2 <- CTMM[[2]]

  sigma <- (CTMM1$sigma + CTMM2$sigma)/2
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]

  D <- as.numeric(mu %*% PDsolve(sigma) %*% mu)/8 + log(det(sigma)/sqrt(det(CTMM1$sigma)*det(CTMM2$sigma)))/2

  return(D)
}

# encounter rate "distance" between stationary Gaussian distributions
RateD <- function(CTMM)
{
  CTMM1 <- CTMM[[1]]
  CTMM2 <- CTMM[[2]]

  sigma <- (CTMM1$sigma + CTMM2$sigma)/2
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]

  D <- as.numeric(mu %*% PDsolve(sigma) %*% mu)/4 + PDlogdet(sigma)/2 + nrow(sigma)/2*log(4*pi)

  return(D)
}

# encounter distance between stationary Gaussian distributions
EncounterD <- function(CTMM)
{
  CTMM1 <- CTMM[[1]]
  CTMM2 <- CTMM[[2]]

  sigma <- (CTMM1$sigma + CTMM2$sigma)/2
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]

  D <- as.numeric(mu %*% PDsolve(sigma) %*% mu)/4 + log(det(sigma)/sqrt(det(CTMM1$sigma)*det(CTMM2$sigma)))/2

  return(D)
}

# square Mahalanobis distance
MahalanobisD <- function(CTMM)
{
  CTMM1 <- CTMM[[1]]
  CTMM2 <- CTMM[[2]]

  sigma <- (CTMM1$sigma + CTMM2$sigma)/2
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]

  D <- as.numeric(mu %*% PDsolve(sigma) %*% mu)
  # D <- sqrt(D)

  return(D)
}


# square Euclidean distance
EuclideanD <- function(CTMM)
{
  CTMM1 <- CTMM[[1]]
  CTMM2 <- CTMM[[2]]

  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]

  D <- sum(mu * mu)
  # D <- sqrt(D)

  return(D)
}


distance <- function(object,method="Mahalanobis",sqrt=FALSE,level=0.95,debias=TRUE,...)
{
  method <- match.arg(method,c("Bhattacharyya","Mahalanobis","Euclidean","Encounter"))
  object <- name.list(object)
  n <- length(object)
  D <- array(NA_real_,c(n,n,3))
  DOF <- array(NA_real_,c(n,n))

  # convert everything to ctmm centroids
  for(i in 1:n)
  {
    CLASS <- class(object[[i]])[1]
    if(CLASS=="UD")
    { object[[i]] <- object[[i]]@CTMM }
    else if(CLASS=="telemetry")
    {
      if(nrow(object[[i]])>1) { stop("Multiple locations in ",names(object)[i]) }
      CTMM <- ctmm(isotropic=TRUE)
      CTMM$mu <- get.telemetry(object[[i]])
      CTMM$COV.mu <- get.error(object[[i]],ctmm(error=uere(object[[1]])$N>0),circle=FALSE,DIM=2)[1,,]
      CTMM$sigma <- covm(0,isotropic=TRUE)
      CTMM$COV <- cbind(0)
      dimnames(CTMM$COV) <- list("major","major")
      object[[i]] <- CTMM
    }
  }

  for(i in 1:(n-1))
  {
    D[i,i,] <- c(0,0,0) # diagonal entries
    DOF[i,i] <- Inf
    for(j in (i+1):n) # off-diagonal entries
    {
      STUFF <- overlap.ctmm(object[c(i,j)],level=FALSE,debias=debias,method=method,distance=TRUE,sqrt=sqrt,...)
      MLE <- STUFF$MLE
      BIAS <- STUFF$BIAS
      VAR <- STUFF$VAR

      if(debias) { MLE <- nant(MLE/BIAS,MLE) }

      dof <- nant(2*MLE^2/VAR,1/VAR)
      CI <- chisq.ci(MLE,DOF=dof,alpha=1-level)

      if(sqrt) # sqrt and debias
      {
        CI <- sqrt(CI)
        if(debias) { CI <- CI/chi.bias(max(dof,1)) }
        VAR <- chi.var(dof,CI[2])
        # effective chi^2 dof (not chi dof)
        dof <- nant(2*CI[2]^2/VAR,1/VAR)
      }

      D[i,j,] <- D[j,i,] <- CI
      DOF[i,j] <- DOF[j,i] <- dof
    }
  }
  D[n,n,] <- c(0,0,0) # diagonal entries
  DOF[n,n] <- Inf

  dimnames(D) <- list(names(object),names(object),NAMES.CI)
  dimnames(DOF) <- list(names(object),names(object))

  R <- list(DOF=DOF,CI=D)
  class(R) <- "distance"
  return(R)
}
