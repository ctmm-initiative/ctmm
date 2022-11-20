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

  D <- as.numeric(mu %*% PDsolve(sigma) %*% mu)/4 + log(det(sigma))/2
  # modulo log(2*pi*2^dim)

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
  n <- length(object)
  D <- array(NA_real_,c(n,n,3))

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
    for(j in (i+1):n) # off-diagonal entries
    { D[i,j,] <- D[j,i,] <- overlap.ctmm(object[c(i,j)],level=level,debias=debias,method=method,distance=TRUE,sqrt=sqrt,...) }
  }
  D[n,n,] <- c(0,0,0) # diagonal entries

  dimnames(D) <- list(names(object),names(object),NAMES.CI)
  return(D)
}
