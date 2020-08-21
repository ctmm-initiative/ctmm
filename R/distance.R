####################
# square distance between stationary Gaussian distributions
BhattacharyyaD <- function(CTMM)
{
  CTMM1 <- CTMM[[1]]
  CTMM2 <- CTMM[[2]]

  sigma <- (CTMM1$sigma + CTMM2$sigma)/2
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]

  D <- as.numeric(mu %*% PDsolve(sigma) %*% mu)/8 + log(det(sigma)/sqrt(det(CTMM1$sigma)*det(CTMM2$sigma)))/2

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


distance <- function(object,method="Mahalanobis",level=0.95,debias=TRUE,...)
{
  method <- match.arg(method,c("Bhattacharyya","Mahalanobis","Euclidean"))
  n <- length(object)
  D <- array(NA_real_,c(n,n,3))

  for(i in 1:(n-1))
  {
    D[i,i,] <- c(0,0,0) # diagonal entries
    for(j in (i+1):n) # off-diagonal entries
    { D[i,j,] <- D[j,i,] <- overlap.ctmm(object[c(i,j)],level=level,debias=debias,method=method,distance=TRUE,...) }
  }
  D[n,n,] <- c(0,0,0) # diagonal entries

  dimnames(D) <- list(names(object),names(object),NAMES.CI)
  return(D)
}
