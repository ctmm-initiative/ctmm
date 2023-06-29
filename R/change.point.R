# functions for estimating a change point that happens at an instant in time

# this function returns the most likely IID change point in a dataset
change.point.iid <- function(data,axes=c('x','y'),IC="AICc")
{
  AXES <- length(axes)
  n <- nrow(data)
  if(n<4) { return(0) }

  z <- get.telemetry(data,axes)

  # forward running mean, covariance
  M1 <- z
  S1 <- array(0,c(n,AXES,AXES))

  # Welford's algorithm
  for(i in 2:n)
  {
    M1[i,] <- ( (i-1)*M1[i-1,] + z[i,] )/i
    S1[i,,] <- S1[i-1,,] + (i-1)/i*outer(z[i,]-M1[i-1,])
  }
  N1 <- 1:n
  S1 <- S1 / (N1 - 1)

  # backward running mean, covariance
  M2 <- z
  S2 <- array(0,c(n,AXES,AXES))

  # Welford's algorithm
  for(i in 2:n)
  {
    j <- n-i+1
    M2[j,] <- ( (i-1)*M2[j+1,] + z[j,] )/i
    S2[j,,] <- S2[j+1,,] + (i-1)/i*outer(z[j,]-M2[j+1,])
  }
  N2 <- n:1
  S2 <- S2 / (N2 - 1)

  # fix 0/0
  S1[1,,] <- S2[n,,] <- diag(Inf,AXES)

  # uncorrelated log-likelihoods -2x and without log(2pi)
  NLLU1 <- AXES*(N1*( rowSums( sapply(1:AXES,function(d){log(S1[,d,d])}) ) + 1 ) - 1 )
  NLLU2 <- AXES*(N2*( rowSums( sapply(1:AXES,function(d){log(S2[,d,d])}) ) + 1 ) - 1 )
  # correlated log-likelihoods
  NLLC1 <- AXES*(N1*( sapply(1:n,function(i){PDlogdet(S1[i,,])}) + 1 ) - 1 )
  NLLC2 <- AXES*(N2*( sapply(1:n,function(i){PDlogdet(S2[i,,])}) + 1 ) - 1 )

  # not enough data for parameters (uncorrelated case)
  NLLU1[1] <- NLLU2[n] <- Inf
  # not enough data for parameters (K parameters - correlated case)
  K <- AXES + (AXES^2+AXES)/2
  MAX <- floor(K/AXES)
  NLLC1[1:MAX] <- NLLC2[n:(n-MAX+1)] <- Inf

  if(IC=="AIC")
  {
    # uncorrelated AICc
    ICU1 <- NLLU1 + 2*( 2*AXES )
    ICU2 <- NLLU2 + 2*( 2*AXES )
    # correlated AICc
    ICC1 <- NLLC1 + 2*( K )
    ICC2 <- NLLC2 + 2*( K )
  }
  else if(IC=="AICc")
  {
    # uncorrelated AICc
    ICU1 <- NLLU1 + AXES*(N1-1)*(2+1+1)/max(N1-1-1-1,0)
    ICU2 <- NLLU2 + AXES*(N2-1)*(2+1+1)/max(N2-1-1-1,0)
    # correlated AICc
    ICC1 <- NLLC1 + AXES*(N1-1)*(2+AXES+1)/max(N1-1-AXES-1,0)
    ICC2 <- NLLC2 + AXES*(N2-1)*(2+AXES+1)/max(N2-1-AXES-1,0)
  }
  else if(IC=="BIC")
  {
    # uncorrelated AICc
    ICU1 <- NLLU1 + log(n)*( 2*AXES )
    ICU2 <- NLLU2 + log(n)*( 2*AXES )
    # correlated AICc
    ICC1 <- NLLC1 + log(n)*( K )
    ICC2 <- NLLC2 + log(n)*( K )
  }

  # selected model's IC
  IC1 <- pmin(ICU1,ICC1)
  IC2 <- pmin(ICU2,ICC2)

  # total IC
  IC <- c(0,IC1) + c(IC2,0)
  I <- which.min(IC) - 1
  IC <- IC[I] - IC[1] # difference from no change point

  R <- list(MIN=I,IC=IC)
  return(R)
}

change.point.guess <- function(data,n=1,axes=c('x','y'),...)
{
  CP <- NULL
  for(i in 1:n)
  {
    # segment delimiters
    AP <- c(0,CP,nrow(data))
    AP <- unique(AP)
    m <- length(AP)-1
    for(j in 1:m) # try every segment
    {
      MIN <- IC <- numeric(m)
      SUB <- (AP[j]+1):AP[j+1]
      STUFF <- change.point.iid(data[SUB,],axes=axes)
      MIN[j] <- STUFF$MIN + AP[j]
      IC[j] <- STUFF$IC
    }
    # !NULL solutions
    SUB <- MIN>0
    MIN <- MIN[SUB]
    IC <- IC[SUB]
    if(length(MIN))
    {
      MIN <- MIN[which.min(IC)]
      CP <- c(CP,MIN)
      CP <- sort(CP)
    }
    else # no solutions
    {
      n <- m-1
      break
    }
  }

  # change points are determined, now form guess objects
  GUESS <- list()
  AP <- c(0,CP,nrow(data))
  m <- length(AP)-1
  for(j in 1:m)
  {
    SUB <- (AP[j]+1):AP[j+1]
    GUESS[[j]] <- ctmm.guess(data[SUB,],ctmm(axes=axes),interactive=FALSE)
  }
  t <- data$t[CP]
  names(GUESS) <- names(t) <- paste0("state.",1:length(t))
  GUESS$commute <- ctmm.commute(GUESS)

  GUESS$axes <- axes
  GUESS$dynamics <- "change.point"
  GUESS$change.point <- t
  GUESS <- new.ctmm(GUESS,info=data@info)
  return(GUESS)
}

ctmm.commute <- function(M)
{
  COM <- sapply(M,M$sigma@par['angle'])
  COM <- abs(diff(COM)) > .Machine$double.eps
  !any(COM)
}

change.point.set <- function(CTMM)
{
  states <- levels(CTMM$change.points$state)
  axes <- CTMM$axes

  get.max <- function(x,ext=identity,norm=max)
  {
    x <- lapply(states,function(s){CTMM[[s]][[x]]})
    K <- sapply(x,length)
    x <- x[K==max(K)]
    x <- sapply(x,ext) # [tau,ind]
    if(length(dim(x))==3)
    { x <- apply(x,1:2,norm) } # not sure about this margins here !!!!!!!!!
    else if(length(dim(x))==2)
    { x <- apply(x,1,norm) }
    else
    { x <- max(x) }
    return(x)
  }

  isotropic <- as.logical( get.max('isotropic',norm=all) )

  CTMM$tau <- get.max('tau')
  CTMM$omega <- get.max('omega')
  CTMM$circle <- get.max('circle')
  CTMM$sigma <- covm( diag(diag(get.max('sigma')),length(axes)), axes=axes,isotropic=isotropic)

  return(CTMM)
}
