corridor <- function(data,CTMM,res.space=10,res.time=100,window=1%#%'day',grid=list(),...)
{
  debug <- FALSE
  axes <- CTMM[[1]]$axes
  info <- mean_info(data)
  # window <- 1 %#% 'day' # minimum window for greedy search
  n <- length(data)
  W <- list()
  dt <- array(0,n)
  for(i in 1:n)
  {
    t <- seq(data[[i]]$t[1],last(data[[i]]$t),length.out=res.time)
    dt[i] <- stats::median(diff(t))
    # use instantaneous speeds as weights for effective distance sampling
    W[[i]] <- speeds(data[[i]],CTMM[[i]],t=t)[,'est']
    # total weight == total distance
    W[[i]] <- W[[i]]*dt[i]
    # W[[i]] <- rep(1,res.time)
    # regularize data
    data[[i]] <- predict(data[[i]],CTMM[[i]],t=t)
  }

  # convert from time to index
  window <- ceiling(window/dt)

  END <- 1:res.time
  DNE <- res.time:1

  ## orient all tracks in the same direction for NN search
  # forward variance minus reverse variance
  VFB <- array(0,c(n,n))
  for(i in 1:n)
  {
    for(j in (i+1)%:%n)
    { VFB[i,j] <- VFB[j,i] <- sum((data[[i]][END,axes]-data[[j]][END,axes])^2) - sum((data[[i]][END,axes]-data[[j]][DNE,axes])^2) }
  }
  # positive VFB is bad

  I <- apply(abs(VFB),1,sum) # importance
  I <- sort(I,index.return=TRUE)$ix # least to most importance
  for(i in I)
  {
    # check if we need to reverse order
    MAX <- +max(VFB[i,]) # worst
    MIN <- -min(VFB[i,]) # best
    if(MAX>MIN) # worst is worse than best
    {
      data[[i]] <- data[[i]][res.time:1,]
      VFB[i,] <- VFB[,i] <- -VFB[i,]
    }
  }

  NN <- array(1,c(res.time,n,n)) # nearest neighbor
  DIST <- array(0,c(res.time,n,n)) # NN pairwise distances^2

  # greedy search algorithm
  search <- function(t,i,j,WIN=NULL)
  {
    if(i==j)
    {
      NN[t,i,j] <<- t
      DIST[t,i,j] <<- 0
      return()
    }

    if(is.null(WIN))
    {
      # initial search point
      s <- NN[t,i,j]
      # one-day window
      WIN <- max(s-window[j],1):min(s+window[j],res.time)
    }
    # all distances^2 in window
    D <- (data[[i]]$x[t]-data[[j]]$x[WIN])^2 + (data[[i]]$y[t]-data[[j]]$y[WIN])^2

    # closest point
    s <- which.min(D)
    FIRST <- s==1
    LAST <- s==length(WIN)

    D <- D[s]
    s <- WIN[s]

    NN[t,i,j] <<- s
    DIST[t,i,j] <<- D

    ## greedy search backwards
    while(FIRST && s>1)
    {
      s <- s-1
      D <- (data[[i]]$x[t]-data[[j]]$x[s])^2 + (data[[i]]$y[t]-data[[j]]$y[s])^2
      if(D<DIST[t,i,j])
      {
        NN[t,i,j] <<- s
        DIST[t,i,j] <<- D
      }
      else
      { break }
    }

    ## greedy search forwards
    while(LAST && s<res.time)
    {
      s <- s+1
      D <- (data[[i]]$x[t]-data[[j]]$x[s])^2 + (data[[i]]$y[t]-data[[j]]$y[s])^2
      if(D<DIST[t,i,j])
      {
        NN[t,i,j] <<- s
        DIST[t,i,j] <<- D
      }
      else
      { break }
    }

    return()
  } # search()

  ## greedy search forward
  for(t in 1:res.time)
  {
    for(i in 1:n) { for(j in 1:n) { search(t,i,j) } }

    # copy to next
    if(t<res.time)
    {
      NN[t+1,,] <- NN[t,,]
      DIST[t+1,,] <- DIST[t,,]
    }
  } # for(t in 1:res)

  # store results
  NN.0 <- NN
  DIST.0 <- DIST

  ## greedy search backwards
  # reset starting point
  NN <- array(res.time,c(res.time,n,n)) # nearest neighbor
  for(t in res.time:1)
  {
    for(i in 1:n) { for(j in 1:n) { search(t,i,j) } }

    # copy to next
    if(t>1)
    {
      NN[t-1,,] <- NN[t,,]
      DIST[t-1,,] <- DIST[t,,]
    }
  } # for(t in res:1)

  # search between any discrepancies
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      SUB <- NN[,i,j]!=NN.0[,i,j]
      SUB <- which(SUB)

      for(t in SUB)
      {
        WIN <- range(NN[t,i,j],NN.0[t,i,j])
        WIN <- WIN[1]:WIN[2]
        search(t,i,j,WIN=WIN)
      }
    }
  }
  rm(NN.0,DIST.0)

  # calculate NN variances
  VAR <- array(0,c(res.time,n))
  VAR <- apply(DIST,1:2,sum) / (4*(n-1))
  rm(DIST)

  # predict missing locations based on NNs
  X <- M <- list()
  E <- list()
  for(i in 1:n)
  {
    X[[i]] <- M[[i]] <- get.telemetry(data[[i]])
    E[[i]] <- get.error(data[[i]],CTMM=ctmm(error=TRUE),DIM=2)
  }

  MU <- array(0,c(n,2))
  SIGMA <- array(0,c(n,2,2))
  for(t in 1:res.time)
  {
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        MU[j,] <- X[[j]][NN[t,i,j],]
        SIGMA[j,,] <- E[[j]][NN[t,i,j],,]
      }

      R <- meta.normal(MU,SIGMA,isotropic=TRUE)
      # population parameters
      mu <- R$mu
      sigma <- R$sigma

      # update variance estimate
      # VAR[t,i] <- sigma[1,1]
      # this can collapse to zero

      # update location prediction
      if(E[[i]][t,1,1]+E[[i]][t,2,2] > 2*.Machine$double.eps)
      { M[[i]][t,] <- mu + sigma %*% PDsolve(sigma+E[[i]][t,,]) %*% (M[[i]][t,]-mu) }

      # conservative variance
      # VAR[t,i] <- sum(sapply(1:n,function(j){(X[[j]][t,]-M[[i]][t,])^2})) / (2*(n-1))
      # too big
    }

    # adjusted variance
    # for(i in 1:n) { VAR[t,i] <- sum(sapply(1:n,function(j){(M[[j]][t,]-M[[i]][t,])^2})) / (4*(n-1)) }
    # too big
  }

  # copy adjusted means for kernel placement
  for(i in 1:n) { data[[i]][,axes] <- M[[i]] }

  H <- list()
  h <- (4/3/n)^(1/5) # Silverman
  for(i in 1:n)
  {
    H[[i]] <- h^2 * VAR[,i] %o% diag(2)

    # velocity integral correction (Fleming et al. Ecoinformatics appendix)
    if(length(CTMM[[i]]$tau)>1)
    {
      v <- c("vx","vy")
      v <- get.telemetry(data[[i]],axes=v)
      v <- vapply(1:nrow(data[[i]]),function(t){v[t,] %o% v[t,]},diag(2))
      v <- aperm(v,c(3,1,2))
      H[[i]] <- H[[i]] + dt[i]^2/12 * v
    }
  } # for(i in 1:n)

  W <- unlist(W)

  COL <- c("t","x","y")
  data <- lapply(data,function(d){d[,COL]})
  data <- do.call(rbind,data)

  H <- lapply(H,function(h){dim(h)=c(dim(h)[1],4);h})
  H <- do.call(rbind,H)
  dim(H) <- c(dim(H)[1],2,2)

  m <- length(W)
  # determine desired (absolute) resolution from smallest of all individual bandwidths
  dr <- sapply(1:m,function(i){sqrt(diag(H[i,,]))/res.space}) # (axes,individuals)
  dim(dr) <- c(2,m)
  dr <- apply(dr,1,min)

  KDE <- kde(data,H=H,W=W,res=res.space,grid=grid,dr=dr)

  KDE <- new.UD(KDE,info=info,type='range',variable="utilization",CTMM=ctmm())

  return(KDE)
}
