# second attempt
ridges2.UD <- function(object,precision=1/8,...)
{
  tol <- .Machine$double.eps^precision
  VOTE.MIN <- 2
  LENGTH.MIN <- 2

  PDF <- object$PDF
  PDF <- log(PDF) # log transform for more accurate derivative estimation

  dx <- object$dr['x']
  dy <- object$dr['y']
  dx2 <- dx^2
  dy2 <- dy^2
  dr2 <- dx^2+dy^2
  dxdydr2 <- dx*dy/dr2

  DIM <- dim(PDF)
  VOTE <- array(0,DIM)

  X <- object$r$x
  Y <- object$r$y
  dR <- c(dx,dy)
  dK <- 1/dR

  if(all(c('grad','hess') %in% names(object)))
  {
    GRAD <- object$grad
    HESS <- object$hess
  }
  else
  {
    # gradient information (under log)
    GRAD <- array(0,c(DIM,2)) # akima won't tolerate -Inf nor NA
    HESS <- array(0,c(DIM,2,2))

    ROWS <- 2:(nrow(PDF)-1)
    for(i in ROWS)
    {
      COLS <- 2:(ncol(PDF)-1)
      for(j in COLS)
      {
        SUB <- PDF[i+(-1):1,j+(-1):1]
        TEST <- all(SUB>-Inf)
        if(TEST) # non-zero density (under logarithm)
        {
          DIFF <- DiffGrid(SUB,dx,dy,dx2,dy2,dr2,dxdydr2)
          # these are better estimates than from splines
          GRAD[i,j,] <- DIFF$GRAD # 1/p D.p # from log
          HESS[i,j,,] <- DIFF$HESS # 1/p D.D.p - 1/p^2 D.p D.p # from log
          HESS[i,j,,] <- HESS[i,j,,] + (GRAD[i,j,] %o% GRAD[i,j,]) # logarithm correction # now 1/p D.D.p
          # last part is unnecessary for modes, but for ridges this looks important
        } # end if non-zero density
      } # end col loop
    } # end row loop
  }

  # where to store ridge metric
  RIDGE <- array(NA,dim(PDF))

  ROWS <- 2:(nrow(PDF)-1)
  for(i in ROWS)
  {
    COLS <- 2:(ncol(PDF)-1)
    for(j in COLS)
    {
      if(PDF[i,j] == -Inf) { next }
      EIGEN <- eigen(HESS[i,j,,])
      if(EIGEN$values[2]>=0) { next }
      # number of cells to ridge
      RIDGE[i,j] <- sqrt(sum((-(dK%*%EIGEN$vectors[,2])/EIGEN$values[2]*(GRAD[i,j,]%*%EIGEN$vectors[,2]))^2))
    } # end col loop
  } # end row loop

  return(RIDGE)

  ## connect ridge points from top to bottom ##
  # TODO !!! enforce minimum connection length - or make new branch !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LINE <- list()
  for(i in 1:nrow(RIDGE))
  {
    if(MODE[i]>0) # start new ridge-line from mode
    { LINE[[length(LINE)+1]] <- RIDGE[i,,drop=FALSE] }
    else # point should connect to existing ridge-line
    {
      r <- RIDGE[i,c('x','y')]
      DIST <- sapply(LINE,function(L){sum((r-tail(L,1))^2)})
      j <- which.min(DIST) # ridge-line to connect to

      # greedy search up tail end
      n <- nrow(LINE[[j]])
      DIST <- DIST[j]
      for(k in n:1)
      {
        if(k==n) { next }
        D.NEXT <- sum((r-LINE[[j]][k,])^2)
        if(DIST<=D.NEXT)
        {
          k <- k+1 # last point was best point
          break
        }
        DIST <- D.NEXT
      }

      DIST <- abs(r-LINE[[j]][k,])
      if(any(DIST>2*dR)) # far from other points
      { LINE[[length(LINE)+1]] <- RIDGE[i,,drop=FALSE] }
      if(k==n) # connect to tail end
      { LINE[[j]] <- rbind(LINE[[j]],RIDGE[i,]) }
      else # branch off and make new ridge
      { LINE[[length(LINE)+1]] <- rbind(LINE[[j]][k,],RIDGE[i,]) }
    } # end connect point to existing ridge
  } # end connect ridge points
  rm(RIDGE)

  ## connect end points into loops
  n <- length(LINE)
  LOOP <- list()

  if(n>1) # connect loops
  {
    DIST2 <- array(Inf,c(n,n,2))
    for(i in 1%:%(n-1))
    {
      for(j in (i+1)%:%n)
      { DIST2[i,j,] <- DIST2[j,i,] <- abs( tail(LINE[[i]],1) - tail(LINE[[j]],1) ) }
    }
    DIST <- DIST2[,,1]^2+DIST2[,,2]^2

    while(dim(DIST)[1]>1)
    {
      IND <- which.min(DIST)
      IND <- arrayInd(IND,dim(DIST))

      if(all(DIST2[IND[1],IND[2],]<=2*dR))
      {
        m <- nrow(LINE[[IND[1]]])
        LOOP[[1+length(LOOP)]] <- rbind( LINE[[IND[2]]] , LINE[[IND[1]]][m:1,] )
        LINE <- LINE[-c(IND)]
        DIST2 <- DIST2[-c(IND),-c(IND),,drop=FALSE]
        DIST <- DIST[-c(IND),-c(IND),drop=FALSE]
      }
      else
      { break }
    }
  } # end multiple ridges

  LINE <- c(LINE,LOOP)
  rm(LOOP)

  # LENGTH <- sapply(LINE,nrow)
  # LINE <- LINE[LENGTH>=LENGTH.MIN]

  LINE <- lapply(LINE,function(L){new.telemetry(data.frame(cbind(L,t=1:nrow(L))))})

  close(pb)
  return(LINE)
}

# first attempt (heavily revised)
ridges.UD <- function(object,...)
{
  # if(all(c('grad','hess') %nin% names(object))) { stop("Please run akde() with grad=TRUE") }

  LOG <- log(object$PDF) # log transform
  POINT <- array(0,dim(LOG))
  CURVE <- array(0,dim(LOG))
  #COST <- array(Inf,dim(LOG))

  dx <- object$dr['x']
  dy <- object$dr['y']
  dx2 <- dx^2
  dy2 <- dy^2
  dr2 <- dx^2+dy^2
  dxdydr2 <- dx*dy/dr2

  ROWS <- 2:(nrow(LOG)-1)
  # ROWS <- 1:nrow(LOG)
  for(i in ROWS)
  {
    COLS <- 2:(ncol(LOG)-1)
    # COLS <- 1:ncol(LOG)
    for(j in COLS)
    {
      SUB <- LOG[i+(-1):1,j+(-1):1]
      TEST <- all(SUB>-Inf)
      # TEST <- LOG[i,j]>0
      if(TEST) # non-zero density (under logarithm)
      {
        # GRAD <- object$grad[i,j,]
        # HESS <- object$hess[i,j,,]

        DIFF <- DiffGrid(SUB,dx,dy,dx2,dy2,dr2,dxdydr2)
        GRAD <- DIFF$GRAD # 1/p D.p # from log
        HESS <- DIFF$HESS # 1/p D.D.p - 1/p^2 D.p D.p # from log
        HESS <- HESS + (GRAD %o% GRAD) # logarithm correction # now 1/p D.D.p

        EIGEN <- eigen(HESS)
        if(EIGEN$values[2]<0) # most negative eigenvalue is negative # a ridge exists somewhere
        {
          # COST[i,j] <- (GRAD %*% EIGEN$vectors[,2])^2 / ( -EIGEN$values[2] )
          # if(EIGEN$values[2]>0) { COST[i,j] <- COST[i,j] + EIGEN$values[1]/EIGEN$values[2] }

          # ridge constraint
          # (GRAD + HESS %*% DL ) %*% EIGEN$vectors[,2] == 0

          # 1. Lagrangian (minimal distance)
          # L == 1/2 * DL %*% DL - lambda * (GRAD + HESS %*% DL ) %*% EIGEN$vectors[,2]
          # 0 == DL - lambda * HESS %*% EIGEN$vectors[,2]
          # 0 == DL - lambda * EIGEN$values[2] * EIGEN$vectors[,2]
          # DL == lambda * EIGEN$values[2] * EIGEN$vectors[,2]

          # 2. climbing directly up to the ridge
          # DL == dl * EIGEN$vectors[,2]
          # GRAD %*% EIGEN$vectors[,2] + dl*(EIGEN$vectors[,2] %*% HESS %*% EIGEN$vectors[,2]) == 0
          # GRAD %*% EIGEN$vectors[,2] + EIGEN$values[2]*dl == 0

          dr <- -c(GRAD %*% EIGEN$vectors[,2])/EIGEN$values[2]
          dr <- dr * EIGEN$vectors[,2]

          # pixel displacement
          dij <- dr/c(dx,dy)
          ip <- i+dij[1]
          jp <- j+dij[2]

          # is point in an adjacent pixel or this pixel
          if(max(abs(dij))<=2)
          {
            # point is in this pixel
            if(max(abs(dij))<=1)
            {
              POINT <- anti.alias(ip,jp,POINT,1)

              # what pixels are next?
              u <- EIGEN$vectors[,1] # bi-direction of ridge
              u <- u/c(dx,dy) # now in units of pixels
              u <- u/sqrt(sum(u^2)) # one pixel length

              POINT <- anti.alias(ip+u[1],jp+u[2],POINT,1/2)
              POINT <- anti.alias(ip-u[1],jp-u[2],POINT,1/2)

              CURVE <- anti.alias(ip+u[1],jp+u[2],CURVE,-EIGEN$values[2]/2)
              CURVE <- anti.alias(ip-u[1],jp-u[2],CURVE,-EIGEN$values[2]/2)
            } # point is in this pixel
            else # point is in adjacent pixel
            {
              POINT <- anti.alias(ip,jp,POINT,1/8)

              CURVE <- anti.alias(ip,jp,CURVE,-EIGEN$values[2]/8)
            }
          } # point is in nearby pixel
        } # end ridge exists
      } # end if non-zero density
    } # end col loop
  } # end row loop

  # CTMM <- object@CTMM
  # DOF <- DOF.mean(CTMM)
  # COST <- PDF * COST # cancel 1/p term
  # COST <- DOF* COST # now like grad^2/std.err^2 ~ F

  # normalize all rasters
  MAX <- 1+2/2+6/8
  POINT <- POINT/MAX
  CURVE <- CURVE/MAX

  MAX <- pmax(POINT,1)
  POINT <- POINT/MAX
  CURVE <- CURVE/MAX

  # average ridge curvature
  SUB <- CURVE>0
  WIDTH <- sum(CURVE[SUB] * object$PDF[SUB]) / sum(object$PDF[SUB])
  # average ridge width
  WIDTH <- 1/sqrt(WIDTH)

  # curvature -> width
  #CURVE[SUB] <- 1/sqrt(CURVE[SUB])
  #CURVE[!SUB] <- 0

  RETURN <- list(Indicator=POINT,Ave.Width=WIDTH)
  return(RETURN)
}


anti.alias <- function(i,j,M,dM=1)
{
  I <- unique( c(floor(i),ceiling(i)) )
  J <- unique( c(floor(j),ceiling(j)) )

  I <- I[0<I & I<nrow(M)]
  J <- J[0<J & J<ncol(M)]

  dM <- array(dM,c(length(I),length(J)))

  if(!length(dM)) { return(M) }

  if(length(I)>1)
  {
    w <- abs(i-I[1])
    dM[1,] <- w*dM[1,]
    dM[2,] <- (1-w)*dM[2,]
  }

  if(length(J)>1)
  {
    w <- abs(j-J[1])
    dM[,1] <- w*dM[,1]
    dM[,2] <- (1-w)*dM[,2]
  }

  M[I,J] <- M[I,J] + dM
  return(M)
}
